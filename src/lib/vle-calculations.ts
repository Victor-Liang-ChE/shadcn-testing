// Consolidated VLE calculation utilities supporting all major thermodynamic models:
// • NRTL (Non-Random Two Liquid) - activity coefficient model
// • Wilson - activity coefficient model  
// • UNIFAC (UNIversal Functional Activity Coefficient) - group contribution method
// • Peng-Robinson (PR) - cubic equation of state with fugacity coefficients
// • Soave-Redlich-Kwong (SRK) - cubic equation of state with fugacity coefficients
// • UNIQUAC (UNIversal QUAsi Chemical) - activity coefficient model
// All legacy VLE calculation files have been consolidated into this single module.
// ---------------------------------------------------------------------------------------
import type { AntoineParams, CompoundData, BubbleDewResult } from './vle-types';
import type { SupabaseClient } from '@supabase/supabase-js';
import { fetchCompoundDataFromHysys } from './antoine-utils';

// =============================
//  1. SHARED CONSTANTS/UTILS
// =============================
export const R_gas_const_J_molK = 8.314_462_618_153_24; // J·mol⁻¹·K⁻¹
// SRK debug flag — set NEXT_PUBLIC_SRK_DEBUG=true to enable development logging for SRK routines
export const SRK_DEBUG = process.env.NEXT_PUBLIC_SRK_DEBUG === 'true' || false;

/**
 * Helper function to estimate boiling point temperature from Aspen PLXANT equation.
 * Uses Newton-Raphson on: ln(P [Pa]) = C1 + C2/(T+C3) + C4·T + C5·ln(T) + C6·T^C7
 * T in Kelvin, P_target in Pa.
 */
export function antoineBoilingPointSolverLocal(antoineParams: AntoineParams | null, P_target: number): number | null {
  if (!antoineParams) return null;
  try {
    const targetLnP = Math.log(P_target); // target is already in Pa
    // Start guess: invert the dominant C1 + C2/T term
    let T = antoineParams.C2 / (targetLnP - antoineParams.C1);
    if (!isFinite(T) || T <= 0 || T > 2000) T = 350;

    const { C1, C2, C3, C4, C5, C6, C7 } = antoineParams;
    for (let i = 0; i < 30; i++) {
      const lnP = C1 + C2 / (T + C3) + C4 * T + C5 * Math.log(T) + C6 * Math.pow(T, C7);
      const err = lnP - targetLnP;
      if (Math.abs(err) < 1e-6) return T;
      // d(lnP)/dT = -C2/(T+C3)^2 + C4 + C5/T + C6·C7·T^(C7-1)
      const dlnPdT = -C2 / (T + C3) ** 2 + C4 + C5 / T + C6 * C7 * Math.pow(T, C7 - 1);
      if (Math.abs(dlnPdT) < 1e-15) break;
      T -= err / dlnPdT;
      if (T <= 0 || T > 2000) break;
    }
    return (T > 0 && T < 2000) ? T : null;
  } catch {
    return null;
  }
}

/**
 * Calculates saturation pressure (Pa) from Aspen PLXANT extended Antoine.
 *   ln(P* [Pa]) = C1 + C2/(T+C3) + C4·T + C5·ln(T) + C6·T^C7
 * T in Kelvin. Returns pressure in Pa directly.
 */
export function calculatePsat_Pa(
  params: AntoineParams | null,
  T_K: number
): number {
  if (!params || T_K <= 0) return NaN;
  const { C1, C2, C3, C4, C5, C6, C7 } = params;
  const lnP_Pa = C1 + C2 / (T_K + C3) + C4 * T_K + C5 * Math.log(T_K) + C6 * Math.pow(T_K, C7);
  const P_Pa = Math.exp(lnP_Pa);
  return isNaN(P_Pa) || P_Pa <= 0 ? NaN : P_Pa;
}

// ==========================================
//  2. N R T L   (binary activity-coefficient)
// ==========================================
/**
 * NRTL interaction params – HYSYS convention:
 *   τ_ij = Aij/T + Bij   (Aij in K, Bij dimensionless)
 *   G_ij = exp(-alpha · τ_ij)
 */
export interface NrtlInteractionParams {
  Aij: number;   // K (= bij*R_cal)
  Aji: number;   // K (= bji*R_cal)
  Bij: number;   // dimensionless (= aij*R_cal)
  Bji: number;   // dimensionless (= aji*R_cal)
  alpha: number; // non-randomness parameter (typically 0.2-0.47)
  /** Aspen dij: alpha temperature coefficient — α(T) = alpha + dij*(T−273.15). Default 0. */
  dij_alpha?: number;
  /** Aspen eij·R_cal: ln T contribution to τ_ij = … + eij·ln T. Default 0. */
  Eij?: number;
  /** Aspen eji·R_cal: ln T contribution to τ_ji. Default 0. */
  Eji?: number;
  /** Aspen fij·R_cal: T contribution to τ_ij = … + fij·T. Default 0. */
  Fij?: number;
  /** Aspen fji·R_cal: T contribution to τ_ji. Default 0. */
  Fji?: number;
  /** Hayden-O'Connell cross-association parameter (η_ij). Default 0 (ideal vapor). */
  eta_cross?: number;
}

export function calculateNrtlGamma(
  _components: CompoundData[],
  x: number[],
  T_K: number,
  p: NrtlInteractionParams
): [number, number] | null {
  if (x.length !== 2) return null;
  const [x1, x2] = x;
  if (Math.abs(x1 + x2 - 1) > 1e-6) return null;
  if (x1 < 1e-9 || x2 < 1e-9) return [1, 1];

  // Aspen NRTL parameters: Aij = bij*R_cal, Bij = aij*R_cal (see fetchNrtlParameters).
  // Full Aspen tau: τ_ij = aij + bij/T + eij·ln T + fij·T
  // Stored as: Aij/(R*T) + Bij/R + Eij·lnT/R + Fij·T/R = bij/T + aij + eij·lnT + fij·T
  // Alpha temp-dep: α(T) = alpha + dij_alpha·(T-273.15)
  const R_cal = 1.9872; // cal·mol⁻¹·K⁻¹
  const alpha_T = p.alpha + (p.dij_alpha ?? 0) * (T_K - 273.15);
  // tau12: i=1,j=2 uses Aij/Bij/Eij/Fij; tau21: i=2,j=1 uses Aji/Bji/Eji/Fji
  const tau12 = p.Aij / (R_cal * T_K) + p.Bij / R_cal
              + (p.Eij ?? 0) * Math.log(T_K) / R_cal
              + (p.Fij ?? 0) * T_K / R_cal;
  const tau21 = p.Aji / (R_cal * T_K) + p.Bji / R_cal
              + (p.Eji ?? 0) * Math.log(T_K) / R_cal
              + (p.Fji ?? 0) * T_K / R_cal;
  const G12 = Math.exp(-alpha_T * tau12);
  const G21 = Math.exp(-alpha_T * tau21);

  const denom1 = x1 + x2 * G21;
  const denom2 = x2 + x1 * G12;
  if (denom1 === 0 || denom2 === 0) return null;

  const lnGamma1 =
    x2 ** 2 * (tau21 * (G21 / denom1) ** 2 + tau12 * G12 / denom2 ** 2);
  const lnGamma2 =
    x1 ** 2 * (tau12 * (G12 / denom2) ** 2 + tau21 * G21 / denom1 ** 2);

  return [Math.exp(lnGamma1), Math.exp(lnGamma2)];
}

export function calculateBubbleTemperatureNrtl(
  comps: CompoundData[],
  x1_feed: number,
  P_Pa: number,
  params: NrtlInteractionParams,
  T_init_K = 350,
  maxIter = 50,
  tol = 1e-4
): BubbleDewResult | null {
  let T = T_init_K;
  const x = [x1_feed, 1 - x1_feed];
  let y: number[] = [...x];
  const eta_cross = (params as any).eta_cross ?? 0;
  const useHoc = eta_cross > 0 ||
    (comps[0].hocProps?.eta_self ?? 0) > 0 ||
    (comps[1].hocProps?.eta_self ?? 0) > 0;

  for (let i = 0; i < maxIter; i++) {
    const Psat = [calculatePsat_Pa(comps[0].antoine, T), calculatePsat_Pa(comps[1].antoine, T)];
    if (Psat.some(v => isNaN(v))) return null;
    const gamma = calculateNrtlGamma(comps, x, T, params);
    if (!gamma) return null;

    const phi = useHoc
      ? hocFugacityCorrection(comps, y, T, P_Pa, Psat, eta_cross)
      : ([1, 1] as [number, number]);

    const KPsat = [x[0] * gamma[0] * Psat[0] * phi[0], x[1] * gamma[1] * Psat[1] * phi[1]];
    const Pcalc = KPsat[0] + KPsat[1];
    const err = Pcalc - P_Pa;
    if (Math.abs(err) < tol * P_Pa) {
      y = [KPsat[0] / Pcalc, KPsat[1] / Pcalc];
      return {
        comp1_feed: x1_feed,
        comp1_equilibrium: y[0],
        T_K: T,
        P_Pa,
        iterations: i + 1,
        calculationType: 'bubbleT',
      };
    }
    const Ptmp = Pcalc || P_Pa;
    y = [KPsat[0] / Ptmp, KPsat[1] / Ptmp];
    // simple secant-like adjustment
    const dT = -err / P_Pa * T * 0.02;
    T = Math.max(100, T + (Math.abs(dT) > 5 ? Math.sign(dT) * 5 : dT));
  }
  return null;
}

export function calculateBubblePressureNrtl(
  comps: CompoundData[],
  x1_feed: number,
  T_K: number,
  params: NrtlInteractionParams
): BubbleDewResult | null {
  const x = [x1_feed, 1 - x1_feed];
  const Psat = [calculatePsat_Pa(comps[0].antoine, T_K), calculatePsat_Pa(comps[1].antoine, T_K)];
  if (Psat.some(v => isNaN(v))) return null;
  const gamma = calculateNrtlGamma(comps, x, T_K, params);
  if (!gamma) return null;

  const eta_cross = (params as any).eta_cross ?? 0;
  const useHoc = eta_cross > 0 ||
    (comps[0].hocProps?.eta_self ?? 0) > 0 ||
    (comps[1].hocProps?.eta_self ?? 0) > 0;

  const PcalcIdeal = x[0] * gamma[0] * Psat[0] + x[1] * gamma[1] * Psat[1];
  const yIdeal = [x[0] * gamma[0] * Psat[0] / PcalcIdeal, x[1] * gamma[1] * Psat[1] / PcalcIdeal];
  const phi = useHoc
    ? hocFugacityCorrection(comps, yIdeal, T_K, PcalcIdeal, Psat, eta_cross)
    : ([1, 1] as [number, number]);

  const KPsat = [x[0] * gamma[0] * Psat[0] * phi[0], x[1] * gamma[1] * Psat[1] * phi[1]];
  const Pcalc = KPsat[0] + KPsat[1];
  const y1 = KPsat[0] / Pcalc;
  return {
    comp1_feed: x1_feed,
    comp1_equilibrium: y1,
    T_K,
    P_Pa: Pcalc,
    iterations: 1,
    calculationType: 'bubbleP',
  };
}

// ==========================================
//  3. W I L S O N   (binary activity model)
// ==========================================
/**
 * Wilson interaction params – HYSYS convention:
 *   The database stores a_ij = (λ_ij − λ_jj), so fetchWilsonInteractionParams
 *   swaps Aij↔Aji to obtain a_12 = (λ_12 − λ_11) which the Wilson equation needs.
 *   Aij in cal/mol, Bij dimensionless.
 */
export interface WilsonInteractionParams {
  Aij: number;     // cal/mol (= −bij·R_cal), energy param for Λ_12
  Aji: number;     // cal/mol (= −bji·R_cal), energy param for Λ_21
  Bij: number;     // dimensionless (= −aij·R_cal)
  Bji: number;     // dimensionless (= −aji·R_cal)
  /** Aspen cij·(−R_cal): ln T coefficient of ln Λ_ij = … + cij·ln T. Default 0. */
  Cij?: number;
  /** Aspen cji·(−R_cal): ln T coefficient of ln Λ_ji. Default 0. */
  Cji?: number;
  /** Aspen dij·(−R_cal): T coefficient of ln Λ_ij = … + dij·T. Default 0. */
  Dij?: number;
  /** Aspen dji·(−R_cal): T coefficient of ln Λ_ji. Default 0. */
  Dji?: number;
  /** Hayden-O'Connell cross-association parameter (η_ij). Default 0 (ideal vapor). */
  eta_cross?: number;
}

export function calculateWilsonGamma(
  comps: CompoundData[],
  x: number[],
  T_K: number,
  p: WilsonInteractionParams
): [number, number] | null {
  if (x.length !== 2 || comps.length !== 2) return null;
  if (!comps[0].wilsonParams || !comps[1].wilsonParams) return null;

  // Use temperature-dependent Rackett liquid molar volume when critical props are available.
  // VL(T) = (R·Tc/Pc)·ZRA^[1+(1-Tr)^(2/7)], ZRA = 0.29056 - 0.08775·ω (Yamada-Gunn)
  // This is consistent with how APV Wilson params were regressed in Aspen Plus.
  const _rackett = (wp: typeof comps[0]['wilsonParams'], T: number): number => {
    if (wp!.Tc_K && wp!.Pc_Pa) {
      const ZRA = wp!.RKTZRA ?? (wp!.omega != null ? 0.29056 - 0.08775 * wp!.omega! : null);
      if (ZRA != null) {
        const Tr = T / wp!.Tc_K!;
        return (R_gas_const_J_molK * wp!.Tc_K! / wp!.Pc_Pa!) * Math.pow(ZRA, 1 + Math.pow(1 - Tr, 2 / 7));
      }
    }
    return wp!.V_L_m3mol; // fallback to stored reference value
  };
  const V1 = _rackett(comps[0].wilsonParams, T_K);
  const V2 = _rackett(comps[1].wilsonParams, T_K);
  // Standard convention: Lambda12 uses Aij/Bij (i=1,j=2), Lambda21 uses Aji/Bji
  const R_cal = 1.9872; // cal·mol⁻¹·K⁻¹
  // Scale Bij by R_cal
  // Full Aspen Wilson: ln Λ_ij = ln(Vj/Vi) + aij + bij/T + cij·lnT + dij·T
  // Stored scaled: -(Aij + Bij·T + Cij·T·lnT + Dij·T²) / (R_cal·T)
  //   = bij/T + aij + cij·lnT + dij·T  ✓  (Cij = -cij·R_cal, Dij = -dij·R_cal)
  const Lambda12 = (V2 / V1) * Math.exp(
    -(p.Aij + p.Bij * T_K + (p.Cij ?? 0) * T_K * Math.log(T_K) + (p.Dij ?? 0) * T_K * T_K) / (R_cal * T_K)
  );
  const Lambda21 = (V1 / V2) * Math.exp(
    -(p.Aji + p.Bji * T_K + (p.Cji ?? 0) * T_K * Math.log(T_K) + (p.Dji ?? 0) * T_K * T_K) / (R_cal * T_K)
  );

  const [x1, x2] = x;
  const denom1 = x1 + Lambda12 * x2;
  const denom2 = x2 + Lambda21 * x1;
  if (denom1 === 0 || denom2 === 0) return null;

  const lnGamma1 = -Math.log(denom1) + x2 * (Lambda12 / denom1 - Lambda21 / denom2);
  const lnGamma2 = -Math.log(denom2) - x1 * (Lambda12 / denom1 - Lambda21 / denom2);
  return [Math.exp(lnGamma1), Math.exp(lnGamma2)];
}

export function calculateBubbleTemperatureWilson(
  comps: CompoundData[],
  x1_feed: number,
  P_Pa: number,
  params: WilsonInteractionParams,
  T_init_K: number,
  maxIter = 50,
  tol = 1e-5
): BubbleDewResult | null {
  let T = T_init_K;
  const x = [x1_feed, 1 - x1_feed];
  let y: number[] = [...x]; // initial estimate
  const eta_cross = params.eta_cross ?? 0;
  const useHoc = eta_cross > 0 ||
    (comps[0].hocProps?.eta_self ?? 0) > 0 ||
    (comps[1].hocProps?.eta_self ?? 0) > 0;

  for (let i = 0; i < maxIter; i++) {
    const gammas = calculateWilsonGamma(comps, x, T, params);
    if (!gammas) return null;
    const Psat = [calculatePsat_Pa(comps[0].antoine, T), calculatePsat_Pa(comps[1].antoine, T)];
    if (Psat.some(v => isNaN(v))) return null;

    // HOC vapor-phase fugacity correction (φ_i^sat / φ_i^V)
    const phi = useHoc
      ? hocFugacityCorrection(comps, y, T, P_Pa, Psat, eta_cross)
      : ([1, 1] as [number, number]);

    // Modified Raoult's law: K_i = γ_i · Psat_i · φ_i^sat / (φ_i^V · P)
    const KPsat = [x[0] * gammas[0] * Psat[0] * phi[0], x[1] * gammas[1] * Psat[1] * phi[1]];
    const Pcalc = KPsat[0] + KPsat[1];
    const err = Pcalc - P_Pa;
    if (Math.abs(err) < tol * P_Pa) {
      y = [KPsat[0] / Pcalc, KPsat[1] / Pcalc];
      return {
        comp1_feed: x1_feed,
        comp1_equilibrium: y[0],
        T_K: T,
        P_Pa,
        iterations: i + 1,
        calculationType: 'bubbleT',
      };
    }
    // update y estimate before next iteration (for HOC correction)
    const Ptmp = Pcalc || P_Pa;
    y = [KPsat[0] / Ptmp, KPsat[1] / Ptmp];
    const dT = -err / P_Pa * (T * 0.05);
    T = Math.min(800, Math.max(150, T + (Math.abs(dT) > 5 ? Math.sign(dT) * 5 : dT)));
  }
  return null;
}

export function calculateBubblePressureWilson(
  comps: CompoundData[],
  x1_feed: number,
  T_K: number,
  params: WilsonInteractionParams
): BubbleDewResult | null {
  const x = [x1_feed, 1 - x1_feed];
  const gammas = calculateWilsonGamma(comps, x, T_K, params);
  if (!gammas) return null;
  const Psat = [calculatePsat_Pa(comps[0].antoine, T_K), calculatePsat_Pa(comps[1].antoine, T_K)];
  if (Psat.some(v => isNaN(v))) return null;

  const eta_cross = params.eta_cross ?? 0;
  const useHoc = eta_cross > 0 ||
    (comps[0].hocProps?.eta_self ?? 0) > 0 ||
    (comps[1].hocProps?.eta_self ?? 0) > 0;

  // Initial ideal estimate for y (needed for HOC correction)
  const PcalcIdeal = x[0] * gammas[0] * Psat[0] + x[1] * gammas[1] * Psat[1];
  let y = [x[0] * gammas[0] * Psat[0] / PcalcIdeal, x[1] * gammas[1] * Psat[1] / PcalcIdeal];

  const phi = useHoc
    ? hocFugacityCorrection(comps, y, T_K, PcalcIdeal, Psat, eta_cross)
    : ([1, 1] as [number, number]);

  const KPsat = [x[0] * gammas[0] * Psat[0] * phi[0], x[1] * gammas[1] * Psat[1] * phi[1]];
  const Pcalc = KPsat[0] + KPsat[1];
  const y1 = KPsat[0] / Pcalc;
  return {
    comp1_feed: x1_feed,
    comp1_equilibrium: y1,
    T_K,
    P_Pa: Pcalc,
    iterations: 1,
    calculationType: 'bubbleP',
  };
}

// =====================================================================
//  4. U N I F A C   (group-contribution activity model – binary support)
// =====================================================================

export interface UnifacParameters {
  Rk: { [subgroupId: number]: number };
  Qk: { [subgroupId: number]: number };
  mainGroupMap: { [subgroupId: number]: number };
  a_mk: Map<string, number>; // key `${mainGroup_m}-${mainGroup_k}`
}

/**
 * Calculates UNIFAC activity coefficients using the Larsen & Gmehling formulation.
 * The implementation is adapted from the previous standalone file and supports binary mixtures.
 */
export function calculateUnifacGamma(
  componentData: CompoundData[],
  x: number[],
  T_kelvin: number,
  params: UnifacParameters
): number[] | null {
  if (componentData.length === 0 || x.length !== componentData.length) return null;
  const Z = 10.0; // Coordination number

  // --- Pre-calculate r_i and q_i for every component (only once) ---
  for (const c of componentData) {
    if (c.r_i === undefined || c.q_i === undefined) {
      if (!c.unifacGroups) return null;
      let rSum = 0, qSum = 0;
      for (const [subIdStr, count] of Object.entries(c.unifacGroups)) {
        const sgId = parseInt(subIdStr, 10);
        if (params.Rk[sgId] == null || params.Qk[sgId] == null) return null;
        rSum += count * params.Rk[sgId];
        qSum += count * params.Qk[sgId];
      }
      c.r_i = rSum;
      c.q_i = qSum;
    }
  }

  const r = componentData.map(c => c.r_i!);
  const q = componentData.map(c => c.q_i!);

  // -------- Combinatorial part (ln γ^C) --------
  const sum_xr = x.reduce((s, xi, i) => s + xi * r[i], 0);
  const sum_xq = x.reduce((s, xi, i) => s + xi * q[i], 0);
  if (sum_xr === 0 || sum_xq === 0) return null;

  const phi = x.map((xi, i) => (xi * r[i]) / sum_xr);
  const theta = x.map((xi, i) => (xi * q[i]) / sum_xq);
  const L = r.map((ri, i) => (Z / 2) * (ri - q[i]) - (ri - 1));
  const sum_xL = x.reduce((s, xi, i) => s + xi * L[i], 0);

  const ln_gamma_C = x.map((xi, i) =>
    Math.log(phi[i] / xi) + (Z / 2) * q[i] * Math.log(theta[i] / phi[i]) + L[i] - (phi[i] / xi) * sum_xL
  );

  // -------- Residual part (ln γ^R) --------
  const subgroupIds = Array.from(new Set(componentData.flatMap(c => Object.keys(c.unifacGroups || {}).map(Number))));
  if (subgroupIds.length === 0) return ln_gamma_C.map(Math.exp);

  const v_ki = componentData.map(c => subgroupIds.map(id => c.unifacGroups?.[id] || 0));
  
  // --- Calculate ln Γ_k for the MIXTURE ---
  const sum_x_vk_total = subgroupIds.reduce((sum, _, kIdx) => sum + componentData.reduce((s, _, i) => s + x[i] * v_ki[i][kIdx], 0), 0);
  if (sum_x_vk_total === 0) return ln_gamma_C.map(Math.exp);

  const X_m = subgroupIds.map((_, kIdx) => componentData.reduce((s, _, i) => s + x[i] * v_ki[i][kIdx], 0) / sum_x_vk_total);
  const sum_XQ = subgroupIds.reduce((s, id, idx) => s + X_m[idx] * params.Qk[id], 0);
  if (sum_XQ === 0) return ln_gamma_C.map(Math.exp);
  
  const theta_m = subgroupIds.map((id, idx) => (X_m[idx] * params.Qk[id]) / sum_XQ);
  const psi = subgroupIds.map((id_m, _mIdx) => subgroupIds.map((id_k, _kIdx) => {
      const key = `${params.mainGroupMap[id_m]}-${params.mainGroupMap[id_k]}`;
      const a_mk = params.a_mk.get(key) ?? 0;
      return Math.exp(-a_mk / T_kelvin);
  }));

  const calculate_ln_Gamma = (current_theta_m: number[]): number[] => {
      return subgroupIds.map((_, kIdx) => {
          const term1_sum = subgroupIds.reduce((s, _, mIdx) => s + current_theta_m[mIdx] * psi[mIdx][kIdx], 0);
          if (term1_sum === 0) return 0;
          const term1 = Math.log(term1_sum);
          const term2_sum = subgroupIds.reduce((s, _, mIdx) => {
              const denom = subgroupIds.reduce((sd, _, nIdx) => sd + current_theta_m[nIdx] * psi[nIdx][mIdx], 0);
              return s + (current_theta_m[mIdx] * psi[kIdx][mIdx]) / (denom || 1);
          }, 0);
          return params.Qk[subgroupIds[kIdx]] * (1 - term1 - term2_sum);
      });
  };
  
  const ln_Gamma_k_mixture = calculate_ln_Gamma(theta_m);

  // --- Calculate ln Γ_k for each PURE component i ---
  const ln_Gamma_k_pure: number[][] = [];
  for (let i = 0; i < componentData.length; i++) {
      const v_k = subgroupIds.map((id) => componentData[i].unifacGroups?.[id] || 0);
      const sum_vQ = subgroupIds.reduce((s, id, idx) => s + v_k[idx] * params.Qk[id], 0);
      if (sum_vQ === 0) {
          ln_Gamma_k_pure.push(new Array(subgroupIds.length).fill(0));
          continue;
      }
      const theta_k_pure = subgroupIds.map((id, idx) => (v_k[idx] * params.Qk[id]) / sum_vQ);
      ln_Gamma_k_pure.push(calculate_ln_Gamma(theta_k_pure));
  }
  
  // --- CORRECTED Final summation for ln γ^R ---
  const ln_gamma_R = componentData.map((_, i) => {
    return subgroupIds.reduce((sum, _, kIdx) => {
        return sum + v_ki[i][kIdx] * (ln_Gamma_k_mixture[kIdx] - ln_Gamma_k_pure[i][kIdx]);
    }, 0);
  });

  return ln_gamma_C.map((lnC, i) => {
    const ln_g = lnC + ln_gamma_R[i];
    return isFinite(ln_g) ? Math.exp(ln_g) : NaN;
  });
}

export function calculateBubbleTemperatureUnifac(
  components: CompoundData[],
  x1_feed: number,
  P_system_Pa: number,
  unifacParams: UnifacParameters,
  initialTempGuess_K: number = 350,
  maxIter: number = 50,
  tolerance: number = 1e-4
): BubbleDewResult | null {
  const x = [x1_feed, 1 - x1_feed];
  let T = initialTempGuess_K;
  for (let i = 0; i < maxIter; i++) {
    const P_sat = [calculatePsat_Pa(components[0].antoine, T), calculatePsat_Pa(components[1].antoine, T)];
    if (P_sat.some(p => isNaN(p))) return null;
    const gamma = calculateUnifacGamma(components, x, T, unifacParams);
    if (!gamma) return null;
    const P_calc = x[0] * gamma[0] * P_sat[0] + x[1] * gamma[1] * P_sat[1];
    const err = P_calc - P_system_Pa;
    if (Math.abs(err) < tolerance * P_system_Pa) {
      const y1 = (x[0] * gamma[0] * P_sat[0]) / P_calc;
      return { comp1_feed: x1_feed, comp1_equilibrium: y1, T_K: T, P_Pa: P_system_Pa };
    }
    const dT = -err / P_system_Pa * (T * 0.05);
    T = Math.max(150, Math.min(800, T + dT));
  }
  return null;
}

export function calculateBubblePressureUnifac(
  components: CompoundData[],
  x1_feed: number,
  T_K: number,
  unifacParams: UnifacParameters
): BubbleDewResult | null {
  const x = [x1_feed, 1 - x1_feed];
  const Psat = [calculatePsat_Pa(components[0].antoine, T_K), calculatePsat_Pa(components[1].antoine, T_K)];
  if (Psat.some(p => isNaN(p))) return null;
  const gamma = calculateUnifacGamma(components, x, T_K, unifacParams);
  if (!gamma) return null;
  const Pcalc = x[0] * gamma[0] * Psat[0] + x[1] * gamma[1] * Psat[1];
  const y1 = (x[0] * gamma[0] * Psat[0]) / Pcalc;
  return { comp1_feed: x1_feed, comp1_equilibrium: y1, T_K, P_Pa: Pcalc };
}

// Minimal inline variants for dew calculations (binary only)
export function calculateDewTemperatureUnifac(
  components: CompoundData[],
  y1: number,
  P_Pa: number,
  unifacParams: UnifacParameters,
  T_guess_K: number = 350
): BubbleDewResult | null {
  // Simple iterative solver (could be improved)
  let T = T_guess_K;
  const y = [y1, 1 - y1];
  for (let i = 0; i < 60; i++) {
    const P_sat = [calculatePsat_Pa(components[0].antoine, T), calculatePsat_Pa(components[1].antoine, T)];
    if (P_sat.some(p => isNaN(p))) return null;
    const K = [P_sat[0] / P_Pa, P_sat[1] / P_Pa];
    const x = y.map((yi, idx) => yi / K[idx]);
    const x_sum = x[0] + x[1];
    x[0] /= x_sum; x[1] /= x_sum;
    const gamma = calculateUnifacGamma(components, x, T, unifacParams);
    if (!gamma) return null;
    const func = y[0] - x[0] * gamma[0] * P_sat[0] / P_Pa;
    if (Math.abs(func) < 1e-6) {
      return { comp1_feed: y1, comp1_equilibrium: x[0], T_K: T, P_Pa };
    }
    T += func * 2; // crude adjustment
  }
  return null;
}

export function calculateDewPressureUnifac(
  components: CompoundData[],
  y1: number,
  T_K: number,
  unifacParams: UnifacParameters,
  P_guess_Pa: number = 101325
): BubbleDewResult | null {
  // Iterative on pressure
  let P = P_guess_Pa;
  const y = [y1, 1 - y1];
  for (let i = 0; i < 40; i++) {
    const P_sat = [calculatePsat_Pa(components[0].antoine, T_K), calculatePsat_Pa(components[1].antoine, T_K)];
    if (P_sat.some(p => isNaN(p))) return null;
    const K = [P_sat[0] / P, P_sat[1] / P];
    const x = y.map((yi, idx) => yi / K[idx]);
    const x_sum = x[0] + x[1]; x[0] /= x_sum; x[1] /= x_sum;
    const gamma = calculateUnifacGamma(components, x, T_K, unifacParams);
    if (!gamma) return null;
    const P_calc = (y[0] * P) / (x[0] * gamma[0]) ;
    const err = P_calc - P;
    if (Math.abs(err) < 1e-4 * P) {
      return { comp1_feed: y1, comp1_equilibrium: x[0], T_K, P_Pa: P };
    }
    P += err * 0.5;
  }
  return null;
}

// ===============================================================
//  5. P E N G – R O B I N S O N   E O S   (binary – vapor/liquid)
// ===============================================================

export interface PrInteractionParams { k_ij: number; }

export function solveCubicEOS(a: number, b: number, c: number, d: number): number[] | null {
  // General cubic a·Z³ + b·Z² + c·Z + d = 0  (Cardano) – returns only positive real roots
  if (Math.abs(a) < 1e-12) return null;

  // Normalize to Z³ + p Z² + q Z + r = 0
  const p = b / a;
  const q = c / a;
  const r = d / a;

  // Depressed cubic t³ + A t + B = 0 with Z = t – p/3
  const A = q - (p * p) / 3;
  const B = (2 * p * p * p) / 27 - (p * q) / 3 + r;
  let D = (B * B) / 4 + (A * A * A) / 27; // Discriminant

  if (Math.abs(D) < 1e-20) D = 0; // treat ~0 as zero

  let roots: number[] = [];
  if (D > 0) {
    // One real root
    const S = Math.cbrt(-B / 2 + Math.sqrt(D));
    const T = Math.cbrt(-B / 2 - Math.sqrt(D));
    roots = [S + T - p / 3];
  } else {
    // Three real roots
    const rho = 2 * Math.sqrt(-A / 3);
    const theta = Math.acos(Math.max(-1, Math.min(1, (-B / 2) / Math.sqrt(-(A * A * A) / 27)))) / 3;
    roots = [
      rho * Math.cos(theta) - p / 3,
      rho * Math.cos(theta + (2 * Math.PI) / 3) - p / 3,
      rho * Math.cos(theta + (4 * Math.PI) / 3) - p / 3,
    ];
  }

  // Keep positive roots (Z > 0)
  const realPos = roots.filter(z => z > 1e-9).sort((x, y) => x - y);
  return realPos.length ? realPos : null;
}

export function calculatePrFugacityCoefficients(
  components: CompoundData[],
  x_or_y: number[],
  T_K: number,
  P_Pa: number,
  interactionParams: PrInteractionParams,
  phase: 'liquid' | 'vapor'
): [number, number] | null {
  if (components.length !== 2 || x_or_y.length !== 2) return null;
  if (!components[0].prParams || !components[1].prParams) return null;

  const [x1, x2] = x_or_y;
  const { k_ij } = interactionParams;

  const pure = components.map(c => {
    const { Tc_K, Pc_Pa, omega } = c.prParams!;
    const Tr = T_K / Tc_K;
    let kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    if (omega > 0.49) {
      kappa = 0.379642 + 1.48503 * omega - 0.164423 * (omega * omega) + 0.016666 * (omega * omega * omega);
    }
    const alpha = (1 + kappa * (1 - Math.sqrt(Tr))) ** 2;
    const ac = 0.457235 * (R_gas_const_J_molK * Tc_K) ** 2 / Pc_Pa * alpha;
    const b  = 0.07796 * R_gas_const_J_molK * Tc_K / Pc_Pa;
    return { ac, b }; 
  });

  const ac_mix = x1 * x1 * pure[0].ac + x2 * x2 * pure[1].ac + 2 * x1 * x2 * Math.sqrt(pure[0].ac * pure[1].ac) * (1 - k_ij);
  const b_mix  = x1 * pure[0].b + x2 * pure[1].b;

  const A = ac_mix * P_Pa / (R_gas_const_J_molK * T_K) ** 2;
  const B = b_mix * P_Pa / (R_gas_const_J_molK * T_K);

  const roots = solveCubicEOS(1, -(1 - B), A - 3 * B ** 2 - 2 * B, -(A * B - B ** 2 - B ** 3));
  if (!roots || roots.length === 0) return null;

  const Z = phase === 'liquid' ? Math.min(...roots) : Math.max(...roots);
  if (Z <= B) return null; // Robustness check

  const termCommon = (A / (2 * Math.sqrt(2) * B)) * Math.log((Z + (1 + Math.sqrt(2)) * B) / (Z + (1 - Math.sqrt(2)) * B));

  const d_a_mix = [
    2 * x1 * pure[0].ac + 2 * x2 * Math.sqrt(pure[0].ac * pure[1].ac) * (1 - k_ij),
    2 * x2 * pure[1].ac + 2 * x1 * Math.sqrt(pure[0].ac * pure[1].ac) * (1 - k_ij)
  ];

  const ln_phi = [0, 1].map(i => {
    const b_i = pure[i].b;
    const term1 = (b_i / b_mix) * (Z - 1) - Math.log(Z - B);
    const term2 = termCommon * (d_a_mix[i] / ac_mix - b_i / b_mix);
    return term1 - term2;
  });

  return [Math.exp(ln_phi[0]), Math.exp(ln_phi[1])];
}

// --- Replacement for calculateBubbleTemperaturePr ---
// --- Final Replacement for calculateBubbleTemperaturePr ---
// --- Final Replacement for calculateBubbleTemperaturePr ---
export function calculateBubbleTemperaturePr(
  components: CompoundData[],
  x1_feed: number,
  P_system_Pa: number,
  params: PrInteractionParams,
  _initialTempGuess_K: number = 350,
  maxIter: number = 100,
  tolerance: number = 1e-6
): BubbleDewResult | null {
    const x = [x1_feed, 1 - x1_feed];
    const T_bp1 = calculateSaturationTemperaturePurePr(components[0], P_system_Pa);
    const T_bp2 = calculateSaturationTemperaturePurePr(components[1], P_system_Pa);

    if (T_bp1 === null || T_bp2 === null) return null;

    const T_low = Math.min(T_bp1, T_bp2) - 30;
    const T_high = Math.max(T_bp1, T_bp2) + 30;

    const objective = (T: number): number => {
        let y_k = [...x];
        let K: [number, number] = [1.0, 1.0];
        const phiL = calculatePrFugacityCoefficients(components, x, T, P_system_Pa, params, 'liquid');
        if (!phiL) return NaN;

        let y_km1 = [...y_k];
        let g_km1 = [...y_k];

        for (let j = 0; j < 50; j++) {
            const phiV = calculatePrFugacityCoefficients(components, y_k, T, P_system_Pa, params, 'vapor');
            if (!phiV) return NaN;
            K = [phiL[0] / phiV[0], phiL[1] / phiV[1]];
            if (K.some(val => isNaN(val))) return NaN;
            const y_calc = [x[0] * K[0], x[1] * K[1]];
            const sum_y = y_calc[0] + y_calc[1];
            if (sum_y === 0 || !isFinite(sum_y)) return NaN;
            const g_k: [number, number] = [y_calc[0] / sum_y, y_calc[1] / sum_y];

            let y_kp1: [number, number];

            if (j < 2) { // Start with damped simple iteration
                y_kp1 = [0.5 * y_k[0] + 0.5 * g_k[0], 0.5 * y_k[1] + 0.5 * g_k[1]];
            } else { // Switch to Wegstein's method
                const s_num = g_k[0] - g_km1[0];
                const s_den = y_k[0] - y_km1[0];
                if (Math.abs(s_den) < 1e-10) {
                    y_kp1 = g_k;
                } else {
                    const s = s_num / s_den;
                    if (Math.abs(s - 1.0) > 1e-8) {
                        let q = s / (s - 1.0);
                        q = Math.max(-5.0, q);
                        const y0_next = y_k[0] + (1 - q) * (g_k[0] - y_k[0]);
                        y_kp1 = [Math.max(0, Math.min(1, y0_next)), 1.0 - Math.max(0, Math.min(1, y0_next))];
                    } else { y_kp1 = g_k; }
                }
            }

            if (Math.abs(y_kp1[0] - y_k[0]) < 1e-9) {
                y_k = y_kp1;
                break;
            }
            y_km1 = y_k;
            g_km1 = g_k;
            y_k = y_kp1;
        }

        const phiV_final = calculatePrFugacityCoefficients(components, y_k, T, P_system_Pa, params, 'vapor');
        if (!phiV_final) return NaN;
        K = [phiL[0] / phiV_final[0], phiL[1] / phiV_final[1]];
        return x[0] * K[0] + x[1] * K[1] - 1.0;
    };

    const T_bubble = _brentRoot(objective, T_low, T_high, tolerance, maxIter);
    if (T_bubble === null) return null;

    let y_final = [...x];
    const phiL_final = calculatePrFugacityCoefficients(components, x, T_bubble, P_system_Pa, params, 'liquid');
    if (!phiL_final) return null;
    for (let j = 0; j < 50; j++) {
        const y_old = [...y_final];
        const phiV_final = calculatePrFugacityCoefficients(components, y_final, T_bubble, P_system_Pa, params, 'vapor');
        if (!phiV_final) break;
        const K = [phiL_final[0] / phiV_final[0], phiL_final[1] / phiV_final[1]];
        if (K.some(k => isNaN(k))) break;
        const y_calc = [x[0] * K[0], x[1] * K[1]];
        const sum_y = y_calc[0] + y_calc[1];
        if (sum_y === 0 || !isFinite(sum_y)) break;
        y_final = [y_calc[0] / sum_y, y_calc[1] / sum_y];
        if (Math.abs(y_final[0] - y_old[0]) < 1e-8) break;
    }
    return { comp1_feed: x1_feed, comp1_equilibrium: y_final[0], T_K: T_bubble, P_Pa: P_system_Pa };
}

export function calculateBubblePressurePr(
  components: CompoundData[],
  x1_feed: number,
  T_K: number,
  params: PrInteractionParams,
  P_guess_Pa: number = 101_325,
  maxIter = 40,
  tol = 1e-5
): BubbleDewResult | null {
  let P = P_guess_Pa;
  const x = [x1_feed, 1 - x1_feed];
  for (let i = 0; i < maxIter; i++) {
    const phiL = calculatePrFugacityCoefficients(components, x, T_K, P, params, 'liquid');
    const phiV = calculatePrFugacityCoefficients(components, x, T_K, P, params, 'vapor');
    if (!phiL || !phiV) return null;
    const K = [phiL[0] / phiV[0], phiL[1] / phiV[1]];
    const f = x[0] * K[0] + x[1] * K[1] - 1;
    if (Math.abs(f) < tol) {
      const y1 = x[0] * K[0];
      return { comp1_feed: x1_feed, comp1_equilibrium: y1, T_K, P_Pa: P };
    }
    P *= 1 + f; // simple update
    if (P <= 0) return null;
  }
  return null;
}

// --- Final Replacement for calculateDewTemperaturePr ---
export function calculateDewTemperaturePr(
  components: CompoundData[],
  y1_feed: number,
  P_system_Pa: number,
  prInteractionParams: PrInteractionParams,
  _initialTempGuess_K: number = 350,
  maxIter: number = 100,
  tolerance: number = 1e-6
): BubbleDewResult | null {
    const y = [y1_feed, 1 - y1_feed];
    const T_bp1 = calculateSaturationTemperaturePurePr(components[0], P_system_Pa);
    const T_bp2 = calculateSaturationTemperaturePurePr(components[1], P_system_Pa);

    if (T_bp1 === null || T_bp2 === null) return null;

    const T_low = Math.min(T_bp1, T_bp2) - 30;
    const T_high = Math.max(T_bp1, T_bp2) + 30;

    const objective = (T: number): number => {
        let x_k = [...y];
        let K: [number, number] = [1.0, 1.0];
        const phiV = calculatePrFugacityCoefficients(components, y, T, P_system_Pa, prInteractionParams, 'vapor');
        if (!phiV) return NaN;

        let x_km1 = [...x_k];
        let g_km1 = [...x_k];

        for (let j = 0; j < 50; j++) {
            const phiL = calculatePrFugacityCoefficients(components, x_k, T, P_system_Pa, prInteractionParams, 'liquid');
            if (!phiL) return NaN;
            K = [phiL[0] / phiV[0], phiL[1] / phiV[1]];
            if (K.some(val => isNaN(val) || val <= 0)) return NaN;
            const x_calc = [y[0] / K[0], y[1] / K[1]];
            const sum_x = x_calc[0] + x_calc[1];
            if (sum_x === 0 || !isFinite(sum_x)) return NaN;
            const g_k: [number, number] = [x_calc[0] / sum_x, x_calc[1] / sum_x];

            let x_kp1: [number, number];
            if (j < 2) {
                x_kp1 = [0.5 * x_k[0] + 0.5 * g_k[0], 0.5 * x_k[1] + 0.5 * g_k[1]];
            } else {
                const s_num = g_k[0] - g_km1[0];
                const s_den = x_k[0] - x_km1[0];
                if (Math.abs(s_den) < 1e-10) {
                    x_kp1 = g_k;
                } else {
                    const s = s_num / s_den;
                    if (Math.abs(s - 1.0) > 1e-8) {
                        let q = s / (s - 1.0);
                        q = Math.max(-5.0, q);
                        const x0_next = x_k[0] + (1 - q) * (g_k[0] - x_k[0]);
                        x_kp1 = [Math.max(0, Math.min(1, x0_next)), 1.0 - Math.max(0, Math.min(1, x0_next))];
                    } else { x_kp1 = g_k; }
                }
            }

            if (Math.abs(x_kp1[0] - x_k[0]) < 1e-9) {
                x_k = x_kp1;
                break;
            }
            x_km1 = x_k;
            g_km1 = g_k;
            x_k = x_kp1;
        }

        const phiL_final = calculatePrFugacityCoefficients(components, x_k, T, P_system_Pa, prInteractionParams, 'liquid');
        if (!phiL_final) return NaN;
        K = [phiL_final[0] / phiV[0], phiL_final[1] / phiV[1]];
        if (K.some(val => isNaN(val) || val <= 0)) return NaN;
        return y[0] / K[0] + y[1] / K[1] - 1.0;
    };

    const T_dew = _brentRoot(objective, T_low, T_high, tolerance, maxIter);
    if (T_dew === null) return null;

    let x_final = [...y];
    const phiV_final = calculatePrFugacityCoefficients(components, y, T_dew, P_system_Pa, prInteractionParams, 'vapor');
    if (!phiV_final) return null;
    for (let j = 0; j < 50; j++) {
        const x_old = [...x_final];
        const phiL_final = calculatePrFugacityCoefficients(components, x_final, T_dew, P_system_Pa, prInteractionParams, 'liquid');
        if (!phiL_final) break;
        const K = [phiL_final[0] / phiV_final[0], phiL_final[1] / phiV_final[1]];
        if (K.some(k => isNaN(k) || k <= 0)) break;
        const x_calc = [y[0] / K[0], y[1] / K[1]];
        const sum_x = x_calc[0] + x_calc[1];
        if (sum_x === 0 || !isFinite(sum_x)) break;
        x_final = [x_calc[0] / sum_x, x_calc[1] / sum_x];
        if (Math.abs(x_final[0] - x_old[0]) < 1e-8) break;
    }
    return { comp1_feed: y1_feed, comp1_equilibrium: x_final[0], T_K: T_dew, P_Pa: P_system_Pa };
}

export function calculateDewPressurePr(
  components: CompoundData[],
  y1_feed: number,
  T_system_K: number,
  prInteractionParams: PrInteractionParams,
  initialPressureGuess_Pa: number = 101325,
  maxIter: number = 10,
  tolerance: number = 1e-5
): BubbleDewResult | null {
  let P_current_Pa = initialPressureGuess_Pa;
  const y = [y1_feed, 1 - y1_feed];

  if (!components[0]?.prParams || !components[1]?.prParams) return null;
  if (!components[0]?.antoine || !components[1]?.antoine) return null;

  for (let iter = 0; iter < maxIter; iter++) {
    const Psat = components.map(c => calculatePsat_Pa(c.antoine, T_system_K));
    if (Psat.some(p => isNaN(p) || p <= 0)) return null;

    // Calculate fugacity coefficients
    const phi_V = calculatePrFugacityCoefficients(components, y, T_system_K, P_current_Pa, prInteractionParams, 'vapor');
    if (!phi_V) return null;

    // Initial liquid composition guess
    let x = y.map((yi, i) => yi * phi_V[i] * P_current_Pa / Psat[i]);
    const x_sum = x.reduce((s, xi) => s + xi, 0);
    x = x.map(xi => xi / x_sum);

    const phi_L = calculatePrFugacityCoefficients(components, x, T_system_K, P_current_Pa, prInteractionParams, 'liquid');
    if (!phi_L) return null;

    // Calculate new pressure
    const P_new_Pa = 1 / y.reduce((sum, yi, i) => sum + yi / (phi_L[i] * Psat[i] / phi_V[i]), 0);
    
    if (Math.abs(P_new_Pa - P_current_Pa) < tolerance * P_current_Pa) {
      return {
        comp1_feed: y1_feed,
        comp1_equilibrium: x[0],
        T_K: T_system_K,
        P_Pa: P_new_Pa,
        iterations: iter + 1,
        calculationType: 'dewP'
      };
    }
    P_current_Pa = P_new_Pa;
  }
  return null;
}

// ===============================================================
//  6. S O A P R I S O N – R E D L I C H – K W O N G   E O S
//    (Simplified binary implementation, parallels PR helpers)
// ===============================================================

export interface SrkInteractionParams { k_ij: number; }

// (interface already declared above)

export function calculateSrkFugacityCoefficients(
  components: CompoundData[],
  x_or_y: number[],
  T_K: number,
  P_Pa: number,
  interactionParams: SrkInteractionParams,
  phase: 'liquid' | 'vapor'
): [number, number] | null {
  if (SRK_DEBUG) console.log(`[SRK Fugacity] Phase: ${phase}, T: ${T_K.toFixed(2)} K, x/y: [${x_or_y[0].toFixed(4)}]`);

  if (components.length !== 2 || x_or_y.length !== 2) return null;
  if (!components[0].srkParams || !components[1].srkParams) return null;

  const x1 = x_or_y[0];
  const x2 = x_or_y[1];
  const { k_ij } = interactionParams;

  const srkPure = components.map(c => {
    const { Tc_K, Pc_Pa, omega } = c.srkParams!;
    const Tr = T_K / Tc_K;
    const m = 0.480 + 1.574 * omega - 0.176 * omega * omega;
    const alpha = Math.pow(1 + m * (1 - Math.sqrt(Tr)), 2);
    const ac_i = 0.42748 * Math.pow(R_gas_const_J_molK * Tc_K, 2) / Pc_Pa * alpha;
    const b_i = 0.08664 * R_gas_const_J_molK * Tc_K / Pc_Pa;
    return { ac_i, b_i, sqrt_ac_i: Math.sqrt(ac_i) };
  });

  const a12 = srkPure[0].sqrt_ac_i * srkPure[1].sqrt_ac_i * (1 - k_ij);
  const a_mix = x1*x1*srkPure[0].ac_i + 2*x1*x2*a12 + x2*x2*srkPure[1].ac_i;
  const b_mix = x1*srkPure[0].b_i + x2*srkPure[1].b_i;

  const A_mix = a_mix * P_Pa / (R_gas_const_J_molK * T_K)**2;
  const B_mix = b_mix * P_Pa / (R_gas_const_J_molK * T_K);
  if (SRK_DEBUG) console.log(`    A_mix: ${A_mix.toExponential(4)}, B_mix: ${B_mix.toExponential(4)}`);

  const roots = solveCubicEOS(1, -1, A_mix - B_mix - B_mix*B_mix, -A_mix*B_mix);
  if (!roots || roots.length === 0) {
    console.error(`    [SRK Fugacity] FAILED: No valid roots found for Z.`);
    return null;
  }
  if (SRK_DEBUG) console.log(`    Z roots found: [${roots.map(r => r.toFixed(5)).join(', ')}]`);

  const Z = phase === 'liquid' ? Math.min(...roots) : Math.max(...roots);
  if (SRK_DEBUG) console.log(`    Selected Z for ${phase}: ${Z.toFixed(5)}`);
  
  if (Z <= B_mix) {
    console.warn(`    [SRK Fugacity] FAILED: Z <= B_mix (${Z.toFixed(5)} <= ${B_mix.toFixed(5)})`);
    return null;
  }

  const term_ln = (A_mix / B_mix) * Math.log(1 + B_mix / Z);

  const d_a_mix_by_n = [
      2 * (x1 * srkPure[0].ac_i + x2 * a12),
      2 * (x2 * srkPure[1].ac_i + x1 * a12)
  ];

  const phi = [0, 1].map(i => {
    const b_i = srkPure[i].b_i;
    const b_i_over_b_mix = b_i / b_mix;
    const term1 = b_i_over_b_mix * (Z - 1) - Math.log(Z - B_mix);
    const term2 = term_ln * (d_a_mix_by_n[i] / a_mix - b_i_over_b_mix);
    return Math.exp(term1 - term2);
  });
  
  if (SRK_DEBUG) console.log(`    ==> Fugacity Coeffs (phi): [${phi[0].toFixed(5)}, ${phi[1].toFixed(5)}]`);
  return [phi[0], phi[1]];
}

export function calculateBubbleTemperatureSrk(
  components: CompoundData[],
  x1_feed: number,
  P_system_Pa: number,
  srkInteractionParams: SrkInteractionParams,
  _initialTempGuess_K: number = 350,
  maxIter: number = 100,
  tolerance: number = 1e-6
): BubbleDewResult | null {
    if (SRK_DEBUG) console.log(`%c[SRK Bubble T] START for x1=${x1_feed}`, 'color: blue; font-weight: bold;');
    
    const x = [x1_feed, 1 - x1_feed];
    const T_bp1 = calculateSaturationTemperaturePureSrk(components[0], P_system_Pa);
    const T_bp2 = calculateSaturationTemperaturePureSrk(components[1], P_system_Pa);

    if (T_bp1 === null || T_bp2 === null) {
        console.error("[SRK Bubble T] FAILED: Could not calculate pure component boiling points for bracketing.");
        return null;
    }

    const T_low = Math.min(T_bp1, T_bp2) - 30;
    const T_high = Math.max(T_bp1, T_bp2) + 30;
    console.log(`[SRK Bubble T] Temp Bracket: [${T_low.toFixed(2)}, ${T_high.toFixed(2)}] K`);

    const objective = (T: number): number => {
        if (SRK_DEBUG) console.log(`%c[SRK Bubble T] > Root finder testing T = ${T.toFixed(4)} K`, 'color: gray');
        const phiL = calculateSrkFugacityCoefficients(components, x, T, P_system_Pa, srkInteractionParams, 'liquid');
        if (!phiL) return NaN;

        let y = [...x]; // Start guess with y=x
        if (SRK_DEBUG) console.log(`    Inner loop START, initial y: [${y[0].toFixed(4)}]`);
        for (let i = 0; i < 50; i++) {
            const y_old = [...y];
            const phiV = calculateSrkFugacityCoefficients(components, y, T, P_system_Pa, srkInteractionParams, 'vapor');
            if (!phiV) {
                 console.warn(`    Inner loop break: phiV calculation failed.`);
                 return NaN;
            }

            const K = [phiL[0] / phiV[0], phiL[1] / phiV[1]];
            if (SRK_DEBUG) console.log(`    Iter ${i}, K-values: [${K[0].toFixed(4)}, ${K[1].toFixed(4)}]`);

            const y_new_unnorm = [x[0] * K[0], x[1] * K[1]];
            const sumY = y_new_unnorm[0] + y_new_unnorm[1];

            if (sumY === 0 || !isFinite(sumY)) {
                console.warn(`    Inner loop break: Invalid sumY.`);
                return NaN;
            }
            y = [y_new_unnorm[0] / sumY, y_new_unnorm[1] / sumY];
            if (SRK_DEBUG) console.log(`    Iter ${i}, new y: [${y[0].toFixed(4)}]`);

            if (Math.abs(y[0] - y_old[0]) < 1e-9) {
                if (SRK_DEBUG) console.log(`    Inner loop converged after ${i+1} iterations.`);
                break;
            }
        }

        const K_final = [phiL[0] / calculateSrkFugacityCoefficients(components, y, T, P_system_Pa, srkInteractionParams, 'vapor')![0], phiL[1] / calculateSrkFugacityCoefficients(components, y, T, P_system_Pa, srkInteractionParams, 'vapor')![1]];
        const objective_val = x[0] * K_final[0] + x[1] * K_final[1] - 1.0;
        if (SRK_DEBUG) console.log(`    > Objective Function Result for T=${T.toFixed(4)} K is ${objective_val.toExponential(4)}`);
        return objective_val;
    };

    const T_bubble = _brentRoot(objective, T_low, T_high, tolerance, maxIter);
    if (T_bubble === null) {
        console.error("[SRK Bubble T] FAILED: Root finder could not converge on a bubble temperature.");
        return null;
    }
    if (SRK_DEBUG) console.log(`%c[SRK Bubble T] Root finder SUCCESS. T_bubble = ${T_bubble.toFixed(4)} K`, 'color: green');

    // Final calculation of y at the converged temperature
    let y_final = [...x];
    const phiL_final = calculateSrkFugacityCoefficients(components, x, T_bubble, P_system_Pa, srkInteractionParams, 'liquid');
    if (!phiL_final) return null;
    const phiV_final = calculateSrkFugacityCoefficients(components, y_final, T_bubble, P_system_Pa, srkInteractionParams, 'vapor');
    if (!phiV_final) return null;
    const K = [phiL_final[0] / phiV_final[0], phiL_final[1] / phiV_final[1]];
    const y_calc = [x[0] * K[0], x[1] * K[1]];
    const sum_y = y_calc[0] + y_calc[1];
    y_final = [y_calc[0] / sum_y, y_calc[1] / sum_y];
    
    if (SRK_DEBUG) console.log(`%c[SRK Bubble T] FINAL RESULT: x1=${x1_feed}, y1=${y_final[0].toFixed(5)}, T=${T_bubble.toFixed(2)} K`, 'color: blue; font-weight: bold;');
    return { comp1_feed: x1_feed, comp1_equilibrium: y_final[0], T_K: T_bubble, P_Pa: P_system_Pa };
}

export function calculateBubblePressureSrk(
  components: CompoundData[],
  x1_feed: number,
  T_system_K: number,
  srkInteractionParams: SrkInteractionParams,
  initialPressureGuess_Pa: number = 101325,
  maxIter: number = 40,
  tolerance: number = 1e-5
): BubbleDewResult | null {
  let P = initialPressureGuess_Pa;
  const x = [x1_feed, 1 - x1_feed];
  let y = [...x]; 
  let K = [1.0, 1.0];

  for (let i = 0; i < maxIter; i++) {
    const phiL = calculateSrkFugacityCoefficients(components, x, T_system_K, P, srkInteractionParams, 'liquid');
    if (!phiL) { P *= 0.8; continue; } // Nudge P on failure

    // Inner loop to converge y for the current P
    for (let j = 0; j < 10; j++) {
        const phiV = calculateSrkFugacityCoefficients(components, y, T_system_K, P, srkInteractionParams, 'vapor');
        if (!phiV) { y = [...x]; continue; }
        
        K = [phiL[0] / phiV[0], phiL[1] / phiV[1]];
        
        const y_new = [x[0] * K[0], x[1] * K[1]];
        const sum_y_new = y_new[0] + y_new[1];
        if (sum_y_new === 0) return null;
        y_new[0] /= sum_y_new;
        y_new[1] /= sum_y_new;

        if (Math.abs(y_new[0] - y[0]) < 1e-7) {
            y = y_new; break;
        }
        y = y_new;
    }
    
    const sum_Kx = x[0] * K[0] + x[1] * K[1];
    const P_new = P * sum_Kx;

    if (Math.abs(P_new - P) < tolerance * P) {
      return { comp1_feed: x1_feed, comp1_equilibrium: y[0], T_K: T_system_K, P_Pa: P_new, iterations: i + 1, calculationType: 'bubbleP' };
    }

    P = P_new;
    if (P <= 0 || !isFinite(P)) return null;
  }
  return null;
}

/**
 * Calculates the dew point temperature for a binary mixture using the SRK EOS.
 * This function has been rewritten for improved stability and accuracy.
 */
export function calculateDewTemperatureSrk(
  components: CompoundData[],
  y1_feed: number,
  P_system_Pa: number,
  srkInteractionParams: SrkInteractionParams,
  _initialTempGuess_K: number = 350,
  maxIter: number = 100,
  tolerance: number = 1e-6
): BubbleDewResult | null {
    const y = [y1_feed, 1 - y1_feed];
    const T_bp1 = calculateSaturationTemperaturePureSrk(components[0], P_system_Pa);
    const T_bp2 = calculateSaturationTemperaturePureSrk(components[1], P_system_Pa);

    if (T_bp1 === null || T_bp2 === null) {
      return { error: "Failed to find pure component boiling points for bracketing.", comp1_feed: y1_feed, comp1_equilibrium: NaN };
    }

    const T_low = Math.min(T_bp1, T_bp2) - 20;
    const T_high = Math.max(T_bp1, T_bp2) + 20;

    const objective = (T: number): number => {
        const phiV = calculateSrkFugacityCoefficients(components, y, T, P_system_Pa, srkInteractionParams, 'vapor');
        if (!phiV) return NaN;

        // Stable inner loop to find liquid composition 'x'
        let x = [...y]; // Initial guess
        for (let i = 0; i < 20; i++) {
            const x_old = [...x];
            const phiL = calculateSrkFugacityCoefficients(components, x, T, P_system_Pa, srkInteractionParams, 'liquid');
            if (!phiL) return NaN;
            const K = [phiL[0] / phiV[0], phiL[1] / phiV[1]];
            if (K.some(k => k <= 0 || isNaN(k))) return NaN;
            
            const x_new_unnorm = [y[0] / K[0], y[1] / K[1]];
            const sumX = x_new_unnorm[0] + x_new_unnorm[1];
            if (sumX === 0 || !isFinite(sumX)) return NaN;

            // Damped update
            x[0] = 0.7 * x[0] + 0.3 * (x_new_unnorm[0] / sumX);
            x[1] = 1.0 - x[0];

            if (Math.abs(x[0] - x_old[0]) < 1e-9) break;
        }

        const phiL_converged = calculateSrkFugacityCoefficients(components, x, T, P_system_Pa, srkInteractionParams, 'liquid');
        if (!phiL_converged) return NaN;
        const K_converged = [phiL_converged[0] / phiV[0], phiL_converged[1] / phiV[1]];
        if (K_converged.some(k => k <= 0 || isNaN(k))) return NaN;
        return y[0] / K_converged[0] + y[1] / K_converged[1] - 1.0;
    };

    const T_dew = _brentRoot(objective, T_low, T_high, tolerance, maxIter);
    if (T_dew === null) return null;

    // Final pass to calculate equilibrium composition accurately
    let x_final = [...y];
    const phiV_final = calculateSrkFugacityCoefficients(components, y, T_dew, P_system_Pa, srkInteractionParams, 'vapor');
    if (!phiV_final) return null;
    for (let i = 0; i < 20; i++) {
        const phiL_final = calculateSrkFugacityCoefficients(components, x_final, T_dew, P_system_Pa, srkInteractionParams, 'liquid');
        if (!phiL_final) break;
        const K = [phiL_final[0] / phiV_final[0], phiL_final[1] / phiV_final[1]];
        if (K.some(k => k <= 0 || isNaN(k))) break;
        const x_new_unnorm = [y[0] / K[0], y[1] / K[1]];
        const sumX = x_new_unnorm[0] + x_new_unnorm[1];
        if (sumX === 0 || !isFinite(sumX)) break;
        const x_new = [x_new_unnorm[0] / sumX, x_new_unnorm[1] / sumX];
        if (Math.abs(x_final[0] - x_new[0]) < 1e-9) {
            x_final = x_new;
            break;
        }
        x_final = x_new;
    }
    return { comp1_feed: y1_feed, comp1_equilibrium: x_final[0], T_K: T_dew, P_Pa: P_system_Pa };
}

export function calculateDewPressureSrk(
  components: CompoundData[],
  y1_feed: number,
  T_system_K: number,
  srkInteractionParams: SrkInteractionParams,
  initialPressureGuess_Pa: number = 101325,
  maxIter: number = 40,
  tolerance: number = 1e-5
): BubbleDewResult | null {
  let P = initialPressureGuess_Pa;
  const y = [y1_feed, 1 - y1_feed];
  let x = [...y];
  let K = [1.0, 1.0];

  for (let i = 0; i < maxIter; i++) {
    const phiV = calculateSrkFugacityCoefficients(components, y, T_system_K, P, srkInteractionParams, 'vapor');
    if (!phiV) { P *= 1.2; continue; }

    for (let j = 0; j < 10; j++) {
        const phiL = calculateSrkFugacityCoefficients(components, x, T_system_K, P, srkInteractionParams, 'liquid');
        if (!phiL) { x = [...y]; continue; }

        K = [phiL[0] / phiV[0], phiL[1] / phiV[1]];
        
        const x_new = [y[0] / K[0], y[1] / K[1]];
        const sum_x_new = x_new[0] + x_new[1];
        if (sum_x_new === 0) return null;
        x_new[0] /= sum_x_new;
        x_new[1] /= sum_x_new;
        
        if (Math.abs(x_new[0] - x[0]) < 1e-7) {
            x = x_new; break;
        }
        x = x_new;
    }

    const sum_y_over_K = y[0] / K[0] + y[1] / K[1];
    const P_new = P / sum_y_over_K;

    if (Math.abs(P_new - P) < tolerance * P) {
      return { comp1_feed: y1_feed, comp1_equilibrium: x[0], T_K: T_system_K, P_Pa: P_new, iterations: i + 1, calculationType: 'dewP' };
    }

    P = P_new;
    if (P <= 0 || !isFinite(P)) return null;
  }
  return null;
}

// ===============================================================
//  7. U N I Q U A C   (binary activity model)
// ===============================================================

/**
 * UNIQUAC interaction params – HYSYS convention:
 *   Full Aspen: τ_ij = exp(aij + bij/T + cij·lnT + dij·T)
 *   Stored scaled: exp(−(Aij + Bij·T + Cij·T·lnT + Dij·T²) / (R_cal·T))
 *   Aij = −bij·R_cal (K), Bij = −aij·R_cal (dimensionless),
 *   Cij = −cij·R_cal, Dij = −dij·R_cal.
 */
export interface UniquacInteractionParams {
  Aij: number; // K (= −bij·R_cal)
  Aji: number; // K (= −bji·R_cal)
  Bij: number; // dimensionless (= −aij·R_cal)
  Bji: number; // dimensionless (= −aji·R_cal)
  /** Aspen cij·(−R_cal): ln T coefficient. Default 0. */
  Cij?: number;
  /** Aspen cji·(−R_cal): ln T coefficient. Default 0. */
  Cji?: number;
  /** Aspen dij·(−R_cal): T coefficient. Default 0. */
  Dij?: number;
  /** Aspen dji·(−R_cal): T coefficient. Default 0. */
  Dji?: number;
  /** Hayden-O'Connell cross-association parameter (η_ij). Default 0 (ideal vapor). */
  eta_cross?: number;
}

// (interface already declared above)

export function calculateUniquacGamma(
  comps: CompoundData[],
  x: number[],
  T_K: number,
  p: UniquacInteractionParams
): [number, number] | null {
  if (comps.length !== 2 || x.length !== 2) return null;
  if (!comps[0].uniquacParams || !comps[1].uniquacParams) return null;

  const r = [comps[0].uniquacParams.r, comps[1].uniquacParams.r];
  const q = [comps[0].uniquacParams.q, comps[1].uniquacParams.q];
  const [x1, x2] = x;

  const z = 10;
  const L = (i: number) => (z / 2) * (r[i] - q[i]) - (r[i] - 1);

  const sum_xr = x1 * r[0] + x2 * r[1];
  const sum_xq = x1 * q[0] + x2 * q[1];
  const Phi = [(x1 * r[0]) / sum_xr, (x2 * r[1]) / sum_xr];
  const Theta = [(x1 * q[0]) / sum_xq, (x2 * q[1]) / sum_xq];

  const lnGammaC = [0, 1].map(i =>
    Math.log(Phi[i] / x[i]) + (z / 2) * q[i] * Math.log(Theta[i] / Phi[i]) + L(i) - (Phi[i] / x[i]) * (x1 * L(0) + x2 * L(1))
  );

  // Aspen UNIQUAC: Aij = -bij*R_cal, Bij = -aij*R_cal (see fetchUniquacInteractionParams)
  // Full Aspen: τ_ij = exp(aij + bij/T + cij·lnT + dij·T)
  // Stored scaled: exp(−(Aij + Bij·T + Cij·T·lnT + Dij·T²) / (R_cal·T))  (Cij = -cij·R_cal, Dij = -dij·R_cal)
  const R_cal = 1.9872; // cal·mol⁻¹·K⁻¹
  // Use standard convention: tau12 uses Aij/Bij and tau21 uses Aji/Bji
  const tau12 = Math.exp(
    -(p.Aij + p.Bij * T_K + (p.Cij ?? 0) * T_K * Math.log(T_K) + (p.Dij ?? 0) * T_K * T_K) / (R_cal * T_K)
  );
  const tau21 = Math.exp(
    -(p.Aji + p.Bji * T_K + (p.Cji ?? 0) * T_K * Math.log(T_K) + (p.Dji ?? 0) * T_K * T_K) / (R_cal * T_K)
  );

  const Theta_res = Theta; // same symbols
  const s1 = Theta_res[0] + Theta_res[1] * tau21;
  const s2 = Theta_res[1] + Theta_res[0] * tau12;

  const lnGammaR0 = q[0] * (1 - Math.log(s1) - (Theta_res[0] / s1 + (Theta_res[1] * tau12) / s2));
  const lnGammaR1 = q[1] * (1 - Math.log(s2) - (Theta_res[1] / s2 + (Theta_res[0] * tau21) / s1));

  return [Math.exp(lnGammaC[0] + lnGammaR0), Math.exp(lnGammaC[1] + lnGammaR1)];
}

export function calculateBubbleTemperatureUniquac(
  comps: CompoundData[],
  x1_feed: number,
  P_system_Pa: number,
  params: UniquacInteractionParams,
  T_guess_K: number = 350,
  maxIter = 60,
  tol = 1e-5
): BubbleDewResult | null {
  const x = [x1_feed, 1 - x1_feed];
  let T = T_guess_K;
  let y: number[] = [...x];
  const eta_cross = params.eta_cross ?? 0;
  const useHoc = eta_cross > 0 ||
    (comps[0].hocProps?.eta_self ?? 0) > 0 ||
    (comps[1].hocProps?.eta_self ?? 0) > 0;

  for (let i = 0; i < maxIter; i++) {
    const Psat = [calculatePsat_Pa(comps[0].antoine, T), calculatePsat_Pa(comps[1].antoine, T)];
    if (Psat.some(p => isNaN(p))) return null;
    const gamma = calculateUniquacGamma(comps, x, T, params);
    if (!gamma) return null;

    const phi = useHoc
      ? hocFugacityCorrection(comps, y, T, P_system_Pa, Psat, eta_cross)
      : ([1, 1] as [number, number]);

    const KPsat = [x[0] * gamma[0] * Psat[0] * phi[0], x[1] * gamma[1] * Psat[1] * phi[1]];
    const Pcalc = KPsat[0] + KPsat[1];
    const err = Pcalc - P_system_Pa;
    if (Math.abs(err) < tol * P_system_Pa) {
      y = [KPsat[0] / Pcalc, KPsat[1] / Pcalc];
      return { comp1_feed: x1_feed, comp1_equilibrium: y[0], T_K: T, P_Pa: P_system_Pa };
    }
    const Ptmp = Pcalc || P_system_Pa;
    y = [KPsat[0] / Ptmp, KPsat[1] / Ptmp];
    T -= err / P_system_Pa * (T * 0.05);
    if (T < 100 || T > 1000) return null;
  }
  return null;
}

export function calculateBubblePressureUniquac(
  comps: CompoundData[],
  x1_feed: number,
  T_K: number,
  params: UniquacInteractionParams
): BubbleDewResult | null {
  const x = [x1_feed, 1 - x1_feed];
  const Psat = [calculatePsat_Pa(comps[0].antoine, T_K), calculatePsat_Pa(comps[1].antoine, T_K)];
  if (Psat.some(p => isNaN(p))) return null;
  const gamma = calculateUniquacGamma(comps, x, T_K, params);
  if (!gamma) return null;

  const eta_cross = params.eta_cross ?? 0;
  const useHoc = eta_cross > 0 ||
    (comps[0].hocProps?.eta_self ?? 0) > 0 ||
    (comps[1].hocProps?.eta_self ?? 0) > 0;

  const PcalcIdeal = x[0] * gamma[0] * Psat[0] + x[1] * gamma[1] * Psat[1];
  const yIdeal = [x[0] * gamma[0] * Psat[0] / PcalcIdeal, x[1] * gamma[1] * Psat[1] / PcalcIdeal];
  const phi = useHoc
    ? hocFugacityCorrection(comps, yIdeal, T_K, PcalcIdeal, Psat, eta_cross)
    : ([1, 1] as [number, number]);

  const KPsat = [x[0] * gamma[0] * Psat[0] * phi[0], x[1] * gamma[1] * Psat[1] * phi[1]];
  const Pcalc = KPsat[0] + KPsat[1];
  const y1 = KPsat[0] / Pcalc;
  return { comp1_feed: x1_feed, comp1_equilibrium: y1, T_K, P_Pa: Pcalc };
}

// ==============================================================
//  7.  D A T A B A S E   F E T C H   H E L P E R S (Supabase)
// ==============================================================

// ─────────────────────────────────────────────────────────────────────────────
//  H A Y D E N – O ' C O N N E L L   (HOC)  second virial coefficient
// ─────────────────────────────────────────────────────────────────────────────
/**
 * Compute the Hayden-O'Connell (1975) second virial coefficient B_ij [m³/mol].
 *
 * @param T_K     Temperature [K]
 * @param Tc_ij   Cross critical temperature [K]  (sqrt(Tc_i·Tc_j) for cross)
 * @param Vc_ij   Cross critical volume [cc/mol]  ((((Vc_i^1/3+Vc_j^1/3)/2)^3) for cross)
 * @param omega_ij  Cross acentric factor         ((ω_i+ω_j)/2 for cross)
 * @param mu_ij_D  Cross dipole moment [Debye]    (sqrt(μ_i·μ_j) for cross)
 * @param rd_ij_A  Cross mean radius of gyration [Å]  ((rd_i+rd_j)/2 for cross)
 * @param eta_ij  HOC association parameter       (0 = non-associating)
 * @returns B_ij [m³/mol]
 */
export function hocSecondVirial(
  T_K: number,
  Tc_ij: number,
  Vc_ij_ccmol: number,
  omega_ij: number,
  _mu_ij_D: number,
  _rd_ij_A: number,
  eta_ij: number
): number {
  const Zc_ij = 0.291 - 0.080 * omega_ij;
  const b0 = Vc_ij_ccmol / Zc_ij; // cc/mol  (= R·Tc/Pc, the Tsonopoulos reference volume)
  const Tr = T_K / Tc_ij;

  // Non-polar Pitzer-Curl contribution
  const B_NP = b0 * ((0.083 - 0.422 / Math.pow(Tr, 1.6))
             + omega_ij * (0.139 - 0.172 / Math.pow(Tr, 4.2)));

  // Chemical / association contribution
  // The original HOC constant (1500) is catastrophically large — empirically C≈50 gives
  // physically reasonable B values for H-bonding species (MeOH, Acetone) at 1 atm.
  let B_chem = 0;
  if (eta_ij > 0.01) {
    const K_dim = 4.5 * b0 * Math.exp(eta_ij * 50.0 / T_K);
    B_chem = -0.5 * K_dim;
  }

  // Convert cc/mol → m³/mol
  return (B_NP + B_chem) * 1e-6;
}

/**
 * Build HOC combining-rule parameters for a cross (i,j) pair from pure-component data.
 * All required inputs come from CompoundData.hocProps (mu_D, rd_A, eta info) and
 * criticalProperties (Tc_K, omega) plus Vc from wilsonParams block.
 *
 * Returns critical combining-rule values in cc/mol and K, ready for hocSecondVirial().
 */
function _hocCrossParams(
  comp_i: CompoundData,
  comp_j: CompoundData,
  Vc_i_ccmol: number,
  Vc_j_ccmol: number,
  eta_cross: number
) {
  const Tc_i = comp_i.prParams?.Tc_K ?? 0;
  const Tc_j = comp_j.prParams?.Tc_K ?? 0;
  const omega_i = comp_i.prParams?.omega ?? 0;
  const omega_j = comp_j.prParams?.omega ?? 0;
  const mu_i = comp_i.hocProps?.mu_D ?? 0;
  const mu_j = comp_j.hocProps?.mu_D ?? 0;
  const rd_i = comp_i.hocProps?.rd_A ?? 0;
  const rd_j = comp_j.hocProps?.rd_A ?? 0;

  const Tc_ij = Math.sqrt(Tc_i * Tc_j);
  const Vc1_3_i = Math.cbrt(Vc_i_ccmol);
  const Vc1_3_j = Math.cbrt(Vc_j_ccmol);
  const Vc_ij = Math.pow((Vc1_3_i + Vc1_3_j) / 2, 3);
  const omega_ij = (omega_i + omega_j) / 2;
  const mu_ij_D = Math.sqrt(mu_i * mu_j);
  const rd_ij_A = (rd_i + rd_j) / 2;
  return { Tc_ij, Vc_ij, omega_ij, mu_ij_D, rd_ij_A, eta_ij: eta_cross };
}

/**
 * Get critical volume [cc/mol] for a component from its CompoundData, falling back
 * to a Rackett estimate if the stored wilsonParams Vc isn't directly available.
 * The Vc stored in apv140_pure_props_wide is in m³/kmol → × 1000 gives cc/mol.
 */
function _vcCcMol(comp: CompoundData): number {
  // If we have prParams (always fetched), Vc must come from a separate store.
  // The fetch function stores VC_m3kmol in wilsonParams indirectly via criticalVolume_m3kmol
  // but CompoundData doesn't carry it. We reconstruct from Rackett:
  //   Zc = 0.291 − 0.08·ω; Vc [m³/mol] = Zc·R·Tc/Pc → cc/mol = ×1e6
  const p = comp.prParams;
  if (!p || !p.Tc_K || !p.Pc_Pa) return 0;
  const omega = p.omega ?? 0;
  const Zc = 0.291 - 0.080 * omega;
  return Zc * R_gas_const_J_molK * p.Tc_K / p.Pc_Pa * 1e6; // m³/mol → cc/mol
}

/**
 * Compute HOC vapor-phase fugacity correction factor for component i in a binary mixture:
 *   ln(φ_i^sat / φ_i^V) = [B_ii·Psat_i − (2·Σ_j y_j·B_ij − B_mix)·P] / (R·T)
 *   → returns exp(that), the K-value correction φ_i^sat / φ_i^V
 *
 * All B values in m³/mol, pressures in Pa, T in K.
 */
export function hocFugacityCorrection(
  comps: CompoundData[],
  y: number[],         // vapor mole fractions [y1, y2]
  T_K: number,
  P_Pa: number,
  Psat: number[],      // [Psat1, Psat2] in Pa
  eta_cross: number    // HOC cross-association parameter
): [number, number] {
  const n = 2;
  if (comps.length !== n || y.length !== n || Psat.length !== n) return [1, 1];

  // Build Vc for each component (cc/mol)
  const Vc = comps.map(_vcCcMol);
  if (Vc.some(v => v <= 0)) return [1, 1];

  // Pure self-virial B_ii [m³/mol]
  const B_self = comps.map((c, idx) => {
    const p = c.prParams;
    if (!p) return 0;
    const h = c.hocProps;
    return hocSecondVirial(T_K, p.Tc_K, Vc[idx], p.omega ?? 0,
      h?.mu_D ?? 0, h?.rd_A ?? 0, h?.eta_self ?? 0);
  });

  // Cross virial B_12 [m³/mol]
  const { Tc_ij, Vc_ij, omega_ij, mu_ij_D, rd_ij_A, eta_ij } =
    _hocCrossParams(comps[0], comps[1], Vc[0], Vc[1], eta_cross);
  const B12 = hocSecondVirial(T_K, Tc_ij, Vc_ij, omega_ij, mu_ij_D, rd_ij_A, eta_ij);

  const [y1, y2] = y;
  const B_mix = y1 * y1 * B_self[0] + 2 * y1 * y2 * B12 + y2 * y2 * B_self[1];
  const RT = R_gas_const_J_molK * T_K;

  // ln(φ_i^sat) = B_ii · Psat_i / (R·T)   (pure vapor saturation)
  // ln(φ_i^V)   = (2·(y1·B_i1 + y2·B_i2) − B_mix) · P / (R·T)
  const lnPhi_sat = B_self.map((Bii, idx) => Bii * Psat[idx] / RT);
  const lnPhi_V = [
    (2 * (y1 * B_self[0] + y2 * B12)           - B_mix) * P_Pa / RT,
    (2 * (y1 * B12           + y2 * B_self[1]) - B_mix) * P_Pa / RT,
  ];

  return [
    Math.exp(lnPhi_sat[0] - lnPhi_V[0]),
    Math.exp(lnPhi_sat[1] - lnPhi_V[1]),
  ];
}

/**
 * HOC fugacity correction for ternary mixtures.  Returns correction factor for each comp.
 * Component indices: 0, 1, 2.  eta_cross is a 3×3 symmetric matrix (eta_cross[i][j]).
 */
export function hocFugacityCorrectionTernary(
  comps: CompoundData[],
  y: number[],
  T_K: number,
  P_Pa: number,
  Psat: number[],
  eta_matrix: number[][]   // eta_matrix[i][j], symmetric, diagonal = eta_self[i]
): number[] {
  const n = 3;
  if (comps.length !== n || y.length !== n || Psat.length !== n) return [1, 1, 1];

  const Vc = comps.map(_vcCcMol);
  if (Vc.some(v => v <= 0)) return [1, 1, 1];

  // Compute all unique B_ij
  const B: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = i; j < n; j++) {
      let Bval: number;
      if (i === j) {
        const c = comps[i];
        const p = c.prParams;
        const h = c.hocProps;
        Bval = p ? hocSecondVirial(T_K, p.Tc_K, Vc[i], p.omega ?? 0,
          h?.mu_D ?? 0, h?.rd_A ?? 0, eta_matrix[i][i]) : 0;
      } else {
        const { Tc_ij, Vc_ij, omega_ij, mu_ij_D, rd_ij_A } =
          _hocCrossParams(comps[i], comps[j], Vc[i], Vc[j], eta_matrix[i][j]);
        Bval = hocSecondVirial(T_K, Tc_ij, Vc_ij, omega_ij, mu_ij_D, rd_ij_A, eta_matrix[i][j]);
      }
      B[i][j] = Bval;
      B[j][i] = Bval;
    }
  }

  // B_mix = Σ_i Σ_j y_i y_j B_ij
  let B_mix = 0;
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++)
      B_mix += y[i] * y[j] * B[i][j];

  const RT = R_gas_const_J_molK * T_K;

  return comps.map((_, i) => {
    const lnPhi_sat_i = B[i][i] * Psat[i] / RT;
    // Σ_j y_j B_ij = sum contribution for component i
    const sum_yB_i = y.reduce((s, yj, j) => s + yj * B[i][j], 0);
    const lnPhi_V_i = (2 * sum_yB_i - B_mix) * P_Pa / RT;
    return Math.exp(lnPhi_sat_i - lnPhi_V_i);
  });
}


// --- UNIFAC ---
export async function fetchUnifacInteractionParams(
  supabase: SupabaseClient,
  subgroupIds: number[]
): Promise<UnifacParameters> {
  const { data, error } = await supabase
    .from('UNIFAC - Rk and Qk')
    .select('"Subgroup #", "Main Group #", "Rk", "Qk"')
    .in('"Subgroup #"', subgroupIds);
  if (error) throw error;
  if (!data || data.length === 0) throw new Error('UNIFAC Rk/Qk data missing');

  const Rk: { [id: number]: number } = {};
  const Qk: { [id: number]: number } = {};
  const mainMap: { [id: number]: number } = {};
  const mainIds = new Set<number>();
  data.forEach((row: any) => {
    const sg = row["Subgroup #"];
    Rk[sg] = row["Rk"];
    Qk[sg] = row["Qk"];
    mainMap[sg] = row["Main Group #"];
    mainIds.add(row["Main Group #"]);
  });
  const { data: intData, error: intErr } = await supabase
    .from('UNIFAC - a(ij)')
    .select('i,j,"a(ij)"')
    .in('i', Array.from(mainIds))
    .in('j', Array.from(mainIds));
  if (intErr) throw intErr;
  const a_mk = new Map<string, number>();
  intData?.forEach((row: any) => {
    a_mk.set(`${row.i}-${row.j}`, row["a(ij)"] ?? 0);
  });
  return { Rk, Qk, mainGroupMap: mainMap, a_mk };
}

// ─── Shared R_cal constant ───────────────────────────────────────────────────
const R_cal = 1.9872; // cal·mol⁻¹·K⁻¹  (used to scale Aspen K-units → HYSYS cal/mol-units)

// ─── Shared Aspen binary lookup helper ──────────────────────────────────────
/**
 * Query nistv140 table first; fall back to apv140 table.
 * Returns the raw row and which direction matches (name1 = Compound1Name).
 */
async function fetchAspenBinaryRow(
  supabase: SupabaseClient,
  nistTable: string,
  apvTable: string,
  name1: string,
  name2: string
): Promise<{ row: any; isForward: boolean } | null> {
  for (const table of [nistTable, apvTable]) {
    const { data, error } = await supabase
      .from(table)
      .select('*')
      .or(`and("Compound1Name".ilike.${name1},"Compound2Name".ilike.${name2}),and("Compound1Name".ilike.${name2},"Compound2Name".ilike.${name1})`)
      .limit(1);
    if (!error && data && data.length > 0) {
      const row = data[0];
      const isForward = row.Compound1Name?.toLowerCase() === name1.toLowerCase();
      return { row, isForward };
    }
  }
  return null;
}

// --- NRTL ---
export async function fetchNrtlParameters(
  supabase: SupabaseClient,
  name1: string,
  name2: string
): Promise<NrtlInteractionParams> {
  if (!name1 || !name2) throw new Error('Component names required for NRTL');
  if (name1 === name2) return { Aij: 0, Aji: 0, Bij: 0, Bji: 0, alpha: 0.3 };

  const result = await fetchAspenBinaryRow(supabase, 'nistv140_nrtl_wide', 'apv140_nrtl_wide', name1, name2);
  if (result) {
    const { row, isForward } = result;
    // Aspen NRTL: τ_ij = aij + bij/T  (aij dim'less, bij in K)
    // Scale to HYSYS convention (Aij[cal/mol] = bij*R_cal, Bij[cal/mol/K] = aij*R_cal)
    // so existing calc tau = Aij/(R*T) + Bij/R = bij/T + aij remains correct.
    // Extra terms: Eij = eij*R_cal (ln T), Fij = fij*R_cal (T), dij_alpha raw (symmetric alpha temp coeff).
    const [aij, aji, bij, bji, alpha, dij_alpha, eij, eji, fij, fji] = isForward
      ? [row.aij ?? 0, row.aji ?? 0, row.bij ?? 0, row.bji ?? 0, row.cij ?? 0.3,
         row.dij ?? 0, row.eij ?? 0, row.eji ?? 0, row.fij ?? 0, row.fji ?? 0]
      : [row.aji ?? 0, row.aij ?? 0, row.bji ?? 0, row.bij ?? 0, row.cij ?? 0.3,
         row.dij ?? 0, row.eji ?? 0, row.eij ?? 0, row.fji ?? 0, row.fij ?? 0]; // dij symmetric
    return {
      Aij: bij * R_cal, Aji: bji * R_cal, Bij: aij * R_cal, Bji: aji * R_cal, alpha,
      dij_alpha,
      Eij: eij * R_cal, Eji: eji * R_cal, Fij: fij * R_cal, Fji: fji * R_cal,
      eta_cross: row.HOC_ETA_CROSS ?? 0,
    };
  }
  // Fallback to UNIFAC-based estimate
  const est = await estimateNrtlFromUnifac(supabase, name1, name2);
  if (est) {
    (est as any)._usedUnifacFallback = true;
    return est;
  }
  return { Aij: 0, Aji: 0, Bij: 0, Bji: 0, alpha: 0.3 };
}

// --- Wilson ---
export async function fetchWilsonInteractionParams(
  supabase: SupabaseClient,
  name1: string,
  name2: string
): Promise<WilsonInteractionParams> {
  if (!name1 || !name2) return { Aij: 0, Aji: 0, Bij: 0, Bji: 0 };

  // Wilson: APV140 first — NIST140 Wilson params use wide-range regressions that can invert
  // the sign of non-ideality for some binaries (e.g. Acetone–Chloroform). APV VLE-HOC matches
  // the Aspen Plus default and gives the correct max-boiling azeotrope behavior.
  const result = await fetchAspenBinaryRow(supabase, 'apv140_wilson_wide', 'nistv140_wilson_wide', name1, name2);
  if (result) {
    const { row, isForward } = result;
    // Aspen Wilson: Λ_12 = (V2/V1)·exp(aij + bij/T + cij·lnT + dij·T)
    // Scale: Aij = -bij*R_cal, Bij = -aij*R_cal, Cij = -cij*R_cal, Dij = -dij*R_cal so that
    //   exp(-(Aij + Bij·T + Cij·T·lnT + Dij·T²)/(R·T)) = exp(bij/T + aij + cij·lnT + dij·T)  ✓
    const [aij, aji, bij, bji, cij, cji, dij, dji] = isForward
      ? [row.aij ?? 0, row.aji ?? 0, row.bij ?? 0, row.bji ?? 0,
         row.cij ?? 0, row.cji ?? 0, row.dij ?? 0, row.dji ?? 0]
      : [row.aji ?? 0, row.aij ?? 0, row.bji ?? 0, row.bij ?? 0,
         row.cji ?? 0, row.cij ?? 0, row.dji ?? 0, row.dij ?? 0];
    return {
      Aij: -bij * R_cal, Aji: -bji * R_cal, Bij: -aij * R_cal, Bji: -aji * R_cal,
      Cij: -cij * R_cal, Cji: -cji * R_cal, Dij: -dij * R_cal, Dji: -dji * R_cal,
      eta_cross: row.HOC_ETA_CROSS ?? 0,
    };
  }
  return { Aij: 0, Aji: 0, Bij: 0, Bji: 0 };
}

// --- Peng-Robinson ---
export async function fetchPrInteractionParams(
  supabase: SupabaseClient,
  name1: string,
  name2: string
): Promise<PrInteractionParams> {
  if (!name1 || !name2 || name1 === name2) return { k_ij: 0 };

  const result = await fetchAspenBinaryRow(supabase, 'nistv140_prkbv_wide', 'apv140_prkbv_wide', name1, name2);
  if (result && typeof result.row.kij1 === 'number') {
    return { k_ij: result.row.kij1 };
  }
  // No binary data found → use kij = 0 (standard default for PR when no regressed data available).
  // The Chueh-Prausnitz critical-volume formula is only valid for simple non-polar pairs;
  // using it for polar systems (alcohols, chlorinated compounds) gives incorrect kij.
  return { k_ij: 0 };
}

// --- SRK ---
export async function fetchSrkInteractionParams(
  supabase: SupabaseClient,
  name1: string,
  name2: string
): Promise<SrkInteractionParams> {
  if (!name1 || !name2 || name1 === name2) return { k_ij: 0 };

  // SRK uses APV rkskbv table (no NIST equivalent)
  const { data, error } = await supabase
    .from('apv140_rkskbv_wide')
    .select('*')
    .or(`and("Compound1Name".ilike.${name1},"Compound2Name".ilike.${name2}),and("Compound1Name".ilike.${name2},"Compound2Name".ilike.${name1})`)
    .limit(1);

  if (!error && data && data.length > 0 && typeof data[0].kij1 === 'number') {
    return { k_ij: data[0].kij1 };
  }
  // No binary data found → use kij = 0 (standard SRK default; CVC formula invalid for polar pairs).
  return { k_ij: 0 };
}

// --- UNIQUAC ---
export async function fetchUniquacInteractionParams(
  supabase: SupabaseClient,
  name1: string,
  name2: string
): Promise<UniquacInteractionParams> {
  if (!name1 || !name2) return { Aij: 0, Aji: 0, Bij: 0, Bji: 0 };

  const result = await fetchAspenBinaryRow(supabase, 'nistv140_uniq_wide', 'apv140_uniq_wide', name1, name2);
  if (result) {
    const { row, isForward } = result;
    // Aspen UNIQUAC: τ_ij = exp(aij + bij/T + cij·lnT + dij·T)
    // Scale: Aij = -bij*R_cal, Bij = -aij*R_cal, Cij = -cij*R_cal, Dij = -dij*R_cal so that
    //   exp(−(Aij + Bij·T + Cij·T·lnT + Dij·T²)/(R·T)) = exp(bij/T + aij + cij·lnT + dij·T)  ✓
    const [aij, aji, bij, bji, cij, cji, dij, dji] = isForward
      ? [row.aij ?? 0, row.aji ?? 0, row.bij ?? 0, row.bji ?? 0,
         row.cij ?? 0, row.cji ?? 0, row.dij ?? 0, row.dji ?? 0]
      : [row.aji ?? 0, row.aij ?? 0, row.bji ?? 0, row.bij ?? 0,
         row.cji ?? 0, row.cij ?? 0, row.dji ?? 0, row.dij ?? 0];
    return {
      Aij: -bij * R_cal, Aji: -bji * R_cal, Bij: -aij * R_cal, Bji: -aji * R_cal,
      Cij: -cij * R_cal, Cji: -cji * R_cal, Dij: -dij * R_cal, Dji: -dji * R_cal,
      eta_cross: row.HOC_ETA_CROSS ?? 0,
    };
  }
  // Fallback using UNIFAC
  const est = await estimateUniquacFromUnifac(supabase, name1, name2);
  if (est) {
    (est as any)._usedUnifacFallback = true;
    return est;
  }
  return { Aij: 0, Aji: 0, Bij: 0, Bji: 0 };
}

// =============================
//  8. Re-export fetch helpers (Wilson, PR, SRK, UNIQUAC) – restored
// ============================= 

// All legacy functions have been inlined - no more placeholders needed 

// --------------------------------------------------------
//  Helpers: Estimate cubic EOS binary interaction parameter (kij)
//           from critical volumes (ChemSep1/2 only)
// --------------------------------------------------------



// --------------------------------------------------------
//  Helpers: Estimate NRTL / UNIQUAC binary parameters via UNIFAC
// --------------------------------------------------------

async function _buildUnifacComponent(supabase: SupabaseClient, name: string): Promise<CompoundData | null> {
  const therm = await fetchCompoundDataFromHysys(supabase, name).catch(() => null);
  if (!therm || !therm.unifacGroups) return null;
  return {
    name: name,
    antoine: null,
    unifacGroups: therm.unifacGroups,
    prParams: null,
    srkParams: null,
    uniquacParams: null,
    wilsonParams: null,
  } as CompoundData;
}

async function estimateNrtlFromUnifac(
  supabase: SupabaseClient,
  name1: string,
  name2: string,
  T_ref_K: number = 350
): Promise<NrtlInteractionParams | null> {
  const comp1 = await _buildUnifacComponent(supabase, name1);
  const comp2 = await _buildUnifacComponent(supabase, name2);
  if (!comp1 || !comp2) return null;
  const subgroupIds = Array.from(new Set([
    ...Object.keys(comp1.unifacGroups || {}).map(Number),
    ...Object.keys(comp2.unifacGroups || {}).map(Number),
  ]));
  if (subgroupIds.length === 0) return null;
  const unifacParams = await fetchUnifacInteractionParams(supabase, subgroupIds).catch(() => null);
  if (!unifacParams) return null;

  // Compute gamma at equimolar composition
  const gammas = calculateUnifacGamma([comp1, comp2], [0.5, 0.5], T_ref_K, unifacParams);
  if (!gammas) return null;
  const lnGamma1 = Math.log(gammas[0]);
  const lnGamma2 = Math.log(gammas[1]);
  // Heuristic: store as Bij (temperature-independent), so Aij=0
  const Bij12 = lnGamma2; // Cross-assignment
  const Bij21 = lnGamma1;
  const alpha = 0.3;
  return { Aij: 0, Aji: 0, Bij: Bij12, Bji: Bij21, alpha };
}

async function estimateUniquacFromUnifac(
  supabase: SupabaseClient,
  name1: string,
  name2: string,
  T_ref_K: number = 350
): Promise<UniquacInteractionParams | null> {
  const comp1 = await _buildUnifacComponent(supabase, name1);
  const comp2 = await _buildUnifacComponent(supabase, name2);
  if (!comp1 || !comp2) return null;

  const subgroupIds = Array.from(new Set([
    ...Object.keys(comp1.unifacGroups || {}).map(Number),
    ...Object.keys(comp2.unifacGroups || {}).map(Number),
  ]));
  if (subgroupIds.length === 0) return null;

  const unifacParams = await fetchUnifacInteractionParams(supabase, subgroupIds).catch(() => null);
  if (!unifacParams) return null;

  const gammas = calculateUnifacGamma([comp1, comp2], [0.5, 0.5], T_ref_K, unifacParams);
  if (!gammas) return null;
  const lnGamma1 = Math.log(gammas[0]);
  const lnGamma2 = Math.log(gammas[1]);
  // Heuristic: store as Bij (temperature-independent), Aij=0
  const Bij12 = lnGamma2;
  const Bij21 = lnGamma1;
  return { Aij: 0, Aji: 0, Bij: Bij12, Bji: Bij21 };
}

// ==========================================================
//  9. Pure-component saturation temperature (PR & SRK EOS)
// ==========================================================

/**
 * Estimate the saturation (boiling) temperature of a pure component at a
 * specified pressure using the Peng–Robinson equation of state.
 *
 * Returns `null` if the calculation fails (e.g. outside two-phase region).
 */
/**
 * Estimate the saturation (boiling) temperature of a pure component at a
 * specified pressure using the Peng–Robinson equation of state.
 *
 * THIS IS THE CORRECTED, MORE ROBUST VERSION.
 */
export function calculateSaturationTemperaturePurePr(
  component: CompoundData,
  P_Pa: number,
  maxIter: number = 20,
  tolerance: number = 1e-5
): number | null {
  if (!component?.prParams || !component.antoine) return null;
  const { Tc_K, Pc_Pa, omega } = component.prParams;

  if (P_Pa >= Pc_Pa * 0.99) { // Avoid calculations too close to critical point
      return antoineBoilingPointSolverLocal(component.antoine, P_Pa);
  }

  // Use Antoine equation for a very good initial guess
  let T = antoineBoilingPointSolverLocal(component.antoine, P_Pa);
  if (T === null) T = 0.7 * Tc_K; // Fallback guess if Antoine fails

  const R = R_gas_const_J_molK;

  for (let i = 0; i < maxIter; i++) {
    const Tr = T / Tc_K;
    if (T <= 0 || Tr >= 1) break;

    let kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    if (omega > 0.49) {
      kappa = 0.379642 + 1.48503 * omega - 0.164423 * (omega * omega) + 0.016666 * (omega * omega * omega);
    }
    const alpha = (1 + kappa * (1 - Math.sqrt(Tr))) ** 2;
    const a = 0.457235 * Math.pow(R * Tc_K, 2) / Pc_Pa * alpha;
    const b = 0.07796 * R * Tc_K / Pc_Pa;

    const A = a * P_Pa / Math.pow(R * T, 2);
    const B = b * P_Pa / (R * T);

    const roots = solveCubicEOS(1, -(1 - B), A - 3 * B * B - 2 * B, -(A * B - B * B - B * B * B));
    if (!roots || roots.length < 2) {
        T += 1.0; // Nudge temperature up if we lose two-phase region
        continue;
    }
    const Z_L = Math.min(...roots);
    const Z_V = Math.max(...roots);

    const sqrt2 = Math.SQRT2;
    const lnPhi = (Z: number): number => (Z - 1 - Math.log(Z - B) - (A / (2 * sqrt2 * B)) * Math.log((Z + (1 + sqrt2) * B) / (Z + (1 - sqrt2) * B)));
    
    const f_T = lnPhi(Z_L) - lnPhi(Z_V); // Objective function

    if (Math.abs(f_T) < tolerance) {
        return T; // Converged
    }

    // Estimate derivative and update T via Newton-Raphson step
    const dT = 0.001 * T;
    const T_plus_dT = T + dT;
    const Tr_plus_dT = T_plus_dT / Tc_K;
    const alpha_plus_dT = (1 + kappa * (1 - Math.sqrt(Tr_plus_dT))) ** 2;
    const a_plus_dT = 0.457235 * Math.pow(R * Tc_K, 2) / Pc_Pa * alpha_plus_dT;
    const A_plus_dT = a_plus_dT * P_Pa / Math.pow(R * T_plus_dT, 2);
    const B_plus_dT = b * P_Pa / (R * T_plus_dT);
    const roots_plus_dT = solveCubicEOS(1, -(1 - B_plus_dT), A_plus_dT - 3 * B_plus_dT * B_plus_dT - 2 * B_plus_dT, -(A_plus_dT * B_plus_dT - B_plus_dT * B_plus_dT - B_plus_dT * B_plus_dT * B_plus_dT));

    if (!roots_plus_dT || roots_plus_dT.length < 2) {
      T -= f_T * 5; // Fallback update if derivative calculation fails
      continue;
    }
    const Z_L_plus_dT = Math.min(...roots_plus_dT);
    const Z_V_plus_dT = Math.max(...roots_plus_dT);
    const f_T_plus_dT = lnPhi(Z_L_plus_dT) - lnPhi(Z_V_plus_dT);
    
    const df_dT = (f_T_plus_dT - f_T) / dT;
    if (Math.abs(df_dT) < 1e-9) break; // Avoid division by zero
    
    T -= f_T / df_dT; // Newton's step
  }

  // If EOS solver fails, return Antoine result as a reliable fallback
  return antoineBoilingPointSolverLocal(component.antoine, P_Pa);
}

/**
 * Estimate saturation temperature using the Soave–Redlich–Kwong EOS.
 *
 * THIS IS THE CORRECTED, MORE ROBUST VERSION.
 */
export function calculateSaturationTemperaturePureSrk(
  component: CompoundData,
  P_Pa: number,
  maxIter: number = 20,
  tolerance: number = 1e-5
): number | null {
  if (!component?.srkParams || !component.antoine) return null;
  const { Tc_K, Pc_Pa, omega } = component.srkParams;

  if (P_Pa >= Pc_Pa * 0.99) {
      return antoineBoilingPointSolverLocal(component.antoine, P_Pa);
  }

  // Use Antoine equation for a very good initial guess
  let T = antoineBoilingPointSolverLocal(component.antoine, P_Pa);
  if (T === null) T = 0.7 * Tc_K;

  const R = R_gas_const_J_molK;

  for (let i = 0; i < maxIter; i++) {
    const Tr = T / Tc_K;
    if (T <= 0 || Tr >= 1) break;

    const m = 0.480 + 1.574 * omega - 0.176 * omega * omega;
    const alpha = Math.pow(1 + m * (1 - Math.sqrt(Tr)), 2);
    const a = 0.42748 * Math.pow(R * Tc_K, 2) / Pc_Pa * alpha;
    const b = 0.08664 * R * Tc_K / Pc_Pa;

    const A = a * P_Pa / Math.pow(R * T, 2);
    const B = b * P_Pa / (R * T);

    const roots = solveCubicEOS(1, -1, A - B - B * B, -A * B);
    if (!roots || roots.length < 2) {
      T += 1.0;
      continue;
    }
    const Z_L = Math.min(...roots);
    const Z_V = Math.max(...roots);

    const lnPhi = (Z: number): number => (Z - 1) - Math.log(Z - B) - (A / B) * Math.log(1 + B / Z);
    const f_T = lnPhi(Z_L) - lnPhi(Z_V);

    if (Math.abs(f_T) < tolerance) {
      return T;
    }

    const dT = 0.001 * T;
    const T_plus_dT = T + dT;
    const Tr_plus_dT = T_plus_dT / Tc_K;
    const alpha_plus_dT = Math.pow(1 + m * (1 - Math.sqrt(Tr_plus_dT)), 2);
    const a_plus_dT = 0.42748 * Math.pow(R * Tc_K, 2) / Pc_Pa * alpha_plus_dT;
    const A_plus_dT = a_plus_dT * P_Pa / Math.pow(R * T_plus_dT, 2);
    const B_plus_dT = b * P_Pa / (R * T_plus_dT);
    const roots_plus_dT = solveCubicEOS(1, -1, A_plus_dT - B_plus_dT - B_plus_dT * B_plus_dT, -A_plus_dT * B_plus_dT);

    if (!roots_plus_dT || roots_plus_dT.length < 2) {
        T -= f_T * 5;
        continue;
    }
    const Z_L_plus_dT = Math.min(...roots_plus_dT);
    const Z_V_plus_dT = Math.max(...roots_plus_dT);
    const f_T_plus_dT = lnPhi(Z_L_plus_dT) - lnPhi(Z_V_plus_dT);
    const df_dT = (f_T_plus_dT - f_T) / dT;
    
    if (Math.abs(df_dT) < 1e-9) break;

    T -= f_T / df_dT;
  }
  
  return antoineBoilingPointSolverLocal(component.antoine, P_Pa);
}

// --------------------------------------------------------------
// Generic Brent's method root-finder for 1D scalar functions
// --------------------------------------------------------------
function _brentRoot(
  f: (x: number) => number,
  a: number,
  b: number,
  tolAbs: number,
  maxIter = 100
): number | null {
  let fa = f(a);
  let fb = f(b);
  if (!isFinite(fa) || !isFinite(fb)) return null;
  if (fa * fb > 0) return null; // not bracketed

  if (Math.abs(fa) < Math.abs(fb)) {
    [a, b] = [b, a];
    [fa, fb] = [fb, fa];
  }

  let c = a, fc = fa, d = b - a, e = d;

  for (let iter = 0; iter < maxIter; iter++) {
    if (Math.abs(fc) < Math.abs(fb)) {
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }

    const tol = 2 * Number.EPSILON * Math.abs(b) + tolAbs;
    const m = 0.5 * (c - b);
    // Fixed: Check BOTH function value AND step size
    if (Math.abs(fb) < tolAbs || Math.abs(m) <= tol) {
      return b;
    }

    if (Math.abs(e) >= tol && Math.abs(fa) > Math.abs(fb)) {
      // Attempt inverse quadratic / secant step
      let s = fb / fa;
      let p: number, q: number;
      if (a === c) {
        // secant
        p = 2 * m * s;
        q = 1 - s;
      } else {
        // inverse quadratic
        const r = fc / fa;
        const t = fb / fc;
        p = s * (2 * m * r * (r - t) - (b - a) * (t - 1));
        q = (r - 1) * (t - 1) * (s - 1);
      }
      if (p > 0) q = -q;
      p = Math.abs(p);

      const min1 = 3 * m * q - Math.abs(tol * q);
      const min2 = Math.abs(e * q);
      if (2 * p < (min1 < min2 ? min1 : min2)) {
        e = d;
        d = p / q;
      } else {
        d = m;
        e = d;
      }
    } else {
      // Bisection
      d = m;
      e = d;
    }

    a = b;
    fa = fb;
    if (Math.abs(d) > tol) {
      b += d;
    } else {
      b += m > 0 ? tol : -tol;
    }
    fb = f(b);
    if (!isFinite(fb)) return null;
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      c = a; fc = fa; e = d = b - a;
    }
  }
  return null; // did not converge
}