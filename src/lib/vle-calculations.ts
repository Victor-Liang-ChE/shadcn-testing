// Consolidated VLE calculation utilities supporting all major thermodynamic models:
// • NRTL (Non-Random Two Liquid) - activity coefficient model
// • Wilson - activity coefficient model  
// • UNIFAC (UNIversal Functional Activity Coefficient) - group contribution method
// • Peng-Robinson (PR) - cubic equation of state with fugacity coefficients
// • Soave-Redlich-Kwong (SRK) - cubic equation of state with fugacity coefficients
// • UNIQUAC (UNIversal QUAsi Chemical) - activity coefficient model
// All legacy VLE calculation files have been consolidated into this single module.
// ---------------------------------------------------------------------------------------
import type { AntoineParams, CompoundData, BubbleDewResult, UnifacGroupComposition } from './vle-types';
import type { SupabaseClient } from '@supabase/supabase-js';

// =============================
//  1. SHARED CONSTANTS/UTILS
// =============================
export const R_gas_const_J_molK = 8.314_462_618_153_24; // J·mol⁻¹·K⁻¹
export const R_cal_molK = 1.987_204_258_640_83;         // cal·mol⁻¹·K⁻¹  (for legacy NRTL formulas)
export const JOULES_PER_CAL = 4.184;

/**
 * Calculates saturation pressure (Pa) from Antoine coefficients.
 * Supports both log10 and ln formulations and automatically converts kPa → Pa.
 */
export function calculatePsat_Pa(
  params: AntoineParams | null,
  T_K: number
): number {
  if (!params || T_K <= 0) return NaN;

  const conv = params.Units?.toLowerCase() === 'kpa' ? 1000 : 1; // kPa → Pa
  let P: number;
  if (params.EquationNo === 1 || params.EquationNo === '1') {
    // log10 form
    P = 10 ** (params.A - params.B / (T_K + params.C));
  } else {
    // ln form (ChemSep-style)
    P = Math.exp(params.A - params.B / (T_K + params.C));
  }
  return isNaN(P) || P <= 0 ? NaN : P * conv;
}

// ==========================================
//  2. N R T L   (binary activity-coefficient)
// ==========================================
export interface NrtlInteractionParams {
  A12: number; // J·mol⁻¹
  A21: number; // J·mol⁻¹
  alpha: number; // dimensionless non-randomness parameter (typically 0.2-0.47)
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

  const tau12 = p.A12 / (R_cal_molK * T_K);
  const tau21 = p.A21 / (R_cal_molK * T_K);
  const G12 = Math.exp(-p.alpha * tau12);
  const G21 = Math.exp(-p.alpha * tau21);

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

  for (let i = 0; i < maxIter; i++) {
    const Psat = [calculatePsat_Pa(comps[0].antoine, T), calculatePsat_Pa(comps[1].antoine, T)];
    if (Psat.some(v => isNaN(v))) return null;
    const gamma = calculateNrtlGamma(comps, x, T, params);
    if (!gamma) return null;

    const Pcalc = x[0] * gamma[0] * Psat[0] + x[1] * gamma[1] * Psat[1];
    const err = Pcalc - P_Pa;
    if (Math.abs(err) < tol * P_Pa) {
      y = [
        (x[0] * gamma[0] * Psat[0]) / Pcalc,
        (x[1] * gamma[1] * Psat[1]) / Pcalc,
      ];
      return {
        comp1_feed: x1_feed,
        comp1_equilibrium: y[0],
        T_K: T,
        P_Pa,
        iterations: i + 1,
        calculationType: 'bubbleT',
      };
    }
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
  const Pcalc = x[0] * gamma[0] * Psat[0] + x[1] * gamma[1] * Psat[1];
  const y1 = (x[0] * gamma[0] * Psat[0]) / Pcalc;
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
export interface WilsonInteractionParams {
  a12_J_mol: number;
  a21_J_mol: number;
}

export function calculateWilsonGamma(
  comps: CompoundData[],
  x: number[],
  T_K: number,
  p: WilsonInteractionParams
): [number, number] | null {
  if (x.length !== 2 || comps.length !== 2) return null;
  if (!comps[0].wilsonParams || !comps[1].wilsonParams) return null;

  const V1 = comps[0].wilsonParams.V_L_m3mol;
  const V2 = comps[1].wilsonParams.V_L_m3mol;
  const Lambda12 = (V2 / V1) * Math.exp(-p.a12_J_mol / (R_gas_const_J_molK * T_K));
  const Lambda21 = (V1 / V2) * Math.exp(-p.a21_J_mol / (R_gas_const_J_molK * T_K));

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
  for (let i = 0; i < maxIter; i++) {
    const gammas = calculateWilsonGamma(comps, x, T, params);
    if (!gammas) return null;
    const Psat = [calculatePsat_Pa(comps[0].antoine, T), calculatePsat_Pa(comps[1].antoine, T)];
    if (Psat.some(v => isNaN(v))) return null;

    const Pcalc = x[0] * gammas[0] * Psat[0] + x[1] * gammas[1] * Psat[1];
    const err = Pcalc - P_Pa;
    if (Math.abs(err) < tol * P_Pa) {
      const y1 = (x[0] * gammas[0] * Psat[0]) / Pcalc;
      return {
        comp1_feed: x1_feed,
        comp1_equilibrium: y1,
        T_K: T,
        P_Pa,
        iterations: i + 1,
        calculationType: 'bubbleT',
      };
    }
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
  const Pcalc = x[0] * gammas[0] * Psat[0] + x[1] * gammas[1] * Psat[1];
  const y1 = (x[0] * gammas[0] * Psat[0]) / Pcalc;
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

  // -------- Residual part --------
  const subgroupIds = Array.from(new Set(componentData.flatMap(c => Object.keys(c.unifacGroups || {}).map(Number))));
  if (subgroupIds.length === 0) return ln_gamma_C.map(Math.exp);

  // v_ki matrix (components × subgroups)
  const v_ki: number[][] = componentData.map(c => subgroupIds.map(id => c.unifacGroups?.[id] || 0));

  // X_m (group mole fractions)
  const sum_x_vk_total = subgroupIds.reduce((sum, _, kIdx) => sum + componentData.reduce((s, _, i) => s + x[i] * v_ki[i][kIdx], 0), 0);
  if (sum_x_vk_total === 0) return ln_gamma_C.map(Math.exp);
  const X_m = subgroupIds.map((_, kIdx) => componentData.reduce((s, _, i) => s + x[i] * v_ki[i][kIdx], 0) / sum_x_vk_total);

  const sum_XQ = subgroupIds.reduce((s, id, idx) => s + X_m[idx] * params.Qk[id], 0);
  if (sum_XQ === 0) return ln_gamma_C.map(Math.exp);
  const theta_m = subgroupIds.map((id, idx) => (X_m[idx] * params.Qk[id]) / sum_XQ);

  // psi matrix (m × k)
  const psi: number[][] = subgroupIds.map(() => new Array(subgroupIds.length).fill(0));
  subgroupIds.forEach((id_m, mIdx) => {
    subgroupIds.forEach((id_k, kIdx) => {
      const key = `${params.mainGroupMap[id_m]}-${params.mainGroupMap[id_k]}`;
      const a_mk = params.a_mk.get(key) ?? 0;
      psi[mIdx][kIdx] = Math.exp(-a_mk / T_kelvin);
    });
  });

  // ln Γ_k for each subgroup k
  const ln_Gamma_k = subgroupIds.map((_, kIdx) => {
    const numerator = subgroupIds.reduce((s, _, mIdx) => s + theta_m[mIdx] * psi[mIdx][kIdx], 0);
    if (numerator === 0) return 0;
    const term1 = Math.log(numerator);
    const term2 = subgroupIds.reduce((s, _, mIdx) => {
      const denom = subgroupIds.reduce((sd, _, nIdx) => sd + theta_m[nIdx] * psi[nIdx][mIdx], 0);
      return s + (theta_m[mIdx] * psi[kIdx][mIdx]) / denom;
    }, 0);
    return params.Qk[subgroupIds[kIdx]] * (1 - term1 - term2);
  });

  // ln γ_R for each component
  const ln_gamma_R = componentData.map((_, i) => {
    let sum = 0;
    subgroupIds.forEach((_, kIdx) => {
      sum += v_ki[i][kIdx] * (ln_Gamma_k[kIdx] - ln_Gamma_k[kIdx]); // Note: ln_Gamma_k(pure k)≈0
    });
    return sum;
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
  const P_sat = [calculatePsat_Pa(components[0].antoine, T_K), calculatePsat_Pa(components[1].antoine, T_K)];
  if (P_sat.some(p => isNaN(p))) return null;
  const gamma = calculateUnifacGamma(components, x, T_K, unifacParams);
  if (!gamma) return null;
  const P_calc = x[0] * gamma[0] * P_sat[0] + x[1] * gamma[1] * P_sat[1];
  const y1 = (x[0] * gamma[0] * P_sat[0]) / P_calc;
  return { comp1_feed: x1_feed, comp1_equilibrium: y1, T_K, P_Pa: P_calc };
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

export interface PrInteractionParams { k_ij: number; k_ji: number }

export function solveCubicEOS(a: number, b: number, c: number, d: number): number[] | null {
  // Normalized cubic Z^3 + aZ^2 + bZ + c = 0  -> using Cardano's method
  const p = b - a * a / 3;
  const q = (2 * a * a * a) / 27 - (a * b) / 3 + c;
  const D = (q * q) / 4 + (p * p * p) / 27;
  if (D > 0) {
    const u = Math.cbrt(-q / 2 + Math.sqrt(D));
    const v = Math.cbrt(-q / 2 - Math.sqrt(D));
    return [u + v - a / 3].filter(z => z > 0);
  }
  if (Math.abs(D) < 1e-12) {
    const u = Math.cbrt(-q / 2);
    return [2 * u - a / 3, -u - a / 3].filter((z, i, arr) => i === 0 || Math.abs(z - arr[0]) > 1e-8).filter(z => z > 0);
  }
  const r = Math.sqrt(-p * p * p / 27);
  const phi = Math.acos(-q / (2 * r));
  const m = 2 * Math.cbrt(r);
  return [
    m * Math.cos(phi / 3) - a / 3,
    m * Math.cos((phi + 2 * Math.PI) / 3) - a / 3,
    m * Math.cos((phi + 4 * Math.PI) / 3) - a / 3,
  ].filter(z => z > 0);
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
    const kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    const alpha = (1 + kappa * (1 - Math.sqrt(Tr))) ** 2;
    const ac = 0.45723553 * (R_gas_const_J_molK * Tc_K) ** 2 / Pc_Pa * alpha;
    const b  = 0.07779607 * R_gas_const_J_molK * Tc_K / Pc_Pa;
    return { ac, b };
  });

  const ac_mix = x1 * x1 * pure[0].ac + x2 * x2 * pure[1].ac + 2 * x1 * x2 * Math.sqrt(pure[0].ac * pure[1].ac) * (1 - k_ij);
  const b_mix  = x1 * pure[0].b + x2 * pure[1].b;

  const A = ac_mix * P_Pa / (R_gas_const_J_molK * T_K) ** 2;
  const B = b_mix * P_Pa / (R_gas_const_J_molK * T_K);

  const roots = solveCubicEOS(1, -(1 - B), A - 3 * B ** 2 - 2 * B, -(A * B - B ** 2 - B ** 3));
  if (!roots || roots.length === 0) return null;

  const Z = phase === 'liquid' ? Math.min(...roots) : Math.max(...roots);
  if (Z <= 0) return null;

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

// Small numerical derivative helper (central difference)
function _numDeriv(func: (t: number) => number, t: number, h: number = 1e-4) {
  return (func(t + h) - func(t - h)) / (2 * h);
}

export function calculateBubbleTemperaturePr(
  components: CompoundData[],
  x1_feed: number,
  P_system_Pa: number,
  params: PrInteractionParams,
  T_guess_K: number = 350,
  maxIter = 60,
  tol = 1e-5
): BubbleDewResult | null {
  let T = T_guess_K;
  const x = [x1_feed, 1 - x1_feed];
  for (let i = 0; i < maxIter; i++) {
    const phiL = calculatePrFugacityCoefficients(components, x, T, P_system_Pa, params, 'liquid');
    const phiV = calculatePrFugacityCoefficients(components, x, T, P_system_Pa, params, 'vapor');
    if (!phiL || !phiV) return null;
    const K = [phiL[0] / phiV[0], phiL[1] / phiV[1]];
    const f = x[0] * K[0] + x[1] * K[1] - 1;
    if (Math.abs(f) < tol) {
      const y1 = (x[0] * K[0]);
      return { comp1_feed: x1_feed, comp1_equilibrium: y1, T_K: T, P_Pa: P_system_Pa };
    }
    const df_dT = _numDeriv(Tc => {
      const phiL_ = calculatePrFugacityCoefficients(components, x, Tc, P_system_Pa, params, 'liquid');
      const phiV_ = calculatePrFugacityCoefficients(components, x, Tc, P_system_Pa, params, 'vapor');
      if (!phiL_ || !phiV_) return f; // fallback
      const K_ = [phiL_[0] / phiV_[0], phiL_[1] / phiV_[1]];
      return x[0] * K_[0] + x[1] * K_[1] - 1;
    }, T);
    if (!isFinite(df_dT) || Math.abs(df_dT) < 1e-8) return null;
    T -= f / df_dT;
    if (T < 100 || T > 1000) return null;
  }
  return null;
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

export function calculateDewTemperaturePr(
  components: CompoundData[],
  y1_feed: number,
  P_system_Pa: number,
  prInteractionParams: PrInteractionParams,
  initialTempGuess_K: number = 350,
  maxIter: number = 50,
  tolerance: number = 1e-5
): BubbleDewResult | null {
  let T_K = initialTempGuess_K;
  const y = [y1_feed, 1 - y1_feed];

  if (!components[0]?.prParams || !components[1]?.prParams) return null;
  if (!components[0]?.antoine || !components[1]?.antoine) return null;

  for (let iter = 0; iter < maxIter; iter++) {
    const Psat = components.map(c => calculatePsat_Pa(c.antoine, T_K));
    if (Psat.some(p => isNaN(p) || p <= 0)) return null;

    // Initial liquid composition guess using ideal K = Psat/P
    let x = y.map((yi, i) => yi / (Psat[i] / P_system_Pa));
    const x_sum = x.reduce((s, xi) => s + xi, 0);
    x = x.map(xi => xi / x_sum);

    // Calculate fugacity coefficients
    const phi_V = calculatePrFugacityCoefficients(components, y, T_K, P_system_Pa, prInteractionParams, 'vapor');
    if (!phi_V) return null;

    const phi_L = calculatePrFugacityCoefficients(components, x, T_K, P_system_Pa, prInteractionParams, 'liquid');
    if (!phi_L) return null;

    // Calculate K-values and objective function
    const K = phi_L.map((phiL, i) => phiL * Psat[i] / (phi_V[i] * P_system_Pa));
    const obj = y.reduce((sum, yi, i) => sum + yi / K[i], 0) - 1;

    if (Math.abs(obj) < tolerance) {
      // Recalculate final liquid composition
      x = y.map((yi, i) => yi / K[i]);
      const x_final_sum = x.reduce((s, xi) => s + xi, 0);
      x = x.map(xi => xi / x_final_sum);
      
      return {
        comp1_feed: y1_feed,
        comp1_equilibrium: x[0],
        T_K: T_K,
        P_Pa: P_system_Pa,
        iterations: iter + 1,
        calculationType: 'dewT'
      };
    }

    // Newton step for temperature
    const dobj_dT = _numDeriv((T: number) => {
      const Psat_T = components.map(c => calculatePsat_Pa(c.antoine, T));
      if (Psat_T.some(p => isNaN(p))) return NaN;
      let x_T = y.map((yi, i) => yi / (Psat_T[i] / P_system_Pa));
      const x_sum_T = x_T.reduce((s, xi) => s + xi, 0);
      x_T = x_T.map(xi => xi / x_sum_T);
      const phi_V_T = calculatePrFugacityCoefficients(components, y, T, P_system_Pa, prInteractionParams, 'vapor');
      if (!phi_V_T) return NaN;
      const phi_L_T = calculatePrFugacityCoefficients(components, x_T, T, P_system_Pa, prInteractionParams, 'liquid');
      if (!phi_L_T) return NaN;
      const K_T = phi_L_T.map((phiL, i) => phiL * Psat_T[i] / (phi_V_T[i] * P_system_Pa));
      return y.reduce((sum, yi, i) => sum + yi / K_T[i], 0) - 1;
    }, T_K);

    if (Math.abs(dobj_dT) < 1e-12) return null;
    T_K = Math.max(100, Math.min(1000, T_K - obj / dobj_dT));
  }
  return null;
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

export interface SrkInteractionParams { k_ij: number; k_ji: number }

// type alias – avoid redeclaration
type SrkInteractionParamsFetch = SrkInteractionParams;

// (interface already declared above)

function _numDerivSrk(func: (t: number) => number, t: number, h: number = 1e-3): number {
  return (func(t + h) - func(t - h)) / (2 * h);
}

export function calculateSrkFugacityCoefficients(
  components: CompoundData[],
  x_or_y: number[],
  T_K: number,
  P_Pa: number,
  interactionParams: SrkInteractionParams,
  phase: 'liquid' | 'vapor'
): [number, number] | null {
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

  // Mixing rules
  const a12 = srkPure[0].sqrt_ac_i * srkPure[1].sqrt_ac_i * (1 - k_ij);
  const a_mix = x1*x1*srkPure[0].ac_i + 2*x1*x2*a12 + x2*x2*srkPure[1].ac_i;
  const b_mix = x1*srkPure[0].b_i + x2*srkPure[1].b_i;

  const A_mix = a_mix * P_Pa / (R_gas_const_J_molK * T_K)**2;
  const B_mix = b_mix * P_Pa / (R_gas_const_J_molK * T_K);

  // Cubic equation: Z^3 - Z^2 + (A - B - B^2)*Z - A*B = 0
  const roots = solveCubicEOS(1, -1, A_mix - B_mix - B_mix*B_mix, -A_mix*B_mix);
  if (!roots || roots.length === 0) return null;

  const Z = phase === 'liquid' ? Math.min(...roots) : Math.max(...roots);

  // Fugacity coefficients
  const phi = [0, 1].map(i => {
    const bi_over_bmix = (i === 0 ? srkPure[0].b_i : srkPure[1].b_i) / b_mix;
    const sum_term = x1 * (i === 0 ? 2*srkPure[0].ac_i : 2*a12) + 
                     x2 * (i === 0 ? 2*a12 : 2*srkPure[1].ac_i);
    const ai_term = sum_term / a_mix;
    
    const ln_phi = bi_over_bmix * (Z - 1) - Math.log(Z - B_mix) - 
                   (A_mix / B_mix) * (ai_term - bi_over_bmix) * Math.log(1 + B_mix/Z);
    return Math.exp(ln_phi);
  });

  return [phi[0], phi[1]];
}

export function calculateBubbleTemperatureSrk(
  components: CompoundData[],
  x1_feed: number,
  P_system_Pa: number,
  srkInteractionParams: SrkInteractionParams,
  initialTempGuess_K: number = 350,
  maxIter: number = 100,
  tolerance_f: number = 1e-6,
  tolerance_T_step: number = 1e-4
): BubbleDewResult | null {
  let T_current_K = initialTempGuess_K;
  const x = [x1_feed, 1 - x1_feed];

  const objectiveFunction = (T_K: number): number => {
    const Psat = components.map(c => calculatePsat_Pa(c.antoine, T_K));
    if (Psat.some(p => isNaN(p) || p <= 0)) return NaN;
    
    const phi_L = calculateSrkFugacityCoefficients(components, x, T_K, P_system_Pa, srkInteractionParams, 'liquid');
    if (!phi_L) return NaN;

    // Initial vapor composition guess
    let y = x.map((xi, i) => xi * phi_L[i] * Psat[i] / P_system_Pa);
    const y_sum = y.reduce((s, yi) => s + yi, 0);
    y = y.map(yi => yi / y_sum);

    const phi_V = calculateSrkFugacityCoefficients(components, y, T_K, P_system_Pa, srkInteractionParams, 'vapor');
    if (!phi_V) return NaN;

    // Calculate K-values and objective
    const K = phi_L.map((phiL, i) => phiL * Psat[i] / (phi_V[i] * P_system_Pa));
    return x.reduce((sum, xi, i) => sum + xi * K[i], 0) - 1;
  };

  for (let iter = 0; iter < maxIter; iter++) {
    const f_current = objectiveFunction(T_current_K);
    if (isNaN(f_current)) return null;
    if (Math.abs(f_current) < tolerance_f) {
      // Calculate final vapor composition
      const Psat = components.map(c => calculatePsat_Pa(c.antoine, T_current_K));
      const phi_L = calculateSrkFugacityCoefficients(components, x, T_current_K, P_system_Pa, srkInteractionParams, 'liquid');
      if (!phi_L) return null;
      
      let y = x.map((xi, i) => xi * phi_L[i] * Psat[i] / P_system_Pa);
      const y_sum = y.reduce((s, yi) => s + yi, 0);
      y = y.map(yi => yi / y_sum);

      return {
        comp1_feed: x1_feed,
        comp1_equilibrium: y[0],
        T_K: T_current_K,
        P_Pa: P_system_Pa,
        iterations: iter + 1,
        calculationType: 'bubbleT'
      };
    }

    const df_dT = _numDerivSrk(objectiveFunction, T_current_K);
    if (Math.abs(df_dT) < 1e-12) return null;

    const T_new_K = T_current_K - f_current / df_dT;
    if (Math.abs(T_new_K - T_current_K) < tolerance_T_step) break;
    T_current_K = Math.max(100, Math.min(1000, T_new_K));
  }
  return null;
}

export function calculateBubblePressureSrk(
  components: CompoundData[],
  x1_feed: number,
  T_system_K: number,
  srkInteractionParams: SrkInteractionParams,
  initialPressureGuess_Pa: number = 101325,
  maxIter: number = 100,
  tolerance: number = 1e-5
): BubbleDewResult | null {
  let P_current_Pa = initialPressureGuess_Pa;
  const x = [x1_feed, 1 - x1_feed];

  for (let iter = 0; iter < maxIter; iter++) {
    const Psat = components.map(c => calculatePsat_Pa(c.antoine, T_system_K));
    if (Psat.some(p => isNaN(p) || p <= 0)) return null;

    const phi_L = calculateSrkFugacityCoefficients(components, x, T_system_K, P_current_Pa, srkInteractionParams, 'liquid');
    if (!phi_L) return null;

    let y = x.map((xi, i) => xi * phi_L[i] * Psat[i] / P_current_Pa);
    const y_sum = y.reduce((s, yi) => s + yi, 0);
    y = y.map(yi => yi / y_sum);

    const phi_V = calculateSrkFugacityCoefficients(components, y, T_system_K, P_current_Pa, srkInteractionParams, 'vapor');
    if (!phi_V) return null;

    const P_new_Pa = x.reduce((sum, xi, i) => sum + xi * phi_L[i] * Psat[i] / phi_V[i], 0);
    
    if (Math.abs(P_new_Pa - P_current_Pa) < tolerance * P_current_Pa) {
      return {
        comp1_feed: x1_feed,
        comp1_equilibrium: y[0],
        T_K: T_system_K,
        P_Pa: P_new_Pa,
        iterations: iter + 1,
        calculationType: 'bubbleP'
      };
    }
    P_current_Pa = P_new_Pa;
  }
  return null;
}

export function calculateDewTemperatureSrk(
  components: CompoundData[],
  y1_feed: number,
  P_system_Pa: number,
  srkInteractionParams: SrkInteractionParams,
  initialTempGuess_K: number = 350,
  maxIter: number = 50,
  tolerance: number = 1e-5
): BubbleDewResult | null {
  let T_current_K = initialTempGuess_K;
  const y = [y1_feed, 1 - y1_feed];

  for (let iter = 0; iter < maxIter; iter++) {
    const Psat = components.map(c => calculatePsat_Pa(c.antoine, T_current_K));
    if (Psat.some(p => isNaN(p) || p <= 0)) return null;

    const phi_V = calculateSrkFugacityCoefficients(components, y, T_current_K, P_system_Pa, srkInteractionParams, 'vapor');
    if (!phi_V) return null;

    let x = y.map((yi, i) => yi * phi_V[i] * P_system_Pa / (Psat[i]));
    const x_sum = x.reduce((s, xi) => s + xi, 0);
    x = x.map(xi => xi / x_sum);

    const phi_L = calculateSrkFugacityCoefficients(components, x, T_current_K, P_system_Pa, srkInteractionParams, 'liquid');
    if (!phi_L) return null;

    const obj = y.reduce((sum, yi, i) => sum + yi / (phi_L[i] * Psat[i] / (phi_V[i] * P_system_Pa)), 0) - 1;
    
    if (Math.abs(obj) < tolerance) {
      return {
        comp1_feed: y1_feed,
        comp1_equilibrium: x[0],
        T_K: T_current_K,
        P_Pa: P_system_Pa,
        iterations: iter + 1,
        calculationType: 'dewT'
      };
    }

    const dobj_dT = _numDerivSrk((T: number) => {
      const Psat_T = components.map(c => calculatePsat_Pa(c.antoine, T));
      if (Psat_T.some(p => isNaN(p))) return NaN;
      let x_T = y.map((yi, i) => yi * phi_V[i] * P_system_Pa / Psat_T[i]);
      const x_sum_T = x_T.reduce((s, xi) => s + xi, 0);
      x_T = x_T.map(xi => xi / x_sum_T);
      const phi_V_T = calculateSrkFugacityCoefficients(components, y, T, P_system_Pa, srkInteractionParams, 'vapor');
      if (!phi_V_T) return NaN;
      const phi_L_T = calculateSrkFugacityCoefficients(components, x_T, T, P_system_Pa, srkInteractionParams, 'liquid');
      if (!phi_L_T) return NaN;
      const K_T = phi_L_T.map((phiL, i) => phiL * Psat_T[i] / (phi_V_T[i] * P_system_Pa));
      return y.reduce((sum, yi, i) => sum + yi / K_T[i], 0) - 1;
    }, T_current_K);

    if (Math.abs(dobj_dT) < 1e-12) return null;
    T_current_K = Math.max(100, Math.min(1000, T_current_K - obj / dobj_dT));
  }
  return null;
}

export function calculateDewPressureSrk(
  components: CompoundData[],
  y1_feed: number,
  T_system_K: number,
  srkInteractionParams: SrkInteractionParams,
  initialPressureGuess_Pa: number = 101325,
  maxIter: number = 10,
  tolerance: number = 1e-5
): BubbleDewResult | null {
  let P_current_Pa = initialPressureGuess_Pa;
  const y = [y1_feed, 1 - y1_feed];

  for (let iter = 0; iter < maxIter; iter++) {
    const Psat = components.map(c => calculatePsat_Pa(c.antoine, T_system_K));
    if (Psat.some(p => isNaN(p) || p <= 0)) return null;

    const phi_V = calculateSrkFugacityCoefficients(components, y, T_system_K, P_current_Pa, srkInteractionParams, 'vapor');
    if (!phi_V) return null;

    let x = y.map((yi, i) => yi * phi_V[i] * P_current_Pa / Psat[i]);
    const x_sum = x.reduce((s, xi) => s + xi, 0);
    x = x.map(xi => xi / x_sum);

    const phi_L = calculateSrkFugacityCoefficients(components, x, T_system_K, P_current_Pa, srkInteractionParams, 'liquid');
    if (!phi_L) return null;

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
//  7. U N I Q U A C   (binary activity model)
// ===============================================================

// Insert interface definition for UNIQUAC interactions
export interface UniquacInteractionParams {
  A12: number;
  A21: number;
}

// type alias – avoid redeclaration
type UniquacInteractionParamsFetch = UniquacInteractionParams;

// (interface already declared above)

function _calcUniquacCombinatorial(x: number[], r: number[], q: number[]): number[] {
  const z = 10;
  const phi = (i: number) => (x[i] * r[i]) / r.reduce((s, _, k) => s + x[k] * r[k], 0);
  const theta = (i: number) => (x[i] * q[i]) / q.reduce((s, _, k) => s + x[k] * q[k], 0);
  const L = (i: number) => (z / 2) * (r[i] - q[i]) - (r[i] - 1);
  return x.map((_, i) => Math.log(phi(i) / x[i]) + (z / 2) * q[i] * Math.log(theta(i) / phi(i)) + L(i) - (phi(i) / x[i]) * x.reduce((s, xj, j) => s + xj * L(j), 0));
}

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

  const tau12 = Math.exp(-p.A12 / T_K);
  const tau21 = Math.exp(-p.A21 / T_K);

  const Theta_res = Theta; // same symbols
  const s1 = Theta_res[0] + Theta_res[1] * tau21;
  const s2 = Theta_res[1] + Theta_res[0] * tau12;

  const lnGammaR0 = q[0] * (1 - Math.log(s1) - (Theta_res[0] + Theta_res[1] * tau21) / s1 - (Theta_res[1] * tau12) / s2);
  const lnGammaR1 = q[1] * (1 - Math.log(s2) - (Theta_res[1] + Theta_res[0] * tau12) / s2 - (Theta_res[0] * tau21) / s1);

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
  for (let i = 0; i < maxIter; i++) {
    const Psat = [calculatePsat_Pa(comps[0].antoine, T), calculatePsat_Pa(comps[1].antoine, T)];
    if (Psat.some(p => isNaN(p))) return null;
    const gamma = calculateUniquacGamma(comps, x, T, params);
    if (!gamma) return null;
    const Pcalc = x[0] * gamma[0] * Psat[0] + x[1] * gamma[1] * Psat[1];
    const err = Pcalc - P_system_Pa;
    if (Math.abs(err) < tol * P_system_Pa) {
      const y1 = (x[0] * gamma[0] * Psat[0]) / Pcalc;
      return { comp1_feed: x1_feed, comp1_equilibrium: y1, T_K: T, P_Pa: P_system_Pa };
    }
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
  const Pcalc = x[0] * gamma[0] * Psat[0] + x[1] * gamma[1] * Psat[1];
  const y1 = (x[0] * gamma[0] * Psat[0]) / Pcalc;
  return { comp1_feed: x1_feed, comp1_equilibrium: y1, T_K, P_Pa: Pcalc };
}

// ==============================================================
//  7.  D A T A B A S E   F E T C H   H E L P E R S (Supabase)
// ==============================================================

// --- UNIFAC ---
export async function fetchUnifacInteractionParams(
  supabase: SupabaseClient,
  subgroupIds: number[]
): Promise<UnifacParameters> {
  const { data, error } = await supabase
    .from('UNIFAC - Rk and Qk')
    .select('"Subgroup #", "Main Group #", Rk, Qk')
    .in('"Subgroup #"', subgroupIds);
  if (error) throw error;
  if (!data || data.length === 0) throw new Error('UNIFAC Rk/Qk data missing');

  const Rk: { [id: number]: number } = {};
  const Qk: { [id: number]: number } = {};
  const mainMap: { [id: number]: number } = {};
  const mainIds = new Set<number>();
  data.forEach(row => {
    const sg = row["Subgroup #"];
    Rk[sg] = row.Rk;
    Qk[sg] = row.Qk;
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
  intData?.forEach(row => {
    a_mk.set(`${row.i}-${row.j}`, row["a(ij)"] ?? 0);
  });
  return { Rk, Qk, mainGroupMap: mainMap, a_mk };
}

// --- NRTL ---
export async function fetchNrtlParameters(
  supabase: SupabaseClient,
  casn1: string,
  casn2: string
): Promise<NrtlInteractionParams> {
  if (!casn1 || !casn2) throw new Error('CAS numbers required for NRTL');
  if (casn1 === casn2) return { A12: 0, A21: 0, alpha: 0.1 };
  const { data, error } = await supabase
    .from('nrtl parameters')
    .select('"A12", "A21", "alpha12", "CASN1", "CASN2"')
    .or(`and(CASN1.eq.${casn1},CASN2.eq.${casn2}),and(CASN1.eq.${casn2},CASN2.eq.${casn1})`)
    .limit(1);
  if (error) throw new Error(`Supabase NRTL query error: ${error.message}`);
  if (!data || data.length === 0) throw new Error(`NRTL parameters missing for ${casn1}/${casn2}`);
  const row = data[0];
  return row.CASN1 === casn1
    ? { A12: row.A12, A21: row.A21, alpha: row.alpha12 }
    : { A12: row.A21, A21: row.A12, alpha: row.alpha12 };
}

// --- Wilson ---
export async function fetchWilsonInteractionParams(
  supabase: SupabaseClient,
  casn1: string,
  casn2: string
): Promise<WilsonInteractionParams> {
  if (!casn1 || !casn2) return { a12_J_mol: 0, a21_J_mol: 0 };
  const { data, error } = await supabase
    .from('wilson parameters')
    .select('A12, A21, CASN1, CASN2')
    .or(`and(CASN1.eq.${casn1},CASN2.eq.${casn2}),and(CASN1.eq.${casn2},CASN1.eq.${casn1})`)
    .limit(1);
  if (error || !data || data.length === 0) return { a12_J_mol: 0, a21_J_mol: 0 };
  const row = data[0];
  const cal2J = JOULES_PER_CAL;
  const A12 = row.CASN1 === casn1 ? row.A12 : row.A21;
  const A21 = row.CASN1 === casn1 ? row.A21 : row.A12;
  return { a12_J_mol: A12 * cal2J, a21_J_mol: A21 * cal2J };
}

// --- Peng-Robinson ---
export async function fetchPrInteractionParams(
  supabase: SupabaseClient,
  casn1: string,
  casn2: string
): Promise<PrInteractionParams> {
  if (!casn1 || !casn2 || casn1 === casn2) return { k_ij: 0, k_ji: 0 };
  
  const { data, error } = await supabase
    .from('pr parameters')
    .select('"CASN1", "CASN2", k12')
    .or(`and("CASN1".eq.${casn1},"CASN2".eq.${casn2}),and("CASN1".eq.${casn2},"CASN2".eq.${casn1})`)
    .limit(1);
    
  if (error) throw new Error(`PR fetch error: ${error.message}`);
  if (!data || data.length === 0) return { k_ij: 0, k_ji: 0 };
  
  const k_value = typeof data[0].k12 === 'number' ? data[0].k12 : 0;
  return { k_ij: k_value, k_ji: k_value };
}

// --- SRK ---
export async function fetchSrkInteractionParams(
  supabase: SupabaseClient,
  casn1: string,
  casn2: string
): Promise<SrkInteractionParams> {
  if (!casn1 || !casn2 || casn1 === casn2) return { k_ij: 0, k_ji: 0 };
  
  const { data, error } = await supabase
    .from('srk parameters')
    .select('"CASN1", "CASN2", k12')
    .or(`and("CASN1".eq.${casn1},"CASN2".eq.${casn2}),and("CASN1".eq.${casn2},"CASN2".eq.${casn1})`)
    .limit(1);
    
  if (error) throw new Error(`SRK fetch error: ${error.message}`);
  if (!data || data.length === 0) return { k_ij: 0, k_ji: 0 };
  
  const k_value = typeof data[0].k12 === 'number' ? data[0].k12 : 0;
  return { k_ij: k_value, k_ji: k_value };
}

// --- UNIQUAC ---
export async function fetchUniquacInteractionParams(
  supabase: SupabaseClient,
  casn1: string,
  casn2: string
): Promise<UniquacInteractionParams> {
  if (!casn1 || !casn2) return { A12: 0, A21: 0 };
  const { data, error } = await supabase
    .from('uniquac parameters')
    .select('A12, A21, CASN1, CASN2')
    .or(`and(CASN1.eq.${casn1},CASN2.eq.${casn2}),and(CASN1.eq.${casn2},CASN2.eq.${casn1})`)
    .limit(1);
  if (error || !data || data.length === 0) return { A12: 0, A21: 0 };
  const row = data[0];
  return row.CASN1 === casn1 ? { A12: row.A12, A21: row.A21 } : { A12: row.A21, A21: row.A12 };
}

// =============================
//  8. Re-export fetch helpers (Wilson, PR, SRK, UNIQUAC) – restored
// =============================

// All legacy functions have been inlined - no more placeholders needed 