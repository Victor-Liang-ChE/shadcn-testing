// =======================================================================================
//  Vapor-Liquid Equilibrium (VLE) Calculation Utilities - MULTICOMPONENT VERSION
// =======================================================================================
//  This module provides functions to calculate VLE properties for multicomponent mixtures
//  using various thermodynamic models. It is an extension of the original binary-only file.
//
//  Supported Models:
//  • NRTL (Non-Random Two Liquid) - Activity coefficient model
//  • Wilson - Activity coefficient model
//  • UNIFAC - Group contribution method (inherently multicomponent)
//  • Peng-Robinson (PR) - Cubic equation of state with fugacity coefficients
//  • Soave-Redlich-Kwong (SRK) - Cubic equation of state with fugacity coefficients
//  • UNIQUAC - Activity coefficient model (inherently multicomponent)
// =======================================================================================

import type { CompoundData, AntoineParams } from './vle-types';
export const R_gas_const_J_molK = 8.31446261815324; // J·mol⁻¹·K⁻¹

/**
 * Calculates saturation pressure (Pa) from HYSYS extended Antoine.
 * ln(P[kPa]) = A + B/(T+C) + D·ln(T) + E·T^F
 * T in Kelvin. Returns pressure in Pa (kPa × 1000).
 */
export function calculatePsat_Pa(
  params: AntoineParams | null,
  T_K: number
): number {
  if (!params || T_K <= 0) return NaN;
  const { A, B, C, D, E, F } = params;
  const lnP_kPa = A + B / (T_K + C) + D * Math.log(T_K) + E * Math.pow(T_K, F);
  const P_Pa = Math.exp(lnP_kPa) * 1000; // kPa → Pa
  return isNaN(P_Pa) || P_Pa <= 0 ? NaN : P_Pa;
}
/**
 * Solves a cubic equation of state: aZ³ + bZ² + cZ + d = 0.
 * Returns only the positive real roots, which is essential for Z-factor calculations.
 */
export function solveCubicEOS(a: number, b: number, c: number, d: number): number[] | null {
  if (Math.abs(a) < 1e-12) return null;
  const p = b / a, q = c / a, r = d / a;
  const A = q - p * p / 3;
  const B = (2 * p * p * p) / 27 - (p * q) / 3 + r;
  let D = (B * B) / 4 + (A * A * A) / 27;
  if (Math.abs(D) < 1e-20) D = 0;
  let roots: number[] = [];
  if (D > 0) {
    const S = Math.cbrt(-B / 2 + Math.sqrt(D));
    const T = Math.cbrt(-B / 2 - Math.sqrt(D));
    roots = [S + T - p / 3];
  } else {
    const rho = 2 * Math.sqrt(-A / 3);
    const theta = Math.acos(Math.max(-1, Math.min(1, (-B / 2) / Math.sqrt(-(A * A * A) / 27)))) / 3;
    roots = [
      rho * Math.cos(theta) - p / 3,
      rho * Math.cos(theta + (2 * Math.PI) / 3) - p / 3,
      rho * Math.cos(theta + (4 * Math.PI) / 3) - p / 3,
    ];
  }
  const realPos = roots.filter(z => z > 1e-9).sort((x, y) => x - y);
  return realPos.length ? realPos : null;
}

// ==================================================================
//  2. N R T L (Multicomponent Activity-Coefficient Model)
// ==================================================================
export interface NrtlInteractionParams {
  Aij: number;  // K (temperature-dependent part: tau = Aij/T + Bij)
  Aji: number;  // K
  Bij: number;  // dimensionless
  Bji: number;  // dimensionless
  alpha: number; // non-randomness parameter
}
export type NrtlParameterMatrix = Map<string, NrtlInteractionParams>; // Key: "i-j"

/**
 * Calculates multicomponent NRTL activity coefficients (gamma).
 */
export function calculateNrtlGammaMulticomponent(
  components: (CompoundData & { phase?: 'Liquid' | 'Gas' })[], // Component metadata array (used for better error messages)
  x: number[], // Mole fractions array [x1, x2, ..., xN]
  T_K: number,
  params: NrtlParameterMatrix
): number[] | null {
  const n = x.length;

  // Quick sanity: ensure components array matches mole-fraction length
  if (n === 0 || n !== components.length) {
    // If mismatch or empty, return ideal gammas to avoid downstream crashes
    return new Array(n).fill(1.0);
  }

  // =======================================================
  // ▼▼▼ SAFEGUARD BLOCK: handle invalid / empty inputs ▼▼▼
  // If input is empty, contains NaN, or sums to (near) zero,
  // return ideal activity coefficients (gamma = 1) to
  // avoid exceptions/crashes in callers (ODE solvers etc.).
  const sumX = x.reduce((a, b) => a + b, 0);
  if (n === 0 || isNaN(sumX) || sumX < 1e-12) {
    return new Array(n).fill(1.0);
  }
  // ▲▲▲ END SAFEGUARD BLOCK ▲▲▲

  const tau = Array(n).fill(0).map(() => Array(n).fill(0));
  const G = Array(n).fill(0).map(() => Array(n).fill(0));
  // NRTL in HYSYS table uses Aij in cal/mol for this dataset
  const R_nrtl_cal = 1.9872; // cal·mol⁻¹·K⁻¹

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i === j) continue;
      // If either component is a gas, treat this pair as ideal by skipping parameter lookup.
      if (components[i]?.phase === 'Gas' || components[j]?.phase === 'Gas') {
        continue;
      }

      const pairParams = params.get(`${i}-${j}`) || params.get(`${j}-${i}`);
      if (!pairParams) {
        const comp1Name = components[i]?.name || `Component ${i}`;
        const comp2Name = components[j]?.name || `Component ${j}`;
        console.error(`❌ [NRTL CALC FAILED] Missing interaction parameters for pair: ${comp1Name} - ${comp2Name}.`);
        return null; // Signal failure
      }

      const p = (params.has(`${i}-${j}`)) 
          ? pairParams 
          : { Aij: pairParams.Aji, Aji: pairParams.Aij, Bij: pairParams.Bji, Bji: pairParams.Bij, alpha: pairParams.alpha };
      
      // Use cal-based R for NRTL interactions from HYSYS, with swap convention:
      // tau[i][j] uses Aji/Bji and tau[j][i] uses Aij/Bij to match HYSYS export format.
      tau[i][j] = p.Aji / (R_nrtl_cal * T_K) + p.Bji;
      tau[j][i] = p.Aij / (R_nrtl_cal * T_K) + p.Bij;
      G[i][j] = Math.exp(-p.alpha * tau[i][j]);
      G[j][i] = Math.exp(-p.alpha * tau[j][i]);
    }
  }

  const lnGamma = new Array(n).fill(0);

  for (let i = 0; i < n; i++) {
    let term1_num = 0;
    let term1_den = 0;
    for (let j = 0; j < n; j++) {
      term1_num += x[j] * tau[j][i] * G[j][i];
      term1_den += x[j] * G[j][i];
    }
    const term1 = term1_den !== 0 ? term1_num / term1_den : 0;

    let term2 = 0;
    for (let j = 0; j < n; j++) {
      const num_j = x[j] * G[i][j];
      let den_k = 0;
      for (let k = 0; k < n; k++) {
        den_k += x[k] * G[k][j];
      }
      if (den_k === 0) continue;

      let inner_sum_m = 0;
      for (let m = 0; m < n; m++) {
        inner_sum_m += x[m] * tau[m][j] * G[m][j];
      }
      term2 += (num_j / den_k) * (tau[i][j] - (inner_sum_m / den_k));
    }
    lnGamma[i] = term1 + term2;
  }

  const gammas = lnGamma.map(val => Math.exp(val));

  // --- START: ADD DETAILED LOGGING ---
  console.log('--- NRTL Gamma Calculation Step ---');
  components.forEach((comp, i) => {
      // Only log for components with significant presence to keep console clean
      if (x[i] > 1e-6) {
          console.log(`  [idx ${i}] ${comp.name} (x = ${x[i].toFixed(4)}): ln(γ) = ${lnGamma[i].toFixed(4)}, γ = ${gammas[i].toFixed(4)}`);
      }
  });
  console.log('---------------------------------');
  // --- END: ADD DETAILED LOGGING ---

  return gammas;
}


// ==================================================================
//  3. W I L S O N (Multicomponent Activity Model)
// ==================================================================
export interface WilsonInteractionParams {
  Aij: number;  // K (Lambda = (Vj/Vi)·exp(-(Aij/T + Bij)))
  Aji: number;  // K
  Bij: number;  // dimensionless
  Bji: number;  // dimensionless
}
export type WilsonParameterMatrix = Map<string, WilsonInteractionParams>; // Key: "i-j"

/**
 * Calculates multicomponent Wilson activity coefficients (gamma).
 */
export function calculateWilsonGammaMulticomponent(
  comps: (CompoundData & { phase?: 'Liquid' | 'Gas' })[],
  x: number[],
  T_K: number,
  params: WilsonParameterMatrix
): number[] | null {
  const n = comps.length;
  if (x.length !== n || n === 0) return null;

  const V_L = comps.map(c => c.wilsonParams?.V_L_m3mol);
  if (V_L.some(v => v === undefined)) return null;
  // HYSYS Wilson Aij are in cal/mol. Use R_cal and swap convention.
  const R_wilson_cal = 1.9872; // cal·mol⁻¹·K⁻¹

  const Lambda = Array(n).fill(0).map(() => Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i === j) {
        Lambda[i][j] = 1;
        continue;
      }
      // If either component is a gas, treat as ideal (Lambda = 1) and skip parameter lookup.
      if (comps[i]?.phase === 'Gas' || comps[j]?.phase === 'Gas') {
        Lambda[i][j] = 1;
        continue;
      }
      const pairParams = params.get(`${i}-${j}`) || params.get(`${j}-${i}`);
      if (!pairParams) {
        const comp1Name = comps[i]?.name || `Component ${i}`;
        const comp2Name = comps[j]?.name || `Component ${j}`;
        console.error(`❌ [Wilson CALC FAILED] Missing interaction parameters for pair: ${comp1Name} - ${comp2Name}.`);
        return null;
      }

      const p = params.has(`${i}-${j}`)
        ? pairParams
        : { Aij: pairParams.Aji, Aji: pairParams.Aij, Bij: pairParams.Bji, Bji: pairParams.Bij };
      // Swap convention: Lambda[i][j] uses Aji/Bji to match HYSYS export format.
      Lambda[i][j] = (V_L[j]! / V_L[i]!) * Math.exp(-(p.Aji / (R_wilson_cal * T_K) + p.Bji));
    }
  }

  const lnGamma = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    let term1_sum = 0;
    for (let j = 0; j < n; j++) {
      term1_sum += x[j] * Lambda[i][j];
    }
    const term1 = -Math.log(term1_sum);

    let term2 = 0;
    for (let k = 0; k < n; k++) {
      let term2_num = 0;
      let term2_den = 0;
      for (let j = 0; j < n; j++) {
        term2_den += x[j] * Lambda[k][j];
      }
      if (term2_den === 0) continue;
      term2_num = x[k] * Lambda[k][i];
      term2 += term2_num / term2_den;
    }
    lnGamma[i] = term1 + (1 - term2);
  }

  return lnGamma.map(val => Math.exp(val));
}


// ==================================================================
//  4. P E N G – R O B I N S O N   E O S (Multicomponent)
// ==================================================================
export interface PrSrkInteractionParams { k_ij: number; }
export type PrSrkParameterMatrix = Map<string, PrSrkInteractionParams>; // Key: "i-j"

/**
 * Calculates multicomponent Peng-Robinson fugacity coefficients (phi).
 */
export function calculatePrFugacityCoefficientsMulticomponent(
  components: CompoundData[],
  x_or_y: number[], // Mole fractions [x1, x2, ..., xN]
  T_K: number,
  P_Pa: number,
  interactionParams: PrSrkParameterMatrix,
  phase: 'liquid' | 'vapor'
): number[] | null {
  const n = components.length;
  if (x_or_y.length !== n || n === 0) return null;

  const pureParams = components.map(c => {
    if (!c.prParams) return null;
    const { Tc_K, Pc_Pa, omega } = c.prParams;
    const Tr = T_K / Tc_K;
    const kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    const alpha = (1 + kappa * (1 - Math.sqrt(Tr))) ** 2;
    const ac = 0.45723553 * (R_gas_const_J_molK * Tc_K) ** 2 / Pc_Pa * alpha;
    const b = 0.07779607 * R_gas_const_J_molK * Tc_K / Pc_Pa;
    return { ac, b };
  });
  if (pureParams.some(p => p === null)) return null;

  // Multicomponent mixing rules
  let b_mix = 0;
  for (let i = 0; i < n; i++) {
    b_mix += x_or_y[i] * pureParams[i]!.b;
  }

  let a_mix = 0;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const k_ij = interactionParams.get(`${i}-${j}`)?.k_ij || interactionParams.get(`${j}-${i}`)?.k_ij || 0;
      const a_ij = Math.sqrt(pureParams[i]!.ac * pureParams[j]!.ac) * (1 - k_ij);
      a_mix += x_or_y[i] * x_or_y[j] * a_ij;
    }
  }

  const A = a_mix * P_Pa / (R_gas_const_J_molK * T_K) ** 2;
  const B = b_mix * P_Pa / (R_gas_const_J_molK * T_K);

  const roots = solveCubicEOS(1, -(1 - B), A - 3 * B ** 2 - 2 * B, -(A * B - B ** 2 - B ** 3));
  if (!roots || roots.length === 0) return null;
  const Z = phase === 'liquid' ? Math.min(...roots) : Math.max(...roots);
  if (Z <= B) return null;

  const termCommon = (A / (2 * Math.sqrt(2) * B)) * Math.log((Z + (1 + Math.sqrt(2)) * B) / (Z + (1 - Math.sqrt(2)) * B));

  const ln_phi = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const b_i = pureParams[i]!.b;
    let sum_a_ij = 0;
    for (let j = 0; j < n; j++) {
        const k_ij = interactionParams.get(`${i}-${j}`)?.k_ij || interactionParams.get(`${j}-${i}`)?.k_ij || 0;
        sum_a_ij += x_or_y[j] * Math.sqrt(pureParams[i]!.ac * pureParams[j]!.ac) * (1 - k_ij);
    }
    const term1 = (b_i / b_mix) * (Z - 1) - Math.log(Z - B);
    const term2 = termCommon * ((2 * sum_a_ij / a_mix) - (b_i / b_mix));
    ln_phi[i] = term1 - term2;
  }

  return ln_phi.map(val => Math.exp(val));
}


// ===============================================================
//  5. S R K   E O S (Multicomponent)
// ===============================================================
/**
 * Calculates multicomponent Soave-Redlich-Kwong fugacity coefficients (phi).
 */
export function calculateSrkFugacityCoefficientsMulticomponent(
  components: CompoundData[],
  x_or_y: number[],
  T_K: number,
  P_Pa: number,
  interactionParams: PrSrkParameterMatrix,
  phase: 'liquid' | 'vapor'
): number[] | null {
  const n = components.length;
  if (x_or_y.length !== n || n === 0) return null;

  const pureParams = components.map(c => {
    if (!c.srkParams) return null;
    const { Tc_K, Pc_Pa, omega } = c.srkParams;
    const Tr = T_K / Tc_K;
    const m = 0.480 + 1.574 * omega - 0.176 * omega * omega;
    const alpha = Math.pow(1 + m * (1 - Math.sqrt(Tr)), 2);
    const ac = 0.42748 * (R_gas_const_J_molK * Tc_K) ** 2 / Pc_Pa * alpha;
    const b = 0.08664 * R_gas_const_J_molK * Tc_K / Pc_Pa;
    return { ac, b };
  });
  if (pureParams.some(p => p === null)) return null;

  let b_mix = 0;
  for (let i = 0; i < n; i++) {
    b_mix += x_or_y[i] * pureParams[i]!.b;
  }

  let a_mix = 0;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const k_ij = interactionParams.get(`${i}-${j}`)?.k_ij || interactionParams.get(`${j}-${i}`)?.k_ij || 0;
      const a_ij = Math.sqrt(pureParams[i]!.ac * pureParams[j]!.ac) * (1 - k_ij);
      a_mix += x_or_y[i] * x_or_y[j] * a_ij;
    }
  }

  const A = a_mix * P_Pa / (R_gas_const_J_molK * T_K) ** 2;
  const B = b_mix * P_Pa / (R_gas_const_J_molK * T_K);
  
  const roots = solveCubicEOS(1, -1, A - B - B * B, -A * B);
  if (!roots || roots.length === 0) return null;
  const Z = phase === 'liquid' ? Math.min(...roots) : Math.max(...roots);
  if (Z <= B) return null;

  const termCommon = (A / B) * Math.log(1 + B / Z);

  const ln_phi = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const b_i = pureParams[i]!.b;
    let sum_a_ij = 0;
    for (let j = 0; j < n; j++) {
        const k_ij = interactionParams.get(`${i}-${j}`)?.k_ij || interactionParams.get(`${j}-${i}`)?.k_ij || 0;
        sum_a_ij += x_or_y[j] * Math.sqrt(pureParams[i]!.ac * pureParams[j]!.ac) * (1 - k_ij);
    }
    const term1 = (b_i / b_mix) * (Z - 1) - Math.log(Z - B);
    const term2 = termCommon * ((2 * sum_a_ij / a_mix) - (b_i / b_mix));
    ln_phi[i] = term1 - term2;
  }

  return ln_phi.map(val => Math.exp(val));
}


// ===============================================================
//  6. U N I Q U A C (Multicomponent - from original file)
// ===============================================================
export interface UniquacInteractionParams {
  Aij: number;  // K (tau = exp(-(Aij/T + Bij)))
  Aji: number;  // K
  Bij: number;  // dimensionless
  Bji: number;  // dimensionless
}
export type UniquacParameterMatrix = Map<string, UniquacInteractionParams>; // Key "i-j"

/**
 * Calculates multicomponent UNIQUAC activity coefficients (gamma).
 */
export function calculateUniquacGammaMulticomponent(
  comps: (CompoundData & { phase?: 'Liquid' | 'Gas' })[],
  x: number[],
  T_K: number,
  params: UniquacParameterMatrix
): number[] | null {
  const n = comps.length;
  if (x.length !== n || n === 0) return null;
  if (comps.some(c => !c.uniquacParams)) return null;

  const r = comps.map(c => c.uniquacParams!.r);
  const q = comps.map(c => c.uniquacParams!.q);
  const z = 10;

  // Combinatorial part
  const sum_xr = x.reduce((s, xi, i) => s + xi * r[i], 0);
  const sum_xq = x.reduce((s, xi, i) => s + xi * q[i], 0);
  if (sum_xr === 0 || sum_xq === 0) return null;

  const phi = x.map((xi, i) => (xi * r[i]) / sum_xr);
  const theta = x.map((xi, i) => (xi * q[i]) / sum_xq);
  const L = r.map((ri, i) => (z / 2) * (ri - q[i]) - (ri - 1));
  const sum_xL = x.reduce((s, xi, i) => s + xi * L[i], 0);

  const lnGammaC = x.map((_, i) =>
    Math.log(phi[i] / x[i]) + (z / 2) * q[i] * Math.log(theta[i] / phi[i]) + L[i] - (phi[i] / x[i]) * sum_xL
  );

  // Residual part
  const R_si = R_gas_const_J_molK; // J·mol⁻¹·K⁻¹ - used in UNIQUAC tau calculation
  const tau = Array(n).fill(0).map(() => Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i === j) { 
        tau[i][j] = 1; 
        continue; 
      }
      // If either component is a gas, treat as ideal (tau = 1) and skip parameter lookup.
      if (comps[i]?.phase === 'Gas' || comps[j]?.phase === 'Gas') {
        tau[i][j] = 1;
        continue;
      }
    const pairParams = params.get(`${i}-${j}`) || params.get(`${j}-${i}`);
    if (!pairParams) {
      const comp1Name = comps[i]?.name || `Component ${i}`;
      const comp2Name = comps[j]?.name || `Component ${j}`;
      console.error(`❌ [UNIQUAC CALC FAILED] Missing interaction parameters for pair: ${comp1Name} - ${comp2Name}.`);
      return null;
    }
    const p = params.has(`${i}-${j}`)
      ? pairParams
      : { Aij: pairParams.Aji, Aji: pairParams.Aij, Bij: pairParams.Bji, Bji: pairParams.Bij };
        // Swap convention: tau[i][j] uses Aji/Bji to match HYSYS export format.
        tau[i][j] = Math.exp(-(p.Aji / (R_si * T_K) + p.Bji));
    }
  }

  const lnGammaR = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
      let term1_sum = 0;
      for (let j = 0; j < n; j++) {
          term1_sum += theta[j] * tau[j][i];
      }
      const term1 = -Math.log(term1_sum);

      let term2 = 0;
      for (let j = 0; j < n; j++) {
          let den = 0;
          for (let k = 0; k < n; k++) {
              den += theta[k] * tau[k][j];
          }
          if (den === 0) continue;
          term2 += (theta[j] * tau[i][j]) / den;
      }
      lnGammaR[i] = q[i] * (1 - term1 - term2);
  }

  return lnGammaC.map((lnC, i) => Math.exp(lnC + lnGammaR[i]));
}

// NOTE: UNIFAC implementation from the original file is already multicomponent-ready and can be used directly.

// =====================================================================
//  UNIFAC (for direct use and parameter estimation)
// =====================================================================

export interface UnifacParameters {
  Rk: { [subgroupId: number]: number };
  Qk: { [subgroupId: number]: number };
  mainGroupMap: { [subgroupId: number]: number };
  a_mk: Map<string, number>; // key `${mainGroup_m}-${mainGroup_k}`
}

/**
 * Calculates UNIFAC activity coefficients. This is needed for the estimation helpers.
 * This function is adapted from your original binary file and supports multicomponent mixtures.
 */
export function calculateUnifacGamma(
  componentData: CompoundData[],
  x: number[],
  T_kelvin: number,
  params: UnifacParameters
): number[] | null {
  if (componentData.length === 0 || x.length !== componentData.length) return null;
  const Z = 10.0;

  for (const c of componentData) {
    if (c.r_i === undefined || c.q_i === undefined) {
      if (!c.unifacGroups) {
        console.error(`[UNIFAC Calc Error] Component ${c.name} is missing UNIFAC group data.`);
        return null;
      }
      let rSum = 0, qSum = 0;
      for (const [subIdStr, count] of Object.entries(c.unifacGroups)) {
        const sgId = parseInt(subIdStr, 10);
        if (params.Rk[sgId] == null || params.Qk[sgId] == null) {
          // This is the new, detailed error log
          console.error(`[UNIFAC Calc Error] Missing Rk/Qk parameters for Subgroup ID #${sgId}, which is required by component ${c.name}. Please check your "UNIFAC - Rk and Qk" table.`);
          return null; // This is the failure point
        }
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

  const subgroupIds = Array.from(new Set(componentData.flatMap(c => Object.keys(c.unifacGroups || {}).map(Number))));
  if (subgroupIds.length === 0) return ln_gamma_C.map(Math.exp);

  const v_ki = componentData.map(c => subgroupIds.map(id => c.unifacGroups?.[id] || 0));
  
  const sum_x_vk_total = subgroupIds.reduce((sum, _, kIdx) => sum + componentData.reduce((s, _, i) => s + x[i] * v_ki[i][kIdx], 0), 0);
  if (sum_x_vk_total === 0) return ln_gamma_C.map(Math.exp);

  const X_m = subgroupIds.map((_, kIdx) => componentData.reduce((s, _, i) => s + x[i] * v_ki[i][kIdx], 0) / sum_x_vk_total);
  const sum_XQ = subgroupIds.reduce((s, id, idx) => s + X_m[idx] * params.Qk[id], 0);
  if (sum_XQ === 0) return ln_gamma_C.map(Math.exp);
  
  const theta_m = subgroupIds.map((id, idx) => (X_m[idx] * params.Qk[id]) / sum_XQ);
  const psi = subgroupIds.map((id_m, mIdx) => subgroupIds.map((id_k, kIdx) => {
      const main_group_m = params.mainGroupMap[id_m];
      const main_group_k = params.mainGroupMap[id_k];
      if (main_group_m === undefined || main_group_k === undefined) return 1.0;
      const key = `${main_group_m}-${main_group_k}`;
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

// Helper to fetch all UNIFAC group parameters (Rk, Qk, a_mk)
export async function fetchUnifacInteractionParams(
  supabase: SupabaseClient,
  subgroupIds: number[]
): Promise<UnifacParameters | null> {
  return _getUnifacParams(supabase);
}

// --- DATABASE FETCHING & ESTIMATION ---

import type { SupabaseClient } from '@supabase/supabase-js';
import { fetchCompoundDataFromHysys } from './antoine-utils';

let unifacParamsCache: UnifacParameters | null = null;
async function _getUnifacParams(supabase: SupabaseClient): Promise<UnifacParameters | null> {
    if (unifacParamsCache) return unifacParamsCache;

    const { data: subGroupsData, error: subGroupsError } = await supabase
        .from('UNIFAC - Rk and Qk')
        .select('id, "Main Group #", "Subgroup #", "Rk", "Qk"');
    if (subGroupsError) {
        console.error("Supabase error fetching from 'UNIFAC - Rk and Qk':", subGroupsError);
        return null;
    }

    const { data: interactionData, error: interactionError } = await supabase
        .from('UNIFAC - a(ij)')
        .select('i, j, "a(ij)"');
    if (interactionError) {
        console.error("Supabase error fetching from 'UNIFAC - a(ij)':", interactionError);
        return null;
    }

    const Rk = Object.fromEntries(subGroupsData.map((row: any) => [row['Subgroup #'], row.Rk]));
    const Qk = Object.fromEntries(subGroupsData.map((row: any) => [row['Subgroup #'], row.Qk]));
    const mainGroupMap = Object.fromEntries(subGroupsData.map((row: any) => [row['Subgroup #'], row['Main Group #']]));
    const a_mk = new Map(interactionData.map((row: any) => [`${row.i}-${row.j}`, row['a(ij)']]));

    unifacParamsCache = { Rk, Qk, mainGroupMap, a_mk: a_mk as Map<string, number> };
    return unifacParamsCache;
}

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
    const unifacParams = await _getUnifacParams(supabase);
    const comp1 = await _buildUnifacComponent(supabase, name1);
    const comp2 = await _buildUnifacComponent(supabase, name2);

    if (!unifacParams || !comp1 || !comp2) return null;

    const gammas_1_in_2 = calculateUnifacGamma([comp1, comp2], [1e-9, 1.0 - 1e-9], T_ref_K, unifacParams);
    const gammas_2_in_1 = calculateUnifacGamma([comp1, comp2], [1.0 - 1e-9, 1e-9], T_ref_K, unifacParams);

    if (!gammas_1_in_2 || !gammas_2_in_1) return null;

    const lnGamma1_inf = Math.log(gammas_2_in_1[0]);
    const lnGamma2_inf = Math.log(gammas_1_in_2[1]);

    // Heuristic: store as Bij (temperature-independent), Aij = 0
    return { Aij: 0, Aji: 0, Bij: lnGamma2_inf, Bji: lnGamma1_inf, alpha: 0.3 };
}

// --- NRTL ---
export async function fetchNrtlParameters(
  supabase: SupabaseClient,
  name1: string,
  name2: string,
): Promise<NrtlInteractionParams | null> {
  if (!name1 || !name2 || name1 === name2) return { Aij: 0, Aji: 0, Bij: 0, Bji: 0, alpha: 0.3 };

  const { data, error } = await supabase
    .from('HYSYS NRTL')
    .select('*')
    .or(`and(Component_i.ilike.${name1},Component_j.ilike.${name2}),and(Component_i.ilike.${name2},Component_j.ilike.${name1})`)
    .limit(1);

  if (error) console.warn(`Supabase HYSYS NRTL query error: ${error.message}`);

  if (data && data.length > 0) {
    const row = data[0];
    const isForward = row.Component_i.toLowerCase() === name1.toLowerCase();
    return isForward
      ? { Aij: row.Aij ?? 0, Aji: row.Aji ?? 0, Bij: row.Bij ?? 0, Bji: row.Bji ?? 0, alpha: row.Cij_Alpha ?? 0.3 }
      : { Aij: row.Aji ?? 0, Aji: row.Aij ?? 0, Bij: row.Bji ?? 0, Bji: row.Bij ?? 0, alpha: row.Cji_Alpha ?? row.Cij_Alpha ?? 0.3 };
  }
  
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
): Promise<WilsonInteractionParams | null> {
  if (!name1 || !name2) return { Aij: 0, Aji: 0, Bij: 0, Bji: 0 };

  const { data, error } = await supabase
    .from('HYSYS WILSON')
    .select('*')
    .or(`and(Component_i.ilike.${name1},Component_j.ilike.${name2}),and(Component_i.ilike.${name2},Component_j.ilike.${name1})`)
    .limit(1);

  if (error || !data || data.length === 0) return { Aij: 0, Aji: 0, Bij: 0, Bji: 0 };
  const row = data[0];
  const isForward = row.Component_i.toLowerCase() === name1.toLowerCase();
  return isForward
    ? { Aij: row.Aij ?? 0, Aji: row.Aji ?? 0, Bij: row.Bij ?? 0, Bji: row.Bji ?? 0 }
    : { Aij: row.Aji ?? 0, Aji: row.Aij ?? 0, Bij: row.Bji ?? 0, Bji: row.Bij ?? 0 };
}

// --- Peng-Robinson ---
export async function fetchPrInteractionParams(
  supabase: SupabaseClient,
  name1: string,
  name2: string
): Promise<PrSrkInteractionParams | null> {
  if (!name1 || !name2 || name1 === name2) return { k_ij: 0 };
  
  const { data, error } = await supabase
    .from('HYSYS PR')
    .select('*')
    .or(`and(Component_i.ilike.${name1},Component_j.ilike.${name2}),and(Component_i.ilike.${name2},Component_j.ilike.${name1})`)
    .limit(1);
    
  if (error) console.warn(`HYSYS PR fetch error: ${error.message}`);
  if (data && data.length > 0 && typeof data[0].Kij === 'number') {
    return { k_ij: data[0].Kij };
  }
  return { k_ij: 0 };
}

// --- SRK ---
export async function fetchSrkInteractionParams(
  supabase: SupabaseClient,
  name1: string,
  name2: string
): Promise<PrSrkInteractionParams | null> {
  if (!name1 || !name2 || name1 === name2) return { k_ij: 0 };
  
  const { data, error } = await supabase
    .from('HYSYS SRK')
    .select('*')
    .or(`and(Component_i.ilike.${name1},Component_j.ilike.${name2}),and(Component_i.ilike.${name2},Component_j.ilike.${name1})`)
    .limit(1);
    
  if (error) console.warn(`HYSYS SRK fetch error: ${error.message}`);
  if (data && data.length > 0 && typeof data[0].Kij === 'number') {
    return { k_ij: data[0].Kij };
  }
  return { k_ij: 0 };
}

// --- UNIQUAC ---
export async function fetchUniquacInteractionParams(
  supabase: SupabaseClient,
  name1: string,
  name2: string
): Promise<UniquacInteractionParams | null> {
  if (!name1 || !name2) return { Aij: 0, Aji: 0, Bij: 0, Bji: 0 };

  const { data, error } = await supabase
    .from('HYSYS UNIQUAC')
    .select('*')
    .or(`and(Component_i.ilike.${name1},Component_j.ilike.${name2}),and(Component_i.ilike.${name2},Component_j.ilike.${name1})`)
    .limit(1);

  if (!error && data && data.length > 0) {
    const row = data[0];
    const isForward = row.Component_i.toLowerCase() === name1.toLowerCase();
    return isForward
      ? { Aij: row.Aij ?? 0, Aji: row.Aji ?? 0, Bij: row.Bij ?? 0, Bji: row.Bji ?? 0 }
      : { Aij: row.Aji ?? 0, Aji: row.Aij ?? 0, Bij: row.Bji ?? 0, Bji: row.Bij ?? 0 };
  }
  if (error) console.warn(`HYSYS UNIQUAC fetch error: ${error.message}`);
  return { Aij: 0, Aji: 0, Bij: 0, Bji: 0 };
}