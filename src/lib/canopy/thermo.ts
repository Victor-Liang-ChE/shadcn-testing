// ─── Canopy Thermodynamic Engine ───────────────────────────────────────────────
// Provides flash calculations, enthalpy computation, and property evaluation
// for the Canopy process simulator.
//
// Phase 1 scope: Ideal (Raoult's Law) fluid package with DIPPR correlations.
// ──────────────────────────────────────────────────────────────────────────────

import type { CanopyCompound, Phase } from './types';

const R = 8.314462618; // J/(mol·K)

// ────────────────────────────────────────────────────────────────
// 1. Vapor Pressure — Extended Antoine (PLXANT / Aspen)
// ────────────────────────────────────────────────────────────────

/** ln(P* [Pa]) = C1 + C2/(T+C3) + C4·T + C5·ln(T) + C6·T^C7 */
export function Psat_Pa(comp: CanopyCompound, T_K: number): number {
  const a = comp.antoine;
  if (!a || T_K <= 0) return NaN;
  const lnP = a.C1 + a.C2 / (T_K + a.C3) + a.C4 * T_K
    + a.C5 * Math.log(T_K) + a.C6 * Math.pow(T_K, a.C7);
  const P = Math.exp(lnP);
  return isFinite(P) && P > 0 ? P : NaN;
}

// ────────────────────────────────────────────────────────────────
// 2. DIPPR Ideal Gas Heat Capacity  (Eq. 16)
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR Eq. 16:  Cp_ig = A + exp(B/T + C + D·T + E·T²)
 * Returns J/(kmol·K). Divide by 1000 for J/(mol·K).
 */
export function CpIG_JkmolK(comp: CanopyCompound, T_K: number): number {
  const c = comp.dipprCoeffs['IdealGasHeatCapacityCp'];
  if (!c) return NaN;
  const { A, B, C, D, E } = c;
  const E_ = E ?? 0;
  return A + Math.exp(B / T_K + C + D * T_K + E_ * T_K * T_K);
}

/** Cp in J/(mol·K) */
export function CpIG_JmolK(comp: CanopyCompound, T_K: number): number {
  return CpIG_JkmolK(comp, T_K) / 1000;
}

// ────────────────────────────────────────────────────────────────
// 3. DIPPR Heat of Vaporization
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR 106: ΔHvap = A · (1 - Tr)^(B + C·Tr + D·Tr² + E·Tr³)
 * where Tr = T/Tc.  Returns J/kmol.  Divide by 1000 for J/mol.
 */
export function Hvap_Jkmol(comp: CanopyCompound, T_K: number): number {
  const c = comp.dipprCoeffs['HeatOfVaporization'];
  if (!c) return NaN;
  const Tr = T_K / comp.Tc_K;
  if (Tr >= 1) return 0; // Above critical — no vaporization
  const { A, B, C, D, E } = c;
  const E_ = E ?? 0;
  return A * Math.pow(1 - Tr, B + C * Tr + D * Tr * Tr + E_ * Tr * Tr * Tr);
}

export function Hvap_Jmol(comp: CanopyCompound, T_K: number): number {
  return Hvap_Jkmol(comp, T_K) / 1000;
}

// ────────────────────────────────────────────────────────────────
// 4. DIPPR Liquid Heat Capacity
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR 114:  Cp_L/R = A²/τ + B - 2Aτ(C+3Dτ²+C(ε²-1)+4Dε(ε²-1)τ)
 * However many compounds use simpler forms (polynomial, etc.).
 * The database stores them as equation code + coeffs A-E.
 * We'll use the raw correlation form from compound_properties:
 *   Cp_L = A + exp(B/T + C + D·T + E·T²)   [J/kmol/K]  (DIPPR code 114 alt form)
 *
 * Actually in this DB, LiquidHeatCapacityCp uses the same template as vapor pressure:
 *   value = A + exp(B/T + C + D·T + E·T²)
 */
export function CpL_JkmolK(comp: CanopyCompound, T_K: number): number {
  const c = comp.dipprCoeffs['LiquidHeatCapacityCp'];
  if (!c) {
    // Fallback: use ideal gas Cp (rough approximation)
    return CpIG_JkmolK(comp, T_K);
  }
  const { A, B, C, D, E } = c;
  const E_ = E ?? 0;
  return A + Math.exp(B / T_K + C + D * T_K + E_ * T_K * T_K);
}

export function CpL_JmolK(comp: CanopyCompound, T_K: number): number {
  return CpL_JkmolK(comp, T_K) / 1000;
}

// ────────────────────────────────────────────────────────────────
// 5. Enthalpy Calculations (Reference: ideal gas at 298.15 K)
// ────────────────────────────────────────────────────────────────

const T_REF = 298.15; // K

/** Integrate Cp_ig from T_ref to T using 50-point Gauss-Legendre (simplified: Simpson's rule) */
function integrateCpIG(comp: CanopyCompound, T_K: number): number {
  if (Math.abs(T_K - T_REF) < 0.01) return 0;
  const n = 100;
  const h = (T_K - T_REF) / n;
  let sum = CpIG_JmolK(comp, T_REF) + CpIG_JmolK(comp, T_K);
  for (let i = 1; i < n; i++) {
    const T = T_REF + i * h;
    sum += (i % 2 === 0 ? 2 : 4) * CpIG_JmolK(comp, T);
  }
  return sum * h / 3;
}

/**
 * Ideal-gas molar enthalpy relative to reference state.
 * H_ig(T) = ∫[T_ref → T] Cp_ig dT   (J/mol)
 */
export function enthalpyIG(comp: CanopyCompound, T_K: number): number {
  return integrateCpIG(comp, T_K);
}

/**
 * Liquid molar enthalpy relative to ideal gas reference.
 * H_L(T) = H_ig(T) - ΔHvap(T)   (J/mol)
 *
 * This is a simplification. For ideal mixtures it's adequate.
 */
export function enthalpyLiquid(comp: CanopyCompound, T_K: number): number {
  return enthalpyIG(comp, T_K) - Hvap_Jmol(comp, T_K);
}

/**
 * Stream molar enthalpy for a mixture.
 * H = Σ(yi·Hi_ig) for vapor fraction + Σ(xi·Hi_L) for liquid fraction
 */
export function streamEnthalpy(
  compounds: CanopyCompound[],
  T_K: number,
  vaporFraction: number,
  x_liquid: number[],
  y_vapor: number[],
): number {
  let H = 0;
  for (let i = 0; i < compounds.length; i++) {
    const Hv = enthalpyIG(compounds[i], T_K);
    const Hl = enthalpyLiquid(compounds[i], T_K);
    H += vaporFraction * y_vapor[i] * Hv + (1 - vaporFraction) * x_liquid[i] * Hl;
  }
  return H;
}

// ────────────────────────────────────────────────────────────────
// 6. Flash Calculation — Isothermal PT Flash (Rachford-Rice)
// ────────────────────────────────────────────────────────────────

export interface FlashResult {
  converged: boolean;
  phase: Phase;
  vaporFraction: number;   // β — moles vapor / total moles
  x: number[];             // Liquid mole fractions
  y: number[];             // Vapor mole fractions
  K: number[];             // K-values (yi/xi)
  T_K: number;
  P_Pa: number;
  H_Jpmol: number;         // Mixture molar enthalpy
}

/**
 * K-values from Raoult's Law: Ki = Psat_i(T) / P
 * For Ideal fluid package.
 */
export function KvaluesIdeal(compounds: CanopyCompound[], T_K: number, P_Pa: number): number[] {
  return compounds.map(c => {
    const ps = Psat_Pa(c, T_K);
    return isFinite(ps) ? ps / P_Pa : 1;
  });
}

/**
 * Rachford-Rice equation: Σ zi(Ki-1) / (1 + β(Ki-1)) = 0
 * Solves for vapor fraction β.
 */
function rachfordRice(z: number[], K: number[]): number {
  const n = z.length;

  // Check: all K > 1 → all vapor
  if (K.every(k => k > 1)) {
    // Verify with bubble point check
    const sumZK = z.reduce((s, zi, i) => s + zi * K[i], 0);
    if (sumZK > 1) return 1; // Superheated vapor
  }

  // Check: all K < 1 → all liquid
  if (K.every(k => k < 1)) {
    const sumZoverK = z.reduce((s, zi, i) => s + zi / K[i], 0);
    if (sumZoverK > 1) return 0; // Subcooled liquid
  }

  // Bounds for β from Whitson & Michelsen
  let betaMin = 0, betaMax = 1;
  for (let i = 0; i < n; i++) {
    if (K[i] > 1) betaMin = Math.max(betaMin, (K[i] * z[i] - 1) / (K[i] - 1));
    if (K[i] < 1) betaMax = Math.min(betaMax, (1 - z[i]) / (1 - K[i]));
  }
  betaMin = Math.max(0, betaMin);
  betaMax = Math.min(1, betaMax);
  if (betaMin > betaMax) return z.reduce((s, zi, i) => s + zi * K[i], 0) > 1 ? 1 : 0;

  // Newton-Raphson on RR equation
  let beta = (betaMin + betaMax) / 2;
  for (let iter = 0; iter < 100; iter++) {
    let f = 0, df = 0;
    for (let i = 0; i < n; i++) {
      const km1 = K[i] - 1;
      const denom = 1 + beta * km1;
      f += z[i] * km1 / denom;
      df -= z[i] * km1 * km1 / (denom * denom);
    }
    if (Math.abs(f) < 1e-12) break;
    let step = -f / df;
    // Ensure beta stays in bounds
    let betaNew = beta + step;
    if (betaNew < betaMin) betaNew = (beta + betaMin) / 2;
    if (betaNew > betaMax) betaNew = (beta + betaMax) / 2;
    beta = betaNew;
  }

  return Math.max(0, Math.min(1, beta));
}

/**
 * Isothermal PT flash using Raoult's Law (Ideal fluid package).
 */
export function flashPT_Ideal(
  compounds: CanopyCompound[],
  z: number[],
  T_K: number,
  P_Pa: number,
): FlashResult {
  const n = compounds.length;
  const K = KvaluesIdeal(compounds, T_K, P_Pa);

  // Bubble point check: Σ zi·Ki
  const sumZK = z.reduce((s, zi, i) => s + zi * K[i], 0);
  // Dew point check: Σ zi/Ki
  const sumZoverK = z.reduce((s, zi, i) => s + zi / K[i], 0);

  let beta: number;
  let phase: Phase;
  const x = new Array(n).fill(0);
  const y = new Array(n).fill(0);

  if (sumZK <= 1) {
    // Subcooled liquid (below bubble point)
    beta = 0;
    phase = 'Liquid';
    for (let i = 0; i < n; i++) { x[i] = z[i]; y[i] = z[i] * K[i]; }
    // Normalize y
    const sy = y.reduce((a, b) => a + b, 0);
    if (sy > 0) for (let i = 0; i < n; i++) y[i] /= sy;
  } else if (sumZoverK <= 1) {
    // Superheated vapor (above dew point)
    beta = 1;
    phase = 'Vapor';
    for (let i = 0; i < n; i++) { y[i] = z[i]; x[i] = z[i] / K[i]; }
    const sx = x.reduce((a, b) => a + b, 0);
    if (sx > 0) for (let i = 0; i < n; i++) x[i] /= sx;
  } else {
    // Two-phase region
    beta = rachfordRice(z, K);
    phase = beta > 0.999 ? 'Vapor' : beta < 0.001 ? 'Liquid' : 'Two-Phase';
    for (let i = 0; i < n; i++) {
      x[i] = z[i] / (1 + beta * (K[i] - 1));
      y[i] = K[i] * x[i];
    }
    // Normalize
    const sx = x.reduce((a, b) => a + b, 0);
    const sy = y.reduce((a, b) => a + b, 0);
    if (sx > 0) for (let i = 0; i < n; i++) x[i] /= sx;
    if (sy > 0) for (let i = 0; i < n; i++) y[i] /= sy;
  }

  const H = streamEnthalpy(compounds, T_K, beta, x, y);

  return { converged: true, phase, vaporFraction: beta, x, y, K, T_K, P_Pa, H_Jpmol: H };
}

// ────────────────────────────────────────────────────────────────
// 7. Bubble & Dew Point Solvers (Ideal)
// ────────────────────────────────────────────────────────────────

/** Bubble point temperature at given P (Raoult's Law). Newton-Raphson. */
export function bubblePointT_Ideal(
  compounds: CanopyCompound[],
  z: number[],
  P_Pa: number,
): number | null {
  // Objective: Σ zi·Psat_i(T) / P = 1
  let T = 350; // Initial guess
  for (let iter = 0; iter < 100; iter++) {
    let f = -1;
    let df = 0;
    for (let i = 0; i < compounds.length; i++) {
      const ps = Psat_Pa(compounds[i], T);
      if (!isFinite(ps)) continue;
      f += z[i] * ps / P_Pa;
      // Numerical derivative
      const dT = 0.01;
      const ps2 = Psat_Pa(compounds[i], T + dT);
      if (isFinite(ps2)) df += z[i] * (ps2 - ps) / (dT * P_Pa);
    }
    if (Math.abs(f) < 1e-8) return T;
    if (Math.abs(df) < 1e-15) break;
    T -= f / df;
    if (T < 100 || T > 1000) break;
  }
  return null;
}

/** Dew point temperature at given P (Raoult's Law). Newton-Raphson. */
export function dewPointT_Ideal(
  compounds: CanopyCompound[],
  z: number[],
  P_Pa: number,
): number | null {
  // Objective: Σ zi·P / Psat_i(T) = 1  ⟹  Σ zi / Ki = 1
  let T = 380;
  for (let iter = 0; iter < 100; iter++) {
    let f = -1;
    let df = 0;
    for (let i = 0; i < compounds.length; i++) {
      const ps = Psat_Pa(compounds[i], T);
      if (!isFinite(ps) || ps < 1e-6) continue;
      f += z[i] * P_Pa / ps;
      const dT = 0.01;
      const ps2 = Psat_Pa(compounds[i], T + dT);
      if (isFinite(ps2) && ps2 > 1e-6) df += z[i] * P_Pa * (1 / ps2 - 1 / ps) / dT;
    }
    if (Math.abs(f) < 1e-8) return T;
    if (Math.abs(df) < 1e-15) break;
    T -= f / df;
    if (T < 100 || T > 1000) break;
  }
  return null;
}
