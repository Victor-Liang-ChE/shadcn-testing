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

/** ln(P* [Pa]) = C1 + C2/(T+C3) + C4·T + C5·ln(T) + C6·T^C7
 *  For T >= Tc, returns Pc (critical pressure) – Psat is undefined above Tc.
 */
export function Psat_Pa(comp: CanopyCompound, T_K: number): number {
  const a = comp.antoine;
  if (!a || T_K <= 0) return NaN;
  // Above Tc: Psat is undefined; return Pc for K-value estimation
  if (comp.Tc_K && comp.Pc_Pa && T_K >= comp.Tc_K) return comp.Pc_Pa;
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
function CpIG_DIPPR16_JkmolK(comp: CanopyCompound, T_K: number): number {
  const c = comp.dipprCoeffs['IdealGasHeatCapacityCp'];
  if (!c) return NaN;
  const { A, B, C, D, E } = c;
  const E_ = E ?? 0;
  return A + Math.exp(B / T_K + C + D * T_K + E_ * T_K * T_K);
}

/**
 * DIPPR Eq. 107 (Aly-Lee):
 *   Cp_ig = A + B·[(C/T)/sinh(C/T)]² + D·[(E/T)/cosh(E/T)]²
 * Returns J/(kmol·K).
 */
function CpIG_DIPPR107_JkmolK(comp: CanopyCompound, T_K: number): number {
  const c = comp.cpigdp;
  if (!c) return NaN;
  const CT = c.C / T_K;
  const ET = c.E / T_K;
  // Numerically stable sinh/cosh terms
  const sinhCT = Math.sinh(CT);
  const coshET = Math.cosh(ET);
  const termB = sinhCT !== 0 ? Math.pow(CT / sinhCT, 2) : 0;
  const termD = coshET !== 0 ? Math.pow(ET / coshET, 2) : 0;
  return c.A + c.B * termB + c.D * termD;
}

/**
 * Ideal gas heat capacity in J/(kmol·K).
 * Prefers DIPPR 107 (Aly-Lee / CPIGDP) if available, falls back to DIPPR 16.
 */
export function CpIG_JkmolK(comp: CanopyCompound, T_K: number): number {
  if (comp.cpigdp) return CpIG_DIPPR107_JkmolK(comp, T_K);
  return CpIG_DIPPR16_JkmolK(comp, T_K);
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
  const Tr = T_K / comp.Tc_K;
  if (Tr >= 1) return 0; // Above critical — no vaporization
  const c = comp.dipprCoeffs['HeatOfVaporization'];
  if (c && (c.A !== 0 || c.B !== 0)) {
    const { A, B, C, D, E } = c;
    const E_ = E ?? 0;
    return A * Math.pow(1 - Tr, B + C * Tr + D * Tr * Tr + E_ * Tr * Tr * Tr);
  }
  // Pitzer correlation fallback when DIPPR coefficients are missing/zero
  // ΔHvap = R·Tc·(7.08·(1-Tr)^0.354 + 10.95·ω·(1-Tr)^0.456)  [J/mol → ×1000 for J/kmol]
  const R = 8.314;
  return 1000 * R * comp.Tc_K * (7.08 * Math.pow(1 - Tr, 0.354) + 10.95 * comp.omega * Math.pow(1 - Tr, 0.456));
}

export function Hvap_Jmol(comp: CanopyCompound, T_K: number): number {
  return Hvap_Jkmol(comp, T_K) / 1000;
}

// ────────────────────────────────────────────────────────────────
// 4. DIPPR Liquid Heat Capacity
// ────────────────────────────────────────────────────────────────

/**
 * Liquid heat capacity [J/kmol/K] from DIPPR correlation.
 *
 * Equation form depends on eqNo stored with the coefficients:
 *   eqNo=100:           DIPPR 100 polynomial  Cp = A + B·T + C·T² + D·T³ + E·T⁴
 *   eqNo=114:           DIPPR 114 Zabransky    Cp = A²/τ + B - 2ACτ - ADτ² - C²τ³/3 - CDτ⁴/2 - D²τ⁵/5
 *                       where τ = 1 - T/Tc
 *   (no eqNo, default): ChemSep/PPDS form      Cp = A + exp(B/T + C + D·T + E·T²)
 *
 * APV140 CPLDIP in PURE40 uses eq 100 for most compounds and eq 114 for some (e.g. propane).
 */
export function CpL_JkmolK(comp: CanopyCompound, T_K: number): number {
  const c = comp.dipprCoeffs['LiquidHeatCapacityCp'];
  if (!c) {
    // Fallback: use ideal gas Cp (rough approximation)
    return CpIG_JkmolK(comp, T_K);
  }
  const { A, B, C, D, E } = c;
  const E_ = E ?? 0;
  const eqNo = c.eqNo;

  if (eqNo === 114) {
    // DIPPR 114 (Zabransky): Cp = A²/τ + B - 2ACτ - ADτ² - C²τ³/3 - CDτ⁴/2 - D²τ⁵/5
    const tau = 1 - T_K / comp.Tc_K;
    if (tau <= 0) return CpIG_JkmolK(comp, T_K); // above critical
    const t2 = tau * tau, t3 = t2 * tau, t4 = t3 * tau, t5 = t4 * tau;
    return A * A / tau + B - 2 * A * C * tau - A * D * t2
      - C * C * t3 / 3 - C * D * t4 / 2 - D * D * t5 / 5;
  }

  if (eqNo === 100) {
    // DIPPR 100 polynomial: Cp = A + B·T + C·T² + D·T³ + E·T⁴
    return A + B * T_K + C * T_K * T_K + D * T_K * T_K * T_K + E_ * Math.pow(T_K, 4);
  }

  // Default (legacy/ChemSep form): Cp = A + exp(B/T + C + D·T + E·T²)
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
 * Ideal-gas molar enthalpy using Aspen elements reference state.
 * H_ig(T) = Hf°(298.15 K) + ∫[T_ref → T] Cp_ig dT   (J/mol)
 *
 * When Hf° is available, enthalpy is referenced to elements in their
 * standard states. This makes heat of reaction emerge automatically from
 * the energy balance (H_out - H_in) without needing a separate ΔH_rxn input.
 */
export function enthalpyIG(comp: CanopyCompound, T_K: number): number {
  const Hf = comp.Hf_298_Jmol ?? 0;
  return Hf + integrateCpIG(comp, T_K);
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
  P_Pa?: number,
  fluidPackage?: FluidPackageType,
  kij?: number[][],
  interactionParams?: InteractionParams,
): number {
  const useEOS = P_Pa != null && (fluidPackage === 'Peng-Robinson' || fluidPackage === 'SRK');
  let H = 0;
  if (useEOS) {
    // EOS enthalpy: H = H_ig + H_dep  (departure replaces ΔHvap approach)
    for (let i = 0; i < compounds.length; i++) {
      const Hig = enthalpyIG(compounds[i], T_K);
      H += vaporFraction * y_vapor[i] * Hig + (1 - vaporFraction) * x_liquid[i] * Hig;
    }
    const fn = fluidPackage === 'Peng-Robinson' ? departureEnthalpyPR : departureEnthalpySRK;
    if (vaporFraction > 0) {
      H += vaporFraction * fn(compounds, y_vapor, T_K, P_Pa!, true, kij);
    }
    if (vaporFraction < 1) {
      H += (1 - vaporFraction) * fn(compounds, x_liquid, T_K, P_Pa!, false, kij);
    }
  } else {
    for (let i = 0; i < compounds.length; i++) {
      const Hv = enthalpyIG(compounds[i], T_K);
      const Hl = enthalpyLiquid(compounds[i], T_K);
      H += vaporFraction * y_vapor[i] * Hv + (1 - vaporFraction) * x_liquid[i] * Hl;
    }
    // Add excess enthalpy for activity coefficient models (liquid only)
    if (vaporFraction < 1 && fluidPackage && fluidPackage !== 'Ideal') {
      let gammaFn: ((x: number[], T: number) => number[]) | null = null;
      if (fluidPackage === 'NRTL' && interactionParams?.nrtl) {
        const p = interactionParams.nrtl;
        gammaFn = (xf, Tf) => nrtlGamma(xf, Tf, p);
      } else if (fluidPackage === 'Wilson' && interactionParams?.wilson) {
        const p = interactionParams.wilson;
        gammaFn = (xf, Tf) => wilsonGamma(xf, Tf, p);
      } else if (fluidPackage === 'UNIQUAC' && interactionParams?.uniquac) {
        const p = interactionParams.uniquac;
        gammaFn = (xf, Tf) => uniquacGamma(xf, Tf, p);
      }
      if (gammaFn) {
        H += (1 - vaporFraction) * excessEnthalpy(x_liquid, T_K, gammaFn);
      }
    }
  }
  return H;
}

// ────────────────────────────────────────────────────────────────
// 5b. EOS Departure Functions (Peng-Robinson & SRK)
// ────────────────────────────────────────────────────────────────

/**
 * Peng-Robinson departure enthalpy for a single phase (J/mol).
 *
 * H_dep = RT(Z-1) + [T·(da_mix/dT) - a_mix] / (2√2·b_mix)
 *         · ln[(Z+(1+√2)B) / (Z+(1-√2)B)]
 */
export function departureEnthalpyPR(
  compounds: CanopyCompound[],
  z: number[],
  T_K: number,
  P_Pa: number,
  isVapor: boolean,
  kij?: number[][],
): number {
  const n = compounds.length;
  const a_i = new Array(n);
  const b_i = new Array(n);
  const da_i_dT = new Array(n);

  for (let i = 0; i < n; i++) {
    const { Tc_K, Pc_Pa } = compounds[i];
    const { alpha, dAlpha_dT } = alphaForCompound(compounds[i], T_K, 'PR');
    const Omega_a = 0.45724;
    a_i[i] = Omega_a * R * R * Tc_K * Tc_K / Pc_Pa * alpha;
    b_i[i] = 0.07780 * R * Tc_K / Pc_Pa;
    da_i_dT[i] = Omega_a * R * R * Tc_K * Tc_K / Pc_Pa * dAlpha_dT;
  }

  // Mixing rules
  let a_mix = 0, b_mix = 0, da_mix_dT = 0;
  for (let i = 0; i < n; i++) {
    b_mix += z[i] * b_i[i];
    for (let j = 0; j < n; j++) {
      const k = kij?.[i]?.[j] ?? 0;
      const aij = (1 - k) * Math.sqrt(a_i[i] * a_i[j]);
      a_mix += z[i] * z[j] * aij;
      const sqrtProd = Math.sqrt(a_i[i] * a_i[j]);
      const daij_dT = sqrtProd > 1e-30
        ? (1 - k) * (da_i_dT[i] * a_i[j] + a_i[i] * da_i_dT[j]) / (2 * sqrtProd)
        : 0;
      da_mix_dT += z[i] * z[j] * daij_dT;
    }
  }

  if (b_mix < 1e-30) return 0;

  const A_ = a_mix * P_Pa / (R * R * T_K * T_K);
  const B_ = b_mix * P_Pa / (R * T_K);

  // Solve PR cubic: Z³ - (1-B)Z² + (A-3B²-2B)Z - (AB-B²-B³) = 0
  const c2 = -(1 - B_);
  const c1 = A_ - 3 * B_ * B_ - 2 * B_;
  const c0 = -(A_ * B_ - B_ * B_ - B_ * B_ * B_);
  const roots = solveCubicReal(1, c2, c1, c0);
  const validRoots = roots.filter(r => r > B_);
  const Z = validRoots.length > 0
    ? (isVapor ? Math.max(...validRoots) : Math.min(...validRoots))
    : roots.length > 0 ? Math.max(...roots) : 1;

  const sqrt2 = Math.SQRT2;
  const d1 = 1 + sqrt2;  // δ₁ = 1+√2
  const d2 = 1 - sqrt2;  // δ₂ = 1-√2
  const arg1 = Z + d1 * B_;
  const arg2 = Z + d2 * B_;
  const logTerm = (arg1 > 0 && arg2 > 0) ? Math.log(arg1 / arg2) : 0;

  return R * T_K * (Z - 1)
    + (T_K * da_mix_dT - a_mix) / (2 * sqrt2 * b_mix) * logTerm;
}

/**
 * SRK departure enthalpy for a single phase (J/mol).
 *
 * H_dep = RT(Z-1) + [T·(da_mix/dT) - a_mix] / b_mix · ln(1 + B/Z)
 */
export function departureEnthalpySRK(
  compounds: CanopyCompound[],
  z: number[],
  T_K: number,
  P_Pa: number,
  isVapor: boolean,
  kij?: number[][],
): number {
  const n = compounds.length;
  const a_i = new Array(n);
  const b_i = new Array(n);
  const da_i_dT = new Array(n);

  for (let i = 0; i < n; i++) {
    const { Tc_K, Pc_Pa } = compounds[i];
    const { alpha, dAlpha_dT } = alphaForCompound(compounds[i], T_K, 'SRK');
    const Omega_a = 0.42748;
    a_i[i] = Omega_a * R * R * Tc_K * Tc_K / Pc_Pa * alpha;
    b_i[i] = 0.08664 * R * Tc_K / Pc_Pa;
    da_i_dT[i] = Omega_a * R * R * Tc_K * Tc_K / Pc_Pa * dAlpha_dT;
  }

  let a_mix = 0, b_mix = 0, da_mix_dT = 0;
  for (let i = 0; i < n; i++) {
    b_mix += z[i] * b_i[i];
    for (let j = 0; j < n; j++) {
      const k = kij?.[i]?.[j] ?? 0;
      const aij = (1 - k) * Math.sqrt(a_i[i] * a_i[j]);
      a_mix += z[i] * z[j] * aij;
      const sqrtProd = Math.sqrt(a_i[i] * a_i[j]);
      const daij_dT = sqrtProd > 1e-30
        ? (1 - k) * (da_i_dT[i] * a_i[j] + a_i[i] * da_i_dT[j]) / (2 * sqrtProd)
        : 0;
      da_mix_dT += z[i] * z[j] * daij_dT;
    }
  }

  if (b_mix < 1e-30) return 0;

  const A_ = a_mix * P_Pa / (R * R * T_K * T_K);
  const B_ = b_mix * P_Pa / (R * T_K);

  // SRK cubic: Z³ - Z² + (A-B-B²)Z - AB = 0
  const c2 = -1;
  const c1 = A_ - B_ - B_ * B_;
  const c0 = -A_ * B_;
  const roots = solveCubicReal(1, c2, c1, c0);
  const validRoots = roots.filter(r => r > B_);
  const Z = validRoots.length > 0
    ? (isVapor ? Math.max(...validRoots) : Math.min(...validRoots))
    : roots.length > 0 ? Math.max(...roots) : 1;

  const logTerm = Z > 0 ? Math.log(1 + B_ / Z) : 0;

  return R * T_K * (Z - 1)
    + (T_K * da_mix_dT - a_mix) / b_mix * logTerm;
}

/**
 * Peng-Robinson departure entropy for a single phase (J/(mol·K)).
 *
 * S_dep = R·ln(Z-B) + (da_mix/dT) / (2√2·b_mix)
 *         · ln[(Z+(1+√2)B) / (Z+(1-√2)B)]
 */
export function departureEntropyPR(
  compounds: CanopyCompound[],
  z: number[],
  T_K: number,
  P_Pa: number,
  isVapor: boolean,
  kij?: number[][],
): number {
  const n = compounds.length;
  const a_i = new Array(n);
  const b_i = new Array(n);
  const da_i_dT = new Array(n);

  for (let i = 0; i < n; i++) {
    const { Tc_K, Pc_Pa } = compounds[i];
    const { alpha, dAlpha_dT } = alphaForCompound(compounds[i], T_K, 'PR');
    const Omega_a = 0.45724;
    a_i[i] = Omega_a * R * R * Tc_K * Tc_K / Pc_Pa * alpha;
    b_i[i] = 0.07780 * R * Tc_K / Pc_Pa;
    da_i_dT[i] = Omega_a * R * R * Tc_K * Tc_K / Pc_Pa * dAlpha_dT;
  }

  let a_mix = 0, b_mix = 0, da_mix_dT = 0;
  for (let i = 0; i < n; i++) {
    b_mix += z[i] * b_i[i];
    for (let j = 0; j < n; j++) {
      const k = kij?.[i]?.[j] ?? 0;
      a_mix += z[i] * z[j] * (1 - k) * Math.sqrt(a_i[i] * a_i[j]);
      const sqrtProd = Math.sqrt(a_i[i] * a_i[j]);
      const daij_dT = sqrtProd > 1e-30
        ? (1 - k) * (da_i_dT[i] * a_i[j] + a_i[i] * da_i_dT[j]) / (2 * sqrtProd)
        : 0;
      da_mix_dT += z[i] * z[j] * daij_dT;
    }
  }

  if (b_mix < 1e-30) return 0;

  const B_ = b_mix * P_Pa / (R * T_K);
  const A_ = a_mix * P_Pa / (R * R * T_K * T_K);

  const c2 = -(1 - B_);
  const c1 = A_ - 3 * B_ * B_ - 2 * B_;
  const c0 = -(A_ * B_ - B_ * B_ - B_ * B_ * B_);
  const roots = solveCubicReal(1, c2, c1, c0);
  const validRoots = roots.filter(r => r > B_);
  const Z = validRoots.length > 0
    ? (isVapor ? Math.max(...validRoots) : Math.min(...validRoots))
    : roots.length > 0 ? Math.max(...roots) : 1;

  const sqrt2 = Math.SQRT2;
  const d1 = 1 + sqrt2;
  const d2 = 1 - sqrt2;
  const arg1 = Z + d1 * B_;
  const arg2 = Z + d2 * B_;
  const logTerm = (arg1 > 0 && arg2 > 0) ? Math.log(arg1 / arg2) : 0;

  return R * Math.log(Math.max(Z - B_, 1e-30))
    + da_mix_dT / (2 * sqrt2 * b_mix) * logTerm;
}

/**
 * SRK departure entropy for a single phase (J/(mol·K)).
 *
 * S_dep = R·ln(Z-B) + (da_mix/dT) / b_mix · ln(1 + B/Z)
 */
export function departureEntropySRK(
  compounds: CanopyCompound[],
  z: number[],
  T_K: number,
  P_Pa: number,
  isVapor: boolean,
  kij?: number[][],
): number {
  const n = compounds.length;
  const a_i = new Array(n);
  const b_i = new Array(n);
  const da_i_dT = new Array(n);

  for (let i = 0; i < n; i++) {
    const { Tc_K, Pc_Pa } = compounds[i];
    const { alpha, dAlpha_dT } = alphaForCompound(compounds[i], T_K, 'SRK');
    const Omega_a = 0.42748;
    a_i[i] = Omega_a * R * R * Tc_K * Tc_K / Pc_Pa * alpha;
    b_i[i] = 0.08664 * R * Tc_K / Pc_Pa;
    da_i_dT[i] = Omega_a * R * R * Tc_K * Tc_K / Pc_Pa * dAlpha_dT;
  }

  let a_mix = 0, b_mix = 0, da_mix_dT = 0;
  for (let i = 0; i < n; i++) {
    b_mix += z[i] * b_i[i];
    for (let j = 0; j < n; j++) {
      const k = kij?.[i]?.[j] ?? 0;
      a_mix += z[i] * z[j] * (1 - k) * Math.sqrt(a_i[i] * a_i[j]);
      const sqrtProd = Math.sqrt(a_i[i] * a_i[j]);
      const daij_dT = sqrtProd > 1e-30
        ? (1 - k) * (da_i_dT[i] * a_i[j] + a_i[i] * da_i_dT[j]) / (2 * sqrtProd)
        : 0;
      da_mix_dT += z[i] * z[j] * daij_dT;
    }
  }

  if (b_mix < 1e-30) return 0;

  const B_ = b_mix * P_Pa / (R * T_K);
  const A_ = a_mix * P_Pa / (R * R * T_K * T_K);

  const c2 = -1;
  const c1 = A_ - B_ - B_ * B_;
  const c0 = -A_ * B_;
  const roots = solveCubicReal(1, c2, c1, c0);
  const validRoots = roots.filter(r => r > B_);
  const Z = validRoots.length > 0
    ? (isVapor ? Math.max(...validRoots) : Math.min(...validRoots))
    : roots.length > 0 ? Math.max(...roots) : 1;

  const logTerm = Z > 0 ? Math.log(1 + B_ / Z) : 0;

  return R * Math.log(Math.max(Z - B_, 1e-30))
    + da_mix_dT / b_mix * logTerm;
}

/**
 * Ideal-gas molar entropy relative to reference state (J/(mol·K)).
 * S_ig(T,P) = ∫[T_ref→T] Cp_ig/T dT − R·ln(P/P_ref)
 */
function integrateCpOverT(comp: CanopyCompound, T_K: number): number {
  if (Math.abs(T_K - T_REF) < 0.01) return 0;
  const n = 100;
  const h = (T_K - T_REF) / n;
  let sum = CpIG_JmolK(comp, T_REF) / T_REF + CpIG_JmolK(comp, T_K) / T_K;
  for (let i = 1; i < n; i++) {
    const T = T_REF + i * h;
    sum += (i % 2 === 0 ? 2 : 4) * CpIG_JmolK(comp, T) / T;
  }
  return sum * h / 3;
}

const P_REF = 101325; // Pa (1 atm reference for entropy)

export function entropyIG(comp: CanopyCompound, T_K: number, P_Pa: number): number {
  return integrateCpOverT(comp, T_K) - R * Math.log(P_Pa / P_REF);
}

/**
 * Stream molar entropy (J/(mol·K)).
 * For EOS packages: S = Σ zi · S_ig,i − R·Σ zi·ln(zi) + S_dep
 * For ideal/activity models: S = Σ (β·yi·Sig,i + (1−β)·xi·Sig,i) − R·Σ zi·ln(zi)
 */
export function streamEntropy(
  compounds: CanopyCompound[],
  T_K: number,
  P_Pa: number,
  vaporFraction: number,
  x_liquid: number[],
  y_vapor: number[],
  fluidPackage?: FluidPackageType,
  kij?: number[][],
): number {
  let S = 0;
  const useEOS = fluidPackage === 'Peng-Robinson' || fluidPackage === 'SRK';
  // Ideal-gas contribution
  for (let i = 0; i < compounds.length; i++) {
    const Sig = entropyIG(compounds[i], T_K, P_Pa);
    S += vaporFraction * y_vapor[i] * Sig + (1 - vaporFraction) * x_liquid[i] * Sig;
  }
  // Ideal mixing entropy: −R·Σ zi·ln(zi)
  for (let i = 0; i < compounds.length; i++) {
    const zi = vaporFraction * y_vapor[i] + (1 - vaporFraction) * x_liquid[i];
    if (zi > 1e-30) S -= R * zi * Math.log(zi);
  }
  if (useEOS) {
    const fn = fluidPackage === 'Peng-Robinson' ? departureEntropyPR : departureEntropySRK;
    if (vaporFraction > 0) {
      S += vaporFraction * fn(compounds, y_vapor, T_K, P_Pa, true, kij);
    }
    if (vaporFraction < 1) {
      S += (1 - vaporFraction) * fn(compounds, x_liquid, T_K, P_Pa, false, kij);
    }
  }
  return S;
}

// ────────────────────────────────────────────────────────────────
// Gibbs Free Energy & Chemical Potential
// ────────────────────────────────────────────────────────────────

/**
 * Ideal gas Gibbs free energy (J/mol) at T and P:
 *   G_ig = H_ig - T·S_ig
 */
export function gibbsIG(comp: CanopyCompound, T_K: number, P_Pa: number): number {
  return enthalpyIG(comp, T_K) - T_K * entropyIG(comp, T_K, P_Pa);
}

/**
 * Mixture Gibbs free energy (J/mol) for an ideal gas mixture:
 *   G_mix = Σ yi·(G_ig_i + R·T·ln(yi))
 */
export function mixtureGibbsIG(
  compounds: CanopyCompound[], y: number[], T_K: number, P_Pa: number,
): number {
  let G = 0;
  for (let i = 0; i < compounds.length; i++) {
    const yi = Math.max(y[i], 1e-30);
    G += yi * (gibbsIG(compounds[i], T_K, P_Pa) + R * T_K * Math.log(yi));
  }
  return G;
}

/**
 * Chemical potential μ_i (J/mol) for component i in an ideal gas mixture:
 *   μ_i = G_ig_i + R·T·ln(y_i·P/P_ref)
 */
export function chemicalPotentialIG(
  comp: CanopyCompound, T_K: number, P_Pa: number, y_i: number,
): number {
  return gibbsIG(comp, T_K, P_Pa) + R * T_K * Math.log(Math.max(y_i, 1e-30));
}

/**
 * Chemical potential μ_i (J/mol) for component i in a liquid with activity coefficients:
 *   μ_i = μ_i^pure_liquid + R·T·ln(γ_i·x_i)
 * where μ_i^pure_liquid = G_ig_i(T,Pref) - ΔG_vap (approx = G_ig - R·T·ln(Psat/P))
 */
export function chemicalPotentialLiquid(
  comp: CanopyCompound, T_K: number, x_i: number, gamma_i: number,
): number {
  const Psat = Psat_Pa(comp, T_K);
  const mu_pure = gibbsIG(comp, T_K, P_REF) - R * T_K * Math.log(Math.max(Psat / P_REF, 1e-30));
  return mu_pure + R * T_K * Math.log(Math.max(gamma_i * x_i, 1e-30));
}

// ────────────────────────────────────────────────────────────────
// Mixture Critical Properties (Li's method)
// ────────────────────────────────────────────────────────────────

/**
 * Mixture pseudo-critical temperature (K) — Kay's rule.
 * Tc_mix = Σ zi·Tc_i
 */
export function mixtureTc(compounds: CanopyCompound[], z: number[]): number {
  return compounds.reduce((s, c, i) => s + z[i] * c.Tc_K, 0);
}

/**
 * Mixture pseudo-critical pressure (Pa) — Kay's rule.
 * Pc_mix = Σ zi·Pc_i
 */
export function mixturePc(compounds: CanopyCompound[], z: number[]): number {
  return compounds.reduce((s, c, i) => s + z[i] * c.Pc_Pa, 0);
}

/**
 * Mixture acentric factor — mole-fraction weighted.
 * ω_mix = Σ zi·ω_i
 */
export function mixtureOmega(compounds: CanopyCompound[], z: number[]): number {
  return compounds.reduce((s, c, i) => s + z[i] * c.omega, 0);
}

/**
 * Li's mixing rule for mixture critical volume (m³/mol):
 * Vc_mix = Σ_i Σ_j zi·zj·Vc_ij, Vc_ij = (Vc_i^(1/3) + Vc_j^(1/3))³/8
 */
export function mixtureVc(
  Vc: number[], z: number[],
): number {
  const n = z.length;
  let Vc_mix = 0;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const Vc_ij = Math.pow((Math.pow(Vc[i], 1 / 3) + Math.pow(Vc[j], 1 / 3)) / 2, 3);
      Vc_mix += z[i] * z[j] * Vc_ij;
    }
  }
  return Vc_mix;
}

/**
 * Isentropic (PS) flash — find T at given P such that S(T,P) = S_target.
 * Used by compressors and turbines for isentropic efficiency calculations.
 */
export function flashPS(
  compounds: CanopyCompound[],
  z: number[],
  P_Pa: number,
  S_target: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
  T_guess?: number,
): FlashResult {
  const evalS = (T: number) => {
    const f = flashPT(compounds, z, T, P_Pa, fluidPackage, interactionParams);
    return { flash: f, S: streamEntropy(compounds, T, P_Pa, f.vaporFraction, f.x, f.y, fluidPackage, interactionParams?.kij) };
  };

  let Tlo = 100, Thi = 1500;
  let T = T_guess ?? 350;
  let bestT = T, bestErr = Infinity;

  // Newton + bracket
  for (let iter = 0; iter < 60; iter++) {
    const { flash, S } = evalS(T);
    const err = S - S_target;
    if (Math.abs(err) < Math.abs(bestErr)) { bestT = T; bestErr = err; }
    if (Math.abs(err) < 1e-4) return { ...flash, T_K: T };

    if (err > 0 && T < Thi) Thi = T;
    if (err < 0 && T > Tlo) Tlo = T;

    const dT = 0.1;
    const { S: S2 } = evalS(T + dT);
    const dSdT = (S2 - S) / dT;

    if (Math.abs(dSdT) > 1e-12) {
      let T_new = T - err / dSdT;
      if (T_new < Tlo || T_new > Thi) T_new = (Tlo + Thi) / 2;
      T = T_new;
    } else {
      T = (Tlo + Thi) / 2;
    }
  }

  // Bisection fallback
  for (let iter = 0; iter < 80; iter++) {
    const mid = (Tlo + Thi) / 2;
    const { flash, S: S_mid } = evalS(mid);
    const err = S_mid - S_target;
    if (Math.abs(err) < Math.abs(bestErr)) { bestT = mid; bestErr = err; }
    if (Math.abs(err) < 1e-4 || (Thi - Tlo) < 0.001) break;
    if (err > 0) Thi = mid; else Tlo = mid;
  }

  return { ...flashPT(compounds, z, bestT, P_Pa, fluidPackage, interactionParams), T_K: bestT };
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
 * Wilson K-value estimate for column initialization.
 * Ki = (Pc_i / P) · exp[5.373·(1 + ω_i)·(1 - Tc_i/T)]
 * Used by Aspen Plus for RadFrac initial estimates (Section 1.4 of Aspen docs).
 * Ref: Wilson, G.M., NGSMA 1968 meeting.
 */
export function KvaluesWilsonEstimate(compounds: CanopyCompound[], T_K: number, P_Pa: number): number[] {
  return compounds.map(c => {
    const K = (c.Pc_Pa / P_Pa) * Math.exp(5.373 * (1 + c.omega) * (1 - c.Tc_K / T_K));
    return isFinite(K) && K > 0 ? K : 1;
  });
}

// ────────────────────────────────────────────────────────────────
// 6b. Non-Ideal K-values — NRTL, Wilson, PR EOS
// ────────────────────────────────────────────────────────────────

/** NRTL binary interaction parameter set for multi-component systems.
 *  Parameters are stored as n×n matrices (row=i, col=j).
 *
 *  Simple form:  τ_ij = A_ij/(R·T) + B_ij
 *  Extended form: τ_ij = a_ij + b_ij/T + e_ij·ln(T) + f_ij·T
 *
 *  When extended parameters (a,b,e,f) are provided, they take precedence.
 */
export interface NRTLParams {
  /** Energy parameter A_ij (cal/mol) — simple form */
  A?: number[][];
  /** Energy parameter B_ij (dimensionless — multiplied by T) — simple form */
  B?: number[][];
  /** Non-randomness factor α_ij = c_ij (constant part) */
  alpha: number[][];
  /** Extended: a_ij (intercept, dimensionless) */
  a_ext?: number[][];
  /** Extended: b_ij (1/T coefficient, K) */
  b_ext?: number[][];
  /** Extended: e_ij (ln(T) coefficient) */
  e_ext?: number[][];
  /** Extended: f_ij (T coefficient, 1/K) */
  f_ext?: number[][];
  /** Temperature-dependent α: d_ij such that α_ij(T) = c_ij + d_ij·(T − 273.15) */
  d_alpha?: number[][];
  /** Temperature validity lower bound (K) */
  Tlower?: number;
  /** Temperature validity upper bound (K) */
  Tupper?: number;
}

/**
 * Multi-component NRTL activity coefficients.
 *
 * τ_ij = A_ij/(R·T) + B_ij    where R=1.98721 cal/(mol·K)
 * G_ij = exp(-α_ij · τ_ij)
 *
 * ln(γ_i) = [Σ_j τ_ji G_ji x_j] / [Σ_l G_li x_l]
 *         + Σ_j { (x_j G_ij) / (Σ_l G_lj x_l) · [τ_ij - Σ_m(x_m τ_mj G_mj) / Σ_l(G_lj x_l)] }
 */
export function nrtlGamma(
  x: number[],
  T_K: number,
  params: NRTLParams,
): number[] {
  const n = x.length;
  const R_cal = 1.98721; // cal/(mol·K)

  // Build τ, α, and G matrices
  // Aspen 12-element: τ_ij = a + b/T + e·ln(T) + f·T;  α_ij = c + d·(T − 273.15)
  const tau: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
  const alphaM: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
  const G: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
  const hasExtended = !!(params.a_ext || params.b_ext || params.e_ext || params.f_ext);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      // Resolve temperature-dependent α_ij = c_ij + d_ij·(T − 273.15)
      const c_ij = params.alpha[i]?.[j] ?? 0.3;
      const d_ij = params.d_alpha?.[i]?.[j] ?? 0;
      alphaM[i][j] = c_ij + d_ij * (T_K - 273.15);

      if (i === j) {
        tau[i][j] = 0;
        G[i][j] = 1;
      } else if (hasExtended) {
        // Extended NRTL: τ_ij = a_ij + b_ij/T + e_ij·ln(T) + f_ij·T
        const a_val = params.a_ext?.[i]?.[j] ?? 0;
        const b_val = params.b_ext?.[i]?.[j] ?? 0;
        const e_val = params.e_ext?.[i]?.[j] ?? 0;
        const f_val = params.f_ext?.[i]?.[j] ?? 0;
        tau[i][j] = a_val + b_val / T_K + e_val * Math.log(T_K) + f_val * T_K;
        G[i][j] = Math.exp(-alphaM[i][j] * tau[i][j]);
      } else {
        tau[i][j] = (params.A?.[i]?.[j] ?? 0) / (R_cal * T_K) + (params.B?.[i]?.[j] ?? 0);
        G[i][j] = Math.exp(-alphaM[i][j] * tau[i][j]);
      }
    }
  }

  // Precompute sums
  const sumGxCol: number[] = new Array(n).fill(0); // Σ_l G_lj x_l for each j
  const sumTauGxCol: number[] = new Array(n).fill(0); // Σ_m τ_mj G_mj x_m for each j
  for (let j = 0; j < n; j++) {
    for (let l = 0; l < n; l++) {
      sumGxCol[j] += G[l][j] * x[l];
      sumTauGxCol[j] += tau[l][j] * G[l][j] * x[l];
    }
  }

  const gamma = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    // Term 1: Σ_j τ_ji G_ji x_j / Σ_l G_li x_l
    let num1 = 0;
    for (let j = 0; j < n; j++) {
      num1 += tau[j][i] * G[j][i] * x[j];
    }
    const term1 = num1 / (sumGxCol[i] || 1e-30);

    // Term 2: Σ_j {(x_j G_ij)/(Σ_l G_lj x_l) · [τ_ij - (Σ_m τ_mj G_mj x_m)/(Σ_l G_lj x_l)]}
    let term2 = 0;
    for (let j = 0; j < n; j++) {
      const denom_j = sumGxCol[j] || 1e-30;
      const bracket = tau[i][j] - sumTauGxCol[j] / denom_j;
      term2 += (x[j] * G[i][j] / denom_j) * bracket;
    }

    const lnGamma = term1 + term2;
    gamma[i] = Math.exp(Math.max(-50, Math.min(50, lnGamma)));
  }

  return gamma;
}

// ────────────────────────────────────────────────────────────────
// Electrolyte-NRTL (simplified Chen model)
// ────────────────────────────────────────────────────────────────

/** Parameters for the Electrolyte-NRTL (Chen) model.
 *  Treats the solution as apparent (molecular) species with long-range
 *  (Debye-Hückel) + short-range (NRTL) contributions.
 */
export interface ElecNRTLParams {
  /** NRTL part (molecular basis) */
  nrtl: NRTLParams;
  /** Charge number for each species (0 = molecular, ±n = ionic) */
  charges: number[];
  /** Dielectric constant of pure solvent (e.g., 78.36 for water at 25°C) */
  epsilon_r: number;
  /** Solvent molecular weight (g/mol), e.g. 18.015 for water */
  solventMW: number;
  /** Optional Born / closest-approach radii for ionic species (m) */
  ionRadii_m?: number[];
}

/**
 * Electrolyte-NRTL activity coefficients (simplified Chen model).
 *
 * Long-range contribution (Pitzer-Debye-Hückel):
 *   ln(γ_i^PDH) = -(Aφ·z_i²·√I) / (1 + ρ·√I)
 *   where Aφ = Debye-Hückel parameter, I = ionic strength, ρ = closest approach (14.9)
 *
 * Short-range contribution: standard NRTL.
 *
 * Total: ln(γ_i) = ln(γ_i^PDH) + ln(γ_i^NRTL)
 */
export function elecNRTLGamma(
  x: number[],
  T_K: number,
  params: ElecNRTLParams,
): number[] {
  const nc = x.length;
  const { charges, epsilon_r, solventMW, ionRadii_m } = params;

  // Short-range (NRTL) contribution
  const gammaNRTL = nrtlGamma(x, T_K, params.nrtl);

  // Ionic strength: I = 0.5 · Σ x_i · z_i²
  let ionicStrength = 0;
  for (let i = 0; i < nc; i++) {
    ionicStrength += 0.5 * x[i] * charges[i] * charges[i];
  }

  if (ionicStrength < 1e-15) return gammaNRTL;  // No ions → pure NRTL

  // Debye-Hückel parameter Aφ
  // Aφ = (1/3)·(2πN_A·d_s/1000)^(1/2) · (e²/(4πε₀εᵣkT))^(3/2)
  // Simplified: Aφ ≈ 1.327757e5 · √(d_s) / (ε_r·T)^(3/2)
  // For dilute aqueous: Aφ ≈ 0.391 at 25°C
  const d_s = 1000 / solventMW; // approximate solvent density proxy (mol/L)
  const Aphi = 1.327757e5 * Math.sqrt(d_s) / Math.pow(epsilon_r * T_K, 1.5);

  const sqrtI = Math.sqrt(ionicStrength);

  const gamma = new Array(nc).fill(0);
  for (let i = 0; i < nc; i++) {
    const radiusAngstrom = Math.max(3, Math.min(25, ((ionRadii_m?.[i] ?? 3.5e-10) * 1e10)));
    const rho_DH = Math.max(4, Math.min(20, 2 * radiusAngstrom));
    // PDH contribution
    let lnGammaPDH: number;
    if (charges[i] !== 0) {
      // Ion: ln(γ±^PDH) = -Aφ·z²·√I / (1 + ρ·√I)
      lnGammaPDH = -Aphi * charges[i] * charges[i] * sqrtI / (1 + rho_DH * sqrtI);
    } else {
      // Molecular species (solvent): Born-type contribution
      // ln(γ_s^PDH) = 2·Aφ·I^(3/2) / (1 + ρ·√I)²  (approximate)
      lnGammaPDH = 2 * Aphi * Math.pow(ionicStrength, 1.5)
        / ((1 + rho_DH * sqrtI) * (1 + rho_DH * sqrtI));
    }

    const lnG = Math.log(Math.max(gammaNRTL[i], 1e-30)) + lnGammaPDH;
    gamma[i] = Math.exp(Math.max(-50, Math.min(50, lnG)));
  }

  return gamma;
}

/**
 * Excess enthalpy H^E (J/mol) from an activity coefficient model.
 * Uses numerical differentiation: H^E = -RT² · ∂(Σxi·ln γi)/∂T
 *
 * This is the heat of mixing for liquid solutions, needed for accurate
 * enthalpy balances with NRTL/Wilson/UNIQUAC fluid packages.
 */
export function excessEnthalpy(
  x: number[],
  T_K: number,
  gammaFn: (x: number[], T: number) => number[],
): number {
  const dT = 0.1; // K
  const gamma = gammaFn(x, T_K);
  const gamma_plus = gammaFn(x, T_K + dT);

  let sumLnG = 0, sumLnG_plus = 0;
  for (let i = 0; i < x.length; i++) {
    sumLnG += x[i] * Math.log(Math.max(gamma[i], 1e-30));
    sumLnG_plus += x[i] * Math.log(Math.max(gamma_plus[i], 1e-30));
  }

  const dSumLnG_dT = (sumLnG_plus - sumLnG) / dT;
  return -R * T_K * T_K * dSumLnG_dT;
}

/**
 * K-values using Modified Raoult's Law with Poynting correction:
 *   K_i = γ_i · Psat_i · exp[V_L_i·(P-Psat_i)/(R·T)] / P
 * For activity coefficient models (NRTL, Wilson, UNIQUAC, UNIFAC).
 * The Poynting factor accounts for pressure effect on liquid fugacity.
 */
export function KvaluesGammaPhi(
  compounds: CanopyCompound[],
  x: number[],
  T_K: number,
  P_Pa: number,
  gamma: number[],
): number[] {
  return compounds.map((c, i) => {
    const ps = Psat_Pa(c, T_K);
    const g = gamma[i] ?? 1;
    if (!isFinite(ps)) return 1;
    // Poynting correction factor: exp[V_L·(P-Psat)/(R·T)]
    const V_L = liquidMolarVolume_m3pmol(c, T_K); // m³/mol
    let poynting = 1;
    if (isFinite(V_L) && V_L > 0) {
      const exponent = V_L * (P_Pa - ps) / (R * T_K);
      // Only apply if correction is meaningful (high P)
      if (Math.abs(exponent) < 50) {
        poynting = Math.exp(exponent);
      }
    }
    return g * ps * poynting / P_Pa;
  });
}

// ────────────────────────────────────────────────────────────────
// Alpha function with Boston-Mathias supercritical extrapolation
// ────────────────────────────────────────────────────────────────

/**
 * Standard alpha function:
 *   T ≤ Tc: α = [1 + κ·(1 - √Tr)]²
 *   T > Tc (Boston-Mathias): α = exp[c₁·(1 - Tr^c₂)]
 *     where c₁ = κ/(1 + κ/2), c₂ = 1 + κ/2
 *
 * Returns { alpha, dAlpha_dT } (dα/dT in 1/K).
 */
export function alphaPR(Tc_K: number, omega: number, T_K: number): { alpha: number; dAlpha_dT: number } {
  const kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
  const Tr = T_K / Tc_K;
  if (Tr <= 1) {
    const sqrtTr = Math.sqrt(Tr);
    const f = 1 + kappa * (1 - sqrtTr);
    const alpha = f * f;
    // dα/dT = 2·f · (-κ/(2·√(T·Tc))) = -f·κ / √(T·Tc)
    const dAlpha_dT = -f * kappa / Math.sqrt(T_K * Tc_K);
    return { alpha, dAlpha_dT };
  }
  // Boston-Mathias extrapolation for T > Tc
  const c2 = 1 + kappa / 2;
  const c1 = kappa / c2;
  const alpha = Math.exp(c1 * (1 - Math.pow(Tr, c2)));
  const dAlpha_dT = alpha * (-c1 * c2 * Math.pow(Tr, c2 - 1) / Tc_K);
  return { alpha, dAlpha_dT };
}

/**
 * Compound-aware alpha selector for PR EOS.
 * Uses Mathias-Copeman if MC parameters are available, otherwise standard Soave/BM.
 */
export function alphaForCompound(
  comp: CanopyCompound, T_K: number, eos: 'PR' | 'SRK' = 'PR',
): { alpha: number; dAlpha_dT: number } {
  if (comp.mathiasCopeman) {
    const { C1, C2, C3 } = comp.mathiasCopeman;
    return alphaMathiasCopeman(comp.Tc_K, T_K, C1, C2, C3);
  }
  if (comp.twuAlpha) {
    const { L, M, N } = comp.twuAlpha;
    return alphaTwu(comp.Tc_K, T_K, L, M, N);
  }
  return eos === 'PR'
    ? alphaPR(comp.Tc_K, comp.omega, T_K)
    : alphaSRK(comp.Tc_K, comp.omega, T_K);
}

export function alphaSRK(Tc_K: number, omega: number, T_K: number): { alpha: number; dAlpha_dT: number } {
  const kappa = 0.48 + 1.574 * omega - 0.176 * omega * omega;
  const Tr = T_K / Tc_K;
  if (Tr <= 1) {
    const sqrtTr = Math.sqrt(Tr);
    const f = 1 + kappa * (1 - sqrtTr);
    const alpha = f * f;
    const dAlpha_dT = -f * kappa / Math.sqrt(T_K * Tc_K);
    return { alpha, dAlpha_dT };
  }
  const c2 = 1 + kappa / 2;
  const c1 = kappa / c2;
  const alpha = Math.exp(c1 * (1 - Math.pow(Tr, c2)));
  const dAlpha_dT = alpha * (-c1 * c2 * Math.pow(Tr, c2 - 1) / Tc_K);
  return { alpha, dAlpha_dT };
}

/**
 * Mathias-Copeman alpha function for PR EOS (polar compounds).
 *
 * Below Tc: α = [1 + c1·(1-√Tr) + c2·(1-√Tr)² + c3·(1-√Tr)³]²
 * Above Tc: α = [1 + c1·(1-√Tr)]²  (only first term)
 *
 * @param c1 Mathias-Copeman parameter 1
 * @param c2 Mathias-Copeman parameter 2
 * @param c3 Mathias-Copeman parameter 3
 */
export function alphaMathiasCopeman(
  Tc_K: number, T_K: number, c1: number, c2: number, c3: number,
): { alpha: number; dAlpha_dT: number } {
  const Tr = T_K / Tc_K;
  const sqrtTr = Math.sqrt(Tr);
  const tau = 1 - sqrtTr;

  if (Tr <= 1) {
    const f = 1 + c1 * tau + c2 * tau * tau + c3 * tau * tau * tau;
    const alpha = f * f;
    // df/dT = (c1 + 2*c2*tau + 3*c3*tau²) · d(tau)/dT
    // d(tau)/dT = -1/(2*sqrt(T_K*Tc_K))
    const dtau_dT = -1 / (2 * Math.sqrt(T_K * Tc_K));
    const df_dT = (c1 + 2 * c2 * tau + 3 * c3 * tau * tau) * dtau_dT;
    return { alpha, dAlpha_dT: 2 * f * df_dT };
  }
  // Above Tc: use only c1
  const f = 1 + c1 * tau;
  const alpha = f * f;
  const dtau_dT = -1 / (2 * Math.sqrt(T_K * Tc_K));
  return { alpha, dAlpha_dT: 2 * f * c1 * dtau_dT };
}

/**
 * Twu alpha function (Twu, Coon, Cunningham, 1995).
 *
 * α = Tr^(N·(M-1)) · exp[L·(1-Tr^(N·M))]
 *
 * This generalized form works for both standard and polar compounds.
 * Aspen Plus allows user-specified L, M, N parameters.
 * Default values for PR: L=0.1253, M=0.8560, N=2.0 (nonpolar)
 */
export function alphaTwu(
  Tc_K: number, T_K: number, L: number, M: number, N: number,
): { alpha: number; dAlpha_dT: number } {
  const Tr = T_K / Tc_K;
  const NM = N * M;
  const Tr_NM_1 = Math.pow(Tr, N * (M - 1));
  const Tr_NM = Math.pow(Tr, NM);
  const alpha = Tr_NM_1 * Math.exp(L * (1 - Tr_NM));

  // dα/dT = α · [N*(M-1)/T_K - L*N*M*Tr^(NM-1)/Tc_K]
  const dAlpha_dT = alpha * (N * (M - 1) / T_K - L * NM * Math.pow(Tr, NM - 1) / Tc_K);
  return { alpha, dAlpha_dT };
}

/**
 * Peng-Robinson K-values: K_i = φ_i^L / φ_i^V
 *
 * a_i = 0.45724 · R²·Tc² / Pc · α(T)
 * b_i = 0.07780 · R·Tc / Pc
 * α(T) = [1 + κ·(1-√(T/Tc))]², κ = 0.37464 + 1.54226ω - 0.26992ω²
 *
 * a_mix = ΣΣ xi·xj·(1-kij)·√(ai·aj)
 * b_mix = Σ xi·bi
 *
 * Z³ - (1-B)Z² + (A-3B²-2B)Z - (AB-B²-B³) = 0
 * A = a_mix·P/(R·T)², B = b_mix·P/(R·T)
 *
 * ln(φ_i) = (bi/b_mix)(Z-1) - ln(Z-B) - A/(2√2·B) · [2Σ_j xj·a_ij/(a_mix) - bi/b_mix] · ln[(Z+2.414B)/(Z-0.414B)]
 */
export function KvaluesPR(
  compounds: CanopyCompound[],
  x_liq: number[],
  y_vap: number[],
  T_K: number,
  P_Pa: number,
  kij?: number[][],
): number[] {
  const n = compounds.length;

  // Component parameters
  const a_i = new Array(n).fill(0);
  const b_i = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const { Tc_K, Pc_Pa } = compounds[i];
    const { alpha } = alphaForCompound(compounds[i], T_K, 'PR');
    a_i[i] = 0.45724 * R * R * Tc_K * Tc_K / Pc_Pa * alpha;
    b_i[i] = 0.07780 * R * Tc_K / Pc_Pa;
  }

  // Solve for both phases
  const solvePR = (z_comp: number[], isVapor: boolean): { Z: number; lnPhi: number[] } => {
    // Mixing rules
    let a_mix = 0, b_mix = 0;
    const a_ij: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      b_mix += z_comp[i] * b_i[i];
      for (let j = 0; j < n; j++) {
        const k = kij?.[i]?.[j] ?? 0;
        a_ij[i][j] = (1 - k) * Math.sqrt(a_i[i] * a_i[j]);
        a_mix += z_comp[i] * z_comp[j] * a_ij[i][j];
      }
    }

    const A_ = a_mix * P_Pa / (R * R * T_K * T_K);
    const B_ = b_mix * P_Pa / (R * T_K);

    // Solve cubic: Z³ - (1-B)Z² + (A-3B²-2B)Z - (AB-B²-B³) = 0
    const c2 = -(1 - B_);
    const c1 = A_ - 3 * B_ * B_ - 2 * B_;
    const c0 = -(A_ * B_ - B_ * B_ - B_ * B_ * B_);

    const roots = solveCubicReal(1, c2, c1, c0);
    // Pick appropriate root: largest for vapor, smallest positive for liquid
    const validRoots = roots.filter(r => r > B_);
    const Z = validRoots.length > 0
      ? (isVapor ? Math.max(...validRoots) : Math.min(...validRoots))
      : roots.length > 0 ? Math.max(...roots) : 1;

    // Fugacity coefficients
    const sqrt2 = Math.SQRT2;
    const lnPhi = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      const sum_a = compounds.reduce((s, _, j) => s + z_comp[j] * a_ij[i][j], 0);
      const biOverBm = b_i[i] / (b_mix || 1e-30);
      const arg1 = Z + (1 + sqrt2) * B_;
      const arg2 = Z - (sqrt2 - 1) * B_;
      const logTerm = (arg1 > 0 && arg2 > 0) ? Math.log(arg1 / arg2) : 0;
      lnPhi[i] = biOverBm * (Z - 1) - Math.log(Math.max(Z - B_, 1e-30))
        - A_ / (2 * sqrt2 * (B_ || 1e-30))
        * (2 * sum_a / (a_mix || 1e-30) - biOverBm) * logTerm;
    }
    return { Z, lnPhi };
  };

  const liq = solvePR(x_liq, false);
  const vap = solvePR(y_vap, true);

  // K_i = φ_i^L / φ_i^V
  return compounds.map((_, i) => {
    const ratio = liq.lnPhi[i] - vap.lnPhi[i];
    return Math.exp(Math.max(-50, Math.min(50, ratio)));
  });
}

/** Solve real roots of cubic: a·x³ + b·x² + c·x + d = 0 */
function solveCubicReal(a: number, b: number, c: number, d: number): number[] {
  const p = (3 * a * c - b * b) / (3 * a * a);
  const q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (27 * a * a * a);
  const disc = q * q / 4 + p * p * p / 27;

  const shift = -b / (3 * a);

  if (disc > 1e-15) {
    // One real root
    const sqrtDisc = Math.sqrt(disc);
    const u = Math.cbrt(-q / 2 + sqrtDisc);
    const v = Math.cbrt(-q / 2 - sqrtDisc);
    return [u + v + shift];
  } else if (disc < -1e-15) {
    // Three real roots
    const r = Math.sqrt(-p * p * p / 27);
    const theta = Math.acos(Math.max(-1, Math.min(1, -q / (2 * r))));
    const m = 2 * Math.cbrt(r);
    return [
      m * Math.cos(theta / 3) + shift,
      m * Math.cos((theta + 2 * Math.PI) / 3) + shift,
      m * Math.cos((theta + 4 * Math.PI) / 3) + shift,
    ];
  } else {
    // Repeated root
    const u = Math.cbrt(-q / 2);
    return [2 * u + shift, -u + shift];
  }
}

/**
 * Unified K-value computation dispatching on fluid package.
 * Falls back to Ideal if non-ideal parameters are not provided.
 */
export type FluidPackageType =
  | 'Ideal'
  | 'NRTL'
  | 'Electrolyte-NRTL'
  | 'Wilson'
  | 'UNIQUAC'
  | 'UNIFAC'
  | 'UNIFAC-DMD'
  | 'Peng-Robinson'
  | 'SRK';

/** Wilson binary interaction parameter set.
 *  Λ_ij = (V_j/V_i) · exp(-a_ij / (R·T))  where a_ij are stored in the A matrix.
 *  A[i][j] = a_ij in cal/mol (energy parameter).
 *  molarVolumes[i] = V_i in cm³/mol (liquid molar volumes).
 *
 *  Aspen 12-element extended form:
 *    Λ_ij = (Vj/Vi)·exp(a_ij + b_ij/T + c_ij·ln(T) + d_ij·T + e_ij·T²)
 */
export interface WilsonParams {
  /** Energy parameter a_ij (cal/mol) — simple form */
  A: number[][];
  /** Liquid molar volumes V_i (cm³/mol) */
  molarVolumes: number[];
  /** Extended: a_ij (dimensionless) — ln(Λ) intercept */
  a_ext?: number[][];
  /** Extended: b_ij (K) — 1/T coefficient */
  b_ext?: number[][];
  /** Extended: c_ij — ln(T) coefficient */
  c_ext?: number[][];
  /** Extended: d_ij (1/K) — T coefficient */
  d_ext?: number[][];
  /** Extended: e_ij (1/K²) — T² coefficient */
  e_ext?: number[][];
}

/**
 * Multi-component Wilson activity coefficients.
 *
 * Λ_ij = (V_j/V_i) · exp(-(a_ij) / (R_cal·T))
 *
 * ln(γ_i) = -ln(Σ_j x_j Λ_ij) + 1 - Σ_k [x_k Λ_ki / Σ_j x_j Λ_kj]
 *
 * Reference: Wilson, G.M. (1964) J. Am. Chem. Soc. 86:127.
 */
export function wilsonGamma(
  x: number[],
  T_K: number,
  params: WilsonParams,
): number[] {
  const n = x.length;
  const R_cal = 1.98721; // cal/(mol·K)

  // Build Λ matrix
  const Lambda: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
  const hasExtended = !!(params.a_ext || params.b_ext || params.c_ext || params.d_ext || params.e_ext);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i === j) {
        Lambda[i][j] = 1;
      } else {
        const Vi = params.molarVolumes[i] || 1;
        const Vj = params.molarVolumes[j] || 1;
        if (hasExtended) {
          // Aspen extended: Λ_ij = (Vj/Vi)·exp(a + b/T + c·ln(T) + d·T + e·T²)
          const a_val = params.a_ext?.[i]?.[j] ?? 0;
          const b_val = params.b_ext?.[i]?.[j] ?? 0;
          const c_val = params.c_ext?.[i]?.[j] ?? 0;
          const d_val = params.d_ext?.[i]?.[j] ?? 0;
          const e_val = params.e_ext?.[i]?.[j] ?? 0;
          const lnLam = a_val + b_val / T_K + c_val * Math.log(T_K)
            + d_val * T_K + e_val * T_K * T_K;
          Lambda[i][j] = (Vj / Vi) * Math.exp(lnLam);
        } else {
          Lambda[i][j] = (Vj / Vi) * Math.exp(-params.A[i][j] / (R_cal * T_K));
        }
      }
    }
  }

  // Precompute Σ_j x_j Λ_ij for each i
  const sumXLambda: number[] = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      sumXLambda[i] += x[j] * Lambda[i][j];
    }
  }

  const gamma = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    // Term 1: -ln(Σ_j x_j Λ_ij)
    const term1 = -Math.log(Math.max(sumXLambda[i], 1e-30));

    // Term 2: 1 - Σ_k [x_k Λ_ki / Σ_j x_j Λ_kj]
    let term2 = 1;
    for (let k = 0; k < n; k++) {
      const denom = sumXLambda[k] || 1e-30;
      term2 -= (x[k] * Lambda[k][i]) / denom;
    }

    const lnGamma = term1 + term2;
    gamma[i] = Math.exp(Math.max(-50, Math.min(50, lnGamma)));
  }

  return gamma;
}

// ────────────────────────────────────────────────────────────────
// 6c. UNIQUAC Activity Coefficient Model
// ────────────────────────────────────────────────────────────────

/** UNIQUAC binary interaction parameter set.
 *  Simple form:  τ_ij = exp(-Δu_ij / (R·T))  where Δu_ij = a_ij + b_ij·T
 *  Aspen 12-element extended form:
 *    τ_ij = exp(a_ij + b_ij/T + c_ij·ln(T) + d_ij·T + e_ij·T²)
 */
export interface UNIQUACParams {
  /** a_ij energy interaction parameter (cal/mol) — n×n matrix */
  a: number[][];
  /** b_ij temperature-dependent parameter (cal/mol/K) — n×n matrix */
  b: number[][];
  /** r_i — van der Waals volume parameter per component */
  r: number[];
  /** q_i — van der Waals surface area parameter per component */
  q: number[];
  /** Extended: a_ij (dimensionless) — ln(τ) intercept */
  a_ext?: number[][];
  /** Extended: b_ij (K) — 1/T coefficient in ln(τ) */
  b_ext?: number[][];
  /** Extended: c_ij — ln(T) coefficient in ln(τ) */
  c_ext?: number[][];
  /** Extended: d_ij (1/K) — T coefficient in ln(τ) */
  d_ext?: number[][];
  /** Extended: e_ij (1/K²) — T² coefficient in ln(τ) */
  e_ext?: number[][];
}

/**
 * Multi-component UNIQUAC activity coefficients.
 *
 * Combinatorial part (entropic, size/shape):
 *   ln(γ_i^C) = ln(Φ_i/x_i) + (z/2)·q_i·ln(θ_i/Φ_i) + l_i - (Φ_i/x_i)·Σ_j x_j·l_j
 *   where Φ_i = x_i·r_i / Σ_j x_j·r_j  (volume fraction)
 *         θ_i = x_i·q_i / Σ_j x_j·q_j  (surface fraction)
 *         l_i = (z/2)·(r_i - q_i) - (r_i - 1)  with z = 10
 *
 * Residual part (enthalpic, interactions):
 *   ln(γ_i^R) = q_i · [1 - ln(Σ_j θ_j·τ_ji) - Σ_j (θ_j·τ_ij / Σ_k θ_k·τ_kj)]
 *   where τ_ij = exp(-(a_ij + b_ij·T) / (R_cal·T))
 *
 * Reference: Abrams & Prausnitz, AIChE J. 21:116 (1975).
 */
export function uniquacGamma(
  x: number[],
  T_K: number,
  params: UNIQUACParams,
): number[] {
  const n = x.length;
  const R_cal = 1.98721; // cal/(mol·K)
  const z = 10; // coordination number

  const { r, q, a, b } = params;

  // Volume and surface fractions
  const sumXR = x.reduce((s, xi, i) => s + xi * r[i], 0) || 1e-30;
  const sumXQ = x.reduce((s, xi, i) => s + xi * q[i], 0) || 1e-30;
  const Phi = x.map((xi, i) => xi * r[i] / sumXR);  // volume fraction
  const theta = x.map((xi, i) => xi * q[i] / sumXQ); // surface fraction

  // l_i = (z/2)(r_i - q_i) - (r_i - 1)
  const l = r.map((ri, i) => (z / 2) * (ri - q[i]) - (ri - 1));

  // τ_ij: simple or extended form
  const tau: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
  const hasExtended = !!(params.a_ext || params.b_ext || params.c_ext || params.d_ext || params.e_ext);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i === j) {
        tau[i][j] = 1;
      } else if (hasExtended) {
        // Aspen extended: τ_ij = exp(a + b/T + c·ln(T) + d·T + e·T²)
        const a_val = params.a_ext?.[i]?.[j] ?? 0;
        const b_val = params.b_ext?.[i]?.[j] ?? 0;
        const c_val = params.c_ext?.[i]?.[j] ?? 0;
        const d_val = params.d_ext?.[i]?.[j] ?? 0;
        const e_val = params.e_ext?.[i]?.[j] ?? 0;
        const lnTau = a_val + b_val / T_K + c_val * Math.log(T_K)
          + d_val * T_K + e_val * T_K * T_K;
        tau[i][j] = Math.exp(Math.max(-50, Math.min(50, lnTau)));
      } else {
        const du = a[i][j] + (b[i]?.[j] ?? 0) * T_K;
        tau[i][j] = Math.exp(-du / (R_cal * T_K));
      }
    }
  }

  // Precompute Σ_j θ_j·τ_ji for each i  and  Σ_k θ_k·τ_kj for each j
  const sumThetaTau_ji: number[] = new Array(n).fill(0); // for residual term1
  const sumThetaTau_kj: number[] = new Array(n).fill(0); // for residual term2 denominator
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      sumThetaTau_ji[i] += theta[j] * tau[j][i];
      sumThetaTau_kj[j] += theta[i] * tau[i][j]; // same sum, indexed by j
    }
  }
  // sumThetaTau_kj computed above is Σ_k θ_k·τ_kj for each j — recalculate properly
  const S_j = new Array(n).fill(0); // S_j = Σ_k θ_k · τ_kj
  for (let j = 0; j < n; j++) {
    for (let k = 0; k < n; k++) {
      S_j[j] += theta[k] * tau[k][j];
    }
  }

  const sumXL = x.reduce((s, xi, i) => s + xi * l[i], 0);

  const gamma = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    // Combinatorial
    const phiOverX = Phi[i] / (x[i] || 1e-30);
    const thetaOverPhi = theta[i] / (Phi[i] || 1e-30);
    const lnGammaC = Math.log(Math.max(phiOverX, 1e-30))
      + (z / 2) * q[i] * Math.log(Math.max(thetaOverPhi, 1e-30))
      + l[i] - phiOverX * sumXL;

    // Residual
    // ln(γ_i^R) = q_i * [1 - ln(Σ_j θ_j·τ_ji) - Σ_j(θ_j·τ_ij / S_j)]
    const term1_R = -Math.log(Math.max(sumThetaTau_ji[i], 1e-30));
    let term2_R = 0;
    for (let j = 0; j < n; j++) {
      term2_R += theta[j] * tau[i][j] / (S_j[j] || 1e-30);
    }
    const lnGammaR = q[i] * (1 + term1_R - term2_R);

    const lnG = lnGammaC + lnGammaR;
    gamma[i] = Math.exp(Math.max(-50, Math.min(50, lnG)));
  }

  return gamma;
}

// ────────────────────────────────────────────────────────────────
// 6d. UNIFAC Group Contribution Activity Coefficient Model
// ────────────────────────────────────────────────────────────────

/** UNIFAC subgroup definition */
export interface UNIFACSubgroup {
  subgroupId: number;
  mainGroupId: number;
  R: number;  // van der Waals volume
  Q: number;  // van der Waals surface area
}

/** UNIFAC interaction parameter a_mn (K) between main groups m and n.
 *  τ_mn = exp(-a_mn / T)  — note: asymmetric a_mn ≠ a_nm
 */
export interface UNIFACData {
  subgroups: UNIFACSubgroup[];
  /** Interaction parameters: key is "m-n" → a_mn (K) */
  interactions: Record<string, number>;
}

/**
 * UNIFAC activity coefficients (original Fredenslund et al., 1975).
 *
 * Combinatorial part (same as UNIQUAC):
 *   ln(γ_i^C) = ln(Φ_i/x_i) + (z/2)·q_i·ln(θ_i/Φ_i) + l_i - (Φ_i/x_i)·Σ x_j·l_j
 *   r_i = Σ_k ν_k^i · R_k,  q_i = Σ_k ν_k^i · Q_k
 *
 * Residual part (group contribution):
 *   ln(γ_i^R) = Σ_k ν_k^i · [ln(Γ_k) - ln(Γ_k^i)]
 *   where Γ_k is the group activity coefficient in the mixture
 *   and Γ_k^i is the group activity coefficient in pure component i
 *
 *   ln(Γ_k) = Q_k · [1 - ln(Σ_m Θ_m·ψ_mk) - Σ_m (Θ_m·ψ_km / Σ_n Θ_n·ψ_nm)]
 *   ψ_mn = exp(-a_mn / T)
 *   Θ_m = Q_m·X_m / Σ_n Q_n·X_n   (group surface fraction)
 *   X_m = (Σ_j x_j·ν_m^j) / (Σ_j x_j·Σ_k ν_k^j)  (group mole fraction)
 *
 * @param x Mole fractions of components
 * @param T_K Temperature in K
 * @param compGroups Array of {subgroupId → count} for each component
 * @param unifacData UNIFAC subgroup definitions and interaction parameters
 */
export function unifacGamma(
  x: number[],
  T_K: number,
  compGroups: Record<number, number>[],
  unifacData: UNIFACData,
): number[] {
  const nc = x.length;
  const z = 10;

  // Build lookup: subgroupId → {R, Q, mainGroupId}
  const sgLookup = new Map<number, UNIFACSubgroup>();
  for (const sg of unifacData.subgroups) sgLookup.set(sg.subgroupId, sg);

  // Compute r_i, q_i for each component
  const r_i = new Array(nc).fill(0);
  const q_i = new Array(nc).fill(0);
  for (let i = 0; i < nc; i++) {
    for (const [sgIdStr, count] of Object.entries(compGroups[i])) {
      const sg = sgLookup.get(Number(sgIdStr));
      if (!sg) continue;
      r_i[i] += count * sg.R;
      q_i[i] += count * sg.Q;
    }
  }

  // Collect all unique subgroup IDs present in the mixture
  const allSgIds = new Set<number>();
  for (const groups of compGroups) {
    for (const sgId of Object.keys(groups)) allSgIds.add(Number(sgId));
  }
  const sgIds = Array.from(allSgIds).sort((a, b) => a - b);

  // ψ_mn = exp(-a_mn/T) between main groups
  const psi = (mainM: number, mainN: number): number => {
    if (mainM === mainN) return 1;
    const key = `${mainM}-${mainN}`;
    const a_mn = unifacData.interactions[key];
    if (a_mn === undefined) return 1; // missing → no interaction
    return Math.exp(-a_mn / T_K);
  };

  // ── Combinatorial Part ──
  const sumXR = r_i.reduce((s, ri, i) => s + x[i] * ri, 0) || 1e-30;
  const sumXQ = q_i.reduce((s, qi, i) => s + x[i] * qi, 0) || 1e-30;
  const Phi = x.map((xi, i) => xi * r_i[i] / sumXR);
  const theta = x.map((xi, i) => xi * q_i[i] / sumXQ);
  const l_i = r_i.map((ri, i) => (z / 2) * (ri - q_i[i]) - (ri - 1));
  const sumXL = x.reduce((s, xi, i) => s + xi * l_i[i], 0);

  const lnGammaC = new Array(nc).fill(0);
  for (let i = 0; i < nc; i++) {
    const phiOverX = Phi[i] / (x[i] || 1e-30);
    const tOverP = theta[i] / (Phi[i] || 1e-30);
    lnGammaC[i] = Math.log(Math.max(phiOverX, 1e-30))
      + (z / 2) * q_i[i] * Math.log(Math.max(tOverP, 1e-30))
      + l_i[i] - phiOverX * sumXL;
  }

  // ── Residual Part ──
  // Helper: compute ln(Γ_k) for each subgroup given group mole fractions
  const computeGroupGamma = (X_sg: Map<number, number>): Map<number, number> => {
    // Group surface fractions: Θ_m = Q_m·X_m / Σ Q_n·X_n
    let sumQX = 0;
    for (const [sgId, Xm] of X_sg) {
      const sg = sgLookup.get(sgId)!;
      sumQX += sg.Q * Xm;
    }
    sumQX = sumQX || 1e-30;

    const Theta = new Map<number, number>();
    for (const [sgId, Xm] of X_sg) {
      const sg = sgLookup.get(sgId)!;
      Theta.set(sgId, sg.Q * Xm / sumQX);
    }

    // For each subgroup k: ln(Γ_k) = Q_k · [1 - ln(Σ_m Θ_m·ψ_mk) - Σ_m(Θ_m·ψ_km / Σ_n Θ_n·ψ_nm)]
    const lnGamma_k = new Map<number, number>();
    for (const k of X_sg.keys()) {
      const sgK = sgLookup.get(k)!;
      const mainK = sgK.mainGroupId;

      // Σ_m Θ_m · ψ_mk
      let sumMK = 0;
      for (const [m, theta_m] of Theta) {
        const mainM = sgLookup.get(m)!.mainGroupId;
        sumMK += theta_m * psi(mainM, mainK);
      }

      // Σ_m (Θ_m · ψ_km / Σ_n Θ_n · ψ_nm)
      let sumTerm2 = 0;
      for (const [m, theta_m] of Theta) {
        const mainM = sgLookup.get(m)!.mainGroupId;
        // Σ_n Θ_n · ψ_nm
        let sumNM = 0;
        for (const [nn, theta_n] of Theta) {
          const mainN = sgLookup.get(nn)!.mainGroupId;
          sumNM += theta_n * psi(mainN, mainM);
        }
        sumTerm2 += theta_m * psi(mainK, mainM) / (sumNM || 1e-30);
      }

      lnGamma_k.set(k, sgK.Q * (1 - Math.log(Math.max(sumMK, 1e-30)) - sumTerm2));
    }
    return lnGamma_k;
  };

  // Group mole fractions in mixture: X_m = Σ_j x_j·ν_m^j / Σ_j x_j·Σ_k ν_k^j
  let totalGroupMoles = 0;
  const groupMoles = new Map<number, number>();
  for (let j = 0; j < nc; j++) {
    for (const [sgIdStr, count] of Object.entries(compGroups[j])) {
      const sgId = Number(sgIdStr);
      groupMoles.set(sgId, (groupMoles.get(sgId) ?? 0) + x[j] * count);
      totalGroupMoles += x[j] * count;
    }
  }
  const X_mix = new Map<number, number>();
  for (const [sgId, moles] of groupMoles) {
    X_mix.set(sgId, moles / (totalGroupMoles || 1e-30));
  }

  const lnGamma_mix = computeGroupGamma(X_mix);

  // Group activity coefficients in pure component i
  const lnGamma_pure: Map<number, number>[] = [];
  for (let i = 0; i < nc; i++) {
    const X_pure_i = new Map<number, number>();
    let totalI = 0;
    for (const [sgIdStr, count] of Object.entries(compGroups[i])) {
      totalI += count;
    }
    for (const [sgIdStr, count] of Object.entries(compGroups[i])) {
      X_pure_i.set(Number(sgIdStr), count / (totalI || 1));
    }
    lnGamma_pure.push(computeGroupGamma(X_pure_i));
  }

  // ln(γ_i^R) = Σ_k ν_k^i · [ln(Γ_k) - ln(Γ_k^i)]
  const lnGammaR = new Array(nc).fill(0);
  for (let i = 0; i < nc; i++) {
    for (const [sgIdStr, count] of Object.entries(compGroups[i])) {
      const sgId = Number(sgIdStr);
      const lnGk_mix = lnGamma_mix.get(sgId) ?? 0;
      const lnGk_pure = lnGamma_pure[i].get(sgId) ?? 0;
      lnGammaR[i] += count * (lnGk_mix - lnGk_pure);
    }
  }

  // Total: γ_i = γ_i^C · γ_i^R
  const gamma = new Array(nc).fill(0);
  for (let i = 0; i < nc; i++) {
    const lnG = lnGammaC[i] + lnGammaR[i];
    gamma[i] = Math.exp(Math.max(-50, Math.min(50, lnG)));
  }

  return gamma;
}

// ────────────────────────────────────────────────────────────────
// UNIFAC-DMD (Dortmund Modified UNIFAC)
// ────────────────────────────────────────────────────────────────

/** UNIFAC-DMD interaction parameters: ψ_mn = exp[-(a_mn + b_mn·T + c_mn·T²)/T] */
export interface UNIFACDMDData {
  subgroups: UNIFACSubgroup[];
  /** Key is "m-n" → { a, b, c } */
  interactions: Record<string, { a: number; b: number; c: number }>;
}

/**
 * UNIFAC-DMD (Dortmund Modified) activity coefficients.
 *
 * Differences from original UNIFAC:
 *  1. Modified combinatorial: Φ'_i = x_i·r_i^(3/4) / Σ x_j·r_j^(3/4)
 *     ln(γ_i^C) = 1 - Φ'_i/x_i + ln(Φ'_i/x_i) - 5·q_i·(1 - Φ_i/(θ_i) + ln(Φ_i/θ_i))
 *     where Φ_i uses r_i (not modified), θ_i uses q_i (standard)
 *  2. Temperature-dependent interactions: ψ_mn = exp[-(a_mn + b_mn·T + c_mn·T²)/T]
 */
export function unifacDMDGamma(
  x: number[],
  T_K: number,
  compGroups: Record<number, number>[],
  dmdData: UNIFACDMDData,
): number[] {
  const nc = x.length;
  const z = 10;

  const sgLookup = new Map<number, UNIFACSubgroup>();
  for (const sg of dmdData.subgroups) sgLookup.set(sg.subgroupId, sg);

  const r_i = new Array(nc).fill(0);
  const q_i = new Array(nc).fill(0);
  for (let i = 0; i < nc; i++) {
    for (const [sgIdStr, count] of Object.entries(compGroups[i])) {
      const sg = sgLookup.get(Number(sgIdStr));
      if (!sg) continue;
      r_i[i] += count * sg.R;
      q_i[i] += count * sg.Q;
    }
  }

  // ψ_mn with T-dependent parameters: exp[-(a + b·T + c·T²)/T]
  const psi = (mainM: number, mainN: number): number => {
    if (mainM === mainN) return 1;
    const key = `${mainM}-${mainN}`;
    const p = dmdData.interactions[key];
    if (!p) return 1;
    return Math.exp(-(p.a + p.b * T_K + p.c * T_K * T_K) / T_K);
  };

  // ── Modified Combinatorial Part (Dortmund) ──
  const sumXR = r_i.reduce((s, ri, i) => s + x[i] * ri, 0) || 1e-30;
  const sumXQ = q_i.reduce((s, qi, i) => s + x[i] * qi, 0) || 1e-30;
  const sumXR34 = r_i.reduce((s, ri, i) => s + x[i] * Math.pow(ri, 0.75), 0) || 1e-30;

  // Standard volume/surface fractions
  const Phi = x.map((xi, i) => xi * r_i[i] / sumXR);
  const theta = x.map((xi, i) => xi * q_i[i] / sumXQ);
  // Modified volume fraction: r_i^(3/4)
  const PhiPrime = x.map((xi, i) => xi * Math.pow(r_i[i], 0.75) / sumXR34);

  const lnGammaC = new Array(nc).fill(0);
  for (let i = 0; i < nc; i++) {
    const xi = x[i] || 1e-30;
    const ppx = PhiPrime[i] / xi;
    const phiT = Phi[i] / (theta[i] || 1e-30);
    lnGammaC[i] = 1 - ppx + Math.log(Math.max(ppx, 1e-30))
      - 5 * q_i[i] * (1 - phiT + Math.log(Math.max(phiT, 1e-30)));
  }

  // ── Residual Part (same form as original, but with T-dependent ψ) ──
  const computeGroupGamma = (X_sg: Map<number, number>): Map<number, number> => {
    let sumQX = 0;
    for (const [sgId, Xm] of X_sg) {
      sumQX += (sgLookup.get(sgId)?.Q ?? 0) * Xm;
    }
    sumQX = sumQX || 1e-30;

    const Theta = new Map<number, number>();
    for (const [sgId, Xm] of X_sg) {
      Theta.set(sgId, (sgLookup.get(sgId)?.Q ?? 0) * Xm / sumQX);
    }

    const lnGamma_k = new Map<number, number>();
    for (const k of X_sg.keys()) {
      const sgK = sgLookup.get(k)!;
      const mainK = sgK.mainGroupId;

      let sumMK = 0;
      for (const [m, theta_m] of Theta) {
        sumMK += theta_m * psi(sgLookup.get(m)!.mainGroupId, mainK);
      }

      let sumTerm2 = 0;
      for (const [m, theta_m] of Theta) {
        const mainM = sgLookup.get(m)!.mainGroupId;
        let sumNM = 0;
        for (const [nn, theta_n] of Theta) {
          sumNM += theta_n * psi(sgLookup.get(nn)!.mainGroupId, mainM);
        }
        sumTerm2 += theta_m * psi(mainK, mainM) / (sumNM || 1e-30);
      }

      lnGamma_k.set(k, sgK.Q * (1 - Math.log(Math.max(sumMK, 1e-30)) - sumTerm2));
    }
    return lnGamma_k;
  };

  // Group mole fractions in mixture
  let totalGroupMoles = 0;
  const groupMoles = new Map<number, number>();
  for (let j = 0; j < nc; j++) {
    for (const [sgIdStr, count] of Object.entries(compGroups[j])) {
      const sgId = Number(sgIdStr);
      groupMoles.set(sgId, (groupMoles.get(sgId) ?? 0) + x[j] * count);
      totalGroupMoles += x[j] * count;
    }
  }
  const X_mix = new Map<number, number>();
  for (const [sgId, moles] of groupMoles) {
    X_mix.set(sgId, moles / (totalGroupMoles || 1e-30));
  }

  const lnGamma_mix = computeGroupGamma(X_mix);

  const lnGamma_pure: Map<number, number>[] = [];
  for (let i = 0; i < nc; i++) {
    const X_pure_i = new Map<number, number>();
    let totalI = 0;
    for (const [, count] of Object.entries(compGroups[i])) totalI += count;
    for (const [sgIdStr, count] of Object.entries(compGroups[i])) {
      X_pure_i.set(Number(sgIdStr), count / (totalI || 1));
    }
    lnGamma_pure.push(computeGroupGamma(X_pure_i));
  }

  const lnGammaR = new Array(nc).fill(0);
  for (let i = 0; i < nc; i++) {
    for (const [sgIdStr, count] of Object.entries(compGroups[i])) {
      const sgId = Number(sgIdStr);
      lnGammaR[i] += count * ((lnGamma_mix.get(sgId) ?? 0) - (lnGamma_pure[i].get(sgId) ?? 0));
    }
  }

  const gamma = new Array(nc).fill(0);
  for (let i = 0; i < nc; i++) {
    const lnG = lnGammaC[i] + lnGammaR[i];
    gamma[i] = Math.exp(Math.max(-50, Math.min(50, lnG)));
  }

  return gamma;
}

/**
 * SRK (Soave-Redlich-Kwong) K-values: K_i = φ_i^L / φ_i^V
 *
 * a_i = 0.42748 · R²·Tc² / Pc · α(T)
 * b_i = 0.08664 · R·Tc / Pc
 * α(T) = [1 + κ·(1-√(T/Tc))]², κ = 0.48 + 1.574ω - 0.176ω²
 *
 * a_mix = ΣΣ xi·xj·(1-kij)·√(ai·aj)
 * b_mix = Σ xi·bi
 *
 * Z³ - Z² + (A-B-B²)Z - AB = 0
 * A = a_mix·P/(R·T)², B = b_mix·P/(R·T)
 *
 * ln(φ_i) = (bi/b_mix)(Z-1) - ln(Z-B) - A/B · [2Σ_j xj·a_ij/a_mix - bi/b_mix] · ln(1+B/Z)
 */
export function KvaluesSRK(
  compounds: CanopyCompound[],
  x_liq: number[],
  y_vap: number[],
  T_K: number,
  P_Pa: number,
  kij?: number[][],
): number[] {
  const n = compounds.length;

  // Component parameters — SRK-specific constants
  const a_i = new Array(n).fill(0);
  const b_i = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const { Tc_K, Pc_Pa } = compounds[i];
    const { alpha } = alphaForCompound(compounds[i], T_K, 'SRK');
    a_i[i] = 0.42748 * R * R * Tc_K * Tc_K / Pc_Pa * alpha;
    b_i[i] = 0.08664 * R * Tc_K / Pc_Pa;
  }

  const solveSRK = (z_comp: number[], isVapor: boolean): { Z: number; lnPhi: number[] } => {
    let a_mix = 0, b_mix = 0;
    const a_ij: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      b_mix += z_comp[i] * b_i[i];
      for (let j = 0; j < n; j++) {
        const k = kij?.[i]?.[j] ?? 0;
        a_ij[i][j] = (1 - k) * Math.sqrt(a_i[i] * a_i[j]);
        a_mix += z_comp[i] * z_comp[j] * a_ij[i][j];
      }
    }

    const A_ = a_mix * P_Pa / (R * R * T_K * T_K);
    const B_ = b_mix * P_Pa / (R * T_K);

    // SRK cubic: Z³ - Z² + (A-B-B²)Z - AB = 0
    const c2 = -1;
    const c1 = A_ - B_ - B_ * B_;
    const c0 = -A_ * B_;

    const roots = solveCubicReal(1, c2, c1, c0);
    const validRoots = roots.filter(r => r > B_);
    const Z = validRoots.length > 0
      ? (isVapor ? Math.max(...validRoots) : Math.min(...validRoots))
      : roots.length > 0 ? Math.max(...roots) : 1;

    // Fugacity coefficients — SRK form
    const lnPhi = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      const sum_a = compounds.reduce((s, _, j) => s + z_comp[j] * a_ij[i][j], 0);
      const biOverBm = b_i[i] / (b_mix || 1e-30);
      const arg = 1 + B_ / Z;
      const logTerm = arg > 0 ? Math.log(arg) : 0;
      lnPhi[i] = biOverBm * (Z - 1) - Math.log(Math.max(Z - B_, 1e-30))
        - (A_ / (B_ || 1e-30))
        * (2 * sum_a / (a_mix || 1e-30) - biOverBm) * logTerm;
    }
    return { Z, lnPhi };
  };

  const liq = solveSRK(x_liq, false);
  const vap = solveSRK(y_vap, true);

  return compounds.map((_, i) => {
    const ratio = liq.lnPhi[i] - vap.lnPhi[i];
    return Math.exp(Math.max(-50, Math.min(50, ratio)));
  });
}

// ────────────────────────────────────────────────────────────────
// Fugacity Coefficients (standalone)
// ────────────────────────────────────────────────────────────────

/**
 * Peng-Robinson fugacity coefficients φ_i for a single phase.
 * Returns { Z, phi } where phi[i] is the fugacity coefficient of component i.
 */
export function fugacityPR(
  compounds: CanopyCompound[],
  z_comp: number[],
  T_K: number,
  P_Pa: number,
  isVapor: boolean,
  kij?: number[][],
): { Z: number; phi: number[] } {
  const n = compounds.length;
  const a_i = new Array(n).fill(0);
  const b_i = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const { Tc_K, Pc_Pa } = compounds[i];
    const { alpha } = alphaForCompound(compounds[i], T_K, 'PR');
    a_i[i] = 0.45724 * R * R * Tc_K * Tc_K / Pc_Pa * alpha;
    b_i[i] = 0.07780 * R * Tc_K / Pc_Pa;
  }

  let a_mix = 0, b_mix = 0;
  const a_ij: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    b_mix += z_comp[i] * b_i[i];
    for (let j = 0; j < n; j++) {
      const k = kij?.[i]?.[j] ?? 0;
      a_ij[i][j] = (1 - k) * Math.sqrt(a_i[i] * a_i[j]);
      a_mix += z_comp[i] * z_comp[j] * a_ij[i][j];
    }
  }

  const A_ = a_mix * P_Pa / (R * R * T_K * T_K);
  const B_ = b_mix * P_Pa / (R * T_K);
  const c2 = -(1 - B_);
  const c1 = A_ - 3 * B_ * B_ - 2 * B_;
  const c0 = -(A_ * B_ - B_ * B_ - B_ * B_ * B_);

  const roots = solveCubicReal(1, c2, c1, c0);
  const validRoots = roots.filter(r => r > B_);
  const Z = validRoots.length > 0
    ? (isVapor ? Math.max(...validRoots) : Math.min(...validRoots))
    : roots.length > 0 ? Math.max(...roots) : 1;

  const sqrt2 = Math.SQRT2;
  const phi = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const sum_a = z_comp.reduce((s, zj, j) => s + zj * a_ij[i][j], 0);
    const biOverBm = b_i[i] / (b_mix || 1e-30);
    const arg1 = Z + (1 + sqrt2) * B_;
    const arg2 = Z - (sqrt2 - 1) * B_;
    const logTerm = (arg1 > 0 && arg2 > 0) ? Math.log(arg1 / arg2) : 0;
    const lnPhi = biOverBm * (Z - 1) - Math.log(Math.max(Z - B_, 1e-30))
      - A_ / (2 * sqrt2 * (B_ || 1e-30))
      * (2 * sum_a / (a_mix || 1e-30) - biOverBm) * logTerm;
    phi[i] = Math.exp(Math.max(-50, Math.min(50, lnPhi)));
  }

  return { Z, phi };
}

/**
 * SRK fugacity coefficients φ_i for a single phase.
 * Returns { Z, phi } where phi[i] is the fugacity coefficient of component i.
 */
export function fugacitySRK(
  compounds: CanopyCompound[],
  z_comp: number[],
  T_K: number,
  P_Pa: number,
  isVapor: boolean,
  kij?: number[][],
): { Z: number; phi: number[] } {
  const n = compounds.length;
  const a_i = new Array(n).fill(0);
  const b_i = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const { Tc_K, Pc_Pa } = compounds[i];
    const { alpha } = alphaForCompound(compounds[i], T_K, 'SRK');
    a_i[i] = 0.42748 * R * R * Tc_K * Tc_K / Pc_Pa * alpha;
    b_i[i] = 0.08664 * R * Tc_K / Pc_Pa;
  }

  let a_mix = 0, b_mix = 0;
  const a_ij: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    b_mix += z_comp[i] * b_i[i];
    for (let j = 0; j < n; j++) {
      const k = kij?.[i]?.[j] ?? 0;
      a_ij[i][j] = (1 - k) * Math.sqrt(a_i[i] * a_i[j]);
      a_mix += z_comp[i] * z_comp[j] * a_ij[i][j];
    }
  }

  const A_ = a_mix * P_Pa / (R * R * T_K * T_K);
  const B_ = b_mix * P_Pa / (R * T_K);
  const c2 = -1;
  const c1 = A_ - B_ - B_ * B_;
  const c0 = -(A_ * B_);

  const roots = solveCubicReal(1, c2, c1, c0);
  const validRoots = roots.filter(r => r > B_);
  const Z = validRoots.length > 0
    ? (isVapor ? Math.max(...validRoots) : Math.min(...validRoots))
    : roots.length > 0 ? Math.max(...roots) : 1;

  const phi = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const sum_a = z_comp.reduce((s, zj, j) => s + zj * a_ij[i][j], 0);
    const biOverBm = b_i[i] / (b_mix || 1e-30);
    const lnPhi = biOverBm * (Z - 1) - Math.log(Math.max(Z - B_, 1e-30))
      - A_ / (B_ || 1e-30) * (2 * sum_a / (a_mix || 1e-30) - biOverBm)
      * Math.log(Math.max(1 + B_ / Z, 1e-30));
    phi[i] = Math.exp(Math.max(-50, Math.min(50, lnPhi)));
  }

  return { Z, phi };
}

/**
 * Compressibility factor Z for a mixture using PR or SRK.
 * Convenience wrapper around fugacity functions.
 */
export function compressibilityFactor(
  compounds: CanopyCompound[],
  z_comp: number[],
  T_K: number,
  P_Pa: number,
  isVapor: boolean,
  eos: 'PR' | 'SRK' = 'PR',
  kij?: number[][],
): number {
  const fn = eos === 'PR' ? fugacityPR : fugacitySRK;
  return fn(compounds, z_comp, T_K, P_Pa, isVapor, kij).Z;
}

/**
 * Vapor density (mol/m³) from equation of state.
 * Uses Z factor: ρ = P / (Z·R·T)
 */
export function vaporDensity_molpm3(
  compounds: CanopyCompound[],
  y: number[],
  T_K: number,
  P_Pa: number,
  eos: 'PR' | 'SRK' = 'PR',
  kij?: number[][],
): number {
  const Z = compressibilityFactor(compounds, y, T_K, P_Pa, true, eos, kij);
  return P_Pa / (Z * R * T_K);
}

// ────────────────────────────────────────────────────────────────
// Henry's Law for Dissolved Gases
// ────────────────────────────────────────────────────────────────

/** Henry's law constant coefficients: ln(H) = A + B/T + C·ln(T) + D·T (H in Pa) */
export interface HenryCoeffs {
  A: number;
  B: number;
  C: number;
  D: number;
}

/**
 * Henry's constant H_i(T) in Pa for a dissolved gas.
 * ln(H) = A + B/T + C·ln(T) + D·T
 * This follows the Aspen Plus HENRY-AP form.
 */
export function henryConstant_Pa(coeffs: HenryCoeffs, T_K: number): number {
  const lnH = coeffs.A + coeffs.B / T_K + coeffs.C * Math.log(T_K) + coeffs.D * T_K;
  return Math.exp(Math.max(-100, Math.min(100, lnH)));
}

/**
 * K-value for a Henry's law component (dissolved gas) in a solvent:
 *   K_i = H_i(T) / (P · φ_i^V)
 *
 * For ideal vapor phase, φ_i^V = 1, so K_i = H_i(T) / P.
 * The activity coefficient γ_i can optionally modify this:
 *   K_i = γ_i · H_i(T) / P
 */
export function KvalueHenry(
  henryCoeffs: HenryCoeffs,
  T_K: number,
  P_Pa: number,
  gamma_i: number = 1,
): number {
  const H = henryConstant_Pa(henryCoeffs, T_K);
  return gamma_i * H / P_Pa;
}

export function computeKvalues(
  compounds: CanopyCompound[],
  x: number[],
  y: number[],
  T_K: number,
  P_Pa: number,
  fluidPackage: FluidPackageType,
  interactionParams?: InteractionParams,
): number[] {
  switch (fluidPackage) {
    case 'NRTL': {
      if (!interactionParams?.nrtl) return KvaluesIdeal(compounds, T_K, P_Pa);
      const gamma = nrtlGamma(x, T_K, interactionParams.nrtl);
      return KvaluesGammaPhi(compounds, x, T_K, P_Pa, gamma);
    }
    case 'Electrolyte-NRTL': {
      if (!interactionParams?.elecNrtl) return KvaluesIdeal(compounds, T_K, P_Pa);
      const gamma = elecNRTLGamma(x, T_K, interactionParams.elecNrtl);
      const K = KvaluesGammaPhi(compounds, x, T_K, P_Pa, gamma);
      return K.map((value, index) => {
        const charge = compounds[index]?.chargeNumber ?? 0;
        if (charge !== 0) return 1e-12;
        return value;
      });
    }
    case 'Wilson': {
      if (!interactionParams?.wilson) return KvaluesIdeal(compounds, T_K, P_Pa);
      const gamma = wilsonGamma(x, T_K, interactionParams.wilson);
      return KvaluesGammaPhi(compounds, x, T_K, P_Pa, gamma);
    }
    case 'UNIQUAC': {
      if (!interactionParams?.uniquac) return KvaluesIdeal(compounds, T_K, P_Pa);
      const gamma = uniquacGamma(x, T_K, interactionParams.uniquac);
      return KvaluesGammaPhi(compounds, x, T_K, P_Pa, gamma);
    }
    case 'UNIFAC': {
      if (!interactionParams?.unifac) return KvaluesIdeal(compounds, T_K, P_Pa);
      const gamma = unifacGamma(x, T_K, interactionParams.unifac.compGroups, interactionParams.unifac.data);
      return KvaluesGammaPhi(compounds, x, T_K, P_Pa, gamma);
    }
    case 'UNIFAC-DMD': {
      if (!interactionParams?.unifacDmd) return KvaluesIdeal(compounds, T_K, P_Pa);
      const gamma = unifacDMDGamma(x, T_K, interactionParams.unifacDmd.compGroups, interactionParams.unifacDmd.data);
      return KvaluesGammaPhi(compounds, x, T_K, P_Pa, gamma);
    }
    case 'Peng-Robinson': {
      return KvaluesPR(compounds, x, y, T_K, P_Pa, resolveKij(interactionParams, T_K));
    }
    case 'SRK': {
      return KvaluesSRK(compounds, x, y, T_K, P_Pa, resolveKij(interactionParams, T_K));
    }
    case 'Ideal':
    default:
      return KvaluesIdeal(compounds, T_K, P_Pa);
  }
}

/**
 * Rachford-Rice equation: Σ zi(Ki-1) / (1 + β(Ki-1)) = 0
 * Solves for vapor fraction β.
 */
export function rachfordRice(z: number[], K: number[]): number {
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

/**
 * General PT flash with successive substitution for non-ideal K-values.
 * For 'Ideal' package, delegates directly to flashPT_Ideal.
 * For activity-coefficient or EOS models, iterates K = f(x,y) until convergence.
 *
 * Convergence criterion: max |ln(K_new/K_old)| < 1e-8  (Aspen-style)
 */
/** Interaction parameter bundle for all supported models */
export type InteractionParams = {
  nrtl?: NRTLParams;
  elecNrtl?: ElecNRTLParams;
  wilson?: WilsonParams;
  uniquac?: UNIQUACParams;
  unifac?: { compGroups: Record<number, number>[]; data: UNIFACData };
  unifacDmd?: { compGroups: Record<number, number>[]; data: UNIFACDMDData };
  kij?: number[][];
  /** Temperature-dependent kij: kij(T) = kij_a + kij_b·T + kij_c/T */
  kij_a?: number[][];
  kij_b?: number[][];
  kij_c?: number[][];
};

/**
 * Compute kij matrix at temperature T. Uses T-dependent form if available.
 */
function resolveKij(params: InteractionParams | undefined, T_K: number): number[][] | undefined {
  if (!params) return undefined;
  if (params.kij_a) {
    const n = params.kij_a.length;
    const kij: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        kij[i][j] = (params.kij_a[i]?.[j] ?? 0)
          + (params.kij_b?.[i]?.[j] ?? 0) * T_K
          + (params.kij_c?.[i]?.[j] ?? 0) / T_K;
      }
    }
    return kij;
  }
  return params.kij;
}

export function flashPT(
  compounds: CanopyCompound[],
  z: number[],
  T_K: number,
  P_Pa: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
): FlashResult {
  if (fluidPackage === 'Ideal') {
    return flashPT_Ideal(compounds, z, T_K, P_Pa);
  }

  const n = compounds.length;

  // Initialize K from Wilson correlation: Ki = (Pci/P)·exp(5.373·(1+ωi)·(1-Tci/T))
  let K = compounds.map((c) => {
    const wilson = (c.Pc_Pa / P_Pa) * Math.exp(5.373 * (1 + c.omega) * (1 - c.Tc_K / T_K));
    return isFinite(wilson) && wilson > 0 ? wilson : Psat_Pa(c, T_K) / P_Pa;
  });

  let x = new Array(n).fill(0);
  let y = new Array(n).fill(0);
  let beta = 0;
  let phase: Phase = 'Liquid';
  let converged = false;

  for (let outer = 0; outer < 100; outer++) {
    // Phase boundary checks with current K
    const sumZK = z.reduce((s, zi, i) => s + zi * K[i], 0);
    const sumZoverK = z.reduce((s, zi, i) => s + zi / K[i], 0);

    if (sumZK <= 1) {
      beta = 0;
      phase = 'Liquid';
      for (let i = 0; i < n; i++) { x[i] = z[i]; y[i] = z[i] * K[i]; }
      const sy = y.reduce((a, b) => a + b, 0);
      if (sy > 0) for (let i = 0; i < n; i++) y[i] /= sy;
    } else if (sumZoverK <= 1) {
      beta = 1;
      phase = 'Vapor';
      for (let i = 0; i < n; i++) { y[i] = z[i]; x[i] = z[i] / K[i]; }
      const sx = x.reduce((a, b) => a + b, 0);
      if (sx > 0) for (let i = 0; i < n; i++) x[i] /= sx;
    } else {
      beta = rachfordRice(z, K);
      phase = beta > 0.999 ? 'Vapor' : beta < 0.001 ? 'Liquid' : 'Two-Phase';
      for (let i = 0; i < n; i++) {
        x[i] = z[i] / (1 + beta * (K[i] - 1));
        y[i] = K[i] * x[i];
      }
      const sx = x.reduce((a, b) => a + b, 0);
      const sy = y.reduce((a, b) => a + b, 0);
      if (sx > 0) for (let i = 0; i < n; i++) x[i] /= sx;
      if (sy > 0) for (let i = 0; i < n; i++) y[i] /= sy;
    }

    // Update K from non-ideal model
    const K_new = computeKvalues(compounds, x, y, T_K, P_Pa, fluidPackage, interactionParams);

    // Check convergence: max |ln(K_new/K_old)| < tol
    let maxChange = 0;
    for (let i = 0; i < n; i++) {
      if (K[i] > 0 && K_new[i] > 0) {
        maxChange = Math.max(maxChange, Math.abs(Math.log(K_new[i] / K[i])));
      }
    }

    K = K_new;
    if (maxChange < 1e-8) { converged = true; break; }
  }

  const H = streamEnthalpy(compounds, T_K, beta, x, y, P_Pa, fluidPackage, resolveKij(interactionParams, T_K), interactionParams);
  return { converged, phase, vaporFraction: beta, x, y, K, T_K, P_Pa, H_Jpmol: H };
}

/**
 * General PQ flash — uses flashPT with non-ideal K-values.
 */
export function flashPQ(
  compounds: CanopyCompound[],
  z: number[],
  P_Pa: number,
  H_in_Jpmol: number,
  Q_W: number,
  totalFlow_molps: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
  T_guess?: number,
): FlashResult {
  if (fluidPackage === 'Ideal') {
    return flashPQ_Ideal(compounds, z, P_Pa, H_in_Jpmol, Q_W, totalFlow_molps, T_guess);
  }

  const H_target = totalFlow_molps > 0 ? H_in_Jpmol + Q_W / totalFlow_molps : H_in_Jpmol;
  const evalH = (T: number) => flashPT(compounds, z, T, P_Pa, fluidPackage, interactionParams).H_Jpmol;

  let Tlo = 100, Thi = 1500;
  let T = T_guess ?? 350;
  let bestT = T, bestErr = Infinity;

  // Newton + bracket
  for (let iter = 0; iter < 60; iter++) {
    const flash = flashPT(compounds, z, T, P_Pa, fluidPackage, interactionParams);
    const err = flash.H_Jpmol - H_target;
    if (Math.abs(err) < Math.abs(bestErr)) { bestT = T; bestErr = err; }
    if (Math.abs(err) < 0.1) return { ...flash, T_K: T };

    if (err > 0 && T < Thi) Thi = T;
    if (err < 0 && T > Tlo) Tlo = T;

    const dT = 0.1;
    const H2 = evalH(T + dT);
    const dHdT = (H2 - flash.H_Jpmol) / dT;

    if (Math.abs(dHdT) > 1e-10) {
      let T_new = T - err / dHdT;
      if (T_new < Tlo || T_new > Thi) T_new = (Tlo + Thi) / 2;
      T = T_new;
    } else {
      T = (Tlo + Thi) / 2;
    }
  }

  // Bisection fallback
  for (let iter = 0; iter < 80; iter++) {
    const mid = (Tlo + Thi) / 2;
    const H_mid = evalH(mid);
    const err = H_mid - H_target;
    if (Math.abs(err) < Math.abs(bestErr)) { bestT = mid; bestErr = err; }
    if (Math.abs(err) < 0.1 || (Thi - Tlo) < 0.001) break;
    if (err > 0) Thi = mid; else Tlo = mid;
  }

  return flashPT(compounds, z, bestT, P_Pa, fluidPackage, interactionParams);
}

// ────────────────────────────────────────────────────────────────
// 7. Bubble & Dew Point Solvers (Ideal)
// ────────────────────────────────────────────────────────────────

/**
 * PQ Flash — Pressure-Duty flash.
 * Given outlet pressure and duty Q (W), find outlet T, phase state.
 * H_out = H_in + Q / F, then iterate on T to match H_out.
 */
export function flashPQ_Ideal(
  compounds: CanopyCompound[],
  z: number[],
  P_Pa: number,
  H_in_Jpmol: number,
  Q_W: number,
  totalFlow_molps: number,
  T_guess?: number,
): FlashResult {
  const H_target = totalFlow_molps > 0 ? H_in_Jpmol + Q_W / totalFlow_molps : H_in_Jpmol;

  // Evaluate H at a given T
  const evalH = (T: number) => flashPT_Ideal(compounds, z, T, P_Pa).H_Jpmol;

  // 1. Find a bracket [Tlo, Thi] where H crosses H_target
  let T = T_guess ?? 350;
  let Tlo = 100, Thi = 1500;
  const H_lo_init = evalH(Tlo);
  const H_hi_init = evalH(Thi);

  // Ensure valid bracket — enthalpy should increase with T
  if (H_target < H_lo_init) Tlo = 50;
  if (H_target > H_hi_init) Thi = 2000;

  // 2. Try Newton-Raphson first (fast when it works)
  let bestT = T;
  let bestErr = Infinity;

  for (let iter = 0; iter < 60; iter++) {
    const flash = flashPT_Ideal(compounds, z, T, P_Pa);
    const err = flash.H_Jpmol - H_target;

    if (Math.abs(err) < Math.abs(bestErr)) { bestT = T; bestErr = err; }
    if (Math.abs(err) < 0.1) return { ...flash, T_K: T };

    // Tighten bracket
    if (err > 0 && T < Thi) Thi = T;
    if (err < 0 && T > Tlo) Tlo = T;

    // Numerical derivative
    const dT = 0.1;
    const H2 = evalH(T + dT);
    const dHdT = (H2 - flash.H_Jpmol) / dT;

    if (Math.abs(dHdT) > 1e-10) {
      let T_new = T - err / dHdT;
      // If Newton step goes outside bracket, use bisection instead
      if (T_new < Tlo || T_new > Thi) {
        T_new = (Tlo + Thi) / 2;
      }
      T = T_new;
    } else {
      // Flat region — bisect
      T = (Tlo + Thi) / 2;
    }
  }

  // 3. Bisection fallback if Newton didn't converge
  for (let iter = 0; iter < 80; iter++) {
    const mid = (Tlo + Thi) / 2;
    const H_mid = evalH(mid);
    const err = H_mid - H_target;

    if (Math.abs(err) < Math.abs(bestErr)) { bestT = mid; bestErr = err; }
    if (Math.abs(err) < 0.1 || (Thi - Tlo) < 0.001) break;

    if (err > 0) Thi = mid; else Tlo = mid;
  }

  return flashPT_Ideal(compounds, z, bestT, P_Pa);
}

/**
 * PV Flash — Pressure-Vapor Fraction flash.
 * Given outlet pressure and desired vapor fraction β, find T.
 * Iterate on T until flash gives the target β.
 */
export function flashPV_Ideal(
  compounds: CanopyCompound[],
  z: number[],
  P_Pa: number,
  targetBeta: number,
  T_guess?: number,
): FlashResult {
  // Clamp targetBeta
  const beta = Math.max(0, Math.min(1, targetBeta));

  // If β=0 → bubble point; β=1 → dew point
  if (beta <= 1e-9) {
    const Tbub = bubblePointT_Ideal(compounds, z, P_Pa);
    if (Tbub) return flashPT_Ideal(compounds, z, Tbub, P_Pa);
  }
  if (beta >= 1 - 1e-9) {
    const Tdew = dewPointT_Ideal(compounds, z, P_Pa);
    if (Tdew) return flashPT_Ideal(compounds, z, Tdew, P_Pa);
  }

  // General case: bisect between bubble and dew point
  let Tlo = bubblePointT_Ideal(compounds, z, P_Pa) ?? 200;
  let Thi = dewPointT_Ideal(compounds, z, P_Pa) ?? 600;
  if (Tlo > Thi) [Tlo, Thi] = [Thi, Tlo];
  Tlo = Math.max(100, Tlo - 5);
  Thi = Math.min(2000, Thi + 5);

  let T = T_guess ?? (Tlo + Thi) / 2;

  // Newton-Raphson on f(T) = β(T) - targetBeta
  for (let iter = 0; iter < 100; iter++) {
    const flash = flashPT_Ideal(compounds, z, T, P_Pa);
    const err = flash.vaporFraction - beta;
    if (Math.abs(err) < 1e-8) return flash;

    const dT = 0.1;
    const flash2 = flashPT_Ideal(compounds, z, T + dT, P_Pa);
    const dBdT = (flash2.vaporFraction - flash.vaporFraction) / dT;

    if (Math.abs(dBdT) < 1e-15) {
      // Fallback: bisection
      const flo = flashPT_Ideal(compounds, z, Tlo, P_Pa).vaporFraction - beta;
      const fhi = flashPT_Ideal(compounds, z, Thi, P_Pa).vaporFraction - beta;
      T = (Tlo + Thi) / 2;
      if (flo * fhi > 0) break;
      for (let j = 0; j < 60; j++) {
        const mid = (Tlo + Thi) / 2;
        const fm = flashPT_Ideal(compounds, z, mid, P_Pa).vaporFraction - beta;
        if (Math.abs(fm) < 1e-8) return flashPT_Ideal(compounds, z, mid, P_Pa);
        if (fm * flo < 0) Thi = mid; else Tlo = mid;
      }
      break;
    }

    let T_new = T - err / dBdT;
    if (T_new < Tlo) T_new = (T + Tlo) / 2;
    if (T_new > Thi) T_new = (T + Thi) / 2;
    T = T_new;
  }

  return flashPT_Ideal(compounds, z, T, P_Pa);
}

/**
 * Non-ideal PV flash: given pressure and target vapor fraction, find T.
 * Uses general flashPT with the specified fluid package.
 */
export function flashPV(
  compounds: CanopyCompound[],
  z: number[],
  P_Pa: number,
  targetBeta: number,
  fluidPackage: FluidPackageType,
  interactionParams?: InteractionParams,
  T_guess?: number,
): FlashResult {
  if (fluidPackage === 'Ideal') return flashPV_Ideal(compounds, z, P_Pa, targetBeta, T_guess);

  const beta = Math.max(0, Math.min(1, targetBeta));

  // Bracket: bubble and dew points
  let Tlo = bubblePointT(compounds, z, P_Pa, fluidPackage, interactionParams) ?? bubblePointT_Ideal(compounds, z, P_Pa) ?? 200;
  let Thi = dewPointT(compounds, z, P_Pa, fluidPackage, interactionParams) ?? dewPointT_Ideal(compounds, z, P_Pa) ?? 600;
  if (Tlo > Thi) [Tlo, Thi] = [Thi, Tlo];
  Tlo = Math.max(100, Tlo - 5);
  Thi = Math.min(2000, Thi + 5);

  let T = T_guess ?? Tlo + beta * (Thi - Tlo);

  // Newton + bisection on f(T) = β(T) - targetBeta
  for (let iter = 0; iter < 100; iter++) {
    const flash = flashPT(compounds, z, T, P_Pa, fluidPackage, interactionParams);
    const err = flash.vaporFraction - beta;
    if (Math.abs(err) < 1e-8) return flash;

    const dT = 0.1;
    const flash2 = flashPT(compounds, z, T + dT, P_Pa, fluidPackage, interactionParams);
    const dBdT = (flash2.vaporFraction - flash.vaporFraction) / dT;

    if (Math.abs(dBdT) < 1e-15) {
      // Fallback: bisection
      for (let j = 0; j < 60; j++) {
        const mid = (Tlo + Thi) / 2;
        const fm = flashPT(compounds, z, mid, P_Pa, fluidPackage, interactionParams).vaporFraction - beta;
        if (Math.abs(fm) < 1e-8) return flashPT(compounds, z, mid, P_Pa, fluidPackage, interactionParams);
        if (fm * (flashPT(compounds, z, Tlo, P_Pa, fluidPackage, interactionParams).vaporFraction - beta) < 0) Thi = mid; else Tlo = mid;
      }
      break;
    }

    let T_new = T - err / dBdT;
    if (T_new < Tlo) T_new = (T + Tlo) / 2;
    if (T_new > Thi) T_new = (T + Thi) / 2;
    T = T_new;
  }

  return flashPT(compounds, z, T, P_Pa, fluidPackage, interactionParams);
}

// ────────────────────────────────────────────────────────────────
// 8. Liquid Density — DIPPR 105 (Rackett) or simplified
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR 105: ρ_L = A / B^(1 + (1-T/C)^D)  [kmol/m³]
 * Fallback: Rackett equation using Tc, Pc, omega.
 */
export function liquidDensity_kmolpm3(comp: CanopyCompound, T_K: number): number {
  const c = comp.dipprCoeffs['LiquidDensity'];
  if (c) {
    const { A, B, C, D } = c;
    const density = A / Math.pow(B, 1 + Math.pow(1 - T_K / C, D));
    if (isFinite(density) && density > 0) return density;
  }
  // Fallback: Rackett equation
  // V_L = (R·Tc/Pc)·Z_RA^(1+(1-Tr)^(2/7))  where Z_RA ≈ 0.29056 - 0.08775·ω
  const Tr = T_K / comp.Tc_K;
  if (Tr >= 1) return NaN; // Above critical
  const Z_RA = 0.29056 - 0.08775 * comp.omega;
  const V_m = (R * comp.Tc_K / comp.Pc_Pa) * Math.pow(Z_RA, 1 + Math.pow(1 - Tr, 2 / 7));
  return V_m > 0 ? 1 / V_m : NaN; // kmol/m³ (= 1/V_m with V_m in m³/kmol)
}

/**
 * Liquid molar volume (m³/mol) — convenience wrapper.
 */
export function liquidMolarVolume_m3pmol(comp: CanopyCompound, T_K: number): number {
  const rho = liquidDensity_kmolpm3(comp, T_K);
  return isFinite(rho) && rho > 0 ? 1 / (rho * 1000) : NaN; // m³/mol
}

/**
 * Mixture liquid molar volume (m³/mol) — ideal mixing.
 */
export function mixtureLiquidMolarVolume(
  compounds: CanopyCompound[], x: number[], T_K: number
): number {
  let V = 0;
  for (let i = 0; i < compounds.length; i++) {
    const Vi = liquidMolarVolume_m3pmol(compounds[i], T_K);
    if (!isFinite(Vi)) return NaN;
    V += x[i] * Vi;
  }
  return V;
}

/**
 * Liquid molar volume from EOS (m³/mol) — PR or SRK.
 * Uses Z_liquid to compute V = Z·R·T/P.
 * More accurate at high pressures than DIPPR correlations.
 */
export function liquidMolarVolume_EOS(
  compounds: CanopyCompound[],
  x: number[],
  T_K: number,
  P_Pa: number,
  eos: 'PR' | 'SRK' = 'PR',
  kij?: number[][],
): number {
  const Z = compressibilityFactor(compounds, x, T_K, P_Pa, false, eos, kij);
  return Z * R * T_K / P_Pa;
}

// ────────────────────────────────────────────────────────────────
// COSTALD Compressed Liquid Density
// ────────────────────────────────────────────────────────────────

/**
 * COSTALD (Corresponding States Liquid Density) — Hankinson & Thomson (1979).
 * Saturated liquid molar volume (m³/mol).
 *
 * V_s = V* · V^(0)(Tr) · [1 - ω_SRK · V^(δ)(Tr)]
 * where V* is a characteristic volume (≈ Vc, but fitted),
 * V^(0) and V^(δ) are universal functions of Tr.
 */
export function COSTALD_Vs_m3pmol(
  Tc_K: number, Vc_m3pmol: number, omega: number, T_K: number,
): number {
  const Tr = T_K / Tc_K;
  if (Tr >= 1) return NaN;

  // V^(0) = 1 - 1.52816·(1-Tr)^(1/3) + 1.43907·(1-Tr)^(2/3) - 0.81446·(1-Tr) + 0.190454·(1-Tr)^(4/3)
  const tau = 1 - Tr;
  const tau13 = Math.pow(tau, 1 / 3);
  const V0 = 1 - 1.52816 * tau13 + 1.43907 * tau13 * tau13
    - 0.81446 * tau + 0.190454 * tau * tau13;

  // V^(δ) = (-0.296123 + 0.386914/Tr - 0.0427258/Tr² - 0.0480645/Tr³) / (Tr - 1.00001)
  const Vd = (-0.296123 + 0.386914 / Tr - 0.0427258 / (Tr * Tr) - 0.0480645 / (Tr * Tr * Tr))
    / (Tr - 1.00001);

  return Vc_m3pmol * V0 * (1 - omega * Vd);
}

/**
 * Compressed liquid density via Tait equation (Thomson-Brobst-Hankinson).
 * Applies pressure correction to COSTALD saturated density.
 *
 * V(T,P) = V_s(T) · [1 - C·ln((B+P)/(B+P_sat))]
 * B/Pc = -1 + a·(1-Tr)^(1/3) + b·(1-Tr)^(2/3) + d·(1-Tr) + e·(1-Tr)^(4/3)
 * C = 0.0861488 + 0.0344483·ω
 */
export function COSTALD_compressed_m3pmol(
  Tc_K: number, Pc_Pa: number, Vc_m3pmol: number, omega: number,
  T_K: number, P_Pa: number,
): number {
  const Vs = COSTALD_Vs_m3pmol(Tc_K, Vc_m3pmol, omega, T_K);
  if (!isFinite(Vs)) return NaN;

  const Tr = T_K / Tc_K;
  const tau = 1 - Tr;
  if (tau <= 0) return NaN;
  const tau13 = Math.pow(tau, 1 / 3);

  // Tait B/Pc
  const BoverPc = -1 + 9.070217 * tau13 - 62.45326 * tau13 * tau13
    + 135.1102 * tau + tau13 * tau * (-4.79594 + 0.250047 * tau13 * tau13);
  const B = BoverPc * Pc_Pa;

  const C = 0.0861488 + 0.0344483 * omega;

  // Saturated pressure (approximate via Antoine if not available)
  const Psat = Pc_Pa * Math.exp(5.373 * (1 + omega) * (1 - Tc_K / T_K));

  const ratio = (B + P_Pa) / (B + Psat);
  if (ratio <= 0) return Vs;

  return Vs * (1 - C * Math.log(ratio));
}

// ────────────────────────────────────────────────────────────────
// Second Virial Coefficient (Pitzer Correlation)
// ────────────────────────────────────────────────────────────────

/**
 * Second virial coefficient B (m³/mol) from Pitzer correlation.
 *
 * B·Pc/(R·Tc) = B^(0) + ω·B^(1)
 * B^(0) = 0.083 - 0.422/Tr^1.6
 * B^(1) = 0.139 - 0.172/Tr^4.2
 *
 * Used as a simpler alternative to Hayden-O'Connell for moderately
 * non-ideal vapor phases.
 */
export function secondVirial_m3pmol(
  Tc_K: number, Pc_Pa: number, omega: number, T_K: number,
): number {
  const Tr = T_K / Tc_K;
  const B0 = 0.083 - 0.422 / Math.pow(Tr, 1.6);
  const B1 = 0.139 - 0.172 / Math.pow(Tr, 4.2);
  return (R * Tc_K / Pc_Pa) * (B0 + omega * B1);
}

/**
 * Fugacity coefficient from truncated virial equation (for pure component):
 *   ln(φ) = B·P/(R·T)
 * Valid at low to moderate pressures where 3rd+ virial coefficients are negligible.
 */
export function fugacityVirial(
  Tc_K: number, Pc_Pa: number, omega: number,
  T_K: number, P_Pa: number,
): number {
  const B = secondVirial_m3pmol(Tc_K, Pc_Pa, omega, T_K);
  return Math.exp(B * P_Pa / (R * T_K));
}

// ────────────────────────────────────────────────────────────────
// Mixture Second Virial Coefficient
// ────────────────────────────────────────────────────────────────

/**
 * Mixture second virial coefficient using classical mixing rules.
 * B_mix = Σ_i Σ_j y_i·y_j·B_ij
 * B_ij computed from pseudo-critical properties:
 *   Tc_ij = √(Tc_i·Tc_j)·(1-kij), Pc_ij = Z_c·R·Tc_ij/Vc_ij,
 *   ω_ij = (ω_i+ω_j)/2, Vc_ij = ((Vc_i^(1/3)+Vc_j^(1/3))/2)³
 */
export function mixtureSecondVirial_m3pmol(
  compounds: CanopyCompound[],
  y: number[],
  T_K: number,
  kij?: number[][],
): number {
  const n = compounds.length;
  let Bmix = 0;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const ci = compounds[i], cj = compounds[j];
      const k = kij?.[i]?.[j] ?? 0;
      const Tc_ij = Math.sqrt(ci.Tc_K * cj.Tc_K) * (1 - k);
      const omega_ij = (ci.omega + cj.omega) / 2;
      // Use geometric mean for Pc_ij (simplified)
      const Pc_ij = Math.sqrt(ci.Pc_Pa * cj.Pc_Pa);
      Bmix += y[i] * y[j] * secondVirial_m3pmol(Tc_ij, Pc_ij, omega_ij, T_K);
    }
  }
  return Bmix;
}

// ────────────────────────────────────────────────────────────────
// 9. Bubble & Dew Point Solvers (Ideal)
// ────────────────────────────────────────────────────────────────

/** Bubble point temperature at given P (Raoult's Law). Newton-Raphson. */
export function bubblePointT_Ideal(
  compounds: CanopyCompound[],
  z: number[],
  P_Pa: number,
  initialGuess_K?: number,
): number | null {
  // Objective: Σ zi·Psat_i(T) / P = 1
  let T = initialGuess_K ?? 350; // Initial guess
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
  initialGuess_K?: number,
): number | null {
  // Objective: Σ zi·P / Psat_i(T) = 1  ⟹  Σ zi / Ki = 1
  let T = initialGuess_K ?? 380;
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

/**
 * General bubble point temperature using any K-value model.
 * Objective: Σ zi·Ki(T) = 1  at bubble point (all liquid, incipient vapor).
 */
export function bubblePointT(
  compounds: CanopyCompound[],
  z: number[],
  P_Pa: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
): number | null {
  if (fluidPackage === 'Ideal') return bubblePointT_Ideal(compounds, z, P_Pa);

  let T = bubblePointT_Ideal(compounds, z, P_Pa) ?? 350;
  const y_est = z.slice(); // initial vapor estimate = feed

  for (let iter = 0; iter < 100; iter++) {
    const K = computeKvalues(compounds, z, y_est, T, P_Pa, fluidPackage, interactionParams);
    let f = -1;
    for (let i = 0; i < compounds.length; i++) {
      f += z[i] * K[i];
      y_est[i] = z[i] * K[i]; // y = z·K at bubble point
    }
    const sy = y_est.reduce((a, b) => a + b, 0);
    if (sy > 0) for (let i = 0; i < compounds.length; i++) y_est[i] /= sy;

    if (Math.abs(f) < 1e-8) return T;

    // Numerical derivative
    const dT = 0.05;
    const K2 = computeKvalues(compounds, z, y_est, T + dT, P_Pa, fluidPackage, interactionParams);
    let f2 = -1;
    for (let i = 0; i < compounds.length; i++) f2 += z[i] * K2[i];
    const df = (f2 - f) / dT;
    if (Math.abs(df) < 1e-15) break;

    T -= f / df;
    if (T < 100 || T > 1000) break;
  }
  return null;
}

/**
 * General dew point temperature using any K-value model.
 * Objective: Σ zi/Ki(T) = 1  at dew point (all vapor, incipient liquid).
 */
export function dewPointT(
  compounds: CanopyCompound[],
  z: number[],
  P_Pa: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
): number | null {
  if (fluidPackage === 'Ideal') return dewPointT_Ideal(compounds, z, P_Pa);

  let T = dewPointT_Ideal(compounds, z, P_Pa) ?? 380;
  const x_est = z.slice(); // initial liquid estimate = feed

  for (let iter = 0; iter < 100; iter++) {
    const K = computeKvalues(compounds, x_est, z, T, P_Pa, fluidPackage, interactionParams);
    let f = -1;
    for (let i = 0; i < compounds.length; i++) {
      const Ki = K[i] > 1e-30 ? K[i] : 1e-30;
      f += z[i] / Ki;
      x_est[i] = z[i] / Ki; // x = z/K at dew point
    }
    const sx = x_est.reduce((a, b) => a + b, 0);
    if (sx > 0) for (let i = 0; i < compounds.length; i++) x_est[i] /= sx;

    if (Math.abs(f) < 1e-8) return T;

    const dT = 0.05;
    const K2 = computeKvalues(compounds, x_est, z, T + dT, P_Pa, fluidPackage, interactionParams);
    let f2 = -1;
    for (let i = 0; i < compounds.length; i++) {
      const Ki = K2[i] > 1e-30 ? K2[i] : 1e-30;
      f2 += z[i] / Ki;
    }
    const df = (f2 - f) / dT;
    if (Math.abs(df) < 1e-15) break;

    T -= f / df;
    if (T < 100 || T > 1000) break;
  }
  return null;
}

export function bubblePointTWaterSaturated(
  compounds: CanopyCompound[],
  z: number[],
  P_Pa: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
): number | null {
  let T = bubblePointT(compounds, z, P_Pa, fluidPackage, interactionParams) ?? 350;
  for (let iter = 0; iter < 60; iter++) {
    const saturated = waterSaturatedComposition(compounds, z, T);
    if (!saturated) return bubblePointT(compounds, z, P_Pa, fluidPackage, interactionParams);
    const bubble = bubblePointT(compounds, saturated.zSaturated, P_Pa, fluidPackage, interactionParams);
    if (bubble == null) return null;
    const err = bubble - T;
    if (Math.abs(err) < 1e-6) return bubble;
    T += 0.6 * err;
    if (T < 100 || T > 1200) break;
  }
  const saturated = waterSaturatedComposition(compounds, z, T);
  return saturated ? bubblePointT(compounds, saturated.zSaturated, P_Pa, fluidPackage, interactionParams) : null;
}

export function dewPointTWaterSaturated(
  compounds: CanopyCompound[],
  z: number[],
  P_Pa: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
): number | null {
  let T = dewPointT(compounds, z, P_Pa, fluidPackage, interactionParams) ?? 380;
  for (let iter = 0; iter < 60; iter++) {
    const saturated = waterSaturatedComposition(compounds, z, T);
    if (!saturated) return dewPointT(compounds, z, P_Pa, fluidPackage, interactionParams);
    const dew = dewPointT(compounds, saturated.zSaturated, P_Pa, fluidPackage, interactionParams);
    if (dew == null) return null;
    const err = dew - T;
    if (Math.abs(err) < 1e-6) return dew;
    T += 0.6 * err;
    if (T < 100 || T > 1200) break;
  }
  const saturated = waterSaturatedComposition(compounds, z, T);
  return saturated ? dewPointT(compounds, saturated.zSaturated, P_Pa, fluidPackage, interactionParams) : null;
}

export function flashPTWaterSaturated(
  compounds: CanopyCompound[],
  z: number[],
  T_K: number,
  P_Pa: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
): FlashResult & { waterAddedFraction: number; excessFreeWaterFraction: number } {
  const saturated = waterSaturatedComposition(compounds, z, T_K);
  const flash = flashPT(
    compounds,
    saturated?.zSaturated ?? z,
    T_K,
    P_Pa,
    fluidPackage,
    interactionParams,
  );
  return {
    ...flash,
    waterAddedFraction: saturated?.waterAddedFraction ?? 0,
    excessFreeWaterFraction: saturated?.excessFreeWaterFraction ?? 0,
  };
}

// ────────────────────────────────────────────────────────────────
// Bubble & Dew Point Pressure Solvers
// ────────────────────────────────────────────────────────────────

/**
 * Bubble point pressure at given T (Raoult's Law).
 * P_bub = Σ x_i · P_sat_i(T)
 */
export function bubblePointP_Ideal(
  compounds: CanopyCompound[], z: number[], T_K: number,
): number | null {
  let Pbub = 0;
  for (let i = 0; i < compounds.length; i++) {
    const ps = Psat_Pa(compounds[i], T_K);
    if (!isFinite(ps) || ps <= 0) return null;
    Pbub += z[i] * ps;
  }
  return Pbub > 0 ? Pbub : null;
}

/**
 * Dew point pressure at given T (Raoult's Law).
 * 1/P_dew = Σ y_i / P_sat_i(T)
 */
export function dewPointP_Ideal(
  compounds: CanopyCompound[], z: number[], T_K: number,
): number | null {
  let invP = 0;
  for (let i = 0; i < compounds.length; i++) {
    const ps = Psat_Pa(compounds[i], T_K);
    if (!isFinite(ps) || ps <= 0) return null;
    invP += z[i] / ps;
  }
  return invP > 0 ? 1 / invP : null;
}

/**
 * Bubble point pressure with activity coefficients / EOS.
 * Iterative Newton on Σ x_i·γ_i·P_sat_i / P - 1 = 0.
 */
export function bubblePointP(
  compounds: CanopyCompound[],
  z: number[],
  T_K: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
): number | null {
  if (fluidPackage === 'Ideal') return bubblePointP_Ideal(compounds, z, T_K);

  // Initial guess from Raoult
  let P = bubblePointP_Ideal(compounds, z, T_K) ?? 101325;
  const y_est = z.slice();

  for (let iter = 0; iter < 100; iter++) {
    const K = computeKvalues(compounds, z, y_est, T_K, P, fluidPackage, interactionParams);
    let sumKz = 0;
    for (let i = 0; i < compounds.length; i++) {
      sumKz += z[i] * K[i];
      y_est[i] = z[i] * K[i];
    }
    // Normalize y
    if (sumKz > 0) for (let i = 0; i < compounds.length; i++) y_est[i] /= sumKz;

    const f = sumKz - 1;
    if (Math.abs(f) < 1e-8) return P;

    // Since K ∝ 1/P for activity models, f ∝ 1/P → df/dP ≈ -f/P
    P *= sumKz; // Direct substitution: P_new = P · Σ(Kz)
    if (P < 10 || P > 1e9) break;
  }
  return null;
}

/**
 * Dew point pressure with activity coefficients / EOS.
 */
export function dewPointP(
  compounds: CanopyCompound[],
  z: number[],
  T_K: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
): number | null {
  if (fluidPackage === 'Ideal') return dewPointP_Ideal(compounds, z, T_K);

  let P = dewPointP_Ideal(compounds, z, T_K) ?? 101325;
  const x_est = z.slice();

  for (let iter = 0; iter < 100; iter++) {
    const K = computeKvalues(compounds, x_est, z, T_K, P, fluidPackage, interactionParams);
    let sumZoverK = 0;
    for (let i = 0; i < compounds.length; i++) {
      const Ki = Math.max(K[i], 1e-30);
      sumZoverK += z[i] / Ki;
      x_est[i] = z[i] / Ki;
    }
    if (sumZoverK > 0) for (let i = 0; i < compounds.length; i++) x_est[i] /= sumZoverK;

    const f = sumZoverK - 1;
    if (Math.abs(f) < 1e-8) return P;

    P /= sumZoverK; // Direct substitution
    if (P < 10 || P > 1e9) break;
  }
  return null;
}

// ────────────────────────────────────────────────────────────────
// 10. Three-Phase VLLE Flash
// ────────────────────────────────────────────────────────────────

export interface VLLEFlashResult {
  converged: boolean;
  T_K: number;
  P_Pa: number;
  phase: 'Liquid' | 'Vapor' | 'Two-Phase' | 'Three-Phase';
  vaporFraction: number;
  liquidIFraction: number;
  liquidIIFraction: number;
  y: number[];
  xI: number[];
  xII: number[];
  K_VI: number[];
  K_LII_LI: number[];
  H_Jpmol: number;
}

/**
 * Three-phase VLLE flash (PT specification).
 * Algorithm:
 * 1. Perform normal VLE flash
 * 2. Stability test on liquid phase (Michelsen TPD criterion)
 * 3. If unstable, converge three-phase split via Newton on Rachford-Rice
 *
 * Reference: Michelsen, M.L. (1982) Fluid Phase Equilibria 9:1-19.
 */
export function flashPT_VLLE(
  compounds: CanopyCompound[],
  z: number[],
  T_K: number,
  P_Pa: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
): VLLEFlashResult {
  const n = compounds.length;

  // Step 1: Normal VLE flash
  const vleResult = flashPT(compounds, z, T_K, P_Pa, fluidPackage, interactionParams);

  const makeVLEReturn = (vle: FlashResult): VLLEFlashResult => ({
    converged: vle.converged,
    T_K, P_Pa,
    phase: vle.phase === 'Two-Phase' ? 'Two-Phase' : vle.phase,
    vaporFraction: vle.vaporFraction,
    liquidIFraction: 1 - vle.vaporFraction,
    liquidIIFraction: 0,
    y: vle.y, xI: vle.x, xII: vle.x,
    K_VI: vle.K,
    K_LII_LI: new Array(n).fill(1),
    H_Jpmol: vle.H_Jpmol,
  });

  // For Ideal or extreme phases, LLE impossible
  if (fluidPackage === 'Ideal' || vleResult.vaporFraction > 0.999 || vleResult.vaporFraction < 0.001) {
    return makeVLEReturn(vleResult);
  }

  // Step 2: Stability test on liquid phase (TPD)
  const gammaLiq = computeGamma(compounds, vleResult.x, T_K, fluidPackage, interactionParams);

  let isUnstable = false;
  let bestW: number[] | null = null;
  let bestTPD = 0;

  for (let trial = 0; trial < n; trial++) {
    const w = new Array(n).fill(1e-6);
    w[trial] = 1 - (n - 1) * 1e-6;

    for (let iter = 0; iter < 50; iter++) {
      const gammaW = computeGamma(compounds, w, T_K, fluidPackage, interactionParams);
      let sumW = 0;
      for (let i = 0; i < n; i++) {
        w[i] = vleResult.x[i] * gammaLiq[i] / (gammaW[i] || 1e-30);
        if (w[i] < 1e-15) w[i] = 1e-15;
        sumW += w[i];
      }
      for (let i = 0; i < n; i++) w[i] /= sumW;

      const gammaW2 = computeGamma(compounds, w, T_K, fluidPackage, interactionParams);
      let tpd = 0;
      for (let i = 0; i < n; i++) {
        tpd += w[i] * (
          Math.log(Math.max(w[i], 1e-30)) + Math.log(Math.max(gammaW2[i], 1e-30))
          - Math.log(Math.max(vleResult.x[i], 1e-30)) - Math.log(Math.max(gammaLiq[i], 1e-30))
        );
      }
      if (tpd < bestTPD - 1e-6) {
        bestTPD = tpd;
        bestW = [...w];
        isUnstable = true;
      }
    }
  }

  if (!isUnstable || !bestW) return makeVLEReturn(vleResult);

  // Step 3: Three-phase convergence
  let xI = [...vleResult.x];
  let xII = [...bestW];
  let y = [...vleResult.y];
  let V = vleResult.vaporFraction;
  let LI = (1 - V) * 0.5;
  let LII = (1 - V) * 0.5;
  let converged = false;

  for (let iter = 0; iter < 200; iter++) {
    const gammaI = computeGamma(compounds, xI, T_K, fluidPackage, interactionParams);
    const gammaII = computeGamma(compounds, xII, T_K, fluidPackage, interactionParams);

    const K_VI = compounds.map((c, i) => {
      const ps = Psat_Pa(c, T_K);
      return isFinite(ps) ? gammaI[i] * ps / P_Pa : 1;
    });
    const K_LL = compounds.map((_, i) => gammaI[i] / (gammaII[i] || 1e-30));

    // Newton on two-variable Rachford-Rice
    for (let inner = 0; inner < 30; inner++) {
      LI = 1 - V - LII;
      if (LI < 0) { LI = 0.01; LII = Math.max(0.01, 1 - V - LI); }

      let f1 = -1, f2 = -1;
      let df1dV = 0, df1dL2 = 0, df2dV = 0, df2dL2 = 0;

      for (let i = 0; i < n; i++) {
        const D = V * K_VI[i] + LI + LII * K_LL[i];
        const D2 = D * D || 1e-60;
        const dDdV = K_VI[i] - 1;
        const dDdL2 = K_LL[i] - 1;
        f1 += z[i] / D;
        f2 += z[i] * K_LL[i] / D;
        df1dV -= z[i] * dDdV / D2;
        df1dL2 -= z[i] * dDdL2 / D2;
        df2dV -= z[i] * K_LL[i] * dDdV / D2;
        df2dL2 -= z[i] * K_LL[i] * dDdL2 / D2;
      }

      if (Math.abs(f1) < 1e-10 && Math.abs(f2) < 1e-10) break;
      const det = df1dV * df2dL2 - df1dL2 * df2dV;
      if (Math.abs(det) < 1e-30) break;

      const dV = -(f1 * df2dL2 - f2 * df1dL2) / det;
      const dL2 = -(df1dV * f2 - df2dV * f1) / det;
      V = Math.max(0, Math.min(1, V + dV * 0.8));
      LII = Math.max(0, Math.min(1, LII + dL2 * 0.8));
      if (V + LII > 1) { const sc = 0.98 / (V + LII); V *= sc; LII *= sc; }
    }

    LI = 1 - V - LII;

    const xI_new = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      const D = V * K_VI[i] + LI + LII * K_LL[i];
      xI_new[i] = z[i] / (D || 1e-30);
    }
    const sumXI = xI_new.reduce((s, v) => s + v, 0) || 1;
    for (let i = 0; i < n; i++) xI_new[i] /= sumXI;

    let maxDx = 0;
    for (let i = 0; i < n; i++) maxDx = Math.max(maxDx, Math.abs(xI_new[i] - xI[i]));

    xI = xI_new;
    y = xI.map((xi, i) => K_VI[i] * xi);
    xII = xI.map((xi, i) => K_LL[i] * xi);
    const sumY = y.reduce((s, v) => s + v, 0) || 1;
    const sumXII = xII.reduce((s, v) => s + v, 0) || 1;
    for (let i = 0; i < n; i++) { y[i] /= sumY; xII[i] /= sumXII; }

    if (maxDx < 1e-8) { converged = true; break; }
  }

  let phase: 'Liquid' | 'Vapor' | 'Two-Phase' | 'Three-Phase' = 'Three-Phase';
  if (V < 0.001 && LII < 0.001) phase = 'Liquid';
  else if (V > 0.999) phase = 'Vapor';
  else if (LII < 0.001) phase = 'Two-Phase';

  return {
    converged, T_K, P_Pa, phase, vaporFraction: V,
    liquidIFraction: LI, liquidIIFraction: LII,
    y, xI, xII,
    K_VI: compounds.map((_, i) => y[i] / (xI[i] || 1e-30)),
    K_LII_LI: compounds.map((_, i) => xII[i] / (xI[i] || 1e-30)),
    H_Jpmol: streamEnthalpy(compounds, T_K, V, xI, y),
  };
}

/** Compute activity coefficients for the given thermodynamic model */
export function computeGamma(
  compounds: CanopyCompound[],
  x: number[],
  T_K: number,
  fluidPackage: FluidPackageType,
  interactionParams?: InteractionParams,
): number[] {
  switch (fluidPackage) {
    case 'NRTL':
      return interactionParams?.nrtl ? nrtlGamma(x, T_K, interactionParams.nrtl) : x.map(() => 1);
    case 'Electrolyte-NRTL':
      return interactionParams?.elecNrtl ? elecNRTLGamma(x, T_K, interactionParams.elecNrtl) : x.map(() => 1);
    case 'Wilson':
      return interactionParams?.wilson ? wilsonGamma(x, T_K, interactionParams.wilson) : x.map(() => 1);
    case 'UNIQUAC':
      return interactionParams?.uniquac ? uniquacGamma(x, T_K, interactionParams.uniquac) : x.map(() => 1);
    case 'UNIFAC':
      return interactionParams?.unifac
        ? unifacGamma(x, T_K, interactionParams.unifac.compGroups, interactionParams.unifac.data)
        : x.map(() => 1);
    case 'UNIFAC-DMD':
      return interactionParams?.unifacDmd
        ? unifacDMDGamma(x, T_K, interactionParams.unifacDmd.compGroups, interactionParams.unifacDmd.data)
        : x.map(() => 1);
    default:
      return x.map(() => 1);
  }
}

// ────────────────────────────────────────────────────────────────
// 11. Transport Properties (DIPPR Correlations)
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR 101: ln(μ) = A + B/T + C·ln(T) + D·T^E
 * Returns liquid viscosity in Pa·s
 */
export function liquidViscosity_Pas(comp: CanopyCompound, T_K: number): number {
  const c = comp.transportCoeffs?.LiquidViscosity;
  if (!c) return NaN;
  const E = c.E ?? 2;
  const lnMu = c.A + c.B / T_K + c.C * Math.log(T_K) + c.D * Math.pow(T_K, E);
  const mu = Math.exp(lnMu);
  return isFinite(mu) && mu > 0 ? mu : NaN;
}

/**
 * DIPPR 102: μ = A·T^B / (1 + C/T + D/T²)
 * Returns vapor viscosity in Pa·s
 */
export function vaporViscosity_Pas(comp: CanopyCompound, T_K: number): number {
  const c = comp.transportCoeffs?.VaporViscosity;
  if (!c) return NaN;
  const mu = c.A * Math.pow(T_K, c.B) / (1 + c.C / T_K + c.D / (T_K * T_K));
  return isFinite(mu) && mu > 0 ? mu : NaN;
}

/**
 * DIPPR 100: λ = A + B·T + C·T² + D·T³ + E·T⁴
 * Returns liquid thermal conductivity in W/(m·K)
 */
export function liquidThermalConductivity_WpmK(comp: CanopyCompound, T_K: number): number {
  const c = comp.transportCoeffs?.LiquidThermalConductivity;
  if (!c) return NaN;

  // Check TRNSWT S3 to determine equation form:
  // 123 = Sato-Riedel: k = A·(1 + B·τ^(1/3) + C·τ^(2/3) + D·τ), τ = 1 - T/Tc
  // 100 = polynomial: k = A + B·T + C·T² + D·T³ + E·T⁴
  const eqn = comp.TRNSWT?.[2] ?? 100;
  let k: number;
  if (eqn === 123) {
    const Tr = T_K / comp.Tc_K;
    if (Tr >= 1) return NaN;
    const tau = 1 - Tr;
    k = c.A * (1 + c.B * Math.cbrt(tau) + c.C * Math.pow(tau, 2 / 3) + c.D * tau);
  } else {
    const E = c.E ?? 0;
    k = c.A + c.B * T_K + c.C * T_K * T_K + c.D * T_K * T_K * T_K + E * Math.pow(T_K, 4);
  }
  return isFinite(k) && k > 0 ? k : NaN;
}

/**
 * DIPPR 102: λ = A·T^B / (1 + C/T + D/T²)
 * Returns vapor thermal conductivity in W/(m·K)
 */
export function vaporThermalConductivity_WpmK(comp: CanopyCompound, T_K: number): number {
  const c = comp.transportCoeffs?.VaporThermalConductivity;
  if (!c) return NaN;
  const k = c.A * Math.pow(T_K, c.B) / (1 + c.C / T_K + c.D / (T_K * T_K));
  return isFinite(k) && k > 0 ? k : NaN;
}

/**
 * DIPPR 106: σ = A · (1-Tr)^(B + C·Tr + D·Tr² + E·Tr³)
 * Returns surface tension in N/m
 */
export function surfaceTension_Npm(comp: CanopyCompound, T_K: number): number {
  const c = comp.transportCoeffs?.SurfaceTension;
  if (!c) return NaN;
  const Tr = T_K / comp.Tc_K;
  if (Tr >= 1) return 0;
  const E = c.E ?? 0;
  const sigma = c.A * Math.pow(1 - Tr, c.B + c.C * Tr + c.D * Tr * Tr + E * Tr * Tr * Tr);
  return isFinite(sigma) && sigma >= 0 ? sigma : NaN;
}

/**
 * Mixture liquid viscosity (Grunberg-Nissan):
 *   ln(μ_mix) = Σ xi · ln(μ_i)
 */
export function mixtureLiquidViscosity_Pas(
  compounds: CanopyCompound[], x: number[], T_K: number
): number {
  let lnMu = 0;
  for (let i = 0; i < compounds.length; i++) {
    const mu_i = liquidViscosity_Pas(compounds[i], T_K);
    if (!isFinite(mu_i) || mu_i <= 0) return NaN;
    lnMu += x[i] * Math.log(mu_i);
  }
  return Math.exp(lnMu);
}

/**
 * Mixture liquid thermal conductivity (simplified linear):
 *   λ_mix = Σ xi · λ_i
 */
export function mixtureLiquidThermalConductivity_WpmK(
  compounds: CanopyCompound[], x: number[], T_K: number
): number {
  let k = 0;
  for (let i = 0; i < compounds.length; i++) {
    const ki = liquidThermalConductivity_WpmK(compounds[i], T_K);
    if (!isFinite(ki)) return NaN;
    k += x[i] * ki;
  }
  return k;
}

/**
 * Mixture surface tension (simplified linear):
 *   σ_mix = Σ xi · σ_i
 */
export function mixtureSurfaceTension_Npm(
  compounds: CanopyCompound[], x: number[], T_K: number
): number {
  let sigma = 0;
  for (let i = 0; i < compounds.length; i++) {
    const si = surfaceTension_Npm(compounds[i], T_K);
    if (!isFinite(si)) return NaN;
    sigma += x[i] * si;
  }
  return sigma;
}

/**
 * Wilke interaction parameter Φ_ij for vapor viscosity/conductivity mixing.
 * Φ_ij = [1 + (μ_i/μ_j)^0.5 · (M_j/M_i)^0.25]^2 / [8·(1 + M_i/M_j)]^0.5
 * Ref: Wilke, C.R., J. Chem. Phys. 18, 517 (1950). Used in Aspen Plus.
 */
function wilkePhi(mu_i: number, mu_j: number, M_i: number, M_j: number): number {
  const ratio = Math.sqrt(mu_i / mu_j) * Math.pow(M_j / M_i, 0.25);
  return (1 + ratio) * (1 + ratio) / Math.sqrt(8 * (1 + M_i / M_j));
}

/**
 * Mixture vapor viscosity (Wilke method):
 *   μ_mix = Σ (y_i · μ_i) / Σ (y_j · Φ_ij)
 * Ref: Wilke, C.R., J. Chem. Phys. 18, 517 (1950).
 * This is the default method in Aspen Plus for low-pressure vapor viscosity.
 */
export function mixtureVaporViscosity_Pas(
  compounds: CanopyCompound[], y: number[], T_K: number
): number {
  const nc = compounds.length;
  const mu = new Array(nc);
  for (let i = 0; i < nc; i++) {
    mu[i] = vaporViscosity_Pas(compounds[i], T_K);
    if (!isFinite(mu[i]) || mu[i] <= 0) return NaN;
  }
  let muMix = 0;
  for (let i = 0; i < nc; i++) {
    if (y[i] < 1e-15) continue;
    let denom = 0;
    for (let j = 0; j < nc; j++) {
      if (y[j] < 1e-15) continue;
      denom += y[j] * wilkePhi(mu[i], mu[j], compounds[i].molecularWeight, compounds[j].molecularWeight);
    }
    muMix += y[i] * mu[i] / denom;
  }
  return muMix;
}

/**
 * Mixture vapor thermal conductivity (Wassilijewa equation with Mason-Saxena modification):
 *   λ_mix = Σ (y_i · λ_i) / Σ (y_j · A_ij)
 * where A_ij uses the Wilke Φ_ij interaction parameter.
 * Ref: Wassilijewa, A., Physik. Z. 5, 737 (1904); Mason & Saxena, Phys. Fluids 1, 361 (1958).
 * This is the default method in Aspen Plus for low-pressure vapor thermal conductivity.
 */
export function mixtureVaporThermalConductivity_WpmK(
  compounds: CanopyCompound[], y: number[], T_K: number
): number {
  const nc = compounds.length;
  const mu = new Array(nc);
  const lam = new Array(nc);
  for (let i = 0; i < nc; i++) {
    mu[i] = vaporViscosity_Pas(compounds[i], T_K);
    lam[i] = vaporThermalConductivity_WpmK(compounds[i], T_K);
    if (!isFinite(mu[i]) || mu[i] <= 0 || !isFinite(lam[i]) || lam[i] <= 0) return NaN;
  }
  let lamMix = 0;
  for (let i = 0; i < nc; i++) {
    if (y[i] < 1e-15) continue;
    let denom = 0;
    for (let j = 0; j < nc; j++) {
      if (y[j] < 1e-15) continue;
      denom += y[j] * wilkePhi(mu[i], mu[j], compounds[i].molecularWeight, compounds[j].molecularWeight);
    }
    lamMix += y[i] * lam[i] / denom;
  }
  return lamMix;
}

// ────────────────────────────────────────────────────────────────
// Phase Envelope (PT diagram)
// ────────────────────────────────────────────────────────────────

/** Point on a phase envelope (bubble or dew curve). */
export interface PhaseEnvelopePoint {
  T_K: number;
  P_Pa: number;
  type: 'bubble' | 'dew';
}

/**
 * Compute the phase envelope (bubble point and dew point curves) for a mixture.
 * Returns arrays of {T, P} points tracing the two-phase boundary.
 *
 * @param nPoints Number of points on each curve
 * @param Tmin_K Starting temperature
 * @param Tmax_K Ending temperature
 */
export function phaseEnvelope(
  compounds: CanopyCompound[],
  z: number[],
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
  nPoints: number = 50,
  Tmin_K: number = 200,
  Tmax_K: number = 600,
): { bubbleCurve: PhaseEnvelopePoint[]; dewCurve: PhaseEnvelopePoint[] } {
  const bubbleCurve: PhaseEnvelopePoint[] = [];
  const dewCurve: PhaseEnvelopePoint[] = [];

  const dT = (Tmax_K - Tmin_K) / (nPoints - 1);

  for (let i = 0; i < nPoints; i++) {
    const T = Tmin_K + i * dT;

    const Pbub = fluidPackage === 'Ideal'
      ? bubblePointP_Ideal(compounds, z, T)
      : bubblePointP(compounds, z, T, fluidPackage, interactionParams);
    if (Pbub != null && isFinite(Pbub) && Pbub > 0) {
      bubbleCurve.push({ T_K: T, P_Pa: Pbub, type: 'bubble' });
    }

    const Pdew = fluidPackage === 'Ideal'
      ? dewPointP_Ideal(compounds, z, T)
      : dewPointP(compounds, z, T, fluidPackage, interactionParams);
    if (Pdew != null && isFinite(Pdew) && Pdew > 0) {
      dewCurve.push({ T_K: T, P_Pa: Pdew, type: 'dew' });
    }
  }

  return { bubbleCurve, dewCurve };
}

// ────────────────────────────────────────────────────────────────
// Txy / Pxy Diagrams
// ────────────────────────────────────────────────────────────────

/** Point on a Txy (or Pxy) diagram. */
export interface TxyPoint {
  x1: number;  // liquid composition of component 1
  y1: number;  // vapor composition of component 1
  T_K: number;
  P_Pa: number;
}

/**
 * Generate Txy diagram data at constant pressure for a binary system.
 * Returns bubble and dew temperatures as a function of x1.
 */
export function txyDiagram(
  compounds: [CanopyCompound, CanopyCompound],
  P_Pa: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
  nPoints: number = 51,
): TxyPoint[] {
  const points: TxyPoint[] = [];

  for (let i = 0; i < nPoints; i++) {
    const x1 = i / (nPoints - 1);
    const z = [x1, 1 - x1];

    const Tbub = fluidPackage === 'Ideal'
      ? bubblePointT_Ideal(compounds, z, P_Pa)
      : bubblePointT(compounds, z, P_Pa, fluidPackage, interactionParams);

    if (Tbub == null) continue;

    const K = computeKvalues(
      compounds, z, z, Tbub, P_Pa, fluidPackage, interactionParams,
    );
    const y1 = z[0] * K[0] / (z.reduce((s, zi, j) => s + zi * K[j], 0) || 1);

    points.push({ x1, y1: Math.max(0, Math.min(1, y1)), T_K: Tbub, P_Pa });
  }

  return points;
}

/**
 * Generate xy diagram data (McCabe-Thiele type) at constant pressure.
 */
export function xyDiagram(
  compounds: [CanopyCompound, CanopyCompound],
  P_Pa: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
  nPoints: number = 51,
): { x1: number; y1: number }[] {
  const txy = txyDiagram(compounds, P_Pa, fluidPackage, interactionParams, nPoints);
  return txy.map(p => ({ x1: p.x1, y1: p.y1 }));
}

// ────────────────────────────────────────────────────────────────
// Simplified Steam Properties (IAPWS-IF97 Region 1 & 2 approximation)
// ────────────────────────────────────────────────────────────────

/**
 * Steam saturation pressure (Pa) at temperature T_K.
 * Wagner equation (IAPWS-IF97 boundary equation).
 */
export function steamPsat_Pa(T_K: number): number {
  const Tc = 647.096;
  const Pc = 22064000; // Pa
  if (T_K >= Tc) return Pc;
  if (T_K < 273.15) return NaN;

  const tau = 1 - T_K / Tc;
  // Wagner equation coefficients for water
  const a1 = -7.85951783, a2 = 1.84408259, a3 = -11.7866497;
  const a4 = 22.6807411, a5 = -15.9618719, a6 = 1.80122502;

  const lnPr = (Tc / T_K) * (a1 * tau + a2 * Math.pow(tau, 1.5)
    + a3 * Math.pow(tau, 3) + a4 * Math.pow(tau, 3.5)
    + a5 * Math.pow(tau, 4) + a6 * Math.pow(tau, 7.5));

  return Pc * Math.exp(lnPr);
}

/**
 * Steam saturation temperature (K) at pressure P_Pa.
 * Inverse of Wagner equation via Newton-Raphson.
 */
export function steamTsat_K(P_Pa: number): number {
  const Pc = 22064000;
  if (P_Pa >= Pc) return 647.096;
  if (P_Pa <= 0) return NaN;

  // Initial guess from simplified Clausius-Clapeyron
  let T = 373.15 * (1 + 0.167 * Math.log(P_Pa / 101325));
  T = Math.max(273.15, Math.min(647, T));

  for (let iter = 0; iter < 50; iter++) {
    const Ps = steamPsat_Pa(T);
    const err = Ps - P_Pa;
    if (Math.abs(err / P_Pa) < 1e-8) return T;

    const dT = 0.01;
    const dPsdT = (steamPsat_Pa(T + dT) - Ps) / dT;
    if (Math.abs(dPsdT) < 1e-10) break;
    T -= err / dPsdT;
    T = Math.max(273.15, Math.min(647, T));
  }
  return T;
}

// ────────────────────────────────────────────────────────────────
// DIPPR 114 Liquid Heat Capacity (Zabransky-Ruzicka)
// ────────────────────────────────────────────────────────────────

/** DIPPR 114 coefficients for liquid heat capacity. */
export interface DIPPR114Coeffs {
  A: number;
  B: number;
  C: number;
  D: number;
  Tmin_K: number;
  Tmax_K: number;
}

/**
 * Liquid heat capacity (J/mol·K) from DIPPR 114 equation:
 *   Cp_L / R = A²/τ + B - 2·A·C·τ - A·D·τ² - C²·τ³/3 - C·D·τ⁴/2 - D²·τ⁵/5
 * where τ = 1 - T/Tc.
 * Note: This uses the form common in DIPPR 801 databases.
 */
export function liquidCp_DIPPR114(
  coeffs: DIPPR114Coeffs, Tc_K: number, T_K: number,
): number {
  const tau = 1 - T_K / Tc_K;
  if (tau <= 0) return NaN; // above critical
  const { A, B, C, D } = coeffs;
  // DIPPR 114 actual form: Cp = A^2/t + B - 2ACt - ADt^2 - C^2*t^3/3 - CD*t^4/2 - D^2*t^5/5
  const Cp = A * A / tau + B
    - 2 * A * C * tau
    - A * D * tau * tau
    - C * C * tau * tau * tau / 3
    - C * D * Math.pow(tau, 4) / 2
    - D * D * Math.pow(tau, 5) / 5;
  return Cp; // J/(mol·K)
}

// ────────────────────────────────────────────────────────────────
// Thermal Expansion & Compressibility
// ────────────────────────────────────────────────────────────────

/**
 * Isobaric thermal expansion coefficient β (1/K) for a liquid:
 *   β = (1/V)(∂V/∂T)_P  ≈ (V(T+δ) - V(T-δ))/(2·δ·V(T))
 * Uses Rackett correlation for molar volume.
 */
export function thermalExpansionCoeff(
  comp: CanopyCompound, T_K: number,
): number {
  const dT = 0.5;
  const V0 = liquidMolarVolume_m3pmol(comp, T_K);
  const Vp = liquidMolarVolume_m3pmol(comp, T_K + dT);
  const Vm = liquidMolarVolume_m3pmol(comp, T_K - dT);
  if (!isFinite(V0) || !isFinite(Vp) || !isFinite(Vm) || V0 <= 0) return NaN;
  return (Vp - Vm) / (2 * dT * V0);
}

/**
 * Isothermal compressibility κ_T (1/Pa) from EOS for a liquid:
 *   κ_T = -(1/V)(∂V/∂P)_T ≈ -(V(P+δ) - V(P-δ))/(2·δ·V(P))
 */
export function isothermalCompressibility(
  compounds: CanopyCompound[], x: number[], T_K: number, P_Pa: number,
  eos: 'PR' | 'SRK' = 'PR', kij?: number[][],
): number {
  const dP = P_Pa * 0.001;
  const V0 = liquidMolarVolume_EOS(compounds, x, T_K, P_Pa, eos, kij);
  const Vp = liquidMolarVolume_EOS(compounds, x, T_K, P_Pa + dP, eos, kij);
  const Vm = liquidMolarVolume_EOS(compounds, x, T_K, P_Pa - dP, eos, kij);
  if (!isFinite(V0) || !isFinite(Vp) || !isFinite(Vm) || V0 <= 0) return NaN;
  return -(Vp - Vm) / (2 * dP * V0);
}

// ────────────────────────────────────────────────────────────────
// Ideal Gas Speed of Sound
// ────────────────────────────────────────────────────────────────

/**
 * Speed of sound in an ideal gas (m/s):
 *   c = √(γ·R·T / M)
 * where γ = Cp/(Cp - R), R = 8.314, M = molecular weight (kg/mol).
 */
export function speedOfSoundIG(
  comp: CanopyCompound, T_K: number,
): number {
  const Cp = CpIG_JmolK(comp, T_K); // J/(mol·K)
  if (!isFinite(Cp) || Cp <= R) return NaN;
  const gamma = Cp / (Cp - R);
  const MW_kgpmol = comp.molecularWeight / 1000;
  return Math.sqrt(gamma * R * T_K / MW_kgpmol);
}

/**
 * Mixture ideal gas speed of sound (m/s):
 *   Uses mixture Cp and mixture MW.
 */
export function mixtureSpeedOfSoundIG(
  compounds: CanopyCompound[], y: number[], T_K: number,
): number {
  let CpMix = 0, MWmix = 0;
  for (let i = 0; i < compounds.length; i++) {
    CpMix += y[i] * CpIG_JmolK(compounds[i], T_K);
    MWmix += y[i] * compounds[i].molecularWeight;
  }
  if (CpMix <= R) return NaN;
  const gamma = CpMix / (CpMix - R);
  return Math.sqrt(gamma * R * T_K / (MWmix / 1000));
}

// ────────────────────────────────────────────────────────────────
// Joule-Thomson Coefficient
// ────────────────────────────────────────────────────────────────

/**
 * Joule-Thomson coefficient μ_JT (K/Pa) for an ideal gas = 0.
 * For real gas with PR EOS:
 *   μ_JT = (T·(∂V/∂T)_P - V) / Cp
 * Computed numerically from EOS.
 */
export function jouleThomson_KperPa(
  compounds: CanopyCompound[], y: number[], T_K: number, P_Pa: number,
  eos: 'PR' | 'SRK' = 'PR', kij?: number[][],
): number {
  const dT = 0.1;
  const fn = eos === 'PR' ? fugacityPR : fugacitySRK;

  const res = fn(compounds, y, T_K, P_Pa, true, kij);
  const V = res.Z * R * T_K / P_Pa;

  const resP = fn(compounds, y, T_K + dT, P_Pa, true, kij);
  const Vp = resP.Z * R * (T_K + dT) / P_Pa;

  const resM = fn(compounds, y, T_K - dT, P_Pa, true, kij);
  const Vm = resM.Z * R * (T_K - dT) / P_Pa;

  const dVdT = (Vp - Vm) / (2 * dT);

  // Mixture Cp
  let CpMix = 0;
  for (let i = 0; i < compounds.length; i++) {
    CpMix += y[i] * CpIG_JmolK(compounds[i], T_K);
  }

  if (CpMix <= 0) return NaN;
  return (T_K * dVdT - V) / CpMix;
}

// ────────────────────────────────────────────────────────────────
// Liquid Activity of Water
// ────────────────────────────────────────────────────────────────

/**
 * Water activity a_w = γ_w · x_w, used in food science and biochemistry.
 * @param gamma_w Activity coefficient of water
 * @param x_w Mole fraction of water
 */
export function waterActivity(gamma_w: number, x_w: number): number {
  return gamma_w * x_w;
}

// ────────────────────────────────────────────────────────────────
// Poynting Correction Factor
// ────────────────────────────────────────────────────────────────

/**
 * Poynting correction for liquid fugacity at elevated pressure:
 *   POY = exp(V_L · (P - Psat) / (R · T))
 * Important for high-pressure VLE.
 */
export function poyntingFactor(
  V_L_m3mol: number, P_Pa: number, Psat_Pa: number, T_K: number,
): number {
  return Math.exp((V_L_m3mol * (P_Pa - Psat_Pa)) / (R * T_K));
}

// ────────────────────────────────────────────────────────────────
// Heat of Mixing (from excess Gibbs energy)
// ────────────────────────────────────────────────────────────────

/**
 * Excess Gibbs energy (J/mol) from activity coefficients:
 *   G_E = R·T·Σ xi·ln(γi)
 */
export function excessGibbs(x: number[], T_K: number, gamma: number[]): number {
  let gE = 0;
  for (let i = 0; i < x.length; i++) {
    if (x[i] > 1e-30) {
      gE += x[i] * Math.log(gamma[i]);
    }
  }
  return R * T_K * gE;
}

/**
 * Excess entropy (J/mol·K) by numerical differentiation:
 *   S_E = -(∂G_E/∂T)_P,x
 */
export function excessEntropy(
  x: number[], T_K: number,
  gammaFn: (T: number) => number[],
): number {
  const dT = 0.1;
  const gE_p = excessGibbs(x, T_K + dT, gammaFn(T_K + dT));
  const gE_m = excessGibbs(x, T_K - dT, gammaFn(T_K - dT));
  return -(gE_p - gE_m) / (2 * dT);
}

// ────────────────────────────────────────────────────────────────
// Peng-Robinson Volume Translation (Peneloux)
// ────────────────────────────────────────────────────────────────

/**
 * Peneloux volume translation parameter for Peng-Robinson EOS (m³/mol):
 *   c_i = 0.40768 · R·Tc / Pc · (0.29441 - ZRA_i)
 * where ZRA = 0.29056 - 0.08775·ω.
 */
export function volumeTranslationPR(
  Tc_K: number, Pc_Pa: number, omega: number,
): number {
  const ZRA = 0.29056 - 0.08775 * omega;
  return 0.40768 * R * Tc_K / Pc_Pa * (0.29441 - ZRA);
}

/**
 * Apply Peneloux volume translation to EOS volume:
 *   V_corrected = V_EOS - Σ zi·ci
 */
export function applyVolumeTranslation(
  V_EOS_m3pmol: number, z: number[], c: number[],
): number {
  let cMix = 0;
  for (let i = 0; i < z.length; i++) cMix += z[i] * c[i];
  return V_EOS_m3pmol - cMix;
}

// ────────────────────────────────────────────────────────────────
// Wong-Sandler Mixing Rule
// ────────────────────────────────────────────────────────────────

/**
 * Wong-Sandler mixing rule: combines EOS with excess Gibbs models.
 * Returns mixture a and b parameters for PR EOS.
 *
 * b_mix = (Σi Σj zi·zj·(b - a/RT)_ij) / (1 - A_E/(C_EOS·RT) - Σi zi·ai/(bi·RT))
 * a_mix = b_mix · [Σi zi·ai/bi + A_E/C_EOS]
 *
 * where C_EOS = ln(√2 - 1)/√2 ≈ -0.6232 for PR
 */
export function wongSandlerMixingPR(
  a: number[], b: number[], z: number[], T_K: number,
  gE_RToverRT: number,  // G_E / (R·T) from activity model
  kij?: number[][],
): { a_mix: number; b_mix: number } {
  const n = z.length;
  const C_PR = Math.log(Math.sqrt(2) - 1) / Math.sqrt(2);

  // Cross terms: (b - a/RT)_ij = 0.5·[(b_i - a_i/RT) + (b_j - a_j/RT)]·(1 - k_ij)
  let sumBij = 0;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const bii = b[i] - a[i] / (R * T_K);
      const bjj = b[j] - a[j] / (R * T_K);
      const kijVal = kij?.[i]?.[j] ?? 0;
      sumBij += z[i] * z[j] * 0.5 * (bii + bjj) * (1 - kijVal);
    }
  }

  let sumAiOverBi = 0;
  for (let i = 0; i < n; i++) {
    sumAiOverBi += z[i] * a[i] / b[i];
  }

  const denom = 1 - gE_RToverRT / C_PR - sumAiOverBi / (R * T_K);
  if (Math.abs(denom) < 1e-15) return { a_mix: 0, b_mix: 0 };

  const b_mix = sumBij / denom;
  const a_mix = b_mix * (sumAiOverBi + gE_RToverRT * R * T_K / C_PR);

  return { a_mix: Math.max(a_mix, 0), b_mix: Math.max(b_mix, 1e-30) };
}

// ────────────────────────────────────────────────────────────────
// Huron-Vidal Mixing Rule
// ────────────────────────────────────────────────────────────────

/**
 * Modified Huron-Vidal first-order (MHV1) mixing rule:
 *   a_mix / (b_mix · RT) = Σi zi·ai/(bi·RT) + (1/q1)·[G_E/(RT) + Σi zi·ln(b_mix/bi)]
 *
 * where q1 = -0.53 for SRK, -0.593 for PR.
 */
export function huronVidalMixingPR(
  a: number[], b: number[], z: number[], T_K: number,
  gE_over_RT: number,
): { a_mix: number; b_mix: number } {
  const n = z.length;
  const q1 = -0.593; // for PR

  // b_mix = Σ zi·bi (linear)
  let b_mix = 0;
  for (let i = 0; i < n; i++) b_mix += z[i] * b[i];

  // Σ zi·ai/(bi·RT)
  let sumAiOverBiRT = 0;
  for (let i = 0; i < n; i++) sumAiOverBiRT += z[i] * a[i] / (b[i] * R * T_K);

  // Σ zi·ln(b_mix/bi)
  let sumLnB = 0;
  for (let i = 0; i < n; i++) sumLnB += z[i] * Math.log(b_mix / b[i]);

  const a_over_bRT = sumAiOverBiRT + (1 / q1) * (gE_over_RT + sumLnB);
  const a_mix = a_over_bRT * b_mix * R * T_K;

  return { a_mix: Math.max(a_mix, 0), b_mix };
}

// ────────────────────────────────────────────────────────────────
// Solid-Liquid Equilibrium (SLE)
// ────────────────────────────────────────────────────────────────

/**
 * Ideal solubility (mole fraction) of solid in liquid — simplified van't Hoff:
 *   ln(x_i) = (ΔH_fus / R) · (1/Tm - 1/T)
 *
 * @param deltaHfus_Jmol Enthalpy of fusion (J/mol)
 * @param Tm_K Melting point (K)
 * @param T_K Solution temperature (K)
 */
export function idealSolubility(
  deltaHfus_Jmol: number, Tm_K: number, T_K: number,
): number {
  if (T_K >= Tm_K) return 1; // fully melted
  const lnx = (deltaHfus_Jmol / R) * (1 / Tm_K - 1 / T_K);
  return Math.exp(lnx);
}

/**
 * SLE with activity coefficients:
 *   x_i = (1 / γ_i) · exp[(ΔH_fus / R) · (1/Tm - 1/T)]
 */
export function solubilitySLE(
  deltaHfus_Jmol: number, Tm_K: number, T_K: number, gamma_i: number,
): number {
  if (T_K >= Tm_K) return 1;
  const xIdeal = idealSolubility(deltaHfus_Jmol, Tm_K, T_K);
  return Math.min(1, xIdeal / gamma_i);
}

// ────────────────────────────────────────────────────────────────
// Enthalpy of Formation & Reaction
// ────────────────────────────────────────────────────────────────

/**
 * Standard enthalpy of reaction (J/mol) at T using Kirchhoff's law:
 *   ΔH_rxn(T) = ΔH_rxn(298.15) + ∫[298.15→T] ΔCp dT
 *
 * @param deltaHf_298 Array of standard enthalpies of formation (J/mol) per compound
 * @param stoich Array of stoichiometric coefficients (negative for reactants)
 * @param compounds Array of CanopyCompound objects
 * @param T_K Temperature (K)
 */
export function enthalpyOfReaction(
  deltaHf_298: number[], stoich: number[],
  compounds: CanopyCompound[], T_K: number,
): number {
  // ΔH_rxn at 298.15 K
  let dH298 = 0;
  for (let i = 0; i < stoich.length; i++) dH298 += stoich[i] * deltaHf_298[i];

  // ΔCp integration from T_REF to T
  let dCpIntegral = 0;
  for (let i = 0; i < compounds.length; i++) {
    const H_T = enthalpyIG(compounds[i], T_K);
    const H_ref = enthalpyIG(compounds[i], T_REF);
    dCpIntegral += stoich[i] * (H_T - H_ref);
  }

  return dH298 + dCpIntegral;
}

/**
 * Equilibrium constant K_eq from standard Gibbs energy of reaction:
 *   K_eq = exp(-ΔG_rxn / (R·T))
 *
 * @param deltaGf_298 Standard Gibbs energies of formation (J/mol)
 * @param stoich Stoichiometric coefficients
 * @param T_K Temperature (K)
 */
export function equilibriumConstant(
  deltaGf_298: number[], stoich: number[], T_K: number,
): number {
  let dG298 = 0;
  for (let i = 0; i < stoich.length; i++) dG298 += stoich[i] * deltaGf_298[i];
  return Math.exp(-dG298 / (R * T_K));
}

/**
 * van't Hoff equation — temperature-dependent equilibrium constant:
 *   ln(K(T2)/K(T1)) = (ΔH_rxn / R) · (1/T1 - 1/T2)
 */
export function vantHoffK(
  K_ref: number, deltaH_rxn: number,
  T_ref_K: number, T_K: number,
): number {
  return K_ref * Math.exp((deltaH_rxn / R) * (1 / T_ref_K - 1 / T_K));
}

// ────────────────────────────────────────────────────────────────
// Liquid-Liquid Equilibrium (LLE) Flash
// ────────────────────────────────────────────────────────────────

/**
 * LLE flash: finds compositions of two liquid phases at equilibrium.
 * Uses the isoactivity condition: γ_I·x_I = γ_II·x_II for each component.
 *
 * @param gammaFn Activity coeff function: (x[], T) → γ[]
 * @param z Feed composition
 * @param T_K Temperature
 * @returns { x1, x2, beta } where beta = fraction of phase II, or null if single phase
 */
export function flashLLE(
  z: number[], T_K: number,
  gammaFn: (x: number[], T: number) => number[],
): { x1: number[]; x2: number[]; beta: number } | null {
  const n = z.length;
  if (n < 2) return null;

  // Initial guess: perturb feed composition
  const x1 = z.map((zi, i) => i === 0 ? Math.min(0.95, zi * 1.5) : zi * 0.5);
  const x2 = z.map((zi, i) => i === 0 ? Math.max(0.05, zi * 0.5) : zi * 1.5);

  // Normalize
  const norm = (arr: number[]) => {
    const s = arr.reduce((a, b) => a + b, 0);
    return arr.map(v => v / s);
  };
  let xI = norm(x1);
  let xII = norm(x2);
  let beta = 0.5;

  for (let iter = 0; iter < 200; iter++) {
    const gI = gammaFn(xI, T_K);
    const gII = gammaFn(xII, T_K);

    // Distribution coefficient: K_i = γ_I_i / γ_II_i
    const K = gI.map((gi, i) => gi / gII[i]);

    // Rachford-Rice for LLE
    // f(β) = Σ zi·(Ki - 1) / (1 + β·(Ki - 1))
    const f = (b: number) => {
      let s = 0;
      for (let i = 0; i < n; i++) s += z[i] * (K[i] - 1) / (1 + b * (K[i] - 1));
      return s;
    };

    // Solve for beta by bisection
    let lo = 0, hi = 1;
    for (let bi = 0; bi < 50; bi++) {
      const mid = (lo + hi) / 2;
      if (f(mid) > 0) lo = mid; else hi = mid;
    }
    beta = (lo + hi) / 2;

    if (beta < 1e-8 || beta > 1 - 1e-8) return null; // single phase

    // Update compositions
    const newXI: number[] = [];
    const newXII: number[] = [];
    for (let i = 0; i < n; i++) {
      newXII[i] = z[i] / (1 + beta * (K[i] - 1));
      newXI[i] = K[i] * newXII[i];
    }

    // Check convergence
    let maxDiff = 0;
    for (let i = 0; i < n; i++) {
      maxDiff = Math.max(maxDiff, Math.abs(newXI[i] - xI[i]), Math.abs(newXII[i] - xII[i]));
    }

    xI = norm(newXI);
    xII = norm(newXII);

    if (maxDiff < 1e-10) return { x1: xI, x2: xII, beta };
  }

  return { x1: xI, x2: xII, beta };
}

// ────────────────────────────────────────────────────────────────
// Acentric Factor Estimation
// ────────────────────────────────────────────────────────────────

/**
 * Estimate acentric factor from Tb and Tc/Pc (Lee-Kesler):
 *   ω = [-ln(Pc/1.01325) - 5.92714 + 6.09648/Tbr + 1.28862·ln(Tbr) - 0.169347·Tbr⁶]
 *       / [15.2518 - 15.6875/Tbr - 13.4721·ln(Tbr) + 0.43577·Tbr⁶]
 * where Tbr = Tb/Tc, Pc in bar.
 */
export function estimateOmega_LeeKesler(
  Tc_K: number, Pc_Pa: number, Tb_K: number,
): number {
  const Pc_bar = Pc_Pa / 1e5;
  const Tbr = Tb_K / Tc_K;
  const Tbr6 = Math.pow(Tbr, 6);
  const num = -Math.log(Pc_bar / 1.01325)
    - 5.92714 + 6.09648 / Tbr + 1.28862 * Math.log(Tbr) - 0.169347 * Tbr6;
  const den = 15.2518 - 15.6875 / Tbr - 13.4721 * Math.log(Tbr) + 0.43577 * Tbr6;
  return num / den;
}

// ────────────────────────────────────────────────────────────────
// Critical Property Estimation (Joback-Reid)
// ────────────────────────────────────────────────────────────────

/** Joback group contribution for critical properties. */
export interface JobackGroup {
  Tc: number;   // contribution to Tc (K)
  Pc: number;   // contribution to Pc (bar)
  Vc: number;   // contribution to Vc (cm³/mol)
  Tb: number;   // contribution to Tb (K)
  nAtoms: number; // number of atoms
}

/**
 * Estimate Tc (K) from Joback method:
 *   Tc = Tb / [0.584 + 0.965·Σ(ΔTc) - (Σ(ΔTc))²]
 */
export function jobeckTc(
  Tb_K: number, groupContribs: number[],
): number {
  const sum = groupContribs.reduce((a, b) => a + b, 0);
  return Tb_K / (0.584 + 0.965 * sum - sum * sum);
}

/**
 * Estimate Pc (Pa) from Joback method:
 *   Pc = 1/(0.113 + 0.0032·n_atoms - Σ(ΔPc))²  (bar)
 */
export function jobeckPc(
  nAtoms: number, groupContribs: number[],
): number {
  const sum = groupContribs.reduce((a, b) => a + b, 0);
  const Pc_bar = 1 / Math.pow(0.113 + 0.0032 * nAtoms - sum, 2);
  return Pc_bar * 1e5;  // convert to Pa
}

/**
 * Estimate Tb (K) from Joback method:
 *   Tb = 198.2 + Σ(ΔTb)
 */
export function jobeckTb(groupContribs: number[]): number {
  return 198.2 + groupContribs.reduce((a, b) => a + b, 0);
}

// ────────────────────────────────────────────────────────────────
// Relative Volatility
// ────────────────────────────────────────────────────────────────

/**
 * Relative volatility of component i vs reference (last component):
 *   α_i = K_i / K_ref
 */
export function relativeVolatility(
  compounds: CanopyCompound[], x: number[], y: number[],
  T_K: number, P_Pa: number,
  fluidPackage: FluidPackageType = 'Ideal',
  interactionParams?: InteractionParams,
  refIndex?: number,
): number[] {
  const K = computeKvalues(compounds, x, y, T_K, P_Pa, fluidPackage, interactionParams);
  const ref = refIndex ?? compounds.length - 1;
  const Kref = K[ref];
  if (Kref === 0) return K.map(() => Infinity);
  return K.map(Ki => Ki / Kref);
}

// ────────────────────────────────────────────────────────────────
// Minimum Reflux (Underwood Method)
// ────────────────────────────────────────────────────────────────

/**
 * Underwood equation for minimum reflux in binary distillation:
 *   R_min = (1/(α-1)) · [xD/xF - α·(1-xD)/(1-xF)]
 *
 * @param alpha Relative volatility
 * @param xF Feed composition of light key
 * @param xD Distillate composition of light key
 * @param q Feed quality (1 = saturated liquid, 0 = saturated vapor)
 */
export function underwoodRmin(
  alpha: number, xF: number, xD: number, q: number,
): number {
  // Underwood root θ: Σ αi·zi/(αi - θ) = 1 - q
  // For binary: α·xF/(α - θ) + (1-xF)/(1 - θ) = 1 - q
  // Solve for θ by bisection between 1 and α
  let lo = 1.001, hi = alpha - 0.001;
  for (let iter = 0; iter < 100; iter++) {
    const mid = (lo + hi) / 2;
    const f = alpha * xF / (alpha - mid) + (1 - xF) / (1 - mid) - (1 - q);
    if (f > 0) lo = mid; else hi = mid;
  }
  const theta = (lo + hi) / 2;

  // R_min + 1 = Σ αi·xD_i / (αi - θ)
  const RplusOne = alpha * xD / (alpha - theta) + (1 - xD) / (1 - theta);
  return RplusOne - 1;
}

/**
 * Fenske equation — minimum number of stages:
 *   N_min = ln[(xD/(1-xD))·((1-xB)/xB)] / ln(α_avg)
 */
export function fenskeNmin(
  alpha_avg: number, xD: number, xB: number,
): number {
  return Math.log((xD / (1 - xD)) * ((1 - xB) / xB)) / Math.log(alpha_avg);
}

/**
 * Gilliland correlation — actual stages from R and R_min:
 *   X = (R - R_min) / (R + 1)
 *   Y = 1 - exp[(1 + 54.4X) / (11 + 117.2X) · (X-1) / √X]
 *   N = (N_min + Y) / (1 - Y)
 */
export function gillilandN(
  Nmin: number, R: number, Rmin: number,
): number {
  const X = (R - Rmin) / (R + 1);
  if (X <= 0) return Infinity;
  const Y = 1 - Math.exp(
    ((1 + 54.4 * X) / (11 + 117.2 * X)) * ((X - 1) / Math.sqrt(X))
  );
  return (Nmin + Y) / (1 - Y);
}

// ────────────────────────────────────────────────────────────────
// Heat Exchanger Design (LMTD & NTU)
// ────────────────────────────────────────────────────────────────

/**
 * Log mean temperature difference (K):
 *   LMTD = (ΔT1 - ΔT2) / ln(ΔT1/ΔT2)
 */
export function LMTD(dT1: number, dT2: number): number {
  if (dT1 <= 0 || dT2 <= 0) return NaN;
  if (Math.abs(dT1 - dT2) < 1e-6) return dT1;
  return (dT1 - dT2) / Math.log(dT1 / dT2);
}

/**
 * LMTD correction factor for multi-pass exchangers (1 shell / 2 tube):
 *   F(P, R) from empirical correlation.
 *   P = (T_co - T_ci) / (T_hi - T_ci)
 *   R = (T_hi - T_ho) / (T_co - T_ci)
 */
export function LMTD_F_factor(
  T_hi: number, T_ho: number, T_ci: number, T_co: number,
): number {
  const P_ = (T_co - T_ci) / (T_hi - T_ci);
  const R_ = (T_hi - T_ho) / (T_co - T_ci);

  if (Math.abs(R_ - 1) < 1e-6) {
    // Special case R = 1
    const S = P_ / (1 - P_);
    if (S <= 0) return NaN;
    return (Math.sqrt(2) * S) / ((1 - S) * Math.log((S * Math.sqrt(2) + S) / (S * Math.sqrt(2) - S + 2)));
  }

  const S = Math.sqrt(R_ * R_ + 1) / (R_ - 1);
  const W = ((1 - P_ * R_) / (1 - P_));
  if (W <= 0) return NaN;
  const num = S * Math.log(W);
  const denom = Math.log((2 - P_ * (R_ + 1 - S)) / (2 - P_ * (R_ + 1 + S)));
  if (Math.abs(denom) < 1e-15) return NaN;
  return num / denom;
}

/**
 * NTU-effectiveness method for heat exchanger sizing.
 * ε = f(NTU, C_min/C_max) for counterflow:
 *   ε = [1 - exp(-NTU·(1-Cr))] / [1 - Cr·exp(-NTU·(1-Cr))]
 */
export function hxEffectivenessCounterflow(
  NTU: number, Cr: number,
): number {
  if (Cr < 1e-10) return 1 - Math.exp(-NTU);
  if (Math.abs(Cr - 1) < 1e-6) return NTU / (1 + NTU);
  const expTerm = Math.exp(-NTU * (1 - Cr));
  return (1 - expTerm) / (1 - Cr * expTerm);
}

/**
 * Required heat exchanger area (m²):
 *   A = Q / (U · F · LMTD)
 */
export function hxArea(
  Q_W: number, U_Wpm2K: number, F: number, lmtd: number,
): number {
  if (U_Wpm2K <= 0 || F <= 0 || lmtd <= 0) return NaN;
  return Q_W / (U_Wpm2K * F * lmtd);
}

// ────────────────────────────────────────────────────────────────
// Pipe Pressure Drop (Darcy-Weisbach)
// ────────────────────────────────────────────────────────────────

/**
 * Reynolds number for pipe flow:
 *   Re = ρ·v·D / μ
 */
export function reynoldsNumber(
  density_kgm3: number, velocity_ms: number,
  diameter_m: number, viscosity_Pas: number,
): number {
  return density_kgm3 * velocity_ms * diameter_m / viscosity_Pas;
}

/**
 * Darcy friction factor from Churchill equation (all flow regimes):
 *   f = 8·[(8/Re)^12 + (A+B)^(-3/2)]^(1/12)
 */
export function frictionFactorChurchill(
  Re: number, roughness_m: number, diameter_m: number,
): number {
  if (Re <= 0) return NaN;
  const epsD = roughness_m / diameter_m;
  const A = Math.pow(
    -2.457 * Math.log(Math.pow(7 / Re, 0.9) + 0.27 * epsD),
    16,
  );
  const B = Math.pow(37530 / Re, 16);
  return 8 * Math.pow(
    Math.pow(8 / Re, 12) + Math.pow(A + B, -1.5),
    1 / 12,
  );
}

/**
 * Darcy-Weisbach pressure drop (Pa):
 *   ΔP = f · (L/D) · (ρ·v²/2)
 */
export function pressureDropPipe(
  f: number, length_m: number, diameter_m: number,
  density_kgm3: number, velocity_ms: number,
): number {
  return f * (length_m / diameter_m) * (density_kgm3 * velocity_ms * velocity_ms / 2);
}

/**
 * Pressure drop across fittings (K-factor method):
 *   ΔP = K · (ρ·v²/2)
 */
export function pressureDropFitting(
  K: number, density_kgm3: number, velocity_ms: number,
): number {
  return K * (density_kgm3 * velocity_ms * velocity_ms / 2);
}

// ────────────────────────────────────────────────────────────────
// Pump & Compressor
// ────────────────────────────────────────────────────────────────

/**
 * Pump power (W) from head and flow:
 *   W = ρ·g·H·Q / η
 *
 * @param head_m Head in meters
 * @param flowRate_m3ps Volumetric flow rate m³/s
 * @param density_kgm3 Liquid density
 * @param efficiency Pump efficiency (0-1)
 */
export function pumpPower_W(
  head_m: number, flowRate_m3ps: number,
  density_kgm3: number, efficiency: number,
): number {
  const g = 9.80665;
  return density_kgm3 * g * head_m * flowRate_m3ps / efficiency;
}

/**
 * NPSH available (m):
 *   NPSH_a = (P_suction - P_vapor) / (ρ·g) + z_suction
 */
export function NPSHavailable(
  P_suction_Pa: number, P_vapor_Pa: number,
  density_kgm3: number, z_suction_m: number = 0,
): number {
  const g = 9.80665;
  return (P_suction_Pa - P_vapor_Pa) / (density_kgm3 * g) + z_suction_m;
}

/**
 * Polytropic compressor work (J/mol):
 *   W_poly = Z·R·T₁/(η_p·(n-1)/n) · [(P₂/P₁)^((n-1)/n) - 1]
 * where n = k/(k - η_p·(k-1)) polytropic exponent, k = Cp/Cv.
 */
export function polytropicWork_Jmol(
  T_in_K: number, P_in_Pa: number, P_out_Pa: number,
  k: number, eta_poly: number, Z: number = 1,
): number {
  if (eta_poly <= 0 || k <= 1) return NaN;
  const n = k / (k - eta_poly * (k - 1));
  const exp_ = (n - 1) / n;
  return Z * R * T_in_K / (eta_poly * exp_) * (Math.pow(P_out_Pa / P_in_Pa, exp_) - 1);
}

/**
 * Adiabatic compressor discharge temperature (K):
 *   T₂ = T₁ · (P₂/P₁)^((k-1)/k)
 */
export function adiabaticDischargeT(
  T_in_K: number, P_in_Pa: number, P_out_Pa: number, k: number,
): number {
  return T_in_K * Math.pow(P_out_Pa / P_in_Pa, (k - 1) / k);
}

// ────────────────────────────────────────────────────────────────
// Orifice Plate & Control Valve
// ────────────────────────────────────────────────────────────────

/**
 * Orifice plate flow rate (kg/s) — ISO 5167:
 *   ṁ = C_d · ε · A_orifice · √(2·ρ·ΔP)
 *
 * @param Cd Discharge coefficient (~0.6)
 * @param epsilon Expansibility factor (~1 for liquids)
 * @param d_orifice Orifice diameter (m)
 * @param density_kgm3 Upstream density
 * @param dP_Pa Differential pressure
 */
export function orificeFlowRate_kgps(
  Cd: number, epsilon: number, d_orifice: number,
  density_kgm3: number, dP_Pa: number,
): number {
  const A = Math.PI * d_orifice * d_orifice / 4;
  return Cd * epsilon * A * Math.sqrt(2 * density_kgm3 * Math.abs(dP_Pa));
}

/**
 * Control valve Cv from flow:
 *   Cv = Q · √(SG / ΔP)
 * where Q in US gpm, ΔP in psi, SG relative to water.
 *
 * This version uses SI units:
 *   Q in m³/h, ΔP in Pa, returns Cv in US gpm/√psi.
 */
export function controlValveCv(
  flowRate_m3ph: number, dP_Pa: number, SG: number = 1,
): number {
  // Convert to US gpm and psi
  const Q_gpm = flowRate_m3ph * 4.40287;   // m³/h to US gpm
  const dP_psi = dP_Pa / 6894.76;
  if (dP_psi <= 0) return NaN;
  return Q_gpm * Math.sqrt(SG / dP_psi);
}

// ────────────────────────────────────────────────────────────────
// Reactor Kinetics
// ────────────────────────────────────────────────────────────────

/**
 * Arrhenius rate constant:
 *   k = A · exp(-Ea / (R·T))
 */
export function arrheniusK(
  A: number, Ea_Jmol: number, T_K: number,
): number {
  return A * Math.exp(-Ea_Jmol / (R * T_K));
}

/**
 * Power-law reaction rate:
 *   r = k · Π(C_i ^ n_i)
 * @param k Rate constant
 * @param concentrations Molar concentrations (mol/m³)
 * @param orders Reaction orders for each species
 * @returns Rate (mol/(m³·s))
 */
export function powerLawRate(
  k: number, concentrations: number[], orders: number[],
): number {
  let rate = k;
  for (let i = 0; i < concentrations.length; i++) {
    rate *= Math.pow(Math.max(concentrations[i], 0), orders[i]);
  }
  return rate;
}

/**
 * CSTR design equation — required volume:
 *   V = F_A0 · X_A / (-r_A)
 *
 * @param FA0 Molar feed rate of A (mol/s)
 * @param XA Conversion
 * @param rateA Reaction rate of A at exit conditions (mol/(m³·s)), positive value
 */
export function cstrVolume_m3(
  FA0: number, XA: number, rateA: number,
): number {
  if (rateA <= 0) return NaN;
  return FA0 * XA / rateA;
}

/**
 * PFR design equation (simple first-order):
 *   V = F_A0 / k · ln(1/(1-X_A))
 */
export function pfrVolume_firstOrder_m3(
  FA0: number, XA: number, k: number, CA0: number,
): number {
  if (XA >= 1 || k <= 0) return NaN;
  return (FA0 / (k * CA0)) * Math.log(1 / (1 - XA));
}

// ────────────────────────────────────────────────────────────────
// Vessel Sizing (API-based)
// ────────────────────────────────────────────────────────────────

/**
 * Vertical separator liquid capacity (m³):
 *   V_liquid = Q_L · t_retention
 *
 * @param flowRate_m3ps Liquid volumetric flow rate
 * @param retentionTime_s Retention time in seconds
 */
export function separatorLiquidVolume(
  flowRate_m3ps: number, retentionTime_s: number,
): number {
  return flowRate_m3ps * retentionTime_s;
}

/**
 * Souders-Brown maximum allowable vapor velocity (m/s):
 *   V_max = K_SB · √((ρ_L - ρ_V) / ρ_V)
 *
 * @param K_SB Souders-Brown factor (typically 0.03-0.07 m/s)
 */
export function soudersBrownVelocity(
  K_SB: number, rhoL_kgm3: number, rhoV_kgm3: number,
): number {
  return K_SB * Math.sqrt((rhoL_kgm3 - rhoV_kgm3) / rhoV_kgm3);
}

/**
 * Minimum vessel diameter for vapor disengagement (m):
 *   D = √(4·Q_V / (π·V_max))
 */
export function minVesselDiameter(
  Q_V_m3ps: number, V_max_ms: number,
): number {
  return Math.sqrt(4 * Q_V_m3ps / (Math.PI * V_max_ms));
}

// ────────────────────────────────────────────────────────────────
// Pressure Vessel Thickness (ASME VIII-1)
// ────────────────────────────────────────────────────────────────

/**
 * Required wall thickness for cylindrical shell under internal pressure:
 *   t = P·R / (S·E - 0.6·P)
 *
 * @param P_Pa Internal gauge pressure (Pa)
 * @param R_m Inside radius (m)
 * @param S_Pa Maximum allowable stress (Pa)
 * @param E Joint efficiency (0-1, typically 0.85-1.0)
 */
export function shellThickness_m(
  P_Pa: number, R_m: number, S_Pa: number, E: number = 1,
): number {
  const denom = S_Pa * E - 0.6 * P_Pa;
  if (denom <= 0) return NaN;
  return P_Pa * R_m / denom;
}

/**
 * 2:1 Ellipsoidal head thickness (m):
 *   t = P·D / (2·S·E - 0.2·P)
 */
export function headThickness_m(
  P_Pa: number, D_m: number, S_Pa: number, E: number = 1,
): number {
  const denom = 2 * S_Pa * E - 0.2 * P_Pa;
  if (denom <= 0) return NaN;
  return P_Pa * D_m / denom;
}

// ────────────────────────────────────────────────────────────────
// Wagner Vapor Pressure (DIPPR 115 / Wagner 25 form)
// ────────────────────────────────────────────────────────────────

/**
 * Wagner 25 equation:  ln(P/Pc) = (a·τ + b·τ^1.5 + c·τ^3 + d·τ^6) / Tr
 * where τ = 1 - Tr, Tr = T/Tc.  Highly accurate up to the critical point.
 * RPN from pplcdefs (WAGNER): [22] [43] / 1 - ...
 */
export function PsatWagner_Pa(comp: CanopyCompound, T_K: number): number {
  const w = comp.wagnerVP;
  if (!w || T_K <= 0) return NaN;
  const Tr = T_K / comp.Tc_K;
  if (Tr >= 1) return comp.Pc_Pa; // At/above critical
  const tau = 1 - Tr;
  const lnPr = (w.a * tau + w.b * Math.pow(tau, 1.5) + w.c * tau * tau * tau
    + w.d * Math.pow(tau, 6)) / Tr;
  const P = comp.Pc_Pa * Math.exp(lnPr);
  return isFinite(P) && P > 0 ? P : NaN;
}

/**
 * Best available vapor pressure: Wagner → Extended Antoine → NaN.
 */
export function PsatBest_Pa(comp: CanopyCompound, T_K: number): number {
  if (comp.wagnerVP) {
    const w = PsatWagner_Pa(comp, T_K);
    if (isFinite(w) && w > 0) return w;
  }
  return Psat_Pa(comp, T_K);
}

// ────────────────────────────────────────────────────────────────
// Watson Heat of Vaporization Extrapolation (DHVLWT)
// ────────────────────────────────────────────────────────────────

/**
 * Watson correlation:  ΔHvap(T) = ΔHvap(T_ref) · [(1 - T/Tc) / (1 - T_ref/Tc)]^n
 * Default n = 0.38 (recommended by Watson; Aspen uses fitted n when available).
 * RPN from pplcdefs: [22] A / 1 - 0 1 B ? - 1 C ? - 1 D ? - * * * [22B] A / 1 - ... / ^ *
 */
export function HvapWatson_Jmol(comp: CanopyCompound, T_K: number): number {
  const w = comp.watsonHvap;
  if (!w) return NaN;
  const Tr = T_K / comp.Tc_K;
  const Tr_ref = w.T_ref_K / comp.Tc_K;
  if (Tr >= 1) return 0;
  if (Tr_ref >= 1) return NaN;
  const n = w.n ?? 0.38;
  return w.Hvap_ref_Jmol * Math.pow((1 - Tr) / (1 - Tr_ref), n);
}

/**
 * Best available Hvap: DIPPR 106 → Watson → Pitzer fallback.
 */
export function HvapBest_Jmol(comp: CanopyCompound, T_K: number): number {
  const v = Hvap_Jmol(comp, T_K);
  if (isFinite(v) && v > 0) return v;
  const w = HvapWatson_Jmol(comp, T_K);
  if (isFinite(w) && w > 0) return w;
  // Pitzer fallback (already in Hvap_Jkmol when coefficients are zero)
  return Hvap_Jkmol(comp, T_K) / 1000;
}

// ────────────────────────────────────────────────────────────────
// DIPPR 100 Universal Polynomial Evaluator
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR Equation 100:  f(T) = A + B·T + C·T² + D·T³ + E·T⁴
 * Used for: liquid density, liquid Cp, solid Cp, solid density, thermal conductivity, etc.
 */
export function DIPPR100(T: number, A: number, B: number, C: number, D: number, E: number = 0): number {
  return A + T * (B + T * (C + T * (D + T * E)));
}

/**
 * DIPPR 100 integral from T1 to T2:  ∫f(T)dT = A·ΔT + B/2·ΔT² + C/3·ΔT³ + D/4·ΔT⁴ + E/5·ΔT⁵
 */
export function DIPPR100_integral(T1: number, T2: number, A: number, B: number, C: number, D: number, E: number = 0): number {
  const f = (T: number) => A * T + B / 2 * T * T + C / 3 * T * T * T + D / 4 * T ** 4 + E / 5 * T ** 5;
  return f(T2) - f(T1);
}

// ────────────────────────────────────────────────────────────────
// DIPPR 116 Liquid Density (Near-Critical)
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR 116:  ρ = A + B·τ^0.35 + C·τ^(2/3) + D·τ + E·τ^(4/3)
 * where τ = 1 - T/Tc.  Accurate near the critical point where DIPPR 105 diverges.
 * Returns kmol/m³.
 */
export function liquidDensity_DIPPR116_kmolpm3(comp: CanopyCompound, T_K: number): number {
  const c = comp.liquidDensityDIPPR116;
  if (!c) return NaN;
  const Tr = T_K / comp.Tc_K;
  if (Tr >= 1) return NaN; // Above critical
  const tau = 1 - Tr;
  return c.A + c.B * Math.pow(tau, 0.35) + c.C * Math.pow(tau, 2 / 3) + c.D * tau + c.E * Math.pow(tau, 4 / 3);
}

/**
 * Best liquid density: DIPPR 116 → DIPPR 105 → Rackett fallback.
 */
export function liquidDensityBest_kmolpm3(comp: CanopyCompound, T_K: number): number {
  // Try DIPPR 116 first (best near critical)
  if (comp.liquidDensityDIPPR116) {
    const rho = liquidDensity_DIPPR116_kmolpm3(comp, T_K);
    if (isFinite(rho) && rho > 0) return rho;
  }
  // Fall back to DIPPR 105 (Rackett)
  return liquidDensity_kmolpm3(comp, T_K);
}

// ────────────────────────────────────────────────────────────────
// Aspen Polynomial Ideal Gas Cp (CPIG / 6-coefficient)
// ────────────────────────────────────────────────────────────────

/**
 * Aspen polynomial:  Cp_ig = A + B·T + C·T² + D·T³ + E·T⁴ + F·T⁵  [J/(kmol·K)]
 * Used by many Aspen databanks when DIPPR 107 (Aly-Lee) is not available.
 * pplcdefs CPIG line 1302: A B [22] * + C [22] 2 ^ * + D [22] 3 ^ * + E [22] 4 ^ * + F [22] 5 ^ * +
 */
export function CpIG_AspenPoly_JkmolK(comp: CanopyCompound, T_K: number): number {
  const c = comp.cpigPoly;
  if (!c) return NaN;
  return c.A + T_K * (c.B + T_K * (c.C + T_K * (c.D + T_K * (c.E + T_K * (c.F ?? 0)))));
}

/**
 * Best available CpIG: DIPPR 107 → Aspen polynomial → DIPPR 16.
 */
export function CpIG_Best_JmolK(comp: CanopyCompound, T_K: number): number {
  if (comp.cpigdp) return CpIG_DIPPR107_JkmolK(comp, T_K) / 1000;
  if (comp.cpigPoly) return CpIG_AspenPoly_JkmolK(comp, T_K) / 1000;
  return CpIG_JkmolK(comp, T_K) / 1000;
}

// ────────────────────────────────────────────────────────────────
// Modified Rackett Liquid Density (MRKZRA)
// ────────────────────────────────────────────────────────────────

/**
 * Modified Rackett equation:
 *   V_L = (R·Tc / Pc) · ZRA^(1 + (1-Tr)^(2/7))
 *
 * When ZRA is compound-specific (fitted to data), this is significantly more accurate
 * than the generic ZRA ≈ 0.29056 - 0.08775·ω used in the default Rackett.
 * pplcdefs MRKZRA line 12954.
 */
export function RackettDensity_kmolpm3(comp: CanopyCompound, T_K: number): number {
  const ZRA = comp.rackett_ZRA ?? (0.29056 - 0.08775 * comp.omega);
  const Tr = T_K / comp.Tc_K;
  if (Tr >= 1) return NaN;
  const Vm = (R * comp.Tc_K / comp.Pc_Pa) * Math.pow(ZRA, 1 + Math.pow(1 - Tr, 2 / 7));
  return 1 / (Vm * 1000); // Convert from m³/mol to kmol/m³
}

// ────────────────────────────────────────────────────────────────
// Clausius-Clapeyron Boiling Point Elevation / Depression
// ────────────────────────────────────────────────────────────────

/**
 * Clausius-Clapeyron:  ln(P2/P1) = ΔHvap/R · (1/T1 - 1/T2)
 * Given a reference (Pref, Tref), estimate boiling T at a new pressure.
 */
export function boilingPointAtP(comp: CanopyCompound, P_Pa: number, P_ref_Pa: number, T_ref_K: number): number {
  const dHvap = Hvap_Jmol(comp, T_ref_K);
  if (!isFinite(dHvap) || dHvap <= 0) return NaN;
  const inv_T = 1 / T_ref_K - R * Math.log(P_Pa / P_ref_Pa) / dHvap;
  return inv_T > 0 ? 1 / inv_T : NaN;
}

// ────────────────────────────────────────────────────────────────
// DIPPR 101 Evaluator (Liquid Viscosity, Vapor Pressure alt form)
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR 101:  f(T) = exp(A + B/T + C·ln(T) + D·T^E)
 * Used for liquid viscosity, some vapor pressure correlations.
 */
export function DIPPR101(T: number, A: number, B: number, C: number, D: number, E: number = 0): number {
  return Math.exp(A + B / T + C * Math.log(T) + D * Math.pow(T, E));
}

// ────────────────────────────────────────────────────────────────
// DIPPR 102 Evaluator (Vapor Viscosity, Vapor Thermal Conductivity)
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR 102:  f(T) = A·T^B / (1 + C/T + D/T²)
 * Used for vapor viscosity, vapor thermal conductivity.
 */
export function DIPPR102(T: number, A: number, B: number, C: number, D: number): number {
  return A * Math.pow(T, B) / (1 + C / T + D / (T * T));
}

// ────────────────────────────────────────────────────────────────
// DIPPR 104 Evaluator (Second Virial Coefficient)
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR 104:  f(T) = A + B/T + C/T³ + D/T⁸ + E/T⁹
 * Used for second virial coefficient B(T).
 */
export function DIPPR104(T: number, A: number, B: number, C: number, D: number, E: number = 0): number {
  const T2 = T * T;
  const T3 = T2 * T;
  return A + B / T + C / T3 + D / (T3 * T3 * T2) + E / (T3 * T3 * T3);
}

// ────────────────────────────────────────────────────────────────
// DIPPR 105 Evaluator (Liquid Density, saturated)
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR 105:  ρ = A / B^(1 + (1-T/C)^D)  [kmol/m³]
 * Standard DIPPR liquid density correlation.
 */
export function DIPPR105(T: number, A: number, B: number, C: number, D: number): number {
  const tau = 1 - T / C;
  if (tau <= 0) return NaN;
  return A / Math.pow(B, 1 + Math.pow(tau, D));
}

// ────────────────────────────────────────────────────────────────
// DIPPR 106 Evaluator (Heat of Vaporization, Surface Tension)
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR 106:  f(T) = A · (1 - Tr)^(B + C·Tr + D·Tr² + E·Tr³)
 * where Tr = T/Tc.  Used for ΔHvap, surface tension.
 */
export function DIPPR106(T: number, Tc: number, A: number, B: number, C: number, D: number, E: number = 0): number {
  const Tr = T / Tc;
  if (Tr >= 1) return 0;
  return A * Math.pow(1 - Tr, B + Tr * (C + Tr * (D + Tr * E)));
}

// ────────────────────────────────────────────────────────────────
// DIPPR 107 Evaluator (Aly-Lee Ideal Gas Cp)
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR 107 (Aly-Lee): Cp = A + B·[(C/T)/sinh(C/T)]² + D·[(E/T)/cosh(E/T)]²
 * Returns J/(kmol·K).
 */
export function DIPPR107(T: number, A: number, B: number, C: number, D: number, E: number): number {
  let val = A;
  if (Math.abs(C) > 1e-30) {
    const u = C / T;
    const sh = Math.sinh(u);
    val += B * (u / sh) * (u / sh);
  }
  if (Math.abs(E) > 1e-30) {
    const u = E / T;
    const ch = Math.cosh(u);
    val += D * (u / ch) * (u / ch);
  }
  return val;
}

// ────────────────────────────────────────────────────────────────
// Corresponding States Liquid Cp (Bondi-Rowlinson)
// ────────────────────────────────────────────────────────────────

/**
 * Bondi-Rowlinson corresponding-states correlation for liquid Cp:
 *   (CpL - CpIG) / R = 1.45 + 0.45/(1-Tr) + 0.25·ω·(17.11 + 25.2·(1-Tr)^(1/3)/Tr + 1.742/(1-Tr))
 *
 * Valid for non-polar and slightly polar liquids at Tr < 0.85.
 * Returns J/(mol·K).
 */
export function CpL_BondiRowlinson_JmolK(comp: CanopyCompound, T_K: number): number {
  const Tr = T_K / comp.Tc_K;
  if (Tr >= 0.99) return NaN;
  const CpIG_val = CpIG_JmolK(comp, T_K);
  const dCpR = 1.45 + 0.45 / (1 - Tr) + 0.25 * comp.omega
    * (17.11 + 25.2 * Math.pow(1 - Tr, 1 / 3) / Tr + 1.742 / (1 - Tr));
  return CpIG_val + R * dCpR;
}

// ────────────────────────────────────────────────────────────────
// Lee-Kesler Vapor Pressure Correlation
// ────────────────────────────────────────────────────────────────

/**
 * Lee-Kesler three-parameter correlation:
 *   ln(Pr) = f⁰(Tr) + ω·f¹(Tr)
 * where:
 *   f⁰ = 5.92714 - 6.09648/Tr - 1.28862·ln(Tr) + 0.169347·Tr⁶
 *   f¹ = 15.2518 - 15.6875/Tr - 13.4721·ln(Tr) + 0.43577·Tr⁶
 *
 * Useful fallback when no Antoine or Wagner data is available.
 */
export function PsatLeeKesler_Pa(comp: CanopyCompound, T_K: number): number {
  const Tr = T_K / comp.Tc_K;
  if (Tr >= 1) return comp.Pc_Pa;
  const Tr6 = Tr * Tr * Tr * Tr * Tr * Tr;
  const f0 = 5.92714 - 6.09648 / Tr - 1.28862 * Math.log(Tr) + 0.169347 * Tr6;
  const f1 = 15.2518 - 15.6875 / Tr - 13.4721 * Math.log(Tr) + 0.43577 * Tr6;
  const lnPr = f0 + comp.omega * f1;
  const P = comp.Pc_Pa * Math.exp(lnPr);
  return isFinite(P) && P > 0 ? P : NaN;
}

// ────────────────────────────────────────────────────────────────
// Compressed Liquid Density (Tait Equation)
// ────────────────────────────────────────────────────────────────

/**
 * Thomson-Brobst-Hankinson (Tait) equation for compressed liquid density:
 *   V(T,P) = Vs(T) · [1 - C·ln((B+P)/(B+Psat))]
 *
 * where Vs is the saturated liquid volume (e.g. from COSTALD or Rackett),
 * B and C are functions of Tr and ω, and Psat is the vapor pressure.
 */
export function taitCompressedVolume_m3pmol(
  comp: CanopyCompound, T_K: number, P_Pa: number, Vs_m3pmol: number,
): number {
  const Tr = T_K / comp.Tc_K;
  const omega = comp.omega;
  // Tait parameters (generalized)
  const e = Math.exp(1);
  const B = comp.Pc_Pa * (-1 - 9.070217 * (1 - Tr) + 62.45326 * Math.pow(1 - Tr, 2 / 3)
    - 135.1102 * (1 - Tr) + e * (4.79594 + 0.250047 * omega + 1.14188 * omega * omega));
  const C_tait = 0.0861488 + 0.0344483 * omega;
  const Ps = Psat_Pa(comp, T_K);
  if (!isFinite(Ps) || !isFinite(B) || B + Ps <= 0 || B + P_Pa <= 0) return Vs_m3pmol;
  return Vs_m3pmol * (1 - C_tait * Math.log((B + P_Pa) / (B + Ps)));
}

// ────────────────────────────────────────────────────────────────
// Thermal de Broglie Wavelength (for quantum corrections)
// ────────────────────────────────────────────────────────────────

/**
 * Thermal de Broglie wavelength: Λ = h / √(2πmkT)
 * Used in quantum corrections for He, H2, Ne at low temperatures.
 */
export function thermalDeBroglieWavelength_m(MW_kgpmol: number, T_K: number): number {
  const h = 6.62607015e-34; // Planck's constant (J·s)
  const kB = 1.380649e-23;  // Boltzmann constant (J/K)
  const NA = 6.02214076e23;
  const m = MW_kgpmol / NA; // mass per molecule (kg)
  return h / Math.sqrt(2 * Math.PI * m * kB * T_K);
}

// ────────────────────────────────────────────────────────────────
// Liquid Thermal Expansion Coefficient (improved)
// ────────────────────────────────────────────────────────────────

/**
 * Thermal expansion coefficient from Rackett density:
 *   α = -(1/V)(dV/dT) = (2/7)·ln(ZRA)·(1-Tr)^(-5/7) / (Tc·(1 + (1-Tr)^(2/7)))
 */
export function thermalExpansionRackett(comp: CanopyCompound, T_K: number): number {
  const ZRA = comp.rackett_ZRA ?? (0.29056 - 0.08775 * comp.omega);
  const Tr = T_K / comp.Tc_K;
  if (Tr >= 1 || ZRA <= 0) return NaN;
  const tau = 1 - Tr;
  return (2 / 7) * Math.log(ZRA) * Math.pow(tau, -5 / 7) / (comp.Tc_K * (1 + Math.pow(tau, 2 / 7)));
}

// ────────────────────────────────────────────────────────────────
// Viscosity Mixing Rules (Grunberg-Nissan with interaction)
// ────────────────────────────────────────────────────────────────

/**
 * Grunberg-Nissan with binary interaction parameters d_ij:
 *   ln(μ_mix) = Σ x_i·ln(μ_i) + ΣΣ x_i·x_j·d_ij
 *
 * When d_ij = 0, this reduces to the simple Grunberg-Nissan rule (no interaction).
 * pplcdefs MULIJ provides the d_ij parameters.
 */
export function mixtureLiquidViscosityGN_Pas(
  mu_i: number[], x: number[], d_ij?: number[][],
): number {
  const n = x.length;
  let lnMu = 0;
  for (let i = 0; i < n; i++) {
    if (x[i] > 0 && mu_i[i] > 0) lnMu += x[i] * Math.log(mu_i[i]);
  }
  if (d_ij) {
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        lnMu += x[i] * x[j] * d_ij[i][j];
      }
    }
  }
  return Math.exp(lnMu);
}

// ────────────────────────────────────────────────────────────────
// Additional property models from pplcdefs
// ────────────────────────────────────────────────────────────────

/**
 * DIPPR 105 saturated liquid density from compound's DNLDIP coefficients → [kmol/m³]
 */
export function liquidDensity_DIPPR105_kmolpm3(comp: CanopyCompound, T_K: number): number {
  const c = comp.dipprCoeffs['LiquidDensity'];
  if (!c) return NaN;
  return DIPPR105(T_K, c.A, c.B, c.C ?? 0, c.D ?? 0);
}

/**
 * Second virial coefficient B(T) from compound's SVRDIP coefficients → [m³/kmol]
 * DIPPR 104: B = A + B₁/T + C/T³ + D/T⁸ + E/T⁹
 */
export function secondVirial_m3pkmol(comp: CanopyCompound, T_K: number): number {
  const c = comp.dipprCoeffs['SecondVirialCoefficient'];
  if (!c) return NaN;
  return DIPPR104(T_K, c.A, c.B, c.C ?? 0, c.D ?? 0, c.E ?? 0);
}

/**
 * Truncated virial equation fugacity coefficient:
 *   ln(φ_i) = (2·Σ_j y_j·B_ij - B_mix)·P/(RT)
 * where B_mix = ΣΣ y_i·y_j·B_ij.
 * Returns φ_i array for each component.
 */
export function fugacityCoeffVirial(
  compounds: CanopyCompound[], T_K: number, P_Pa: number, y: number[]
): number[] {
  const n = compounds.length;
  const R = 8314; // J/(kmol·K)

  // Compute B_ii for each compound
  const B_ii: number[] = compounds.map(c => {
    const val = secondVirial_m3pkmol(c, T_K);
    return isFinite(val) ? val : 0;
  });

  // B_mix = ΣΣ y_i y_j B_ij  (arithmetic combining rule for cross-coefficients)
  let B_mix = 0;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const B_ij = (i === j) ? B_ii[i] : 0.5 * (B_ii[i] + B_ii[j]);
      B_mix += y[i] * y[j] * B_ij;
    }
  }

  // Fugacity coefficient: ln(φ_i) = (2·Σ_j y_j·B_ij - B_mix)·P/(RT)
  const phi = new Array(n);
  for (let i = 0; i < n; i++) {
    let sum_yB = 0;
    for (let j = 0; j < n; j++) {
      const B_ij = (i === j) ? B_ii[i] : 0.5 * (B_ii[i] + B_ii[j]);
      sum_yB += y[j] * B_ij;
    }
    phi[i] = Math.exp((2 * sum_yB - B_mix) * P_Pa / (R * T_K));
  }
  return phi;
}

/**
 * Li mixing rule for mixture liquid thermal conductivity → [W/(m·K)]
 *   1/λ_mix = Σ φ_i/λ_i  where φ_i = x_i·Vc_i / Σ x_j·Vc_j
 */
export function mixtureLiquidThermalConductivityLi_WpmK(
  compounds: CanopyCompound[], x: number[], T_K: number
): number {
  const n = x.length;
  const lambda_i = compounds.map(c => liquidThermalConductivity_WpmK(c, T_K));
  const Vc = compounds.map(c => c.Vc_m3pkmol ?? 0);
  let sumXV = 0;
  for (let i = 0; i < n; i++) sumXV += x[i] * Vc[i];
  if (sumXV <= 0) return NaN;

  let invLambda = 0;
  for (let i = 0; i < n; i++) {
    if (lambda_i[i] <= 0 || !isFinite(lambda_i[i])) continue;
    const phi_i = x[i] * Vc[i] / sumXV;
    invLambda += phi_i / lambda_i[i];
  }
  return invLambda > 0 ? 1 / invLambda : NaN;
}

/**
 * Wilke mixing rule for mixture vapor viscosity → [Pa·s]
 *   μ_mix = Σ (x_i·μ_i) / Σ_j x_j·Φ_ij
 * Φ_ij = [1 + (μ_i/μ_j)^0.5·(Mj/Mi)^0.25]² / [8·(1+Mi/Mj)]^0.5
 */
export function mixtureVaporViscosityWilke_Pas(
  compounds: CanopyCompound[], y: number[], T_K: number
): number {
  const n = y.length;
  const mu = compounds.map(c => vaporViscosity_Pas(c, T_K));
  const MW = compounds.map(c => c.molecularWeight);

  let mixMu = 0;
  for (let i = 0; i < n; i++) {
    if (!isFinite(mu[i]) || mu[i] <= 0 || y[i] <= 0) continue;
    let denom = 0;
    for (let j = 0; j < n; j++) {
      if (y[j] <= 0 || !isFinite(mu[j]) || mu[j] <= 0) continue;
      const phiIJ =
        Math.pow(1 + Math.sqrt(mu[i] / mu[j]) * Math.pow(MW[j] / MW[i], 0.25), 2)
        / Math.sqrt(8 * (1 + MW[i] / MW[j]));
      denom += y[j] * phiIJ;
    }
    mixMu += y[i] * mu[i] / denom;
  }
  return mixMu;
}

// ════════════════════════════════════════════════════════════════════════════
//  NRTL Activity Coefficient Model
// ════════════════════════════════════════════════════════════════════════════
//
// The NRTL (Non-Random Two-Liquid) model computes activity coefficients γ_i
// for liquid-phase non-ideality in VLE.
//
// Aspen convention (same as Renon & Prausnitz, 1968):
//   τ_ij = aij + bij/T + eij·ln(T) + fij·T
//   G_ij = exp(−α_ij · τ_ij)
//   α_ij = cij  (symmetric: α_ij = α_ji)
//
//   ln γ_i = (Σ_j x_j·τ_ji·G_ji) / (Σ_k x_k·G_ki)
//          + Σ_j [ x_j·G_ij / (Σ_k x_k·G_kj) ] · [ τ_ij − (Σ_m x_m·τ_mj·G_mj) / (Σ_k x_k·G_kj) ]
//

import type { NrtlBinaryParams } from './types';

/**
 * Look up NRTL binary parameters for compounds named `nameI` and `nameJ`.
 * Keys in the params record are "A|B" with A < B alphabetically.
 * Returns { aij, aji, bij, bji, cij, ... } oriented so "i" corresponds to `nameI`.
 */
export function lookupNrtlParams(
  nameI: string, nameJ: string,
  paramDb: Record<string, NrtlBinaryParams>
): NrtlBinaryParams | null {
  const [a, b] = nameI < nameJ ? [nameI, nameJ] : [nameJ, nameI];
  const key = `${a}|${b}`;
  const p = paramDb[key];
  if (!p) return null;
  // If nameI is the first in sorted order → use as-is
  if (nameI <= nameJ) return p;
  // Otherwise swap ij ↔ ji
  return {
    aij: p.aji, aji: p.aij,
    bij: p.bji, bji: p.bij,
    cij: p.cij, dij: p.dij,
    eij: p.eji ?? 0, eji: p.eij ?? 0,
    fij: p.fji ?? 0, fji: p.fij ?? 0,
    Tlower_K: p.Tlower_K, Tupper_K: p.Tupper_K,
    source: p.source,
  };
}

/**
 * Compute NRTL τ_ij at temperature T (K).
 *   τ_ij = aij + bij/T + eij·ln(T) + fij·T
 */
function nrtlTau(aij: number, bij: number, eij: number, fij: number, T_K: number): number {
  return aij + bij / T_K + eij * Math.log(T_K) + fij * T_K;
}

/**
 * NRTL activity coefficients for an N-component liquid mixture.
 *
 * @param names      Compound names (same order as `x`)
 * @param x          Liquid mole fractions
 * @param T_K        Temperature (K)
 * @param paramDb    Record of NrtlBinaryParams keyed by "A|B" (sorted alphabetically)
 * @returns          Array of activity coefficients γ_i
 */
export function gammaNRTL(
  names: string[], x: number[], T_K: number,
  paramDb: Record<string, NrtlBinaryParams>
): number[] {
  const n = names.length;

  // Build τ and G matrices (τ_ii = 0, G_ii = 1)
  const tau: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
  const G: number[][] = Array.from({ length: n }, () => new Array(n).fill(1));

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i === j) continue;
      const p = lookupNrtlParams(names[i], names[j], paramDb);
      if (!p) continue; // missing pair → treat as ideal (τ=0, G=1)
      const tauIJ = nrtlTau(p.aij, p.bij, p.eij ?? 0, p.fij ?? 0, T_K);
      const alphaIJ = p.cij;
      tau[i][j] = tauIJ;
      G[i][j] = Math.exp(-alphaIJ * tauIJ);
    }
  }

  // Compute ln(γ_i)
  const gamma = new Array(n);
  for (let i = 0; i < n; i++) {
    // Term 1: Σ_j x_j·τ_ji·G_ji / Σ_k x_k·G_ki
    let num1 = 0, den1 = 0;
    for (let j = 0; j < n; j++) {
      num1 += x[j] * tau[j][i] * G[j][i];
      den1 += x[j] * G[j][i];
    }
    const term1 = den1 > 0 ? num1 / den1 : 0;

    // Term 2: Σ_j [ x_j·G_ij / Σ_k x_k·G_kj ] · [ τ_ij − Σ_m x_m·τ_mj·G_mj / Σ_k x_k·G_kj ]
    let term2 = 0;
    for (let j = 0; j < n; j++) {
      let denJ = 0, numMJ = 0;
      for (let k = 0; k < n; k++) {
        denJ += x[k] * G[k][j];
        numMJ += x[k] * tau[k][j] * G[k][j];
      }
      if (denJ <= 0) continue;
      term2 += (x[j] * G[i][j] / denJ) * (tau[i][j] - numMJ / denJ);
    }

    gamma[i] = Math.exp(term1 + term2);
  }
  return gamma;
}

/**
 * UNIQUAC activity coefficient model for an N-component system.
 * Uses structural parameters r (volume) and q (surface area) from CanopyCompound.
 *
 * Note: This is a combinatorial + residual model. The binary interaction
 * parameters (τ_ij = exp(-Δu_ij / RT)) should be stored separately.
 * This function computes only the combinatorial contribution when
 * no binary params are provided (useful for athermal mixtures).
 *
 * Combinatorial:
 *   ln γ_i^C = ln(Φ_i/x_i) + (z/2)·q_i·ln(θ_i/Φ_i) + l_i − (Φ_i/x_i)·Σ_j x_j·l_j
 * where Φ_i = x_i·r_i / Σ x_j·r_j,  θ_i = x_i·q_i / Σ x_j·q_j
 *       l_i = (z/2)·(r_i - q_i) - (r_i - 1),  z = 10 (coordination number)
 */
export function gammaUNIQUAC_combinatorial(
  compounds: { uniquac_r?: number; uniquac_q?: number }[],
  x: number[]
): number[] {
  const n = x.length;
  const z = 10;  // coordination number

  const r = compounds.map(c => c.uniquac_r ?? 0);
  const q = compounds.map(c => c.uniquac_q ?? 0);
  const l = r.map((ri, i) => (z / 2) * (ri - q[i]) - (ri - 1));

  // Volume fraction Φ_i and surface fraction θ_i
  let sumXR = 0, sumXQ = 0;
  for (let i = 0; i < n; i++) { sumXR += x[i] * r[i]; sumXQ += x[i] * q[i]; }
  if (sumXR <= 0 || sumXQ <= 0) return new Array(n).fill(1);

  const Phi = x.map((xi, i) => xi * r[i] / sumXR);
  const theta = x.map((xi, i) => xi * q[i] / sumXQ);

  let sumXL = 0;
  for (let i = 0; i < n; i++) sumXL += x[i] * l[i];

  const gamma = new Array(n);
  for (let i = 0; i < n; i++) {
    if (x[i] <= 0) { gamma[i] = 1; continue; }
    const lnGammaC = Math.log(Phi[i] / x[i])
      + (z / 2) * q[i] * Math.log(theta[i] / Phi[i])
      + l[i] - (Phi[i] / x[i]) * sumXL;
    gamma[i] = Math.exp(lnGammaC);
  }
  return gamma;
}


// ═══════════════════════════════════════════════════════════════════════════════
// Hildebrand solubility parameter from Hvap and molar volume
// δ = √((ΔHvap − RT) / Vm)   [(J/m³)^0.5]
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Estimate Hildebrand solubility parameter from heat of vaporization and molar volume.
 * Uses DIPPR106 Hvap and liquid density at 25°C.
 */
export function hildebrandSolubilityParam(comp: CanopyCompound): number {
  const T = 298.15;
  const R = 8.314;
  // Get Hvap at 298 K
  const hvap = HvapBest_Jmol(comp, T);
  if (Number.isNaN(hvap) || hvap <= 0) return NaN;
  // Get liquid molar volume at 298 K (m³/mol)
  const rho_kmolpm3 = liquidDensityBest_kmolpm3(comp, T);
  if (Number.isNaN(rho_kmolpm3) || rho_kmolpm3 <= 0) return NaN;
  const Vm_m3pmol = 1 / (rho_kmolpm3 * 1000); // m³/mol
  // δ = √((ΔHvap − RT) / Vm)
  return Math.sqrt((hvap - R * T) / Vm_m3pmol);
}

// ═══════════════════════════════════════════════════════════════════════════════
// PC-SAFT EOS: Helmholtz free energy and pressure
// Gross & Sadowski (2001), Ind. Eng. Chem. Res. 40, 1244–1260
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Compute the PC-SAFT hard-sphere compressibility factor Z_hs for a pure component.
 * Uses Carnahan-Starling expression with segment packing fraction η.
 * η = (π/6)·ρ·m·σ³   where ρ = number density (1/ų)
 *
 * Z_hs = (1 + η + η² − η³) / (1 − η)³
 */
export function pcSaftZhs(eta: number): number {
  const e1 = 1 - eta;
  return (1 + eta + eta * eta - eta * eta * eta) / (e1 * e1 * e1);
}

/**
 * PC-SAFT reduced Helmholtz free energy for hard-chain reference (pure component).
 * ā_hc = m·ā_hs − (m−1)·ln(g_hs(σ))
 * where g_hs(σ) = (1−η/2)/(1−η)³
 */
export function pcSaftAhc(m: number, eta: number): number {
  const e1 = 1 - eta;
  // Hard-sphere reduced free energy per segment (Carnahan-Starling)
  const a_hs = (4 * eta - 3 * eta * eta) / (e1 * e1);
  // Radial distribution function at contact
  const g_hs = (1 - eta / 2) / (e1 * e1 * e1);
  return m * a_hs - (m - 1) * Math.log(g_hs);
}

// ═══════════════════════════════════════════════════════════════════════════════
// UNIFAC / UNIFAC-DMD bridge using stored APV140 data
// ═══════════════════════════════════════════════════════════════════════════════

// Standard UNIFAC subgroup table (Fredenslund 1975, groups 1–4)
// subgroupId: 1=CH3, 2=CH2, 3=CH, 4=C, 5=CH2=CH, 6=CH=CH, 7=CH2=C, 8=CH=C,
// 9=ACH, 10=AC, 11=ACCH3, 12=ACCH2, 13=ACCH
const UNIFAC_SUBGROUPS: UNIFACSubgroup[] = [
  { subgroupId: 1, mainGroupId: 1, R: 0.9011, Q: 0.848 },   // CH3
  { subgroupId: 2, mainGroupId: 1, R: 0.6744, Q: 0.540 },   // CH2
  { subgroupId: 3, mainGroupId: 1, R: 0.4469, Q: 0.228 },   // CH
  { subgroupId: 4, mainGroupId: 1, R: 0.2195, Q: 0.000 },   // C
  { subgroupId: 5, mainGroupId: 2, R: 1.3454, Q: 1.176 },   // CH2=CH
  { subgroupId: 6, mainGroupId: 2, R: 1.1167, Q: 0.867 },   // CH=CH
  { subgroupId: 7, mainGroupId: 2, R: 1.1173, Q: 0.988 },   // CH2=C
  { subgroupId: 8, mainGroupId: 2, R: 0.8886, Q: 0.676 },   // CH=C
  { subgroupId: 9, mainGroupId: 3, R: 0.5313, Q: 0.400 },   // ACH
  { subgroupId: 10, mainGroupId: 3, R: 0.3652, Q: 0.120 },  // AC
  { subgroupId: 11, mainGroupId: 4, R: 1.2663, Q: 0.968 },  // ACCH3
  { subgroupId: 12, mainGroupId: 4, R: 1.0396, Q: 0.660 },  // ACCH2
  { subgroupId: 13, mainGroupId: 4, R: 0.8121, Q: 0.348 },  // ACCH
];

/**
 * Compute UNIFAC activity coefficients using the stored APV140 interaction parameters.
 * Builds the UNIFACData structure from UNIFAC_SUBGROUPS and the given interaction
 * parameter map (key format "m|n" → a_mn).
 */
export function unifacGammaFromStore(
  x: number[],
  T_K: number,
  compounds: CanopyCompound[],
  interactionParams: Record<string, number>,
): number[] {
  // Build compGroups array from compounds
  const compGroups = compounds.map(c => c.unifacGroups ?? {});
  // Convert interaction params from "m|n" to "m-n" format
  const interactions: Record<string, number> = {};
  for (const [key, val] of Object.entries(interactionParams)) {
    interactions[key.replace('|', '-')] = val;
  }
  const data: UNIFACData = { subgroups: UNIFAC_SUBGROUPS, interactions };
  return unifacGamma(x, T_K, compGroups, data);
}

export function unifacDMDGammaFromStore(
  x: number[],
  T_K: number,
  compounds: CanopyCompound[],
  interactionParams: Record<string, { a: number; b: number; c: number }>,
): number[] {
  const compGroups = compounds.map(c => c.unifacGroups ?? {});
  const interactions: Record<string, { a: number; b: number; c: number }> = {};
  for (const [key, val] of Object.entries(interactionParams)) {
    interactions[key.replace('|', '-')] = val;
  }
  const data: UNIFACDMDData = { subgroups: UNIFAC_SUBGROUPS, interactions };
  return unifacDMDGamma(x, T_K, compGroups, data);
}

// ═══════════════════════════════════════════════════════════════════════════════
// DIPPR 123 — Sato-Riedel liquid thermal conductivity
// k = A·(1 + B·τ^(1/3) + C·τ^(2/3) + D·τ)  where τ = 1 − Tr
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Evaluate DIPPR equation 123 (Sato-Riedel) for liquid thermal conductivity.
 * @returns Thermal conductivity in W/(m·K)
 */
export function DIPPR123(T: number, Tc: number, A: number, B: number, C: number, D: number): number {
  const Tr = T / Tc;
  if (Tr >= 1) return NaN;
  const tau = 1 - Tr;
  return A * (1 + B * Math.cbrt(tau) + C * Math.pow(tau, 2 / 3) + D * tau);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Sublimation pressure from PSXANT (solid-phase extended Antoine)
// ln(P/Pa) = C1 + C2/(T+C3) + C4·T + C5·ln(T) + C6·T^C7
// Valid from Tmin to Tf (freezing point)
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Compute sublimation (solid→vapor) pressure using solid-phase extended Antoine.
 * @returns Pressure in Pa, or NaN if outside valid range or no data
 */
export function sublimationPressure_Pa(comp: CanopyCompound, T_K: number): number {
  const ps = comp.psxant;
  if (!ps) return NaN;
  if (T_K > ps.Tmax_K || T_K < ps.Tmin_K) return NaN;
  const lnP = ps.C1 + ps.C2 / (T_K + ps.C3) + ps.C4 * T_K
    + ps.C5 * Math.log(T_K) + ps.C6 * Math.pow(T_K, ps.C7);
  return Math.exp(lnP);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Water solubility — mole fraction of water in hydrocarbon phase
// ln(x_w) = A + B/T + C·ln(T)   [WATSOL from APV140]
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Mole fraction of water dissolved in hydrocarbon at temperature T.
 * Uses WATSOL correlation: ln(x_w) = A + B/T + C·ln(T)
 */
export function waterSolubilityInHC(comp: CanopyCompound, T_K: number): number {
  const w = comp.waterSolubility;
  if (!w) return NaN;
  if (T_K < w.Tmin_K || T_K > w.Tmax_K) return NaN;
  const lnXw = w.A + w.B / T_K + w.C * Math.log(T_K);
  return Math.exp(lnXw);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Hydrocarbon solubility in water — mole fraction of HC in aqueous phase
// ln(x_hc) = A + B/T + C·ln(T)   [HCSOL from APV140]
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Mole fraction of hydrocarbon dissolved in water at temperature T.
 * Uses HCSOL correlation: ln(x_hc) = A + B/T + C·ln(T)
 */
export function hcSolubilityInWater(comp: CanopyCompound, T_K: number): number {
  const h = comp.hcSolubility;
  if (!h) return NaN;
  if (T_K < h.Tmin_K || T_K > h.Tmax_K) return NaN;
  const lnXhc = h.A + h.B / T_K + h.C * Math.log(T_K);
  return Math.exp(lnXhc);
}

export interface WaterSaturationResult {
  zSaturated: number[];
  waterIndex: number;
  requiredWaterFraction: number;
  waterAddedFraction: number;
  excessFreeWaterFraction: number;
}

export function waterSaturatedComposition(
  compounds: CanopyCompound[],
  z: number[],
  T_K: number,
): WaterSaturationResult | null {
  const waterIndex = compounds.findIndex(compound => {
    const name = compound.name.trim().toUpperCase();
    return name === 'WATER' || name === 'H2O';
  });
  if (waterIndex < 0 || waterIndex >= z.length) return null;

  const nonWaterFraction = Math.max(0, 1 - (z[waterIndex] ?? 0));
  if (nonWaterFraction <= 1e-12) {
    return {
      zSaturated: [...z],
      waterIndex,
      requiredWaterFraction: z[waterIndex] ?? 0,
      waterAddedFraction: 0,
      excessFreeWaterFraction: 0,
    };
  }

  const normalizedOrganic = z.map((value, index) => index === waterIndex ? 0 : value / nonWaterFraction);
  let waterPerOrganicMol = 0;
  for (let i = 0; i < compounds.length; i++) {
    if (i === waterIndex) continue;
    const xi = normalizedOrganic[i] ?? 0;
    if (xi <= 0) continue;
    const xw = waterSolubilityInHC(compounds[i], T_K);
    if (!Number.isFinite(xw) || xw <= 0 || xw >= 0.5) continue;
    waterPerOrganicMol += xi * xw / Math.max(1 - xw, 1e-9);
  }

  const requiredWaterFraction = waterPerOrganicMol / (1 + waterPerOrganicMol);
  const currentWaterFraction = z[waterIndex] ?? 0;
  const targetWaterFraction = Math.max(currentWaterFraction, requiredWaterFraction);
  const zSaturated = z.map((value, index) => {
    if (index === waterIndex) return targetWaterFraction;
    return (value / Math.max(nonWaterFraction, 1e-12)) * (1 - targetWaterFraction);
  });

  return {
    zSaturated,
    waterIndex,
    requiredWaterFraction,
    waterAddedFraction: Math.max(0, requiredWaterFraction - currentWaterFraction),
    excessFreeWaterFraction: Math.max(0, currentWaterFraction - requiredWaterFraction),
  };
}

// ═══════════════════════════════════════════════════════════════════════════════
// Wilson activity coefficient from stored APV140 binary parameters
// ln(Λ_ij) = aij + bij/T + cij·ln(T) + dij·T
// ln(γ_i) = 1 − ln(Σ_j x_j·Λ_ij) − Σ_k x_k·Λ_ki / (Σ_j x_j·Λ_kj)
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Compute Wilson Λ_ij from Aspen-style parameters.
 * ln(Λ_ij) = aij + bij/T + cij·ln(T) + dij·T
 */
function wilsonLambda(aij: number, bij: number, cij: number, dij: number, T: number): number {
  return Math.exp(aij + bij / T + cij * Math.log(T) + dij * T);
}

/**
 * Wilson activity coefficients from stored binary parameters.
 * Handles 2-component systems only (for now).
 */
export function wilsonGammaFromStore(
  x: number[],
  T_K: number,
  compNames: string[],
  params: Record<string, { aij: number; aji: number; bij: number; bji: number; cij?: number; cji?: number; dij?: number; dji?: number }>,
): number[] {
  const n = x.length;
  // Build Λ matrix
  const Lambda: number[][] = Array.from({ length: n }, () => Array(n).fill(1));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i === j) continue;
      const key = [compNames[i], compNames[j]].sort().join('|');
      const p = params[key];
      if (!p) continue;
      // Determine if stored order matches i→j or j→i
      const iFirst = compNames[i] <= compNames[j];
      const a = iFirst ? p.aij : p.aji;
      const b = iFirst ? p.bij : p.bji;
      const c = iFirst ? (p.cij ?? 0) : (p.cji ?? 0);
      const d = iFirst ? (p.dij ?? 0) : (p.dji ?? 0);
      Lambda[i][j] = wilsonLambda(a, b, c, d, T_K);
    }
  }
  // Compute γ
  const gamma = new Array(n);
  for (let i = 0; i < n; i++) {
    let sum1 = 0;
    for (let j = 0; j < n; j++) sum1 += x[j] * Lambda[i][j];
    let sum2 = 0;
    for (let k = 0; k < n; k++) {
      let sumDenom = 0;
      for (let j = 0; j < n; j++) sumDenom += x[j] * Lambda[k][j];
      sum2 += x[k] * Lambda[k][i] / sumDenom;
    }
    gamma[i] = Math.exp(1 - Math.log(sum1) - sum2);
  }
  return gamma;
}

// ═══════════════════════════════════════════════════════════════════════════════
// COSTALD liquid density from stored APV140 data (VSTCTD + OMGCTD)
// Uses Hankinson-Thomson corresponding-states equation with dedicated V* and ω
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * COSTALD saturated liquid molar volume using APV140 characteristic volume (VSTCTD)
 * and COSTALD acentric factor (OMGCTD).
 * @returns Molar volume in m³/mol, or NaN if missing data / above Tc
 */
export function COSTALD_fromStore(comp: CanopyCompound, T_K: number): number {
  const Vstar = comp.VSTCTD;
  const omega = comp.omegaCostald;
  if (Vstar == null || omega == null) return NaN;
  // VSTCTD is in m³/kmol — convert to m³/mol
  return COSTALD_Vs_m3pmol(comp.Tc_K, Vstar / 1000, omega, T_K);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Rackett liquid density from ZRA (Spencer-Danner modification)
// V_sat = (R·Tc / Pc) · ZRA^[1 + (1-Tr)^(2/7)]
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Rackett equation for saturated liquid molar volume.
 * Uses the modified Rackett parameter ZRA from APV140 (RKTZRA).
 * @returns Molar volume in m³/mol, or NaN if missing data / above Tc
 */
export function rackett_Vm_m3pmol(comp: CanopyCompound, T_K: number): number {
  const ZRA = comp.rackett_ZRA;
  if (ZRA == null || T_K >= comp.Tc_K) return NaN;
  const Tr = T_K / comp.Tc_K;
  const exp = 1 + Math.pow(1 - Tr, 2 / 7);
  return (R * comp.Tc_K / comp.Pc_Pa) * Math.pow(ZRA, exp);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Gibbs energy of ideal mixing
// ΔG_mix = R·T · Σ x_i · ln(x_i)
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Ideal Gibbs energy of mixing (J/mol).
 * ΔG_mix^ideal = R·T · Σ x_i · ln(x_i)
 * Always negative for non-trivial mixtures.
 */
export function idealGibbsMixing_Jpmol(T_K: number, x: number[]): number {
  let sum = 0;
  for (let i = 0; i < x.length; i++) {
    if (x[i] > 0) sum += x[i] * Math.log(x[i]);
  }
  return R * T_K * sum;
}

/**
 * Ideal entropy of mixing (J/mol·K).
 * ΔS_mix^ideal = -R · Σ x_i · ln(x_i)
 */
export function idealEntropyMixing_JpmolK(x: number[]): number {
  let sum = 0;
  for (let i = 0; i < x.length; i++) {
    if (x[i] > 0) sum += x[i] * Math.log(x[i]);
  }
  return -R * sum;
}

/**
 * Excess Gibbs energy from activity coefficients (J/mol).
 * G^E = R·T · Σ x_i · ln(γ_i)
 */
export function excessGibbs_Jpmol(T_K: number, x: number[], gamma: number[]): number {
  let sum = 0;
  for (const i of x.keys()) {
    if (x[i] > 0) sum += x[i] * Math.log(gamma[i]);
  }
  return R * T_K * sum;
}

// ═══════════════════════════════════════════════════════════════════════════════
// UNIQUAC γ from stored APV140 binary parameters
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * UNIQUAC activity coefficients from stored APV140 binary parameters.
 * Looks up r, q from compounds and τ_ij params from the UNIQUAC_BINARY_PARAMS table.
 * APV140 UNIQUAC params use extended form: τ_ij = exp(aij + bij/T).
 */
export function uniquacGammaFromStore(
  x: number[],
  T_K: number,
  compounds: CanopyCompound[],
  params: Record<string, { aij: number; aji: number; bij: number; bji: number; cij?: number; cji?: number; dij?: number; dji?: number }>,
): number[] {
  const n = x.length;
  const r = compounds.map(c => c.uniquac_r ?? 1);
  const q = compounds.map(c => c.uniquac_q ?? 1);
  // Build a_ext and b_ext matrices from stored params
  const a_ext: number[][] = Array.from({ length: n }, () => Array(n).fill(0));
  const b_ext: number[][] = Array.from({ length: n }, () => Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i === j) continue;
      const ni = compounds[i].name;
      const nj = compounds[j].name;
      // Lookup uses CSV names: ISOPROPYLBENZENE for CUMENE
      const lookupName = (name: string) =>
        name === 'CUMENE' ? 'ISOPROPYLBENZENE' : name;
      const key = [lookupName(ni), lookupName(nj)].sort().join('|');
      const p = params[key];
      if (!p) continue;
      const iFirst = lookupName(ni) <= lookupName(nj);
      a_ext[i][j] = iFirst ? p.aij : p.aji;
      b_ext[i][j] = iFirst ? p.bij : p.bji;
    }
  }
  const zeroMatrix = Array.from({ length: compounds.length }, () => new Array(compounds.length).fill(0));
  return uniquacGamma(x, T_K, { r, q, a: zeroMatrix, b: zeroMatrix, a_ext, b_ext });
}

// ═══════════════════════════════════════════════════════════════════════════════
// Henry's constant from stored APV140 binary parameters
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Henry's constant from stored APV140 HENRY-AP parameters.
 * ln(H) = A + B/T + C·ln(T) + D·T + E·T²
 * @returns Henry's constant in Pa, or NaN if outside valid range or not found.
 */
export function henryConstantFromStore(
  soluteName: string,
  solventName: string,
  T_K: number,
  params: Record<string, { A: number; B: number; C: number; D: number; E: number; Tmin_K: number; Tmax_K: number }>,
): number {
  const key = `${soluteName}|${solventName}`;
  const p = params[key];
  if (!p) return NaN;
  if (T_K < p.Tmin_K || T_K > p.Tmax_K) return NaN;
  const lnH = p.A + p.B / T_K + p.C * Math.log(T_K) + p.D * T_K + p.E * T_K * T_K;
  return Math.exp(Math.max(-100, Math.min(100, lnH)));
}

// ═══════════════════════════════════════════════════════════════════════════════
// NRTL γ from stored binary parameters
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * NRTL activity coefficients from stored binary parameters (APV140 or NIST).
 * Builds a_ext, b_ext, alpha matrices from the stored NrtlBinaryParams table,
 * then delegates to nrtlGamma().
 */
export function nrtlGammaFromStore(
  x: number[],
  T_K: number,
  compNames: string[],
  params: Record<string, { aij: number; aji: number; bij: number; bji: number; cij: number }>,
): number[] {
  const n = x.length;
  const a_ext: number[][] = Array.from({ length: n }, () => Array(n).fill(0));
  const b_ext: number[][] = Array.from({ length: n }, () => Array(n).fill(0));
  const alpha: number[][] = Array.from({ length: n }, () => Array(n).fill(0.3));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i === j) continue;
      const key = [compNames[i], compNames[j]].sort().join('|');
      const p = params[key];
      if (!p) continue;
      const iFirst = compNames[i] <= compNames[j];
      a_ext[i][j] = iFirst ? p.aij : p.aji;
      b_ext[i][j] = iFirst ? p.bij : p.bji;
      alpha[i][j] = p.cij;
    }
  }
  return nrtlGamma(x, T_K, { A: [[]], B: [[]], alpha, a_ext, b_ext });
}

// ═══════════════════════════════════════════════════════════════════════════════
// Solid density from stored DIPPR 100 coefficients
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Solid density from DIPPR 100: ρ = A + B·T + C·T² + D·T³ + E·T⁴  (kmol/m³)
 * Uses SolidDensity from dipprCoeffs in store.
 */
export function solidDensity_DIPPR100_kmolpm3(comp: CanopyCompound, T_K: number): number {
  const c = comp.dipprCoeffs['SolidDensity'];
  if (!c) return NaN;
  return DIPPR100(T_K, c.A, c.B, c.C ?? 0, c.D ?? 0, c.E ?? 0);
}

// ═══════════════════════════════════════════════════════════════════════════════
// PR/SRK kij lookup from stored EOS binary params
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Compute kij(T) = kij1 + kij2·T + kij3·T² from stored EOS binary params.
 * Returns 0 if the pair is not found.
 */
export function eosKijFromStore(
  comp1Name: string,
  comp2Name: string,
  T_K: number,
  params: Record<string, { kij1: number; kij2?: number; kij3?: number }>,
): number {
  const key = [comp1Name, comp2Name].sort().join('|');
  const p = params[key];
  if (!p) return 0;
  return p.kij1 + (p.kij2 ?? 0) * T_K + (p.kij3 ?? 0) * T_K * T_K;
}

// ═══════════════════════════════════════════════════════════════════════════════
// Water solubility in hydrocarbon phase
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Mole fraction of water dissolved in a hydrocarbon phase.
 * ln(x_w) = A + B/T + C·ln(T)
 * Returns NaN if out of valid temperature range or compound has no data.
 */
export function waterSolubility_moleFrac(comp: CanopyCompound, T_K: number): number {
  const p = comp.waterSolubility;
  if (!p) return NaN;
  if (T_K < p.Tmin_K || T_K > p.Tmax_K) return NaN;
  const lnx = p.A + p.B / T_K + p.C * Math.log(T_K);
  return Math.exp(lnx);
}

/**
 * Mole fraction of hydrocarbon dissolved in the water phase.
 * ln(x_hc) = A + B/T + C·ln(T)
 * Returns NaN if out of valid temperature range or compound has no data.
 */
export function hcSolubility_moleFrac(comp: CanopyCompound, T_K: number): number {
  const p = comp.hcSolubility;
  if (!p) return NaN;
  if (T_K < p.Tmin_K || T_K > p.Tmax_K) return NaN;
  const lnx = p.A + p.B / T_K + p.C * Math.log(T_K);
  return Math.exp(lnx);
}

// ═══════════════════════════════════════════════════════════════════════════════
// VTPR volume-translated Peng-Robinson
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * VTPR (volume-translated PR) molar volume correction.
 * V_VTPR = V_PR + c  where c = volumeTranslation_PR (m³/kmol).
 * The VTPR model uses its own Tc, Pc, and Twu alpha parameters (VTPRTC, VTPRPC, VTPRTW)
 * that may differ from standard values.
 */
export function vtprVolumeCorrection_m3pkmol(comp: CanopyCompound): number {
  return comp.volumeTranslation_PR ?? 0;
}

// ─── RKS T-dependent kij evaluator ──────────────────────────────────────────
/**
 * Evaluate RKS binary interaction parameter: kij = kij1 + kij2·T + kij3/T
 * Note the different form from PR (which uses kij3·T²).
 */
export function rksKijFromStore(
  comp1Name: string, comp2Name: string, T_K: number,
  params: Record<string, { kij1: number; kij2: number; kij3: number }>
): number {
  const lookupName = (n: string) => n === 'CUMENE' ? 'ISOPROPYLBENZENE' : n;
  const n1 = lookupName(comp1Name), n2 = lookupName(comp2Name);
  const key = [n1, n2].sort().join('|');
  const p = params[key] ?? params[[n2, n1].sort().join('|')];
  if (!p) return 0;
  return p.kij1 + p.kij2 * T_K + p.kij3 / T_K;
}

// ─── Extended Antoine vapor pressure (PLXANT) ───────────────────────────────
/**
 * Extended Antoine equation (PLXANT):
 *   ln(P/Pa) = C1 + C2/(T+C3) + C4·T + C5·ln(T) + C6·T^C7
 * Returns vapor pressure in Pa.
 */
export function extendedAntoinePsat_Pa(comp: CanopyCompound, T_K: number): number {
  const a = comp.antoine;
  if (!a) return NaN;
  const lnP = a.C1 + a.C2 / (T_K + a.C3) + a.C4 * T_K + a.C5 * Math.log(T_K) + a.C6 * Math.pow(T_K, a.C7);
  return Math.exp(lnP);
}

// ─── Dielectric constant evaluator ──────────────────────────────────────────
/**
 * Dielectric constant: ε(T) = e1 + e2·(1/T - 1/e3)
 * If e2=0, returns e1 (constant).
 */
export function dielectricConstant(comp: CanopyCompound, T_K: number): number {
  const d = comp.dielectricConst;
  if (!d) return NaN;
  return d.e1 + d.e2 * (1 / T_K - 1 / d.e3_K);
}
