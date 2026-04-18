// ─── Canopy Unit Operation Solver ──────────────────────────────────────────────
// Solves individual unit operations given their inlet streams and parameters.
// ──────────────────────────────────────────────────────────────────────────────

import type { CanopyCompound, MaterialStream, UnitOperation, StoichiometricReaction, ColumnSpec, Phase } from './types';
import {
  flashPT, flashPQ,
  streamEnthalpy, enthalpyIG, enthalpyLiquid, entropyIG, gibbsIG, streamEntropy,
  Psat_Pa, liquidMolarVolume_m3pmol, mixtureLiquidMolarVolume,
  liquidDensity_kmolpm3, computeGamma,
  CpIG_JmolK, bubblePointT_Ideal, dewPointT_Ideal, dewPointT,
  KvaluesIdeal, flashPT_VLLE,
  computeKvalues, bubblePointT, KvaluesWilsonEstimate,
  waterSolubilityInHC, hcSolubilityInWater,
  bubblePointTWaterSaturated, dewPointTWaterSaturated, flashPTWaterSaturated,
} from './thermo';
import type { FluidPackageType, NRTLParams, WilsonParams, UNIQUACParams, UNIFACData, InteractionParams } from './thermo';

const R = 8.314462618; // J/(mol·K)

/** Options for non-ideal thermodynamics, threaded through solvers. */
export interface ThermoOptions {
  fluidPackage: FluidPackageType;
  interactionParams?: InteractionParams;
}

const DEFAULT_THERMO: ThermoOptions = { fluidPackage: 'Ideal' };

// ────────────────────────────────────────────────────────────────
// Mixer: combines N inlet streams into 1 outlet stream
// ────────────────────────────────────────────────────────────────

export function solveMixer(
  compounds: CanopyCompound[],
  inlets: MaterialStream[],
  outletPressure_Pa: number,
  thermo: ThermoOptions = DEFAULT_THERMO,
): { outlet: Partial<MaterialStream>; duty_W: number } {
  const n = compounds.length;
  let totalMol = 0;
  const componentMols = new Array(n).fill(0);

  for (const inlet of inlets) {
    for (let i = 0; i < n; i++) {
      componentMols[i] += inlet.totalFlow_molps * inlet.moleFractions[i];
    }
    totalMol += inlet.totalFlow_molps;
  }

  const z = componentMols.map(m => totalMol > 0 ? m / totalMol : 0);

  // Enthalpy balance: H_out must equal sum of H_in (adiabatic mixer)
  let H_in_total_W = 0; // J/s = W
  for (const inlet of inlets) {
    H_in_total_W += inlet.totalFlow_molps * inlet.H_Jpmol;
  }
  const H_mix_Jpmol = totalMol > 0 ? H_in_total_W / totalMol : 0;

  // PQ flash at outlet P with Q=0 to find correct T (proper energy balance)
  const T_avg = totalMol > 0
    ? inlets.reduce((s, inp) => s + inp.T_K * inp.totalFlow_molps, 0) / totalMol
    : 298.15;
  const flash = flashPQ(compounds, z, outletPressure_Pa, H_mix_Jpmol, 0, totalMol,
    thermo.fluidPackage, thermo.interactionParams, T_avg);

  return {
    outlet: {
      T_K: flash.T_K,
      P_Pa: outletPressure_Pa,
      totalFlow_molps: totalMol,
      moleFractions: z,
      phase: flash.phase,
      vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    duty_W: 0, // Adiabatic mixer
  };
}

// ────────────────────────────────────────────────────────────────
// Heater: heats/cools a stream to a target temperature
// ────────────────────────────────────────────────────────────────

export function solveHeater(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  params: { targetT_K: number; outletP_Pa: number; waterSaturatedFeed?: boolean },
  thermo: ThermoOptions = DEFAULT_THERMO,
): { outlet: Partial<MaterialStream>; duty_W: number; waterAddedForSaturation_molps: number; excessFreeWater_molps: number } {
  const z = [...inlet.moleFractions];
  const flash = params.waterSaturatedFeed
    ? flashPTWaterSaturated(compounds, z, params.targetT_K, params.outletP_Pa, thermo.fluidPackage, thermo.interactionParams)
    : { ...flashPT(compounds, z, params.targetT_K, params.outletP_Pa, thermo.fluidPackage, thermo.interactionParams), waterAddedFraction: 0, excessFreeWaterFraction: 0 };

  const H_in = inlet.H_Jpmol;
  const H_out = flash.H_Jpmol;
  const saturatedFeedFlow = inlet.totalFlow_molps * (1 + flash.waterAddedFraction);
  const duty_W = saturatedFeedFlow * H_out - inlet.totalFlow_molps * H_in; // W

  return {
    outlet: {
      T_K: params.targetT_K,
      P_Pa: params.outletP_Pa,
      totalFlow_molps: saturatedFeedFlow,
      moleFractions: z,
      phase: flash.phase,
      vaporFraction: flash.vaporFraction,
      H_Jpmol: H_out,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    duty_W,
    waterAddedForSaturation_molps: inlet.totalFlow_molps * flash.waterAddedFraction,
    excessFreeWater_molps: inlet.totalFlow_molps * flash.excessFreeWaterFraction,
  };
}

// ────────────────────────────────────────────────────────────────
// Flash Drum: separates a two-phase stream into V and L
// ────────────────────────────────────────────────────────────────

export interface FlashDrumResult {
  vapor: Partial<MaterialStream>;
  liquid: Partial<MaterialStream>;
  duty_W: number;
  waterAddedForSaturation_molps: number;
  excessFreeWater_molps: number;
}

export function solveFlashDrum(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  params: { T_K: number; P_Pa: number; waterSaturatedFeed?: boolean },
  thermo: ThermoOptions = DEFAULT_THERMO,
): FlashDrumResult {
  const z = [...inlet.moleFractions];
  const flash = params.waterSaturatedFeed
    ? flashPTWaterSaturated(compounds, z, params.T_K, params.P_Pa, thermo.fluidPackage, thermo.interactionParams)
    : { ...flashPT(compounds, z, params.T_K, params.P_Pa, thermo.fluidPackage, thermo.interactionParams), waterAddedFraction: 0, excessFreeWaterFraction: 0 };

  const F = inlet.totalFlow_molps;
  const saturatedFeedFlow = F * (1 + flash.waterAddedFraction);
  const V_flow = saturatedFeedFlow * flash.vaporFraction;
  const L_flow = saturatedFeedFlow * (1 - flash.vaporFraction);

  // Enthalpy of outlet streams
  const H_V = streamEnthalpy(compounds, params.T_K, 1, flash.x, flash.y);
  const H_L = streamEnthalpy(compounds, params.T_K, 0, flash.x, flash.y);

  // Energy balance: Q = V·H_V + L·H_L - F·H_in
  const duty_W = V_flow * H_V + L_flow * H_L - F * inlet.H_Jpmol;

  return {
    vapor: {
      T_K: params.T_K,
      P_Pa: params.P_Pa,
      totalFlow_molps: V_flow,
      moleFractions: flash.y,
      phase: 'Vapor' as const,
      vaporFraction: 1,
      H_Jpmol: H_V,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    liquid: {
      T_K: params.T_K,
      P_Pa: params.P_Pa,
      totalFlow_molps: L_flow,
      moleFractions: flash.x,
      phase: 'Liquid' as const,
      vaporFraction: 0,
      H_Jpmol: H_L,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    duty_W,
    waterAddedForSaturation_molps: F * flash.waterAddedFraction,
    excessFreeWater_molps: F * flash.excessFreeWaterFraction,
  };
}

// ────────────────────────────────────────────────────────────────
// Valve: isenthalpic pressure drop
// ────────────────────────────────────────────────────────────────

export function solveValve(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  params: { outletP_Pa: number; waterSaturatedFeed?: boolean },
  thermo: ThermoOptions = DEFAULT_THERMO,
): { outlet: Partial<MaterialStream>; duty_W: number; waterAddedForSaturation_molps: number; excessFreeWater_molps: number } {
  // Isenthalpic: find T at new P that gives same enthalpy
  const z = [...inlet.moleFractions];
  const H_target = inlet.H_Jpmol;
  let waterAddedForSaturation_molps = 0;
  let excessFreeWater_molps = 0;

  // Newton-Raphson on T to match enthalpy
  let T = inlet.T_K;
  for (let iter = 0; iter < 50; iter++) {
    const flash = params.waterSaturatedFeed
      ? flashPTWaterSaturated(compounds, z, T, params.outletP_Pa, thermo.fluidPackage, thermo.interactionParams)
      : { ...flashPT(compounds, z, T, params.outletP_Pa, thermo.fluidPackage, thermo.interactionParams), waterAddedFraction: 0, excessFreeWaterFraction: 0 };
    const H = flash.H_Jpmol;
    const saturatedFeedFlow = inlet.totalFlow_molps * (1 + flash.waterAddedFraction);
    const err = saturatedFeedFlow * H - inlet.totalFlow_molps * H_target;
    if (Math.abs(err) < 0.01) {
      waterAddedForSaturation_molps = inlet.totalFlow_molps * flash.waterAddedFraction;
      excessFreeWater_molps = inlet.totalFlow_molps * flash.excessFreeWaterFraction;
      return {
        outlet: {
          T_K: T,
          P_Pa: params.outletP_Pa,
          totalFlow_molps: saturatedFeedFlow,
          moleFractions: z,
          phase: flash.phase,
          vaporFraction: flash.vaporFraction,
          H_Jpmol: H,
          x_liquid: flash.x,
          y_vapor: flash.y,
          solved: true,
        },
        duty_W: 0,
        waterAddedForSaturation_molps,
        excessFreeWater_molps,
      };
    }
    // Numerical derivative
    const dT = 0.1;
    const flash2 = params.waterSaturatedFeed
      ? flashPTWaterSaturated(compounds, z, T + dT, params.outletP_Pa, thermo.fluidPackage, thermo.interactionParams)
      : { ...flashPT(compounds, z, T + dT, params.outletP_Pa, thermo.fluidPackage, thermo.interactionParams), waterAddedFraction: 0, excessFreeWaterFraction: 0 };
    const saturatedFeedFlow2 = inlet.totalFlow_molps * (1 + flash2.waterAddedFraction);
    const dHdT = (saturatedFeedFlow2 * flash2.H_Jpmol - saturatedFeedFlow * H) / dT;
    if (Math.abs(dHdT) < 1e-15) break;
    T -= err / dHdT;
    if (T < 100) T = 100;
    if (T > 1000) T = 1000;
  }

  // Fallback: just do flash at inlet T
  const flash = params.waterSaturatedFeed
    ? flashPTWaterSaturated(compounds, z, T, params.outletP_Pa, thermo.fluidPackage, thermo.interactionParams)
    : { ...flashPT(compounds, z, T, params.outletP_Pa, thermo.fluidPackage, thermo.interactionParams), waterAddedFraction: 0, excessFreeWaterFraction: 0 };
  waterAddedForSaturation_molps = inlet.totalFlow_molps * flash.waterAddedFraction;
  excessFreeWater_molps = inlet.totalFlow_molps * flash.excessFreeWaterFraction;
  return {
    outlet: {
      T_K: T,
      P_Pa: params.outletP_Pa,
      totalFlow_molps: inlet.totalFlow_molps * (1 + flash.waterAddedFraction),
      moleFractions: z,
      phase: flash.phase,
      vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    duty_W: 0,
    waterAddedForSaturation_molps,
    excessFreeWater_molps,
  };
}

// ────────────────────────────────────────────────────────────────
// Splitter: splits a stream into multiple outlets at same T,P
// ────────────────────────────────────────────────────────────────

export function solveSplitter(
  inlet: MaterialStream,
  splitRatios: number[], // e.g. [0.5, 0.5] for two equal outlets
): Partial<MaterialStream>[] {
  return splitRatios.map(ratio => ({
    T_K: inlet.T_K,
    P_Pa: inlet.P_Pa,
    totalFlow_molps: inlet.totalFlow_molps * ratio,
    moleFractions: [...inlet.moleFractions],
    phase: inlet.phase,
    vaporFraction: inlet.vaporFraction,
    H_Jpmol: inlet.H_Jpmol,
    x_liquid: [...inlet.x_liquid],
    y_vapor: [...inlet.y_vapor],
    solved: true,
  }));
}

// ────────────────────────────────────────────────────────────────
// Pump: liquid pressure increase with efficiency
// ────────────────────────────────────────────────────────────────
//
// Equations from Aspen Plus reference & aspen_unit_ops_math_reference.md:
//   W_ideal = F · V_L · (P2 - P1)          [W]
//   W_actual = W_ideal / η                  [W]
//   H_out = H_in + W_actual / F             [J/mol]
//   Head = ΔP / (ρ_L · g)                  [m]
//   NPSH_A = (P_suction - Pv) / (ρ_L · g)  [m]
//
// Spec options: outlet pressure, pressure increase, or pressure ratio.

export interface PumpResult {
  outlet: Partial<MaterialStream>;
  work_W: number;        // Actual shaft work (W) — positive = work input
  idealWork_W: number;   // Isentropic work (W)
  head_m: number;        // Developed head (m)
  NPSH_avail_m: number;  // NPSH available (m)
  efficiency: number;    // Used efficiency
}

export function solvePump(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  params: {
    outletP_Pa: number;
    efficiency: number; // 0 < η ≤ 1
    driverEfficiency?: number; // motor efficiency (default 1)
  },
  thermo: ThermoOptions = DEFAULT_THERMO,
): PumpResult {
  const { outletP_Pa, efficiency } = params;
  const eta = Math.max(0.01, Math.min(1, efficiency));
  const z = [...inlet.moleFractions];

  // Liquid molar volume at inlet conditions (m³/mol)
  const V_L = mixtureLiquidMolarVolume(compounds, z, inlet.T_K);
  const V_L_safe = isFinite(V_L) && V_L > 0 ? V_L : 1e-4; // fallback ~100 cm³/mol

  const deltaP = outletP_Pa - inlet.P_Pa;

  // Ideal (isentropic) work: W = F · V_L · ΔP
  const idealWork_W = inlet.totalFlow_molps * V_L_safe * deltaP;

  // Actual work considering efficiency
  const actualWork_W = idealWork_W / eta;

  // Outlet enthalpy = inlet + work/F
  const H_out = inlet.H_Jpmol + (inlet.totalFlow_molps > 0 ? actualWork_W / inlet.totalFlow_molps : 0);

  // PQ flash to find outlet T (work adds heat to liquid)
  const flash = flashPQ(
    compounds, z, outletP_Pa, inlet.H_Jpmol,
    actualWork_W, inlet.totalFlow_molps,
    thermo.fluidPackage, thermo.interactionParams, inlet.T_K
  );

  // Pump head: H = ΔP / (ρ·g)
  const rhoMolar = V_L_safe > 0 ? 1 / V_L_safe : NaN; // mol/m³
  const MW_avg = compounds.reduce((s, c, i) => s + z[i] * c.molecularWeight, 0); // g/mol
  const rhoMass = rhoMolar * MW_avg / 1000; // kg/m³
  const g = 9.80665;
  const head_m = rhoMass > 0 ? deltaP / (rhoMass * g) : 0;

  // NPSH available (simplified — no elevation/velocity terms)
  const Pv_avg = compounds.reduce((s, c, i) => {
    const pv = Psat_Pa(c, inlet.T_K);
    return s + z[i] * (isFinite(pv) ? pv : 0);
  }, 0);
  const NPSH_avail_m = rhoMass > 0
    ? (inlet.P_Pa - Pv_avg) / (rhoMass * g)
    : 0;

  return {
    outlet: {
      T_K: flash.T_K,
      P_Pa: outletP_Pa,
      totalFlow_molps: inlet.totalFlow_molps,
      moleFractions: z,
      phase: flash.phase,
      vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    work_W: actualWork_W,
    idealWork_W,
    head_m,
    NPSH_avail_m,
    efficiency: eta,
  };
}

// ────────────────────────────────────────────────────────────────
// Compressor: isentropic or polytropic gas compression
// ────────────────────────────────────────────────────────────────
//
// From aspen_unit_ops_math_reference.md:
//   Isentropic:
//     1. Flash inlet → get H1, then ideal gas Cp/Cv ratio κ
//     2. T2s = T1 · (P2/P1)^((κ-1)/κ)  (ideal gas isentropic)
//     3. H2s from flash at (T2s, P2)
//     4. Ws = F · (H2s - H1) / η_s
//     5. H2 = H1 + Ws/F → flash at (P2, H2)
//   Polytropic:
//     n/(n-1) = κ/(κ-1) · 1/η_p    (ASME method)
//     T2 = T1 · (P2/P1)^((n-1)/n)
//     Wp = F · Cp · T1 · [(P2/P1)^((n-1)/n) - 1] / η_p

export interface CompressorResult {
  outlet: Partial<MaterialStream>;
  work_W: number;          // Actual shaft work (W)
  idealWork_W: number;     // Isentropic work (W)
  polytropicHead_J_kg: number;
  dischargeT_K: number;
  efficiency: number;
}

export function solveCompressor(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  params: {
    outletP_Pa: number;
    efficiency: number;       // isentropic efficiency (0..1)
    model?: 'isentropic' | 'polytropic'; // default: isentropic
  },
  thermo: ThermoOptions = DEFAULT_THERMO,
): CompressorResult {
  const { outletP_Pa, efficiency } = params;
  const model = params.model ?? 'isentropic';
  const eta = Math.max(0.01, Math.min(1, efficiency));
  const z = [...inlet.moleFractions];

  // Mean Cp and Cv (ideal gas) — κ = Cp/Cv
  const Cp_avg = compounds.reduce((s, c, i) => s + z[i] * CpIG_JmolK(c, inlet.T_K), 0);
  const Cv_avg = Cp_avg - R; // ideal gas: Cp - Cv = R
  const kappa = Cv_avg > 0 ? Cp_avg / Cv_avg : 1.4;

  const P_ratio = outletP_Pa / inlet.P_Pa;
  const MW_avg = compounds.reduce((s, c, i) => s + z[i] * c.molecularWeight, 0); // g/mol

  let T2: number;
  let W_ideal: number;
  let W_actual: number;

  if (model === 'polytropic') {
    // Polytropic: n/(n-1) = κ/(κ-1) · 1/η_p
    const n = 1 / (1 - (kappa - 1) / (kappa * eta));
    T2 = inlet.T_K * Math.pow(P_ratio, (n - 1) / n);
    // Polytropic work
    W_ideal = inlet.totalFlow_molps * Cp_avg * inlet.T_K * (Math.pow(P_ratio, (kappa - 1) / kappa) - 1);
    W_actual = inlet.totalFlow_molps * Cp_avg * inlet.T_K * (Math.pow(P_ratio, (n - 1) / n) - 1);
  } else {
    // Isentropic
    const T2s = inlet.T_K * Math.pow(P_ratio, (kappa - 1) / kappa);
    // Isentropic enthalpy change
    const H1 = inlet.H_Jpmol;
    const flash2s = flashPT(compounds, z, T2s, outletP_Pa,
      thermo.fluidPackage, thermo.interactionParams);
    const H2s = flash2s.H_Jpmol;
    W_ideal = inlet.totalFlow_molps * (H2s - H1); // W (ideal)
    W_actual = W_ideal / eta;
    // Actual outlet enthalpy
    const H2_actual = H1 + (inlet.totalFlow_molps > 0 ? W_actual / inlet.totalFlow_molps : 0);
    // Find actual T2 via PQ flash
    const flashActual = flashPQ(
      compounds, z, outletP_Pa, H1, W_actual, inlet.totalFlow_molps,
      thermo.fluidPackage, thermo.interactionParams, T2s
    );
    T2 = flashActual.T_K;
  }

  const flash = flashPT(compounds, z, T2, outletP_Pa,
    thermo.fluidPackage, thermo.interactionParams);

  // Polytropic head (J/kg)
  const polytropicHead = MW_avg > 0 ? W_actual / (inlet.totalFlow_molps * MW_avg / 1000) : 0;

  return {
    outlet: {
      T_K: T2,
      P_Pa: outletP_Pa,
      totalFlow_molps: inlet.totalFlow_molps,
      moleFractions: z,
      phase: flash.phase,
      vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    work_W: W_actual,
    idealWork_W: W_ideal,
    polytropicHead_J_kg: polytropicHead,
    dischargeT_K: T2,
    efficiency: eta,
  };
}

// ────────────────────────────────────────────────────────────────
// HeatX: Two-stream heat exchanger (LMTD method)
// ────────────────────────────────────────────────────────────────
//
// Shortcut mode uses LMTD with correction factor F_T:
//   Q = U · A · F_T · ΔTLM
//
// Spec options:
//   - Hot outlet T specified → calculate Q, then cold outlet
//   - Cold outlet T specified → calculate Q, then hot outlet
//   - Duty Q specified → calculate both outlet temps
//   - UA specified → iterative convergence

export interface HeatXResult {
  hotOutlet: Partial<MaterialStream>;
  coldOutlet: Partial<MaterialStream>;
  duty_W: number;            // Heat transferred (W), positive = hot→cold
  LMTD_K: number;            // Log-mean temperature difference
  UA_WperK: number;          // UA product (W/K)
  Ft: number;                // LMTD correction factor
}

/**
 * Calculate LMTD for counter-current flow.
 */
function calcLMTD(Th_in: number, Th_out: number, Tc_in: number, Tc_out: number): number {
  const dT1 = Th_in - Tc_out;
  const dT2 = Th_out - Tc_in;
  if (dT1 <= 0 || dT2 <= 0) return NaN; // Temperature cross — infeasible
  if (Math.abs(dT1 - dT2) < 0.01) return dT1; // Equal — avoid ln(1)=0 division
  return (dT1 - dT2) / Math.log(dT1 / dT2);
}

/**
 * LMTD correction factor for 1-shell, 2-tube pass exchanger.
 * R = (Th_in-Th_out)/(Tc_out-Tc_in), P = (Tc_out-Tc_in)/(Th_in-Tc_in)
 */
function calcFt_1shell2tube(R: number, P: number): number {
  if (P <= 0 || P >= 1) return 1;
  if (Math.abs(R - 1) < 1e-6) {
    // Special case R=1
    const sqrt2 = Math.SQRT2;
    const num = P * sqrt2;
    const d1 = 2 - P * (2 - sqrt2);
    const d2 = 2 - P * (2 + sqrt2);
    if (d1 <= 0 || d2 <= 0) return NaN;
    const denom = (1 - P) * Math.log(d1 / d2);
    return denom !== 0 ? num / denom : 1;
  }
  const S = Math.sqrt(R * R + 1);
  const num = S * Math.log((1 - P) / (1 - R * P));
  const d1 = 2 - P * (R + 1 - S);
  const d2 = 2 - P * (R + 1 + S);
  if (d1 <= 0 || d2 <= 0) return NaN;
  const denom = (R - 1) * Math.log(d1 / d2);
  return denom !== 0 ? num / denom : 1;
}

export function solveHeatX(
  compounds: CanopyCompound[],
  hotInlet: MaterialStream,
  coldInlet: MaterialStream,
  params: {
    spec: 'hotOutletT' | 'coldOutletT' | 'duty' | 'UA';
    specValue: number;  // T in K, duty in W, UA in W/K
    hotOutletP_Pa?: number;
    coldOutletP_Pa?: number;
    flowArrangement?: 'counter' | 'co' | '1shell2tube';
    Ft_override?: number; // manual correction factor
  },
  thermo: ThermoOptions = DEFAULT_THERMO,
): HeatXResult {
  const hotP = params.hotOutletP_Pa ?? hotInlet.P_Pa;
  const coldP = params.coldOutletP_Pa ?? coldInlet.P_Pa;
  const z_hot = [...hotInlet.moleFractions];
  const z_cold = [...coldInlet.moleFractions];

  let duty_W: number;
  let Th_out: number;
  let Tc_out: number;

  switch (params.spec) {
    case 'hotOutletT': {
      Th_out = params.specValue;
      const flashHot = flashPT(compounds, z_hot, Th_out, hotP,
        thermo.fluidPackage, thermo.interactionParams);
      duty_W = hotInlet.totalFlow_molps * (hotInlet.H_Jpmol - flashHot.H_Jpmol);
      const flashCold = flashPQ(
        compounds, z_cold, coldP, coldInlet.H_Jpmol,
        duty_W, coldInlet.totalFlow_molps,
        thermo.fluidPackage, thermo.interactionParams, coldInlet.T_K
      );
      Tc_out = flashCold.T_K;
      break;
    }
    case 'coldOutletT': {
      Tc_out = params.specValue;
      const flashCold = flashPT(compounds, z_cold, Tc_out, coldP,
        thermo.fluidPackage, thermo.interactionParams);
      duty_W = coldInlet.totalFlow_molps * (flashCold.H_Jpmol - coldInlet.H_Jpmol);
      const flashHot = flashPQ(
        compounds, z_hot, hotP, hotInlet.H_Jpmol,
        -duty_W, hotInlet.totalFlow_molps,
        thermo.fluidPackage, thermo.interactionParams, hotInlet.T_K
      );
      Th_out = flashHot.T_K;
      break;
    }
    case 'duty': {
      duty_W = params.specValue;
      const flashHot = flashPQ(
        compounds, z_hot, hotP, hotInlet.H_Jpmol,
        -duty_W, hotInlet.totalFlow_molps,
        thermo.fluidPackage, thermo.interactionParams, hotInlet.T_K
      );
      Th_out = flashHot.T_K;
      const flashCold = flashPQ(
        compounds, z_cold, coldP, coldInlet.H_Jpmol,
        duty_W, coldInlet.totalFlow_molps,
        thermo.fluidPackage, thermo.interactionParams, coldInlet.T_K
      );
      Tc_out = flashCold.T_K;
      break;
    }
    case 'UA': {
      const UA = params.specValue;
      Th_out = (hotInlet.T_K + coldInlet.T_K) / 2;
      Tc_out = coldInlet.T_K;

      for (let iter = 0; iter < 50; iter++) {
        const flashHot = flashPT(compounds, z_hot, Th_out, hotP,
          thermo.fluidPackage, thermo.interactionParams);
        duty_W = hotInlet.totalFlow_molps * (hotInlet.H_Jpmol - flashHot.H_Jpmol);
        const flashCold = flashPQ(
          compounds, z_cold, coldP, coldInlet.H_Jpmol,
          duty_W, coldInlet.totalFlow_molps,
          thermo.fluidPackage, thermo.interactionParams, coldInlet.T_K
        );
        Tc_out = flashCold.T_K;

        const lmtd = calcLMTD(hotInlet.T_K, Th_out, coldInlet.T_K, Tc_out);
        if (!isFinite(lmtd) || lmtd <= 0) break;
        const Ft = params.Ft_override ?? 1;
        const Q_calc = UA * Ft * lmtd;
        const err = duty_W - Q_calc;
        if (Math.abs(err) < Math.abs(duty_W) * 1e-6 + 1) break;

        Th_out += err / (hotInlet.totalFlow_molps * 50 + 1e-10);
        Th_out = Math.max(coldInlet.T_K + 1, Math.min(hotInlet.T_K, Th_out));
      }
      duty_W = hotInlet.totalFlow_molps * (hotInlet.H_Jpmol - flashPT(compounds, z_hot, Th_out, hotP,
        thermo.fluidPackage, thermo.interactionParams).H_Jpmol);
      break;
    }
    default:
      duty_W = 0;
      Th_out = hotInlet.T_K;
      Tc_out = coldInlet.T_K;
  }

  // Final flashes
  const flashHotFinal = flashPT(compounds, z_hot, Th_out, hotP,
    thermo.fluidPackage, thermo.interactionParams);
  const flashColdFinal = flashPT(compounds, z_cold, Tc_out, coldP,
    thermo.fluidPackage, thermo.interactionParams);

  // LMTD and UA
  const LMTD = calcLMTD(hotInlet.T_K, Th_out, coldInlet.T_K, Tc_out);
  let Ft = params.Ft_override ?? 1;
  if (!params.Ft_override && params.flowArrangement === '1shell2tube') {
    const dTc = Tc_out - coldInlet.T_K;
    const dTh = hotInlet.T_K - Th_out;
    if (dTc > 0 && dTh > 0) {
      const R_hx = dTh / dTc;
      const P_hx = dTc / (hotInlet.T_K - coldInlet.T_K);
      const ft = calcFt_1shell2tube(R_hx, P_hx);
      if (isFinite(ft) && ft > 0) Ft = ft;
    }
  }
  const UA = (isFinite(LMTD) && LMTD > 0 && Ft > 0) ? Math.abs(duty_W) / (Ft * LMTD) : 0;

  return {
    hotOutlet: {
      T_K: Th_out,
      P_Pa: hotP,
      totalFlow_molps: hotInlet.totalFlow_molps,
      moleFractions: z_hot,
      phase: flashHotFinal.phase,
      vaporFraction: flashHotFinal.vaporFraction,
      H_Jpmol: flashHotFinal.H_Jpmol,
      x_liquid: flashHotFinal.x,
      y_vapor: flashHotFinal.y,
      solved: true,
    },
    coldOutlet: {
      T_K: Tc_out,
      P_Pa: coldP,
      totalFlow_molps: coldInlet.totalFlow_molps,
      moleFractions: z_cold,
      phase: flashColdFinal.phase,
      vaporFraction: flashColdFinal.vaporFraction,
      H_Jpmol: flashColdFinal.H_Jpmol,
      x_liquid: flashColdFinal.x,
      y_vapor: flashColdFinal.y,
      solved: true,
    },
    duty_W: Math.abs(duty_W),
    LMTD_K: isFinite(LMTD) ? LMTD : 0,
    UA_WperK: UA,
    Ft,
  };
}

// ────────────────────────────────────────────────────────────────
// Pipe: pressure drop via Darcy-Weisbach + Moody friction factor
// ────────────────────────────────────────────────────────────────
//
// ΔP_fric = f · (L/D) · (ρ·v²/2)            [Pa]
// ΔP_elev = ρ · g · Δz                       [Pa]
// ΔP_total = ΔP_fric + ΔP_elev + ΔP_fittings
//
// Moody friction factor (Colebrook-White):
//   1/√f = -2·log10(ε/(3.7·D) + 2.51/(Re·√f))
// Approximation (Swamee-Jain):
//   f = 0.25 / [log10(ε/(3.7·D) + 5.74/Re^0.9)]²

export interface PipeResult {
  outlet: Partial<MaterialStream>;
  deltaP_Pa: number;        // Total pressure drop
  deltaP_fric_Pa: number;   // Friction contribution
  deltaP_elev_Pa: number;   // Elevation contribution
  velocity_mps: number;     // Flow velocity
  reynoldsNumber: number;
  frictionFactor: number;
}

export function solvePipe(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  params: {
    length_m: number;
    diameter_m: number;
    roughness_m?: number;   // absolute roughness (default: 0.046e-3 for commercial steel)
    elevation_m?: number;   // Δz (positive = uphill)
    K_fittings?: number;    // Sum of fitting resistance coefficients
  },
  thermo: ThermoOptions = DEFAULT_THERMO,
): PipeResult {
  const { length_m, diameter_m } = params;
  const roughness = params.roughness_m ?? 0.046e-3;
  const elevation = params.elevation_m ?? 0;
  const K_fittings = params.K_fittings ?? 0;
  const z = [...inlet.moleFractions];

  // Get density and viscosity at inlet conditions
  const MW_avg = compounds.reduce((s, c, i) => s + z[i] * c.molecularWeight, 0); // g/mol
  const V_L = mixtureLiquidMolarVolume(compounds, z, inlet.T_K);
  const V_L_safe = isFinite(V_L) && V_L > 0 ? V_L : 1e-4;
  const rhoMolar = 1 / V_L_safe; // mol/m³
  const rhoMass = rhoMolar * MW_avg / 1000; // kg/m³

  // Pipe cross-section area
  const A_pipe = Math.PI * diameter_m * diameter_m / 4;

  // Volumetric flow (m³/s)
  const Q_vol = inlet.totalFlow_molps * V_L_safe; // mol/s · m³/mol

  // Velocity
  const velocity = A_pipe > 0 ? Q_vol / A_pipe : 0;

  // Reynolds number (assume liquid viscosity ~1e-3 Pa·s if not available)
  const mu = 1e-3; // Pa·s (approximate — should come from DIPPR)
  const Re = rhoMass * velocity * diameter_m / mu;

  // Friction factor — Swamee-Jain approximation (valid for Re > 4000)
  let f: number;
  if (Re < 2300) {
    // Laminar
    f = Re > 0 ? 64 / Re : 0;
  } else {
    // Turbulent — Swamee-Jain
    const term = roughness / (3.7 * diameter_m) + 5.74 / Math.pow(Re, 0.9);
    f = term > 0 ? 0.25 / Math.pow(Math.log10(term), 2) : 0.02;
  }

  const g = 9.80665;

  // Darcy-Weisbach: ΔP_fric = f · (L/D) · ρv²/2
  const dP_fric = f * (length_m / diameter_m) * rhoMass * velocity * velocity / 2;

  // Fittings: ΔP_fit = K · ρv²/2
  const dP_fit = K_fittings * rhoMass * velocity * velocity / 2;

  // Elevation: ΔP_elev = ρ·g·Δz
  const dP_elev = rhoMass * g * elevation;

  const dP_total = dP_fric + dP_fit + dP_elev;
  const P_out = Math.max(1000, inlet.P_Pa - dP_total); // Floor at 1 kPa

  // Isenthalpic flash at new pressure
  const flash = flashPQ(
    compounds, z, P_out, inlet.H_Jpmol, 0, inlet.totalFlow_molps,
    thermo.fluidPackage, thermo.interactionParams, inlet.T_K
  );

  return {
    outlet: {
      T_K: flash.T_K,
      P_Pa: P_out,
      totalFlow_molps: inlet.totalFlow_molps,
      moleFractions: z,
      phase: flash.phase,
      vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    deltaP_Pa: dP_total,
    deltaP_fric_Pa: dP_fric + dP_fit,
    deltaP_elev_Pa: dP_elev,
    velocity_mps: velocity,
    reynoldsNumber: Re,
    frictionFactor: f,
  };
}

// ────────────────────────────────────────────────────────────────
// RStoic: Stoichiometric reactor
// ────────────────────────────────────────────────────────────────
//
// Given reactions with coefficients ν_i,r and conversion X_r:
//   ξ_r = X_r · F_{key,in} / |ν_{key,r}|    (extent of reaction)
//   F_i,out = F_i,in + Σ_r ν_i,r · ξ_r
//
// Then flash the product at specified (T, P).

export interface RStoicResult {
  outlet: Partial<MaterialStream>;
  duty_W: number;
  extents: number[];         // Extent of each reaction (mol/s)
  componentFlows: number[];  // Outlet component flows (mol/s)
}

export function solveRStoic(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  params: {
    reactions: StoichiometricReaction[];
    outletT_K: number;
    outletP_Pa: number;
  },
  thermo: ThermoOptions = DEFAULT_THERMO,
): RStoicResult {
  const n = compounds.length;
  const F_in = inlet.moleFractions.map(x => x * inlet.totalFlow_molps);
  const F_out = [...F_in];
  const extents: number[] = [];

  for (const rxn of params.reactions) {
    // Extent = conversion × current flow of key / |ν_key|
    const F_key = F_out[rxn.keyComponentIndex];
    const nu_key = rxn.coefficients[rxn.keyComponentIndex];
    const extent = Math.abs(nu_key) > 0
      ? rxn.conversion * F_key / Math.abs(nu_key)
      : 0;
    extents.push(extent);

    // Apply to all components
    for (let i = 0; i < n; i++) {
      F_out[i] += rxn.coefficients[i] * extent;
      if (F_out[i] < 0) F_out[i] = 0; // Prevent negative flows
    }
  }

  const totalFlow = F_out.reduce((a, b) => a + b, 0);
  const z_out = totalFlow > 0 ? F_out.map(f => f / totalFlow) : new Array(n).fill(1 / n);

  // Flash product
  const flash = flashPT(compounds, z_out, params.outletT_K, params.outletP_Pa,
    thermo.fluidPackage, thermo.interactionParams);

  // Duty = H_out - H_in (energy balance)
  const H_in_W = inlet.totalFlow_molps * inlet.H_Jpmol;
  const H_out_W = totalFlow * flash.H_Jpmol;
  const duty_W = H_out_W - H_in_W;

  return {
    outlet: {
      T_K: params.outletT_K,
      P_Pa: params.outletP_Pa,
      totalFlow_molps: totalFlow,
      moleFractions: z_out,
      phase: flash.phase,
      vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    duty_W,
    extents,
    componentFlows: F_out,
  };
}

// ────────────────────────────────────────────────────────────────
// Column (Shortcut): Fenske-Underwood-Gilliland distillation
// ────────────────────────────────────────────────────────────────
//
// Fenske: N_min = ln[(xLK_D/xHK_D)·(xHK_B/xLK_B)] / ln(α_avg)
// Underwood: Σ αi·zi / (αi - θ) = 1-q  → find θ → R_min
// Gilliland: (N-Nmin)/(N+1) = f((R-Rmin)/(R+1))
//
// This gives distillate and bottoms compositions assuming sharp separation
// of light and heavy keys with specified recoveries.

export interface ColumnResult {
  distillate: Partial<MaterialStream>;
  bottoms: Partial<MaterialStream>;
  duty_condenser_W: number;
  duty_reboiler_W: number;
  N_min: number;           // Minimum stages (Fenske)
  N_actual: number;        // Actual stages (Gilliland)
  R_min: number;           // Minimum reflux ratio (Underwood)
  R_actual: number;        // Actual reflux ratio
  feedStage: number;       // Optimal feed stage (Kirkbride)
}

export function solveColumn(
  compounds: CanopyCompound[],
  feed: MaterialStream,
  spec: ColumnSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): ColumnResult {
  const n = compounds.length;
  const F_feed = feed.moleFractions.map(x => x * feed.totalFlow_molps); // component flows
  const z = [...feed.moleFractions];

  const { lightKeyIndex: lk, heavyKeyIndex: hk } = spec;
  const { lightKeyRecovery: recLK, heavyKeyRecovery: recHK } = spec;

  // Relative volatilities at feed, condenser, reboiler temperatures
  // Use average of condenser and reboiler conditions
  const T_feed = feed.T_K;
  const P_avg = (spec.condenserP_Pa + spec.reboilerP_Pa) / 2;

  const K_feed = KvaluesIdeal(compounds, T_feed, P_avg);
  const alpha_feed = K_feed.map(k => k / K_feed[hk]); // α relative to heavy key

  // Estimate condenser/reboiler T from bubble points (crude)
  const T_cond = bubblePointT_Ideal(compounds, z, spec.condenserP_Pa) ?? T_feed - 20;
  const T_reb = bubblePointT_Ideal(compounds, z, spec.reboilerP_Pa) ?? T_feed + 20;

  const K_cond = KvaluesIdeal(compounds, T_cond, spec.condenserP_Pa);
  const K_reb = KvaluesIdeal(compounds, T_reb, spec.reboilerP_Pa);
  const alpha_cond = K_cond.map(k => k / K_cond[hk]);
  const alpha_reb = K_reb.map(k => k / K_reb[hk]);

  // Geometric mean relative volatility
  const alpha = alpha_cond.map((_, i) =>
    Math.pow(alpha_cond[i] * alpha_feed[i] * alpha_reb[i], 1 / 3)
  );
  const alphaLK = alpha[lk];

  // --- Fenske: Minimum stages ---
  const xLK_D = recLK * z[lk]; // approx mole fraction in distillate
  const xHK_D = (1 - recHK) * z[hk];
  const xLK_B = (1 - recLK) * z[lk];
  const xHK_B = recHK * z[hk];

  let N_min: number;
  if (alphaLK > 1 && xLK_D > 0 && xHK_D > 0 && xLK_B > 0 && xHK_B > 0) {
    N_min = Math.log((xLK_D / xHK_D) * (xHK_B / xLK_B)) / Math.log(alphaLK);
  } else {
    N_min = 10; // fallback
  }
  N_min = Math.max(1, N_min);

  // --- Underwood: Minimum reflux ---
  // Solve Σ αi·zi/(αi-θ) = 1-q for θ (q=1 for saturated liquid feed)
  const q = feed.vaporFraction < 0.01 ? 1 : (1 - feed.vaporFraction);

  // Find θ between α_HK and α_LK via bisection
  let thetaLo = 1.001; // just above α_HK (=1)
  let thetaHi = alphaLK - 0.001;
  if (thetaLo >= thetaHi) thetaHi = thetaLo + 1;

  const underwoodF = (theta: number) =>
    z.reduce((s, zi, i) => s + alpha[i] * zi / (alpha[i] - theta), 0) - (1 - q);

  // Bisection
  for (let iter = 0; iter < 100; iter++) {
    const mid = (thetaLo + thetaHi) / 2;
    const fmid = underwoodF(mid);
    const flo = underwoodF(thetaLo);
    if (Math.abs(fmid) < 1e-8) { thetaLo = mid; break; }
    if (fmid * flo < 0) thetaHi = mid; else thetaLo = mid;
    if (thetaHi - thetaLo < 1e-8) break;
  }
  const theta = (thetaLo + thetaHi) / 2;

  // R_min+1 from Underwood (distillate composition)
  // Build approximate distillate composition using Fenske distribution
  // For non-key component i: d_i/b_i = α_i^N_min
  // → fraction to distillate = α_i^N_min / (1 + α_i^N_min)
  const fenskeDistFrac = (alphaI: number) => {
    const aN = Math.pow(alphaI, N_min);
    return aN / (1 + aN); // fraction of component i going to distillate
  };

  const D_total = F_feed.reduce((s, fi, i) => {
    if (i === lk) return s + recLK * fi;
    if (i === hk) return s + (1 - recHK) * fi;
    return s + fenskeDistFrac(alpha[i]) * fi;
  }, 0);

  const xD = F_feed.map((fi, i) => {
    if (i === lk) return recLK * fi / D_total;
    if (i === hk) return (1 - recHK) * fi / D_total;
    return fenskeDistFrac(alpha[i]) * fi / D_total;
  });
  // Normalize xD
  const xD_sum = xD.reduce((a, b) => a + b, 0);
  for (let i = 0; i < n; i++) xD[i] /= xD_sum;

  let Rmin_p1 = xD.reduce((s, xdi, i) => s + alpha[i] * xdi / (alpha[i] - theta), 0);
  let R_min = Math.max(0.1, Rmin_p1 - 1);

  // --- Gilliland correlation: (N-Nmin)/(N+1) = f(X) ---
  const R_actual = spec.refluxRatioMultiplier * R_min;
  const X = (R_actual - R_min) / (R_actual + 1);

  // Molokanov correlation for Gilliland:
  //   Y = 1 - exp[(1+54.4X)/(11+117.2X) · (X-1)/√X]
  let Y: number;
  if (X <= 0) {
    Y = 1; // operating at minimum reflux
  } else if (X >= 1) {
    Y = 0; // total reflux
  } else {
    Y = 1 - Math.exp(((1 + 54.4 * X) / (11 + 117.2 * X)) * ((X - 1) / Math.sqrt(X)));
  }
  const N_actual = Math.max(Math.ceil((Y * (N_min + 1) + N_min) / (1 - Y + 1e-10)), Math.ceil(N_min) + 1);

  // --- Kirkbride: feed stage location ---
  const B_total = feed.totalFlow_molps - D_total;
  let feedStage: number;
  if (B_total > 0 && D_total > 0 && xD[hk] > 0 && z[lk] > 0) {
    const ratio = Math.pow(
      (B_total / D_total) * (z[hk] / z[lk]) * Math.pow(xD[hk] / (xHK_B / Math.max(B_total, 1e-10)), 2),
      0.206
    );
    feedStage = Math.max(1, Math.min(N_actual - 1, Math.round(N_actual / (1 + ratio))));
  } else {
    feedStage = Math.round(N_actual / 2);
  }

  // --- Material balance ---
  const F_D = new Array(n).fill(0);
  const F_B = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    F_D[i] = xD[i] * D_total;
    F_B[i] = F_feed[i] - F_D[i];
    if (F_B[i] < 0) { F_D[i] = F_feed[i]; F_B[i] = 0; }
  }
  const D_flow = F_D.reduce((a, b) => a + b, 0);
  const B_flow = F_B.reduce((a, b) => a + b, 0);
  const zD = D_flow > 0 ? F_D.map(f => f / D_flow) : new Array(n).fill(1 / n);
  const zB = B_flow > 0 ? F_B.map(f => f / B_flow) : new Array(n).fill(1 / n);

  // Flash distillate and bottoms
  const flashD = flashPT(compounds, zD, T_cond, spec.condenserP_Pa,
    thermo.fluidPackage, thermo.interactionParams);
  const flashB = flashPT(compounds, zB, T_reb, spec.reboilerP_Pa,
    thermo.fluidPackage, thermo.interactionParams);

  // Energy balance (simplified)
  const H_feed = feed.totalFlow_molps * feed.H_Jpmol;
  const H_dist = D_flow * flashD.H_Jpmol;
  const H_bot = B_flow * flashB.H_Jpmol;
  // Q_cond + Q_reb = H_dist + H_bot - H_feed
  // Approximate: Q_cond ≈ -R·D·(H_V - H_L) at condenser, Q_reb = H_dist + H_bot - H_feed + Q_cond
  const H_V_cond = streamEnthalpy(compounds, T_cond, 1, flashD.x, flashD.y);
  const H_L_cond = streamEnthalpy(compounds, T_cond, 0, flashD.x, flashD.y);
  const Q_cond = -(R_actual + (spec.condenserType === 'Total' ? 1 : 0)) * D_flow * (H_V_cond - H_L_cond);
  const Q_reb = H_dist + H_bot - H_feed - Q_cond;

  return {
    distillate: {
      T_K: T_cond,
      P_Pa: spec.condenserP_Pa,
      totalFlow_molps: D_flow,
      moleFractions: zD,
      phase: flashD.phase,
      vaporFraction: spec.condenserType === 'Total' ? 0 : flashD.vaporFraction,
      H_Jpmol: flashD.H_Jpmol,
      x_liquid: flashD.x,
      y_vapor: flashD.y,
      solved: true,
    },
    bottoms: {
      T_K: T_reb,
      P_Pa: spec.reboilerP_Pa,
      totalFlow_molps: B_flow,
      moleFractions: zB,
      phase: flashB.phase,
      vaporFraction: 0,
      H_Jpmol: flashB.H_Jpmol,
      x_liquid: flashB.x,
      y_vapor: flashB.y,
      solved: true,
    },
    duty_condenser_W: Q_cond,
    duty_reboiler_W: Q_reb,
    N_min,
    N_actual,
    R_min,
    R_actual,
    feedStage,
  };
}

// ────────────────────────────────────────────────────────────────
// ComponentSeparator: split feed into products by specified split fractions
// ────────────────────────────────────────────────────────────────

export function solveComponentSeparator(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  /** splitFractions[i] = fraction of component i going to outlet1 (0-1) */
  splitFractions: number[],
  thermo: ThermoOptions = DEFAULT_THERMO,
): { outlet1: MaterialStream; outlet2: MaterialStream } {
  const n = compounds.length;
  const F = inlet.totalFlow_molps;

  const z1 = new Array(n).fill(0);
  const z2 = new Array(n).fill(0);
  let f1 = 0, f2 = 0;

  for (let i = 0; i < n; i++) {
    const fi = F * inlet.moleFractions[i];
    const split = splitFractions[i] ?? 0.5;
    z1[i] = fi * split;
    z2[i] = fi * (1 - split);
    f1 += z1[i];
    f2 += z2[i];
  }

  // Normalize mole fractions
  for (let i = 0; i < n; i++) {
    z1[i] = f1 > 0 ? z1[i] / f1 : 0;
    z2[i] = f2 > 0 ? z2[i] / f2 : 0;
  }

  const flash1 = flashPT(compounds, z1, inlet.T_K, inlet.P_Pa, thermo.fluidPackage, thermo.interactionParams);
  const flash2 = flashPT(compounds, z2, inlet.T_K, inlet.P_Pa, thermo.fluidPackage, thermo.interactionParams);

  return {
    outlet1: {
      ...inlet,
      T_K: inlet.T_K,
      P_Pa: inlet.P_Pa,
      totalFlow_molps: f1,
      moleFractions: z1,
      phase: flash1.phase,
      vaporFraction: flash1.vaporFraction,
      H_Jpmol: flash1.H_Jpmol,
      x_liquid: flash1.x,
      y_vapor: flash1.y,
      solved: true,
    },
    outlet2: {
      ...inlet,
      T_K: inlet.T_K,
      P_Pa: inlet.P_Pa,
      totalFlow_molps: f2,
      moleFractions: z2,
      phase: flash2.phase,
      vaporFraction: flash2.vaporFraction,
      H_Jpmol: flash2.H_Jpmol,
      x_liquid: flash2.x,
      y_vapor: flash2.y,
      solved: true,
    },
  };
}

// ────────────────────────────────────────────────────────────────
// CSTR: Continuous Stirred Tank Reactor
// Supports two modes:
//   (a) Conversion specification (stoichiometric, like RStoic)
//   (b) Kinetic rate law: r = k₀·T^n·exp(-Ea/RT)·Π[Ci]^{αi}
// ────────────────────────────────────────────────────────────────

export interface KineticReaction {
  coefficients: number[];       // stoichiometric coefficients (negative = reactant)
  keyComponentIndex: number;
  conversion: number;           // used only in conversion mode
  // Kinetic parameters (if provided, uses kinetic mode)
  preExponentialFactor?: number; // k₀ (units depend on reaction order)
  activationEnergy_Jpmol?: number; // Ea (J/mol)
  temperatureExponent?: number;  // n in T^n
  reactionOrders?: number[];     // order w.r.t. each component (default = |ν_i| for reactants)
  /** Reaction phase: 'Vapor' uses C_i = y_i·P/(RT), 'Liquid' uses C_i = x_i·ρ_L/M_w */
  phase?: 'Vapor' | 'Liquid';
}

export function solveCSTR(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  params: {
    reactions: KineticReaction[];
    outletT_K: number;
    outletP_Pa: number;
    volume_m3?: number;          // reactor volume (needed for kinetic mode)
  },
  /** Heat duty supplied to reactor (W), 0 for adiabatic */
  duty_W: number = 0,
  thermo: ThermoOptions = DEFAULT_THERMO,
): { outlet: Partial<MaterialStream>; duty_W: number } {
  const nc = compounds.length;
  const F_in = inlet.moleFractions.map(x => x * inlet.totalFlow_molps);
  const F_out = [...F_in];
  const R_gas = 8.314462618;

  const hasKinetics = params.reactions.some(r => r.preExponentialFactor !== undefined);

  if (hasKinetics && params.volume_m3) {
    // ── Kinetic mode: solve steady-state CSTR equations ──
    // 0 = F_in,i − F_out,i + V·Σⱼ(νᵢⱼ·rⱼ)
    // Use successive substitution with damping
    const V = params.volume_m3;
    const T = params.outletT_K ?? inlet.T_K;
    const P = params.outletP_Pa ?? inlet.P_Pa;

    // Estimate total volumetric flow (ideal gas or liquid)
    // Determine reaction phase: defaults to Vapor (ideal gas) unless Liquid specified
    const rxnPhase = params.reactions[0]?.phase ?? 'Vapor';
    for (let iter = 0; iter < 500; iter++) {
      const F_total = F_out.reduce((a, b) => a + b, 0);
      const z_iter = F_out.map(f => F_total > 0 ? f / F_total : 0);

      let C: number[];
      if (rxnPhase === 'Liquid') {
        // Liquid: C_i = x_i · ρ_L / M_w  (mol/m³)
        // ρ_L in kmol/m³ from DIPPR, average Mw for mixture
        const Mw_mix = z_iter.reduce((s, xi, i) => s + xi * compounds[i].molecularWeight, 0); // g/mol
        let rho_kmolm3 = 0;
        for (let i = 0; i < nc; i++) {
          const rho_i = liquidDensity_kmolpm3(compounds[i], T);
          rho_kmolm3 += z_iter[i] * (isFinite(rho_i) && rho_i > 0 ? rho_i : 10);
        }
        // C_i = x_i * rho_mol (mol/m³) where rho_mol = rho_kmolm3 * 1000
        C = z_iter.map(xi => xi * rho_kmolm3 * 1000);
      } else {
        // Vapor (ideal gas): C_i = F_out,i / Q  where Q = F_total·RT/P
        const Q = F_total * R_gas * T / P; // m³/s
        C = F_out.map(f => Q > 0 ? f / Q : 0); // mol/m³
      }

      const F_new = [...F_in];
      for (const rxn of params.reactions) {
        const k0 = rxn.preExponentialFactor ?? 1;
        const Ea = rxn.activationEnergy_Jpmol ?? 0;
        const n_T = rxn.temperatureExponent ?? 0;
        const k = k0 * Math.pow(T, n_T) * Math.exp(-Ea / (R_gas * T));

        // Rate: r = k · Π[C_i]^{order_i}
        let rate = k;
        for (let i = 0; i < nc; i++) {
          const order = rxn.reactionOrders?.[i]
            ?? (rxn.coefficients[i] < 0 ? Math.abs(rxn.coefficients[i]) : 0);
          if (order > 0) {
            rate *= Math.pow(Math.max(C[i], 0), order);
          }
        }

        for (let i = 0; i < nc; i++) {
          F_new[i] += rxn.coefficients[i] * rate * V;
        }
      }

      // Damped update and enforce non-negativity
      let maxChange = 0;
      for (let i = 0; i < nc; i++) {
        const newVal = Math.max(0, 0.5 * F_out[i] + 0.5 * F_new[i]);
        maxChange = Math.max(maxChange, Math.abs(newVal - F_out[i]) / (F_out[i] + 1e-30));
        F_out[i] = newVal;
      }
      if (maxChange < 1e-10) break;
    }
  } else {
    // ── Conversion specification mode (same as RStoic) ──
    for (const rxn of params.reactions) {
      const keyMoles = F_in[rxn.keyComponentIndex];
      const extent = keyMoles * rxn.conversion;
      for (let i = 0; i < nc; i++) {
        const stoichRatio = rxn.coefficients[i] / Math.abs(rxn.coefficients[rxn.keyComponentIndex]);
        F_out[i] += stoichRatio * extent;
        if (F_out[i] < 0) F_out[i] = 0;
      }
    }
  }

  const totalFlow = F_out.reduce((s, m) => s + m, 0);
  const z = F_out.map(m => totalFlow > 0 ? m / totalFlow : 0);
  const T_out = params.outletT_K ?? inlet.T_K;
  const P_out = params.outletP_Pa ?? inlet.P_Pa;

  // Energy balance
  const H_in = inlet.totalFlow_molps * inlet.H_Jpmol;

  if (Math.abs(duty_W) < 1e-10) {
    // Adiabatic: solve for T such that H_out = H_in
    const H_target = totalFlow > 0 ? H_in / totalFlow : inlet.H_Jpmol;
    const flash = flashPQ(compounds, z, P_out, H_target, 0, totalFlow,
      thermo.fluidPackage, thermo.interactionParams, inlet.T_K);
    return {
      outlet: {
        T_K: flash.T_K, P_Pa: P_out, totalFlow_molps: totalFlow,
        moleFractions: z, phase: flash.phase, vaporFraction: flash.vaporFraction,
        H_Jpmol: flash.H_Jpmol, x_liquid: flash.x, y_vapor: flash.y, solved: flash.converged,
      },
      duty_W: 0,
    };
  }

  // With specified outlet T: compute duty
  const flash = flashPT(compounds, z, T_out, P_out, thermo.fluidPackage, thermo.interactionParams);
  const H_out = totalFlow * flash.H_Jpmol;

  return {
    outlet: {
      T_K: T_out, P_Pa: P_out, totalFlow_molps: totalFlow,
      moleFractions: z, phase: flash.phase, vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol, x_liquid: flash.x, y_vapor: flash.y, solved: flash.converged,
    },
    duty_W: H_out - H_in,
  };
}

// ────────────────────────────────────────────────────────────────
// PFR: Plug Flow Reactor
// Supports two modes:
//   (a) Conversion specification (stoichiometric)
//   (b) Kinetic rate law with RK4 integration along reactor volume
// ────────────────────────────────────────────────────────────────

export function solvePFR(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  params: {
    reactions: KineticReaction[];
    outletT_K: number;
    outletP_Pa: number;
    volume_m3?: number;
    /** Tube diameter (m) for heat transfer */
    diameter_m?: number;
    /** Overall heat transfer coefficient U (W/m²/K) */
    U_Wpm2K?: number;
    /** Coolant temperature (K) — constant or initial */
    T_coolant_K?: number;
    /** Coolant flow · Cp product (W/K) for non-isothermal coolant, undefined = isothermal coolant */
    coolantMdotCp_WpK?: number;
    /** Coolant flow direction: true = co-current, false = counter-current */
    coCurrent?: boolean;
  },
  /** Number of integration steps along reactor length */
  nSteps: number = 50,
  /** Heat duty supplied to reactor (W), 0 for adiabatic */
  duty_W: number = 0,
  thermo: ThermoOptions = DEFAULT_THERMO,
): { outlet: Partial<MaterialStream>; duty_W: number } {
  const nc = compounds.length;
  const F_in = inlet.moleFractions.map(x => x * inlet.totalFlow_molps);
  const R_gas = 8.314462618;
  const F_out = [...F_in];

  const hasKinetics = params.reactions.some(r => r.preExponentialFactor !== undefined);

  if (hasKinetics && params.volume_m3) {
    // ── Kinetic mode: integrate dFi/dV = Σ(νij · rj) via RK4 ──
    // Aspen RPlug: coupled material + energy balance ODEs
    // dF_i/dV = r_i
    // dT/dV = [−Σ(ΔH_r · R_r) + U·πD·(Tc − T)/Ac] / (F_total · Cp_mix)
    const V_total = params.volume_m3;
    const dV = V_total / nSteps;
    const P = params.outletP_Pa ?? inlet.P_Pa;

    // Heat transfer parameters
    const D = params.diameter_m ?? 0;
    const U = params.U_Wpm2K ?? 0;
    const hasHeatTransfer = D > 0 && U > 0 && params.T_coolant_K !== undefined;
    const Ac = D > 0 ? Math.PI * D * D / 4 : 1; // cross-sectional area

    // State vector: [F_0, F_1, ..., F_nc-1, T, T_coolant]
    const stateSize = nc + 2;

    const computeDerivs = (state: number[]): number[] => {
      const F = state.slice(0, nc);
      const T_loc = state[nc];
      const Tc_loc = state[nc + 1];
      const F_total = F.reduce((a, b) => a + b, 0);
      const z_loc = F.map(f => F_total > 0 ? f / F_total : 0);
      const rxnPhase = params.reactions[0]?.phase ?? 'Vapor';

      let C: number[];
      if (rxnPhase === 'Liquid') {
        let rho_kmolm3 = 0;
        for (let i = 0; i < nc; i++) {
          const rho_i = liquidDensity_kmolpm3(compounds[i], T_loc);
          rho_kmolm3 += z_loc[i] * (isFinite(rho_i) && rho_i > 0 ? rho_i : 10);
        }
        C = z_loc.map(xi => xi * rho_kmolm3 * 1000);
      } else {
        const Q = F_total * R_gas * T_loc / P;
        C = F.map(f => Q > 0 ? f / Q : 0);
      }

      const dFdV = new Array(nc).fill(0);
      let sumHrRr = 0; // Σ(ΔH_r · R_r) for energy balance

      for (const rxn of params.reactions) {
        const k0 = rxn.preExponentialFactor ?? 1;
        const Ea = rxn.activationEnergy_Jpmol ?? 0;
        const n_T = rxn.temperatureExponent ?? 0;
        const k = k0 * Math.pow(T_loc, n_T) * Math.exp(-Ea / (R_gas * T_loc));

        let rate = k;
        for (let i = 0; i < nc; i++) {
          const order = rxn.reactionOrders?.[i]
            ?? (rxn.coefficients[i] < 0 ? Math.abs(rxn.coefficients[i]) : 0);
          if (order > 0) rate *= Math.pow(Math.max(C[i], 0), order);
        }
        for (let i = 0; i < nc; i++) {
          dFdV[i] += rxn.coefficients[i] * rate;
        }

        // Heat of reaction: ΔH_rxn = Σ ν_i · H_ig,i(T) (from formation enthalpies)
        let deltaHrxn = 0;
        for (let i = 0; i < nc; i++) {
          if (rxn.coefficients[i] !== 0) {
            deltaHrxn += rxn.coefficients[i] * enthalpyIG(compounds[i], T_loc);
          }
        }
        sumHrRr += deltaHrxn * rate;
      }

      // Energy balance: dT/dV
      let dTdV = 0;
      if (hasHeatTransfer || Math.abs(sumHrRr) > 1e-20) {
        // Mixture Cp (J/mol/K)
        let Cp_mix = 0;
        for (let i = 0; i < nc; i++) {
          Cp_mix += z_loc[i] * CpIG_JmolK(compounds[i], T_loc);
        }
        Cp_mix = Math.max(Cp_mix, 1); // safety

        const Qheat = hasHeatTransfer
          ? U * Math.PI * D * (Tc_loc - T_loc) / Ac
          : 0;
        dTdV = (-sumHrRr + Qheat) / (F_total * Cp_mix);
      }

      // Coolant temperature derivative
      let dTcdV = 0;
      if (hasHeatTransfer && params.coolantMdotCp_WpK) {
        const Qwall = U * Math.PI * D * (T_loc - Tc_loc);
        const sign = params.coCurrent !== false ? 1 : -1; // co-current default
        dTcdV = sign * Qwall / params.coolantMdotCp_WpK / Ac;
      }

      const derivs = new Array(stateSize).fill(0);
      for (let i = 0; i < nc; i++) derivs[i] = dFdV[i];
      derivs[nc] = dTdV;
      derivs[nc + 1] = dTcdV;
      return derivs;
    };

    // Initial state
    const state = new Array(stateSize);
    for (let i = 0; i < nc; i++) state[i] = F_in[i];
    state[nc] = params.outletT_K ?? inlet.T_K; // initial T = inlet T
    state[nc + 1] = params.T_coolant_K ?? inlet.T_K;

    // RK4 integration
    for (let step = 0; step < nSteps; step++) {
      const k1 = computeDerivs(state);
      const s1 = state.map((s, i) => s + 0.5 * dV * k1[i]);
      const k2 = computeDerivs(s1);
      const s2 = state.map((s, i) => s + 0.5 * dV * k2[i]);
      const k3 = computeDerivs(s2);
      const s3 = state.map((s, i) => s + dV * k3[i]);
      const k4 = computeDerivs(s3);

      for (let i = 0; i < stateSize; i++) {
        state[i] += dV * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
      }
      // Enforce non-negative flows
      for (let i = 0; i < nc; i++) state[i] = Math.max(0, state[i]);
    }
    for (let i = 0; i < nc; i++) F_out[i] = state[i];

    // Use integrated outlet temperature
    const T_integrated = state[nc];
    const totalFlow = F_out.reduce((s, m) => s + m, 0);
    const z = F_out.map(m => totalFlow > 0 ? m / totalFlow : 0);
    const P_out = params.outletP_Pa ?? inlet.P_Pa;

    const flash = flashPT(compounds, z, T_integrated, P_out, thermo.fluidPackage, thermo.interactionParams);
    const H_in = inlet.totalFlow_molps * inlet.H_Jpmol;
    const H_out = totalFlow * flash.H_Jpmol;

    return {
      outlet: {
        T_K: T_integrated, P_Pa: P_out, totalFlow_molps: totalFlow,
        moleFractions: z, phase: flash.phase, vaporFraction: flash.vaporFraction,
        H_Jpmol: flash.H_Jpmol, x_liquid: flash.x, y_vapor: flash.y, solved: flash.converged,
      },
      duty_W: H_out - H_in,
    };
  } else {
    // ── Conversion specification mode ──
    for (const rxn of params.reactions) {
      const keyMolesInitial = F_in[rxn.keyComponentIndex];
      const extent = keyMolesInitial * rxn.conversion;
      for (let i = 0; i < nc; i++) {
        const stoichRatio = rxn.coefficients[i] / Math.abs(rxn.coefficients[rxn.keyComponentIndex]);
        F_out[i] += stoichRatio * extent;
        if (F_out[i] < 0) F_out[i] = 0;
      }
    }
  }

  const totalFlow = F_out.reduce((s, m) => s + m, 0);
  const z = F_out.map(m => totalFlow > 0 ? m / totalFlow : 0);
  const T_out = params.outletT_K ?? inlet.T_K;
  const P_out = params.outletP_Pa ?? inlet.P_Pa;

  const flash = flashPT(compounds, z, T_out, P_out, thermo.fluidPackage, thermo.interactionParams);
  const H_in = inlet.totalFlow_molps * inlet.H_Jpmol;
  const H_out = totalFlow * flash.H_Jpmol;

  return {
    outlet: {
      T_K: T_out, P_Pa: P_out, totalFlow_molps: totalFlow,
      moleFractions: z, phase: flash.phase, vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol, x_liquid: flash.x, y_vapor: flash.y, solved: flash.converged,
    },
    duty_W: H_out - H_in,
  };
}

// ────────────────────────────────────────────────────────────────
// ThreePhaseFlash: VLLE flash drum
// ────────────────────────────────────────────────────────────────

export function solveThreePhaseFlash(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  T_out_K: number | null,
  P_out_Pa: number | null,
  duty_W: number,
  thermo: ThermoOptions = DEFAULT_THERMO,
): {
  vapor: MaterialStream;
  liquidI: MaterialStream;
  liquidII: MaterialStream;
  duty_W: number;
} {
  const T = T_out_K ?? inlet.T_K;
  const P = P_out_Pa ?? inlet.P_Pa;
  const F = inlet.totalFlow_molps;

  const result = flashPT_VLLE(compounds, inlet.moleFractions, T, P,
    thermo.fluidPackage, thermo.interactionParams);

  const V_flow = F * result.vaporFraction;
  const LI_flow = F * result.liquidIFraction;
  const LII_flow = F * result.liquidIIFraction;

  const H_out = result.H_Jpmol;
  const Q = (H_out - inlet.H_Jpmol) * F - duty_W;

  return {
    vapor: {
      ...inlet,
      T_K: T, P_Pa: P, totalFlow_molps: V_flow,
      moleFractions: result.y,
      phase: 'Vapor', vaporFraction: 1,
      H_Jpmol: streamEnthalpy(compounds, T, 1, result.y, result.y),
      x_liquid: result.y, y_vapor: result.y, solved: true,
    },
    liquidI: {
      ...inlet,
      T_K: T, P_Pa: P, totalFlow_molps: LI_flow,
      moleFractions: result.xI,
      phase: 'Liquid', vaporFraction: 0,
      H_Jpmol: streamEnthalpy(compounds, T, 0, result.xI, result.xI),
      x_liquid: result.xI, y_vapor: result.xI, solved: true,
    },
    liquidII: {
      ...inlet,
      T_K: T, P_Pa: P, totalFlow_molps: LII_flow,
      moleFractions: result.xII,
      phase: 'Liquid', vaporFraction: 0,
      H_Jpmol: streamEnthalpy(compounds, T, 0, result.xII, result.xII),
      x_liquid: result.xII, y_vapor: result.xII, solved: true,
    },
    duty_W: Q,
  };
}

// ────────────────────────────────────────────────────────────────
// Absorber: multi-stage gas absorption (counter-current)
// Simplified model using Kremser-Brown-Souders equation
// ────────────────────────────────────────────────────────────────

export function solveAbsorber(
  compounds: CanopyCompound[],
  gasInlet: MaterialStream,
  liquidInlet: MaterialStream,
  /** Number of theoretical stages */
  nStages: number,
  P_Pa: number,
  thermo: ThermoOptions = DEFAULT_THERMO,
): { gasOutlet: MaterialStream; liquidOutlet: MaterialStream } {
  const n = compounds.length;

  // Average temperature
  const T_avg = (gasInlet.T_K + liquidInlet.T_K) / 2;

  // K-values at average conditions
  let x_est = [...liquidInlet.moleFractions];
  let y_est = [...gasInlet.moleFractions];
  let K = thermo.fluidPackage === 'Ideal'
    ? KvaluesIdeal(compounds, T_avg, P_Pa)
    : computeKvalues(compounds, x_est, y_est, T_avg, P_Pa, thermo.fluidPackage, thermo.interactionParams);

  // Absorption factors A_i = L / (K_i · V)
  const L = liquidInlet.totalFlow_molps;
  const V = gasInlet.totalFlow_molps;

  const computeSplit = (kValues: number[]) => {
    const gasOut_moles = new Array(n).fill(0);
    const liqOut_moles = new Array(n).fill(0);

    for (let i = 0; i < n; i++) {
      const Ki = isFinite(kValues[i]) && kValues[i] > 1e-30 ? kValues[i] : 1;
      const A_i = (L > 1e-30 && V > 1e-30) ? L / (Ki * V) : 0;
      const y_in = gasInlet.moleFractions[i] * V;
      const x_in = liquidInlet.moleFractions[i] * L;

      let E_abs: number;
      if (L <= 1e-30 || V <= 1e-30) {
        E_abs = 0;
      } else if (Math.abs(A_i - 1) < 1e-6) {
        E_abs = nStages / (nStages + 1);
      } else {
        const AN1 = Math.pow(A_i, nStages + 1);
        E_abs = (AN1 - A_i) / (AN1 - 1);
      }
      E_abs = Math.max(0, Math.min(1, E_abs));

      const absorbed = y_in * E_abs;
      gasOut_moles[i] = y_in - absorbed;
      liqOut_moles[i] = x_in + absorbed;
    }

    const gasTotal = gasOut_moles.reduce((s, v) => s + v, 0);
    const liqTotal = liqOut_moles.reduce((s, v) => s + v, 0);
    return {
      gasTotal,
      liqTotal,
      yOut: gasOut_moles.map(m => gasTotal > 0 ? m / gasTotal : 0),
      xOut: liqOut_moles.map(m => liqTotal > 0 ? m / liqTotal : 0),
    };
  };

  let split = computeSplit(K);
  for (let iter = 0; iter < 25; iter++) {
    const nextX = [...split.xOut];
    const nextY = [...split.yOut];
    const maxDelta = Math.max(
      ...nextX.map((xi, i) => Math.abs(xi - x_est[i])),
      ...nextY.map((yi, i) => Math.abs(yi - y_est[i]))
    );
    x_est = nextX;
    y_est = nextY;
    if (maxDelta < 1e-8) break;
    K = thermo.fluidPackage === 'Ideal'
      ? KvaluesIdeal(compounds, T_avg, P_Pa)
      : computeKvalues(compounds, x_est, y_est, T_avg, P_Pa, thermo.fluidPackage, thermo.interactionParams);
    split = computeSplit(K);
  }

  const { gasTotal, liqTotal, yOut, xOut } = split;

  const flashGas = flashPT(compounds, yOut, T_avg, P_Pa, thermo.fluidPackage, thermo.interactionParams);
  const flashLiq = flashPT(compounds, xOut, T_avg, P_Pa, thermo.fluidPackage, thermo.interactionParams);

  return {
    gasOutlet: {
      ...gasInlet,
      T_K: T_avg, P_Pa, totalFlow_molps: gasTotal,
      moleFractions: yOut,
      phase: flashGas.phase, vaporFraction: flashGas.vaporFraction,
      H_Jpmol: flashGas.H_Jpmol,
      x_liquid: flashGas.x, y_vapor: flashGas.y, solved: true,
    },
    liquidOutlet: {
      ...liquidInlet,
      T_K: T_avg, P_Pa, totalFlow_molps: liqTotal,
      moleFractions: xOut,
      phase: flashLiq.phase, vaporFraction: flashLiq.vaporFraction,
      H_Jpmol: flashLiq.H_Jpmol,
      x_liquid: flashLiq.x, y_vapor: flashLiq.y, solved: true,
    },
  };
}

// ────────────────────────────────────────────────────────────────
// Rigorous Distillation — Bubble-Point MESH Method
// Stage-by-stage tray-to-tray calculation with energy balance
// ────────────────────────────────────────────────────────────────

function compositionFromMoles(componentMoles: number[]): { total: number; composition: number[] } {
  const total = componentMoles.reduce((sum, value) => sum + value, 0);
  if (total <= 1e-30) {
    const n = componentMoles.length;
    return { total: 0, composition: new Array(n).fill(n > 0 ? 1 / n : 0) };
  }
  return {
    total,
    composition: componentMoles.map(value => Math.max(0, value) / total),
  };
}

function buildLiquidStreamFromMoles(
  compounds: CanopyCompound[],
  componentMoles: number[],
  T_K: number,
  P_Pa: number,
  thermo: ThermoOptions,
): MaterialStream {
  const { total, composition } = compositionFromMoles(componentMoles);
  const flash = flashPT(compounds, composition, T_K, P_Pa, thermo.fluidPackage, thermo.interactionParams);
  return {
    id: '',
    name: '',
    T_K,
    P_Pa,
    totalFlow_molps: total,
    moleFractions: composition,
    phase: 'Liquid',
    vaporFraction: 0,
    H_Jpmol: flash.H_Jpmol,
    x_liquid: composition,
    y_vapor: composition,
    sourceUnitId: null,
    sourcePortId: null,
    targetUnitId: null,
    targetPortId: null,
    solved: true,
  };
}

function buildStreamFromMoles(
  compounds: CanopyCompound[],
  componentMoles: number[],
  T_K: number,
  P_Pa: number,
  thermo: ThermoOptions,
  phaseHint?: Phase,
): MaterialStream {
  const { total, composition } = compositionFromMoles(componentMoles);
  const flash = flashPT(compounds, composition, T_K, P_Pa, thermo.fluidPackage, thermo.interactionParams);
  const phase = phaseHint ?? flash.phase;
  const vaporFraction = phaseHint === 'Liquid'
    ? 0
    : phaseHint === 'Vapor'
      ? 1
      : flash.vaporFraction;
  return {
    id: '',
    name: '',
    T_K,
    P_Pa,
    totalFlow_molps: total,
    moleFractions: composition,
    phase,
    vaporFraction,
    H_Jpmol: flash.H_Jpmol,
    x_liquid: flash.x,
    y_vapor: flash.y,
    sourceUnitId: null,
    sourcePortId: null,
    targetUnitId: null,
    targetPortId: null,
    solved: true,
  };
}

function maxVectorDelta(a: number[][], b: number[][]): number {
  let maxDelta = 0;
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < a[i].length; j++) {
      maxDelta = Math.max(maxDelta, Math.abs(a[i][j] - b[i][j]));
    }
  }
  return maxDelta;
}

function contactExtractionStage(
  compounds: CanopyCompound[],
  raffinateInMoles: number[],
  extractInMoles: number[],
  T_K: number,
  P_Pa: number,
  solventKeyIndex: number,
  thermo: ThermoOptions,
): { raffinateOut: number[]; extractOut: number[] } {
  const n = compounds.length;
  const totalMoles = raffinateInMoles.map((value, i) => value + extractInMoles[i]);
  const mixed = compositionFromMoles(totalMoles);
  if (mixed.total <= 1e-30) {
    return {
      raffinateOut: [...raffinateInMoles],
      extractOut: [...extractInMoles],
    };
  }

  const vlle = flashPT_VLLE(compounds, mixed.composition, T_K, P_Pa, thermo.fluidPackage, thermo.interactionParams);
  const totalLiquidFraction = Math.max(vlle.liquidIFraction + vlle.liquidIIFraction, 1e-30);
  if (vlle.liquidIIFraction > 1e-8 && totalLiquidFraction > 1e-12) {
    const liquidTotal = mixed.total * Math.max(1 - vlle.vaporFraction, 0);
    const phaseITotal = liquidTotal * vlle.liquidIFraction / totalLiquidFraction;
    const phaseIITotal = liquidTotal * vlle.liquidIIFraction / totalLiquidFraction;
    const phaseIMoles = vlle.xI.map(x => x * phaseITotal);
    const phaseIIMoles = vlle.xII.map(x => x * phaseIITotal);
    const phaseIIsExtract = vlle.xII[solventKeyIndex] >= vlle.xI[solventKeyIndex];

    return phaseIIsExtract
      ? { raffinateOut: phaseIMoles, extractOut: phaseIIMoles }
      : { raffinateOut: phaseIIMoles, extractOut: phaseIMoles };
  }

  if (thermo.fluidPackage === 'Ideal') {
    return {
      raffinateOut: [...raffinateInMoles],
      extractOut: [...extractInMoles],
    };
  }

  const raffinateComp = compositionFromMoles(raffinateInMoles).composition;
  const extractComp = compositionFromMoles(extractInMoles).composition;
  const gammaR = computeGamma(compounds, raffinateComp, T_K, thermo.fluidPackage, thermo.interactionParams);
  const gammaE = computeGamma(compounds, extractComp, T_K, thermo.fluidPackage, thermo.interactionParams);
  const raffinateOut = new Array(n).fill(0);
  const extractOut = new Array(n).fill(0);

  for (let i = 0; i < n; i++) {
    const distribution = Math.max(1e-6, gammaR[i] / Math.max(gammaE[i], 1e-6));
    raffinateOut[i] = totalMoles[i] / (1 + distribution);
    extractOut[i] = totalMoles[i] - raffinateOut[i];
  }

  return { raffinateOut, extractOut };
}

export function solveExtractor(
  compounds: CanopyCompound[],
  feedInlet: MaterialStream,
  solventInlet: MaterialStream,
  nStages: number,
  P_Pa: number,
  thermo: ThermoOptions = DEFAULT_THERMO,
  solventKeyIndex?: number,
): { extractOutlet: MaterialStream; raffinateOutlet: MaterialStream; converged: boolean } {
  const n = compounds.length;
  const stages = Math.max(1, Math.round(nStages));
  const T_avg = (feedInlet.T_K + solventInlet.T_K) / 2;
  const feedMoles = feedInlet.moleFractions.map(x => x * feedInlet.totalFlow_molps);
  const solventMoles = solventInlet.moleFractions.map(x => x * solventInlet.totalFlow_molps);
  const feedFlow = Math.max(feedInlet.totalFlow_molps, 1e-30);
  const solventFlow = Math.max(solventInlet.totalFlow_molps, 1e-30);
  const keyIndex = typeof solventKeyIndex === 'number' && solventKeyIndex >= 0 && solventKeyIndex < n
    ? solventKeyIndex
    : solventInlet.moleFractions.reduce((best, value, index, arr) => value > arr[best] ? index : best, 0);
  const mixedComp = compositionFromMoles(feedMoles.map((value, i) => value + solventMoles[i])).composition;
  const vlle = flashPT_VLLE(compounds, mixedComp, T_avg, P_Pa, thermo.fluidPackage, thermo.interactionParams);
  let distribution = new Array(n).fill(0);

  if (vlle.liquidIIFraction > 1e-8) {
    const phaseIIsExtract = vlle.xII[keyIndex] >= vlle.xI[keyIndex];
    const xExtract = phaseIIsExtract ? vlle.xII : vlle.xI;
    const xRaffinate = phaseIIsExtract ? vlle.xI : vlle.xII;
    distribution = xExtract.map((value, i) => value / Math.max(xRaffinate[i], 1e-9));
  } else if (thermo.fluidPackage !== 'Ideal') {
    const feedComp = compositionFromMoles(feedMoles).composition;
    const solventComp = compositionFromMoles(solventMoles).composition;
    const gammaR = computeGamma(compounds, feedComp, T_avg, thermo.fluidPackage, thermo.interactionParams);
    const gammaE = computeGamma(compounds, solventComp, T_avg, thermo.fluidPackage, thermo.interactionParams);
    distribution = gammaR.map((value, i) => Math.max(1e-6, value / Math.max(gammaE[i], 1e-6)));
  }

  const finalExtract = new Array(n).fill(0);
  const finalRaffinate = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const extractionFactor = (solventFlow * Math.max(distribution[i], 0)) / feedFlow;
    let fractionExtracted = 0;
    if (thermo.fluidPackage !== 'Ideal' && feedMoles[i] > 1e-30) {
      if (Math.abs(extractionFactor - 1) < 1e-6) {
        fractionExtracted = stages / (stages + 1);
      } else {
        const eN1 = Math.pow(extractionFactor, stages + 1);
        fractionExtracted = (eN1 - extractionFactor) / (eN1 - 1);
      }
      fractionExtracted = Math.max(0, Math.min(1, fractionExtracted));
    }
    const transferred = feedMoles[i] * fractionExtracted;
    finalExtract[i] = solventMoles[i] + transferred;
    finalRaffinate[i] = feedMoles[i] - transferred;
  }

  return {
    extractOutlet: buildLiquidStreamFromMoles(compounds, finalExtract, T_avg, P_Pa, thermo),
    raffinateOutlet: buildLiquidStreamFromMoles(compounds, finalRaffinate, T_avg, P_Pa, thermo),
    converged: true,
  };
}

export interface RigorousColumnSpec {
  nStages: number;
  feedStage: number;          // 1-based (1 = condenser top)
  refluxRatio: number;
  distillateRate_molps: number;
  condenserP_Pa: number;
  reboilerP_Pa: number;
  condenserType: 'Total' | 'Partial';
  murphreeEfficiency?: number;
  liquidSideDrawStage?: number;
  liquidSideDrawFraction?: number;
  vaporSideDrawStage?: number;
  vaporSideDrawFraction?: number;
  spec1Mode?: RigorousColumnSpecMode;
  spec1Value?: number;
  spec1ComponentIndex?: number;
  spec2Mode?: RigorousColumnSpecMode;
  spec2Value?: number;
  spec2ComponentIndex?: number;
}

export type RigorousColumnSpecMode =
  | 'None'
  | 'DistillateRate'
  | 'DistillateFraction'
  | 'BottomsRate'
  | 'RefluxRatio'
  | 'BoilupRatio'
  | 'DistillateComposition';

export interface RigorousColumnSpecDiagnostics {
  status: 'ok' | 'underspecified' | 'overspecified' | 'invalid';
  activeSpecCount: number;
  messages: string[];
}

export interface RigorousColumnResult {
  converged: boolean;
  distillate: Partial<MaterialStream>;
  bottoms: Partial<MaterialStream>;
  liquidSideDraw?: Partial<MaterialStream>;
  vaporSideDraw?: Partial<MaterialStream>;
  duty_condenser_W: number;
  duty_reboiler_W: number;
  stageDuties_W: number[];
  stageTemperatures: number[];
  stageCompositions: number[][];  // x[stage][comp]
  stageVaporCompositions: number[][];
}

export interface RateFracSpec extends RigorousColumnSpec {
  columnDiameter_m?: number;
  traySpacing_m?: number;
  weirHeight_m?: number;
  interfacialArea_m2pm3?: number;
  kL_mps?: number;
  kV_mps?: number;
  liquidHoldup_s?: number;
  vaporHoldup_s?: number;
  internalsMode?: 'Tray' | 'Packing';
  pressureRelaxation?: number;
}

export interface RateFracResult extends RigorousColumnResult {
  effectiveMurphreeEfficiency: number;
  stageEfficiencies: number[];
  floodingFractions: number[];
  stagePressureDrops_Pa: number[];
  interfacialLiquidCompositions: number[][];
  interfacialVaporCompositions: number[][];
  liquidMassTransferCoefficients_mps: number[][];
  vaporMassTransferCoefficients_mps: number[][];
  componentMolarFluxes_molpsm2: number[][];
  stageHeatTransferRates_W: number[];
  couplingIterations: number;
  pressureCouplingResidual_Pa: number;
  hydraulicRegimes: string[];
}

const FLOW_SPEC_MODES = new Set<RigorousColumnSpecMode>(['DistillateRate', 'DistillateFraction', 'BottomsRate']);
const ENERGY_SPEC_MODES = new Set<RigorousColumnSpecMode>(['RefluxRatio', 'BoilupRatio']);
const COMPOSITION_SPEC_MODES = new Set<RigorousColumnSpecMode>(['DistillateComposition']);

interface ColumnSpecSelection {
  mode: RigorousColumnSpecMode;
  value: number;
  componentIndex?: number;
  slot: 'spec1' | 'spec2';
}

function getActiveRigorousColumnSpecs(spec: RigorousColumnSpec): ColumnSpecSelection[] {
  const selections: ColumnSpecSelection[] = [];
  if (spec.spec1Mode && spec.spec1Mode !== 'None' && Number.isFinite(spec.spec1Value)) {
    selections.push({
      mode: spec.spec1Mode,
      value: spec.spec1Value as number,
      componentIndex: spec.spec1ComponentIndex,
      slot: 'spec1',
    });
  }
  if (spec.spec2Mode && spec.spec2Mode !== 'None' && Number.isFinite(spec.spec2Value)) {
    selections.push({
      mode: spec.spec2Mode,
      value: spec.spec2Value as number,
      componentIndex: spec.spec2ComponentIndex,
      slot: 'spec2',
    });
  }
  return selections;
}

export function getRigorousColumnSpecDiagnostics(spec: RigorousColumnSpec): RigorousColumnSpecDiagnostics {
  const activeSpecs = getActiveRigorousColumnSpecs(spec);
  if (activeSpecs.length === 0) {
    return {
      status: 'ok',
      activeSpecCount: 0,
      messages: ['Using legacy manual inputs for reflux ratio and distillate rate.'],
    };
  }
  if (activeSpecs.length < 2) {
    return {
      status: 'underspecified',
      activeSpecCount: activeSpecs.length,
      messages: ['Rigorous columns need two active operating specifications.'],
    };
  }
  const categories = activeSpecs.map(selection =>
    FLOW_SPEC_MODES.has(selection.mode)
      ? 'flow'
      : ENERGY_SPEC_MODES.has(selection.mode)
        ? 'energy'
        : COMPOSITION_SPEC_MODES.has(selection.mode)
          ? 'composition'
          : 'other',
  );
  const flowCount = categories.filter(category => category === 'flow').length;
  const energyCount = categories.filter(category => category === 'energy').length;
  const compositionCount = categories.filter(category => category === 'composition').length;
  if (flowCount > 1 || energyCount > 1 || compositionCount > 1) {
    return {
      status: 'overspecified',
      activeSpecCount: activeSpecs.length,
      messages: ['Selected column specs overconstrain the same degree of freedom.'],
    };
  }
  if (flowCount === 0 && compositionCount === 0) {
    return {
      status: 'underspecified',
      activeSpecCount: activeSpecs.length,
      messages: ['At least one flow or composition spec is needed to determine product split.'],
    };
  }
  if (energyCount === 0 && compositionCount === 0) {
    return {
      status: 'underspecified',
      activeSpecCount: activeSpecs.length,
      messages: ['At least one reflux or boilup-style spec is needed to determine column internal traffic.'],
    };
  }
  if (flowCount === 1 && energyCount === 1) {
    return {
      status: 'ok',
      activeSpecCount: activeSpecs.length,
      messages: ['Flow and energy specs are sufficient for a rigorous column solve.'],
    };
  }
  if (compositionCount === 1 && (flowCount === 1 || energyCount === 1)) {
    return {
      status: 'ok',
      activeSpecCount: activeSpecs.length,
      messages: ['Composition spec will be paired with one operating spec to solve the remaining degree of freedom.'],
    };
  }
  return {
    status: 'invalid',
    activeSpecCount: activeSpecs.length,
    messages: ['Column spec combination is not currently supported.'],
  };
}

function evaluateRootCandidate(
  lower: number,
  upper: number,
  evaluator: (value: number) => number,
): { value: number; converged: boolean } {
  const samples = 8;
  let bestValue = lower;
  let bestError = Number.POSITIVE_INFINITY;
  let previousValue = lower;
  let previousError = evaluator(lower);
  bestError = Math.abs(previousError);
  for (let i = 1; i <= samples; i++) {
    const value = lower + (upper - lower) * i / samples;
    const error = evaluator(value);
    if (Math.abs(error) < bestError) {
      bestError = Math.abs(error);
      bestValue = value;
    }
    if (previousError === 0 || error === 0 || previousError * error < 0) {
      let lo = previousValue;
      let hi = value;
      let fLo = previousError;
      let fHi = error;
      for (let iter = 0; iter < 18; iter++) {
        const mid = 0.5 * (lo + hi);
        const fMid = evaluator(mid);
        if (Math.abs(fMid) < 1e-4) return { value: mid, converged: true };
        if (fLo * fMid <= 0) {
          hi = mid;
          fHi = fMid;
        } else {
          lo = mid;
          fLo = fMid;
        }
        if (Math.abs(hi - lo) < 1e-5 * Math.max(1, Math.abs(mid))) {
          return { value: 0.5 * (lo + hi), converged: true };
        }
        void fHi;
      }
      return { value: 0.5 * (lo + hi), converged: true };
    }
    previousValue = value;
    previousError = error;
  }
  return { value: bestValue, converged: false };
}

/**
 * Rigorous distillation column using the Bubble-Point (BP) method.
 *
 * MESH equations solved iteratively:
 *   M: Material balance on each stage
 *   E: Equilibrium (K-values from bubble-point)
 *   S: Summation (Σ x_i = 1 on each stage)
 *   H: Heat balance (energy balance for stage temperatures)
 *
 * Algorithm (Wang-Henke Bubble-Point method):
 * 1. Assume linear temperature profile
 * 2. Calculate K-values at each stage temperature
 * 3. Solve tridiagonal system for component liquid flows
 * 4. Normalize compositions, update T via bubble-point
 * 5. Repeat until converged
 *
 * Reference: Wang & Henke, Hydrocarbon Processing 45(8):155 (1966).
 */
function solveColumnRigorousCore(
  compounds: CanopyCompound[],
  feed: MaterialStream,
  spec: RigorousColumnSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): RigorousColumnResult {
  const nc = compounds.length;
  const ns = spec.nStages;
  const fs = spec.feedStage - 1; // 0-based feed stage

  const F = feed.totalFlow_molps;
  const z = [...feed.moleFractions];
  const D = spec.distillateRate_molps;
  const RR = spec.refluxRatio;
  const L0 = RR * D; // reflux from condenser
  const q = feed.vaporFraction < 0.01 ? 1 : (1 - feed.vaporFraction);
  const murphreeEfficiency = Math.max(0.01, Math.min(1, spec.murphreeEfficiency ?? 1));
  const liquidSideStage = spec.liquidSideDrawStage != null ? Math.max(0, Math.min(ns - 1, spec.liquidSideDrawStage - 1)) : null;
  const vaporSideStage = spec.vaporSideDrawStage != null ? Math.max(0, Math.min(ns - 1, spec.vaporSideDrawStage - 1)) : null;
  const liquidSideFraction = Math.max(0, Math.min(0.5, spec.liquidSideDrawFraction ?? 0));
  const vaporSideFraction = Math.max(0, Math.min(0.5, spec.vaporSideDrawFraction ?? 0));

  // Pressure profile (linear)
  const P = new Array(ns).fill(0);
  for (let j = 0; j < ns; j++) {
    P[j] = spec.condenserP_Pa + (spec.reboilerP_Pa - spec.condenserP_Pa) * j / Math.max(ns - 1, 1);
  }

  // Initial temperature profile (linear between estimated bubble points)
  // Aspen approach: bubble point at top, dew point at bottom, linear interpolation
  const T_top = bubblePointT_Ideal(compounds, z, spec.condenserP_Pa) ?? feed.T_K - 20;
  const T_bot = dewPointT_Ideal(compounds, z, spec.reboilerP_Pa)
    ?? bubblePointT_Ideal(compounds, z, spec.reboilerP_Pa)
    ?? feed.T_K + 20;
  const T = new Array(ns).fill(0);
  for (let j = 0; j < ns; j++) {
    T[j] = T_top + (T_bot - T_top) * j / Math.max(ns - 1, 1);
  }

  // Initial compositions using Wilson K-value estimates (Aspen Plus Section 1.4)
  // At each stage, estimate K-values from Wilson correlation, then
  // normalize x and y from isothermal flash assumption
  const x: number[][] = Array.from({ length: ns }, () => [...z]);
  const y: number[][] = Array.from({ length: ns }, () => [...z]);
  for (let j = 0; j < ns; j++) {
    const Kw = KvaluesWilsonEstimate(compounds, T[j], P[j]);
    // Simple Rachford-Rice: assume feed-like composition, compute y = K*x
    // then normalize
    let sumKz = 0;
    for (let i = 0; i < nc; i++) sumKz += Kw[i] * z[i];
    for (let i = 0; i < nc; i++) {
      y[j][i] = Kw[i] * z[i] / (sumKz > 0 ? sumKz : 1);
      x[j][i] = z[i]; // start with feed composition for liquid
    }
  }

  // Feed component flows
  const fj = z.map(zi => zi * F);

  // ── Liquid and vapor flows (constant molar overflow) ──
  // Convention: L[j] = total liquid leaving stage j going down
  //             V[j] = total vapor leaving stage j going up
  // Stage 0 = condenser (total), Stage ns-1 = reboiler
  const Lbase = new Array(ns).fill(0);
  const Vbase = new Array(ns).fill(0);

  // Rectifying section (above feed): L = L0, V = L0 + D
  // Stripping section (feed & below): L = L0 + q*F, V = L0 + D - (1-q)*F
  Vbase[0] = 0; // total condenser: no vapor product
  Lbase[0] = L0;
  for (let j = 1; j < ns; j++) {
    if (j <= fs) {
      // Rectifying section (between condenser and feed)
      Lbase[j] = L0;
      Vbase[j] = L0 + D;
    } else {
      // Stripping section (feed stage and below)
      Lbase[j] = L0 + q * F;
      Vbase[j] = L0 + D - (1 - q) * F;
    }
  }
  const U = new Array(ns).fill(0);
  const W = new Array(ns).fill(0);
  if (liquidSideStage != null && liquidSideStage > 0 && liquidSideStage < ns - 1) {
    U[liquidSideStage] = liquidSideFraction * Math.max(Lbase[liquidSideStage], 0);
  }
  if (vaporSideStage != null && vaporSideStage > 0 && vaporSideStage < ns - 1) {
    W[vaporSideStage] = vaporSideFraction * Math.max(Vbase[vaporSideStage], 0);
  }

  const totalLiquidSideFlow = U.reduce((sum, value) => sum + value, 0);
  const totalVaporSideFlow = W.reduce((sum, value) => sum + value, 0);
  const B = Math.max(F - D - totalLiquidSideFlow - totalVaporSideFlow, 1e-9);

  const L = [...Lbase];
  const V = [...Vbase];
  for (let j = 0; j < ns; j++) {
    const liquidDrawAbove = U.slice(0, j + 1).reduce((sum, value) => sum + value, 0);
    L[j] = Math.max(Lbase[j] - liquidDrawAbove, 1e-9);
  }
  for (let j = 1; j < ns; j++) {
    const vaporDrawBelow = W.slice(j, ns).reduce((sum, value) => sum + value, 0);
    V[j] = Math.max(Vbase[j] - vaporDrawBelow, 1e-9);
  }
  L[ns - 1] = B; // bottoms leaves reboiler after side draws

  let converged = false;

  for (let iter = 0; iter < 200; iter++) {
    // Step 1: K-values at each stage temperature (respects fluid package)
    const K: number[][] = Array.from({ length: ns }, (_, j) =>
      thermo.fluidPackage === 'Ideal'
        ? KvaluesIdeal(compounds, T[j], P[j])
        : computeKvalues(compounds, x[j], y[j], T[j], P[j], thermo.fluidPackage, thermo.interactionParams)
    );
    const Keff: number[][] = K.map(row =>
      row.map(k => Math.max(1e-6, 1 + murphreeEfficiency * (k - 1)))
    );

    // Step 2: Solve tridiagonal system for each component
    // Variable: l_{j,i} = L_j * x_{j,i}  (component liquid molar flow)
    //
    // Stage j material balance (1 ≤ j ≤ ns-2):
    //   l_{j-1,i} + V_{j+1}·K_{j+1,i}·l_{j+1,i}/L_{j+1} + F_j·z_i
    //       = l_{j,i} + V_j·K_{j,i}·l_{j,i}/L_j
    //
    // Tridiagonal:  a_j·l_{j-1} + b_j·l_j + c_j·l_{j+1} = d_j
    //   a_j = 1
    //   b_j = -(1 + V_j·K_{j,i} / L_j)
    //   c_j = V_{j+1}·K_{j+1,i} / L_{j+1}
    //   d_j = -F_j·z_i
    //
    // Condenser (j=0, total):
    //   V_1·K_{1,i}·l_{1,i}/L_1 = (L_0+D)/L_0 · l_{0,i}
    //   b_0 = -(1 + D/L_0),  c_0 = V_1·K_{1,i}/L_1,  d_0 = 0
    //
    // Reboiler (j=ns-1):
    //   a = 1,  b = -(1 + V_{ns-1}·K_{ns-1,i}/L_{ns-1}),  c = 0

    for (let i = 0; i < nc; i++) {
      const a = new Array(ns).fill(0);
      const b = new Array(ns).fill(0);
      const c = new Array(ns).fill(0);
      const d = new Array(ns).fill(0);

      // Condenser (j=0)
      a[0] = 0;
        b[0] = -(1 + D / Math.max(L[0], 1e-9) + U[0] / Math.max(L[0], 1e-9));
        c[0] = ns > 1 ? V[1] * Keff[1][i] / L[1] : 0;
        d[0] = 0;

      // Interior stages (j=1 to ns-2)
        for (let j = 1; j < ns - 1; j++) {
          a[j] = 1;
          b[j] = -(1 + U[j] / Math.max(L[j], 1e-9) + (V[j] + W[j]) * Keff[j][i] / Math.max(L[j], 1e-9));
          c[j] = V[j + 1] * Keff[j + 1][i] / L[j + 1];
          d[j] = j === fs ? -fj[i] : 0;
        }

      // Reboiler (j=ns-1)
      a[ns - 1] = 1;
      b[ns - 1] = -(1 + U[ns - 1] / Math.max(L[ns - 1], 1e-9) + (V[ns - 1] + W[ns - 1]) * Keff[ns - 1][i] / Math.max(L[ns - 1], 1e-9));
      c[ns - 1] = 0;
      d[ns - 1] = fs === ns - 1 ? -fj[i] : 0;

      // Thomas algorithm (tridiagonal solver)
      const cp = new Array(ns).fill(0);
      const dp = new Array(ns).fill(0);

      cp[0] = c[0] / b[0];
      dp[0] = d[0] / b[0];

      for (let j = 1; j < ns; j++) {
        const m = b[j] - a[j] * cp[j - 1];
        if (Math.abs(m) < 1e-30) { cp[j] = 0; dp[j] = 0; continue; }
        cp[j] = c[j] / m;
        dp[j] = (d[j] - a[j] * dp[j - 1]) / m;
      }

      // Back substitution
      const liq = new Array(ns).fill(0);
      liq[ns - 1] = dp[ns - 1];
      for (let j = ns - 2; j >= 0; j--) {
        liq[j] = dp[j] - cp[j] * liq[j + 1];
      }

      for (let j = 0; j < ns; j++) {
        x[j][i] = Math.max(liq[j], 1e-15);
      }
    }

    // Step 3: Normalize compositions and compute vapor.
    // Murphree vapor efficiency uses the incoming vapor from the stage above
    // as the reference instead of forcing full equilibrium on every tray.
    for (let j = 0; j < ns; j++) {
      const sumX = x[j].reduce((s, v) => s + v, 0) || 1;
      for (let i = 0; i < nc; i++) x[j][i] /= sumX;

      const yEq = new Array(nc).fill(0);
      let sumYEq = 0;
      for (let i = 0; i < nc; i++) {
        yEq[i] = Keff[j][i] * x[j][i];
        sumYEq += yEq[i];
      }
      if (sumYEq > 0) {
        for (let i = 0; i < nc; i++) yEq[i] /= sumYEq;
      }

      let sumY = 0;
      for (let i = 0; i < nc; i++) {
        const yIn = j > 0 ? y[j - 1][i] : yEq[i];
        y[j][i] = yIn + murphreeEfficiency * (yEq[i] - yIn);
        sumY += y[j][i];
      }
      if (sumY > 0) {
        for (let i = 0; i < nc; i++) y[j][i] /= sumY;
      }
    }

    // Step 4: Update temperatures via bubble-point (respects fluid package)
    let maxDeltaT = 0;
    for (let j = 0; j < ns; j++) {
      const T_new = thermo.fluidPackage === 'Ideal'
        ? bubblePointT_Ideal(compounds, x[j], P[j])
        : bubblePointT(compounds, x[j], P[j], thermo.fluidPackage, thermo.interactionParams);
      if (T_new !== null) {
        const dT = Math.abs(T_new - T[j]);
        if (dT > maxDeltaT) maxDeltaT = dT;
        T[j] = T_new;
      }
    }

    if (maxDeltaT < 0.01) { converged = true; break; }
  }

  // Build outputs
  const distComps = spec.condenserType === 'Total' ? [...x[0]] : [...y[0]];
  const botComps = [...x[ns - 1]];
  const liquidSideFlow = totalLiquidSideFlow;
  const vaporSideFlow = totalVaporSideFlow;

  const flashD = flashPT(compounds, distComps, T[0], P[0], thermo.fluidPackage, thermo.interactionParams);
  const flashB = flashPT(compounds, botComps, T[ns - 1], P[ns - 1], thermo.fluidPackage, thermo.interactionParams);
  const liquidSideDraw = liquidSideStage != null && liquidSideFlow > 1e-12
    ? buildStreamFromMoles(
      compounds,
      x[liquidSideStage].map(value => value * liquidSideFlow),
      T[liquidSideStage],
      P[liquidSideStage],
      thermo,
      'Liquid',
    )
    : undefined;
  const vaporSideDraw = vaporSideStage != null && vaporSideFlow > 1e-12
    ? buildStreamFromMoles(
      compounds,
      y[vaporSideStage].map(value => value * vaporSideFlow),
      T[vaporSideStage],
      P[vaporSideStage],
      thermo,
      'Vapor',
    )
    : undefined;

  const stageLiquidEnthalpy = x.map((stageX, j) =>
    streamEnthalpy(compounds, T[j], 0, stageX, stageX, P[j], thermo.fluidPackage, undefined, thermo.interactionParams),
  );
  const stageVaporEnthalpy = y.map((stageY, j) =>
    streamEnthalpy(compounds, T[j], 1, stageY, stageY, P[j], thermo.fluidPackage, undefined, thermo.interactionParams),
  );
  const stageDuties_W = new Array(ns).fill(0);
  for (let j = 0; j < ns; j++) {
    const liquidIn = j > 0 ? L[j - 1] * stageLiquidEnthalpy[j - 1] : 0;
    const vaporIn = j < ns - 1 ? V[j + 1] * stageVaporEnthalpy[j + 1] : 0;
    const feedIn = j === fs ? F * feed.H_Jpmol : 0;
    const liquidProduct =
      (j === 0 && spec.condenserType === 'Total' ? D : 0)
      + (j === ns - 1 ? B : 0)
      + (liquidSideStage === j ? liquidSideFlow : 0);
    const vaporProduct =
      (j === 0 && spec.condenserType === 'Partial' ? D : 0)
      + (vaporSideStage === j ? vaporSideFlow : 0);
    const liquidOut = L[j] * stageLiquidEnthalpy[j] + liquidProduct * stageLiquidEnthalpy[j];
    const vaporOut = Math.max(V[j], 0) * stageVaporEnthalpy[j] + vaporProduct * stageVaporEnthalpy[j];
    stageDuties_W[j] = liquidIn + vaporIn + feedIn - liquidOut - vaporOut;
  }

  // Energy balance
  const H_D = D * streamEnthalpy(compounds, T[0], 0, distComps, distComps);
  const H_B = B * streamEnthalpy(compounds, T[ns - 1], 0, botComps, botComps);
  const H_F = F * feed.H_Jpmol;
  const H_V1 = (L0 + D) * streamEnthalpy(compounds, T[1] ?? T[0], 1, y[1] ?? y[0], y[1] ?? y[0]);
  const H_L0 = L0 * streamEnthalpy(compounds, T[0], 0, x[0], x[0]);
  const Q_cond = H_L0 + H_D - H_V1;
  const Q_reb = H_D + H_B - H_F - Q_cond;

  return {
    converged,
    distillate: {
      T_K: T[0], P_Pa: P[0], totalFlow_molps: D,
      moleFractions: distComps,
      phase: spec.condenserType === 'Total' ? 'Liquid' : flashD.phase,
      vaporFraction: spec.condenserType === 'Total' ? 0 : flashD.vaporFraction,
      H_Jpmol: flashD.H_Jpmol,
      x_liquid: flashD.x, y_vapor: flashD.y, solved: true,
    },
    bottoms: {
      T_K: T[ns - 1], P_Pa: P[ns - 1], totalFlow_molps: B,
      moleFractions: botComps,
      phase: 'Liquid', vaporFraction: 0,
      H_Jpmol: flashB.H_Jpmol,
      x_liquid: flashB.x, y_vapor: flashB.y, solved: true,
    },
    liquidSideDraw,
    vaporSideDraw,
    duty_condenser_W: Q_cond,
    duty_reboiler_W: Q_reb,
    stageDuties_W,
    stageTemperatures: [...T],
    stageCompositions: x.map(row => [...row]),
    stageVaporCompositions: y.map(row => [...row]),
  };
}

function resolveRigorousColumnOperatingSpec(
  compounds: CanopyCompound[],
  feed: MaterialStream,
  spec: RigorousColumnSpec,
  thermo: ThermoOptions,
): {
  diagnostics: RigorousColumnSpecDiagnostics;
  resolvedDistillateRate_molps: number;
  resolvedRefluxRatio: number;
} {
  const diagnostics = getRigorousColumnSpecDiagnostics(spec);
  if (diagnostics.status !== 'ok') {
    return {
      diagnostics,
      resolvedDistillateRate_molps: spec.distillateRate_molps,
      resolvedRefluxRatio: spec.refluxRatio,
    };
  }

  const activeSpecs = getActiveRigorousColumnSpecs(spec);
  if (activeSpecs.length === 0) {
    return {
      diagnostics,
      resolvedDistillateRate_molps: spec.distillateRate_molps,
      resolvedRefluxRatio: spec.refluxRatio,
    };
  }

  const feedFlow = Math.max(feed.totalFlow_molps, 1e-6);
  const q = feed.vaporFraction < 0.01 ? 1 : (1 - feed.vaporFraction);
  const flowSpec = activeSpecs.find(selection => FLOW_SPEC_MODES.has(selection.mode));
  const energySpec = activeSpecs.find(selection => ENERGY_SPEC_MODES.has(selection.mode));
  const compositionSpec = activeSpecs.find(selection => COMPOSITION_SPEC_MODES.has(selection.mode));

  const flowToDistillate = (selection: ColumnSpecSelection): number => {
    switch (selection.mode) {
      case 'DistillateRate':
        return selection.value;
      case 'DistillateFraction':
        return selection.value * feedFlow;
      case 'BottomsRate':
        return feedFlow - selection.value;
      default:
        return spec.distillateRate_molps;
    }
  };

  const boilupToReflux = (boilupRatio: number, distillateRate_molps: number): number => {
    const bottomsRate = Math.max(feedFlow - distillateRate_molps, 1e-6);
    const vaporRate = Math.max(boilupRatio, 0.01) * bottomsRate;
    return Math.max((vaporRate - distillateRate_molps + (1 - q) * feedFlow) / Math.max(distillateRate_molps, 1e-6), 0.05);
  };

  let resolvedDistillateRate_molps = flowSpec
    ? flowToDistillate(flowSpec)
    : spec.distillateRate_molps;
  resolvedDistillateRate_molps = Math.max(1e-6, Math.min(feedFlow - 1e-6, resolvedDistillateRate_molps));

  let resolvedRefluxRatio = energySpec
    ? (energySpec.mode === 'BoilupRatio'
        ? boilupToReflux(energySpec.value, resolvedDistillateRate_molps)
        : energySpec.value)
    : spec.refluxRatio;
  resolvedRefluxRatio = Math.max(0.05, resolvedRefluxRatio);

  if (compositionSpec) {
    const componentIndex = compositionSpec.componentIndex ?? 0;
    const target = Math.max(1e-6, Math.min(0.999999, compositionSpec.value));
    const evaluateDistillateComposition = (trialD: number, trialRR: number) => {
      const trial = solveColumnRigorousCore(compounds, feed, {
        ...spec,
        distillateRate_molps: trialD,
        refluxRatio: trialRR,
        spec1Mode: 'None',
        spec2Mode: 'None',
      }, thermo);
      return (trial.distillate.moleFractions?.[componentIndex] ?? 0) - target;
    };

    if (flowSpec && !energySpec) {
      const solved = evaluateRootCandidate(0.05, 20, value =>
        evaluateDistillateComposition(resolvedDistillateRate_molps, value),
      );
      resolvedRefluxRatio = solved.value;
    } else if (!flowSpec && energySpec) {
      const solved = evaluateRootCandidate(1e-6, feedFlow - 1e-6, value =>
        evaluateDistillateComposition(value, resolvedRefluxRatio),
      );
      resolvedDistillateRate_molps = solved.value;
    }
  }

  return {
    diagnostics,
    resolvedDistillateRate_molps,
    resolvedRefluxRatio,
  };
}

export function solveColumnRigorous(
  compounds: CanopyCompound[],
  feed: MaterialStream,
  spec: RigorousColumnSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): RigorousColumnResult {
  const resolved = resolveRigorousColumnOperatingSpec(compounds, feed, spec, thermo);
  if (resolved.diagnostics.status !== 'ok') {
    return {
      converged: false,
      distillate: {},
      bottoms: {},
      duty_condenser_W: 0,
      duty_reboiler_W: 0,
      stageDuties_W: [],
      stageTemperatures: [],
      stageCompositions: [],
      stageVaporCompositions: [],
    };
  }

  return solveColumnRigorousCore(compounds, feed, {
    ...spec,
    distillateRate_molps: resolved.resolvedDistillateRate_molps,
    refluxRatio: resolved.resolvedRefluxRatio,
    spec1Mode: 'None',
    spec2Mode: 'None',
  }, thermo);
}

export function solveRateFrac(
  compounds: CanopyCompound[],
  feed: MaterialStream,
  spec: RateFracSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): RateFracResult {
  const ns = Math.max(spec.nStages, 2);
  const diameter = Math.max(spec.columnDiameter_m ?? 1.5, 0.1);
  const traySpacing = Math.max(spec.traySpacing_m ?? 0.6, 0.05);
  const weirHeight = Math.max(spec.weirHeight_m ?? 0.05, 0.005);
  const internalsMode = spec.internalsMode === 'Packing' ? 'Packing' : 'Tray';
  const pressureRelaxation = Math.max(0.1, Math.min(1, spec.pressureRelaxation ?? 0.6));
  const area = Math.PI * diameter * diameter / 4;
  const kL = Math.max(spec.kL_mps ?? 2e-4, 1e-7);
  const kV = Math.max(spec.kV_mps ?? 5e-3, 1e-6);
  const aInt = Math.max(spec.interfacialArea_m2pm3 ?? 150, 1);
  const liquidResidence = Math.max(spec.liquidHoldup_s ?? 5, 0.1);
  const vaporResidence = Math.max(spec.vaporHoldup_s ?? 0.5, 0.01);
  const Ttop = bubblePointT_Ideal(compounds, feed.moleFractions, spec.condenserP_Pa) ?? Math.max(feed.T_K - 20, 250);
  const Tbot = dewPointT_Ideal(compounds, feed.moleFractions, spec.reboilerP_Pa)
    ?? bubblePointT_Ideal(compounds, feed.moleFractions, spec.reboilerP_Pa)
    ?? (feed.T_K + 20);
  const D = Math.max(spec.distillateRate_molps, 1e-6);
  const RR = Math.max(spec.refluxRatio, 0.01);
  const L0 = RR * D;
  const q = feed.vaporFraction < 0.01 ? 1 : (1 - feed.vaporFraction);
  const fs = Math.max(0, Math.min(ns - 1, spec.feedStage - 1));
  const initialTprofile = new Array(ns).fill(0).map((_, j) => Ttop + (Tbot - Ttop) * j / Math.max(ns - 1, 1));
  const stageInterfacialArea_m2 = Math.max(area * traySpacing * aInt, 1e-6);
  const maxOuterIterations = 4;

  let temperatureProfile = [...initialTprofile];
  let pressureProfile = new Array(ns).fill(0).map((_, j) =>
    spec.condenserP_Pa + (spec.reboilerP_Pa - spec.condenserP_Pa) * j / Math.max(ns - 1, 1),
  );
  let stageEfficiencies = new Array(ns).fill(Math.max(0.1, Math.min(0.95, spec.murphreeEfficiency ?? 0.65)));
  let floodingFractions = new Array(ns).fill(0);
  let stagePressureDrops_Pa = new Array(ns).fill(0);
  let liquidMassTransferCoefficients_mps = new Array(ns).fill(null).map(() => compounds.map(() => kL));
  let vaporMassTransferCoefficients_mps = new Array(ns).fill(null).map(() => compounds.map(() => kV));
  let interfacialLiquidCompositions = new Array(ns).fill(null).map(() => compounds.map(() => 1 / Math.max(compounds.length, 1)));
  let interfacialVaporCompositions = new Array(ns).fill(null).map(() => compounds.map(() => 1 / Math.max(compounds.length, 1)));
  let componentMolarFluxes_molpsm2 = new Array(ns).fill(null).map(() => compounds.map(() => 0));
  let stageHeatTransferRates_W = new Array(ns).fill(0);
  let hydraulicRegimes = new Array(ns).fill('normal');
  let pressureCouplingResidual_Pa = 0;
  let rigorous: RigorousColumnResult = {
    converged: false,
    distillate: {},
    bottoms: {},
    duty_condenser_W: 0,
    duty_reboiler_W: 0,
    stageDuties_W: [],
    stageTemperatures: [],
    stageCompositions: [],
    stageVaporCompositions: [],
  };
  let couplingIterations = 0;

  for (let outer = 0; outer < maxOuterIterations; outer++) {
    couplingIterations = outer + 1;
    const baseEfficiencies: number[] = [];
    floodingFractions = [];
    stagePressureDrops_Pa = [];
    liquidMassTransferCoefficients_mps = [];
    vaporMassTransferCoefficients_mps = [];
    hydraulicRegimes = [];

    for (let j = 0; j < ns; j++) {
      const L = j === 0 ? L0 : (j <= fs ? L0 : L0 + q * feed.totalFlow_molps);
      const V = j === 0 ? Math.max(L0 + D, 1e-6) : (j <= fs ? L0 + D : L0 + D - (1 - q) * feed.totalFlow_molps);
      const z = rigorous.stageCompositions[j] ?? feed.moleFractions;
      const Tstage = rigorous.stageTemperatures[j] ?? temperatureProfile[j] ?? feed.T_K;
      const Pstage = pressureProfile[j] ?? spec.condenserP_Pa;
      const MW = compounds.reduce((sum, compound, i) => sum + (z[i] ?? 0) * compound.molecularWeight, 0) / 1000;
      const rhoV = Math.max(Pstage * Math.max(MW, 1e-6) / (R * Math.max(Tstage, 150)), 0.05);
      const Vliq = Math.max(mixtureLiquidMolarVolume(compounds, z, Tstage), 1e-5);
      const rhoL = Math.max(MW / Vliq, rhoV + 1);
      const vaporVolFlow = Math.max(V, 1e-6) * R * Math.max(Tstage, 150) / Math.max(Pstage, 1000);
      const liquidVolFlow = Math.max(L, 1e-6) * Vliq;
      const uV = vaporVolFlow / area;
      const uL = liquidVolFlow / area;
      const cSB = internalsMode === 'Packing'
        ? 0.09 * Math.max(0.7, Math.pow(traySpacing / 0.6, 0.35)) * Math.pow(aInt / 150, 0.1)
        : 0.12 * Math.max(0.75, Math.sqrt(traySpacing / 0.6));
      const floodVelocity = Math.max(cSB * Math.sqrt(Math.max((rhoL - rhoV) / rhoV, 0.1)), 1e-4);
      const floodFraction = Math.max(0, Math.min(1.5, uV / floodVelocity));
      const hydraulicsBoost = internalsMode === 'Packing'
        ? Math.max(0.6, Math.min(1.8, Math.pow(aInt / 150, 0.35) * Math.pow(Math.max(uL, 1e-6) / 0.01, 0.08)))
        : Math.max(0.7, Math.min(1.4, Math.pow(Math.max(uL, 1e-6) / 0.01, 0.04)));
      const effectiveKL = kL * hydraulicsBoost;
      const effectiveKV = kV * (internalsMode === 'Packing' ? 1.15 : 1);
      const kOverall = (effectiveKL * liquidResidence + effectiveKV * vaporResidence) * aInt;
      const massTransferEfficiency = 1 - Math.exp(-Math.max(kOverall, 0));
      const hydraulicPenalty = floodFraction <= 0.75
        ? 1
        : Math.max(0.2, 1 - 1.6 * (floodFraction - 0.75));
      const baseEfficiency = Math.max(0.1, Math.min(0.95, massTransferEfficiency * hydraulicPenalty));
      const dryDp = internalsMode === 'Packing'
        ? Math.max(40, 95 * traySpacing * rhoV * uV * uV / Math.max(diameter, 0.15))
        : 250 + 0.5 * rhoV * uV * uV;
      const liquidHeadDp = internalsMode === 'Packing'
        ? rhoL * 9.80665 * Math.max(0.15 * weirHeight, 0.002)
        : rhoL * 9.80665 * weirHeight;
      const hydraulicRegime = floodFraction < 0.45
        ? (internalsMode === 'Packing' ? 'wetted' : 'spray')
        : floodFraction < 0.85
          ? (internalsMode === 'Packing' ? 'loading' : 'froth')
          : 'flooding';

      floodingFractions.push(floodFraction);
      stagePressureDrops_Pa.push(dryDp + liquidHeadDp + 150 * uL);
      baseEfficiencies.push(baseEfficiency);
      hydraulicRegimes.push(hydraulicRegime);
      liquidMassTransferCoefficients_mps.push(
        compounds.map(compound => effectiveKL * Math.pow(100 / Math.max(compound.molecularWeight, 1), 0.2)),
      );
      vaporMassTransferCoefficients_mps.push(
        compounds.map(compound => effectiveKV * Math.pow(100 / Math.max(compound.molecularWeight, 1), 0.1)),
      );
    }

    const effectiveMurphreeEfficiency = stageEfficiencies.reduce((sum, value) => sum + value, 0) / Math.max(stageEfficiencies.length, 1);
    const totalDp = stagePressureDrops_Pa.reduce((sum, value) => sum + value, 0);
    rigorous = solveColumnRigorous(compounds, feed, {
      ...spec,
      murphreeEfficiency: effectiveMurphreeEfficiency,
      condenserP_Pa: spec.condenserP_Pa,
      reboilerP_Pa: Math.max(spec.reboilerP_Pa, spec.condenserP_Pa + totalDp),
    }, thermo);

    temperatureProfile = rigorous.stageTemperatures.length === ns ? [...rigorous.stageTemperatures] : [...initialTprofile];
    const rawPressureProfile = stagePressureDrops_Pa.reduce<number[]>((profile, dp, index) => {
      if (index === 0) {
        profile.push(spec.condenserP_Pa);
      } else {
        profile.push(profile[index - 1] + stagePressureDrops_Pa[index - 1]);
      }
      return profile;
    }, []);
    const fallbackPressureProfile = new Array(ns).fill(0).map((_, j) => spec.condenserP_Pa + totalDp * j / Math.max(ns - 1, 1));
    const unconvergedPressureProfile = rawPressureProfile.length === ns ? rawPressureProfile : fallbackPressureProfile;
    pressureCouplingResidual_Pa = unconvergedPressureProfile.reduce((maxResidual, value, index) => {
      const previous = pressureProfile[index] ?? unconvergedPressureProfile[index];
      return Math.max(maxResidual, Math.abs(value - previous));
    }, 0);
    pressureProfile = unconvergedPressureProfile.map((value, index) => {
      const previous = pressureProfile[index] ?? value;
      return previous + pressureRelaxation * (value - previous);
    });
    if (pressureProfile.length < ns) {
      pressureProfile = [...fallbackPressureProfile];
    }

    interfacialLiquidCompositions = rigorous.stageCompositions.map((stageX, j) => {
      const kLStage = liquidMassTransferCoefficients_mps[j] ?? compounds.map(() => kL);
      const kVStage = vaporMassTransferCoefficients_mps[j] ?? compounds.map(() => kV);
      const stageY = rigorous.stageVaporCompositions[j] ?? stageX;
      const Pstage = pressureProfile[j] ?? spec.condenserP_Pa;
      const K = thermo.fluidPackage === 'Ideal'
        ? KvaluesIdeal(compounds, rigorous.stageTemperatures[j] ?? temperatureProfile[j], Pstage)
        : computeKvalues(
          compounds,
          stageX,
          stageY,
          rigorous.stageTemperatures[j] ?? temperatureProfile[j],
          Pstage,
          thermo.fluidPackage,
          thermo.interactionParams,
        );
      const weights = stageX.map((x, i) => {
        const backCalculatedXi = Math.max((stageY[i] ?? 0) / Math.max(K[i] ?? 1, 1e-12), 1e-12);
        return Math.max(kLStage[i] * x + kVStage[i] * backCalculatedXi, 1e-12);
      });
      const sum = weights.reduce((acc, value) => acc + value, 0) || 1;
      return weights.map(value => value / sum);
    });

    interfacialVaporCompositions = interfacialLiquidCompositions.map((stageXi, j) => {
      const Pstage = pressureProfile[j] ?? spec.condenserP_Pa;
      const K = thermo.fluidPackage === 'Ideal'
        ? KvaluesIdeal(compounds, rigorous.stageTemperatures[j] ?? temperatureProfile[j], Pstage)
        : computeKvalues(
          compounds,
          rigorous.stageCompositions[j] ?? stageXi,
          rigorous.stageVaporCompositions[j] ?? stageXi,
          rigorous.stageTemperatures[j] ?? temperatureProfile[j],
          Pstage,
          thermo.fluidPackage,
          thermo.interactionParams,
        );
      const yStar = stageXi.map((x, i) => Math.max(x * (K[i] ?? 1), 1e-12));
      const sum = yStar.reduce((acc, value) => acc + value, 0) || 1;
      return yStar.map(value => value / sum);
    });

    componentMolarFluxes_molpsm2 = rigorous.stageCompositions.map((stageX, j) => {
      const stageY = rigorous.stageVaporCompositions[j] ?? stageX;
      const stageXi = interfacialLiquidCompositions[j] ?? stageX;
      const stageYi = interfacialVaporCompositions[j] ?? stageY;
      const Tstage = rigorous.stageTemperatures[j] ?? temperatureProfile[j];
      const Pstage = pressureProfile[j] ?? spec.condenserP_Pa;
      const C_L = 1 / Math.max(mixtureLiquidMolarVolume(compounds, stageX, Tstage), 1e-6);
      const C_V = Math.max(Pstage / (R * Math.max(Tstage, 150)), 1e-6);
      return stageX.map((x, i) => {
        const liquidFlux = (liquidMassTransferCoefficients_mps[j]?.[i] ?? kL) * C_L * Math.abs(x - (stageXi[i] ?? x));
        const vaporFlux = (vaporMassTransferCoefficients_mps[j]?.[i] ?? kV) * C_V * Math.abs((stageYi[i] ?? stageY[i] ?? 0) - (stageY[i] ?? 0));
        return 0.5 * (liquidFlux + vaporFlux);
      });
    });

    stageHeatTransferRates_W = componentMolarFluxes_molpsm2.map((fluxes, j) => {
      const stageX = rigorous.stageCompositions[j] ?? feed.moleFractions;
      const stageY = rigorous.stageVaporCompositions[j] ?? feed.moleFractions;
      const Tstage = rigorous.stageTemperatures[j] ?? temperatureProfile[j];
      const Hv = streamEnthalpy(compounds, Tstage, 1, stageX, stageY);
      const Hl = streamEnthalpy(compounds, Tstage, 0, stageX, stageY);
      const deltaH = Math.abs(Hv - Hl);
      const totalFlux = fluxes.reduce((sum, value) => sum + value, 0);
      return totalFlux * stageInterfacialArea_m2 * deltaH;
    });

    const nextStageEfficiencies = baseEfficiencies.map((baseEfficiency, j) => {
      const stageInternalTraffic = Math.max(
        (j === 0 ? L0 : (j <= fs ? L0 : L0 + q * feed.totalFlow_molps))
        + (j === 0 ? L0 + D : (j <= fs ? L0 + D : L0 + D - (1 - q) * feed.totalFlow_molps)),
        1e-6,
      );
      const transferredMoles = componentMolarFluxes_molpsm2[j].reduce((sum, value) => sum + value, 0) * stageInterfacialArea_m2;
      const fluxUtilization = Math.tanh(3 * transferredMoles / stageInternalTraffic);
      const heatPenalty = 1 / (1 + stageHeatTransferRates_W[j] / Math.max(stageInternalTraffic * 5e4, 1));
      return Math.max(0.1, Math.min(0.98, baseEfficiency * (0.75 + 0.5 * fluxUtilization) * heatPenalty));
    });

    const maxEfficiencyChange = nextStageEfficiencies.reduce((acc, value, index) => Math.max(acc, Math.abs(value - (stageEfficiencies[index] ?? 0))), 0);
    stageEfficiencies = nextStageEfficiencies;
    if (maxEfficiencyChange < 1e-3 && pressureCouplingResidual_Pa < 25) break;
  }

  return {
    ...rigorous,
    effectiveMurphreeEfficiency: stageEfficiencies.reduce((sum, value) => sum + value, 0) / Math.max(stageEfficiencies.length, 1),
    stageEfficiencies,
    floodingFractions,
    stagePressureDrops_Pa,
    interfacialLiquidCompositions,
    interfacialVaporCompositions,
    liquidMassTransferCoefficients_mps,
    vaporMassTransferCoefficients_mps,
    componentMolarFluxes_molpsm2,
    stageHeatTransferRates_W,
    couplingIterations,
    pressureCouplingResidual_Pa,
    hydraulicRegimes,
  };
}

// ────────────────────────────────────────────────────────────────
// RYield (Yield Reactor)
// User specifies mass/mole yield vector. Normalize & flash product.
// ────────────────────────────────────────────────────────────────

export function solveRYield(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  yields: number[],               // mole yields per component (will be normalized)
  outletT_K: number | undefined,  // if undefined, use inlet T
  outletP_Pa: number | undefined, // if undefined, use inlet P
  thermo: ThermoOptions = DEFAULT_THERMO,
): { outlet: Partial<MaterialStream>; duty_W: number } {
  const nc = compounds.length;
  const T = outletT_K ?? inlet.T_K;
  const P = outletP_Pa ?? inlet.P_Pa;

  // Normalize yields
  let yieldSum = yields.reduce((a, b) => a + b, 0);
  if (yieldSum <= 0) yieldSum = 1;
  const normYields = yields.map(y => y / yieldSum);

  // Outlet flows: total feed flow redistributed by yields
  const F_total = inlet.totalFlow_molps;
  const flows = normYields.map(y => y * F_total);
  const totalOut = flows.reduce((a, b) => a + b, 0);
  const z_out = totalOut > 0 ? flows.map(f => f / totalOut) : new Array(nc).fill(1 / nc);

  const flash = flashPT(compounds, z_out, T, P, thermo.fluidPackage, thermo.interactionParams);
  const H_in = inlet.totalFlow_molps * inlet.H_Jpmol;
  const H_out = totalOut * flash.H_Jpmol;
  const duty = H_out - H_in;

  return {
    outlet: {
      T_K: T, P_Pa: P, totalFlow_molps: totalOut,
      moleFractions: z_out,
      phase: flash.phase, vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x, y_vapor: flash.y, solved: true,
    },
    duty_W: duty,
  };
}

// ────────────────────────────────────────────────────────────────
// Decanter (Liquid-Liquid Equilibrium separator)
// Performs LLE flash to split feed into 2 liquid phases.
// Uses activity coefficients to determine phase split.
export interface DecanterSpec {
  dirtyWaterMode?: boolean;
  freeWaterSplitFraction?: number;
  waterCarryoverFraction?: number;
  hydrocarbonToAqueousFraction?: number;
}

function buildDecanterLiquidOutlet(
  compounds: CanopyCompound[],
  componentFlows_molps: number[],
  T_K: number,
  P_Pa: number,
): Partial<MaterialStream> {
  const totalFlow_molps = componentFlows_molps.reduce((sum, value) => sum + Math.max(0, value), 0);
  const nc = compounds.length;
  const moleFractions = totalFlow_molps > 0
    ? componentFlows_molps.map(value => Math.max(0, value) / totalFlow_molps)
    : new Array(nc).fill(1 / Math.max(1, nc));
  const flash = totalFlow_molps > 0
    ? flashPT(compounds, moleFractions, T_K, P_Pa)
    : { H_Jpmol: 0, x: moleFractions, y: moleFractions };

  return {
    T_K,
    P_Pa,
    totalFlow_molps,
    moleFractions,
    phase: 'Liquid' as Phase,
    vaporFraction: 0,
    H_Jpmol: flash.H_Jpmol,
    x_liquid: flash.x,
    y_vapor: flash.y,
    solved: true,
  };
}

function dirtyWaterDecanterSplit(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  T_K: number,
  P_Pa: number,
  spec: DecanterSpec,
): { liquidI: Partial<MaterialStream>; liquidII: Partial<MaterialStream>; duty_W: number; freeWaterSeparated_molps: number; aqueousPhaseIndex: 1 | 2; waterSaturationFraction: number } | null {
  const waterIndex = compounds.findIndex(compound => {
    const name = compound.name.trim().toUpperCase();
    return name === 'WATER' || name === 'H2O';
  });
  if (waterIndex < 0) return null;

  const componentFlows = inlet.moleFractions.map(fraction => fraction * inlet.totalFlow_molps);
  const waterFeed = componentFlows[waterIndex] ?? 0;
  if (waterFeed <= 1e-12) return null;

  const aqueousFlows = new Array(compounds.length).fill(0);
  const organicFlows = new Array(compounds.length).fill(0);
  const hydrocarbonToAqueousFraction = Math.min(0.2, Math.max(0, spec.hydrocarbonToAqueousFraction ?? 0.002));
  const waterCarryoverFraction = Math.min(0.25, Math.max(0, spec.waterCarryoverFraction ?? 0.002));
  const freeWaterSplitFraction = Math.min(0.999, Math.max(0, spec.freeWaterSplitFraction ?? 0.985));

  let waterCapacityInOrganic = 0;

  for (let i = 0; i < compounds.length; i++) {
    if (i === waterIndex) continue;
    const flow = Math.max(0, componentFlows[i] ?? 0);
    if (flow <= 0) continue;

    const compound = compounds[i];
    const electrolyteLike = (compound.chargeNumber ?? 0) !== 0 || compound.isElectrolyteSolvent;

    if (electrolyteLike) {
      aqueousFlows[i] = flow;
      continue;
    }

    const hcSolubility = hcSolubilityInWater(compound, T_K);
    const dissolvedFraction = Number.isFinite(hcSolubility)
      ? Math.min(0.15, Math.max(hydrocarbonToAqueousFraction, hcSolubility * 20))
      : hydrocarbonToAqueousFraction;

    aqueousFlows[i] = flow * dissolvedFraction;
    organicFlows[i] = flow - aqueousFlows[i];

    const waterSolubility = waterSolubilityInHC(compound, T_K);
    if (Number.isFinite(waterSolubility) && waterSolubility > 0 && waterSolubility < 0.5) {
      waterCapacityInOrganic += organicFlows[i] * waterSolubility / Math.max(1 - waterSolubility, 1e-9);
    }
  }

  const targetOrganicWater = Math.min(
    waterFeed,
    Math.max(waterFeed * waterCarryoverFraction, waterCapacityInOrganic),
  );
  organicFlows[waterIndex] = targetOrganicWater;
  aqueousFlows[waterIndex] = Math.max(0, waterFeed - targetOrganicWater);

  const aqueousTotal = aqueousFlows.reduce((sum, value) => sum + value, 0);
  const organicTotal = organicFlows.reduce((sum, value) => sum + value, 0);
  if (aqueousTotal <= 1e-10 || organicTotal <= 1e-10) return null;

  const liquidI = buildDecanterLiquidOutlet(compounds, organicFlows, T_K, P_Pa);
  const liquidII = buildDecanterLiquidOutlet(compounds, aqueousFlows, T_K, P_Pa);
  const H_in = inlet.totalFlow_molps * inlet.H_Jpmol;
  const H_out = (liquidI.totalFlow_molps ?? 0) * (liquidI.H_Jpmol ?? 0)
    + (liquidII.totalFlow_molps ?? 0) * (liquidII.H_Jpmol ?? 0);

  return {
    liquidI,
    liquidII,
    duty_W: H_out - H_in,
    freeWaterSeparated_molps: Math.max(0, waterFeed * freeWaterSplitFraction),
    aqueousPhaseIndex: 2,
    waterSaturationFraction: waterCapacityInOrganic > 1e-12
      ? waterFeed / waterCapacityInOrganic
      : waterFeed > 0 ? Number.POSITIVE_INFINITY : 0,
  };
}
// ────────────────────────────────────────────────────────────────

export function solveDecanter(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  T_K: number | undefined,  // if undefined, use inlet T
  P_Pa: number | undefined, // if undefined, use inlet P
  spec: DecanterSpec = {},
  thermo: ThermoOptions = DEFAULT_THERMO,
): { liquidI: Partial<MaterialStream>; liquidII: Partial<MaterialStream>; duty_W: number; freeWaterSeparated_molps: number; aqueousPhaseIndex: 1 | 2 | null; waterSaturationFraction: number } {
  const T = T_K ?? inlet.T_K;
  const P = P_Pa ?? inlet.P_Pa;
  const nc = compounds.length;

  if (spec.dirtyWaterMode) {
    const dirtyWaterResult = dirtyWaterDecanterSplit(compounds, inlet, T, P, spec);
    if (dirtyWaterResult) return dirtyWaterResult;
  }

  // Use VLLE flash — it will detect liquid-liquid split
  const vlle = flashPT_VLLE(compounds, inlet.moleFractions, T, P,
    thermo.fluidPackage, thermo.interactionParams);

  const F = inlet.totalFlow_molps;

  if (vlle.liquidIIFraction > 1e-10) {
    // Actual liquid-liquid split
    const F_LI = F * vlle.liquidIFraction / (1 - vlle.vaporFraction || 1);
    const F_LII = F * vlle.liquidIIFraction / (1 - vlle.vaporFraction || 1);

    const flashI = flashPT(compounds, vlle.xI, T, P, thermo.fluidPackage, thermo.interactionParams);
    const flashII = flashPT(compounds, vlle.xII, T, P, thermo.fluidPackage, thermo.interactionParams);

    const H_in = F * inlet.H_Jpmol;
    const H_out = F_LI * flashI.H_Jpmol + F_LII * flashII.H_Jpmol;

    return {
      liquidI: {
        T_K: T, P_Pa: P, totalFlow_molps: F_LI,
        moleFractions: vlle.xI,
        phase: 'Liquid' as Phase, vaporFraction: 0,
        H_Jpmol: flashI.H_Jpmol,
        x_liquid: flashI.x, y_vapor: flashI.y, solved: true,
      },
      liquidII: {
        T_K: T, P_Pa: P, totalFlow_molps: F_LII,
        moleFractions: vlle.xII,
        phase: 'Liquid' as Phase, vaporFraction: 0,
        H_Jpmol: flashII.H_Jpmol,
        x_liquid: flashII.x, y_vapor: flashII.y, solved: true,
      },
      duty_W: H_out - H_in,
      freeWaterSeparated_molps: 0,
      aqueousPhaseIndex: null,
      waterSaturationFraction: 0,
    };
  }

  // No split — all goes to liquid I, liquid II is empty
  const flashAll = flashPT(compounds, inlet.moleFractions, T, P,
    thermo.fluidPackage, thermo.interactionParams);
  return {
    liquidI: {
      T_K: T, P_Pa: P, totalFlow_molps: F,
      moleFractions: [...inlet.moleFractions],
      phase: 'Liquid' as Phase, vaporFraction: 0,
      H_Jpmol: flashAll.H_Jpmol,
      x_liquid: flashAll.x, y_vapor: flashAll.y, solved: true,
    },
    liquidII: {
      T_K: T, P_Pa: P, totalFlow_molps: 0,
      moleFractions: new Array(nc).fill(1 / nc),
      phase: 'Liquid' as Phase, vaporFraction: 0,
      H_Jpmol: 0, x_liquid: new Array(nc).fill(1 / nc),
      y_vapor: new Array(nc).fill(1 / nc), solved: true,
    },
    duty_W: F * (flashAll.H_Jpmol - inlet.H_Jpmol),
    freeWaterSeparated_molps: 0,
    aqueousPhaseIndex: null,
    waterSaturationFraction: 0,
  };
}

// ────────────────────────────────────────────────────────────────
// REquil (Chemical Equilibrium Reactor)
// Finds equilibrium composition by minimizing Gibbs free energy
// via the equilibrium constant approach.
// ────────────────────────────────────────────────────────────────

export function solveREquil(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  reactions: StoichiometricReaction[],
  outletT_K: number | undefined,
  outletP_Pa: number | undefined,
  thermo: ThermoOptions = DEFAULT_THERMO,
): { outlet: Partial<MaterialStream>; duty_W: number } {
  const nc = compounds.length;
  const T = outletT_K ?? inlet.T_K;
  const P = outletP_Pa ?? inlet.P_Pa;
  const R_gas = 8.314462618; // J/(mol·K)
  const P_ref = 101325; // Pa

  // Inlet component flows
  const n_in = inlet.moleFractions.map(x => x * inlet.totalFlow_molps);
  const n_out = [...n_in];

  // For each reaction, compute equilibrium constant K(T)
  // Method: ΔG°_rxn(T) = Σ νᵢ·G_ig,i(T) → K = exp(-ΔG°/(RT))
  // With formation enthalpies in H_ig, this gives the proper K from
  // ΔG° = ΔH° - TΔS° at the reaction temperature.

  for (const rxn of reactions) {
    const stoich = rxn.coefficients; // negative for reactants, positive for products

    // Compute ΔG°_rxn at T from individual Gibbs energies
    let deltaG = 0;
    for (let i = 0; i < nc; i++) {
      if (stoich[i] !== 0) {
        // If Gf_298 is available, use proper Gibbs via gibbsIG
        // Otherwise fall back to enthalpy approximation
        if (compounds[i].Gf_298_Jmol !== undefined) {
          deltaG += stoich[i] * gibbsIG(compounds[i], T, P_ref);
        } else {
          // Fallback: assume ΔG ≈ ΔH (no entropy data)
          deltaG += stoich[i] * enthalpyIG(compounds[i], T);
        }
      }
    }

    // K_eq from ΔG° = -RT·ln(K)
    const lnK = -deltaG / (R_gas * T);
    const K_eq = Math.exp(Math.max(-50, Math.min(50, lnK)));

    // Find equilibrium extent ξ via bisection
    // n_i = n_i0 + ν_i * ξ
    // K = Π (n_i/n_total)^ν_i  (for ideal solution approximation)

    const maxExtent = (() => {
      let maxE = Infinity;
      for (let i = 0; i < nc; i++) {
        if (stoich[i] < 0 && n_out[i] > 0) {
          maxE = Math.min(maxE, n_out[i] / (-stoich[i]));
        }
      }
      return maxE === Infinity ? inlet.totalFlow_molps : maxE;
    })();

    let lo = 0, hi = maxExtent * 0.9999;
    for (let iter = 0; iter < 100; iter++) {
      const mid = (lo + hi) / 2;
      const n_try = n_out.map((ni, i) => Math.max(1e-20, ni + stoich[i] * mid));
      const total = n_try.reduce((a, b) => a + b, 0);
      const x_try = n_try.map(ni => ni / total);

      // Compute activity coefficients for non-ideal K
      // Aspen: K = Π(a_i^ν_i) where a_i = γ_i·x_i (liquid) or φ_i·y_i·P/P° (gas)
      // For liquid-phase reactions with activity model:
      const gamma = computeGamma(compounds, x_try, T, thermo.fluidPackage, thermo.interactionParams);

      let lnQ = 0;
      for (let i = 0; i < nc; i++) {
        if (stoich[i] !== 0) {
          // a_i = γ_i · x_i for liquid species
          const activity = gamma[i] * x_try[i];
          lnQ += stoich[i] * Math.log(Math.max(activity, 1e-30));
        }
      }

      if (lnQ < lnK) {
        lo = mid;  // reaction needs to proceed further
      } else {
        hi = mid;  // reaction has gone too far
      }
      if (hi - lo < 1e-12 * maxExtent) break;
    }

    const xi = (lo + hi) / 2;
    for (let i = 0; i < nc; i++) {
      n_out[i] = Math.max(0, n_out[i] + stoich[i] * xi);
    }
  }

  const totalOut = n_out.reduce((a, b) => a + b, 0);
  const z_out = totalOut > 0 ? n_out.map(n => n / totalOut) : new Array(nc).fill(1 / nc);

  const flash = flashPT(compounds, z_out, T, P, thermo.fluidPackage, thermo.interactionParams);
  const H_in = inlet.totalFlow_molps * inlet.H_Jpmol;
  const H_out = totalOut * flash.H_Jpmol;

  return {
    outlet: {
      T_K: T, P_Pa: P, totalFlow_molps: totalOut,
      moleFractions: z_out,
      phase: flash.phase, vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x, y_vapor: flash.y, solved: true,
    },
    duty_W: H_out - H_in,
  };
}

// ────────────────────────────────────────────────────────────────
// 20b. RGibbs — Gibbs Free Energy Minimization Reactor
// ────────────────────────────────────────────────────────────────
/**
 * RGibbs reactor: finds equilibrium composition by minimizing total Gibbs energy.
 * Does NOT require user-specified reactions — all possible products are considered.
 *
 * Method: Steepest-descent + Lagrange multiplier for element balance.
 * For each component: μ_i = μ_i° + RT·ln(y_i·P/P°)  (ideal gas approximation)
 *   where μ_i° = H_ig(T) - T·S_ig(T, P°)
 *
 * Element conservation is enforced: A·n = b (where A is the formula matrix, b is element totals).
 */
export function solveRGibbs(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  outletT_K: number | undefined,
  outletP_Pa: number | undefined,
  /** Element formula matrix: elementMatrix[i][e] = number of atoms of element e in compound i */
  elementMatrix: number[][],
  thermo: ThermoOptions = DEFAULT_THERMO,
): { outlet: Partial<MaterialStream>; duty_W: number } {
  const nc = compounds.length;
  const T = outletT_K ?? inlet.T_K;
  const P = outletP_Pa ?? inlet.P_Pa;
  const R_gas = 8.314462618;
  const P_ref = 101325;

  // Total inlet moles
  const n_total_in = inlet.totalFlow_molps;
  const n_in = inlet.moleFractions.map(x => x * n_total_in);

  // Element totals from inlet (conservation constraint)
  const nElements = elementMatrix[0]?.length ?? 0;
  const b_elem = new Array(nElements).fill(0);
  for (let i = 0; i < nc; i++) {
    for (let e = 0; e < nElements; e++) {
      b_elem[e] += n_in[i] * (elementMatrix[i]?.[e] ?? 0);
    }
  }

  // Compute standard chemical potential for each species at T
  // μ°_i = H_ig(T) - T·S_ig(T, P_ref) ≈ ∫Cp dT - T·∫Cp/T dT + T·R·ln(1)
  const mu0 = new Array(nc).fill(0);
  for (let i = 0; i < nc; i++) {
    const H_i = enthalpyIG(compounds[i], T);         // J/mol
    const S_i = entropyIG(compounds[i], T, P_ref);    // J/(mol·K)
    mu0[i] = H_i - T * S_i;
  }

  // Initialize: start from inlet composition
  let n_out = [...n_in];
  // Ensure all moles are positive
  for (let i = 0; i < nc; i++) {
    if (n_out[i] < 1e-15) n_out[i] = 1e-15;
  }

  // Gibbs minimization by successive approximation with Lagrange multipliers
  for (let outer = 0; outer < 200; outer++) {
    const n_total = n_out.reduce((a, b) => a + b, 0);
    if (n_total < 1e-30) break;

    // Chemical potential of each species
    const mu = new Array(nc).fill(0);
    for (let i = 0; i < nc; i++) {
      const y_i = n_out[i] / n_total;
      mu[i] = mu0[i] + R_gas * T * Math.log(Math.max(y_i, 1e-30) * P / P_ref);
    }

    // Solve Lagrange multiplier system: A^T · λ = -μ (least squares)
    // where λ are element multipliers
    // This minimizes G subject to element balance
    // Using normal equations: (A^T A) λ = A^T (-μ)

    // Build A^T A and A^T (-μ)
    const ATA = Array.from({ length: nElements }, () => new Array(nElements).fill(0));
    const ATmu = new Array(nElements).fill(0);
    for (let e1 = 0; e1 < nElements; e1++) {
      for (let e2 = 0; e2 < nElements; e2++) {
        for (let i = 0; i < nc; i++) {
          ATA[e1][e2] += (elementMatrix[i]?.[e1] ?? 0) * (elementMatrix[i]?.[e2] ?? 0) * n_out[i];
        }
      }
      for (let i = 0; i < nc; i++) {
        ATmu[e1] -= (elementMatrix[i]?.[e1] ?? 0) * mu[i] * n_out[i];
      }
    }

    // Solve ATA·λ = ATmu by Gaussian elimination
    const lambda = solveLinearSystem(ATA, ATmu);
    if (!lambda) break;

    // Update moles: n_i_new = n_i · exp[-{μ_i + Σ_e λ_e·a_ie}/(RT)]
    let maxChange = 0;
    const n_new = new Array(nc).fill(0);
    for (let i = 0; i < nc; i++) {
      let correction = mu[i];
      for (let e = 0; e < nElements; e++) {
        correction += lambda[e] * (elementMatrix[i]?.[e] ?? 0);
      }
      // Damped update to prevent oscillation
      const factor = Math.exp(-correction / (R_gas * T));
      n_new[i] = n_out[i] * Math.min(Math.max(factor, 0.1), 10);
      if (n_new[i] < 1e-20) n_new[i] = 1e-20;
      maxChange = Math.max(maxChange, Math.abs(n_new[i] - n_out[i]) / (n_out[i] + 1e-30));
    }

    // Project onto element balance constraint (scale to match element totals)
    const n_total_new = n_new.reduce((a, b) => a + b, 0);
    // Check element balance
    for (let e = 0; e < nElements; e++) {
      let b_curr = 0;
      for (let i = 0; i < nc; i++) {
        b_curr += n_new[i] * (elementMatrix[i]?.[e] ?? 0);
      }
      if (Math.abs(b_curr) > 1e-30 && Math.abs(b_elem[e]) > 1e-30) {
        const scale = b_elem[e] / b_curr;
        // Only scale if departure is significant
        if (Math.abs(scale - 1) > 0.01) {
          for (let i = 0; i < nc; i++) {
            if ((elementMatrix[i]?.[e] ?? 0) > 0) {
              n_new[i] *= scale;
            }
          }
        }
      }
    }

    n_out = n_new;
    if (maxChange < 1e-8) break;
  }

  const totalOut = n_out.reduce((a, b) => a + b, 0);
  const z_out = totalOut > 0 ? n_out.map(n => n / totalOut) : new Array(nc).fill(1 / nc);

  const flash = flashPT(compounds, z_out, T, P, thermo.fluidPackage, thermo.interactionParams);
  const H_in = n_total_in * inlet.H_Jpmol;
  const H_out = totalOut * flash.H_Jpmol;

  return {
    outlet: {
      T_K: T, P_Pa: P, totalFlow_molps: totalOut,
      moleFractions: z_out,
      phase: flash.phase, vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x, y_vapor: flash.y, solved: true,
    },
    duty_W: H_out - H_in,
  };
}

/** Solve n×n linear system Ax=b by Gaussian elimination with partial pivoting */
function solveLinearSystem(A: number[][], b: number[]): number[] | null {
  const n = A.length;
  // Augmented matrix
  const M = A.map((row, i) => [...row, b[i]]);

  for (let col = 0; col < n; col++) {
    // Partial pivoting
    let maxRow = col;
    let maxVal = Math.abs(M[col][col]);
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(M[row][col]) > maxVal) {
        maxVal = Math.abs(M[row][col]);
        maxRow = row;
      }
    }
    if (maxVal < 1e-30) return null;
    if (maxRow !== col) [M[col], M[maxRow]] = [M[maxRow], M[col]];

    // Eliminate below
    for (let row = col + 1; row < n; row++) {
      const factor = M[row][col] / M[col][col];
      for (let j = col; j <= n; j++) {
        M[row][j] -= factor * M[col][j];
      }
    }
  }

  // Back substitution
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = M[i][n];
    for (let j = i + 1; j < n; j++) {
      x[i] -= M[i][j] * x[j];
    }
    x[i] /= M[i][i];
  }
  return x;
}

// ────────────────────────────────────────────────────────────────
// 21. MHeatX — Multi-Stream Heat Exchanger
// ────────────────────────────────────────────────────────────────
export interface MHeatXSpec {
  hotStreamIndices: number[];
  coldStreamIndices: number[];
  deltaT_min_K: number;
  UA_WpK?: number;
}

export function solveMHeatX(
  compounds: CanopyCompound[],
  inlets: MaterialStream[],
  spec: MHeatXSpec,
  thermo: ThermoOptions,
): { outlets: MaterialStream[]; duty_W: number } {
  const { hotStreamIndices, coldStreamIndices, deltaT_min_K } = spec;

  const coldInletT_min = Math.min(...coldStreamIndices.map(i => inlets[i].T_K));
  const hotInletT_max = Math.max(...hotStreamIndices.map(i => inlets[i].T_K));

  let Q_hot_max = 0;
  let Q_cold_max = 0;

  for (const hi of hotStreamIndices) {
    const s = inlets[hi];
    const T_out_min = coldInletT_min + deltaT_min_K;
    const H_hot = s.totalFlow_molps * s.H_Jpmol;
    const flash_cold = flashPT(compounds, s.moleFractions,
      Math.max(T_out_min, coldInletT_min + deltaT_min_K), s.P_Pa,
      thermo.fluidPackage, thermo.interactionParams);
    const H_cold = s.totalFlow_molps * flash_cold.H_Jpmol;
    Q_hot_max += H_hot - H_cold;
  }

  for (const ci of coldStreamIndices) {
    const s = inlets[ci];
    const T_out_max = hotInletT_max - deltaT_min_K;
    const H_cold_in = s.totalFlow_molps * s.H_Jpmol;
    const flash_hot = flashPT(compounds, s.moleFractions,
      Math.min(T_out_max, hotInletT_max - deltaT_min_K), s.P_Pa,
      thermo.fluidPackage, thermo.interactionParams);
    const H_hot_out = s.totalFlow_molps * flash_hot.H_Jpmol;
    Q_cold_max += H_hot_out - H_cold_in;
  }

  const Q = Math.min(Q_hot_max, Q_cold_max);

  const outlets: MaterialStream[] = [];
  for (let idx = 0; idx < inlets.length; idx++) {
    const s = inlets[idx];
    const isHot = hotStreamIndices.includes(idx);
    const isCold = coldStreamIndices.includes(idx);

    if (isHot) {
      const Q_stream = Q / hotStreamIndices.length;
      const H_out_per_mol = (s.totalFlow_molps * s.H_Jpmol - Q_stream) / s.totalFlow_molps;
      const result = flashPQ(compounds, s.moleFractions, s.P_Pa, H_out_per_mol, 0,
        s.totalFlow_molps, thermo.fluidPackage, thermo.interactionParams, s.T_K);
      outlets.push({
        ...s, T_K: result.T_K, phase: result.phase, vaporFraction: result.vaporFraction,
        H_Jpmol: result.H_Jpmol, x_liquid: result.x, y_vapor: result.y, solved: true,
      });
    } else if (isCold) {
      const Q_stream = Q / coldStreamIndices.length;
      const H_out_per_mol = (s.totalFlow_molps * s.H_Jpmol + Q_stream) / s.totalFlow_molps;
      const result = flashPQ(compounds, s.moleFractions, s.P_Pa, H_out_per_mol, 0,
        s.totalFlow_molps, thermo.fluidPackage, thermo.interactionParams, s.T_K);
      outlets.push({
        ...s, T_K: result.T_K, phase: result.phase, vaporFraction: result.vaporFraction,
        H_Jpmol: result.H_Jpmol, x_liquid: result.x, y_vapor: result.y, solved: true,
      });
    } else {
      outlets.push({ ...s, solved: true });
    }
  }

  return { outlets, duty_W: Q };
}

// ────────────────────────────────────────────────────────────────
// 22. RBatch — Batch Reactor (ODE integration in time)
// ────────────────────────────────────────────────────────────────
export interface RBatchSpec {
  reactions: StoichiometricReaction[];
  rateConstants_1ps: number[];
  Ea_Jpmol: number[];
  T_ref_K: number;
  batchTime_s: number;
  nSteps?: number;
  isothermal?: boolean;
  duty_W?: number;
}

export function solveRBatch(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: RBatchSpec,
  thermo: ThermoOptions,
): { outlet: MaterialStream; duty_W: number; conversionProfile: number[] } {
  const { reactions, rateConstants_1ps, Ea_Jpmol, T_ref_K, batchTime_s } = spec;
  const nSteps = spec.nSteps ?? 200;
  const isothermal = spec.isothermal ?? true;
  const nc = compounds.length;
  const R_gas = 8.314462618;

  const N = inlet.moleFractions.map(x => x * inlet.totalFlow_molps);
  let T = inlet.T_K;
  const dt = batchTime_s / nSteps;
  const conversionProfile: number[] = [];
  const N0_key = N[reactions[0]?.keyComponentIndex ?? 0];
  let totalDuty = 0;

  for (let step = 0; step < nSteps; step++) {
    const keyIdx = reactions[0]?.keyComponentIndex ?? 0;
    conversionProfile.push(N0_key > 0 ? 1 - N[keyIdx] / N0_key : 0);

    const totalMols = N.reduce((s, n) => s + n, 0);
    const x = N.map(n => n / totalMols);
    const dN = new Array(nc).fill(0);
    let dH_rxn = 0;

    for (let r = 0; r < reactions.length; r++) {
      const rxn = reactions[r];
      const k0 = rateConstants_1ps[r] ?? 0.01;
      const Ea = Ea_Jpmol[r] ?? 0;
      const k = k0 * Math.exp(-Ea / R_gas * (1 / T - 1 / T_ref_K));

      let minRatio = Infinity;
      for (let i = 0; i < nc && i < rxn.coefficients.length; i++) {
        if (rxn.coefficients[i] < 0) {
          const ratio = N[i] / Math.abs(rxn.coefficients[i]);
          if (ratio < minRatio) minRatio = ratio;
        }
      }
      const extent = k * Math.max(minRatio, 0);

      for (let i = 0; i < nc && i < rxn.coefficients.length; i++) {
        dN[i] += rxn.coefficients[i] * extent;
      }

      let dH = 0;
      for (let i = 0; i < nc && i < rxn.coefficients.length; i++) {
        dH += rxn.coefficients[i] * enthalpyIG(compounds[i], T);
      }
      dH_rxn += dH * extent;
    }

    for (let i = 0; i < nc; i++) {
      N[i] = Math.max(0, N[i] + dN[i] * dt);
    }

    if (!isothermal) {
      const totalMolsNow = N.reduce((s, n) => s + n, 0);
      const Cp_total = x.reduce((s, xi, i) => s + xi * CpIG_JmolK(compounds[i], T), 0) * totalMolsNow;
      if (Cp_total > 0) {
        T += (-dH_rxn + (spec.duty_W ?? 0)) * dt / Cp_total;
      }
    } else {
      totalDuty += dH_rxn * dt;
    }
  }

  const keyIdx = reactions[0]?.keyComponentIndex ?? 0;
  conversionProfile.push(N0_key > 0 ? 1 - N[keyIdx] / N0_key : 0);

  const totalMolsOut = N.reduce((s, n) => s + n, 0);
  const zOut = N.map(n => n / (totalMolsOut || 1));
  const flash = flashPT(compounds, zOut, T, inlet.P_Pa, thermo.fluidPackage, thermo.interactionParams);

  return {
    outlet: {
      ...inlet,
      T_K: T, P_Pa: inlet.P_Pa, totalFlow_molps: totalMolsOut,
      moleFractions: zOut, phase: flash.phase, vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol, x_liquid: flash.x, y_vapor: flash.y, solved: true,
    },
    duty_W: isothermal ? totalDuty / batchTime_s : (spec.duty_W ?? 0),
    conversionProfile,
  };
}

// ────────────────────────────────────────────────────────────────
// 23. Crystallizer (MSMPR — Mixed-Suspension Mixed-Product Removal)
// ────────────────────────────────────────────────────────────────
export interface CrystallizerSpec {
  residenceTime_s: number;
  soluteIndex: number;
  C_sat_molpm3: number;
  k_g: number;
  g_exponent: number;
  k_b: number;
  b_exponent: number;
  j_exponent?: number;
  k_v?: number;
  rho_crystal_kgpm3?: number;
  liquidDensity_kgpm3?: number;
  T_K?: number;
}

export function solveCrystallizer(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: CrystallizerSpec,
  thermo: ThermoOptions,
): {
  outlet: MaterialStream;
  crystalFlow_kgps: number;
  meanCrystalSize_m: number;
  duty_W: number;
  supersaturation_molpm3: number;
  growthRate_mps: number;
  nucleationRate_1pm3ps: number;
  moments: { mu0: number; mu1: number; mu2: number; mu3: number };
} {
  const {
    residenceTime_s: tau, soluteIndex, C_sat_molpm3,
    k_g, g_exponent: g, k_b, b_exponent: b,
    j_exponent: j = 0.5,
    k_v = Math.PI / 6, rho_crystal_kgpm3: rho_c = 2000,
  } = spec;
  const T = spec.T_K ?? inlet.T_K;
  const nc = compounds.length;

  const avgMW = compounds.reduce((s, c, i) => s + c.molecularWeight * inlet.moleFractions[i], 0);
  const liquidDensity_kgpm3 = spec.liquidDensity_kgpm3 ?? Math.max(50, avgMW / Math.max(mixtureLiquidMolarVolume(compounds, inlet.moleFractions, T), 1e-9));
  const V_feed = inlet.totalFlow_molps * avgMW / (liquidDensity_kgpm3 * 1000);
  const C_feed = inlet.moleFractions[soluteIndex] * inlet.totalFlow_molps / (V_feed || 1e-10);

  const deltaC = Math.max(C_feed - C_sat_molpm3, 0);
  const G = k_g * Math.pow(Math.max(deltaC, 0), g);
  const Gtau = G * tau;

  // Analytical solution for MSMPR with secondary nucleation B0 = kb * ΔC^b * μ3^0.5:
  //   μ3 = [kb * ΔC^b * (Gτ)^4 * 6 / G]^2
  const coeff = (G > 0) ? 6 * (k_b * Math.pow(deltaC, b) / G) * Math.pow(Gtau, 4) : 0;
  const mu3 = coeff > 0 && j < 1 ? Math.pow(coeff, 1 / (1 - j)) : (j === 0 ? coeff : 0);
  const B0 = (G > 0) ? k_b * Math.pow(deltaC, b) * Math.pow(Math.max(mu3, 0), j) : 0;
  const n0 = (G > 0) ? B0 / G : 0;
  const mu0 = n0 * Gtau;
  const mu1 = n0 * Math.pow(Gtau, 2);
  const mu2 = 2 * n0 * Math.pow(Gtau, 3);
  const crystalVolFrac = k_v * mu3;
  let crystalMassRate = crystalVolFrac * rho_c * (V_feed || 0);
  const L_mean = mu0 > 0 ? 4 * Gtau : 0;

  // Cap crystal removal: cannot remove more solute than available (convert g→kg)
  const maxSoluteMass = inlet.moleFractions[soluteIndex] * inlet.totalFlow_molps * (compounds[soluteIndex]?.molecularWeight ?? 100) / 1000;
  crystalMassRate = Math.min(crystalMassRate, maxSoluteMass * 0.99); // leave at least 1%
  if (!isFinite(crystalMassRate)) crystalMassRate = maxSoluteMass * 0.99;

  const molesRemoved = crystalMassRate * 1000 / (compounds[soluteIndex]?.molecularWeight ?? 100);
  const outMoleFracs = [...inlet.moleFractions];
  const totalOut = Math.max(inlet.totalFlow_molps - molesRemoved, 1e-20);
  for (let i = 0; i < nc; i++) {
    if (i === soluteIndex) {
      outMoleFracs[i] = Math.max(0, inlet.moleFractions[i] * inlet.totalFlow_molps - molesRemoved) / totalOut;
    } else {
      outMoleFracs[i] = inlet.moleFractions[i] * inlet.totalFlow_molps / totalOut;
    }
  }
  const sum = outMoleFracs.reduce((s, x) => s + x, 0);
  for (let i = 0; i < nc; i++) outMoleFracs[i] /= (sum || 1);

  const flash = flashPT(compounds, outMoleFracs, T, inlet.P_Pa, thermo.fluidPackage, thermo.interactionParams);
  const H_in = inlet.totalFlow_molps * inlet.H_Jpmol;
  const H_out = totalOut * flash.H_Jpmol;

  return {
    outlet: {
      ...inlet,
      T_K: T, P_Pa: inlet.P_Pa, totalFlow_molps: totalOut,
      moleFractions: outMoleFracs, phase: flash.phase, vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol, x_liquid: flash.x, y_vapor: flash.y, solved: true,
    },
    crystalFlow_kgps: crystalMassRate,
    meanCrystalSize_m: L_mean,
    duty_W: H_out - H_in,
    supersaturation_molpm3: deltaC,
    growthRate_mps: G,
    nucleationRate_1pm3ps: B0,
    moments: { mu0, mu1, mu2, mu3 },
  };
}

// ────────────────────────────────────────────────────────────────
// 24. Crusher — Size Reduction (Bond's Law)
// ────────────────────────────────────────────────────────────────
export interface CrusherSpec {
  workIndex_kWhpton: number;
  feedSize_um: number;
  productSize_um: number;
  solidFlow_kgps: number;
}

export function solveCrusher(
  _compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: CrusherSpec,
): { outlet: MaterialStream; power_W: number } {
  const { workIndex_kWhpton, feedSize_um, productSize_um, solidFlow_kgps } = spec;

  // Bond's Law: W = Wi * 10 * (1/sqrt(P80) - 1/sqrt(F80))  [kWh/ton]
  const W_kWhpton = workIndex_kWhpton * 10 * (
    1 / Math.sqrt(productSize_um) - 1 / Math.sqrt(feedSize_um)
  );

  const solidFlow_tonps = solidFlow_kgps / 1000;
  const power_W = W_kWhpton * solidFlow_tonps * 3600 * 1000;

  return {
    outlet: { ...inlet, solved: true },
    power_W,
  };
}

// ────────────────────────────────────────────────────────────────
// 25. Dryer — Convective Drying (Constant + Falling Rate)
// ────────────────────────────────────────────────────────────────
export interface DryerSpec {
  X_in_kgpkg: number;
  X_out_kgpkg: number;
  X_c_kgpkg: number;
  X_eq_kgpkg?: number;
  dryFlow_kgps: number;
  gasT_in_K: number;
  gasFlow_kgps: number;
  Lv_Jpkg?: number;
  h_Wpm2K?: number;
  area_m2?: number;
  wetBulbT_K?: number;
  gasHumidityIn_kgpkg?: number;
  Cp_gas_JpkgK?: number;
}

export function solveDryer(
  _compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: DryerSpec,
): {
  outlet: MaterialStream;
  waterRemoved_kgps: number;
  gasT_out_K: number;
  duty_W: number;
  gasHumidityOut_kgpkg: number;
  dryingTime_s: number;
  constantRateTime_s: number;
  fallingRateTime_s: number;
  constantRate_kgpm2s: number;
} {
  const {
    X_in_kgpkg: X_in, X_out_kgpkg: X_out, X_c_kgpkg: X_c,
    X_eq_kgpkg: X_eq = 0,
    dryFlow_kgps, gasT_in_K,
    gasFlow_kgps, Lv_Jpkg = 2260e3,
    h_Wpm2K = 35,
    area_m2 = 1,
    wetBulbT_K = inlet.T_K,
    gasHumidityIn_kgpkg = 0,
    Cp_gas_JpkgK = 1006,
  } = spec;

  const X_target = Math.max(X_eq, Math.min(X_in, X_out));
  const waterRemoved_kgps = Math.max(0, dryFlow_kgps * (X_in - X_target));
  const duty_W = waterRemoved_kgps * Lv_Jpkg;

  const Cp_air = 1006; // J/(kg·K)
  const constantRate_kgpm2s = Math.max(0, h_Wpm2K * Math.max(gasT_in_K - wetBulbT_K, 0) / Lv_Jpkg);
  const gasT_out_K = gasT_in_K - duty_W / Math.max(gasFlow_kgps * Cp_gas_JpkgK, 1e-9);
  const gasHumidityOut_kgpkg = gasHumidityIn_kgpkg + waterRemoved_kgps / Math.max(gasFlow_kgps, 1e-9);

  let constantRateTime_s = 0;
  let fallingRateTime_s = 0;
  if (constantRate_kgpm2s > 0 && area_m2 > 0) {
    if (X_in > X_c) {
      const constDrop = Math.max(X_in - Math.max(X_target, X_c), 0);
      constantRateTime_s = dryFlow_kgps * constDrop / (area_m2 * constantRate_kgpm2s);
    }
    if (X_target < Math.min(X_in, X_c)) {
      const xStart = Math.min(X_in, X_c);
      const numerator = Math.max(xStart - X_eq, 1e-9);
      const denominator = Math.max(X_target - X_eq, 1e-9);
      const scale = Math.max(X_c - X_eq, 1e-9);
      fallingRateTime_s = dryFlow_kgps * scale * Math.log(numerator / denominator) / (area_m2 * constantRate_kgpm2s);
    }
  }

  return {
    outlet: { ...inlet, solved: true },
    waterRemoved_kgps,
    gasT_out_K,
    duty_W,
    gasHumidityOut_kgpkg,
    dryingTime_s: constantRateTime_s + fallingRateTime_s,
    constantRateTime_s,
    fallingRateTime_s,
    constantRate_kgpm2s,
  };
}

export interface MembraneSpec {
  model?: 'stage-cut' | 'rejection' | 'solution-diffusion';
  stageCut: number;
  selectivities?: number[];
  rejectionCoefficients?: number[];
  membraneArea_m2?: number;
  permeanceGPU?: number[];
  foulingFactor?: number;
  permeateP_Pa?: number;
  retentateP_Pa?: number;
}

export function solveMembrane(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: MembraneSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): { retentate: MaterialStream; permeate: MaterialStream } {
  const stageCut = Math.max(0.01, Math.min(0.95, spec.stageCut));
  const feedMoles = inlet.moleFractions.map(value => value * inlet.totalFlow_molps);
  const model = spec.model ?? 'stage-cut';
  const permeatePressure = spec.permeateP_Pa ?? Math.max(1000, inlet.P_Pa * 0.5);
  const retentatePressure = spec.retentateP_Pa ?? inlet.P_Pa;
  let permeateMoles = new Array(feedMoles.length).fill(0);

  if (model === 'solution-diffusion') {
    const area = Math.max(spec.membraneArea_m2 ?? 50, 1e-6);
    const fouling = Math.max(0.05, Math.min(1, spec.foulingFactor ?? 1));
    const gpuToMolPerM2sPa = 3.348e-10;
    const permeance = compounds.map((_, i) => Math.max(spec.permeanceGPU?.[i] ?? 100, 1e-9) * gpuToMolPerM2sPa);
    let trialPermeate = [...inlet.moleFractions];
    for (let iter = 0; iter < 4; iter++) {
      const fluxes = compounds.map((_, i) =>
        Math.max(0, permeance[i] * area * fouling * (inlet.P_Pa * inlet.moleFractions[i] - permeatePressure * trialPermeate[i])),
      );
      const totalFlux = fluxes.reduce((sum, value) => sum + value, 0);
      if (totalFlux <= 1e-18) break;
      trialPermeate = fluxes.map(value => value / totalFlux);
      permeateMoles = fluxes.map(value => value * stageCut * inlet.totalFlow_molps / Math.max(totalFlux, 1e-18));
    }
  } else if (model === 'rejection') {
    const rejections = compounds.map((_, i) => Math.max(0, Math.min(0.999, spec.rejectionCoefficients?.[i] ?? 0)));
    const weights = inlet.moleFractions.map((zi, i) => zi * Math.max(1 - rejections[i], 1e-6));
    const avgWeight = weights.reduce((sum, value) => sum + value, 0) || 1;
    permeateMoles = feedMoles.map((mol, i) => {
      const split = Math.max(0, Math.min(0.999, stageCut * weights[i] / avgWeight));
      return mol * split;
    });
  } else {
    const selectivities = compounds.map((_, i) => Math.max(spec.selectivities?.[i] ?? 1, 1e-6));
    const averageSelectivity = inlet.moleFractions.reduce(
      (sum, zi, i) => sum + zi * selectivities[i],
      0,
    ) || 1;
    permeateMoles = feedMoles.map((mol, i) => {
      const split = Math.max(0, Math.min(0.999, stageCut * selectivities[i] / averageSelectivity));
      return mol * split;
    });
  }
  const retentateMoles = feedMoles.map((mol, i) => Math.max(0, mol - permeateMoles[i]));

  return {
    retentate: buildStreamFromMoles(
      compounds,
      retentateMoles,
      inlet.T_K,
      retentatePressure,
      thermo,
    ),
    permeate: buildStreamFromMoles(
      compounds,
      permeateMoles,
      inlet.T_K,
      permeatePressure,
      thermo,
    ),
  };
}

export interface CycloneSpec {
  cutSize_um: number;
  particleSizes_um?: number[];
  pressureDrop_Pa?: number;
}

export function solveCyclone(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: CycloneSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): { overflow: MaterialStream; underflow: MaterialStream; pressureDrop_Pa: number } {
  const cutSize = Math.max(spec.cutSize_um, 1e-6);
  const pressureDrop_Pa = Math.max(spec.pressureDrop_Pa ?? 1500, 0);
  const feedMoles = inlet.moleFractions.map(value => value * inlet.totalFlow_molps);
  const underflowMoles = feedMoles.map((mol, i) => {
    const dp = Math.max(spec.particleSizes_um?.[i] ?? 0, 0);
    const efficiency = dp > 0 ? 1 - Math.exp(-Math.pow(dp / cutSize, 2)) : 0;
    return mol * Math.max(0, Math.min(0.999, efficiency));
  });
  const overflowMoles = feedMoles.map((mol, i) => Math.max(0, mol - underflowMoles[i]));
  const outletP = Math.max(1000, inlet.P_Pa - pressureDrop_Pa);

  return {
    overflow: buildStreamFromMoles(compounds, overflowMoles, inlet.T_K, outletP, thermo),
    underflow: buildStreamFromMoles(compounds, underflowMoles, inlet.T_K, outletP, thermo, 'Liquid'),
    pressureDrop_Pa,
  };
}

export interface FilterSpec {
  retainedFractions: number[];
  pressureDrop_Pa?: number;
}

export function solveFilter(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: FilterSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): { filtrate: MaterialStream; cake: MaterialStream; pressureDrop_Pa: number } {
  const pressureDrop_Pa = Math.max(spec.pressureDrop_Pa ?? 5000, 0);
  const feedMoles = inlet.moleFractions.map(value => value * inlet.totalFlow_molps);
  const cakeMoles = feedMoles.map((mol, i) => mol * Math.max(0, Math.min(1, spec.retainedFractions[i] ?? 0)));
  const filtrateMoles = feedMoles.map((mol, i) => Math.max(0, mol - cakeMoles[i]));
  const outletP = Math.max(1000, inlet.P_Pa - pressureDrop_Pa);

  return {
    filtrate: buildStreamFromMoles(compounds, filtrateMoles, inlet.T_K, outletP, thermo),
    cake: buildStreamFromMoles(compounds, cakeMoles, inlet.T_K, outletP, thermo, 'Liquid'),
    pressureDrop_Pa,
  };
}

export interface ScreenSpec {
  cutSize_um: number;
  sharpness?: number;
  particleSizes_um?: number[];
}

export function solveScreen(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: ScreenSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): { oversize: MaterialStream; undersize: MaterialStream } {
  const cutSize = Math.max(spec.cutSize_um, 1e-6);
  const sharpness = Math.max(spec.sharpness ?? 4, 0.1);
  const feedMoles = inlet.moleFractions.map(value => value * inlet.totalFlow_molps);
  const oversizeMoles = feedMoles.map((mol, i) => {
    const dp = Math.max(spec.particleSizes_um?.[i] ?? 0, 0);
    const partition = dp > 0 ? 1 / (1 + Math.exp(-sharpness * (dp / cutSize - 1))) : 0;
    return mol * partition;
  });
  const undersizeMoles = feedMoles.map((mol, i) => Math.max(0, mol - oversizeMoles[i]));

  return {
    oversize: buildStreamFromMoles(compounds, oversizeMoles, inlet.T_K, inlet.P_Pa, thermo, 'Liquid'),
    undersize: buildStreamFromMoles(compounds, undersizeMoles, inlet.T_K, inlet.P_Pa, thermo),
  };
}

export interface CentrifugeSpec {
  cutSize_um: number;
  sharpness?: number;
  particleSizes_um?: number[];
  liquidSplitToCake?: number;
  pressureDrop_Pa?: number;
}

export function solveCentrifuge(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: CentrifugeSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): { centrate: MaterialStream; cake: MaterialStream; pressureDrop_Pa: number } {
  const cutSize = Math.max(spec.cutSize_um, 1e-6);
  const sharpness = Math.max(spec.sharpness ?? 7, 0.1);
  const liquidSplitToCake = Math.max(0, Math.min(spec.liquidSplitToCake ?? 0.08, 0.5));
  const pressureDrop_Pa = Math.max(spec.pressureDrop_Pa ?? 5000, 0);
  const feedMoles = inlet.moleFractions.map(value => value * inlet.totalFlow_molps);
  const cakeMoles = feedMoles.map((mol, i) => {
    const dp = Math.max(spec.particleSizes_um?.[i] ?? 0, 0);
    const solidsCapture = dp > 0
      ? 1 / (1 + Math.exp(-sharpness * (dp / cutSize - 1)))
      : 0;
    const entrainment = dp > 0
      ? liquidSplitToCake * Math.exp(-dp / cutSize)
      : liquidSplitToCake;
    const split = Math.max(0, Math.min(0.999, solidsCapture + (1 - solidsCapture) * entrainment));
    return mol * split;
  });
  const centrateMoles = feedMoles.map((mol, i) => Math.max(0, mol - cakeMoles[i]));
  const outletP = Math.max(1000, inlet.P_Pa - pressureDrop_Pa);

  return {
    centrate: buildStreamFromMoles(compounds, centrateMoles, inlet.T_K, outletP, thermo),
    cake: buildStreamFromMoles(compounds, cakeMoles, inlet.T_K, outletP, thermo, 'Liquid'),
    pressureDrop_Pa,
  };
}

export interface SteamDrumSpec {
  outletP_Pa: number;
}

export function solveSteamDrum(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: SteamDrumSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): { vapor: MaterialStream; liquid: MaterialStream; duty_W: number } {
  const flash = flashPT(compounds, inlet.moleFractions, inlet.T_K, spec.outletP_Pa, thermo.fluidPackage, thermo.interactionParams);
  const vaporMoles = flash.y.map(value => value * inlet.totalFlow_molps * flash.vaporFraction);
  const liquidMoles = flash.x.map(value => value * inlet.totalFlow_molps * (1 - flash.vaporFraction));
  return {
    vapor: buildStreamFromMoles(compounds, vaporMoles, inlet.T_K, spec.outletP_Pa, thermo, 'Vapor'),
    liquid: buildStreamFromMoles(compounds, liquidMoles, inlet.T_K, spec.outletP_Pa, thermo, 'Liquid'),
    duty_W: inlet.totalFlow_molps * (flash.H_Jpmol - inlet.H_Jpmol),
  };
}

export interface SteamHeaterSpec {
  targetT_K: number;
  outletP_Pa?: number;
  steamLevel?: 'LP' | 'MP' | 'HP' | 'Auto';
}

export function solveSteamHeater(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: SteamHeaterSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): {
  outlet: Partial<MaterialStream>;
  duty_W: number;
  steamFlow_kgps: number;
  utilityType: 'lpSteam' | 'mpSteam' | 'hpSteam';
} {
  const result = solveHeater(compounds, inlet, {
    targetT_K: spec.targetT_K,
    outletP_Pa: spec.outletP_Pa ?? inlet.P_Pa,
  }, thermo);
  const utilityType = spec.steamLevel === 'HP'
    ? 'hpSteam'
    : spec.steamLevel === 'MP'
      ? 'mpSteam'
      : spec.steamLevel === 'LP'
        ? 'lpSteam'
        : spec.targetT_K > 523
          ? 'hpSteam'
          : spec.targetT_K > 423
            ? 'mpSteam'
            : 'lpSteam';
  const latentHeat_Jpkg = utilityType === 'hpSteam'
    ? 1.90e6
    : utilityType === 'mpSteam'
      ? 2.00e6
      : 2.10e6;
  const steamFlow_kgps = Math.max(result.duty_W, 0) / latentHeat_Jpkg;

  return {
    outlet: result.outlet,
    duty_W: result.duty_W,
    steamFlow_kgps,
    utilityType,
  };
}

export interface SteamTurbineSpec {
  outletP_Pa: number;
  efficiency?: number;
  generatorEfficiency?: number;
}

export function solveSteamTurbine(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: SteamTurbineSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): {
  outlet: Partial<MaterialStream>;
  shaftWork_W: number;
  electricPower_W: number;
  idealWork_W: number;
  exhaustT_K: number;
} {
  const outletP_Pa = Math.max(1000, spec.outletP_Pa);
  if (outletP_Pa >= inlet.P_Pa || inlet.totalFlow_molps <= 0) {
    return {
      outlet: { ...inlet, P_Pa: outletP_Pa, solved: true },
      shaftWork_W: 0,
      electricPower_W: 0,
      idealWork_W: 0,
      exhaustT_K: inlet.T_K,
    };
  }

  const eta = Math.max(0.01, Math.min(spec.efficiency ?? 0.75, 1));
  const generatorEfficiency = Math.max(0.01, Math.min(spec.generatorEfficiency ?? 0.96, 1));
  const z = [...inlet.moleFractions];
  const Cp_avg = compounds.reduce((sum, compound, i) => sum + z[i] * CpIG_JmolK(compound, inlet.T_K), 0);
  const Cv_avg = Math.max(Cp_avg - R, 1e-6);
  const kappa = Math.max(Cp_avg / Cv_avg, 1.01);
  const P_ratio = outletP_Pa / Math.max(inlet.P_Pa, 1);
  const T2s = Math.max(120, inlet.T_K * Math.pow(P_ratio, (kappa - 1) / kappa));
  const flash2s = flashPT(compounds, z, T2s, outletP_Pa, thermo.fluidPackage, thermo.interactionParams);
  const idealDrop_Jpmol = Math.max(0, inlet.H_Jpmol - flash2s.H_Jpmol);
  const shaftWork_W = inlet.totalFlow_molps * idealDrop_Jpmol * eta;
  const electricPower_W = shaftWork_W * generatorEfficiency;
  const flash = flashPQ(
    compounds,
    z,
    outletP_Pa,
    inlet.H_Jpmol,
    -shaftWork_W,
    inlet.totalFlow_molps,
    thermo.fluidPackage,
    thermo.interactionParams,
    T2s,
  );

  return {
    outlet: {
      T_K: flash.T_K,
      P_Pa: outletP_Pa,
      totalFlow_molps: inlet.totalFlow_molps,
      moleFractions: z,
      phase: flash.phase,
      vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    shaftWork_W,
    electricPower_W,
    idealWork_W: inlet.totalFlow_molps * idealDrop_Jpmol,
    exhaustT_K: flash.T_K,
  };
}

export interface SteamValveSpec {
  outletP_Pa: number;
}

export function solveSteamValve(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: SteamValveSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): { outlet: Partial<MaterialStream>; pressureDrop_Pa: number } {
  const result = solveValve(compounds, inlet, { outletP_Pa: spec.outletP_Pa }, thermo);
  return {
    outlet: result.outlet,
    pressureDrop_Pa: Math.max(inlet.P_Pa - spec.outletP_Pa, 0),
  };
}

export interface SteamHeaderSpec {
  outletP_Pa?: number;
}

export function solveSteamHeader(
  compounds: CanopyCompound[],
  inlets: MaterialStream[],
  spec: SteamHeaderSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): { outlet: Partial<MaterialStream>; duty_W: number } {
  const outletP_Pa = Math.max(
    spec.outletP_Pa ?? Math.min(...inlets.map(inlet => inlet.P_Pa)),
    1000,
  );
  return solveMixer(compounds, inlets, outletP_Pa, thermo);
}

export interface SteamTrapSpec {
  outletP_Pa: number;
  maxFlashFraction?: number;
}

export function solveSteamTrap(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: SteamTrapSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): {
  vaporOutlet: Partial<MaterialStream>;
  liquidOutlet: Partial<MaterialStream>;
  flashFraction: number;
  pressureDrop_Pa: number;
} {
  const z = [...inlet.moleFractions];
  const flash = flashPQ(
    compounds,
    z,
    spec.outletP_Pa,
    inlet.H_Jpmol,
    0,
    inlet.totalFlow_molps,
    thermo.fluidPackage,
    thermo.interactionParams,
    inlet.T_K,
  );
  const flashFraction = Math.max(0, Math.min(flash.vaporFraction, spec.maxFlashFraction ?? 1));
  const vaporFlow = inlet.totalFlow_molps * flashFraction;
  const liquidFlow = Math.max(0, inlet.totalFlow_molps - vaporFlow);

  return {
    vaporOutlet: {
      T_K: flash.T_K,
      P_Pa: spec.outletP_Pa,
      totalFlow_molps: vaporFlow,
      moleFractions: vaporFlow > 0 ? flash.y : z,
      phase: 'Vapor',
      vaporFraction: 1,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    liquidOutlet: {
      T_K: flash.T_K,
      P_Pa: spec.outletP_Pa,
      totalFlow_molps: liquidFlow,
      moleFractions: liquidFlow > 0 ? flash.x : z,
      phase: 'Liquid',
      vaporFraction: 0,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    flashFraction,
    pressureDrop_Pa: Math.max(0, inlet.P_Pa - spec.outletP_Pa),
  };
}

export interface PIDControllerSpec {
  setpoint: number;
  measuredProperty: 'T_K' | 'P_Pa' | 'totalFlow_molps' | 'vaporFraction' | 'signalValue' | string;
  gain: number;
  integralTime_s?: number;
  derivativeTime_s?: number;
  timestep_s?: number;
  measurementLag_s?: number;
  outputRateLimit_per_s?: number;
  bias?: number;
  lowLimit?: number;
  highLimit?: number;
  reverseActing?: boolean;
  previousError?: number;
  previousPreviousError?: number;
  previousMeasuredValue?: number;
  filteredMeasuredValue?: number;
  lastOutput?: number;
}

export function solvePIDControllerSignal(
  measuredValue: number,
  spec: PIDControllerSpec,
): {
  controllerOutput: number;
  measuredValue: number;
  filteredMeasuredValue: number;
  error: number;
  previousError: number;
  previousPreviousError: number;
  previousMeasuredValue: number;
} {
  const dt = Math.max(spec.timestep_s ?? 1, 1e-9);
  const lag = Math.max(spec.measurementLag_s ?? 0, 0);
  const prevFilteredMeasured = spec.filteredMeasuredValue ?? measuredValue;
  const filterAlpha = lag > 1e-9 ? Math.min(1, dt / (lag + dt)) : 1;
  const filteredMeasuredValue = prevFilteredMeasured + filterAlpha * (measuredValue - prevFilteredMeasured);
  const rawError = spec.setpoint - filteredMeasuredValue;
  const error = spec.reverseActing ? -rawError : rawError;
  const Ti = spec.integralTime_s;
  const Td = Math.max(spec.derivativeTime_s ?? 0, 0);
  const prevError = spec.previousError ?? 0;
  const prevPrevError = spec.previousPreviousError ?? prevError;
  const prevMeasured = spec.previousMeasuredValue ?? filteredMeasuredValue;
  const baseOutput = spec.lastOutput ?? spec.bias ?? 0;
  const proportionalIncrement = spec.gain * (error - prevError);
  const integralIncrement = Ti && isFinite(Ti) && Ti > 1e-9
    ? spec.gain * error * dt / Ti
    : 0;
  const derivativeIncrement = Td > 0
    ? -spec.gain * Td * (filteredMeasuredValue - 2 * prevMeasured + prevFilteredMeasured) / dt
    : 0;
  const unconstrained = baseOutput + proportionalIncrement + integralIncrement + derivativeIncrement;
  let controllerOutput = Math.max(
    spec.lowLimit ?? -Infinity,
    Math.min(spec.highLimit ?? Infinity, unconstrained),
  );
  const outputRateLimit = spec.outputRateLimit_per_s;
  if (outputRateLimit != null && isFinite(outputRateLimit) && outputRateLimit > 0) {
    const maxDelta = outputRateLimit * dt;
    controllerOutput = Math.max(baseOutput - maxDelta, Math.min(baseOutput + maxDelta, controllerOutput));
  }

  return {
    controllerOutput,
    measuredValue,
    filteredMeasuredValue,
    error,
    previousError: error,
    previousPreviousError: prevError,
    previousMeasuredValue: filteredMeasuredValue,
  };
}

export function solvePIDController(
  inlet: MaterialStream,
  spec: PIDControllerSpec,
): {
  outlet: MaterialStream;
  controllerOutput: number;
  measuredValue: number;
  filteredMeasuredValue: number;
  error: number;
  previousError: number;
  previousPreviousError: number;
  previousMeasuredValue: number;
} {
  const measuredValue = typeof (inlet as unknown as Record<string, unknown>)[spec.measuredProperty] === 'number'
    ? ((inlet as unknown as Record<string, number>)[spec.measuredProperty] as number)
    : inlet.T_K;
  const result = solvePIDControllerSignal(measuredValue, spec);
  return {
    outlet: { ...inlet, solved: true },
    ...result,
  };
}

export interface LeadLagSpec {
  leadTime_s?: number;
  lagTime_s?: number;
  timestep_s?: number;
  previousInput?: number;
  lastOutput?: number;
}

export function solveLeadLagSignal(
  inputValue: number,
  spec: LeadLagSpec,
): { signalValue: number; previousInput: number; lastOutput: number } {
  const dt = Math.max(spec.timestep_s ?? 1, 1e-9);
  const lead = Math.max(spec.leadTime_s ?? 0, 0);
  const lag = Math.max(spec.lagTime_s ?? 0, 0);
  const prevInput = spec.previousInput ?? inputValue;
  const prevOutput = spec.lastOutput ?? inputValue;
  const derivative = (inputValue - prevInput) / dt;
  const leadTarget = inputValue + lead * derivative;
  const alpha = lag > 1e-9 ? Math.min(1, dt / (lag + dt)) : 1;
  const signalValue = prevOutput + alpha * (leadTarget - prevOutput);
  return { signalValue, previousInput: inputValue, lastOutput: signalValue };
}

export interface DeadTimeSpec {
  deadTime_s?: number;
  timestep_s?: number;
  history?: number[];
}

export function solveDeadTimeSignal(
  inputValue: number,
  spec: DeadTimeSpec,
): { signalValue: number; history: number[] } {
  const dt = Math.max(spec.timestep_s ?? 1, 1e-9);
  const deadTime = Math.max(spec.deadTime_s ?? 0, 0);
  const steps = Math.max(1, Math.round(deadTime / dt));
  const history = Array.isArray(spec.history) ? [...spec.history] : [];
  history.push(inputValue);
  while (history.length > steps + 1) history.shift();
  const signalValue = history.length > steps ? history[0] : inputValue;
  return { signalValue, history };
}

export interface SignalSelectorSpec {
  mode: 'High' | 'Low' | 'Average';
}

export function solveSignalSelector(
  inputValues: number[],
  spec: SignalSelectorSpec,
): { signalValue: number } {
  const values = inputValues.filter(value => Number.isFinite(value));
  if (values.length === 0) return { signalValue: 0 };
  if (spec.mode === 'Low') return { signalValue: Math.min(...values) };
  if (spec.mode === 'Average') return { signalValue: values.reduce((sum, value) => sum + value, 0) / values.length };
  return { signalValue: Math.max(...values) };
}

export interface AnalyzerSpec {
  propertyCode: string;
  componentIndex?: number;
}

export interface AnalyzerResult {
  signalValue: number;
  propertyLabel: string;
  warnings: string[];
}

function normalizeAnalyzerCode(propertyCode: string): string {
  return propertyCode.trim().toUpperCase().replace(/[^A-Z0-9]+/g, '_');
}

function streamMassFractions(compounds: CanopyCompound[], moleFractions: number[]): number[] {
  const numerators = moleFractions.map((z, index) => z * (compounds[index]?.molecularWeight ?? 0));
  const total = numerators.reduce((sum, value) => sum + value, 0);
  if (total <= 0) return moleFractions.map(() => 0);
  return numerators.map(value => value / total);
}

function streamMixtureMW(compounds: CanopyCompound[], moleFractions: number[]): number {
  return moleFractions.reduce((sum, z, index) => sum + z * (compounds[index]?.molecularWeight ?? 0), 0);
}

function streamSpecificGravity(compounds: CanopyCompound[], moleFractions: number[]): number | null {
  const massFractions = streamMassFractions(compounds, moleFractions);
  let reciprocal = 0;
  for (let i = 0; i < compounds.length; i++) {
    const sg = compounds[i]?.specificGravity;
    if (sg == null || sg <= 0) return null;
    reciprocal += massFractions[i] / sg;
  }
  return reciprocal > 0 ? 1 / reciprocal : null;
}

function weightedCompoundProperty(
  compounds: CanopyCompound[],
  moleFractions: number[],
  getter: (compound: CanopyCompound) => number | undefined,
): number | null {
  let totalWeight = 0;
  let weighted = 0;
  for (let i = 0; i < compounds.length; i++) {
    const value = getter(compounds[i]);
    if (value == null || !Number.isFinite(value)) continue;
    const weight = moleFractions[i] ?? 0;
    weighted += weight * value;
    totalWeight += weight;
  }
  return totalWeight > 0 ? weighted / totalWeight : null;
}

function weightedCompoundPropertyByVolume(
  compounds: CanopyCompound[],
  moleFractions: number[],
  getter: (compound: CanopyCompound) => number | undefined,
): number | null {
  let totalWeight = 0;
  let weighted = 0;
  for (let i = 0; i < compounds.length; i++) {
    const value = getter(compounds[i]);
    if (value == null || !Number.isFinite(value)) continue;
    const sg = compounds[i]?.specificGravity;
    const mw = compounds[i]?.molecularWeight ?? 0;
    const pseudoVolumeWeight = (moleFractions[i] ?? 0) * mw / Math.max(sg ?? 1, 1e-6);
    if (pseudoVolumeWeight <= 0) continue;
    weighted += pseudoVolumeWeight * value;
    totalWeight += pseudoVolumeWeight;
  }
  return totalWeight > 0 ? weighted / totalWeight : null;
}

function estimatedDistillationTemperature(
  compounds: CanopyCompound[],
  moleFractions: number[],
  recoveryFraction: number,
): number | null {
  const ranked = compounds
    .map((compound, index) => ({ Tb_K: compound.Tb_K, fraction: moleFractions[index] ?? 0 }))
    .filter(item => item.Tb_K != null && Number.isFinite(item.Tb_K) && item.fraction > 0)
    .sort((a, b) => (a.Tb_K ?? 0) - (b.Tb_K ?? 0));

  const total = ranked.reduce((sum, item) => sum + item.fraction, 0);
  if (total <= 0) return null;

  const target = Math.min(0.999, Math.max(0.001, recoveryFraction)) * total;
  let cumulative = 0;
  for (let i = 0; i < ranked.length; i++) {
    const previous = cumulative;
    cumulative += ranked[i].fraction;
    if (cumulative >= target) {
      if (i === 0) return ranked[i].Tb_K ?? null;
      const prevTb = ranked[i - 1].Tb_K ?? ranked[i].Tb_K ?? null;
      const currTb = ranked[i].Tb_K ?? prevTb;
      if (prevTb == null || currTb == null) return currTb;
      const span = Math.max(cumulative - previous, 1e-12);
      const blend = (target - previous) / span;
      return prevTb + (currTb - prevTb) * blend;
    }
  }
  return ranked[ranked.length - 1]?.Tb_K ?? null;
}

function streamWaterSaturationFraction(
  compounds: CanopyCompound[],
  stream: MaterialStream,
): number | null {
  const waterIndex = compounds.findIndex(compound => {
    const name = compound.name.trim().toUpperCase();
    return name === 'WATER' || name === 'H2O';
  });
  if (waterIndex < 0) return null;

  const componentFlows = stream.moleFractions.map(fraction => fraction * stream.totalFlow_molps);
  const waterFlow = componentFlows[waterIndex] ?? 0;
  if (waterFlow <= 0) return 0;

  let waterCapacity = 0;
  for (let i = 0; i < compounds.length; i++) {
    if (i === waterIndex) continue;
    const compound = compounds[i];
    const solubility = waterSolubilityInHC(compound, stream.T_K);
    const flow = componentFlows[i] ?? 0;
    if (!Number.isFinite(solubility) || solubility <= 0 || solubility >= 0.5 || flow <= 0) continue;
    waterCapacity += flow * solubility / Math.max(1 - solubility, 1e-9);
  }

  if (waterCapacity <= 1e-12) {
    return waterFlow > 0 ? 1 : 0;
  }
  return waterFlow / waterCapacity;
}

export function solveAnalyzer(
  compounds: CanopyCompound[],
  stream: MaterialStream,
  spec: AnalyzerSpec,
): AnalyzerResult {
  const code = normalizeAnalyzerCode(spec.propertyCode || 'TEMP');
  const componentIndex = Math.max(0, Math.floor(spec.componentIndex ?? 0));
  const massFractions = streamMassFractions(compounds, stream.moleFractions);
  const specificGravity = streamSpecificGravity(compounds, stream.moleFractions);

  switch (code) {
    case 'TEMP':
    case 'TEMPERATURE':
      return { signalValue: stream.T_K, propertyLabel: 'Temperature (K)', warnings: [] };
    case 'PRES':
    case 'PRESSURE':
      return { signalValue: stream.P_Pa, propertyLabel: 'Pressure (Pa)', warnings: [] };
    case 'ENTH':
    case 'MOLE_ENTH':
      return { signalValue: stream.H_Jpmol, propertyLabel: 'Molar Enthalpy (J/mol)', warnings: [] };
    case 'MOLE_FLOW':
    case 'MOLEFLOW':
      return { signalValue: stream.totalFlow_molps, propertyLabel: 'Molar Flow (mol/s)', warnings: [] };
    case 'MASS_FLOW':
    case 'MASSFLOW': {
      const massFlow_kgps = stream.totalFlow_molps * streamMixtureMW(compounds, stream.moleFractions) / 1000;
      return { signalValue: massFlow_kgps, propertyLabel: 'Mass Flow (kg/s)', warnings: [] };
    }
    case 'VFRA':
    case 'VAPOR_FRACTION':
      return { signalValue: stream.vaporFraction, propertyLabel: 'Vapor Fraction', warnings: [] };
    case 'LFRA':
    case 'LIQUID_FRACTION':
      return { signalValue: 1 - stream.vaporFraction, propertyLabel: 'Liquid Fraction', warnings: [] };
    case 'MOLE_FRAC':
      return {
        signalValue: stream.moleFractions[componentIndex] ?? 0,
        propertyLabel: `Mole Fraction x[${componentIndex}]`,
        warnings: [],
      };
    case 'MASS_FRAC':
      return {
        signalValue: massFractions[componentIndex] ?? 0,
        propertyLabel: `Mass Fraction w[${componentIndex}]`,
        warnings: [],
      };
    case 'TBPT':
    case 'TBPWT': {
      const tbp = weightedCompoundProperty(compounds, stream.moleFractions, compound => compound.Tb_K);
      return tbp != null
        ? { signalValue: tbp, propertyLabel: 'Average TBP (K)', warnings: [] }
        : { signalValue: 0, propertyLabel: 'Average TBP (K)', warnings: ['TBP could not be estimated because the stream compounds do not all have normal boiling points.'] };
    }
    case 'CETA': {
      const cetane = weightedCompoundProperty(compounds, stream.moleFractions, compound => compound.cetaneNumber);
      return cetane != null
        ? {
            signalValue: cetane,
            propertyLabel: 'Cetane Number Approx.',
            warnings: ['CETA is estimated from component cetane-number metadata and is not a fitted ASTM engine-test prediction for the mixture.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'Cetane Number Approx.',
            warnings: ['CETA could not be estimated because compound cetane-number metadata is unavailable for this stream.'],
          };
    }
    case 'AIT':
    case 'AUTOIGNITION':
    case 'AUTOIGNITION_TEMP': {
      const autoignition = weightedCompoundProperty(compounds, stream.moleFractions, compound => compound.autoignitionTemp_K);
      return autoignition != null
        ? {
            signalValue: autoignition,
            propertyLabel: 'Autoignition Temperature Approx. (K)',
            warnings: ['AIT is estimated from component autoignition metadata and is not a full ASTM E659-style mixture prediction.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'Autoignition Temperature Approx. (K)',
            warnings: ['AIT could not be estimated because compound autoignition metadata is unavailable.'],
          };
    }
    case 'VABP': {
      const vabp = weightedCompoundPropertyByVolume(compounds, stream.moleFractions, compound => compound.Tb_K);
      return vabp != null
        ? {
            signalValue: vabp,
            propertyLabel: 'Volume Average Boiling Point Approx. (K)',
            warnings: ['VABP is estimated from normal boiling points weighted by pseudo-liquid volume from molecular weight and specific gravity.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'Volume Average Boiling Point Approx. (K)',
            warnings: ['VABP could not be estimated because boiling-point or specific-gravity data is missing.'],
          };
    }
    case 'D86T': {
      const d86 = estimatedDistillationTemperature(compounds, stream.moleFractions, 0.5);
      return d86 != null
        ? {
            signalValue: d86,
            propertyLabel: 'ASTM D86 Equivalent Temperature (K)',
            warnings: ['D86 is estimated from the stream boiling-point distribution, not from a full ASTM D86 curve reconstruction.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'ASTM D86 Equivalent Temperature (K)',
            warnings: ['D86 could not be estimated because the stream compounds are missing normal boiling points.'],
          };
    }
    case 'D1160T': {
      const d1160 = estimatedDistillationTemperature(compounds, stream.moleFractions, 0.7);
      return d1160 != null
        ? {
            signalValue: d1160,
            propertyLabel: 'ASTM D1160 Equivalent Temperature (K)',
            warnings: ['D1160 is estimated from the stream boiling-point distribution, not from a full vacuum-distillation reconstruction.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'ASTM D1160 Equivalent Temperature (K)',
            warnings: ['D1160 could not be estimated because the stream compounds are missing normal boiling points.'],
          };
    }
    case 'D2887T': {
      const d2887 = estimatedDistillationTemperature(compounds, stream.moleFractions, 0.9);
      return d2887 != null
        ? {
            signalValue: d2887,
            propertyLabel: 'ASTM D2887 Equivalent Temperature (K)',
            warnings: ['D2887 is estimated from the stream boiling-point distribution, not from a full simulated distillation calibration.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'ASTM D2887 Equivalent Temperature (K)',
            warnings: ['D2887 could not be estimated because the stream compounds are missing normal boiling points.'],
          };
    }
    case 'BUBPT': {
      const Tbub = bubblePointT(compounds, stream.moleFractions, stream.P_Pa);
      return Tbub != null
        ? { signalValue: Tbub, propertyLabel: 'Bubble Point Temperature (K)', warnings: [] }
        : { signalValue: 0, propertyLabel: 'Bubble Point Temperature (K)', warnings: ['Bubble point could not be estimated for this stream.'] };
    }
    case 'DEWPT': {
      const Tdew = dewPointT(compounds, stream.moleFractions, stream.P_Pa);
      return Tdew != null
        ? { signalValue: Tdew, propertyLabel: 'Dew Point Temperature (K)', warnings: [] }
        : { signalValue: 0, propertyLabel: 'Dew Point Temperature (K)', warnings: ['Dew point could not be estimated for this stream.'] };
    }
    case 'WATBUB':
    case 'WAT_BUBPT': {
      const Tbub = bubblePointTWaterSaturated(compounds, stream.moleFractions, stream.P_Pa);
      return Tbub != null
        ? {
            signalValue: Tbub,
            propertyLabel: 'Water-Saturated Bubble Point (K)',
            warnings: ['Water-saturated bubble point is estimated by enriching the stream to hydrocarbon water solubility at the solved temperature.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'Water-Saturated Bubble Point (K)',
            warnings: ['Water-saturated bubble point could not be estimated because no water-saturation path was available.'],
          };
    }
    case 'WATDEW':
    case 'WAT_DEWPT': {
      const Tdew = dewPointTWaterSaturated(compounds, stream.moleFractions, stream.P_Pa);
      return Tdew != null
        ? {
            signalValue: Tdew,
            propertyLabel: 'Water-Saturated Dew Point (K)',
            warnings: ['Water-saturated dew point is estimated by enriching the stream to hydrocarbon water solubility at the solved temperature.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'Water-Saturated Dew Point (K)',
            warnings: ['Water-saturated dew point could not be estimated because no water-saturation path was available.'],
          };
    }
    case 'SG':
    case 'SPECIFIC_GRAVITY':
      return specificGravity != null
        ? { signalValue: specificGravity, propertyLabel: 'Specific Gravity', warnings: [] }
        : { signalValue: 0, propertyLabel: 'Specific Gravity', warnings: ['Specific gravity could not be estimated because one or more compounds are missing standard-density metadata.'] };
    case 'API':
      return specificGravity != null
        ? { signalValue: 141.5 / specificGravity - 131.5, propertyLabel: 'API Gravity', warnings: [] }
        : { signalValue: 0, propertyLabel: 'API Gravity', warnings: ['API gravity could not be estimated because specific gravity is unavailable for this stream.'] };
    case 'RI':
    case 'REFI':
    case 'NDEX':
    case 'REFINDEX':
    case 'REFRACTIVE_INDEX': {
      const refractiveIndex = weightedCompoundProperty(compounds, stream.moleFractions, compound => compound.refractiveIndex);
      return refractiveIndex != null
        ? { signalValue: refractiveIndex, propertyLabel: 'Refractive Index', warnings: [] }
        : { signalValue: 0, propertyLabel: 'Refractive Index', warnings: ['Refractive index could not be estimated because one or more compounds are missing refractive-index metadata.'] };
    }
    case 'WAT': {
      const waterIndex = compounds.findIndex(compound => compound.name.toUpperCase() === 'WATER');
      const waterMassPct = waterIndex >= 0 ? (massFractions[waterIndex] ?? 0) * 100 : 0;
      return {
        signalValue: waterMassPct,
        propertyLabel: 'Water Content (mass %)',
        warnings: waterIndex >= 0 ? [] : ['No WATER component is present, so water content reports as zero.'],
      };
    }
    case 'WATSAT':
    case 'WATER_SATURATION': {
      const waterSaturation = streamWaterSaturationFraction(compounds, stream);
      return waterSaturation != null
        ? {
            signalValue: waterSaturation * 100,
            propertyLabel: 'Water Saturation (%)',
            warnings: ['Water saturation is estimated from APV140-style water-in-hydrocarbon solubility correlations where available.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'Water Saturation (%)',
            warnings: ['Water saturation could not be estimated because the stream has no WATER component metadata.'],
          };
    }
    case 'GRS':
    case 'QVAL':
    case 'HHV': {
      const grossHeatingValue = weightedCompoundProperty(compounds, stream.moleFractions, compound =>
        compound.Hcomb_Jmol != null
          ? Math.abs(compound.Hcomb_Jmol)
          : compound.hCombustion_Jpkmol != null
            ? Math.abs(compound.hCombustion_Jpkmol) / 1000
            : undefined,
      );
      return grossHeatingValue != null
        ? {
            signalValue: grossHeatingValue,
            propertyLabel: 'Gross Heating Value Approx. (J/mol)',
            warnings: ['GRS/QVAL is estimated from component combustion metadata and is not a bomb-calorimetry mixture reconstruction.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'Gross Heating Value Approx. (J/mol)',
            warnings: ['GRS/QVAL could not be estimated because combustion-heating metadata is unavailable.'],
          };
    }
    case 'NET':
    case 'LHV': {
      const netHeatingValue = weightedCompoundProperty(compounds, stream.moleFractions, compound =>
        compound.hCombustion_Jpkmol != null
          ? Math.abs(compound.hCombustion_Jpkmol) / 1000
          : compound.Hcomb_Jmol != null
            ? Math.abs(compound.Hcomb_Jmol)
            : undefined,
      );
      return netHeatingValue != null
        ? {
            signalValue: netHeatingValue,
            propertyLabel: 'Net Heating Value Approx. (J/mol)',
            warnings: ['NET is estimated from component lower-heating-value metadata where available, otherwise it falls back to gross combustion values.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'Net Heating Value Approx. (J/mol)',
            warnings: ['NET could not be estimated because combustion-heating metadata is unavailable.'],
          };
    }
    case 'FLASHPOINT': {
      const flashPoint = weightedCompoundProperty(compounds, stream.moleFractions, compound => compound.flashPoint_K);
      return flashPoint != null
        ? {
            signalValue: flashPoint,
            propertyLabel: 'Flash Point Approx. (K)',
            warnings: ['Flash point is estimated from compound flash-point metadata and is not a rigorous closed-cup mixture prediction.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'Flash Point Approx. (K)',
            warnings: ['Flash point could not be estimated because compound flash-point metadata is unavailable.'],
          };
    }
    case 'ANILPT':
    case 'ANILINEPT':
    case 'ANILINE_POINT': {
      const anilinePoint = weightedCompoundProperty(compounds, stream.moleFractions, compound => compound.anilinePoint_K);
      return anilinePoint != null
        ? {
            signalValue: anilinePoint,
            propertyLabel: 'Aniline Point Approx. (K)',
            warnings: ['Aniline point is estimated from compound petroleum-characterization metadata and is not a full mixture correlation.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'Aniline Point Approx. (K)',
            warnings: ['Aniline point could not be estimated because compound petroleum-characterization metadata is unavailable.'],
          };
    }
    case 'REID':
    case 'REIDVP': {
      const T_rvp = 310.95;
      let totalPressure = 0;
      let contributingFraction = 0;
      for (let i = 0; i < compounds.length; i++) {
        const name = compounds[i].name.trim().toUpperCase();
        if (name === 'NITROGEN' || name === 'N2' || name === 'OXYGEN' || name === 'O2') continue;
        const psat = Psat_Pa(compounds[i], T_rvp);
        if (!Number.isFinite(psat) || psat <= 0) continue;
        const z = stream.moleFractions[i] ?? 0;
        totalPressure += z * psat;
        contributingFraction += z;
      }
      return contributingFraction > 0
        ? {
            signalValue: totalPressure,
            propertyLabel: 'Reid Vapor Pressure Approx. (Pa)',
            warnings: ['REIDVP is approximated from component vapor pressures at 37.8°C and does not reproduce the full ASTM D323 procedure.'],
          }
        : {
            signalValue: 0,
            propertyLabel: 'Reid Vapor Pressure Approx. (Pa)',
            warnings: ['REIDVP could not be estimated because no valid vapor-pressure data was available at 37.8°C.'],
          };
    }
    case 'SMX': {
      const entropy = streamEntropy(
        compounds,
        stream.T_K,
        stream.P_Pa,
        stream.vaporFraction,
        stream.x_liquid,
        stream.y_vapor,
      );
      return {
        signalValue: entropy,
        propertyLabel: 'Mixture Entropy (J/mol-K)',
        warnings: [],
      };
    }
    default:
      return {
        signalValue: 0,
        propertyLabel: spec.propertyCode || 'Analyzer',
        warnings: [`Analyzer property "${spec.propertyCode}" is not implemented yet.`],
      };
  }
}

export interface ChargeBalanceSpec {
  adjustComponentIndex: number;
  targetCharge?: number;
}

export interface ChargeBalanceResult {
  outlet: Partial<MaterialStream>;
  residualCharge: number;
  adjustedComponentFlow_molps: number;
}

export function solveChargeBalance(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: ChargeBalanceSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): ChargeBalanceResult {
  const charges = compounds.map(compound => compound.chargeNumber ?? 0);
  const adjustIndex = Math.max(0, Math.min(compounds.length - 1, Math.floor(spec.adjustComponentIndex)));
  const adjustCharge = charges[adjustIndex] ?? 0;
  if (adjustCharge === 0) {
    throw new Error(`Selected adjustment species "${compounds[adjustIndex]?.displayName ?? adjustIndex}" is neutral.`);
  }

  const componentFlows = inlet.moleFractions.map(fraction => inlet.totalFlow_molps * fraction);
  const currentCharge = componentFlows.reduce((sum, flow, index) => sum + flow * charges[index], 0);
  const targetCharge = spec.targetCharge ?? 0;
  const deltaFlow = (targetCharge - currentCharge) / adjustCharge;
  const newAdjustFlow = (componentFlows[adjustIndex] ?? 0) + deltaFlow;
  if (newAdjustFlow < -1e-12) {
    throw new Error('Charge balance would require a negative flow for the selected adjustment species.');
  }

  componentFlows[adjustIndex] = Math.max(0, newAdjustFlow);
  const totalFlow = componentFlows.reduce((sum, flow) => sum + flow, 0);
  if (totalFlow <= 0) {
    throw new Error('Charge-balanced stream flow became non-positive.');
  }

  const moleFractions = componentFlows.map(flow => flow / totalFlow);
  const residualCharge = componentFlows.reduce((sum, flow, index) => sum + flow * charges[index], 0) - targetCharge;
  const flash = flashPT(
    compounds,
    moleFractions,
    inlet.T_K,
    inlet.P_Pa,
    thermo.fluidPackage,
    thermo.interactionParams,
  );

  return {
    outlet: {
      T_K: inlet.T_K,
      P_Pa: inlet.P_Pa,
      totalFlow_molps: totalFlow,
      moleFractions,
      phase: flash.phase,
      vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    residualCharge,
    adjustedComponentFlow_molps: componentFlows[adjustIndex],
  };
}

export interface ElectrolyteEquilibriumSpec {
  apparentSpeciesIndex: number;
  cationIndex: number;
  anionIndex: number;
  cationStoich?: number;
  anionStoich?: number;
  log10K?: number;
  maxDissociationFraction?: number;
  reactionSet?: ElectrolyteReactionEntry[];
  maxPasses?: number;
}

export interface ElectrolyteReactionEntry {
  name?: string;
  active?: boolean;
  apparentSpeciesIndex?: number;
  cationIndex?: number;
  anionIndex?: number;
  cationStoich?: number;
  anionStoich?: number;
  log10K?: number;
  maxDissociationFraction?: number;
  stoichiometry?: number[];
  referenceT_K?: number;
  deltaH_Jmol?: number;
  deltaS_JmolK?: number;
}

export interface ElectrolyteEquilibriumResult {
  outlet: Partial<MaterialStream>;
  dissociationFraction: number;
  ionicStrength: number;
  equilibriumConstant: number;
  activeReactionCount: number;
  equilibriumPasses: number;
  reactionExtents_molps: number[];
}

function electrolyteReactionStoichiometry(
  compounds: CanopyCompound[],
  reaction: ElectrolyteReactionEntry,
): number[] {
  if (Array.isArray(reaction.stoichiometry) && reaction.stoichiometry.length > 0) {
    const stoichiometry = new Array(compounds.length).fill(0);
    for (let i = 0; i < Math.min(compounds.length, reaction.stoichiometry.length); i++) {
      stoichiometry[i] = Number(reaction.stoichiometry[i] ?? 0);
    }
    if (stoichiometry.some(value => value < 0) && stoichiometry.some(value => value > 0)) {
      return stoichiometry;
    }
  }

  const apparentSpeciesIndex = Math.max(0, Math.min(compounds.length - 1, Math.floor(reaction.apparentSpeciesIndex ?? 0)));
  const cationIndex = Math.max(0, Math.min(compounds.length - 1, Math.floor(reaction.cationIndex ?? 0)));
  const anionIndex = Math.max(0, Math.min(compounds.length - 1, Math.floor(reaction.anionIndex ?? 0)));
  const stoichiometry = new Array(compounds.length).fill(0);
  stoichiometry[apparentSpeciesIndex] -= 1;
  stoichiometry[cationIndex] += Math.max(1, reaction.cationStoich ?? 1);
  stoichiometry[anionIndex] += Math.max(1, reaction.anionStoich ?? 1);
  return stoichiometry;
}

function electrolyteReactionExtentBounds(baseFlows: number[], stoichiometry: number[]): { lo: number; hi: number } {
  let lo = Number.NEGATIVE_INFINITY;
  let hi = Number.POSITIVE_INFINITY;
  for (let i = 0; i < stoichiometry.length; i++) {
    const nu = stoichiometry[i] ?? 0;
    if (Math.abs(nu) < 1e-14) continue;
    const flow = Math.max(baseFlows[i] ?? 0, 0);
    if (nu > 0) {
      lo = Math.max(lo, -flow / nu);
    } else {
      hi = Math.min(hi, flow / -nu);
    }
  }
  if (!Number.isFinite(lo)) lo = 0;
  if (!Number.isFinite(hi)) hi = 0;
  return { lo, hi };
}

function electrolyteReactionEquilibriumConstant(
  reaction: ElectrolyteReactionEntry,
  T_K: number,
): number {
  if (Number.isFinite(reaction.deltaH_Jmol) && Number.isFinite(reaction.deltaS_JmolK)) {
    return Math.exp((reaction.deltaS_JmolK as number) / R - (reaction.deltaH_Jmol as number) / (R * Math.max(T_K, 1)));
  }
  const baseK = Math.pow(10, reaction.log10K ?? 2);
  if (Number.isFinite(reaction.deltaH_Jmol) && Number.isFinite(reaction.referenceT_K) && (reaction.referenceT_K as number) > 0) {
    const lnK = Math.log(Math.max(baseK, 1e-30))
      - (reaction.deltaH_Jmol as number) / R * (1 / Math.max(T_K, 1) - 1 / (reaction.referenceT_K as number));
    return Math.exp(Math.max(-80, Math.min(80, lnK)));
  }
  return baseK;
}

function electrolyteReactionResidual(
  compounds: CanopyCompound[],
  baseFlows: number[],
  extent: number,
  stoichiometry: number[],
  T_K: number,
  thermo: ThermoOptions,
  K_eq: number,
): number {
  const flows = baseFlows.map((flow, index) => Math.max(0, flow + (stoichiometry[index] ?? 0) * extent));
  const total = flows.reduce((sum, value) => sum + value, 0);
  if (total <= 0) return Number.POSITIVE_INFINITY;
  const x = flows.map(value => value / total);
  const gamma = computeGamma(compounds, x, T_K, thermo.fluidPackage, thermo.interactionParams);

  let lnQ = 0;
  for (let i = 0; i < stoichiometry.length; i++) {
    const nu = stoichiometry[i] ?? 0;
    if (Math.abs(nu) < 1e-14) continue;
    const activity = Math.max(1e-30, (x[i] ?? 0) * (gamma[i] ?? 1));
    lnQ += nu * Math.log(activity);
  }
  return lnQ - Math.log(Math.max(K_eq, 1e-30));
}

function materializeStream(stream: Partial<MaterialStream>, fallback: MaterialStream): MaterialStream {
  return {
    id: stream.id ?? fallback.id,
    name: stream.name ?? fallback.name,
    T_K: stream.T_K ?? fallback.T_K,
    P_Pa: stream.P_Pa ?? fallback.P_Pa,
    totalFlow_molps: stream.totalFlow_molps ?? fallback.totalFlow_molps,
    moleFractions: stream.moleFractions ?? fallback.moleFractions,
    phase: stream.phase ?? fallback.phase,
    vaporFraction: stream.vaporFraction ?? fallback.vaporFraction,
    H_Jpmol: stream.H_Jpmol ?? fallback.H_Jpmol,
    x_liquid: stream.x_liquid ?? fallback.x_liquid,
    y_vapor: stream.y_vapor ?? fallback.y_vapor,
    sourceUnitId: stream.sourceUnitId ?? fallback.sourceUnitId,
    sourcePortId: stream.sourcePortId ?? fallback.sourcePortId,
    targetUnitId: stream.targetUnitId ?? fallback.targetUnitId,
    targetPortId: stream.targetPortId ?? fallback.targetPortId,
    solved: stream.solved ?? fallback.solved,
  };
}

function solveSingleElectrolyteReaction(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  reaction: ElectrolyteReactionEntry,
  thermo: ThermoOptions,
): { outlet: MaterialStream; dissociationFraction: number; equilibriumConstant: number; extent_molps: number } {
  const stoichiometry = electrolyteReactionStoichiometry(compounds, reaction);
  const K_eq = electrolyteReactionEquilibriumConstant(reaction, inlet.T_K);

  const baseFlows = inlet.moleFractions.map(fraction => inlet.totalFlow_molps * fraction);
  const { lo: extentLoBound, hi: extentHiBound } = electrolyteReactionExtentBounds(baseFlows, stoichiometry);
  let lo = extentLoBound;
  let hi = extentHiBound;
  if (reaction.maxDissociationFraction != null) {
    const apparentIndex = reaction.apparentSpeciesIndex != null ? Math.floor(reaction.apparentSpeciesIndex) : null;
    if (apparentIndex != null && apparentIndex >= 0 && apparentIndex < compounds.length) {
      const apparentFeed = baseFlows[apparentIndex] ?? 0;
      hi = Math.min(hi, apparentFeed * Math.max(0, Math.min(0.999999, reaction.maxDissociationFraction)));
    }
  }
  if (!(hi > lo)) {
    return {
      outlet: inlet,
      dissociationFraction: 0,
      equilibriumConstant: K_eq,
      extent_molps: 0,
    };
  }

  let extent = Math.max(lo, Math.min(0, hi));
  let fLo = electrolyteReactionResidual(compounds, baseFlows, lo, stoichiometry, inlet.T_K, thermo, K_eq);
  let fHi = electrolyteReactionResidual(compounds, baseFlows, hi, stoichiometry, inlet.T_K, thermo, K_eq);
  if (Number.isFinite(fLo) && Number.isFinite(fHi) && fLo * fHi < 0) {
    for (let iter = 0; iter < 80; iter++) {
      const mid = 0.5 * (lo + hi);
      const fMid = electrolyteReactionResidual(compounds, baseFlows, mid, stoichiometry, inlet.T_K, thermo, K_eq);
      extent = mid;
      if (Math.abs(fMid) < 1e-10 || Math.abs(hi - lo) < 1e-10) break;
      if (fLo * fMid <= 0) {
        hi = mid;
        fHi = fMid;
      } else {
        lo = mid;
        fLo = fMid;
      }
    }
  } else {
    extent = Math.abs(fHi) < Math.abs(fLo) ? hi : lo;
  }

  const outletFlows = baseFlows.map((flow, index) => Math.max(0, flow + (stoichiometry[index] ?? 0) * extent));
  const totalFlow = outletFlows.reduce((sum, value) => sum + value, 0);
  const moleFractions = outletFlows.map(value => value / Math.max(totalFlow, 1e-30));
  const flash = flashPT(compounds, moleFractions, inlet.T_K, inlet.P_Pa, thermo.fluidPackage, thermo.interactionParams);
  const reactantReference = stoichiometry.reduce((sum, nu, index) => {
    if (nu >= 0) return sum;
    return sum + (-nu) * Math.max(baseFlows[index] ?? 0, 0);
  }, 0);

  return {
    outlet: materializeStream({
      ...inlet,
      totalFlow_molps: totalFlow,
      moleFractions,
      phase: flash.phase,
      vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    }, inlet),
    dissociationFraction: Math.abs(extent) / Math.max(reactantReference, 1e-30),
    equilibriumConstant: K_eq,
    extent_molps: extent,
  };
}

export function solveElectrolyteEquilibrium(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  spec: ElectrolyteEquilibriumSpec,
  thermo: ThermoOptions = DEFAULT_THERMO,
): ElectrolyteEquilibriumResult {
  const configuredReactionSet = (spec.reactionSet ?? []).filter(reaction => reaction.active !== false);
  const reactionSet = configuredReactionSet.length > 0
    ? configuredReactionSet
    : [{
        apparentSpeciesIndex: spec.apparentSpeciesIndex,
        cationIndex: spec.cationIndex,
        anionIndex: spec.anionIndex,
        cationStoich: spec.cationStoich,
        anionStoich: spec.anionStoich,
        log10K: spec.log10K,
        maxDissociationFraction: spec.maxDissociationFraction,
      }];

  const reactionStoichiometries = reactionSet.map(reaction => electrolyteReactionStoichiometry(compounds, reaction));
  const uniqueApparentSpecies = [...new Set(reactionSet
    .map(reaction => reaction.apparentSpeciesIndex)
    .filter((index): index is number => Number.isFinite(index)))];
  const initialFlows = inlet.moleFractions.map(fraction => inlet.totalFlow_molps * fraction);
  if (uniqueApparentSpecies.length > 0 && uniqueApparentSpecies.every(index => (initialFlows[index] ?? 0) <= 0)) {
    const label = compounds[uniqueApparentSpecies[0] ?? 0]?.displayName ?? String(uniqueApparentSpecies[0] ?? 0);
    throw new Error(`Electrolyte equilibrium requires positive flow of "${label}".`);
  }

  let currentStream = inlet;
  let lastK = Math.pow(10, spec.log10K ?? 2);
  let passesUsed = 0;
  const cumulativeExtents = reactionSet.map(() => 0);
  const maxPasses = Math.max(1, Math.min(12, Math.floor(spec.maxPasses ?? Math.max(3, reactionSet.length * 2))));
  for (let pass = 0; pass < maxPasses; pass++) {
    const before = [...currentStream.moleFractions];
    for (let reactionIndex = 0; reactionIndex < reactionSet.length; reactionIndex++) {
      const reaction = reactionSet[reactionIndex];
      const step = solveSingleElectrolyteReaction(compounds, currentStream, reaction, thermo);
      currentStream = step.outlet;
      lastK = step.equilibriumConstant;
      cumulativeExtents[reactionIndex] += step.extent_molps;
    }
    passesUsed = pass + 1;
    const maxDelta = currentStream.moleFractions.reduce((acc, value, index) => Math.max(acc, Math.abs(value - (before[index] ?? 0))), 0);
    if (maxDelta < 1e-8) break;
  }

  const finalFlows = currentStream.moleFractions.map(fraction => currentStream.totalFlow_molps * fraction);
  const totalApparentFeed = uniqueApparentSpecies.reduce((sum, index) => sum + Math.max(initialFlows[index] ?? 0, 0), 0);
  const totalApparentRemaining = uniqueApparentSpecies.reduce((sum, index) => sum + Math.max(finalFlows[index] ?? 0, 0), 0);
  const totalReferenceReactants = reactionStoichiometries.reduce((sum, stoichiometry) => {
    return sum + stoichiometry.reduce((inner, nu, index) => {
      if (nu >= 0) return inner;
      return inner + (-nu) * Math.max(initialFlows[index] ?? 0, 0);
    }, 0);
  }, 0);
  const totalReacted = reactionStoichiometries.reduce((sum, stoichiometry, index) => {
    const reactantMultiplier = stoichiometry.reduce((inner, nu) => inner + (nu < 0 ? -nu : 0), 0);
    return sum + Math.abs(cumulativeExtents[index] ?? 0) * reactantMultiplier;
  }, 0);
  const ionicStrength = 0.5 * currentStream.moleFractions.reduce((sum, value, index) => {
    const charge = compounds[index]?.chargeNumber ?? 0;
    return sum + value * charge * charge;
  }, 0);

  return {
    outlet: currentStream,
    dissociationFraction: uniqueApparentSpecies.length > 0
      ? Math.max(0, totalApparentFeed - totalApparentRemaining) / Math.max(totalApparentFeed, 1e-30)
      : totalReacted / Math.max(totalReferenceReactants, 1e-30),
    ionicStrength,
    equilibriumConstant: lastK,
    activeReactionCount: reactionSet.length,
    equilibriumPasses: passesUsed,
    reactionExtents_molps: cumulativeExtents,
  };
}
