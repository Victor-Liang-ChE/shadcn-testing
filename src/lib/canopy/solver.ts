// ─── Canopy Unit Operation Solver ──────────────────────────────────────────────
// Solves individual unit operations given their inlet streams and parameters.
// ──────────────────────────────────────────────────────────────────────────────

import type { CanopyCompound, MaterialStream, UnitOperation } from './types';
import { flashPT_Ideal, streamEnthalpy, enthalpyIG, enthalpyLiquid } from './thermo';

// ────────────────────────────────────────────────────────────────
// Mixer: combines N inlet streams into 1 outlet stream
// ────────────────────────────────────────────────────────────────

export function solveMixer(
  compounds: CanopyCompound[],
  inlets: MaterialStream[],
  outletPressure_Pa: number,
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
  let H_in_total = 0;
  for (const inlet of inlets) {
    H_in_total += inlet.totalFlow_molps * inlet.H_Jpmol;
  }

  // Flash at outlet P to find T (energy balance)
  // For simplicity in Phase 1: weighted average temperature as starting guess
  let T_avg = 0;
  if (totalMol > 0) {
    for (const inlet of inlets) {
      T_avg += inlet.T_K * inlet.totalFlow_molps / totalMol;
    }
  } else {
    T_avg = 298.15;
  }

  // Flash at T_avg, P_out
  const flash = flashPT_Ideal(compounds, z, T_avg, outletPressure_Pa);

  return {
    outlet: {
      T_K: T_avg,
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
  params: { targetT_K: number; outletP_Pa: number },
): { outlet: Partial<MaterialStream>; duty_W: number } {
  const z = [...inlet.moleFractions];
  const flash = flashPT_Ideal(compounds, z, params.targetT_K, params.outletP_Pa);

  const H_in = inlet.H_Jpmol;
  const H_out = flash.H_Jpmol;
  const duty_W = inlet.totalFlow_molps * (H_out - H_in); // W

  return {
    outlet: {
      T_K: params.targetT_K,
      P_Pa: params.outletP_Pa,
      totalFlow_molps: inlet.totalFlow_molps,
      moleFractions: z,
      phase: flash.phase,
      vaporFraction: flash.vaporFraction,
      H_Jpmol: H_out,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    duty_W,
  };
}

// ────────────────────────────────────────────────────────────────
// Flash Drum: separates a two-phase stream into V and L
// ────────────────────────────────────────────────────────────────

export interface FlashDrumResult {
  vapor: Partial<MaterialStream>;
  liquid: Partial<MaterialStream>;
  duty_W: number;
}

export function solveFlashDrum(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  params: { T_K: number; P_Pa: number },
): FlashDrumResult {
  const z = [...inlet.moleFractions];
  const flash = flashPT_Ideal(compounds, z, params.T_K, params.P_Pa);

  const F = inlet.totalFlow_molps;
  const V_flow = F * flash.vaporFraction;
  const L_flow = F * (1 - flash.vaporFraction);

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
  };
}

// ────────────────────────────────────────────────────────────────
// Valve: isenthalpic pressure drop
// ────────────────────────────────────────────────────────────────

export function solveValve(
  compounds: CanopyCompound[],
  inlet: MaterialStream,
  params: { outletP_Pa: number },
): { outlet: Partial<MaterialStream>; duty_W: number } {
  // Isenthalpic: find T at new P that gives same enthalpy
  const z = [...inlet.moleFractions];
  const H_target = inlet.H_Jpmol;

  // Newton-Raphson on T to match enthalpy
  let T = inlet.T_K;
  for (let iter = 0; iter < 50; iter++) {
    const flash = flashPT_Ideal(compounds, z, T, params.outletP_Pa);
    const H = flash.H_Jpmol;
    const err = H - H_target;
    if (Math.abs(err) < 0.01) {
      return {
        outlet: {
          T_K: T,
          P_Pa: params.outletP_Pa,
          totalFlow_molps: inlet.totalFlow_molps,
          moleFractions: z,
          phase: flash.phase,
          vaporFraction: flash.vaporFraction,
          H_Jpmol: H,
          x_liquid: flash.x,
          y_vapor: flash.y,
          solved: true,
        },
        duty_W: 0,
      };
    }
    // Numerical derivative
    const dT = 0.1;
    const flash2 = flashPT_Ideal(compounds, z, T + dT, params.outletP_Pa);
    const dHdT = (flash2.H_Jpmol - H) / dT;
    if (Math.abs(dHdT) < 1e-15) break;
    T -= err / dHdT;
    if (T < 100) T = 100;
    if (T > 1000) T = 1000;
  }

  // Fallback: just do flash at inlet T
  const flash = flashPT_Ideal(compounds, z, T, params.outletP_Pa);
  return {
    outlet: {
      T_K: T,
      P_Pa: params.outletP_Pa,
      totalFlow_molps: inlet.totalFlow_molps,
      moleFractions: z,
      phase: flash.phase,
      vaporFraction: flash.vaporFraction,
      H_Jpmol: flash.H_Jpmol,
      x_liquid: flash.x,
      y_vapor: flash.y,
      solved: true,
    },
    duty_W: 0,
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
