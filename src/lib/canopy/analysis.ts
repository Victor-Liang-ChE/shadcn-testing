// ─── Canopy Process Simulator: Property Analysis Module ─────────────────────
//
// Generates T-xy, P-xy, and residue curve map data for binary/ternary systems.
// Uses the same thermo engine (flashPT, bubble/dew point iteration).
// ────────────────────────────────────────────────────────────────────────────────

import type { CanopyCompound } from './types';
import { flashPT } from './thermo';
import type { InteractionParams } from './thermo';

// ─── T-xy Diagram ───────────────────────────────────────────────

export interface TxyPoint {
  x1: number;          // Liquid mole fraction of component 1
  y1: number;          // Vapor mole fraction of component 1
  T_K: number;         // Bubble/dew temperature at this composition
}

export interface TxyData {
  comp1: string;
  comp2: string;
  P_Pa: number;
  bubbleCurve: TxyPoint[];     // x1 vs T (bubble line)
  dewCurve: TxyPoint[];        // y1 vs T (dew line)
}

/**
 * Generate T-xy diagram data for a binary system at constant pressure.
 * Iterates over x1 from 0 to 1 and finds bubble temperature at each point.
 */
export function generateTxy(
  compounds: CanopyCompound[],
  comp1Index: number,
  comp2Index: number,
  P_Pa: number,
  fluidPackage: string,
  interactionParams?: InteractionParams,
  nPoints: number = 51,
): TxyData {
  const c1 = compounds[comp1Index];
  const c2 = compounds[comp2Index];
  const twoCompounds = [c1, c2];

  const bubbleCurve: TxyPoint[] = [];
  const dewCurve: TxyPoint[] = [];

  for (let i = 0; i <= nPoints; i++) {
    const x1 = i / nPoints;
    const z = [x1, 1 - x1];

    // Find bubble temperature by bisection
    let Tlow = 200, Thigh = 600;
    // Bracket: at Tlow should be subcooled (all liquid), at Thigh should be superheated (all vapor)
    for (let iter = 0; iter < 60; iter++) {
      const Tmid = (Tlow + Thigh) / 2;
      const flash = flashPT(twoCompounds, z, Tmid, P_Pa,
        fluidPackage as any, interactionParams);
      if (flash.vaporFraction > 1e-8) {
        Thigh = Tmid;
      } else {
        Tlow = Tmid;
      }
      if (Thigh - Tlow < 0.01) break;
    }
    const Tbubble = (Tlow + Thigh) / 2;

    // Get equilibrium y at bubble point
    const flashBub = flashPT(twoCompounds, z, Tbubble + 0.01, P_Pa,
      fluidPackage as any, interactionParams);
    const y1 = flashBub.y[0];

    bubbleCurve.push({ x1, y1, T_K: Tbubble });
    dewCurve.push({ x1: y1, y1: x1, T_K: Tbubble });
  }

  // Sort dew curve by x1 for clean plotting
  dewCurve.sort((a, b) => a.x1 - b.x1);

  return {
    comp1: c1.displayName,
    comp2: c2.displayName,
    P_Pa,
    bubbleCurve,
    dewCurve,
  };
}

// ─── P-xy Diagram ───────────────────────────────────────────────

export interface PxyPoint {
  x1: number;
  y1: number;
  P_Pa: number;
}

export interface PxyData {
  comp1: string;
  comp2: string;
  T_K: number;
  bubbleCurve: PxyPoint[];
  dewCurve: PxyPoint[];
}

/**
 * Generate P-xy diagram data for a binary system at constant temperature.
 */
export function generatePxy(
  compounds: CanopyCompound[],
  comp1Index: number,
  comp2Index: number,
  T_K: number,
  fluidPackage: string,
  interactionParams?: InteractionParams,
  nPoints: number = 51,
): PxyData {
  const c1 = compounds[comp1Index];
  const c2 = compounds[comp2Index];
  const twoCompounds = [c1, c2];

  const bubbleCurve: PxyPoint[] = [];
  const dewCurve: PxyPoint[] = [];

  for (let i = 0; i <= nPoints; i++) {
    const x1 = i / nPoints;
    const z = [x1, 1 - x1];

    // Find bubble pressure by bisection
    let Plow = 100, Phigh = 1e7; // 100 Pa to 100 bar
    for (let iter = 0; iter < 60; iter++) {
      const Pmid = (Plow + Phigh) / 2;
      const flash = flashPT(twoCompounds, z, T_K, Pmid,
        fluidPackage as any, interactionParams);
      if (flash.vaporFraction > 1e-8) {
        Plow = Pmid;
      } else {
        Phigh = Pmid;
      }
      if (Phigh - Plow < 10) break;
    }
    const Pbubble = (Plow + Phigh) / 2;

    const flashBub = flashPT(twoCompounds, z, T_K, Pbubble - 10,
      fluidPackage as any, interactionParams);
    const y1 = flashBub.y[0];

    bubbleCurve.push({ x1, y1, P_Pa: Pbubble });
    dewCurve.push({ x1: y1, y1: x1, P_Pa: Pbubble });
  }

  dewCurve.sort((a, b) => a.x1 - b.x1);

  return {
    comp1: c1.displayName,
    comp2: c2.displayName,
    T_K,
    bubbleCurve,
    dewCurve,
  };
}

// ─── x-y Diagram ────────────────────────────────────────────────

export interface XYPoint {
  x1: number;    // Liquid mole fraction
  y1: number;    // Vapor mole fraction
}

export interface XYData {
  comp1: string;
  comp2: string;
  P_Pa: number;
  points: XYPoint[];
}

/**
 * Generate x-y equilibrium diagram for a binary system at constant P.
 */
export function generateXY(
  compounds: CanopyCompound[],
  comp1Index: number,
  comp2Index: number,
  P_Pa: number,
  fluidPackage: string,
  interactionParams?: InteractionParams,
  nPoints: number = 51,
): XYData {
  const c1 = compounds[comp1Index];
  const c2 = compounds[comp2Index];
  const twoCompounds = [c1, c2];

  const points: XYPoint[] = [];

  for (let i = 0; i <= nPoints; i++) {
    const x1 = i / nPoints;
    const z = [x1, 1 - x1];

    // Find bubble T, then get y
    let Tlow = 200, Thigh = 600;
    for (let iter = 0; iter < 60; iter++) {
      const Tmid = (Tlow + Thigh) / 2;
      const flash = flashPT(twoCompounds, z, Tmid, P_Pa,
        fluidPackage as any, interactionParams);
      if (flash.vaporFraction > 1e-8) {
        Thigh = Tmid;
      } else {
        Tlow = Tmid;
      }
      if (Thigh - Tlow < 0.01) break;
    }
    const Tbub = (Tlow + Thigh) / 2;
    const flash = flashPT(twoCompounds, z, Tbub + 0.01, P_Pa,
      fluidPackage as any, interactionParams);

    points.push({ x1, y1: flash.y[0] });
  }

  return {
    comp1: c1.displayName,
    comp2: c2.displayName,
    P_Pa,
    points,
  };
}

// ─── Relative Volatility ────────────────────────────────────────

export interface AlphaPoint {
  x1: number;
  alpha: number;
  T_K: number;
}

/**
 * Calculate relative volatility α₁₂ = (y₁/x₁) / (y₂/x₂) across composition range.
 */
export function generateAlpha(
  compounds: CanopyCompound[],
  comp1Index: number,
  comp2Index: number,
  P_Pa: number,
  fluidPackage: string,
  interactionParams?: InteractionParams,
  nPoints: number = 51,
): AlphaPoint[] {
  const c1 = compounds[comp1Index];
  const c2 = compounds[comp2Index];
  const twoCompounds = [c1, c2];
  const result: AlphaPoint[] = [];

  for (let i = 1; i < nPoints; i++) {
    const x1 = i / nPoints;
    const z = [x1, 1 - x1];

    let Tlow = 200, Thigh = 600;
    for (let iter = 0; iter < 60; iter++) {
      const Tmid = (Tlow + Thigh) / 2;
      const flash = flashPT(twoCompounds, z, Tmid, P_Pa,
        fluidPackage as any, interactionParams);
      if (flash.vaporFraction > 1e-8) Thigh = Tmid; else Tlow = Tmid;
      if (Thigh - Tlow < 0.01) break;
    }
    const T = (Tlow + Thigh) / 2;
    const flash = flashPT(twoCompounds, z, T + 0.01, P_Pa,
      fluidPackage as any, interactionParams);
    const y1 = flash.y[0];
    const x2 = 1 - x1;
    const y2 = 1 - y1;
    const alpha = x2 > 1e-10 && y2 > 1e-10
      ? (y1 / x1) / (y2 / x2)
      : NaN;

    result.push({ x1, alpha, T_K: T });
  }

  return result;
}
