// ─── Canopy Process Simulator: Core Types ──────────────────────────────────────

export type FluidPackage = 'NRTL' | 'Wilson' | 'UNIQUAC' | 'UNIFAC' | 'Peng-Robinson' | 'SRK' | 'Ideal';

export type Phase = 'Liquid' | 'Vapor' | 'Two-Phase';

/** Compound as stored in the simulation — resolved from the DB. */
export interface CanopyCompound {
  name: string;               // Canonical Aspen name (e.g. "BENZENE")
  displayName: string;         // Formatted display name (e.g. "Benzene")
  molecularWeight: number;     // g/mol
  Tc_K: number;              // Critical temperature (K)
  Pc_Pa: number;             // Critical pressure (Pa)
  omega: number;             // Acentric factor
  // DIPPR property coefficients (keyed by property name)
  dipprCoeffs: Record<string, { A: number; B: number; C: number; D: number; E?: number }>;
  // Extended Antoine (PLXANT) parameters
  antoine: {
    C1: number; C2: number; C3: number; C4: number;
    C5: number; C6: number; C7: number;
    Tmin_K: number; Tmax_K: number;
  } | null;
}

/** A material stream connecting two unit operations. */
export interface MaterialStream {
  id: string;
  name: string;
  // Thermodynamic state
  T_K: number;                 // Temperature (K)
  P_Pa: number;               // Pressure (Pa)
  totalFlow_molps: number;     // Total molar flow (mol/s)
  moleFractions: number[];     // Mole fractions (same order as simulation compounds)
  // Calculated results (filled after solving)
  phase: Phase;
  vaporFraction: number;       // 0 = all liquid, 1 = all vapor
  H_Jpmol: number;            // Stream molar enthalpy (J/mol)
  // Per-phase compositions (after flash)
  x_liquid: number[];          // Liquid mole fractions
  y_vapor: number[];           // Vapor mole fractions
  // Connection info
  sourceUnitId: string | null;
  sourcePortId: string | null;
  targetUnitId: string | null;
  targetPortId: string | null;
  // Status
  solved: boolean;
}

/** Types of unit operations supported */
export type UnitOpType = 'Mixer' | 'Heater' | 'Flash' | 'Splitter' | 'Valve';

/** Port definition on a unit operation */
export interface UnitPort {
  id: string;
  label: string;
  type: 'inlet' | 'outlet';
  position: 'left' | 'right' | 'top' | 'bottom';
  streamId: string | null;
}

/** Base unit operation interface */
export interface UnitOperation {
  id: string;
  name: string;
  type: UnitOpType;
  ports: UnitPort[];
  // Unit-specific parameters
  params: Record<string, number | string | boolean>;
  // Results
  solved: boolean;
  duty_W: number;             // Heat duty (W) — positive = heat added
  errors: string[];
}

/** Simulation case — the top-level container */
export interface SimulationCase {
  id: string;
  name: string;
  compounds: CanopyCompound[];
  fluidPackage: FluidPackage;
  streams: Record<string, MaterialStream>;
  units: Record<string, UnitOperation>;
  solved: boolean;
}

/** Default empty stream */
export function createDefaultStream(id: string, name: string, nCompounds: number): MaterialStream {
  return {
    id,
    name,
    T_K: 298.15,
    P_Pa: 101325,
    totalFlow_molps: 0,
    moleFractions: new Array(nCompounds).fill(1 / Math.max(nCompounds, 1)),
    phase: 'Liquid',
    vaporFraction: 0,
    H_Jpmol: 0,
    x_liquid: new Array(nCompounds).fill(1 / Math.max(nCompounds, 1)),
    y_vapor: new Array(nCompounds).fill(1 / Math.max(nCompounds, 1)),
    sourceUnitId: null,
    sourcePortId: null,
    targetUnitId: null,
    targetPortId: null,
    solved: false,
  };
}
