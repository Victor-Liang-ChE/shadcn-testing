// ─── Canopy Process Simulator: Core Types ──────────────────────────────────────

export type FluidPackage =
  | 'NRTL'
  | 'Electrolyte-NRTL'
  | 'Wilson'
  | 'UNIQUAC'
  | 'UNIFAC'
  | 'UNIFAC-DMD'
  | 'Peng-Robinson'
  | 'SRK'
  | 'Ideal';

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
  // eqNo selects the DIPPR equation form (100 = polynomial, 106 = Watson, 114 = Zabransky, etc.)
  dipprCoeffs: Record<string, { A: number; B: number; C: number; D: number; E?: number; eqNo?: number }>;
  // DIPPR 107 (Aly-Lee) ideal gas Cp coefficients: Cp = A + B·[(C/T)/sinh(C/T)]² + D·[(E/T)/cosh(E/T)]²
  cpigdp?: { A: number; B: number; C: number; D: number; E: number; Tmin_K?: number; Tmax_K?: number };
  // Extended Antoine (PLXANT) parameters
  antoine: {
    C1: number; C2: number; C3: number; C4: number;
    C5: number; C6: number; C7: number;
    Tmin_K: number; Tmax_K: number;
  } | null;
  // UNIQUAC structural parameters (from GMUQR / GMUQQ)
  uniquac_r?: number;         // Volume parameter (van der Waals volume)
  uniquac_q?: number;         // Surface area parameter (van der Waals area)
  // UNIFAC subgroup decomposition: map of subgroupID → count
  unifacGroups?: Record<number, number>;
  // Electrolyte metadata for electrolyte property packages
  chargeNumber?: number;
  isElectrolyteSolvent?: boolean;
  // Standard heat of formation at 298.15 K (J/mol, ideal gas)
  Hf_298_Jmol?: number;
  // Standard Gibbs energy of formation at 298.15 K (J/mol, ideal gas)
  Gf_298_Jmol?: number;
  // Mathias-Copeman alpha function parameters for PR/SRK EOS (better Psat fit)
  // α = [1 + C1(1-√Tr) + C2(1-√Tr)² + C3(1-√Tr)³]²
  mathiasCopeman?: { C1: number; C2: number; C3: number };
  twuAlpha?: { L: number; M: number; N: number };
  // Wagner vapor pressure: ln(Pr) = (a·τ + b·τ^1.5 + c·τ^3 + d·τ^6) / Tr  where τ = 1 - Tr
  wagnerVP?: { a: number; b: number; c: number; d: number; Tmin_K?: number; Tmax_K?: number };
  // Watson Hvap reference point for extrapolation: ΔHvap(T) = ΔHvap(T1)·((1-Tr)/(1-Tr1))^n
  watsonHvap?: { Hvap_ref_Jmol: number; T_ref_K: number; n?: number };
  // Aspen polynomial Cp_ig: A + B·T + C·T² + D·T³ + E·T⁴ + F·T⁵ [J/(kmol·K)]
  cpigPoly?: { A: number; B: number; C: number; D: number; E: number; F?: number; Tmin_K?: number; Tmax_K?: number };
  // DIPPR 116 liquid density: ρ = A + B·τ^0.35 + C·τ^(2/3) + D·τ + E·τ^(4/3) [kmol/m³]
  liquidDensityDIPPR116?: { A: number; B: number; C: number; D: number; E: number };
  // Critical volume (m³/kmol) – used for mixing rules and liquid density
  Vc_m3pkmol?: number;
  // Critical compressibility factor Zc = Pc·Vc/(R·Tc)
  Zc?: number;
  // Rackett equation parameter for liquid density (ZRA)
  rackett_ZRA?: number;
  // Molar volume at normal boiling point (m³/mol) – used in Le Bas and transport correlations
  Vb_m3pmol?: number;
  // Normal boiling point (K) – used by many correlations
  Tb_K?: number;
  // Freezing/melting point (K) – used for SLE
  Tf_K?: number;
  // Heat of fusion at melting point (J/mol) – used for SLE
  Hfus_Jmol?: number;
  // Heat of combustion at 298.15 K (J/mol, negative = exothermic)
  Hcomb_Jmol?: number;
  // Heat of vaporization at normal boiling point (J/mol)
  Hvap_nb_Jmol?: number;
  // Standard liquid molar volume (m³/kmol) at 60°F / 15.56°C
  VLSTD_m3pkmol?: number;
  // Mathias-Copeman alpha function for SRK EOS
  mathiasCopeman_SRK?: { C1: number; C2: number; C3: number };
  // PR EOS Peneloux volume translation (m³/kmol)
  volumeTranslation_PR?: number;
  // Hildebrand solubility parameter at 25°C ((J/m³)^0.5)
  delta_Jm3_05?: number;
  // Radius of gyration (m)
  radiusOfGyration_m?: number;
  // PC-SAFT pure-component parameters (Gross & Sadowski 2001)
  pcSaft?: {
    m: number;       // Segment number
    sigma: number;   // Segment diameter (Å)
    epsilon_k: number; // Segment energy parameter ε/k (K)
  };
  // CPA (Cubic-Plus-Association) EOS parameters (Kontogeorgis et al.)
  cpa?: {
    Tc_K: number;              // CPA-fitted critical temperature
    Pc_Pa: number;             // CPA-fitted critical pressure
    m: number;                 // alpha-function parameter
    epsilon_Pa_m3: number;     // association energy (0 for non-associating)
    beta: number;              // effective association volume (0 for non-associating)
  };
  // Brelvi-O'Connell characteristic volume (m³/kmol) for Henry's Poynting correction
  brelviOConnellV?: number;
  // Transport property DIPPR coefficients
  transportCoeffs?: {
    /** DIPPR 101: ln(μ) = A + B/T + C*ln(T) + D*T^E  [Pa·s] */
    LiquidViscosity?: { A: number; B: number; C: number; D: number; E?: number };
    /** DIPPR 102: μ = A*T^B / (1 + C/T + D/T²)  [Pa·s] */
    VaporViscosity?: { A: number; B: number; C: number; D: number; E?: number };
    /** DIPPR 100/123: λ evaluated as polynomial (100) or Sato-Riedel (123); check TRNSWT[2] */
    LiquidThermalConductivity?: { A: number; B: number; C: number; D: number; E?: number };
    /** DIPPR 102: λ = A*T^B / (1 + C/T + D/T²)  [W/(m·K)] */
    VaporThermalConductivity?: { A: number; B: number; C: number; D: number; E?: number };
    /** DIPPR 106: σ = A*(1-Tr)^(B + C*Tr + D*Tr² + E*Tr³)  [N/m] */
    SurfaceTension?: { A: number; B: number; C: number; D: number; E?: number };
  };
  // Specific gravity at 60°F/60°F (water = 1)
  specificGravity?: number;
  // Flash point (K) — closed cup
  flashPoint_K?: number;
  // Autoignition temperature (K)
  autoignitionTemp_K?: number;
  // Thermodynamic equation switches [S1..S8]: S1=ρ_solid, S2=ρ_liquid, S3=Psat,
  // S4=ΔHvap, S5=Cp_solid, S6=Cp_liquid, S7=Cp_ig, S8=B_virial
  THRSWT?: [number, number, number, number, number, number, number, number];
  // Transport equation switches [S1..S5]: S1=μ_liq, S2=μ_vap, S3=λ_liq, S4=λ_vap, S5=σ
  TRNSWT?: [number, number, number, number, number];
  // Solid-phase extended Antoine (sublimation pressure)
  // ln(P/Pa) = C1 + C2/(T+C3) + C4·T + C5·ln(T) + C6·T^C7
  psxant?: {
    C1: number; C2: number; C3: number; C4: number;
    C5: number; C6: number; C7: number;
    Tmin_K: number; Tmax_K: number;
  };
  // Water solubility in hydrocarbon phase: ln(x_w) = A + B/T + C·ln(T)
  waterSolubility?: { A: number; B: number; C: number; Tmin_K: number; Tmax_K: number };
  // Hydrocarbon solubility in water: ln(x_hc) = A + B/T + C·ln(T)
  hcSolubility?: { A: number; B: number; C: number; Tmin_K: number; Tmax_K: number };
  // Schwartentruber-Renon polar parameters for PR EOS
  schwartentruberRenon_PR?: { p1: number; p2: number; p3: number };
  // Schwartentruber-Renon polar parameters for SRK EOS
  schwartentruberRenon_SRK?: { p1: number; p2: number; p3: number };
  // COSTALD liquid volume constant (m³/kmol)
  VLCVT1?: number;
  // COSTALD characteristic volume V* (m³/kmol) — Hankinson-Thomson
  VSTCTD?: number;
  // Corrected acentric factor for COSTALD
  omegaCostald?: number;

  // ── Thermochemical formation/combustion data (J/kmol unless noted) ──
  /** Standard enthalpy of formation at 298.15 K, ideal gas (J/kmol) */
  dhForm_Jpkmol?: number;
  /** Standard Gibbs energy of formation at 298.15 K, ideal gas (J/kmol) */
  dgForm_Jpkmol?: number;
  /** Standard ideal-gas entropy at 298.15 K (J/kmol·K) */
  entropy_JpkmolK?: number;
  /** Standard heat of combustion, net (J/kmol, negative=exothermic) */
  hCombustion_Jpkmol?: number;
  /** Heat of fusion at triple point (J/kmol) */
  hFusion_Jpkmol?: number;
  /** Heat of vaporization at normal boiling point (J/kmol) */
  dhVapAtNBP_Jpkmol?: number;
  /** Liquid Cp departure from ideal-gas Cp at 298.15 K (J/kmol·K) */
  cpDepartureLiq_JpkmolK?: number;

  // ── Physical constants ──
  /** Standard liquid molar volume at 60°F (m³/kmol) */
  vLiqStd_m3pkmol?: number;
  /** Liquid molar volume at normal boiling point (m³/kmol) */
  vLiqAtNBP_m3pkmol?: number;
  /** Dipole moment (C·m) */
  dipoleMoment_Cm?: number;
  /** Refractive index at 20°C (sodium D-line) */
  refractiveIndex?: number;
  /** Lower flammability limit in air (vol%) */
  flammabilityLower_pct?: number;
  /** Upper flammability limit in air (vol%) */
  flammabilityUpper_pct?: number;
  /** Triple point temperature (K) */
  triplePointT_K?: number;
  /** Triple point pressure (Pa) */
  triplePointP_Pa?: number;
  /** Dielectric constant ε(T) = e1 + e2·(1/T - 1/e3), or simply e1 at reference temp */
  dielectricConst?: { e1: number; e2: number; e3_K: number };
  /** COSMO-SAC cavity volume (Å³, dimensionless in pplcdefs) */
  cosmoVolume?: number;
  /** Solid-liquid Cp difference at triple point (J/kmol·K) */
  solidLiquidCpDiff_JpkmolK?: number;
  /** Aniline point (K) — petroleum characterization */
  anilinePoint_K?: number;
  /** Cetane number — petroleum ignition quality */
  cetaneNumber?: number;
  /** Hayden-O'Connell solvation parameter η for second virial cross-coefficients */
  hocEta?: Record<string, number>;
}

/**
 * NRTL binary interaction parameters for a pair (i,j).
 * τ_ij = aij + bij/T + eij·ln(T) + fij·T
 * G_ij = exp(-αij·τ_ij)  where αij = cij + dij·(T - 273.15)
 */
export interface NrtlBinaryParams {
  aij: number; aji: number;
  bij: number; bji: number;  // K
  cij: number;                // alpha (non-randomness), symmetric
  dij?: number;               // temperature-dependent alpha
  eij?: number; eji?: number;
  fij?: number; fji?: number;
  Tlower_K?: number; Tupper_K?: number;
  source?: string;            // e.g. "APV140 VLE-IG", "NIST-IG"
}

/**
 * Peng-Robinson / SRK kij: kij(T) = kij1 + kij2·T + kij3·T²
 */
export interface EosKij {
  kij1: number;
  kij2?: number;
  kij3?: number;
  source?: string;
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
export type UnitOpType =
  | 'Mixer' | 'Heater' | 'Flash' | 'Splitter' | 'Valve'
  | 'Pump' | 'Compressor' | 'HeatX' | 'Pipe' | 'RStoic' | 'Column' | 'RadFrac' | 'RateFrac'
  | 'ThreePhaseFlash' | 'CSTR' | 'PFR' | 'Absorber' | 'Extractor' | 'ComponentSeparator'
  | 'RYield' | 'Decanter' | 'REquil' | 'RGibbs'
  | 'MHeatX' | 'RBatch' | 'Crystallizer' | 'Crusher' | 'Dryer'
  | 'Membrane' | 'Cyclone' | 'Filter' | 'Screen' | 'Centrifuge'
  | 'SteamDrum' | 'SteamHeater' | 'SteamTurbine' | 'SteamValve' | 'SteamHeader' | 'SteamTrap'
  | 'PIDController' | 'LeadLag' | 'DeadTime' | 'SignalSelector'
  | 'Analyzer' | 'ChargeBalance' | 'ElectrolyteEquilibrium';

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
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  params: Record<string, any>;
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

// ─── Reaction / Stoichiometry Types ───────────────────────────────────────────

/** A single stoichiometric reaction */
export interface StoichiometricReaction {
  /** Stoichiometric coefficients per component (negative = reactant, positive = product) */
  coefficients: number[];
  /** Index of the key component for conversion specification */
  keyComponentIndex: number;
  /** Fractional conversion of the key component (0..1) */
  conversion: number;
}

/** Column specification for shortcut distillation */
export interface ColumnSpec {
  /** Light key component index */
  lightKeyIndex: number;
  /** Heavy key component index */
  heavyKeyIndex: number;
  /** Recovery of light key in distillate (0..1) */
  lightKeyRecovery: number;
  /** Recovery of heavy key in bottoms (0..1) */
  heavyKeyRecovery: number;
  /** Reflux ratio multiplier over minimum (e.g. 1.2 = 1.2× Rmin) */
  refluxRatioMultiplier: number;
  /** Condenser pressure (Pa) */
  condenserP_Pa: number;
  /** Reboiler pressure (Pa) */
  reboilerP_Pa: number;
  /** Condenser type */
  condenserType: 'Total' | 'Partial';
}

/** Default empty stream */
export function createDefaultStream(
  id: string = 'stream',
  name: string = 'Stream',
  nCompounds: number = 2,
): MaterialStream {
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

// ────────────────────────────────────────────────────────────────
// Design Specification
// ────────────────────────────────────────────────────────────────

/** A design spec defines a target-variable ↔ manipulated-variable pair.
 *  The solver adjusts the manipulated variable until the target is met.
 *
 *  Example: set stream S3 temperature to 350 K by adjusting Heater-1 targetT_K.
 */
export interface DesignSpec {
  id: string;
  name: string;

  /** What we want to achieve */
  target: {
    /** 'stream' or 'unit' */
    objectType: 'stream' | 'unit';
    /** ID of the stream or unit */
    objectId: string;
    /** Property path — e.g. 'T_K', 'P_Pa', 'moleFractions[0]', 'vaporFraction' */
    property: string;
    /** Desired value */
    value: number;
  };

  /** What the solver manipulates */
  manipulated: {
    /** 'unit' parameter */
    unitId: string;
    /** Parameter name in unit.params — e.g. 'targetT_K', 'outletP_Pa' */
    paramName: string;
    /** Lower bound for the manipulated variable */
    lowerBound: number;
    /** Upper bound for the manipulated variable */
    upperBound: number;
  };

  /** Tolerance (absolute) for convergence */
  tolerance: number;

  /** Whether this spec is active */
  active: boolean;
}

/**
 * Sensitivity analysis: vary one (or two) parameters and record outputs.
 * Equivalent to Aspen's Sensitivity Analysis model.
 */
export interface SensitivityAnalysis {
  id: string;
  name: string;

  /** Variable to vary */
  variable: {
    unitId: string;
    paramName: string;
    values: number[]; // explicit list of values, or generate from range
  };

  /** Optional second variable (2D sweep) */
  variable2?: {
    unitId: string;
    paramName: string;
    values: number[];
  };

  /** Outputs to record at each sweep point */
  outputs: Array<{
    objectType: 'stream' | 'unit';
    objectId: string;
    property: string;
    label: string;
  }>;

  /** Results matrix: rows = sweep points, columns = output values */
  results: number[][];

  active: boolean;
}

/**
 * Calculator block: evaluate an expression before or after solving.
 * Similar to Aspen Plus Calculator blocks.
 */
export interface CalculatorBlock {
  id: string;
  name: string;

  /** When to execute: 'before' or 'after' solving */
  executionOrder: 'before' | 'after';

  /** Import variables: read from streams/units */
  imports: Array<{
    varName: string;
    objectType: 'stream' | 'unit';
    objectId: string;
    property: string;
  }>;

  /** Export variables: write to unit params */
  exports: Array<{
    varName: string;
    unitId: string;
    paramName: string;
  }>;

  /** Simple arithmetic expression (evaluated with imported vars as context) */
  expression: string;

  active: boolean;
}
