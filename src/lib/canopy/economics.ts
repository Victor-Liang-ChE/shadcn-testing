// ─── Canopy Process Simulator: Economics / Costing Module ──────────────────────
//
// Equipment cost correlations from:
//   Turton, R., Shaeiwitz, J.A., Bhattacharyya, D., Whiting, W.B.
//   "Analysis, Synthesis and Design of Chemical Processes", 5th ed., Pearson, 2018.
//
//   Purchase cost: log₁₀(Cp°) = K₁ + K₂·log₁₀(A) + K₃·[log₁₀(A)]²  [Table A.1]
//   Bare module:   CBM = Cp° · FBM                                     [Table A.4]
//   Utility rates: from Table 8.3 (2013 USD, CEPCI 567)
//   CEPCI scaling: Turton Eq. (A.16): Cost₂ = Cost₁ · (CEPCI₂/CEPCI₁)
//
// These are the same correlations used by Aspen Process Economic Analyzer (APEA)
// as default/screening estimates. All costs in USD (2024 basis, CEPCI ≈ 800).
// ────────────────────────────────────────────────────────────────────────────────

import type { UnitOperation, MaterialStream, CanopyCompound } from './types';

// ─── Constants ──────────────────────────────────────────────────

/** Chemical Engineering Plant Cost Index — base year and current */
const CEPCI_BASE = 567;   // 2013 (Turton 5th ed. basis)
const CEPCI_CURRENT = 800; // 2024 estimate
const CEPCI_RATIO = CEPCI_CURRENT / CEPCI_BASE;

/** Operating hours per year */
const HOURS_PER_YEAR = 8400;

// ─── Utility cost rates (Turton Table 8.3, 2013 USD) ────────────

export interface UtilityRates {
  /** Electricity cost ($/kWh) */
  electricity: number;
  /** Cooling water cost ($/GJ) */
  coolingWater: number;
  /** Low-pressure steam cost ($/GJ) — 150°C */
  lpSteam: number;
  /** Medium-pressure steam cost ($/GJ) — 250°C */
  mpSteam: number;
  /** High-pressure steam cost ($/GJ) — 400°C */
  hpSteam: number;
  /** Refrigeration cost ($/GJ) */
  refrigeration: number;
  /** Natural gas / fuel cost ($/GJ) */
  fuel: number;
}

export const DEFAULT_UTILITY_RATES: UtilityRates = {
  electricity: 0.07,     // $/kWh
  coolingWater: 0.354,   // $/GJ
  lpSteam: 13.28,        // $/GJ
  mpSteam: 14.19,        // $/GJ
  hpSteam: 17.70,        // $/GJ
  refrigeration: 4.43,   // $/GJ
  fuel: 3.68,            // $/GJ
};

// ─── Equipment purchase cost correlations (Turton Table A.1) ────
//
// log₁₀(Cp°) = K₁ + K₂·log₁₀(A) + K₃·(log₁₀(A))²
// where A = sizing parameter, Cp° = purchase cost at base CEPCI
//
// Bare module cost: CBM = Cp° × FBM  (Turton Table A.4)
// where FBM = B₁ + B₂·Fp·FM (pressure & material factors)

interface CostCorrelation {
  K1: number; K2: number; K3: number;
  /** Bare-module factor (simple) */
  fbm: number;
  /** Sizing parameter description */
  sizingUnit: string;
  /** Min/Max valid range for A */
  Amin: number; Amax: number;
}

const COST_DATA: Record<string, CostCorrelation> = {
  // Heat exchangers (shell & tube) — A in m²
  heatExchanger: {
    K1: 4.3247, K2: -0.3030, K3: 0.1634,
    fbm: 3.17, sizingUnit: 'm²', Amin: 10, Amax: 1000,
  },
  // Centrifugal pump — A in kW
  pump: {
    K1: 3.3892, K2: 0.0536, K3: 0.1538,
    fbm: 3.30, sizingUnit: 'kW', Amin: 1, Amax: 300,
  },
  // Centrifugal compressor — A in kW
  compressor: {
    K1: 2.2897, K2: 1.3604, K3: -0.1027,
    fbm: 2.15, sizingUnit: 'kW', Amin: 450, Amax: 3000,
  },
  // Fired heater (furnace) — A in kW
  firedHeater: {
    K1: 7.3488, K2: -1.1666, K3: 0.2028,
    fbm: 2.19, sizingUnit: 'kW', Amin: 1000, Amax: 100000,
  },
  // Pressure vessel (vertical) — A in m³
  vessel: {
    K1: 3.4974, K2: 0.4485, K3: 0.1074,
    fbm: 4.16, sizingUnit: 'm³', Amin: 0.3, Amax: 520,
  },
  // Tray column (sieve trays) — A per tray in m²
  tray: {
    K1: 2.9949, K2: 0.4465, K3: 0.3961,
    fbm: 1.0, sizingUnit: 'm²', Amin: 0.07, Amax: 12.3,
  },
  // Tower / column shell (vertical vessel) — A in m³
  column: {
    K1: 3.4974, K2: 0.4485, K3: 0.1074,
    fbm: 4.16, sizingUnit: 'm³', Amin: 0.3, Amax: 520,
  },
  // Reactor (jacketed, agitated) — A in m³
  reactor: {
    K1: 4.1052, K2: 0.5320, K3: -0.0005,
    fbm: 4.16, sizingUnit: 'm³', Amin: 0.5, Amax: 100,
  },
  // Valve — minimal cost
  valve: {
    K1: 2.5, K2: 0.0, K3: 0.0,
    fbm: 1.0, sizingUnit: '-', Amin: 0, Amax: 1e9,
  },
  // Mixer (static) — minimal cost
  mixer: {
    K1: 2.0, K2: 0.0, K3: 0.0,
    fbm: 1.0, sizingUnit: '-', Amin: 0, Amax: 1e9,
  },
  // Generic / fallback
  generic: {
    K1: 3.5, K2: 0.0, K3: 0.0,
    fbm: 3.0, sizingUnit: '-', Amin: 0, Amax: 1e9,
  },
};

/** Map unit operation type to cost correlation key */
function getCostKey(type: string): string {
  const map: Record<string, string> = {
    Heater: 'firedHeater',
    Flash: 'vessel',
    ThreePhaseFlash: 'vessel',
    Mixer: 'mixer',
    Splitter: 'mixer',
    Valve: 'valve',
    Pump: 'pump',
    Compressor: 'compressor',
    HeatX: 'heatExchanger',
    MHeatX: 'heatExchanger',
    Pipe: 'valve',
    RStoic: 'reactor',
    CSTR: 'reactor',
    PFR: 'reactor',
    REquil: 'reactor',
    RGibbs: 'reactor',
    RYield: 'reactor',
    RBatch: 'reactor',
    Column: 'column',
    RadFrac: 'column',
    RateFrac: 'column',
    Absorber: 'column',
    Extractor: 'column',
    Decanter: 'vessel',
    ComponentSeparator: 'vessel',
    Crystallizer: 'reactor',
    Crusher: 'generic',
    Dryer: 'generic',
    Membrane: 'vessel',
    Cyclone: 'vessel',
    Filter: 'vessel',
    Screen: 'vessel',
    Centrifuge: 'vessel',
    SteamDrum: 'vessel',
    SteamHeater: 'heatExchanger',
    SteamTurbine: 'compressor',
    SteamValve: 'valve',
    SteamHeader: 'mixer',
    SteamTrap: 'valve',
    PIDController: 'valve',
    LeadLag: 'valve',
    DeadTime: 'valve',
    SignalSelector: 'valve',
    Analyzer: 'valve',
    ChargeBalance: 'vessel',
    ElectrolyteEquilibrium: 'reactor',
  };
  return map[type] ?? 'generic';
}

/** Estimate the sizing parameter A for a given unit */
function estimateSizingParam(unit: UnitOperation, streams: Record<string, MaterialStream>): number {
  const costKey = getCostKey(unit.type);
  const duty = Math.abs(unit.duty_W);

  switch (costKey) {
    case 'pump':
    case 'compressor':
      // Use duty in kW
      return Math.max(duty / 1000, 1);
    case 'firedHeater':
      return Math.max(duty / 1000, 1000);
    case 'heatExchanger': {
      // Estimate area from Q = U·A·ΔT, assume U=500 W/m²K, ΔTlm=30K
      const U = 500;
      const dT = 30;
      return Math.max(duty / (U * dT), 10);
    }
    case 'vessel':
    case 'column': {
      // Estimate volume from total flow
      const inletPorts = unit.ports.filter(p => p.type === 'inlet' && p.streamId);
      let totalFlow = 0;
      for (const port of inletPorts) {
        const s = streams[port.streamId!];
        if (s) totalFlow += s.totalFlow_molps;
      }
      // Rough: 1 mol/s ≈ 0.1 m³ vessel volume (5 min residence)
      return Math.max(totalFlow * 0.1, 0.3);
    }
    case 'reactor': {
      const vol = unit.params.volume_m3 ?? unit.params.reactorVolume_m3 ?? 1;
      return Math.max(vol, 0.5);
    }
    default:
      return 1;
  }
}

// ─── Equipment cost calculation ─────────────────────────────────

export interface EquipmentCost {
  unitId: string;
  unitName: string;
  unitType: string;
  sizingParam: number;
  sizingUnit: string;
  purchaseCost: number;       // Cp° (base CEPCI)
  bareModuleCost: number;     // CBM (installed)
  /** The cost correlation used */
  correlationKey: string;
}

export function calculateEquipmentCost(
  unit: UnitOperation,
  streams: Record<string, MaterialStream>,
): EquipmentCost {
  const key = getCostKey(unit.type);
  const corr = COST_DATA[key];
  const A = estimateSizingParam(unit, streams);
  const Aclamped = Math.max(corr.Amin, Math.min(A, corr.Amax));
  const logA = Math.log10(Aclamped);
  const logCp = corr.K1 + corr.K2 * logA + corr.K3 * logA * logA;
  const Cp_base = Math.pow(10, logCp);
  const Cp = Cp_base * CEPCI_RATIO;
  const CBM = Cp * corr.fbm;

  return {
    unitId: unit.id,
    unitName: unit.name,
    unitType: unit.type,
    sizingParam: A,
    sizingUnit: corr.sizingUnit,
    purchaseCost: Cp,
    bareModuleCost: CBM,
    correlationKey: key,
  };
}

// ─── Utility cost calculation ───────────────────────────────────

export interface UtilityCost {
  unitId: string;
  unitName: string;
  duty_W: number;
  utilityType: 'electricity' | 'coolingWater' | 'lpSteam' | 'mpSteam' | 'hpSteam' | 'refrigeration' | 'fuel';
  annualCost: number;        // $/year
  rateUsed: number;          // rate in $/unit
}

/** Classify a unit's duty into a utility type */
function classifyUtility(unit: UnitOperation): UtilityCost['utilityType'] {
  const type = unit.type;
  const duty = unit.duty_W;

  // Pumps and compressors use electricity
  if (type === 'Pump' || type === 'Compressor') return 'electricity';

  // Positive duty = heating needed
  if (duty > 0) {
    // Check heater temperature to select steam grade
    const targetT = unit.params.targetT_K ?? unit.params.T_K ?? 373;
    if (targetT > 523) return 'hpSteam';      // > 250°C
    if (targetT > 423) return 'mpSteam';       // > 150°C
    return 'lpSteam';
  }

  // Negative duty = cooling needed
  if (duty < 0) {
    const targetT = unit.params.targetT_K ?? unit.params.T_K ?? 298;
    if (targetT < 273) return 'refrigeration'; // sub-zero
    return 'coolingWater';
  }

  return 'coolingWater';
}

export function calculateUtilityCost(
  unit: UnitOperation,
  rates: UtilityRates,
): UtilityCost | null {
  const duty = unit.duty_W;
  if (Math.abs(duty) < 1) return null; // negligible duty

  const utilityType = classifyUtility(unit);
  let annualCost: number;
  let rate: number;

  if (utilityType === 'electricity') {
    // duty in W → kWh/year
    const kWh_year = (Math.abs(duty) / 1000) * HOURS_PER_YEAR;
    rate = rates.electricity;
    annualCost = kWh_year * rate;
  } else {
    // duty in W → GJ/year
    const GJ_year = (Math.abs(duty) / 1e9) * 3600 * HOURS_PER_YEAR;
    rate = rates[utilityType];
    annualCost = GJ_year * rate;
  }

  return {
    unitId: unit.id,
    unitName: unit.name,
    duty_W: duty,
    utilityType,
    annualCost,
    rateUsed: rate,
  };
}

// ─── Full economics summary ─────────────────────────────────────

export interface EconomicsSummary {
  /** Per-unit equipment costs */
  equipmentCosts: EquipmentCost[];
  /** Per-unit utility costs */
  utilityCosts: UtilityCost[];
  /** Total equipment purchase cost ($) */
  totalPurchaseCost: number;
  /** Total bare module (installed) cost ($) */
  totalBareModuleCost: number;
  /** Total grass-roots capital (CBM × 1.18 for contingency + fees) */
  totalGrassRootsCost: number;
  /** Total annual utility operating cost ($/year) */
  totalAnnualUtilityCost: number;
  /** CEPCI used */
  cepci: number;
  /** Operating hours/year */
  operatingHours: number;
}

export function calculateEconomics(
  units: Record<string, UnitOperation>,
  streams: Record<string, MaterialStream>,
  rates: UtilityRates = DEFAULT_UTILITY_RATES,
): EconomicsSummary {
  const equipmentCosts: EquipmentCost[] = [];
  const utilityCosts: UtilityCost[] = [];

  for (const unit of Object.values(units)) {
    if (!unit.solved) continue;

    equipmentCosts.push(calculateEquipmentCost(unit, streams));

    const uc = calculateUtilityCost(unit, rates);
    if (uc) utilityCosts.push(uc);
  }

  const totalPurchaseCost = equipmentCosts.reduce((s, c) => s + c.purchaseCost, 0);
  const totalBareModuleCost = equipmentCosts.reduce((s, c) => s + c.bareModuleCost, 0);
  const totalGrassRootsCost = totalBareModuleCost * 1.18; // 18% contingency + fees
  const totalAnnualUtilityCost = utilityCosts.reduce((s, c) => s + c.annualCost, 0);

  return {
    equipmentCosts,
    utilityCosts,
    totalPurchaseCost,
    totalBareModuleCost,
    totalGrassRootsCost,
    totalAnnualUtilityCost,
    cepci: CEPCI_CURRENT,
    operatingHours: HOURS_PER_YEAR,
  };
}
