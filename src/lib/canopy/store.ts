// ─── Canopy Zustand Store ────────────────────────────────────────────────────
// Manages all simulation state: compounds, streams, units, and solver results.
// ─────────────────────────────────────────────────────────────────────────────

import { create } from 'zustand';
import type {
  CanopyCompound,
  FluidPackage,
  MaterialStream,
  UnitOperation,
  UnitOpType,
  UnitPort,
  DesignSpec,
  SensitivityAnalysis,
  CalculatorBlock,
} from './types';
import { createDefaultStream } from './types';
import { flashPT } from './thermo';
import type { NRTLParams, WilsonParams, InteractionParams } from './thermo';
import {
  solveMixer, solveHeater, solveFlashDrum, solveValve, solveSplitter,
  solvePump, solveCompressor, solveHeatX, solvePipe,
  solveRStoic, solveColumn, solveColumnRigorous, solveRateFrac,
  solveComponentSeparator, solveCSTR, solvePFR, solveThreePhaseFlash, solveAbsorber, solveExtractor,
  solveRYield, solveDecanter, solveREquil, solveRGibbs,
  solveMHeatX, solveRBatch, solveCrystallizer, solveCrusher, solveDryer,
  solveMembrane, solveCyclone, solveFilter, solveScreen, solveCentrifuge,
  solveSteamDrum, solveSteamHeater, solveSteamTurbine, solveSteamValve, solveSteamHeader, solveSteamTrap,
  solvePIDControllerSignal, solveLeadLagSignal, solveDeadTimeSignal, solveSignalSelector, getRigorousColumnSpecDiagnostics,
  solveAnalyzer, solveChargeBalance, solveElectrolyteEquilibrium,
} from './solver';
import type { ThermoOptions } from './solver';

// ────────────────────────────────────────────────────────────────
// Store Interface
// ────────────────────────────────────────────────────────────────

export type SimScreen = 'setup' | 'flowsheet';

/** Unit system for display (all internal calcs use SI) */
export interface UnitSystem {
  temperature: 'C' | 'K' | 'F';
  pressure: 'bar' | 'atm' | 'kPa' | 'Pa' | 'psi';
  flow: 'mol/s' | 'kmol/h' | 'kg/s' | 'kg/h';
  energy: 'kW' | 'W' | 'BTU/h' | 'kcal/h';
}

export const DEFAULT_UNIT_SYSTEM: UnitSystem = {
  temperature: 'C',
  pressure: 'bar',
  flow: 'mol/s',
  energy: 'kW',
};

/** Tab can be a unit or a stream */
export type TabItem = { type: 'unit'; id: string } | { type: 'stream'; id: string };

interface CanopyState {
  // Navigation
  currentScreen: SimScreen;
  setScreen: (screen: SimScreen) => void;

  // Setup
  compounds: CanopyCompound[];
  fluidPackage: FluidPackage;
  unitSystem: UnitSystem;
  /** Interaction parameters for non-ideal fluid packages */
  interactionParams: InteractionParams;
  setCompounds: (compounds: CanopyCompound[]) => void;
  addCompound: (compound: CanopyCompound) => void;
  removeCompound: (name: string) => void;
  setFluidPackage: (pkg: FluidPackage) => void;
  setUnitSystem: (patch: Partial<UnitSystem>) => void;
  setInteractionParams: (params: Partial<InteractionParams>) => void;

  // Flowsheet
  streams: Record<string, MaterialStream>;
  units: Record<string, UnitOperation>;
  addStream: (stream: MaterialStream) => void;
  updateStream: (id: string, patch: Partial<MaterialStream>) => void;
  removeStream: (id: string) => void;
  addUnit: (unit: UnitOperation) => void;
  updateUnit: (id: string, patch: Partial<UnitOperation>) => void;
  removeUnit: (id: string) => void;
  /** Add a new unit operation to the flowsheet from the object palette */
  addUnitFromPalette: (type: UnitOpType, position: { x: number; y: number }) => void;
  /** Add a new feed stream to the flowsheet */
  addFeedStream: (position?: { x: number; y: number }) => void;
  /** Get the next available ID for a stream or unit */
  nextStreamId: () => string;
  nextUnitId: () => string;
  /** Rename a unit operation */
  renameUnit: (id: string, newName: string) => void;
  /** Rename a stream */
  renameStream: (id: string, newName: string) => void;
  /** Disconnect a stream from its connected units (removes the linking, keeps the stream) */
  disconnectStream: (streamId: string) => void;
  /** Connect two units by creating a stream between an outlet port and an inlet port */
  connectPorts: (sourceUnitId: string, sourceHandleId: string, targetUnitId: string, targetHandleId: string) => void;

  // Tabs (open unit/stream details)
  openTabs: TabItem[];
  activeTab: TabItem | null;
  openTab: (item: TabItem) => void;
  closeTab: (item: TabItem) => void;
  setActiveTab: (item: TabItem | null) => void;

  // Solver
  solved: boolean;
  solverErrors: string[];
  solveAll: () => void;
  /** When true, solver does not auto-run on changes */
  solverPaused: boolean;
  setSolverPaused: (paused: boolean) => void;

  // Warnings (DB failures, missing data, etc.)
  dbWarnings: string[];
  addDbWarning: (msg: string) => void;
  clearDbWarnings: () => void;

  // Design Specs
  designSpecs: DesignSpec[];
  addDesignSpec: (spec: DesignSpec) => void;
  updateDesignSpec: (id: string, patch: Partial<DesignSpec>) => void;
  removeDesignSpec: (id: string) => void;

  // Sensitivity Analysis
  sensitivities: SensitivityAnalysis[];
  addSensitivity: (sa: SensitivityAnalysis) => void;
  removeSensitivity: (id: string) => void;
  runSensitivity: (id: string) => void;

  // Calculator Blocks
  calculatorBlocks: CalculatorBlock[];
  addCalculatorBlock: (cb: CalculatorBlock) => void;
  removeCalculatorBlock: (id: string) => void;
  updateCalculatorBlock: (id: string, patch: Partial<CalculatorBlock>) => void;

  // Preload
  loadPreset: (preset: 'cumene-production') => void;

  // Export / Import
  exportFlowsheet: () => string;
  importFlowsheet: (json: string) => boolean;

  // Undo / Redo
  undo: () => void;
  redo: () => void;
  canUndo: boolean;
  canRedo: boolean;

  // Node positions (for React Flow drag persistence)
  nodePositions: Record<string, { x: number; y: number }>;
  setNodePosition: (nodeId: string, pos: { x: number; y: number }) => void;
}

// ────────────────────────────────────────────────────────────────
// Preloaded Data — Cumene Production (Friedel-Crafts Alkylation)
// ────────────────────────────────────────────────────────────────
//
// All constants from Aspen Plus PURE25 databank (APV140).
// Compounds: PROPYLENE, PROPANE (inert), BENZENE, CUMENE (ISOPROPYLBENZENE), P-DIISOPROPYLBENZENE (DIPB)
//
// Reactions:
//   (1) C₆H₆ + C₃H₆  → C₉H₁₂    (benzene + propylene → cumene)
//   (2) C₉H₁₂ + C₃H₆ → C₁₂H₁₈   (cumene + propylene → p-DIPB, undesired)
//

const PROPYLENE: CanopyCompound = {
  name: 'PROPYLENE',
  displayName: 'Propylene',
  molecularWeight: 42.0797,
  Tc_K: 364.85,
  Pc_Pa: 4600000,
  omega: 0.137588,
  Vc_m3pkmol: 0.185,
  uniquac_r: 2.24654,
  uniquac_q: 2.024,
  rackett_ZRA: 0.27753,    // APV140 PURE40
  Zc: 0.281,
  Vb_m3pmol: 6.88009e-5,   // 0.0688009 m³/kmol
  Tf_K: 87.9,
  Hfus_Jmol: 2936,         // 2936000 J/kmol
  Hcomb_Jmol: -1926200,    // -1926200000 J/kmol
  Hvap_nb_Jmol: 18731.7,   // 18731700 J/kmol at Tb
  VLSTD_m3pkmol: 0.0808566,
  mathiasCopeman: { C1: 0.587594, C2: -0.0686087, C3: 0.36052 },
  mathiasCopeman_SRK: { C1: 0.737737, C2: -0.347116, C3: 0.67512 },
  volumeTranslation_PR: -0.00446284,  // m³/kmol
  Hf_298_Jmol: 20230,   // J/mol (PURE25: 20230000 J/kmol)
  Gf_298_Jmol: 62640,   // J/mol (PURE25: 62640000 J/kmol)
  delta_Jm3_05: 12702.6,  // Hildebrand solubility parameter (J/m³)^0.5
  radiusOfGyration_m: 2.2311e-10,  // m
  pcSaft: { m: 1.9597, sigma: 3.5356, epsilon_k: 207.19 },
  cpa: { Tc_K: 369.031, Pc_Pa: 5025191, m: 0.68414, epsilon_Pa_m3: 0, beta: 0 },
  dielectricConst: { e1: 1.739026, e2: 371.95, e3_K: 298.15 },
  cosmoVolume: 75.60203,
  solidLiquidCpDiff_JpkmolK: 29301.9,
  Tb_K: 225.45,
  unifacGroups: { 1: 1, 5: 1 },  // 1×CH3 (subgroup 1) + 1×CH2=CH (subgroup 5)
  dipprCoeffs: {
    LiquidHeatCapacityCp: { A: 114140, B: -343.72, C: 1.0905, D: 0, E: 0, eqNo: 100 },
    HeatOfVaporization: { A: 25216000, B: 0.33721, C: -0.18399, D: 0.22377, E: 0 },
    IdealGasHeatCapacityCp: { A: 0, B: 0, C: 0, D: 0, E: 0 },
    LiquidDensity: { A: 1.4403, B: 0.26852, C: 364.85, D: 0.28775, E: 0 },
    SecondVirialCoefficient: { A: 0.11727, B: -94.686, C: -4033300, D: 2.0649e18, E: -4.9773e20 },
    SolidDensity: { A: 20.926, B: -0.021754, C: 0, D: 0, E: 0, eqNo: 100 },
    SolidHeatCapacity: { A: -217.22, B: -134.39, C: 46.474, D: -0.79918, E: 0.0043349, eqNo: 100 },
  },
  cpigdp: { A: 43852, B: 150600, C: 1398.8, D: 74754, E: 616.46, Tmin_K: 130, Tmax_K: 1500 },
  antoine: {
    C1: 43.905, C2: -3097.8, C3: 0, C4: 0,
    C5: -3.4425, C6: 9.9989e-17, C7: 6,
    Tmin_K: 87.89, Tmax_K: 364.85,
  },
  transportCoeffs: {
    LiquidViscosity: { A: -92.082, B: 1907.3, C: 15.639, D: -0.043098, E: 1 },
    VaporViscosity: { A: 7.3919e-7, B: 0.5423, C: 263.73, D: 0 },
    LiquidThermalConductivity: { A: 0.099874, B: -1.3409, C: 1.6545, D: 1.3334 },
    VaporThermalConductivity: { A: 4.49e-5, B: 1.2018, C: 421, D: 0 },
    SurfaceTension: { A: 0.053118, B: 1.1993, C: 0, D: 0 },
  },
  specificGravity: 0.522,
  flashPoint_K: 169,
  autoignitionTemp_K: 728.15,
  THRSWT: [100, 105, 101, 106, 100, 100, 107, 104],
  TRNSWT: [101, 102, 123, 102, 106],
  psxant: { C1: 31.977, C2: -3403.7, C3: 0, C4: 0, C5: 0, C6: 0, C7: 0, Tmin_K: 37.89, Tmax_K: 87.89 },
  waterSolubility: { A: 3.47154, B: -3028.13, C: 0, Tmin_K: 295, Tmax_K: 359 },
  schwartentruberRenon_PR: { p1: -0.0658764, p2: 0.170939, p3: -0.110312 },
  schwartentruberRenon_SRK: { p1: -0.041136, p2: 0.159201, p3: -0.142792 },
  VLCVT1: 0.00969,
  omegaCostald: 0.34922,
  twuAlpha: { L: 0.807472, M: 0.930959, N: 0.693585 },
  VSTCTD: 0.191,
  dhForm_Jpkmol: 20230000,
  dgForm_Jpkmol: 62640000,
  entropy_JpkmolK: 267000,
  hCombustion_Jpkmol: -1926200000,
  hFusion_Jpkmol: 2936000,
  dhVapAtNBP_Jpkmol: 18731700,
  cpDepartureLiq_JpkmolK: 29301.9,
  vLiqStd_m3pkmol: 0.0808566,
  vLiqAtNBP_m3pkmol: 0.0688009,
  dipoleMoment_Cm: 1.1566e-25,
  refractiveIndex: 1.305,
  flammabilityLower_pct: 2.15,
  flammabilityUpper_pct: 11.2,
  triplePointT_K: 87.89,
  triplePointP_Pa: 0.00117,
};

const PROPANE: CanopyCompound = {
  name: 'PROPANE',
  displayName: 'Propane',
  molecularWeight: 44.0956,
  Tc_K: 369.83,
  Pc_Pa: 4248000,
  omega: 0.152291,
  Vc_m3pkmol: 0.200,
  uniquac_r: 2.4766,
  uniquac_q: 2.236,
  rackett_ZRA: 0.27657,    // APV140 PURE40
  Zc: 0.276,
  Vb_m3pmol: 7.56892e-5,   // 0.0756892 m³/kmol
  Tf_K: 85.47,
  Hfus_Jmol: 3524,         // 3524000 J/kmol
  Hcomb_Jmol: -2043110,    // -2043110000 J/kmol
  Hvap_nb_Jmol: 18746,     // 18746000 J/kmol at Tb
  VLSTD_m3pkmol: 0.0871442,
  mathiasCopeman: { C1: 0.610032, C2: -0.096182, C3: 0.357249 },
  mathiasCopeman_SRK: { C1: 0.76154, C2: -0.376913, C3: 0.671243 },
  volumeTranslation_PR: -0.00481468,  // m³/kmol
  Hf_298_Jmol: -104680,   // J/mol
  Gf_298_Jmol: -24390,    // J/mol
  delta_Jm3_05: 12281.1,  // Hildebrand solubility parameter (J/m³)^0.5
  radiusOfGyration_m: 2.431e-10,  // m
  pcSaft: { m: 2.002, sigma: 3.6184, epsilon_k: 208.11 },
  cpa: { Tc_K: 378.635, Pc_Pa: 4716099, m: 0.64296, epsilon_Pa_m3: 0, beta: 0 },
  brelviOConnellV: 0.23695,
  cosmoVolume: 80.70296,
  solidLiquidCpDiff_JpkmolK: 31887.1,
  Tb_K: 231.11,
  unifacGroups: { 1: 2, 2: 1 },  // 2×CH3 (subgroup 1) + 1×CH2 (subgroup 2)
  dipprCoeffs: {
    LiquidHeatCapacityCp: { A: 62.983, B: 113630, C: 633.21, D: -873.46, E: 0, eqNo: 114 },
    HeatOfVaporization: { A: 29209000, B: 0.78237, C: -0.77319, D: 0.39246, E: 0 },
    IdealGasHeatCapacityCp: { A: 0, B: 0, C: 0, D: 0, E: 0 },
    LiquidDensity: { A: 1.3757, B: 0.27453, C: 369.83, D: 0.29359, E: 0 },
    SecondVirialCoefficient: { A: 0.1127, B: -99.2, C: -4510000, D: 3.09e17, E: -7.05e19 },
    SolidDensity: { A: 18.861, B: -0.020332, C: 0, D: 0, E: 0, eqNo: 100 },
    SolidHeatCapacity: { A: -11230, B: 1059, C: -3.6, D: 0, E: 0, eqNo: 100 },
  },
  cpigdp: { A: 59474, B: 126610, C: 844.31, D: 86165, E: 2482.7, Tmin_K: 298.15, Tmax_K: 1500 },
  antoine: {
    C1: 59.078, C2: -3492.6, C3: 0, C4: 0,
    C5: -6.0669, C6: 1.0919e-5, C7: 2,
    Tmin_K: 85.47, Tmax_K: 369.83,
  },
  transportCoeffs: {
    LiquidViscosity: { A: -17.156, B: 646.25, C: 1.1101, D: -7.3439e-11, E: 4 },
    VaporViscosity: { A: 4.9054e-8, B: 0.90125, C: 0, D: 0 },
    LiquidThermalConductivity: { A: 0.34489, B: -3.9485, C: 5.8521, D: -2.1706 },
    VaporThermalConductivity: { A: -1.12, B: 0.10972, C: -9834.6, D: -7535800 },
    SurfaceTension: { A: 0.05092, B: 1.2197, C: 0, D: 0 },
  },
  specificGravity: 0.5077,
  flashPoint_K: 171,
  autoignitionTemp_K: 723,
  THRSWT: [100, 105, 101, 106, 100, 114, 107, 104],
  TRNSWT: [101, 102, 123, 102, 106],
  psxant: { C1: 31.429, C2: -3429.2, C3: 0, C4: 0, C5: 0, C6: 0, C7: 0, Tmin_K: 35.47, Tmax_K: 85.47 },
  waterSolubility: { A: 8.37463, B: -4893.61, C: 0, Tmin_K: 285, Tmax_K: 366 },
  schwartentruberRenon_PR: { p1: -0.0557464, p2: 0.152839, p3: -0.104749 },
  schwartentruberRenon_SRK: { p1: -0.0307137, p2: 0.140524, p3: -0.13752 },
  VLCVT1: 0.010332,
  omegaCostald: 0.21454,
  twuAlpha: { L: 0.585635, M: 0.878009, N: 0.926249 },
  VSTCTD: 0.20277,
  dhForm_Jpkmol: -104680000,
  dgForm_Jpkmol: -24390000,
  entropy_JpkmolK: 270200,
  hCombustion_Jpkmol: -2043110000,
  hFusion_Jpkmol: 3524000,
  dhVapAtNBP_Jpkmol: 18746000,
  cpDepartureLiq_JpkmolK: 31887.1,
  vLiqStd_m3pkmol: 0.0871442,
  vLiqAtNBP_m3pkmol: 0.0756892,
  dipoleMoment_Cm: 0,
  refractiveIndex: 1.28614,
  flammabilityLower_pct: 2.1,
  flammabilityUpper_pct: 9.5,
  triplePointT_K: 85.47,
  triplePointP_Pa: 0.0001685,
};

const BENZENE: CanopyCompound = {
  name: 'BENZENE',
  displayName: 'Benzene',
  molecularWeight: 78.1118,
  Tc_K: 562.05,
  Pc_Pa: 4895000,
  omega: 0.2103,
  Vc_m3pkmol: 0.256,
  uniquac_r: 3.19051,
  uniquac_q: 2.4,
  rackett_ZRA: 0.2697,     // APV140 PURE40
  Zc: 0.268,
  Vb_m3pmol: 9.58294e-5,   // 0.0958294 m³/kmol
  Tf_K: 278.68,
  Hfus_Jmol: 9866,         // 9866000 J/kmol
  Hcomb_Jmol: -3136000,    // -3136000000 J/kmol
  Hvap_nb_Jmol: 30736.9,   // 30736900 J/kmol at Tb
  VLSTD_m3pkmol: 0.0885091,
  mathiasCopeman: { C1: 0.700636, C2: -0.254116, C3: 0.986886 },
  mathiasCopeman_SRK: { C1: 0.860821, C2: -0.577298, C3: 1.40125 },
  volumeTranslation_PR: -0.0029029,  // m³/kmol
  Hf_298_Jmol: 82880,   // J/mol
  Gf_298_Jmol: 129600,  // J/mol
  delta_Jm3_05: 18730,  // Hildebrand solubility parameter (J/m³)^0.5
  radiusOfGyration_m: 3.004e-10,  // m
  pcSaft: { m: 2.4653, sigma: 3.6478, epsilon_k: 287.35 },
  cpa: { Tc_K: 572.814, Pc_Pa: 5502449, m: 0.77025, epsilon_Pa_m3: 0, beta: 0 },
  brelviOConnellV: 0.35,
  dielectricConst: { e1: 2.27385, e2: 244.744, e3_K: 298.15 },
  cosmoVolume: 110.22176,
  solidLiquidCpDiff_JpkmolK: 1337.33,
  hocEta: { 'PROPYLENE': 0.1 },
  Tb_K: 353.24,
  unifacGroups: { 9: 6 },  // 6×ACH (subgroup 9, main group 3)
  dipprCoeffs: {
    IdealGasHeatCapacityCp: { A: 34010.24, B: -588.0978, C: 12.81777, D: -0.000197306, E: 5.142899e-8 },
    HeatOfVaporization: { A: 50007000, B: 0.65393, C: -0.27698, D: 0.029569, E: 0 },
    LiquidHeatCapacityCp: { A: 129440, B: -169.5, C: 0.64781, D: 0, E: 0, eqNo: 100 },
    LiquidDensity: { A: 1.0259, B: 0.26666, C: 562.05, D: 0.28394, E: 0 },
    SecondVirialCoefficient: { A: 0.15059, B: -186.94, C: -23146000, D: -7.0493e18, E: -6.8786e20 },
    SolidDensity: { A: 13.061, B: -0.00035714, C: 0, D: 0, E: 0, eqNo: 100 },
    SolidHeatCapacity: { A: 7400, B: 624.9, C: -2.6874, D: 0.007316, E: 0, eqNo: 100 },
  },
  cpigdp: { A: 55238, B: 173380, C: 764.25, D: 72545, E: 2445.7, Tmin_K: 298.15, Tmax_K: 1500 },
  antoine: {
    C1: 83.107, C2: -6486.2, C3: 0, C4: 0,
    C5: -9.2194, C6: 6.9844e-6, C7: 2,
    Tmin_K: 278.68, Tmax_K: 562.05,
  },
  transportCoeffs: {
    LiquidViscosity: { A: 7.5117, B: 294.68, C: -2.794, D: 0 },
    VaporViscosity: { A: 3.134e-8, B: 0.9676, C: 7.9, D: 0 },
    LiquidThermalConductivity: { A: 0.054252, B: 2.7419, C: -7.2256, D: 8.2256 },
    VaporThermalConductivity: { A: 1.652e-5, B: 1.3117, C: 491, D: 0 },
    SurfaceTension: { A: 0.071815, B: 1.2362, C: 0, D: 0 },
  },
  specificGravity: 0.8844,
  flashPoint_K: 262,
  autoignitionTemp_K: 833.15,
  THRSWT: [100, 105, 101, 106, 100, 100, 107, 104],
  TRNSWT: [101, 102, 123, 102, 106],
  psxant: { C1: 72.829, C2: -7042.3, C3: 0, C4: 0, C5: -7.061, C6: 8.6915e-6, C7: 2, Tmin_K: 178.25, Tmax_K: 278.68 },
  waterSolubility: { A: 4.8818, B: -3181.37, C: 0, Tmin_K: 275, Tmax_K: 533 },
  hcSolubility: { A: -192.4824, B: 8053.106, C: 27.67472, Tmin_K: 283.15, Tmax_K: 700 },
  schwartentruberRenon_PR: { p1: -0.138563, p2: 0.374156, p3: -0.24971 },
  schwartentruberRenon_SRK: { p1: -0.120911, p2: 0.382863, p3: -0.299139 },
  VLCVT1: 0.01226,
  omegaCostald: 0.33714,
  twuAlpha: { L: 0.085723, M: 0.871128, N: 3.39897 },
  VSTCTD: 0.2631,
  dhForm_Jpkmol: 82880000,
  dgForm_Jpkmol: 129600000,
  entropy_JpkmolK: 269300,
  hCombustion_Jpkmol: -3136000000,
  hFusion_Jpkmol: 9866000,
  dhVapAtNBP_Jpkmol: 30736900,
  cpDepartureLiq_JpkmolK: 1337.33,
  vLiqStd_m3pkmol: 0.0885091,
  vLiqAtNBP_m3pkmol: 0.0958294,
  dipoleMoment_Cm: 0,
  refractiveIndex: 1.49792,
  flammabilityLower_pct: 1.2,
  flammabilityUpper_pct: 8,
  triplePointT_K: 278.68,
  triplePointP_Pa: 4764.22,
};

const CUMENE: CanopyCompound = {
  name: 'CUMENE',
  displayName: 'Cumene',
  molecularWeight: 120.192,
  Tc_K: 631,
  Pc_Pa: 3209000,
  omega: 0.327406,
  Vc_m3pkmol: 0.434,
  uniquac_r: 5.27093,
  uniquac_q: 4.056,
  rackett_ZRA: 0.26176,    // APV140 PURE40
  Zc: 0.265,
  Vb_m3pmol: 1.61892e-4,   // 0.161892 m³/kmol
  Tf_K: 177.14,
  Hfus_Jmol: 7326,         // 7326000 J/kmol
  Hcomb_Jmol: -4951000,    // -4951000000 J/kmol
  Hvap_nb_Jmol: 37118.5,   // 37118500 J/kmol at Tb
  VLSTD_m3pkmol: 0.139057,
  mathiasCopeman: { C1: 0.8578, C2: -0.179221, C3: 0.741113 },
  mathiasCopeman_SRK: { C1: 1.02898, C2: -0.515535, C3: 1.15396 },
  volumeTranslation_PR: 0.00012767,  // m³/kmol
  Hf_298_Jmol: 4000,    // J/mol (PURE25: 4000000 J/kmol)
  Gf_298_Jmol: 137900,  // J/mol
  delta_Jm3_05: 17490,  // Hildebrand solubility parameter (J/m³)^0.5
  radiusOfGyration_m: 4.322e-10,  // m
  cpa: { Tc_K: 631.15, Pc_Pa: 3208810, m: 0.99147, epsilon_Pa_m3: 0, beta: 0 },
  dielectricConst: { e1: 2.38, e2: 0, e3_K: 293.15 },
  cosmoVolume: 173.99542,
  solidLiquidCpDiff_JpkmolK: 18826.5,
  anilinePoint_K: 258.15,
  Tb_K: 425.56,
  unifacGroups: { 1: 2, 9: 5, 13: 1 },  // 2×CH3 + 5×ACH + 1×ACCH (subgroup 13, main group 4)
  dipprCoeffs: {
    IdealGasHeatCapacityCp: { A: 0, B: 0, C: 0, D: 0, E: 0 },
    HeatOfVaporization: { A: 75255000, B: 1.3714, C: -1.5024, D: 0.59731, E: 0 },
    LiquidHeatCapacityCp: { A: 61723, B: 494.81, C: 0, D: 0, E: 0, eqNo: 100 },
    LiquidDensity: { A: 0.58711, B: 0.25583, C: 631, D: 0.28498, E: 0 },
    SecondVirialCoefficient: { A: 0.25065, B: -340.58, C: -78182000, D: -2.3663e20, E: 3.1813e22 },
    SolidDensity: { A: 10.162, B: -0.0071706, C: 0, D: 0, E: 0, eqNo: 100 },
    SolidHeatCapacity: { A: -10000, B: 1493.8, C: -8.879, D: 0.027804, E: 0, eqNo: 100 },
  },
  cpigdp: { A: 108100, B: 379320, C: 1750.5, D: 300270, E: 794.8, Tmin_K: 200, Tmax_K: 1500 },
  antoine: {
    C1: 102.81, C2: -8674.6, C3: 0, C4: 0,
    C5: -11.922, C6: 7.0048e-6, C7: 2,
    Tmin_K: 177.14, Tmax_K: 631,
  },
  transportCoeffs: {
    LiquidViscosity: { A: -24.988, B: 1807.9, C: 2.0556, D: 0 },
    VaporViscosity: { A: 3.3699e-7, B: 0.60751, C: 221.17, D: 0 },
    LiquidThermalConductivity: { A: 0.06254, B: 1.0399, C: -3.3579, D: 4.3921 },
    VaporThermalConductivity: { A: 1.6743e-7, B: 1.8369, C: -449.46, D: 112760 },
    SurfaceTension: { A: 0.063695, B: 1.3027, C: 0, D: 0 },
  },
  specificGravity: 0.8663,
  flashPoint_K: 309.15,
  autoignitionTemp_K: 697,
  THRSWT: [100, 105, 101, 106, 100, 100, 107, 104],
  TRNSWT: [101, 102, 123, 102, 106],
  psxant: { C1: 33.893, C2: -7360.7, C3: 0, C4: 0, C5: 0, C6: 0, C7: 0, Tmin_K: 127.14, Tmax_K: 177.14 },
  waterSolubility: { A: 2.96, B: -2685.3, C: 0, Tmin_K: 273, Tmax_K: 323 },
  hcSolubility: { A: -304.679, B: 12774.71, C: 43.8994, Tmin_K: 283.15, Tmax_K: 700 },
  schwartentruberRenon_PR: { p1: -0.116462, p2: 0.315099, p3: -0.211926 },
  schwartentruberRenon_SRK: { p1: -0.0961385, p2: 0.320417, p3: -0.262758 },
  VLCVT1: 0.020948,
  omegaCostald: 0.43808,
  twuAlpha: { L: 0.374974, M: 0.800115, N: 1.657 },
  VSTCTD: 0.43932,
  dhForm_Jpkmol: 4000000,
  dgForm_Jpkmol: 137900000,
  entropy_JpkmolK: 386000,
  hCombustion_Jpkmol: -4951000000,
  hFusion_Jpkmol: 7326000,
  dhVapAtNBP_Jpkmol: 37118500,
  cpDepartureLiq_JpkmolK: 18826.5,
  vLiqStd_m3pkmol: 0.139057,
  vLiqAtNBP_m3pkmol: 0.161892,
  dipoleMoment_Cm: 1.23244e-25,
  refractiveIndex: 1.4889,
  flammabilityLower_pct: 0.88,
  flammabilityUpper_pct: 6.5,
  triplePointT_K: 177.14,
  triplePointP_Pa: 0.000471313,
};

const P_DIISOPROPYLBENZENE: CanopyCompound = {
  name: 'P-DIISOPROPYLBENZENE',
  displayName: 'p-Diisopropylbenzene',
  molecularWeight: 162.271,
  Tc_K: 689,
  Pc_Pa: 2450000,
  omega: 0.390023,
  Vc_m3pkmol: 0.598,
  uniquac_r: 6.65788,
  uniquac_q: 5.616,
  rackett_ZRA: 0.25776,    // APV140 PURE40
  Zc: 0.256,
  Vb_m3pmol: 2.30889e-4,   // 0.230889 m³/kmol
  Tf_K: 256.08,
  Hcomb_Jmol: -6770000,    // -6770000000 J/kmol
  Hvap_nb_Jmol: 43219.6,   // 43219600 J/kmol at Tb
  VLSTD_m3pkmol: 0.190236,
  mathiasCopeman: { C1: 0.851522, C2: 0.5151, C3: 0.0699018 },
  mathiasCopeman_SRK: { C1: 1.02094, C2: 0.230081, C3: 0.400452 },
  volumeTranslation_PR: 0.00392245,  // m³/kmol
  Hf_298_Jmol: -77600,  // J/mol (PURE25: -77600000 J/kmol)
  Gf_298_Jmol: 147804,  // J/mol
  delta_Jm3_05: 16880,  // Hildebrand solubility parameter (J/m³)^0.5
  radiusOfGyration_m: 5.178e-10,  // m
  Tb_K: 483.65,
  unifacGroups: { 1: 4, 9: 4, 13: 2 },  // 4×CH3 + 4×ACH + 2×ACCH
  dipprCoeffs: {
    IdealGasHeatCapacityCp: { A: 0, B: 0, C: 0, D: 0, E: 0 },
    HeatOfVaporization: { A: 71697000, B: 0.46446, C: -0.19741, D: 0.1872, E: 0 },
    LiquidHeatCapacityCp: { A: 109600, B: 608.4, C: 0, D: 0, E: 0, eqNo: 100 },
    LiquidDensity: { A: 0.42768, B: 0.25774, C: 689, D: 0.2857, E: 0 },
    SecondVirialCoefficient: { A: 0.34302, B: -577.4, C: -146990000, D: -1.4104e21, E: 3.3952e23 },
    SolidDensity: { A: 6.9624, B: -0.0033986, C: 0, D: 0, E: 0, eqNo: 100 },
    SolidHeatCapacity: { A: 2980, B: 0.79267, C: 0, D: 0, E: 0, eqNo: 100 },
  },
  cpigdp: { A: 157150, B: 534060, C: 1666.5, D: 394600, E: 760.64, Tmin_K: 300, Tmax_K: 1500 },
  antoine: {
    C1: 76.617, C2: -9057, C3: 0, C4: 0,
    C5: -7.5033, C6: 2.5709e-18, C7: 6,
    Tmin_K: 256.08, Tmax_K: 689,
  },
  transportCoeffs: {
    LiquidViscosity: { A: -11.129, B: 1224.8, C: -0.020812, D: 0 },
    VaporViscosity: { A: 8.0509e-7, B: 0.4804, C: 392.06, D: 0 },
    LiquidThermalConductivity: { A: 0.1698, B: -0.0001761, C: 0, D: 0 },
    VaporThermalConductivity: { A: 0.3061, B: -0.06924, C: -303.6, D: 1917500 },
    SurfaceTension: { A: 0.05856, B: 1.2503, C: 0, D: 0 },
  },
  specificGravity: 0.855153,
  flashPoint_K: 349.82,
  autoignitionTemp_K: 722,
  THRSWT: [100, 105, 101, 106, 102, 100, 107, 104],
  TRNSWT: [101, 102, 100, 102, 106],
  psxant: { C1: 34.865, C2: -9020.4, C3: 0, C4: 0, C5: 0, C6: 0, C7: 0, Tmin_K: 206.08, Tmax_K: 256.08 },
  // No WATSOL data in PURE40
  hcSolubility: { A: -429.463, B: 18024.6, C: 61.9402, Tmin_K: 283.15, Tmax_K: 700 },
  schwartentruberRenon_PR: { p1: -0.215409, p2: 0.443208, p3: -0.196636 },
  schwartentruberRenon_SRK: { p1: -0.19338, p2: 0.438702, p3: -0.235485 },
  VLCVT1: 0.0271836,
  omegaCostald: 0.54105,
  twuAlpha: { L: 1.1164, M: 0.853157, N: 0.817496 },
  VSTCTD: 0.63041,
  dhForm_Jpkmol: -77600000,
  dgForm_Jpkmol: 147804000,
  entropy_JpkmolK: 488000,
  hCombustion_Jpkmol: -6770000000,
  dhVapAtNBP_Jpkmol: 43219600,
  cpDepartureLiq_JpkmolK: 23706.4,
  vLiqStd_m3pkmol: 0.190236,
  vLiqAtNBP_m3pkmol: 0.230889,
  dipoleMoment_Cm: 2.18047e-26,
  refractiveIndex: 1.48758,
  flammabilityLower_pct: 0.7,
  flammabilityUpper_pct: 4.3,
  triplePointT_K: 256.08,
  triplePointP_Pa: 0.697953,
  cosmoVolume: 237.22471,
  solidLiquidCpDiff_JpkmolK: 23706.4,
};

// ─── NRTL Binary Interaction Parameters ──────────────────────────────────────
// Merged from APV140 VLE-IG and NIST-IG databanks.
// τ_ij = aij + bij/T + eij·ln(T) + fij·T
// G_ij = exp(-α_ij · τ_ij)   where α_ij = cij
//
// Key: "COMPOUND1|COMPOUND2" (alphabetical order; function handles lookup both ways)
import type { NrtlBinaryParams, EosKij } from './types';

/** The 5 preset compounds for the cumene production flowsheet */
export const PRESET_COMPOUNDS: CanopyCompound[] = [PROPYLENE, PROPANE, BENZENE, CUMENE, P_DIISOPROPYLBENZENE];

export const NRTL_BINARY_PARAMS: Record<string, NrtlBinaryParams> = {
  // ── APV140 VLE-IG pairs ──
  'BENZENE|ISOPROPYLBENZENE': {
    aij: 0, aji: 0, bij: 54.4802, bji: -44.669, cij: 0.3,
    Tlower_K: 354.75, Tupper_K: 409.55, source: 'APV140 VLE-IG',
  },
  'BENZENE|PROPYLENE': {
    aij: 0, aji: 0, bij: -3.7947, bji: 151.4452, cij: 0.3,
    Tlower_K: 298.15, Tupper_K: 298.15, source: 'APV140 VLE-IG',
  },
  // ── NIST-IG pairs (fills gaps not in APV140) ──
  'BENZENE|PROPANE': {
    aij: -0.997877, aji: 0.392697, bij: 301.744, bji: 155.513, cij: 0.220653,
    Tlower_K: 160.014, Tupper_K: 344.274, source: 'NIST-IG',
  },
  'BENZENE|P-DIISOPROPYLBENZENE': {
    // stored as DIPB→BENZENE in NIST; swap to alphabetical: BENZENE is "i", DIPB is "j"
    aij: 0.661867, aji: -1.02131, bij: -45.6104, bji: 143.881, cij: 0.5,
    Tlower_K: 334.8, Tupper_K: 430.8, source: 'NIST-IG (swapped)',
  },
  'ISOPROPYLBENZENE|PROPYLENE': {
    // stored as PROPYLENE→ISOPROPYLBENZENE in NIST; swap to alphabetical
    aij: 8.50646, aji: -3.54578, bij: -2030.88, bji: 867.174, cij: 0.5,
    Tlower_K: 273.15, Tupper_K: 353.123, source: 'NIST-IG (swapped)',
  },
  'PROPANE|PROPYLENE': {
    // stored as PROPYLENE→PROPANE in NIST; swap to alphabetical
    aij: -0.0981501, aji: 0.00125927, bij: 46.9617, bji: 1.66539, cij: 0.499999,
    Tlower_K: 230.008, Tupper_K: 273.15, source: 'NIST-IG (swapped)',
  },
};

// ─── EOS Binary Interaction Parameters (PR / SRK kij) ──────────────────────
// From APV140 PRKBV and RKSKBV (EOS-LIT databank).
// kij(T) = kij1 + kij2·T + kij3·T²

export const PR_KIJ: Record<string, EosKij> = {
  'BENZENE|PROPANE': { kij1: 0.0233, source: 'APV140 EOS-LIT' },
  'PROPANE|PROPYLENE': { kij1: 0.0074, source: 'APV140 EOS-LIT' },
  // HYSYS values for reference:
  // PROPANE+PROPYLENE: 0.0079, PROPANE+BENZENE: 0.019997
};

export const SRK_KIJ: Record<string, EosKij> = {
  'BENZENE|PROPANE': { kij1: 0.02, source: 'APV140 EOS-LIT' },
  'PROPANE|PROPYLENE': { kij1: 0.008, source: 'APV140 EOS-LIT' },
};

// ─── UNIFAC Interaction Parameters (Original UNIFAC) ─────────────────────────
// amn in K.  Ψ_mn = exp(-a_mn / T)
// Source: Aspen APV140, Hansen et al. (1991), Fredenslund et al. (1975)
// Only main groups relevant to cumene system: 1=CH2, 2=C=C, 3=ACH, 4=ACCH2
export const UNIFAC_INTERACTION_PARAMS: Record<string, number> = {
  '1|2': 86.02,   '2|1': -35.36,
  '1|3': 61.13,   '3|1': -11.12,
  '1|4': 76.5,    '4|1': -69.7,
  '2|3': 38.81,   '3|2': 3.446,
  '2|4': 74.15,   '4|2': -113.6,
  '3|4': 167.0,   '4|3': -146.8,
};

// ─── Modified UNIFAC (Dortmund) Interaction Parameters ───────────────────────
// Temperature-dependent: Ψ_mn = exp(-(amn + bmn·T + cmn·T²) / T)
// Source: Aspen APV140 (Gmehling et al.)
export const UNIFAC_DORTMUND_PARAMS: Record<string, { a: number; b: number; c: number }> = {
  '1|2': { a: 189.66, b: -0.2723, c: 0 },    '2|1': { a: -95.42, b: 0.1852, c: 0 },
  '1|3': { a: 114.2,  b: -0.0908, c: 0 },     '3|1': { a: -34.38, b: 0.0425, c: 0 },
  '1|4': { a: 145.2,  b: -0.1951, c: 0 },     '4|1': { a: -41.52, b: -0.0563, c: 0 },
  '2|3': { a: -28.25, b: -0.0174, c: 0 },     '3|2': { a: -8.479, b: -0.0067, c: 0 },
  '3|4': { a: 305.7,  b: -0.855, c: 0.00077 },'4|3': { a: -149.5, b: 0.3345, c: 0 },
  '4|2': { a: -53.46, b: -0.0508, c: 0 },
};

// ─── Wilson Binary Interaction Parameters ────────────────────────────────────
// ln(Λ_ij) = aij + bij/T + cij·ln(T) + dij·T  [Aspen form]
// Source: APV140 VLE-IG
export interface WilsonBinaryParams {
  aij: number; aji: number;
  bij: number; bji: number;  // K
  cij?: number; cji?: number;
  dij?: number; dji?: number;
  Tlower_K?: number; Tupper_K?: number;
}

export const WILSON_BINARY_PARAMS: Record<string, WilsonBinaryParams> = {
  'BENZENE|PROPYLENE': {
    aij: 0, aji: 0,
    bij: -123.2166, bji: -26.1093,
    Tlower_K: 298.15, Tupper_K: 298.15,
  },
  'BENZENE|ISOPROPYLBENZENE': {
    aij: 0, aji: 0,
    bij: 54.6667, bji: -68.1399,
    Tlower_K: 354.75, Tupper_K: 409.55,
  },
};

// ─── UNIQUAC Binary Interaction Parameters (APV140 VLE-IG) ───────────────────
// τ_ij = exp(aij + bij/T + cij·ln(T) + dij·T)
// Key: "COMPOUND1|COMPOUND2" (alphabetical)
export interface UniquacBinaryParams {
  aij: number; aji: number;
  bij: number; bji: number;  // K
  cij?: number; cji?: number;
  dij?: number; dji?: number;
  Tlower_K?: number; Tupper_K?: number;
  source?: string;
}

export const UNIQUAC_BINARY_PARAMS: Record<string, UniquacBinaryParams> = {
  'BENZENE|PROPYLENE': {
    aij: 0, aji: 0,
    bij: 87.5607, bji: -160.6023,
    Tlower_K: 298.15, Tupper_K: 298.15,
    source: 'APV140 VLE-IG',
  },
  'BENZENE|ISOPROPYLBENZENE': {
    aij: 0, aji: 0,
    bij: -105.2781, bji: 74.938,
    Tlower_K: 354.75, Tupper_K: 409.55,
    source: 'APV140 VLE-IG',
  },
};

// ─── Henry's Law Constants (APV140 HENRY-AP) ────────────────────────────────
// ln(H) = A + B/T + C·ln(T) + D·T + E·T²
// H in Pa, T in K. Solute dissolved in the listed solvent.
// Key: "SOLUTE|SOLVENT"
export interface HenryParams {
  A: number; B: number; C: number; D: number; E: number;
  Tmin_K: number; Tmax_K: number;
}

export const HENRY_BINARY_PARAMS: Record<string, HenryParams> = {
  'PROPYLENE|BENZENE': {
    A: 103.945, B: -6196, C: -12.073, D: 0, E: 0,
    Tmin_K: 298.15, Tmax_K: 343.15,
  },
  'PROPANE|BENZENE': {
    A: 28.5346, B: -2200.7, C: -1.1886, D: 0, E: 0,
    Tmin_K: 283.15, Tmax_K: 343.15,
  },
};

// ─── NIST UNIQUAC Binary Params (NIST-IG) ─────────────────────────────────────
// τ_ij = exp(aij + bij/T), same extended form as APV140
export const NIST_UNIQUAC_BINARY_PARAMS: Record<string, UniquacBinaryParams> = {
  'BENZENE|ISOPROPYLBENZENE': {
    aij: 0.0512383, aji: -0.0261228,
    bij: 37.9295, bji: -63.9315,
    Tlower_K: 298.15, Tupper_K: 409.531,
    source: 'NIST-IG',
  },
  'BENZENE|PROPYLENE': {
    aij: 2.6705, aji: -3.4161,
    bij: -727.419, bji: 872.61,
    Tlower_K: 279.99, Tupper_K: 343.132,
    source: 'NIST-IG',
  },
  'BENZENE|PROPANE': {
    aij: 0.430517, aji: -0.0910563,
    bij: -136.3, bji: -58.1696,
    Tlower_K: 160.014, Tupper_K: 344.274,
    source: 'NIST-IG',
  },
  'BENZENE|P-DIISOPROPYLBENZENE': {
    aij: -2.00939, aji: 3.17037,
    bij: -485.372, bji: -798.078,
    Tlower_K: 334.8, Tupper_K: 430.8,
    source: 'NIST-IG',
  },
  'ISOPROPYLBENZENE|PROPYLENE': {
    aij: -3.54578, aji: 8.50646,
    bij: 867.174, bji: -2030.88,
    Tlower_K: 273.15, Tupper_K: 353.123,
    source: 'NIST-IG',
  },
  'PROPANE|PROPYLENE': {
    aij: 1.61905, aji: -1.81887,
    bij: -383.022, bji: 419.087,
    Tlower_K: 230.008, Tupper_K: 273.15,
    source: 'NIST-IG',
  },
};

// ─── NIST Wilson Binary Params (NIST-IG) ──────────────────────────────────────
export const NIST_WILSON_BINARY_PARAMS: Record<string, WilsonBinaryParams> = {
  'BENZENE|ISOPROPYLBENZENE': {
    aij: 0.0754498, aji: 0.0823931,
    bij: -40.5283, bji: -8.48864,
    Tlower_K: 298.15, Tupper_K: 409.531,
  },
  'BENZENE|PROPYLENE': {
    aij: -3.66307, aji: 4.02386,
    bij: 973.333, bji: -1253.8,
    Tlower_K: 279.99, Tupper_K: 343.132,
  },
  'BENZENE|PROPANE': {
    aij: -0.0622453, aji: 0.912073,
    bij: -211.84, bji: -341.663,
    Tlower_K: 160.014, Tupper_K: 344.274,
  },
  'BENZENE|P-DIISOPROPYLBENZENE': {
    aij: -1.15597, aji: 1.30451,
    bij: 160.058, bji: -213.162,
    Tlower_K: 334.8, Tupper_K: 430.8,
  },
  'ISOPROPYLBENZENE|PROPYLENE': {
    aij: -8.94489, aji: 3.08418,
    bij: 2058.27, bji: -702.912,
    Tlower_K: 273.15, Tupper_K: 353.123,
  },
  'PROPANE|PROPYLENE': {
    aij: 1.58833, aji: -1.37373,
    bij: -437.509, bji: 357.612,
    Tlower_K: 230.008, Tupper_K: 273.15,
  },
};

// ─── NIST NRTL Binary Params (NIST-IG) ───────────────────────────────────────
export interface NistNrtlBinaryParams {
  aij: number; aji: number;
  bij: number; bji: number;
  cij: number;  // alpha_ij (NRTL non-randomness)
  Tlower_K: number; Tupper_K: number;
}

export const NIST_NRTL_BINARY_PARAMS: Record<string, NistNrtlBinaryParams> = {
  'BENZENE|ISOPROPYLBENZENE': {
    aij: -0.471595, aji: 0.410255,
    bij: 37.9866, bji: 3.88272,
    cij: 0.5,
    Tlower_K: 298.15, Tupper_K: 409.531,
  },
  'BENZENE|P-DIISOPROPYLBENZENE': {
    aij: -1.02131, aji: 0.661867,
    bij: 143.881, bji: -45.6104,
    cij: 0.5,
    Tlower_K: 334.8, Tupper_K: 430.8,
  },
  'BENZENE|PROPANE': {
    aij: -0.997877, aji: 0.392697,
    bij: 301.744, bji: 155.513,
    cij: 0.220653,
    Tlower_K: 160.014, Tupper_K: 344.274,
  },
  'BENZENE|PROPYLENE': {
    aij: 0, aji: 0,  // no NIST-IG for this pair, use RK
    bij: 0, bji: 0,
    cij: 0.3,
    Tlower_K: 279.99, Tupper_K: 343.132,
  },
  'ISOPROPYLBENZENE|PROPYLENE': {
    aij: -3.54578, aji: 8.50646,
    bij: 867.174, bji: -2030.88,
    cij: 0.5,
    Tlower_K: 273.15, Tupper_K: 353.123,
  },
  'PROPANE|PROPYLENE': {
    aij: 0.00125927, aji: -0.0981501,
    bij: 1.66539, bji: 46.9617,
    cij: 0.499999,
    Tlower_K: 230.008, Tupper_K: 273.15,
  },
};

// ─── NIST PR kij Binary Params (NIST-EOS) ────────────────────────────────────
// kij = kij1 + kij2·T + kij3·T² (T in K)
export interface NistPRKijParams {
  kij1: number; kij2: number; kij3: number;
  Tmin_K: number; Tmax_K: number;
}

export const NIST_PRKBV: Record<string, NistPRKijParams> = {
  'BENZENE|PROPYLENE': { kij1: 0.0105975, kij2: 0, kij3: 0, Tmin_K: 283.098, Tmax_K: 543.15 },
  'PROPANE|PROPYLENE': { kij1: 0.00852572, kij2: 0, kij3: 0, Tmin_K: 230.008, Tmax_K: 360.873 },
  'BENZENE|PROPANE': { kij1: 0.0328503, kij2: 0, kij3: 0, Tmin_K: 310.881, Tmax_K: 510.917 },
  'BENZENE|P-DIISOPROPYLBENZENE': { kij1: -0.0234875, kij2: 0, kij3: 0, Tmin_K: 334.8, Tmax_K: 430.8 },
  'BENZENE|ISOPROPYLBENZENE': { kij1: -0.0112698, kij2: 0, kij3: 0, Tmin_K: 335.19, Tmax_K: 409.531 },
  'ISOPROPYLBENZENE|PROPANE': { kij1: 0.0109642, kij2: 0, kij3: 0, Tmin_K: 376.85, Tmax_K: 377.35 },
};

// ─── CPA (Cubic-Plus-Association) kij Binary Params (AP-EOS) ─────────────────
export const CPA_KIJ: Record<string, EosKij> = {
  'BENZENE|PROPANE': { kij1: 0.030854, source: 'APV140 AP-EOS' },
  'PROPANE|PROPYLENE': { kij1: 0.008, source: 'APV140 AP-EOS' },
};

// ─── RKS T-dependent kij (RKSKBV): kij = kij1 + kij2·T + kij3/T ─────────────
export interface RksKijParams {
  kij1: number; kij2: number; kij3: number;
  Tmin_K: number; Tmax_K: number;
}
export const RKS_KBV: Record<string, RksKijParams> = {
  'PROPANE|PROPYLENE':  { kij1: 0.008, kij2: 0, kij3: 0, Tmin_K: 0, Tmax_K: 1000 },
  'BENZENE|PROPANE':    { kij1: 0.02,  kij2: 0, kij3: 0, Tmin_K: 0, Tmax_K: 1000 },
};

// ─── Lee-Kesler-Plocker EOS kij (LKPKIJ) ─────────────────────────────────────
export const LKP_KIJ: Record<string, { kij: number }> = {
  'PROPANE|PROPYLENE': { kij: -0.0081 },
  'PROPANE|BENZENE':   { kij: 0.0133 },
};

function createPresetFlowsheet(): {
  streams: Record<string, MaterialStream>;
  units: Record<string, UnitOperation>;
} {
  const n = 5; // propylene, propane, benzene, cumene, DIPB

  // ═══════════════════════════════════════════════════════════════
  // Cumene Production via Friedel-Crafts Alkylation
  // ═══════════════════════════════════════════════════════════════
  //
  //  C₆H₆ + C₃H₆ → C₉H₁₂  (main: benzene + propylene → cumene)
  //  C₉H₁₂ + C₃H₆ → C₁₂H₁₈ (side: cumene + propylene → p-DIPB)
  //
  //  Row 1 (y=150):
  //    S1 (Propylene Feed) → U1 (Mixer) → S3 → U2 (Feed Heater) → S4 → U3 (Reactor) → S5 → U4 (Cooler)
  //
  //  Row 2 (y=350):
  //    U4 → S6 → U5 (Flash) ─┬─ S7 (Vapor = lights to fuel) out
  //                           └─ S8 (Liquid) → U6 (Column 1: Benzene/Cumene)
  //
  //  Row 3 (y=550):
  //    U6 ─┬─ S9 (Distillate = unreacted benzene → recycle → U1)
  //        └─ S10 (Bottoms = cumene + DIPB) → U7 (Column 2: Cumene/DIPB)
  //            ├─ S11 (Distillate = Cumene Product)
  //            └─ S12 (Bottoms = DIPB byproduct)
  //
  //  S2 (Fresh Benzene Feed) → U1
  //  S9 (Benzene recycle) → U1
  //

  // ── S1: Fresh Propylene Feed (5% propane inert) ──
  const S1: MaterialStream = {
    ...createDefaultStream('S1', 'Propylene Feed', n),
    T_K: 313.15, P_Pa: 2500000, totalFlow_molps: 105,
    // propylene=0.95, propane=0.05, benzene=0, cumene=0, DIPB=0
    moleFractions: [0.95, 0.05, 0, 0, 0],
    targetUnitId: 'U1',
  };

  // ── S2: Fresh Benzene Feed ──
  const S2: MaterialStream = {
    ...createDefaultStream('S2', 'Benzene Feed', n),
    T_K: 313.15, P_Pa: 2500000, totalFlow_molps: 115,
    // pure benzene with slight excess (molar ratio ~1.1:1 vs propylene)
    moleFractions: [0, 0, 1, 0, 0],
    targetUnitId: 'U1',
  };

  // ── S3: Mixed Feed (Mixer outlet) ──
  const S3: MaterialStream = {
    ...createDefaultStream('S3', 'Mixed Feed', n),
    T_K: 313, P_Pa: 2500000, totalFlow_molps: 220,
    moleFractions: [0.454, 0.024, 0.522, 0, 0],
    sourceUnitId: 'U1', targetUnitId: 'U2',
  };

  // ── S4: Heated Feed (to reactor) ──
  const S4: MaterialStream = {
    ...createDefaultStream('S4', 'Heated Feed', n),
    T_K: 623.15, P_Pa: 2500000, totalFlow_molps: 220,
    moleFractions: [0.454, 0.024, 0.522, 0, 0],
    sourceUnitId: 'U2', targetUnitId: 'U3',
  };

  // ── S5: Reactor Effluent ──
  const S5: MaterialStream = {
    ...createDefaultStream('S5', 'Reactor Effluent', n),
    T_K: 623.15, P_Pa: 2500000, totalFlow_molps: 220,
    // After reaction: propylene mostly consumed, cumene formed, small DIPB
    moleFractions: [0.02, 0.024, 0.10, 0.83, 0.026],
    sourceUnitId: 'U3', targetUnitId: 'U4',
  };

  // ── S6: Cooled Effluent ──
  const S6: MaterialStream = {
    ...createDefaultStream('S6', 'Cooled Effluent', n),
    T_K: 323.15, P_Pa: 2500000, totalFlow_molps: 220,
    moleFractions: [0.02, 0.024, 0.10, 0.83, 0.026],
    sourceUnitId: 'U4', targetUnitId: 'U5',
  };

  // ── S7: Flash Vapor (light gases — fuel gas) ──
  const S7: MaterialStream = {
    ...createDefaultStream('S7', 'Fuel Gas', n),
    sourceUnitId: 'U5',
  };

  // ── S8: Flash Liquid (to Column 1) ──
  const S8: MaterialStream = {
    ...createDefaultStream('S8', 'Flash Liquid', n),
    sourceUnitId: 'U5', targetUnitId: 'U6',
  };

  // ── S9: Column 1 Distillate (benzene recycle → mixer) — TEAR STREAM ──
  const S9: MaterialStream = {
    ...createDefaultStream('S9', 'Benzene Recycle', n),
    sourceUnitId: 'U6', targetUnitId: 'U1',
  };

  // ── S10: Column 1 Bottoms (cumene + DIPB → Column 2) ──
  const S10: MaterialStream = {
    ...createDefaultStream('S10', 'Cumene + DIPB', n),
    sourceUnitId: 'U6', targetUnitId: 'U7',
  };

  // ── S11: Column 2 Distillate (Cumene product) ──
  const S11: MaterialStream = {
    ...createDefaultStream('S11', 'Cumene Product', n),
    sourceUnitId: 'U7',
  };

  // ── S12: Column 2 Bottoms (DIPB byproduct) ──
  const S12: MaterialStream = {
    ...createDefaultStream('S12', 'DIPB Byproduct', n),
    sourceUnitId: 'U7',
  };

  // ═══════════════════════════════════════════════════════════════
  // Unit Operations
  // ═══════════════════════════════════════════════════════════════

  const U1: UnitOperation = {
    id: 'U1', name: 'Feed Mixer', type: 'Mixer',
    ports: [
      { id: 'U1-in1', label: 'Propylene', type: 'inlet', position: 'left', streamId: 'S1' },
      { id: 'U1-in2', label: 'Benzene', type: 'inlet', position: 'left', streamId: 'S2' },
      { id: 'U1-in3', label: 'Recycle', type: 'inlet', position: 'left', streamId: 'S9' },
      { id: 'U1-out', label: 'Outlet', type: 'outlet', position: 'right', streamId: 'S3' },
    ],
    params: { outletP_Pa: 2500000 },
    solved: false, duty_W: 0, errors: [],
  };

  const U2: UnitOperation = {
    id: 'U2', name: 'Feed Heater', type: 'Heater',
    ports: [
      { id: 'U2-in', label: 'Inlet', type: 'inlet', position: 'left', streamId: 'S3' },
      { id: 'U2-out', label: 'Outlet', type: 'outlet', position: 'right', streamId: 'S4' },
    ],
    params: { targetT_K: 623.15, outletP_Pa: 2500000 },
    solved: false, duty_W: 0, errors: [],
  };

  const U3: UnitOperation = {
    id: 'U3', name: 'Alkylation Reactor', type: 'RStoic',
    ports: [
      { id: 'U3-in', label: 'Inlet', type: 'inlet', position: 'left', streamId: 'S4' },
      { id: 'U3-out', label: 'Outlet', type: 'outlet', position: 'right', streamId: 'S5' },
    ],
    params: {
      reactions: [
        {
          // Main: Benzene + Propylene → Cumene
          // coefficients: [propylene, propane, benzene, cumene, DIPB]
          coefficients: [-1, 0, -1, 1, 0],
          keyComponentIndex: 0,     // key = propylene
          conversion: 0.95,         // 95% propylene conversion
        },
        {
          // Side: Cumene + Propylene → DIPB
          coefficients: [-1, 0, 0, -1, 1],
          keyComponentIndex: 3,     // key = cumene
          conversion: 0.05,         // 5% of cumene goes to DIPB
        },
      ],
      outletT_K: 623.15,
      outletP_Pa: 2500000,
    },
    solved: false, duty_W: 0, errors: [],
  };

  const U4: UnitOperation = {
    id: 'U4', name: 'Product Cooler', type: 'Heater',
    ports: [
      { id: 'U4-in', label: 'Inlet', type: 'inlet', position: 'left', streamId: 'S5' },
      { id: 'U4-out', label: 'Outlet', type: 'outlet', position: 'right', streamId: 'S6' },
    ],
    params: { targetT_K: 323.15, outletP_Pa: 2500000 },
    solved: false, duty_W: 0, errors: [],
  };

  const U5: UnitOperation = {
    id: 'U5', name: 'Gas Separator', type: 'Flash',
    ports: [
      { id: 'U5-in', label: 'Inlet', type: 'inlet', position: 'left', streamId: 'S6' },
      { id: 'U5-vap', label: 'Vapor', type: 'outlet', position: 'top', streamId: 'S7' },
      { id: 'U5-liq', label: 'Liquid', type: 'outlet', position: 'bottom', streamId: 'S8' },
    ],
    params: { T_K: 323.15, P_Pa: 101325 },
    solved: false, duty_W: 0, errors: [],
  };

  const U6: UnitOperation = {
    id: 'U6', name: 'Benzene Column', type: 'Column',
    ports: [
      { id: 'U6-feed', label: 'Feed', type: 'inlet', position: 'left', streamId: 'S8' },
      { id: 'U6-dist', label: 'Distillate', type: 'outlet', position: 'top', streamId: 'S9' },
      { id: 'U6-bot', label: 'Bottoms', type: 'outlet', position: 'bottom', streamId: 'S10' },
    ],
    params: {
      rigorousMode: false,
      lightKeyIndex: 2, heavyKeyIndex: 3,   // benzene (light) / cumene (heavy)
      lightKeyRecovery: 0.99, heavyKeyRecovery: 0.99,
      refluxRatioMultiplier: 1.3,
      nStages: 30, feedStage: 15, refluxRatio: 2.0, distillateRate_molps: 0,
      condenserP_Pa: 101325, reboilerP_Pa: 101325, condenserType: 'Total',
    },
    solved: false, duty_W: 0, errors: [],
  };

  const U7: UnitOperation = {
    id: 'U7', name: 'Product Column', type: 'Column',
    ports: [
      { id: 'U7-feed', label: 'Feed', type: 'inlet', position: 'left', streamId: 'S10' },
      { id: 'U7-dist', label: 'Distillate', type: 'outlet', position: 'top', streamId: 'S11' },
      { id: 'U7-bot', label: 'Bottoms', type: 'outlet', position: 'bottom', streamId: 'S12' },
    ],
    params: {
      rigorousMode: false,
      lightKeyIndex: 3, heavyKeyIndex: 4,   // cumene (light) / DIPB (heavy)
      lightKeyRecovery: 0.995, heavyKeyRecovery: 0.995,
      refluxRatioMultiplier: 1.3,
      nStages: 20, feedStage: 10, refluxRatio: 1.5, distillateRate_molps: 0,
      condenserP_Pa: 101325, reboilerP_Pa: 101325, condenserType: 'Total',
    },
    solved: false, duty_W: 0, errors: [],
  };

  return {
    streams: { S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12 },
    units: { U1, U2, U3, U4, U5, U6, U7 },
  };
}

// ────────────────────────────────────────────────────────────────
// Display names for unit operation types
// ────────────────────────────────────────────────────────────────

const TYPE_DISPLAY_NAMES: Record<UnitOpType, string> = {
  Heater: 'Heater',
  Flash: 'Flash Drum',
  Mixer: 'Mixer',
  Valve: 'Valve',
  Splitter: 'Splitter',
  Pump: 'Pump',
  Compressor: 'Compressor',
  HeatX: 'Heat Exchanger',
  Pipe: 'Pipe',
  RStoic: 'RStoic',
  Column: 'Column',
  RadFrac: 'RadFrac',
  RateFrac: 'RateFrac',
  ThreePhaseFlash: 'Flash3',
  CSTR: 'CSTR',
  PFR: 'PFR',
  Absorber: 'Absorber',
  Extractor: 'Extractor',
  ComponentSeparator: 'Sep',
  RYield: 'RYield',
  Decanter: 'Decanter',
  REquil: 'REquil',
  RGibbs: 'RGibbs',
  MHeatX: 'MHeatX',
  RBatch: 'RBatch',
  Crystallizer: 'Crystallizer',
  Crusher: 'Crusher',
  Dryer: 'Dryer',
  Membrane: 'Membrane',
  Cyclone: 'Cyclone',
  Filter: 'Filter',
  Screen: 'Screen',
  Centrifuge: 'Centrifuge',
  SteamDrum: 'Steam Drum',
  SteamHeater: 'Steam Heater',
  SteamTurbine: 'Steam Turbine',
  SteamValve: 'Steam Valve',
  SteamHeader: 'Steam Header',
  SteamTrap: 'Steam Trap',
  PIDController: 'PID Controller',
  LeadLag: 'Lead-Lag',
  DeadTime: 'Dead Time',
  SignalSelector: 'Signal Selector',
  Analyzer: 'Analyzer',
  ChargeBalance: 'Charge Balance',
  ElectrolyteEquilibrium: 'Electrolyte Equilibrium',
};

// ────────────────────────────────────────────────────────────────
// Store
// ────────────────────────────────────────────────────────────────

type UndoSnapshot = {
  streams: Record<string, MaterialStream>;
  units: Record<string, UnitOperation>;
  compounds: CanopyCompound[];
  nodePositions: Record<string, { x: number; y: number }>;
};

const MAX_UNDO = 50;
const _undoStack: UndoSnapshot[] = [];
const _redoStack: UndoSnapshot[] = [];

/** Push current state to undo stack before mutating. */
function _pushUndo(s: CanopyState) {
  _undoStack.push({
    streams: s.streams,
    units: s.units,
    compounds: s.compounds,
    nodePositions: s.nodePositions,
  });
  if (_undoStack.length > MAX_UNDO) _undoStack.shift();
  _redoStack.length = 0; // clear redo on new action
}

export const useCanopyStore = create<CanopyState>((set, get) => ({
  // Navigation
  currentScreen: 'setup',
  setScreen: (screen) => set({ currentScreen: screen }),

  // Setup
  compounds: [],
  fluidPackage: 'Ideal',
  unitSystem: { ...DEFAULT_UNIT_SYSTEM },
  interactionParams: {},
  setCompounds: (compounds) => set({ compounds }),
  addCompound: (compound) => set(s => ({
    compounds: s.compounds.some(c => c.name === compound.name)
      ? s.compounds
      : [...s.compounds, compound],
  })),
  removeCompound: (name) => set(s => ({
    compounds: s.compounds.filter(c => c.name !== name),
  })),
  setFluidPackage: (pkg) => set({ fluidPackage: pkg }),
  setUnitSystem: (patch) => set(s => ({ unitSystem: { ...s.unitSystem, ...patch } })),
  setInteractionParams: (params) => set(s => ({
    interactionParams: { ...s.interactionParams, ...params },
  })),

  // Flowsheet
  streams: {},
  units: {},
  addStream: (stream) => { _pushUndo(get()); set(s => ({ streams: { ...s.streams, [stream.id]: stream }, canUndo: true, canRedo: false })); },
  updateStream: (id, patch) => { _pushUndo(get()); set(s => ({
    streams: { ...s.streams, [id]: { ...s.streams[id], ...patch } }, canUndo: true, canRedo: false,
  })); },
  removeStream: (id) => { _pushUndo(get()); set(s => {
    const { [id]: _, ...rest } = s.streams;
    return { streams: rest, canUndo: true, canRedo: false };
  }); },
  addUnit: (unit) => { _pushUndo(get()); set(s => ({ units: { ...s.units, [unit.id]: unit }, canUndo: true, canRedo: false })); },
  updateUnit: (id, patch) => { _pushUndo(get()); set(s => ({
    units: { ...s.units, [id]: { ...s.units[id], ...patch } }, canUndo: true, canRedo: false,
  })); },
  removeUnit: (id) => { _pushUndo(get()); set(s => {
    const { [id]: _, ...rest } = s.units;
    return { units: rest, canUndo: true, canRedo: false };
  }); },
  renameUnit: (id, newName) => set(s => ({
    units: { ...s.units, [id]: { ...s.units[id], name: newName } },
  })),
  renameStream: (id, newName) => set(s => ({
    streams: { ...s.streams, [id]: { ...s.streams[id], name: newName } },
  })),

  disconnectStream: (streamId) => set(s => {
    const stream = s.streams[streamId];
    if (!stream) return s;
    const updatedUnits = { ...s.units };
    // Clear port references on source unit
    if (stream.sourceUnitId && updatedUnits[stream.sourceUnitId]) {
      updatedUnits[stream.sourceUnitId] = {
        ...updatedUnits[stream.sourceUnitId],
        ports: updatedUnits[stream.sourceUnitId].ports.map(p =>
          p.streamId === streamId ? { ...p, streamId: null } : p
        ),
      };
    }
    // Clear port references on target unit
    if (stream.targetUnitId && updatedUnits[stream.targetUnitId]) {
      updatedUnits[stream.targetUnitId] = {
        ...updatedUnits[stream.targetUnitId],
        ports: updatedUnits[stream.targetUnitId].ports.map(p =>
          p.streamId === streamId ? { ...p, streamId: null } : p
        ),
      };
    }
    // Remove the stream entirely
    const { [streamId]: _, ...restStreams } = s.streams;
    return { streams: restStreams, units: updatedUnits, solved: false };
  }),

  connectPorts: (sourceUnitId, sourceHandleId, targetUnitId, targetHandleId) => {
    const state = get();
    const srcUnit = state.units[sourceUnitId];
    const tgtUnit = state.units[targetUnitId];
    if (!srcUnit || !tgtUnit) return;

    // Find the outlet port on the source unit
    const srcPort = srcUnit.ports.find(p => {
      const handleId = p.id.slice(sourceUnitId.length + 1) || '';
      return p.type === 'outlet' && (handleId === sourceHandleId || p.id === sourceHandleId);
    });
    // Find the inlet port on the target unit
    const tgtPort = tgtUnit.ports.find(p => {
      const handleId = p.id.slice(targetUnitId.length + 1) || '';
      return p.type === 'inlet' && (handleId === targetHandleId || p.id === targetHandleId);
    });
    if (!srcPort || !tgtPort) return;

    // If source port already has a stream, disconnect it
    if (srcPort.streamId) {
      state.disconnectStream(srcPort.streamId);
    }
    // If target port already has a stream, disconnect it
    if (tgtPort.streamId) {
      state.disconnectStream(tgtPort.streamId);
    }

    // Re-read state after disconnections
    const freshState = get();

    // Create a new stream connecting them
    const streamId = freshState.nextStreamId();
    const n = freshState.compounds.length;
    const newStream: MaterialStream = {
      ...createDefaultStream(streamId, streamId, n),
      T_K: 298.15,
      P_Pa: 101325,
      totalFlow_molps: 100,
      moleFractions: n > 0 ? new Array(n).fill(1 / n) : [],
      sourceUnitId,
      sourcePortId: srcPort.id,
      targetUnitId,
      targetPortId: tgtPort.id,
    };

    // Update ports on both units
    const updatedUnits = { ...freshState.units };
    updatedUnits[sourceUnitId] = {
      ...updatedUnits[sourceUnitId],
      ports: updatedUnits[sourceUnitId].ports.map(p =>
        p.id === srcPort.id ? { ...p, streamId } : p
      ),
    };
    updatedUnits[targetUnitId] = {
      ...updatedUnits[targetUnitId],
      ports: updatedUnits[targetUnitId].ports.map(p =>
        p.id === tgtPort.id ? { ...p, streamId } : p
      ),
    };

    // Compute an explicit label position so the label never auto-follows units
    const srcPos = freshState.nodePositions[sourceUnitId] ?? { x: 300, y: 200 };
    const tgtPos = freshState.nodePositions[targetUnitId] ?? { x: 600, y: 200 };
    const labelPos = {
      x: (srcPos.x + tgtPos.x) / 2,
      y: (srcPos.y + tgtPos.y) / 2 - 30,
    };

    set({
      streams: { ...freshState.streams, [streamId]: newStream },
      units: updatedUnits,
      nodePositions: { ...freshState.nodePositions, [`label-${streamId}`]: labelPos },
      solved: false,
    });
  },

  // Add unit from palette
  addUnitFromPalette: (type, position) => {
    const state = get();
    const unitId = state.nextUnitId();
    const portDefs: Record<UnitOpType, () => UnitPort[]> = {
      Heater: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Flash: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-vap`, label: 'Vapor', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-liq`, label: 'Liquid', type: 'outlet', position: 'bottom', streamId: null },
      ],
      Mixer: () => [
        { id: `${unitId}-in1`, label: 'Inlet 1', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-in2`, label: 'Inlet 2', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Valve: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Splitter: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out1`, label: 'Outlet 1', type: 'outlet', position: 'right', streamId: null },
        { id: `${unitId}-out2`, label: 'Outlet 2', type: 'outlet', position: 'right', streamId: null },
      ],
      Pump: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Compressor: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      HeatX: () => [
        { id: `${unitId}-hot-in`, label: 'Hot Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-cold-in`, label: 'Cold Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-hot-out`, label: 'Hot Outlet', type: 'outlet', position: 'right', streamId: null },
        { id: `${unitId}-cold-out`, label: 'Cold Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Pipe: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      RStoic: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Column: () => [
        { id: `${unitId}-feed`, label: 'Feed', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-dist`, label: 'Distillate', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-bot`, label: 'Bottoms', type: 'outlet', position: 'bottom', streamId: null },
      ],
      RadFrac: () => [
        { id: `${unitId}-feed`, label: 'Feed', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-dist`, label: 'Distillate', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-side-vap`, label: 'Side Vapor', type: 'outlet', position: 'right', streamId: null },
        { id: `${unitId}-side-liq`, label: 'Side Liquid', type: 'outlet', position: 'right', streamId: null },
        { id: `${unitId}-bot`, label: 'Bottoms', type: 'outlet', position: 'bottom', streamId: null },
      ],
      RateFrac: () => [
        { id: `${unitId}-feed`, label: 'Feed', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-dist`, label: 'Distillate', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-side-vap`, label: 'Side Vapor', type: 'outlet', position: 'right', streamId: null },
        { id: `${unitId}-side-liq`, label: 'Side Liquid', type: 'outlet', position: 'right', streamId: null },
        { id: `${unitId}-bot`, label: 'Bottoms', type: 'outlet', position: 'bottom', streamId: null },
      ],
      ThreePhaseFlash: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-vap`, label: 'Vapor', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-liq1`, label: 'Liquid I', type: 'outlet', position: 'right', streamId: null },
        { id: `${unitId}-liq2`, label: 'Liquid II', type: 'outlet', position: 'bottom', streamId: null },
      ],
      CSTR: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      PFR: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Absorber: () => [
        { id: `${unitId}-gas-in`, label: 'Gas Inlet', type: 'inlet', position: 'bottom', streamId: null },
        { id: `${unitId}-liq-in`, label: 'Liquid Inlet', type: 'inlet', position: 'top', streamId: null },
        { id: `${unitId}-gas-out`, label: 'Gas Outlet', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-liq-out`, label: 'Liquid Outlet', type: 'outlet', position: 'bottom', streamId: null },
      ],
      Extractor: () => [
        { id: `${unitId}-feed-in`, label: 'Feed Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-solvent-in`, label: 'Solvent Inlet', type: 'inlet', position: 'top', streamId: null },
        { id: `${unitId}-extract-out`, label: 'Extract', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-raffinate-out`, label: 'Raffinate', type: 'outlet', position: 'right', streamId: null },
      ],
      ComponentSeparator: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-top`, label: 'Top', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-bot`, label: 'Bottom', type: 'outlet', position: 'bottom', streamId: null },
      ],
      RYield: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Decanter: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-liq1`, label: 'Liquid I', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-liq2`, label: 'Liquid II', type: 'outlet', position: 'bottom', streamId: null },
      ],
      REquil: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      RGibbs: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      MHeatX: () => [
        { id: `${unitId}-hot-in`, label: 'Hot Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-cold-in`, label: 'Cold Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-hot-out`, label: 'Hot Outlet', type: 'outlet', position: 'right', streamId: null },
        { id: `${unitId}-cold-out`, label: 'Cold Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      RBatch: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Crystallizer: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Crusher: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Dryer: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Membrane: () => [
        { id: `${unitId}-in`, label: 'Feed', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-perm`, label: 'Permeate', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-ret`, label: 'Retentate', type: 'outlet', position: 'right', streamId: null },
      ],
      Cyclone: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-ovhd`, label: 'Overflow', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-uflow`, label: 'Underflow', type: 'outlet', position: 'bottom', streamId: null },
      ],
      Filter: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-filtrate`, label: 'Filtrate', type: 'outlet', position: 'right', streamId: null },
        { id: `${unitId}-cake`, label: 'Cake', type: 'outlet', position: 'bottom', streamId: null },
      ],
      Screen: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-oversize`, label: 'Oversize', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-undersize`, label: 'Undersize', type: 'outlet', position: 'right', streamId: null },
      ],
      Centrifuge: () => [
        { id: `${unitId}-in`, label: 'Feed', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-centrate`, label: 'Centrate', type: 'outlet', position: 'right', streamId: null },
        { id: `${unitId}-cake`, label: 'Cake', type: 'outlet', position: 'bottom', streamId: null },
      ],
      SteamDrum: () => [
        { id: `${unitId}-in`, label: 'Feed', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-vap`, label: 'Vapor', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-liq`, label: 'Liquid', type: 'outlet', position: 'bottom', streamId: null },
      ],
      SteamHeater: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      SteamTurbine: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      SteamValve: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      SteamHeader: () => [
        { id: `${unitId}-in1`, label: 'Inlet 1', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-in2`, label: 'Inlet 2', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      SteamTrap: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-vap`, label: 'Flash Vapor', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-liq`, label: 'Condensate', type: 'outlet', position: 'right', streamId: null },
      ],
      PIDController: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      LeadLag: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      DeadTime: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      SignalSelector: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Analyzer: () => [],
      ChargeBalance: () => [
        { id: `${unitId}-in`, label: 'Feed', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Balanced Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      ElectrolyteEquilibrium: () => [
        { id: `${unitId}-in`, label: 'Feed', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Equilibrated Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
    };
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const defaultParams: Record<UnitOpType, Record<string, any>> = {
      Heater: { targetT_K: 373.15, outletP_Pa: 101325, waterSaturatedFeed: false, waterAddedForSaturation_molps: 0, excessFreeWater_molps: 0 },
      Flash: { T_K: 373.15, P_Pa: 101325, waterSaturatedFeed: false, waterAddedForSaturation_molps: 0, excessFreeWater_molps: 0 },
      Mixer: { outletP_Pa: 101325 },
      Valve: { outletP_Pa: 50000, waterSaturatedFeed: false, waterAddedForSaturation_molps: 0, excessFreeWater_molps: 0 },
      Splitter: { splitRatio: 0.5, outletP_Pa: 101325 },
      Pump: { outletP_Pa: 500000, efficiency: 0.75, driverEfficiency: 0.95 },
      Compressor: { outletP_Pa: 500000, efficiency: 0.72, model: 'isentropic' },
      HeatX: { spec: 'hotOutletT', specValue: 323.15, flowArrangement: 'counter', hotOutletP_Pa: 0, coldOutletP_Pa: 0 },
      Pipe: { length_m: 100, diameter_m: 0.1, roughness_m: 0.000046, elevation_m: 0, K_fittings: 0 },
      RStoic: { reactions: [], outletT_K: 373.15, outletP_Pa: 101325 },
      Column: {
        rigorousMode: false,
        // Shortcut (FUG) parameters
        lightKeyIndex: 0, heavyKeyIndex: 1,
        lightKeyRecovery: 0.95, heavyKeyRecovery: 0.95,
        refluxRatioMultiplier: 1.3,
        // Rigorous parameters
        nStages: 20, feedStage: 10, refluxRatio: 2.0, distillateRate_molps: 0, murphreeEfficiency: 1.0,
        spec1Mode: 'None', spec1Value: 0, spec1ComponentIndex: 0,
        spec2Mode: 'None', spec2Value: 0, spec2ComponentIndex: 0,
        // Common
        condenserP_Pa: 101325, reboilerP_Pa: 101325, condenserType: 'Total',
      },
      RadFrac: {
        rigorousMode: true,
        nStages: 20, feedStage: 10, refluxRatio: 2.0, distillateRate_molps: 0,
        murphreeEfficiency: 1.0,
        spec1Mode: 'None', spec1Value: 0, spec1ComponentIndex: 0,
        spec2Mode: 'None', spec2Value: 0, spec2ComponentIndex: 0,
        liquidSideDrawStage: undefined, liquidSideDrawFraction: 0,
        vaporSideDrawStage: undefined, vaporSideDrawFraction: 0,
        condenserP_Pa: 101325, reboilerP_Pa: 101325, condenserType: 'Total',
      },
      RateFrac: {
        rigorousMode: true,
        nStages: 20, feedStage: 10, refluxRatio: 2.0, distillateRate_molps: 0,
        spec1Mode: 'None', spec1Value: 0, spec1ComponentIndex: 0,
        spec2Mode: 'None', spec2Value: 0, spec2ComponentIndex: 0,
        liquidSideDrawStage: undefined, liquidSideDrawFraction: 0,
        vaporSideDrawStage: undefined, vaporSideDrawFraction: 0,
        condenserP_Pa: 101325, reboilerP_Pa: 150000, condenserType: 'Total',
        columnDiameter_m: 1.5, traySpacing_m: 0.6, weirHeight_m: 0.05,
        interfacialArea_m2pm3: 150, kL_mps: 2e-4, kV_mps: 5e-3,
        liquidHoldup_s: 5, vaporHoldup_s: 0.5,
        internalsMode: 'Tray', pressureRelaxation: 0.6,
      },
      ThreePhaseFlash: { T_K: 373.15, P_Pa: 101325 },
      CSTR: { reactions: [], outletT_K: 373.15, outletP_Pa: 101325, volume_m3: 1.0 },
      PFR: { reactions: [], outletT_K: 373.15, outletP_Pa: 101325, volume_m3: 1.0, nSteps: 100 },
      Absorber: { nStages: 10, P_Pa: 101325 },
      Extractor: { nStages: 6, P_Pa: 101325, solventKeyIndex: 0 },
      ComponentSeparator: { splitFractions: [] },
      RYield: { yields: [], outletT_K: 373.15, outletP_Pa: 101325 },
      Decanter: {
        T_K: 298.15,
        P_Pa: 101325,
        dirtyWaterMode: false,
        freeWaterSplitFraction: 0.985,
        waterCarryoverFraction: 0.002,
        hydrocarbonToAqueousFraction: 0.002,
      },
      REquil: { reactions: [], outletT_K: 373.15, outletP_Pa: 101325 },
      RGibbs: { outletT_K: 373.15, outletP_Pa: 101325, elementMatrix: [] },
      MHeatX: { hotStreamIndices: [0], coldStreamIndices: [1], deltaT_min_K: 10 },
      RBatch: { reactions: [], batchTime_s: 3600, T_K: 373.15, P_Pa: 101325, volume_m3: 1.0 },
      Crystallizer: { residenceTime_s: 3600, soluteIndex: 0, C_sat_molpm3: 100, k_g: 1e-7, g_exponent: 1, k_b: 1e8, b_exponent: 1, j_exponent: 0.5 },
      Crusher: { workIndex_kWhpton: 15, feedSize_um: 25000, productSize_um: 5000, solidFlow_kgps: 1 },
      Dryer: { X_in_kgpkg: 0.5, X_out_kgpkg: 0.05, X_c_kgpkg: 0.2, X_eq_kgpkg: 0.02, dryFlow_kgps: 1, gasT_in_K: 423, gasFlow_kgps: 5, h_Wpm2K: 35, area_m2: 1, wetBulbT_K: 313.15, gasHumidityIn_kgpkg: 0 },
      Membrane: {
        model: 'stage-cut',
        stageCut: 0.25,
        selectivities: [],
        rejectionCoefficients: [],
        membraneArea_m2: 50,
        permeanceGPU: [],
        foulingFactor: 1,
        permeateP_Pa: 50662.5,
        retentateP_Pa: 101325,
      },
      Cyclone: { cutSize_um: 10, particleSizes_um: [], pressureDrop_Pa: 1500 },
      Filter: { retainedFractions: [], pressureDrop_Pa: 5000 },
      Screen: { cutSize_um: 5, sharpness: 4, particleSizes_um: [] },
      Centrifuge: { cutSize_um: 3, sharpness: 7, particleSizes_um: [], liquidSplitToCake: 0.08, pressureDrop_Pa: 5000 },
      SteamDrum: { outletP_Pa: 101325 },
      SteamHeater: { targetT_K: 423.15, outletP_Pa: 101325, steamLevel: 'Auto' },
      SteamTurbine: { outletP_Pa: 300000, efficiency: 0.75, generatorEfficiency: 0.96 },
      SteamValve: { outletP_Pa: 250000 },
      SteamHeader: { outletP_Pa: 250000 },
      SteamTrap: { outletP_Pa: 200000, maxFlashFraction: 0.2 },
      PIDController: {
        setpoint: 350,
        measuredProperty: 'T_K',
        measuredObjectType: 'stream',
        measuredObjectId: '',
        manipulatedUnitId: '',
        manipulatedParamName: 'targetT_K',
        gain: 1,
        integralTime_s: 5,
        derivativeTime_s: 0,
        controllerTimestep_s: 1,
        measurementLag_s: 0,
        outputRateLimit_per_s: undefined,
        bias: undefined,
        lowLimit: 0,
        highLimit: 1000,
        reverseActing: false,
      },
      LeadLag: {
        measuredObjectType: 'stream',
        measuredObjectId: '',
        measuredProperty: 'T_K',
        leadTime_s: 1,
        lagTime_s: 5,
        timestep_s: 1,
      },
      DeadTime: {
        measuredObjectType: 'stream',
        measuredObjectId: '',
        measuredProperty: 'T_K',
        deadTime_s: 5,
        timestep_s: 1,
      },
      SignalSelector: {
        mode: 'High',
        measuredObjectTypeA: 'stream',
        measuredObjectIdA: '',
        measuredPropertyA: 'T_K',
        measuredObjectTypeB: 'stream',
        measuredObjectIdB: '',
        measuredPropertyB: 'T_K',
      },
      Analyzer: {
        measuredObjectType: 'stream',
        measuredObjectId: '',
        propertyCode: 'TEMP',
        componentIndex: 0,
        signalValue: 0,
        propertyLabel: 'Temperature (K)',
      },
      ChargeBalance: {
        adjustComponentIndex: 0,
        targetCharge: 0,
        residualCharge: 0,
        adjustedComponentFlow_molps: 0,
      },
      ElectrolyteEquilibrium: {
        apparentSpeciesIndex: 0,
        cationIndex: 0,
        anionIndex: 0,
        cationStoich: 1,
        anionStoich: 1,
        log10K: 2,
        maxDissociationFraction: 0.999,
        reactionSetJson: '[]',
        dissociationFraction: 0,
        ionicStrength: 0,
        equilibriumConstant: 100,
        activeReactionCount: 1,
        equilibriumPasses: 1,
      },
    };
    const newUnit: UnitOperation = {
      id: unitId,
      name: `${TYPE_DISPLAY_NAMES[type] ?? type} ${Object.values(state.units).filter(u => u.type === type).length + 1}`,
      type,
      ports: portDefs[type](),
      params: defaultParams[type],
      solved: false,
      duty_W: 0,
      errors: [],
    };
    set(s => ({
      units: { ...s.units, [unitId]: newUnit },
      nodePositions: { ...s.nodePositions, [unitId]: position },
      solved: false,
    }));
  },

  // Add feed stream
  addFeedStream: (position) => {
    const state = get();
    const streamId = state.nextStreamId();
    const n = state.compounds.length;
    const newStream: MaterialStream = {
      ...createDefaultStream(streamId, streamId, n),
      T_K: 298.15,
      P_Pa: 101325,
      totalFlow_molps: 100,
      moleFractions: n > 0 ? new Array(n).fill(1 / n) : [],
    };
    set(s => ({
      streams: { ...s.streams, [streamId]: newStream },
      nodePositions: { ...s.nodePositions, [`label-${streamId}`]: position ?? { x: 50, y: 200 } },
      solved: false,
    }));
  },

  nextStreamId: () => {
    const ids = Object.keys(get().streams).map(id => parseInt(id.slice(1))).filter(n => !isNaN(n));
    return `S${Math.max(0, ...ids) + 1}`;
  },
  nextUnitId: () => {
    const ids = Object.keys(get().units).map(id => parseInt(id.slice(1))).filter(n => !isNaN(n));
    return `U${Math.max(0, ...ids) + 1}`;
  },

  // Tabs (supports both units and streams)
  openTabs: [],
  activeTab: null,
  openTab: (item) => set(s => ({
    openTabs: s.openTabs.some(t => t.type === item.type && t.id === item.id)
      ? s.openTabs
      : [...s.openTabs, item],
    activeTab: item,
  })),
  closeTab: (item) => set(s => {
    const newTabs = s.openTabs.filter(t => !(t.type === item.type && t.id === item.id));
    const isActive = s.activeTab?.type === item.type && s.activeTab?.id === item.id;
    return {
      openTabs: newTabs,
      activeTab: isActive ? (newTabs[newTabs.length - 1] ?? null) : s.activeTab,
    };
  }),
  setActiveTab: (item) => set({ activeTab: item }),

  // Solver
  solved: false,
  solverErrors: [],
  solverPaused: false,
  setSolverPaused: (paused) => set({ solverPaused: paused }),

  // Warnings
  dbWarnings: [],
  addDbWarning: (msg) => set((s) => ({
    dbWarnings: s.dbWarnings.includes(msg) ? s.dbWarnings : [...s.dbWarnings, msg],
  })),
  clearDbWarnings: () => set({ dbWarnings: [] }),

  // Design Specs
  designSpecs: [],
  addDesignSpec: (spec) => set((s) => ({ designSpecs: [...s.designSpecs, spec] })),
  updateDesignSpec: (id, patch) => set((s) => ({
    designSpecs: s.designSpecs.map(ds => ds.id === id ? { ...ds, ...patch } : ds),
  })),
  removeDesignSpec: (id) => set((s) => ({
    designSpecs: s.designSpecs.filter(ds => ds.id !== id),
  })),

  // Sensitivity Analysis
  sensitivities: [],
  addSensitivity: (sa) => set((s) => ({ sensitivities: [...s.sensitivities, sa] })),
  removeSensitivity: (id) => set((s) => ({
    sensitivities: s.sensitivities.filter(sa => sa.id !== id),
  })),
  runSensitivity: (id) => {
    const state = get();
    const sa = state.sensitivities.find(s => s.id === id);
    if (!sa || !sa.active) return;

    const results: number[][] = [];
    const vals1 = sa.variable.values;
    const vals2 = sa.variable2?.values ?? [0];
    const is2D = !!sa.variable2;

    for (const v1 of vals1) {
      for (const v2 of vals2) {
        // Set the variable value
        const unit1 = state.units[sa.variable.unitId];
        if (unit1) {
          state.units[sa.variable.unitId] = {
            ...unit1,
            params: { ...unit1.params, [sa.variable.paramName]: v1 },
          };
        }
        if (is2D && sa.variable2) {
          const unit2 = state.units[sa.variable2.unitId];
          if (unit2) {
            state.units[sa.variable2.unitId] = {
              ...unit2,
              params: { ...unit2.params, [sa.variable2.paramName]: v2 },
            };
          }
        }

        // Run solver
        state.solveAll();

        // Collect outputs
        const row: number[] = is2D ? [v1, v2] : [v1];
        for (const out of sa.outputs) {
          let val = NaN;
          if (out.objectType === 'stream') {
            const s = state.streams[out.objectId];
            if (s) val = (s as any)[out.property] ?? NaN;
          } else {
            const u = state.units[out.objectId];
            if (u) val = (u as any)[out.property] ?? (u.params as any)?.[out.property] ?? NaN;
          }
          row.push(val);
        }
        results.push(row);
      }
    }

    set((s) => ({
      sensitivities: s.sensitivities.map(si =>
        si.id === id ? { ...si, results } : si
      ),
    }));
  },

  // Calculator Blocks
  calculatorBlocks: [],
  addCalculatorBlock: (cb) => set((s) => ({ calculatorBlocks: [...s.calculatorBlocks, cb] })),
  removeCalculatorBlock: (id) => set((s) => ({
    calculatorBlocks: s.calculatorBlocks.filter(cb => cb.id !== id),
  })),
  updateCalculatorBlock: (id, patch) => set((s) => ({
    calculatorBlocks: s.calculatorBlocks.map(cb => cb.id === id ? { ...cb, ...patch } : cb),
  })),

  solveAll: () => {
    if (get().solverPaused) return; // Solver frozen by user
    const { compounds, streams, units, designSpecs, fluidPackage, interactionParams } = get();
    if (compounds.length === 0) {
      set({ solverErrors: ['No compounds selected.'], solved: false });
      return;
    }

    // ── Validate flowsheet before solving ──
    const warnings: string[] = [];
    const unitList = Object.values(units);
    const streamList = Object.values(streams);

    // Check for disconnected ports
    for (const u of unitList) {
      const unconnected = u.ports.filter(p => !p.streamId);
      if (unconnected.length > 0) {
        warnings.push(`${u.name}: ${unconnected.length} unconnected port${unconnected.length > 1 ? 's' : ''} (${unconnected.map(p => p.label).join(', ')})`);
      }
    }

    // Check feed streams have valid data
    for (const s of streamList) {
      if (!s.sourceUnitId) {
        if (!s.totalFlow_molps || s.totalFlow_molps <= 0) {
          warnings.push(`Feed stream "${s.name}": zero or missing flow rate.`);
        }
        if (!s.moleFractions || s.moleFractions.length !== compounds.length) {
          warnings.push(`Feed stream "${s.name}": composition does not match compound count.`);
        }
      }
    }

    // Check RStoic reactions match compound count
    for (const u of unitList) {
      if (u.type === 'RStoic' && u.params.reactions) {
        for (const rxn of u.params.reactions as any[]) {
          if (rxn.coefficients && rxn.coefficients.length !== compounds.length) {
            warnings.push(`${u.name}: reaction has ${rxn.coefficients.length} coefficients but ${compounds.length} compounds.`);
          }
        }
      }
    }

    if (warnings.length > 0) {
      set((s) => ({
        dbWarnings: [...s.dbWarnings.filter(w => !w.startsWith('Validation:')),
          ...warnings.map(w => `Validation: ${w}`)],
      }));
    }

    // ── Execute Calculator Blocks (before) ──
    const calcBlocks = get().calculatorBlocks;
    const execCalcBlock = (cb: CalculatorBlock, currentStreams: Record<string, MaterialStream>, currentUnits: Record<string, UnitOperation>) => {
      if (!cb.active) return;
      try {
        // Build context from imports
        const ctx: Record<string, number> = {};
        for (const imp of cb.imports) {
          if (imp.objectType === 'stream') {
            const s = currentStreams[imp.objectId];
            if (s) ctx[imp.varName] = (s as any)[imp.property] ?? 0;
          } else {
            const u = currentUnits[imp.objectId];
            if (u) ctx[imp.varName] = (u as any)[imp.property] ?? (u.params as any)?.[imp.property] ?? 0;
          }
        }
        // Evaluate expression (safe subset: only math operations + imported vars)
        const varNames = Object.keys(ctx);
        const varValues = Object.values(ctx);
        // eslint-disable-next-line no-new-func
        const fn = new Function(...varNames, `"use strict"; return (${cb.expression});`);
        const result = fn(...varValues);
        // Write exports
        for (const exp of cb.exports) {
          if (currentUnits[exp.unitId]) {
            currentUnits[exp.unitId] = {
              ...currentUnits[exp.unitId],
              params: { ...currentUnits[exp.unitId].params, [exp.paramName]: typeof result === 'number' ? result : 0 },
            };
          }
        }
      } catch {
        warnings.push(`Calculator "${cb.name}": evaluation error.`);
      }
    };
    for (const cb of calcBlocks.filter(c => c.executionOrder === 'before')) {
      execCalcBlock(cb, streams, units);
    }

    // ── Build ThermoOptions from store state ──
    const thermo: ThermoOptions = {
      fluidPackage: fluidPackage as ThermoOptions['fluidPackage'],
      interactionParams: Object.keys(interactionParams).length > 0 ? interactionParams : undefined,
    };

    const updatedStreams = { ...streams };
    const updatedUnits = { ...units };
    const errors: string[] = [];

    // ── Helper: apply outlet to stream preserving identity fields ──
    const applyOutlet = (
      streamId: string,
      outlet: Partial<MaterialStream>,
      unitId: string,
      portId: string,
    ) => {
      updatedStreams[streamId] = {
        ...updatedStreams[streamId],
        ...outlet,
        id: streamId,
        name: updatedStreams[streamId].name,
        sourceUnitId: unitId,
        sourcePortId: portId,
        targetUnitId: updatedStreams[streamId].targetUnitId,
        targetPortId: updatedStreams[streamId].targetPortId,
      } as MaterialStream;
    };

    const readObjectValue = (
      objectType: 'stream' | 'unit',
      objectId: string,
      property: string,
    ): number | null => {
      if (objectType === 'stream') {
        const stream = updatedStreams[objectId];
        if (!stream) return null;
        const value = (stream as any)[property];
        return typeof value === 'number' ? value : null;
      }
      const unit = updatedUnits[objectId];
      if (!unit) return null;
      const unitValue = (unit as any)[property];
      if (typeof unitValue === 'number') return unitValue;
      const paramValue = (unit.params as any)?.[property];
      return typeof paramValue === 'number' ? paramValue : null;
    };

    // ── Flash feed streams ──
    for (const [id, stream] of Object.entries(updatedStreams)) {
      if (!stream.sourceUnitId) {
        const flash = flashPT(compounds, stream.moleFractions, stream.T_K, stream.P_Pa,
          thermo.fluidPackage, thermo.interactionParams);
        updatedStreams[id] = {
          ...stream,
          phase: flash.phase,
          vaporFraction: flash.vaporFraction,
          H_Jpmol: flash.H_Jpmol,
          x_liquid: flash.x,
          y_vapor: flash.y,
          solved: true,
        };
      }
    }

    // ── Tear stream detection (DFS cycle detection) ──
    const tearStreamIds = new Set<string>();
    const adj = new Map<string, string[]>();
    for (const [, stream] of Object.entries(updatedStreams)) {
      const from = stream.sourceUnitId;
      const to = stream.targetUnitId;
      if (from && to) {
        if (!adj.has(from)) adj.set(from, []);
        adj.get(from)!.push(to);
      }
    }
    const visited = new Set<string>();
    const inStack = new Set<string>();
    const dfs = (u: string) => {
      visited.add(u);
      inStack.add(u);
      for (const v of adj.get(u) ?? []) {
        if (inStack.has(v)) {
          for (const [sid, s] of Object.entries(updatedStreams)) {
            if (s.sourceUnitId === u && s.targetUnitId === v) tearStreamIds.add(sid);
          }
        } else if (!visited.has(v)) {
          dfs(v);
        }
      }
      inStack.delete(u);
    };
    for (const unitId of Object.keys(updatedUnits)) {
      if (!visited.has(unitId)) dfs(unitId);
    }

    // Initialize tear streams (Aspen Plus approach):
    // If a tear stream has no flow, estimate from feed streams entering
    // the recycle loop. This avoids hardcoding initial guesses.
    for (const tearId of tearStreamIds) {
      const s = updatedStreams[tearId];
      if (!s.solved || s.totalFlow_molps < 1e-12) {
        // Estimate tear stream initial values from total feed to the loop.
        // Find all feed streams (streams without sourceUnitId, or from units outside the loop).
        const loopUnits = new Set<string>();
        // Trace the recycle loop from tear target back to tear source
        const target = s.targetUnitId;
        const source = s.sourceUnitId;
        if (target && source) {
          // BFS from target to source to find loop units
          const queue = [target];
          const bfsVisited = new Set<string>([target]);
          while (queue.length > 0) {
            const u = queue.shift()!;
            loopUnits.add(u);
            if (u === source) continue; // don't traverse past the source
            for (const v of adj.get(u) ?? []) {
              if (!bfsVisited.has(v)) {
                bfsVisited.add(v);
                queue.push(v);
              }
            }
          }
        }

        // Sum all external feeds entering the loop
        let totalFeed = 0;
        const feedComp = new Array(compounds.length).fill(0);
        let feedT = 0;
        let feedP = 0;
        let feedCount = 0;
        for (const [sid, str] of Object.entries(updatedStreams)) {
          if (sid === tearId) continue;
          if (str.targetUnitId && loopUnits.has(str.targetUnitId)) {
            const fromOutside = !str.sourceUnitId || !loopUnits.has(str.sourceUnitId);
            if (fromOutside && str.totalFlow_molps > 0) {
              totalFeed += str.totalFlow_molps;
              for (let ci = 0; ci < compounds.length; ci++) {
                feedComp[ci] += str.totalFlow_molps * str.moleFractions[ci];
              }
              feedT += str.T_K;
              feedP += str.P_Pa;
              feedCount++;
            }
          }
        }

        // Initial guess: fraction of total feed (Aspen uses ~20-50% as rule of thumb)
        const estFlow = Math.max(totalFeed * 0.2, 1e-6);
        const estComp = totalFeed > 0
          ? feedComp.map(c => c / totalFeed)
          : new Array(compounds.length).fill(1 / compounds.length);
        const estT = feedCount > 0 ? feedT / feedCount : 298.15;
        const estP = feedCount > 0 ? feedP / feedCount : 101325;

        // Flash the estimated tear stream
        const flash = flashPT(compounds, estComp, estT, estP,
          thermo.fluidPackage, thermo.interactionParams);
        updatedStreams[tearId] = {
          ...s,
          T_K: estT,
          P_Pa: estP,
          totalFlow_molps: estFlow,
          moleFractions: estComp,
          phase: flash.phase,
          vaporFraction: flash.vaporFraction,
          H_Jpmol: flash.H_Jpmol,
          x_liquid: flash.x,
          y_vapor: flash.y,
          solved: true,
        };
      }
    }

    // ── Forward sweep: solve all units in topological order ──
    const runForwardSweep = () => {
      const sweepErrors: string[] = [];
      const solvedUnits = new Set<string>();
      // Mark all units unsolved before sweep
      for (const uid of Object.keys(updatedUnits)) {
        updatedUnits[uid] = { ...updatedUnits[uid], solved: false };
      }
      let changed = true;
      let maxPasses = 10;

      while (changed && maxPasses-- > 0) {
        changed = false;
        for (const [unitId, unit] of Object.entries(updatedUnits)) {
          if (solvedUnits.has(unitId)) continue;

          const inletPorts = unit.ports.filter(p => p.type === 'inlet');
          const inletStreams = inletPorts
            .map(p => p.streamId ? updatedStreams[p.streamId] : null)
            .filter((s): s is MaterialStream => s !== null && s.solved);

          if (inletStreams.length < inletPorts.length) continue;

          try {
            switch (unit.type) {
              case 'Heater': {
                const result = solveHeater(compounds, inletStreams[0], {
                  targetT_K: unit.params.targetT_K as number,
                  outletP_Pa: unit.params.outletP_Pa as number,
                  waterSaturatedFeed: Boolean(unit.params.waterSaturatedFeed),
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = {
                  ...unit,
                  duty_W: result.duty_W,
                  solved: true,
                  errors: [],
                  params: {
                    ...unit.params,
                    waterAddedForSaturation_molps: result.waterAddedForSaturation_molps,
                    excessFreeWater_molps: result.excessFreeWater_molps,
                  },
                };
                break;
              }

              case 'Flash': {
                const result = solveFlashDrum(compounds, inletStreams[0], {
                  T_K: unit.params.T_K as number,
                  P_Pa: unit.params.P_Pa as number,
                  waterSaturatedFeed: Boolean(unit.params.waterSaturatedFeed),
                }, thermo);
                const vapPort = unit.ports.find(p => p.id.includes('vap') || p.label === 'Vapor');
                const liqPort = unit.ports.find(p => p.id.includes('liq') || p.label === 'Liquid');
                if (vapPort?.streamId) applyOutlet(vapPort.streamId, result.vapor, unitId, vapPort.id);
                if (liqPort?.streamId) applyOutlet(liqPort.streamId, result.liquid, unitId, liqPort.id);
                updatedUnits[unitId] = {
                  ...unit,
                  duty_W: result.duty_W,
                  solved: true,
                  errors: [],
                  params: {
                    ...unit.params,
                    waterAddedForSaturation_molps: result.waterAddedForSaturation_molps,
                    excessFreeWater_molps: result.excessFreeWater_molps,
                  },
                };
                break;
              }

              case 'Mixer': {
                const result = solveMixer(compounds, inletStreams,
                  (unit.params.outletP_Pa as number) ?? inletStreams[0].P_Pa, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: 0, solved: true, errors: [] };
                break;
              }

              case 'Valve': {
                const result = solveValve(compounds, inletStreams[0], {
                  outletP_Pa: unit.params.outletP_Pa as number,
                  waterSaturatedFeed: Boolean(unit.params.waterSaturatedFeed),
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = {
                  ...unit,
                  duty_W: 0,
                  solved: true,
                  errors: [],
                  params: {
                    ...unit.params,
                    waterAddedForSaturation_molps: result.waterAddedForSaturation_molps,
                    excessFreeWater_molps: result.excessFreeWater_molps,
                  },
                };
                break;
              }

              case 'Pump': {
                const result = solvePump(compounds, inletStreams[0], {
                  outletP_Pa: unit.params.outletP_Pa as number,
                  efficiency: (unit.params.efficiency as number) ?? 0.75,
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.work_W, solved: true, errors: [] };
                break;
              }

              case 'Compressor': {
                const result = solveCompressor(compounds, inletStreams[0], {
                  outletP_Pa: unit.params.outletP_Pa as number,
                  efficiency: (unit.params.efficiency as number) ?? 0.72,
                  model: (unit.params.model as 'isentropic' | 'polytropic') ?? 'isentropic',
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.work_W, solved: true, errors: [] };
                break;
              }

              case 'Splitter': {
                const ratio = (unit.params.splitRatio as number) ?? 0.5;
                const result = solveSplitter(inletStreams[0], [ratio, 1 - ratio]);
                const outPorts = unit.ports.filter(p => p.type === 'outlet');
                for (let oi = 0; oi < outPorts.length && oi < result.length; oi++) {
                  const streamId = outPorts[oi].streamId;
                  if (streamId) applyOutlet(streamId, result[oi], unitId, outPorts[oi].id);
                }
                updatedUnits[unitId] = { ...unit, duty_W: 0, solved: true, errors: [] };
                break;
              }

              case 'HeatX': {
                const hotPort = unit.ports.find(p => p.id.includes('hot-in') || p.label === 'Hot Inlet');
                const coldPort = unit.ports.find(p => p.id.includes('cold-in') || p.label === 'Cold Inlet');
                const hotStream = hotPort?.streamId ? updatedStreams[hotPort.streamId] : null;
                const coldStream = coldPort?.streamId ? updatedStreams[coldPort.streamId] : null;
                if (!hotStream?.solved || !coldStream?.solved) continue;

                const result = solveHeatX(compounds, hotStream, coldStream, {
                  spec: (unit.params.spec as 'hotOutletT' | 'coldOutletT' | 'duty' | 'UA') ?? 'hotOutletT',
                  specValue: unit.params.specValue as number,
                  flowArrangement: (unit.params.flowArrangement as 'counter' | 'co' | '1shell2tube') ?? 'counter',
                }, thermo);
                const hotOutPort = unit.ports.find(p => p.id.includes('hot-out') || p.label === 'Hot Outlet');
                const coldOutPort = unit.ports.find(p => p.id.includes('cold-out') || p.label === 'Cold Outlet');
                if (hotOutPort?.streamId) applyOutlet(hotOutPort.streamId, result.hotOutlet, unitId, hotOutPort.id);
                if (coldOutPort?.streamId) applyOutlet(coldOutPort.streamId, result.coldOutlet, unitId, coldOutPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'Pipe': {
                const result = solvePipe(compounds, inletStreams[0], {
                  length_m: (unit.params.length_m as number) ?? 100,
                  diameter_m: (unit.params.diameter_m as number) ?? 0.1,
                  roughness_m: unit.params.roughness_m as number | undefined,
                  elevation_m: unit.params.elevation_m as number | undefined,
                  K_fittings: unit.params.K_fittings as number | undefined,
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: 0, solved: true, errors: [] };
                break;
              }

              case 'RStoic': {
                const result = solveRStoic(compounds, inletStreams[0], {
                  reactions: (unit.params.reactions as any[]) ?? [],
                  outletT_K: unit.params.outletT_K as number,
                  outletP_Pa: unit.params.outletP_Pa as number,
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'Column':
              case 'RadFrac':
              case 'RateFrac': {
              if (unit.type === 'RateFrac') {
                  const specDiagnostics = getRigorousColumnSpecDiagnostics({
                    nStages: (unit.params.nStages as number) ?? 20,
                    feedStage: (unit.params.feedStage as number) ?? 10,
                    refluxRatio: (unit.params.refluxRatio as number) ?? 2,
                    distillateRate_molps: (unit.params.distillateRate_molps as number) ?? inletStreams[0].totalFlow_molps / 2,
                    condenserP_Pa: (unit.params.condenserP_Pa as number) ?? 101325,
                    reboilerP_Pa: (unit.params.reboilerP_Pa as number) ?? 101325,
                    condenserType: (unit.params.condenserType as 'Total' | 'Partial') ?? 'Total',
                    spec1Mode: unit.params.spec1Mode as any,
                    spec1Value: unit.params.spec1Value as number | undefined,
                    spec1ComponentIndex: unit.params.spec1ComponentIndex as number | undefined,
                    spec2Mode: unit.params.spec2Mode as any,
                    spec2Value: unit.params.spec2Value as number | undefined,
                    spec2ComponentIndex: unit.params.spec2ComponentIndex as number | undefined,
                  });
                  if (specDiagnostics.status !== 'ok') {
                    updatedUnits[unitId] = { ...unit, solved: false, errors: specDiagnostics.messages };
                    solvedUnits.add(unitId);
                    changed = true;
                    break;
                  }
                  const rateResult = solveRateFrac(compounds, inletStreams[0], {
                    nStages: (unit.params.nStages as number) ?? 20,
                    feedStage: (unit.params.feedStage as number) ?? 10,
                    refluxRatio: (unit.params.refluxRatio as number) ?? 2.0,
                    distillateRate_molps: (unit.params.distillateRate_molps as number) ?? inletStreams[0].totalFlow_molps / 2,
                    liquidSideDrawStage: unit.params.liquidSideDrawStage as number | undefined,
                    liquidSideDrawFraction: (unit.params.liquidSideDrawFraction as number) ?? 0,
                    vaporSideDrawStage: unit.params.vaporSideDrawStage as number | undefined,
                    vaporSideDrawFraction: (unit.params.vaporSideDrawFraction as number) ?? 0,
                    spec1Mode: unit.params.spec1Mode as any,
                    spec1Value: unit.params.spec1Value as number | undefined,
                    spec1ComponentIndex: unit.params.spec1ComponentIndex as number | undefined,
                    spec2Mode: unit.params.spec2Mode as any,
                    spec2Value: unit.params.spec2Value as number | undefined,
                    spec2ComponentIndex: unit.params.spec2ComponentIndex as number | undefined,
                    condenserP_Pa: (unit.params.condenserP_Pa as number) ?? 101325,
                    reboilerP_Pa: (unit.params.reboilerP_Pa as number) ?? 101325,
                    condenserType: (unit.params.condenserType as 'Total' | 'Partial') ?? 'Total',
                    columnDiameter_m: unit.params.columnDiameter_m as number | undefined,
                    traySpacing_m: unit.params.traySpacing_m as number | undefined,
                    weirHeight_m: unit.params.weirHeight_m as number | undefined,
                    interfacialArea_m2pm3: unit.params.interfacialArea_m2pm3 as number | undefined,
                    kL_mps: unit.params.kL_mps as number | undefined,
                    kV_mps: unit.params.kV_mps as number | undefined,
                    liquidHoldup_s: unit.params.liquidHoldup_s as number | undefined,
                    vaporHoldup_s: unit.params.vaporHoldup_s as number | undefined,
                    internalsMode: (unit.params.internalsMode as 'Tray' | 'Packing' | undefined) ?? 'Tray',
                    pressureRelaxation: unit.params.pressureRelaxation as number | undefined,
                  }, thermo);
                  const distPort = unit.ports.find(p => p.id.includes('dist') || p.label === 'Distillate');
                  const sideVapPort = unit.ports.find(p => p.id.includes('side-vap') || p.label === 'Side Vapor');
                  const sideLiqPort = unit.ports.find(p => p.id.includes('side-liq') || p.label === 'Side Liquid');
                  const botPort = unit.ports.find(p => p.id.includes('bot') || p.label === 'Bottoms');
                  if (distPort?.streamId) applyOutlet(distPort.streamId, rateResult.distillate, unitId, distPort.id);
                  if (sideVapPort?.streamId && rateResult.vaporSideDraw) applyOutlet(sideVapPort.streamId, rateResult.vaporSideDraw, unitId, sideVapPort.id);
                  if (sideLiqPort?.streamId && rateResult.liquidSideDraw) applyOutlet(sideLiqPort.streamId, rateResult.liquidSideDraw, unitId, sideLiqPort.id);
                  if (botPort?.streamId) applyOutlet(botPort.streamId, rateResult.bottoms, unitId, botPort.id);
                  updatedUnits[unitId] = {
                    ...unit,
                    params: {
                      ...unit.params,
                      effectiveMurphreeEfficiency: rateResult.effectiveMurphreeEfficiency,
                      floodingFractions: rateResult.floodingFractions,
                      stagePressureDrops_Pa: rateResult.stagePressureDrops_Pa,
                      stageEfficiencies: rateResult.stageEfficiencies,
                      interfacialLiquidCompositions: rateResult.interfacialLiquidCompositions,
                      interfacialVaporCompositions: rateResult.interfacialVaporCompositions,
                      liquidMassTransferCoefficients_mps: rateResult.liquidMassTransferCoefficients_mps,
                      vaporMassTransferCoefficients_mps: rateResult.vaporMassTransferCoefficients_mps,
                      pressureCouplingResidual_Pa: rateResult.pressureCouplingResidual_Pa,
                      hydraulicRegimes: rateResult.hydraulicRegimes,
                    },
                    duty_W: rateResult.duty_reboiler_W + rateResult.duty_condenser_W,
                    solved: true,
                    errors: rateResult.converged
                      ? []
                      : specDiagnostics.status !== 'ok'
                        ? specDiagnostics.messages
                        : ['RateFrac did not converge'],
                  };
                } else if (unit.type === 'RadFrac' || unit.params.rigorousMode) {
                  // Rigorous MESH column (Wang-Henke Bubble-Point method)
                  const specDiagnostics = getRigorousColumnSpecDiagnostics({
                    nStages: (unit.params.nStages as number) ?? 20,
                    feedStage: (unit.params.feedStage as number) ?? 10,
                    refluxRatio: (unit.params.refluxRatio as number) ?? 2,
                    distillateRate_molps: (unit.params.distillateRate_molps as number) ?? inletStreams[0].totalFlow_molps / 2,
                    condenserP_Pa: (unit.params.condenserP_Pa as number) ?? 101325,
                    reboilerP_Pa: (unit.params.reboilerP_Pa as number) ?? 101325,
                    condenserType: (unit.params.condenserType as 'Total' | 'Partial') ?? 'Total',
                    spec1Mode: unit.params.spec1Mode as any,
                    spec1Value: unit.params.spec1Value as number | undefined,
                    spec1ComponentIndex: unit.params.spec1ComponentIndex as number | undefined,
                    spec2Mode: unit.params.spec2Mode as any,
                    spec2Value: unit.params.spec2Value as number | undefined,
                    spec2ComponentIndex: unit.params.spec2ComponentIndex as number | undefined,
                  });
                  if (specDiagnostics.status !== 'ok') {
                    updatedUnits[unitId] = { ...unit, solved: false, errors: specDiagnostics.messages };
                    solvedUnits.add(unitId);
                    changed = true;
                    break;
                  }
                  const rigResult = solveColumnRigorous(compounds, inletStreams[0], {
                    nStages: (unit.params.nStages as number) ?? 20,
                    feedStage: (unit.params.feedStage as number) ?? 10,
                    refluxRatio: (unit.params.refluxRatio as number) ?? 2.0,
                    distillateRate_molps: (unit.params.distillateRate_molps as number) ?? inletStreams[0].totalFlow_molps / 2,
                    murphreeEfficiency: (unit.params.murphreeEfficiency as number) ?? 1.0,
                    liquidSideDrawStage: unit.params.liquidSideDrawStage as number | undefined,
                    liquidSideDrawFraction: (unit.params.liquidSideDrawFraction as number) ?? 0,
                    vaporSideDrawStage: unit.params.vaporSideDrawStage as number | undefined,
                    vaporSideDrawFraction: (unit.params.vaporSideDrawFraction as number) ?? 0,
                    spec1Mode: unit.params.spec1Mode as any,
                    spec1Value: unit.params.spec1Value as number | undefined,
                    spec1ComponentIndex: unit.params.spec1ComponentIndex as number | undefined,
                    spec2Mode: unit.params.spec2Mode as any,
                    spec2Value: unit.params.spec2Value as number | undefined,
                    spec2ComponentIndex: unit.params.spec2ComponentIndex as number | undefined,
                    condenserP_Pa: (unit.params.condenserP_Pa as number) ?? 101325,
                    reboilerP_Pa: (unit.params.reboilerP_Pa as number) ?? 101325,
                    condenserType: (unit.params.condenserType as 'Total' | 'Partial') ?? 'Total',
                  }, thermo);
                  const distPort = unit.ports.find(p => p.id.includes('dist') || p.label === 'Distillate');
                  const sideVapPort = unit.ports.find(p => p.id.includes('side-vap') || p.label === 'Side Vapor');
                  const sideLiqPort = unit.ports.find(p => p.id.includes('side-liq') || p.label === 'Side Liquid');
                  const botPort = unit.ports.find(p => p.id.includes('bot') || p.label === 'Bottoms');
                  if (distPort?.streamId) applyOutlet(distPort.streamId, rigResult.distillate, unitId, distPort.id);
                  if (sideVapPort?.streamId && rigResult.vaporSideDraw) applyOutlet(sideVapPort.streamId, rigResult.vaporSideDraw, unitId, sideVapPort.id);
                  if (sideLiqPort?.streamId && rigResult.liquidSideDraw) applyOutlet(sideLiqPort.streamId, rigResult.liquidSideDraw, unitId, sideLiqPort.id);
                  if (botPort?.streamId) applyOutlet(botPort.streamId, rigResult.bottoms, unitId, botPort.id);
                  updatedUnits[unitId] = {
                    ...unit,
                    duty_W: rigResult.duty_reboiler_W + rigResult.duty_condenser_W,
                    solved: true,
                    errors: rigResult.converged
                      ? []
                      : specDiagnostics.status !== 'ok'
                        ? specDiagnostics.messages
                        : ['MESH did not converge'],
                  };
                } else {
                  // Shortcut (Fenske-Underwood-Gilliland)
                  const result = solveColumn(compounds, inletStreams[0], {
                    lightKeyIndex: (unit.params.lightKeyIndex as number) ?? 0,
                    heavyKeyIndex: (unit.params.heavyKeyIndex as number) ?? 1,
                    lightKeyRecovery: (unit.params.lightKeyRecovery as number) ?? 0.95,
                    heavyKeyRecovery: (unit.params.heavyKeyRecovery as number) ?? 0.95,
                    refluxRatioMultiplier: (unit.params.refluxRatioMultiplier as number) ?? 1.3,
                    condenserP_Pa: (unit.params.condenserP_Pa as number) ?? 101325,
                    reboilerP_Pa: (unit.params.reboilerP_Pa as number) ?? 101325,
                    condenserType: (unit.params.condenserType as 'Total' | 'Partial') ?? 'Total',
                  }, thermo);
                  const distPort = unit.ports.find(p => p.id.includes('dist') || p.label === 'Distillate');
                  const botPort = unit.ports.find(p => p.id.includes('bot') || p.label === 'Bottoms');
                  if (distPort?.streamId) applyOutlet(distPort.streamId, result.distillate, unitId, distPort.id);
                  if (botPort?.streamId) applyOutlet(botPort.streamId, result.bottoms, unitId, botPort.id);
                  updatedUnits[unitId] = {
                    ...unit,
                    duty_W: result.duty_reboiler_W + result.duty_condenser_W,
                    solved: true,
                    errors: [],
                  };
                }
                break;
              }

              case 'ComponentSeparator': {
                const splitFracs = (unit.params.splitFractions as number[])
                  ?? new Array(compounds.length).fill(0.5);
                const result = solveComponentSeparator(compounds, inletStreams[0], splitFracs, thermo);
                const topPort = unit.ports.find(p => p.id.includes('top') || p.label === 'Top');
                const botPort2 = unit.ports.find(p => p.id.includes('bot') || p.label === 'Bottom');
                if (topPort?.streamId) applyOutlet(topPort.streamId, result.outlet1, unitId, topPort.id);
                if (botPort2?.streamId) applyOutlet(botPort2.streamId, result.outlet2, unitId, botPort2.id);
                updatedUnits[unitId] = { ...unit, duty_W: 0, solved: true, errors: [] };
                break;
              }

              case 'CSTR': {
                const result = solveCSTR(compounds, inletStreams[0], {
                  reactions: (unit.params.reactions as any[]) ?? [],
                  outletT_K: unit.params.outletT_K as number,
                  outletP_Pa: unit.params.outletP_Pa as number,
                }, (unit.params.duty_W as number) ?? 0, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'PFR': {
                const result = solvePFR(compounds, inletStreams[0], {
                  reactions: (unit.params.reactions as any[]) ?? [],
                  outletT_K: unit.params.outletT_K as number,
                  outletP_Pa: unit.params.outletP_Pa as number,
                }, (unit.params.nSteps as number) ?? 20, (unit.params.duty_W as number) ?? 0, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'ThreePhaseFlash': {
                const result = solveThreePhaseFlash(compounds, inletStreams[0],
                  (unit.params.T_K as number) ?? inletStreams[0].T_K,
                  (unit.params.P_Pa as number) ?? inletStreams[0].P_Pa,
                  (unit.params.duty_W as number) ?? 0, thermo);
                const vapPort = unit.ports.find(p => p.id.includes('vap') || p.label === 'Vapor');
                const liq1Port = unit.ports.find(p => p.id.includes('liq1') || p.label === 'Liquid I');
                const liq2Port = unit.ports.find(p => p.id.includes('liq2') || p.label === 'Liquid II');
                if (vapPort?.streamId) applyOutlet(vapPort.streamId, result.vapor, unitId, vapPort.id);
                if (liq1Port?.streamId) applyOutlet(liq1Port.streamId, result.liquidI, unitId, liq1Port.id);
                if (liq2Port?.streamId) applyOutlet(liq2Port.streamId, result.liquidII, unitId, liq2Port.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'Absorber': {
                const gasPort = unit.ports.find(p => p.id.includes('gas-in') || p.label === 'Gas Inlet');
                const liqPort = unit.ports.find(p => p.id.includes('liq-in') || p.label === 'Liquid Inlet');
                const gasStream = gasPort?.streamId ? updatedStreams[gasPort.streamId] : null;
                const liqStream = liqPort?.streamId ? updatedStreams[liqPort.streamId] : null;
                if (!gasStream?.solved || !liqStream?.solved) continue;

                const result = solveAbsorber(compounds, gasStream, liqStream,
                  (unit.params.nStages as number) ?? 10,
                  (unit.params.P_Pa as number) ?? gasStream.P_Pa, thermo);
                const gasOutPort = unit.ports.find(p => p.id.includes('gas-out') || p.label === 'Gas Outlet');
                const liqOutPort = unit.ports.find(p => p.id.includes('liq-out') || p.label === 'Liquid Outlet');
                if (gasOutPort?.streamId) applyOutlet(gasOutPort.streamId, result.gasOutlet, unitId, gasOutPort.id);
                if (liqOutPort?.streamId) applyOutlet(liqOutPort.streamId, result.liquidOutlet, unitId, liqOutPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: 0, solved: true, errors: [] };
                break;
              }

              case 'Extractor': {
                const feedPort = unit.ports.find(p => p.id.includes('feed-in') || p.label === 'Feed Inlet');
                const solventPort = unit.ports.find(p => p.id.includes('solvent-in') || p.label === 'Solvent Inlet');
                const feedStream = feedPort?.streamId ? updatedStreams[feedPort.streamId] : null;
                const solventStream = solventPort?.streamId ? updatedStreams[solventPort.streamId] : null;
                if (!feedStream?.solved || !solventStream?.solved) continue;

                const result = solveExtractor(
                  compounds,
                  feedStream,
                  solventStream,
                  (unit.params.nStages as number) ?? 6,
                  (unit.params.P_Pa as number) ?? feedStream.P_Pa,
                  thermo,
                  unit.params.solventKeyIndex as number | undefined,
                );
                const extractPort = unit.ports.find(p => p.id.includes('extract-out') || p.label === 'Extract');
                const raffinatePort = unit.ports.find(p => p.id.includes('raffinate-out') || p.label === 'Raffinate');
                if (extractPort?.streamId) applyOutlet(extractPort.streamId, result.extractOutlet, unitId, extractPort.id);
                if (raffinatePort?.streamId) applyOutlet(raffinatePort.streamId, result.raffinateOutlet, unitId, raffinatePort.id);
                updatedUnits[unitId] = {
                  ...unit,
                  duty_W: 0,
                  solved: true,
                  errors: result.converged ? [] : ['Extractor stage iteration did not fully converge'],
                };
                break;
              }

              case 'RYield': {
                const result = solveRYield(compounds, inletStreams[0],
                  (unit.params.yields as number[]) ?? compounds.map(() => 1 / compounds.length),
                  unit.params.outletT_K as number | undefined,
                  unit.params.outletP_Pa as number | undefined, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'Decanter': {
                const result = solveDecanter(compounds, inletStreams[0],
                  unit.params.T_K as number | undefined,
                  unit.params.P_Pa as number | undefined,
                  {
                    dirtyWaterMode: Boolean(unit.params.dirtyWaterMode),
                    freeWaterSplitFraction: unit.params.freeWaterSplitFraction as number | undefined,
                    waterCarryoverFraction: unit.params.waterCarryoverFraction as number | undefined,
                    hydrocarbonToAqueousFraction: unit.params.hydrocarbonToAqueousFraction as number | undefined,
                  },
                  thermo);
                const liq1Port = unit.ports.find(p => p.id.includes('liq1') || p.label === 'Liquid I');
                const liq2Port = unit.ports.find(p => p.id.includes('liq2') || p.label === 'Liquid II');
                if (liq1Port?.streamId) applyOutlet(liq1Port.streamId, result.liquidI, unitId, liq1Port.id);
                if (liq2Port?.streamId) applyOutlet(liq2Port.streamId, result.liquidII, unitId, liq2Port.id);
                updatedUnits[unitId] = {
                  ...unit,
                  duty_W: result.duty_W,
                  solved: true,
                  errors: [],
                  params: {
                    ...unit.params,
                    freeWaterSeparated_molps: result.freeWaterSeparated_molps,
                    aqueousPhaseIndex: result.aqueousPhaseIndex,
                    waterSaturationFraction: result.waterSaturationFraction,
                  },
                };
                break;
              }

              case 'REquil': {
                const result = solveREquil(compounds, inletStreams[0],
                  (unit.params.reactions as any[]) ?? [],
                  unit.params.outletT_K as number | undefined,
                  unit.params.outletP_Pa as number | undefined, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'RGibbs': {
                const result = solveRGibbs(compounds, inletStreams[0],
                  unit.params.outletT_K as number | undefined,
                  unit.params.outletP_Pa as number | undefined,
                  (unit.params.elementMatrix as number[][]) ?? [],
                  thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'MHeatX': {
                const result = solveMHeatX(compounds, inletStreams, {
                  hotStreamIndices: (unit.params.hotStreamIndices as number[]) ?? [0],
                  coldStreamIndices: (unit.params.coldStreamIndices as number[]) ?? [1],
                  deltaT_min_K: (unit.params.deltaT_min_K as number) ?? 10,
                  UA_WpK: unit.params.UA_WpK as number | undefined,
                }, thermo);
                const outPorts = unit.ports.filter(p => p.type === 'outlet');
                for (let oi = 0; oi < outPorts.length && oi < result.outlets.length; oi++) {
                  const streamId = outPorts[oi].streamId;
                  if (streamId) applyOutlet(streamId, result.outlets[oi], unitId, outPorts[oi].id);
                }
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'RBatch': {
                const result = solveRBatch(compounds, inletStreams[0], {
                  reactions: (unit.params.reactions as any[]) ?? [],
                  rateConstants_1ps: (unit.params.rateConstants_1ps as number[]) ?? [0.01],
                  Ea_Jpmol: (unit.params.Ea_Jpmol as number[]) ?? [0],
                  T_ref_K: (unit.params.T_ref_K as number) ?? 300,
                  batchTime_s: (unit.params.batchTime_s as number) ?? 3600,
                  nSteps: unit.params.nSteps as number | undefined,
                  isothermal: unit.params.isothermal as boolean | undefined,
                  duty_W: unit.params.duty_W as number | undefined,
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'Crystallizer': {
                const result = solveCrystallizer(compounds, inletStreams[0], {
                  residenceTime_s: (unit.params.residenceTime_s as number) ?? 3600,
                  soluteIndex: (unit.params.soluteIndex as number) ?? 0,
                  C_sat_molpm3: (unit.params.C_sat_molpm3 as number) ?? 100,
                  k_g: (unit.params.k_g as number) ?? 1e-8,
                  g_exponent: (unit.params.g_exponent as number) ?? 1,
                  k_b: (unit.params.k_b as number) ?? 1e10,
                  b_exponent: (unit.params.b_exponent as number) ?? 1,
                  k_v: unit.params.k_v as number | undefined,
                  rho_crystal_kgpm3: unit.params.rho_crystal_kgpm3 as number | undefined,
                  T_K: unit.params.T_K as number | undefined,
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'Crusher': {
                const result = solveCrusher(compounds, inletStreams[0], {
                  workIndex_kWhpton: (unit.params.workIndex_kWhpton as number) ?? 15,
                  feedSize_um: (unit.params.feedSize_um as number) ?? 25000,
                  productSize_um: (unit.params.productSize_um as number) ?? 5000,
                  solidFlow_kgps: (unit.params.solidFlow_kgps as number) ?? 1,
                });
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.power_W, solved: true, errors: [] };
                break;
              }

              case 'Dryer': {
                const result = solveDryer(compounds, inletStreams[0], {
                  X_in_kgpkg: (unit.params.X_in_kgpkg as number) ?? 0.5,
                  X_out_kgpkg: (unit.params.X_out_kgpkg as number) ?? 0.05,
                  X_c_kgpkg: (unit.params.X_c_kgpkg as number) ?? 0.2,
                  X_eq_kgpkg: unit.params.X_eq_kgpkg as number | undefined,
                  dryFlow_kgps: (unit.params.dryFlow_kgps as number) ?? 1,
                  gasT_in_K: (unit.params.gasT_in_K as number) ?? 423,
                  gasFlow_kgps: (unit.params.gasFlow_kgps as number) ?? 5,
                  Lv_Jpkg: unit.params.Lv_Jpkg as number | undefined,
                });
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'Membrane': {
                const result = solveMembrane(compounds, inletStreams[0], {
                  model: (unit.params.model as 'stage-cut' | 'rejection' | 'solution-diffusion' | undefined) ?? 'stage-cut',
                  stageCut: (unit.params.stageCut as number) ?? 0.25,
                  selectivities: (unit.params.selectivities as number[]) ?? [],
                  rejectionCoefficients: (unit.params.rejectionCoefficients as number[]) ?? [],
                  membraneArea_m2: unit.params.membraneArea_m2 as number | undefined,
                  permeanceGPU: (unit.params.permeanceGPU as number[]) ?? [],
                  foulingFactor: unit.params.foulingFactor as number | undefined,
                  permeateP_Pa: unit.params.permeateP_Pa as number | undefined,
                  retentateP_Pa: unit.params.retentateP_Pa as number | undefined,
                }, thermo);
                const permPort = unit.ports.find(p => p.id.includes('perm') || p.label === 'Permeate');
                const retPort = unit.ports.find(p => p.id.includes('ret') || p.label === 'Retentate');
                if (permPort?.streamId) applyOutlet(permPort.streamId, result.permeate, unitId, permPort.id);
                if (retPort?.streamId) applyOutlet(retPort.streamId, result.retentate, unitId, retPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: 0, solved: true, errors: [] };
                break;
              }

              case 'Cyclone': {
                const result = solveCyclone(compounds, inletStreams[0], {
                  cutSize_um: (unit.params.cutSize_um as number) ?? 10,
                  particleSizes_um: (unit.params.particleSizes_um as number[]) ?? [],
                  pressureDrop_Pa: unit.params.pressureDrop_Pa as number | undefined,
                }, thermo);
                const ovhdPort = unit.ports.find(p => p.id.includes('ovhd') || p.label === 'Overflow');
                const uflowPort = unit.ports.find(p => p.id.includes('uflow') || p.label === 'Underflow');
                if (ovhdPort?.streamId) applyOutlet(ovhdPort.streamId, result.overflow, unitId, ovhdPort.id);
                if (uflowPort?.streamId) applyOutlet(uflowPort.streamId, result.underflow, unitId, uflowPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: 0, solved: true, errors: [] };
                break;
              }

              case 'Filter': {
                const result = solveFilter(compounds, inletStreams[0], {
                  retainedFractions: (unit.params.retainedFractions as number[]) ?? [],
                  pressureDrop_Pa: unit.params.pressureDrop_Pa as number | undefined,
                }, thermo);
                const filtratePort = unit.ports.find(p => p.id.includes('filtrate') || p.label === 'Filtrate');
                const cakePort = unit.ports.find(p => p.id.includes('cake') || p.label === 'Cake');
                if (filtratePort?.streamId) applyOutlet(filtratePort.streamId, result.filtrate, unitId, filtratePort.id);
                if (cakePort?.streamId) applyOutlet(cakePort.streamId, result.cake, unitId, cakePort.id);
                updatedUnits[unitId] = { ...unit, duty_W: 0, solved: true, errors: [] };
                break;
              }

              case 'Screen': {
                const result = solveScreen(compounds, inletStreams[0], {
                  cutSize_um: (unit.params.cutSize_um as number) ?? 5,
                  sharpness: unit.params.sharpness as number | undefined,
                  particleSizes_um: (unit.params.particleSizes_um as number[]) ?? [],
                }, thermo);
                const oversizePort = unit.ports.find(p => p.id.includes('oversize') || p.label === 'Oversize');
                const undersizePort = unit.ports.find(p => p.id.includes('undersize') || p.label === 'Undersize');
                if (oversizePort?.streamId) applyOutlet(oversizePort.streamId, result.oversize, unitId, oversizePort.id);
                if (undersizePort?.streamId) applyOutlet(undersizePort.streamId, result.undersize, unitId, undersizePort.id);
                updatedUnits[unitId] = { ...unit, duty_W: 0, solved: true, errors: [] };
                break;
              }

              case 'Centrifuge': {
                const result = solveCentrifuge(compounds, inletStreams[0], {
                  cutSize_um: (unit.params.cutSize_um as number) ?? 3,
                  sharpness: unit.params.sharpness as number | undefined,
                  particleSizes_um: (unit.params.particleSizes_um as number[]) ?? [],
                  liquidSplitToCake: unit.params.liquidSplitToCake as number | undefined,
                  pressureDrop_Pa: unit.params.pressureDrop_Pa as number | undefined,
                }, thermo);
                const centratePort = unit.ports.find(p => p.id.includes('centrate') || p.label === 'Centrate');
                const cakePort = unit.ports.find(p => p.id.includes('cake') || p.label === 'Cake');
                if (centratePort?.streamId) applyOutlet(centratePort.streamId, result.centrate, unitId, centratePort.id);
                if (cakePort?.streamId) applyOutlet(cakePort.streamId, result.cake, unitId, cakePort.id);
                updatedUnits[unitId] = { ...unit, duty_W: 0, solved: true, errors: [] };
                break;
              }

              case 'SteamDrum': {
                const result = solveSteamDrum(compounds, inletStreams[0], {
                  outletP_Pa: (unit.params.outletP_Pa as number) ?? inletStreams[0].P_Pa,
                }, thermo);
                const vaporPort = unit.ports.find(p => p.id.includes('vap') || p.label === 'Vapor');
                const liquidPort = unit.ports.find(p => p.id.includes('liq') || p.label === 'Liquid');
                if (vaporPort?.streamId) applyOutlet(vaporPort.streamId, result.vapor, unitId, vaporPort.id);
                if (liquidPort?.streamId) applyOutlet(liquidPort.streamId, result.liquid, unitId, liquidPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'SteamHeater': {
                const result = solveSteamHeater(compounds, inletStreams[0], {
                  targetT_K: (unit.params.targetT_K as number) ?? 423.15,
                  outletP_Pa: unit.params.outletP_Pa as number | undefined,
                  steamLevel: (unit.params.steamLevel as 'LP' | 'MP' | 'HP' | 'Auto' | undefined) ?? 'Auto',
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = {
                  ...unit,
                  params: { ...unit.params, steamFlow_kgps: result.steamFlow_kgps, utilityType: result.utilityType },
                  duty_W: result.duty_W,
                  solved: true,
                  errors: [],
                };
                break;
              }

              case 'SteamTurbine': {
                const result = solveSteamTurbine(compounds, inletStreams[0], {
                  outletP_Pa: (unit.params.outletP_Pa as number) ?? Math.max(1000, inletStreams[0].P_Pa * 0.3),
                  efficiency: unit.params.efficiency as number | undefined,
                  generatorEfficiency: unit.params.generatorEfficiency as number | undefined,
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = {
                  ...unit,
                  params: {
                    ...unit.params,
                    shaftWork_W: result.shaftWork_W,
                    electricPower_W: result.electricPower_W,
                    exhaustT_K: result.exhaustT_K,
                  },
                  duty_W: -result.electricPower_W,
                  solved: true,
                  errors: [],
                };
                break;
              }

              case 'SteamValve': {
                const result = solveSteamValve(compounds, inletStreams[0], {
                  outletP_Pa: (unit.params.outletP_Pa as number) ?? Math.max(1000, inletStreams[0].P_Pa * 0.5),
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = {
                  ...unit,
                  params: { ...unit.params, pressureDrop_Pa: result.pressureDrop_Pa },
                  duty_W: 0,
                  solved: true,
                  errors: [],
                };
                break;
              }

              case 'SteamHeader': {
                const result = solveSteamHeader(compounds, inletStreams, {
                  outletP_Pa: unit.params.outletP_Pa as number | undefined,
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
                break;
              }

              case 'SteamTrap': {
                const result = solveSteamTrap(compounds, inletStreams[0], {
                  outletP_Pa: (unit.params.outletP_Pa as number) ?? Math.max(1000, inletStreams[0].P_Pa * 0.5),
                  maxFlashFraction: unit.params.maxFlashFraction as number | undefined,
                }, thermo);
                const vaporPort = unit.ports.find(p => p.id.includes('vap') || p.label === 'Flash Vapor');
                const liquidPort = unit.ports.find(p => p.id.includes('liq') || p.label === 'Condensate');
                if (vaporPort?.streamId) applyOutlet(vaporPort.streamId, result.vaporOutlet, unitId, vaporPort.id);
                if (liquidPort?.streamId) applyOutlet(liquidPort.streamId, result.liquidOutlet, unitId, liquidPort.id);
                updatedUnits[unitId] = {
                  ...unit,
                  params: {
                    ...unit.params,
                    flashFraction: result.flashFraction,
                    pressureDrop_Pa: result.pressureDrop_Pa,
                  },
                  duty_W: 0,
                  solved: true,
                  errors: [],
                };
                break;
              }

              case 'Analyzer': {
                const measuredObjectType = (unit.params.measuredObjectType as 'stream' | 'unit' | undefined) ?? 'stream';
                const measuredObjectId = String(unit.params.measuredObjectId ?? '');
                const sourceStream = inletStreams[0]
                  ?? (measuredObjectType === 'stream' && measuredObjectId ? updatedStreams[measuredObjectId] : undefined);

                if (!inletStreams[0] && measuredObjectType === 'stream' && measuredObjectId && (!sourceStream || !sourceStream.solved)) {
                  continue;
                }

                if (sourceStream) {
                  const result = solveAnalyzer(compounds, sourceStream, {
                    propertyCode: String(unit.params.propertyCode ?? 'TEMP'),
                    componentIndex: Number(unit.params.componentIndex ?? 0),
                  });
                  updatedUnits[unitId] = {
                    ...unit,
                    duty_W: 0,
                    solved: true,
                    errors: result.warnings,
                    params: {
                      ...unit.params,
                      signalValue: result.signalValue,
                      propertyLabel: result.propertyLabel,
                    },
                  };
                  break;
                }

                if (measuredObjectType === 'unit' && measuredObjectId) {
                  if (!updatedUnits[measuredObjectId]?.solved) continue;
                  const result = readObjectValue('unit', measuredObjectId, String(unit.params.propertyCode ?? 'signalValue'));
                  if (result == null) {
                    updatedUnits[unitId] = {
                      ...unit,
                      duty_W: 0,
                      solved: false,
                      errors: ['Analyzer could not read the selected unit property.'],
                    };
                    break;
                  }
                  updatedUnits[unitId] = {
                    ...unit,
                    duty_W: 0,
                    solved: true,
                    errors: [],
                    params: {
                      ...unit.params,
                      signalValue: result,
                      propertyLabel: String(unit.params.propertyCode ?? 'signalValue'),
                    },
                  };
                  break;
                }

                updatedUnits[unitId] = {
                  ...unit,
                  duty_W: 0,
                  solved: false,
                  errors: ['Analyzer requires a connected stream or a valid measured object reference.'],
                };
                break;
              }

              case 'ChargeBalance': {
                const result = solveChargeBalance(compounds, inletStreams[0], {
                  adjustComponentIndex: Number(unit.params.adjustComponentIndex ?? 0),
                  targetCharge: Number(unit.params.targetCharge ?? 0),
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = {
                  ...unit,
                  duty_W: 0,
                  solved: true,
                  errors: [],
                  params: {
                    ...unit.params,
                    residualCharge: result.residualCharge,
                    adjustedComponentFlow_molps: result.adjustedComponentFlow_molps,
                  },
                };
                break;
              }

              case 'ElectrolyteEquilibrium': {
                let reactionSet: Array<{
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
                }> | undefined;
                try {
                  const rawReactionSet = typeof unit.params.reactionSetJson === 'string'
                    ? JSON.parse(unit.params.reactionSetJson)
                    : [];
                  reactionSet = Array.isArray(rawReactionSet) ? rawReactionSet : undefined;
                } catch {
                  reactionSet = undefined;
                }
                const result = solveElectrolyteEquilibrium(compounds, inletStreams[0], {
                  apparentSpeciesIndex: Number(unit.params.apparentSpeciesIndex ?? 0),
                  cationIndex: Number(unit.params.cationIndex ?? 0),
                  anionIndex: Number(unit.params.anionIndex ?? 0),
                  cationStoich: Number(unit.params.cationStoich ?? 1),
                  anionStoich: Number(unit.params.anionStoich ?? 1),
                  log10K: Number(unit.params.log10K ?? 2),
                  maxDissociationFraction: Number(unit.params.maxDissociationFraction ?? 0.999),
                  reactionSet,
                }, thermo);
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId) applyOutlet(outPort.streamId, result.outlet, unitId, outPort.id);
                updatedUnits[unitId] = {
                  ...unit,
                  duty_W: 0,
                  solved: true,
                  errors: [],
                  params: {
                    ...unit.params,
                    dissociationFraction: result.dissociationFraction,
                    ionicStrength: result.ionicStrength,
                    equilibriumConstant: result.equilibriumConstant,
                    activeReactionCount: result.activeReactionCount,
                    equilibriumPasses: result.equilibriumPasses,
                    reactionExtents_molps: result.reactionExtents_molps,
                  },
                };
                break;
              }

              case 'PIDController':
              case 'LeadLag':
              case 'DeadTime':
              case 'SignalSelector': {
                const passthrough = inletStreams[0]
                  ? { ...inletStreams[0], solved: true }
                  : undefined;
                const outPort = unit.ports.find(p => p.type === 'outlet');
                if (outPort?.streamId && passthrough) applyOutlet(outPort.streamId, passthrough, unitId, outPort.id);
                updatedUnits[unitId] = {
                  ...unit,
                  duty_W: 0,
                  solved: true,
                  errors: [],
                };
                break;
              }

              default:
                sweepErrors.push(`Unknown unit type: ${unit.type}`);
                continue;
            }
            solvedUnits.add(unitId);
            changed = true;
          } catch (e: any) {
            sweepErrors.push(`Error solving ${unit.name}: ${e.message}`);
            updatedUnits[unitId] = { ...unit, solved: false, errors: [e.message] };
          }
        }
      }
      return sweepErrors;
    };

    // ── Recycle convergence loop (Wegstein + Aitken Δ² acceleration) ──
    const tearState = new Map<string, {
      x_prev?: number[];
      g_prev?: number[];
      x_hist: number[][];  // last 3 x values for Aitken
      iter: number;
    }>();
    for (const tid of tearStreamIds) {
      tearState.set(tid, { iter: 0, x_hist: [] });
    }

    const maxRecycleIter = tearStreamIds.size > 0 ? 300 : 1;
    let convergedRecycle = tearStreamIds.size === 0;

    for (let recycleIter = 0; recycleIter < maxRecycleIter; recycleIter++) {
      const tearBefore = new Map<string, number[]>();
      for (const tid of tearStreamIds) {
        const s = updatedStreams[tid];
        tearBefore.set(tid, [s.T_K, s.P_Pa, s.totalFlow_molps, ...s.moleFractions]);
      }

      const sweepErrors = runForwardSweep();
      errors.push(...sweepErrors);

      if (tearStreamIds.size === 0) break;

      // Check tear stream convergence
      let maxRelErr = 0;
      for (const tid of tearStreamIds) {
        const before = tearBefore.get(tid)!;
        const s = updatedStreams[tid];
        const after = [s.T_K, s.P_Pa, s.totalFlow_molps, ...s.moleFractions];
        for (let i = 0; i < before.length; i++) {
          const ref = Math.max(Math.abs(before[i]), 1e-10);
          maxRelErr = Math.max(maxRelErr, Math.abs(after[i] - before[i]) / ref);
        }
      }

      if (maxRelErr < 1e-6) {
        convergedRecycle = true;
        break;
      }

      // Convergence acceleration strategy: Iterated Aitken Δ²
      // Cycles of 3 direct substitutions + 1 Aitken extrapolation.
      // For near-unity loop gains (high-recycle systems), this converges
      // quadratically vs Wegstein's linear convergence.
      // Fallback to Wegstein q ∈ [-5, 0] when Aitken produces bad estimates.
      for (const tid of tearStreamIds) {
        const state = tearState.get(tid)!;
        const s = updatedStreams[tid];
        const x_curr = tearBefore.get(tid)!;
        const g_curr = [s.T_K, s.P_Pa, s.totalFlow_molps, ...s.moleFractions];

        // Track g(x) history for Aitken (collect 3 consecutive DS outputs)
        state.x_hist.push([...g_curr]);

        // Aitken cycle: every 4th iteration (after 3 DS steps), apply Aitken
        const cyclePos = state.iter % 4; // 0,1,2 = DS, 3 = Aitken
        if (cyclePos === 3 && state.x_hist.length >= 3) {
          const h = state.x_hist;
          const n0 = h[h.length - 3];
          const n1 = h[h.length - 2];
          const n2 = h[h.length - 1];
          const x_next = new Array(x_curr.length).fill(0);
          let aitkenOk = true;
          for (let i = 0; i < x_curr.length; i++) {
            const d2 = n2[i] - 2 * n1[i] + n0[i];
            if (Math.abs(d2) > 1e-15) {
              const d1 = n1[i] - n0[i];
              const x_aitken = n0[i] - (d1 * d1) / d2;
              if (isFinite(x_aitken) && (i < 3 || x_aitken >= 0)) {
                x_next[i] = x_aitken;
              } else {
                x_next[i] = g_curr[i];
                aitkenOk = false;
              }
            } else {
              x_next[i] = g_curr[i];
            }
          }

          // Enforce physical bounds
          x_next[0] = Math.max(100, x_next[0]); // T_K >= 100
          x_next[1] = Math.max(1000, x_next[1]); // P_Pa >= 1000
          x_next[2] = Math.max(1e-10, x_next[2]); // flow > 0

          // Normalize mole fractions
          const mf = x_next.slice(3);
          let sumMf = 0;
          for (let i = 0; i < mf.length; i++) {
            mf[i] = Math.max(0, mf[i]);
            sumMf += mf[i];
          }
          if (sumMf > 1e-15) {
            for (let i = 0; i < mf.length; i++) mf[i] /= sumMf;
          }

          updatedStreams[tid] = {
            ...s,
            T_K: x_next[0],
            P_Pa: x_next[1],
            totalFlow_molps: x_next[2],
            moleFractions: mf,
          };

          // Re-flash for consistency
          const tearFlash = flashPT(compounds, mf, x_next[0], x_next[1],
            thermo.fluidPackage, thermo.interactionParams);
          updatedStreams[tid].phase = tearFlash.phase;
          updatedStreams[tid].vaporFraction = tearFlash.vaporFraction;
          updatedStreams[tid].H_Jpmol = tearFlash.H_Jpmol;
          updatedStreams[tid].x_liquid = tearFlash.x;
          updatedStreams[tid].y_vapor = tearFlash.y;

          // Clear history for next Aitken cycle
          state.x_hist.length = 0;
        }
        // Direct substitution iterations (cyclePos 0,1,2): g(x) is already
        // in updatedStreams from the forward sweep, no action needed.

        state.x_prev = x_curr;
        state.g_prev = g_curr;
        state.iter++;
      }
    }

    if (!convergedRecycle && tearStreamIds.size > 0) {
      errors.push('Recycle loop did not converge within 300 iterations.');
    }

    // ── Design Spec Convergence (Secant Method) ──
    const signalBlockTypes = new Set(['LeadLag', 'DeadTime', 'SignalSelector']);
    const evaluateSignalBlocks = (): boolean => {
      let changed = false;
      for (const unit of Object.values(updatedUnits)) {
        if (!signalBlockTypes.has(unit.type) || !(((unit.params.active as boolean | undefined) ?? true))) continue;
        if (unit.type === 'SignalSelector') {
          const inputA = readObjectValue(
            (unit.params.measuredObjectTypeA as 'stream' | 'unit' | undefined) ?? 'stream',
            (unit.params.measuredObjectIdA as string | undefined) ?? '',
            (unit.params.measuredPropertyA as string | undefined) ?? 'T_K',
          );
          const inputB = readObjectValue(
            (unit.params.measuredObjectTypeB as 'stream' | 'unit' | undefined) ?? 'stream',
            (unit.params.measuredObjectIdB as string | undefined) ?? '',
            (unit.params.measuredPropertyB as string | undefined) ?? 'T_K',
          );
          const result = solveSignalSelector(
            [inputA ?? Number.NaN, inputB ?? Number.NaN],
            { mode: (unit.params.mode as 'High' | 'Low' | 'Average' | undefined) ?? 'High' },
          );
          if (Math.abs(((unit.params.signalValue as number | undefined) ?? 0) - result.signalValue) > 1e-9) changed = true;
          updatedUnits[unit.id] = {
            ...unit,
            solved: true,
            duty_W: 0,
            errors: [],
            params: { ...unit.params, signalValue: result.signalValue },
          };
          continue;
        }

        const measuredObjectType = (unit.params.measuredObjectType as 'stream' | 'unit' | undefined) ?? 'stream';
        const measuredObjectId = (unit.params.measuredObjectId as string | undefined) ?? '';
        const measuredProperty = (unit.params.measuredProperty as string | undefined) ?? 'T_K';
        if (!measuredObjectId) continue;
        const measuredValue = readObjectValue(measuredObjectType, measuredObjectId, measuredProperty);
        if (measuredValue == null) continue;

        if (unit.type === 'LeadLag') {
          const result = solveLeadLagSignal(measuredValue, {
            leadTime_s: unit.params.leadTime_s as number | undefined,
            lagTime_s: unit.params.lagTime_s as number | undefined,
            timestep_s: unit.params.timestep_s as number | undefined,
            previousInput: unit.params.previousInput as number | undefined,
            lastOutput: unit.params.lastOutput as number | undefined,
          });
          if (Math.abs(((unit.params.signalValue as number | undefined) ?? 0) - result.signalValue) > 1e-9) changed = true;
          updatedUnits[unit.id] = {
            ...unit,
            solved: true,
            duty_W: 0,
            errors: [],
            params: { ...unit.params, signalValue: result.signalValue, previousInput: result.previousInput, lastOutput: result.lastOutput },
          };
          continue;
        }

        const result = solveDeadTimeSignal(measuredValue, {
          deadTime_s: unit.params.deadTime_s as number | undefined,
          timestep_s: unit.params.timestep_s as number | undefined,
          history: unit.params.history as number[] | undefined,
        });
        if (Math.abs(((unit.params.signalValue as number | undefined) ?? 0) - result.signalValue) > 1e-9) changed = true;
        updatedUnits[unit.id] = {
          ...unit,
          solved: true,
          duty_W: 0,
          errors: [],
          params: { ...unit.params, signalValue: result.signalValue, history: result.history },
        };
      }
      return changed;
    };

    evaluateSignalBlocks();

    const activeControllers = Object.values(updatedUnits).filter(
      unit => unit.type === 'PIDController' && ((unit.params.active as boolean | undefined) ?? true),
    );
    if (activeControllers.length > 0) {
      for (let ctrlIter = 0; ctrlIter < 12; ctrlIter++) {
        evaluateSignalBlocks();
        let changedControllerOutput = false;

        for (const controller of activeControllers) {
          const measuredObjectType = (controller.params.measuredObjectType as 'stream' | 'unit' | undefined) ?? 'stream';
          const measuredObjectId = (controller.params.measuredObjectId as string | undefined) ?? '';
          const measuredProperty = (controller.params.measuredProperty as string | undefined) ?? 'T_K';
          const manipulatedUnitId = (controller.params.manipulatedUnitId as string | undefined) ?? '';
          const manipulatedParamName = (controller.params.manipulatedParamName as string | undefined) ?? '';
          if (!measuredObjectId || !manipulatedUnitId || !manipulatedParamName) continue;

          const measuredValue = readObjectValue(measuredObjectType, measuredObjectId, measuredProperty);
          const manipulatedUnit = updatedUnits[manipulatedUnitId];
          if (measuredValue == null || !manipulatedUnit) continue;

          const currentManipulated = (manipulatedUnit.params[manipulatedParamName] as number | undefined) ?? 0;
          const controlResult = solvePIDControllerSignal(measuredValue, {
            setpoint: (controller.params.setpoint as number) ?? 0,
            measuredProperty,
            gain: (controller.params.gain as number) ?? 1,
            integralTime_s: controller.params.integralTime_s as number | undefined,
            derivativeTime_s: controller.params.derivativeTime_s as number | undefined,
            timestep_s: (controller.params.controllerTimestep_s as number | undefined) ?? 1,
            measurementLag_s: controller.params.measurementLag_s as number | undefined,
            outputRateLimit_per_s: controller.params.outputRateLimit_per_s as number | undefined,
            bias: controller.params.bias as number | undefined,
            lowLimit: controller.params.lowLimit as number | undefined,
            highLimit: controller.params.highLimit as number | undefined,
            reverseActing: controller.params.reverseActing as boolean | undefined,
            previousError: controller.params.previousError as number | undefined,
            previousPreviousError: controller.params.previousPreviousError as number | undefined,
            previousMeasuredValue: controller.params.previousMeasuredValue as number | undefined,
            filteredMeasuredValue: controller.params.filteredMeasuredValue as number | undefined,
            lastOutput: (controller.params.lastOutput as number | undefined) ?? currentManipulated,
          });

          if (Math.abs(controlResult.controllerOutput - currentManipulated) > 1e-9) {
            changedControllerOutput = true;
          }

          updatedUnits[manipulatedUnitId] = {
            ...manipulatedUnit,
            params: {
              ...manipulatedUnit.params,
              [manipulatedParamName]: controlResult.controllerOutput,
            },
          };

          updatedUnits[controller.id] = {
            ...updatedUnits[controller.id],
            solved: true,
            duty_W: 0,
            errors: [],
            params: {
              ...updatedUnits[controller.id].params,
              controllerOutput: controlResult.controllerOutput,
              measuredValue: controlResult.measuredValue,
              filteredMeasuredValue: controlResult.filteredMeasuredValue,
              controlError: controlResult.error,
              lastOutput: controlResult.controllerOutput,
              previousError: controlResult.previousError,
              previousPreviousError: controlResult.previousPreviousError,
              previousMeasuredValue: controlResult.previousMeasuredValue,
            },
          };
        }

        if (!changedControllerOutput) break;
        runForwardSweep();
      }
    }

    const activeSpecs = designSpecs.filter(ds => ds.active);

    for (const spec of activeSpecs) {
      const readTarget = (): number | null => {
        if (spec.target.objectType === 'stream') {
          const s = updatedStreams[spec.target.objectId];
          if (s) return (s as any)[spec.target.property] ?? null;
        } else {
          const u = updatedUnits[spec.target.objectId];
          if (u) return (u as any)[spec.target.property] ?? (u.params as any)?.[spec.target.property] ?? null;
        }
        return null;
      };

      const setManipulated = (val: number) => {
        const u = updatedUnits[spec.manipulated.unitId];
        if (u) {
          updatedUnits[spec.manipulated.unitId] = {
            ...u,
            params: { ...u.params, [spec.manipulated.paramName]: val },
          };
        }
      };

      const getManipulated = (): number => {
        const u = updatedUnits[spec.manipulated.unitId];
        return (u?.params[spec.manipulated.paramName] as number) ?? 0;
      };

      let targetVal = readTarget();
      if (targetVal === null || typeof targetVal !== 'number') continue;

      let err = targetVal - spec.target.value;
      if (Math.abs(err) <= spec.tolerance) continue;

      // Secant method: need two initial points
      let x0 = getManipulated();
      const perturbFrac = 0.01;
      let x1 = x0 * (1 + perturbFrac) + (Math.abs(x0) < 1e-10 ? 1 : 0);
      // Clamp to bounds
      x1 = Math.max(spec.manipulated.lowerBound, Math.min(spec.manipulated.upperBound, x1));

      let f0 = err;

      // Evaluate at x1
      setManipulated(x1);
      // Re-flash feeds
      for (const [id, stream] of Object.entries(updatedStreams)) {
        if (!stream.sourceUnitId) {
          const flash = flashPT(compounds, stream.moleFractions, stream.T_K, stream.P_Pa,
            thermo.fluidPackage, thermo.interactionParams);
          updatedStreams[id] = { ...stream, phase: flash.phase, vaporFraction: flash.vaporFraction,
            H_Jpmol: flash.H_Jpmol, x_liquid: flash.x, y_vapor: flash.y, solved: true };
        }
      }
      runForwardSweep();
      let f1 = (readTarget() ?? 0) - spec.target.value;

      for (let dsIter = 0; dsIter < 20; dsIter++) {
        if (Math.abs(f1) <= spec.tolerance) break;

        // Secant step
        const dfdx = (f1 - f0);
        if (Math.abs(dfdx) < 1e-30) break;
        let x2 = x1 - f1 * (x1 - x0) / dfdx;
        // Clamp to bounds
        x2 = Math.max(spec.manipulated.lowerBound, Math.min(spec.manipulated.upperBound, x2));

        x0 = x1;
        f0 = f1;
        x1 = x2;

        setManipulated(x1);
        for (const [id, stream] of Object.entries(updatedStreams)) {
          if (!stream.sourceUnitId) {
            const flash = flashPT(compounds, stream.moleFractions, stream.T_K, stream.P_Pa,
              thermo.fluidPackage, thermo.interactionParams);
            updatedStreams[id] = { ...stream, phase: flash.phase, vaporFraction: flash.vaporFraction,
              H_Jpmol: flash.H_Jpmol, x_liquid: flash.x, y_vapor: flash.y, solved: true };
          }
        }
        runForwardSweep();
        f1 = (readTarget() ?? 0) - spec.target.value;
      }

      if (Math.abs(f1) > spec.tolerance) {
        errors.push(
          `Design spec "${spec.name}": target not met after 20 secant iterations ` +
          `(current=${(readTarget() ?? 0).toFixed(4)}, desired=${spec.target.value.toFixed(4)}).`
        );
      }
    }

    // ── Execute Calculator Blocks (after) ──
    for (const cb of calcBlocks.filter(c => c.executionOrder === 'after')) {
      execCalcBlock(cb, updatedStreams, updatedUnits);
    }

    set({
      streams: updatedStreams,
      units: updatedUnits,
      solved: errors.length === 0 && Object.values(updatedUnits).every(u => u.solved),
      solverErrors: errors,
    });
  },

  // Preload
  loadPreset: (preset) => {
    if (preset === 'cumene-production') {
      const { streams, units } = createPresetFlowsheet();
      set({
        compounds: [PROPYLENE, PROPANE, BENZENE, CUMENE, P_DIISOPROPYLBENZENE],
        fluidPackage: 'Ideal',
        unitSystem: { ...DEFAULT_UNIT_SYSTEM },
        streams,
        units,
        currentScreen: 'setup',
        solved: false,
        solverErrors: [],
        dbWarnings: [],
        openTabs: [],
        activeTab: null,
        nodePositions: {
          // Row 1: Feed → Reactor → Cooler
          'label-S1':  { x: 40,  y: 120 },
          'label-S2':  { x: 40,  y: 300 },
          'U1':        { x: 280, y: 210 },
          'U2':        { x: 580, y: 210 },
          'U3':        { x: 900, y: 210 },
          'U4':        { x: 1220, y: 210 },
          // Row 2: Flash → Column 1
          'U5':        { x: 1540, y: 210 },
          'label-S7':  { x: 1860, y: 60  },
          'U6':        { x: 1540, y: 560 },
          // Row 3: Column 2 + products
          'label-S9':  { x: 1120, y: 560 },
          'U7':        { x: 1540, y: 930 },
          'label-S11': { x: 1860, y: 860 },
          'label-S12': { x: 1860, y: 1040 },
        },
      });

      // Auto-run the flowsheet so user sees results immediately
      setTimeout(() => {
        get().solveAll();
      }, 0);
    }
  },

  // Node positions
  nodePositions: {},
  setNodePosition: (nodeId, pos) => set(s => ({
    nodePositions: { ...s.nodePositions, [nodeId]: pos },
  })),

  // ── Export / Import ──────────────────────────────────────────
  exportFlowsheet: () => {
    const s = get();
    return JSON.stringify({
      version: 1,
      compounds: s.compounds,
      fluidPackage: s.fluidPackage,
      interactionParams: s.interactionParams,
      unitSystem: s.unitSystem,
      streams: s.streams,
      units: s.units,
      nodePositions: s.nodePositions,
      designSpecs: s.designSpecs,
      sensitivities: s.sensitivities,
      calculatorBlocks: s.calculatorBlocks,
    }, null, 2);
  },

  importFlowsheet: (json: string) => {
    try {
      const data = JSON.parse(json);
      if (!data.version || !data.compounds || !data.streams || !data.units) return false;
      set({
        compounds: data.compounds,
        fluidPackage: data.fluidPackage ?? 'Ideal',
        interactionParams: data.interactionParams ?? {},
        unitSystem: data.unitSystem ?? DEFAULT_UNIT_SYSTEM,
        streams: data.streams,
        units: data.units,
        nodePositions: data.nodePositions ?? {},
        designSpecs: data.designSpecs ?? [],
        sensitivities: data.sensitivities ?? [],
        calculatorBlocks: data.calculatorBlocks ?? [],
        solved: false,
        solverErrors: [],
        dbWarnings: [],
        openTabs: [],
        activeTab: null,
      });
      return true;
    } catch {
      return false;
    }
  },

  // ── Undo / Redo ──────────────────────────────────────────────
  canUndo: false,
  canRedo: false,
  undo: () => {
    if (_undoStack.length === 0) return;
    const current = get();
    _redoStack.push({
      streams: current.streams,
      units: current.units,
      compounds: current.compounds,
      nodePositions: current.nodePositions,
    });
    const prev = _undoStack.pop()!;
    set({
      streams: prev.streams,
      units: prev.units,
      compounds: prev.compounds,
      nodePositions: prev.nodePositions,
      canUndo: _undoStack.length > 0,
      canRedo: true,
      solved: false,
      solverErrors: [],
      dbWarnings: [],
    });
  },
  redo: () => {
    if (_redoStack.length === 0) return;
    const current = get();
    _undoStack.push({
      streams: current.streams,
      units: current.units,
      compounds: current.compounds,
      nodePositions: current.nodePositions,
    });
    const next = _redoStack.pop()!;
    set({
      streams: next.streams,
      units: next.units,
      compounds: next.compounds,
      nodePositions: next.nodePositions,
      canUndo: true,
      canRedo: _redoStack.length > 0,
      solved: false,
      solverErrors: [],
      dbWarnings: [],
    });
  },
}));
