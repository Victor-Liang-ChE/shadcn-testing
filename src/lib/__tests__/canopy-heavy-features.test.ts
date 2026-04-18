import { describe, expect, it } from 'vitest';

import { PRESET_COMPOUNDS, useCanopyStore } from '../canopy/store';
import type { MaterialStream } from '../canopy/types';
import {
  getRigorousColumnSpecDiagnostics,
  solveColumnRigorous,
  solveRateFrac,
  solveSteamHeader,
  solveSteamTrap,
  solveSteamValve,
  solveLeadLagSignal,
  solveDeadTimeSignal,
  solveSignalSelector,
} from '../canopy/solver';
import { flashPT } from '../canopy/thermo';
import { createExampleCrudeAssay, generatePseudoComponentsFromAssay } from '../canopy/petro';

const benzene = PRESET_COMPOUNDS.find(c => c.name === 'BENZENE')!;
const cumene = PRESET_COMPOUNDS.find(c => c.name === 'CUMENE')!;
const dipb = PRESET_COMPOUNDS.find(c => c.name === 'P-DIISOPROPYLBENZENE')!;
const compounds = [benzene, cumene, dipb];

function makeStream(name: string, totalFlow_molps: number, moleFractions: number[], T_K = 360, P_Pa = 101325): MaterialStream {
  const flash = flashPT(compounds, moleFractions, T_K, P_Pa);
  return {
    id: name,
    name,
    T_K,
    P_Pa,
    totalFlow_molps,
    moleFractions,
    phase: flash.phase,
    vaporFraction: flash.vaporFraction,
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

describe('heavy canopy Aspen-style features', () => {
  it('RateFrac returns hydraulics and pressure-drop diagnostics on top of staged separation', () => {
    const feed = makeStream('ratefrac-feed', 120, [0.52, 0.33, 0.15], 400, 140000);
    const trayResult = solveRateFrac(compounds, feed, {
      nStages: 16,
      feedStage: 8,
      refluxRatio: 2.3,
      distillateRate_molps: 45,
      condenserP_Pa: 120000,
      reboilerP_Pa: 180000,
      condenserType: 'Total',
      columnDiameter_m: 1.8,
      traySpacing_m: 0.6,
      interfacialArea_m2pm3: 180,
      kL_mps: 3e-4,
      kV_mps: 8e-3,
      internalsMode: 'Tray',
      pressureRelaxation: 0.5,
    });
    const packingResult = solveRateFrac(compounds, feed, {
      nStages: 16,
      feedStage: 8,
      refluxRatio: 2.3,
      distillateRate_molps: 45,
      condenserP_Pa: 120000,
      reboilerP_Pa: 180000,
      condenserType: 'Total',
      columnDiameter_m: 1.8,
      traySpacing_m: 0.6,
      interfacialArea_m2pm3: 220,
      kL_mps: 3e-4,
      kV_mps: 8e-3,
      internalsMode: 'Packing',
      pressureRelaxation: 0.5,
    });

    expect(trayResult.stageEfficiencies).toHaveLength(16);
    expect(trayResult.floodingFractions).toHaveLength(16);
    expect(trayResult.stagePressureDrops_Pa).toHaveLength(16);
    expect(trayResult.effectiveMurphreeEfficiency).toBeGreaterThan(0.1);
    expect(trayResult.effectiveMurphreeEfficiency).toBeLessThan(0.95);
    expect(trayResult.stagePressureDrops_Pa.every(value => value > 0)).toBe(true);
    expect(trayResult.componentMolarFluxes_molpsm2).toHaveLength(16);
    expect(trayResult.componentMolarFluxes_molpsm2[0]).toHaveLength(compounds.length);
    expect(trayResult.stageHeatTransferRates_W).toHaveLength(16);
    expect(trayResult.stageHeatTransferRates_W.some(value => value >= 0)).toBe(true);
    expect(trayResult.couplingIterations).toBeGreaterThanOrEqual(1);
    expect(packingResult.hydraulicRegimes).toHaveLength(16);
    expect(packingResult.pressureCouplingResidual_Pa).toBeGreaterThanOrEqual(0);
    expect(packingResult.stagePressureDrops_Pa.reduce((sum, value) => sum + value, 0))
      .toBeLessThan(trayResult.stagePressureDrops_Pa.reduce((sum, value) => sum + value, 0));
  });

  it('tracks rigorous column under/over-specification and supports composition-driven specs', () => {
    const diagnostics = getRigorousColumnSpecDiagnostics({
      nStages: 16,
      feedStage: 8,
      refluxRatio: 2.3,
      distillateRate_molps: 45,
      condenserP_Pa: 120000,
      reboilerP_Pa: 180000,
      condenserType: 'Total',
      spec1Mode: 'DistillateRate',
      spec1Value: 45,
    });
    expect(diagnostics.status).toBe('underspecified');

    const feed = makeStream('column-feed', 120, [0.52, 0.33, 0.15], 400, 140000);
    const result = solveColumnRigorous(compounds, feed, {
      nStages: 16,
      feedStage: 8,
      refluxRatio: 2.3,
      distillateRate_molps: 45,
      condenserP_Pa: 120000,
      reboilerP_Pa: 180000,
      condenserType: 'Total',
      spec1Mode: 'DistillateRate',
      spec1Value: 45,
      spec2Mode: 'DistillateComposition',
      spec2Value: 0.55,
      spec2ComponentIndex: 0,
    });

    expect(result.converged).toBe(true);
    expect(result.distillate.totalFlow_molps).toBeCloseTo(45, 1);
    expect(result.distillate.moleFractions?.[0] ?? 0).toBeGreaterThan(0.52);
  });

  it('steam header mixes utility streams and steam valve performs letdown', () => {
    const steam1 = makeStream('steam-1', 20, [1, 0, 0], 520, 600000);
    const steam2 = makeStream('steam-2', 10, [1, 0, 0], 430, 350000);
    const header = solveSteamHeader(compounds, [steam1, steam2], { outletP_Pa: 300000 });
    const letdown = solveSteamValve(compounds, header.outlet as MaterialStream, { outletP_Pa: 180000 });

    expect(header.outlet.totalFlow_molps).toBeCloseTo(30, 8);
    expect(header.outlet.P_Pa).toBe(300000);
    expect(letdown.outlet.P_Pa).toBe(180000);
    expect(letdown.pressureDrop_Pa).toBeCloseTo(120000, 8);
  });

  it('steam trap produces flash vapor and condensate outlets during letdown', () => {
    const condensate = makeStream('condensate', 12, [1, 0, 0], 430, 500000);
    const result = solveSteamTrap(compounds, condensate, { outletP_Pa: 180000, maxFlashFraction: 0.25 });

    expect(result.flashFraction).toBeGreaterThanOrEqual(0);
    expect(result.flashFraction).toBeLessThanOrEqual(0.25);
    expect((result.vaporOutlet.totalFlow_molps ?? 0) + (result.liquidOutlet.totalFlow_molps ?? 0)).toBeCloseTo(12, 8);
    expect(result.pressureDrop_Pa).toBeCloseTo(320000, 8);
  });

  it('lead-lag, dead-time, and selector signals produce usable secondary measurements', () => {
    const leadLag = solveLeadLagSignal(120, {
      leadTime_s: 2,
      lagTime_s: 8,
      timestep_s: 1,
      previousInput: 100,
      lastOutput: 90,
    });
    const deadTime = solveDeadTimeSignal(55, {
      deadTime_s: 3,
      timestep_s: 1,
      history: [40, 42, 44, 46],
    });
    const selector = solveSignalSelector([leadLag.signalValue, deadTime.signalValue], { mode: 'High' });

    expect(leadLag.signalValue).toBeGreaterThan(90);
    expect(deadTime.signalValue).toBeCloseTo(42, 8);
    expect(selector.signalValue).toBeCloseTo(Math.max(leadLag.signalValue, deadTime.signalValue), 8);
  });

  it('generates petroleum pseudo-components from a crude assay without remote tables', () => {
    const assay = createExampleCrudeAssay();
    const pseudos = generatePseudoComponentsFromAssay(assay);

    expect(pseudos).toHaveLength(5);
    expect(pseudos[0].Tb_K!).toBeLessThan(pseudos[4].Tb_K!);
    expect(pseudos[0].molecularWeight).toBeLessThan(pseudos[4].molecularWeight);
    expect(pseudos.every(comp => comp.Pc_Pa > 0 && comp.Tc_K > comp.Tb_K!)).toBe(true);
  });

  it('adds RateFrac and signal family units to the palette', () => {
    useCanopyStore.setState(state => ({ ...state, units: {}, streams: {}, nodePositions: {} }));
    const store = useCanopyStore.getState();
    store.addUnitFromPalette('RateFrac', { x: 100, y: 100 });
    store.addUnitFromPalette('SteamValve', { x: 140, y: 120 });
    store.addUnitFromPalette('SteamHeader', { x: 180, y: 140 });
    store.addUnitFromPalette('SteamTrap', { x: 220, y: 160 });
    store.addUnitFromPalette('LeadLag', { x: 260, y: 180 });
    store.addUnitFromPalette('DeadTime', { x: 300, y: 200 });
    store.addUnitFromPalette('SignalSelector', { x: 340, y: 220 });

    const types = Object.values(useCanopyStore.getState().units).map(unit => unit.type);
    expect(types).toEqual(['RateFrac', 'SteamValve', 'SteamHeader', 'SteamTrap', 'LeadLag', 'DeadTime', 'SignalSelector']);
  });
});
