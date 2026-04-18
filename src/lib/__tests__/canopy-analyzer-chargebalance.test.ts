import { describe, expect, it } from 'vitest';

import { PRESET_COMPOUNDS, useCanopyStore } from '../canopy/store';
import { solveAnalyzer, solveChargeBalance, solveDecanter } from '../canopy/solver';
import { flashPT } from '../canopy/thermo';
import type { CanopyCompound, MaterialStream } from '../canopy/types';

const benzene = PRESET_COMPOUNDS.find(c => c.name === 'BENZENE')!;
const cumene = PRESET_COMPOUNDS.find(c => c.name === 'CUMENE')!;
const dipb = PRESET_COMPOUNDS.find(c => c.name === 'P-DIISOPROPYLBENZENE')!;
const water: CanopyCompound = {
  ...benzene,
  name: 'WATER',
  displayName: 'Water',
  molecularWeight: 18.015,
  Tc_K: 647.096,
  Pc_Pa: 22064000,
  omega: 0.3443,
  Tb_K: 373.15,
  specificGravity: 1,
  isElectrolyteSolvent: true,
  dielectricConst: { e1: 78.36, e2: 0, e3_K: 298.15 },
  antoine: {
    C1: 73.649,
    C2: -7258.2,
    C3: 0,
    C4: 0,
    C5: -7.3037,
    C6: 4.1653e-6,
    C7: 2,
    Tmin_K: 273.16,
    Tmax_K: 647.096,
  },
  waterSolubility: undefined,
  hcSolubility: undefined,
};

function makeStream(
  compounds: CanopyCompound[],
  totalFlow_molps: number,
  moleFractions: number[],
  T_K = 360,
  P_Pa = 101325,
): MaterialStream {
  const flash = flashPT(compounds, moleFractions, T_K, P_Pa);
  return {
    id: 'S1',
    name: 'S1',
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

describe('canopy analyzer and charge-balance blocks', () => {
  it('computes Aspen-style analyzer properties from a stream', () => {
    const compounds = [benzene, cumene];
    const stream = makeStream(compounds, 20, [0.65, 0.35], 365, 150000);

    const temp = solveAnalyzer(compounds, stream, { propertyCode: 'TEMP' });
    const massFlow = solveAnalyzer(compounds, stream, { propertyCode: 'MASS_FLOW' });
    const tbpt = solveAnalyzer(compounds, stream, { propertyCode: 'TBPT' });
    const vabp = solveAnalyzer(compounds, stream, { propertyCode: 'VABP' });
    const smx = solveAnalyzer(compounds, stream, { propertyCode: 'SMX' });
    const api = solveAnalyzer(compounds, stream, { propertyCode: 'API' });
    const refractiveIndex = solveAnalyzer(compounds, stream, { propertyCode: 'REFINDEX' });

    expect(temp.signalValue).toBeCloseTo(365, 8);
    expect(massFlow.signalValue).toBeGreaterThan(1);
    expect(tbpt.signalValue).toBeGreaterThan(benzene.Tb_K!);
    expect(tbpt.signalValue).toBeLessThan(cumene.Tb_K!);
    expect(vabp.signalValue).toBeGreaterThan(benzene.Tb_K!);
    expect(vabp.signalValue).toBeLessThan(cumene.Tb_K!);
    expect(smx.signalValue).toBeGreaterThan(0);
    expect(api.signalValue).toBeGreaterThan(0);
    expect(refractiveIndex.signalValue).toBeGreaterThan(1);
    expect(refractiveIndex.warnings).toEqual([]);
  });

  it('computes cetane-number analyzer output when component metadata exists', () => {
    const cetaneCompounds: CanopyCompound[] = [
      { ...benzene, name: 'PARAFFIN-A', displayName: 'Paraffin A', cetaneNumber: 55 },
      { ...cumene, name: 'PARAFFIN-B', displayName: 'Paraffin B', cetaneNumber: 35 },
    ];
    const stream = makeStream(cetaneCompounds, 15, [0.6, 0.4], 340, 101325);

    const ceta = solveAnalyzer(cetaneCompounds, stream, { propertyCode: 'CETA' });

    expect(ceta.signalValue).toBeCloseTo(47, 8);
    expect(ceta.warnings[0]).toMatch(/metadata/i);
  });

  it('computes extended petroleum-style analyzer outputs with explicit approximation warnings', () => {
    const compounds = [benzene, cumene, dipb];
    const stream = makeStream(compounds, 20, [0.45, 0.35, 0.2], 365, 150000);

    const d86 = solveAnalyzer(compounds, stream, { propertyCode: 'D86T' });
    const d1160 = solveAnalyzer(compounds, stream, { propertyCode: 'D1160T' });
    const d2887 = solveAnalyzer(compounds, stream, { propertyCode: 'D2887T' });
    const reidvp = solveAnalyzer(compounds, stream, { propertyCode: 'REIDVP' });
    const ri = solveAnalyzer(compounds, stream, { propertyCode: 'RI' });

    expect(d86.signalValue).toBeGreaterThan(benzene.Tb_K!);
    expect(d1160.signalValue).toBeGreaterThan(d86.signalValue);
    expect(d2887.signalValue).toBeGreaterThan(d1160.signalValue);
    expect(reidvp.signalValue).toBeGreaterThan(0);
    expect(ri.signalValue).toBeGreaterThan(1);
    expect(d86.warnings[0]).toMatch(/estimated/i);
    expect(d1160.warnings[0]).toMatch(/estimated/i);
    expect(d2887.warnings[0]).toMatch(/estimated/i);
    expect(reidvp.warnings[0]).toMatch(/approximated/i);
  });

  it('computes water saturation from hydrocarbon-water solubility data', () => {
    const compounds = [water, benzene];
    const stream = makeStream(compounds, 10, [0.02, 0.98], 298.15, 101325);

    const waterSaturation = solveAnalyzer(compounds, stream, { propertyCode: 'WATSAT' });

    expect(waterSaturation.signalValue).toBeGreaterThan(100);
    expect(waterSaturation.warnings[0]).toMatch(/solubility/i);
  });

  it('adjusts one ionic component to achieve charge balance', () => {
    const ionicCompounds: CanopyCompound[] = [
      { ...benzene, name: 'SOLVENT', displayName: 'Solvent', chargeNumber: 0 },
      { ...cumene, name: 'CAT+', displayName: 'Cation', chargeNumber: 1 },
      { ...dipb, name: 'AN-', displayName: 'Anion', chargeNumber: -1 },
    ];
    const inlet = makeStream(ionicCompounds, 10, [0.7, 0.2, 0.1], 330, 101325);

    const result = solveChargeBalance(ionicCompounds, inlet, {
      adjustComponentIndex: 2,
      targetCharge: 0,
    });

    expect(result.residualCharge).toBeCloseTo(0, 8);
    expect(result.adjustedComponentFlow_molps).toBeCloseTo(2, 8);
    expect(result.outlet.totalFlow_molps).toBeCloseTo(11, 8);
    expect(result.outlet.moleFractions?.reduce((sum, value) => sum + value, 0)).toBeCloseTo(1, 8);
  });

  it('adds analyzer and charge-balance blocks to the palette', () => {
    useCanopyStore.setState(state => ({ ...state, units: {}, streams: {}, nodePositions: {} }));
    const store = useCanopyStore.getState();
    store.addUnitFromPalette('Analyzer', { x: 100, y: 100 });
    store.addUnitFromPalette('ChargeBalance', { x: 140, y: 140 });

    const types = Object.values(useCanopyStore.getState().units).map(unit => unit.type);
    expect(types).toEqual(['Analyzer', 'ChargeBalance']);
  });

  it('splits dirty water into aqueous and hydrocarbon phases using solubility-guided carryover', () => {
    const compounds = [water, benzene];
    const feed = makeStream(compounds, 100, [0.4, 0.6], 298.15, 101325);

    const result = solveDecanter(compounds, feed, 298.15, 101325, {
      dirtyWaterMode: true,
      freeWaterSplitFraction: 0.99,
      waterCarryoverFraction: 0.002,
      hydrocarbonToAqueousFraction: 0.003,
    });

    expect(result.liquidI.totalFlow_molps).toBeGreaterThan(0);
    expect(result.liquidII.totalFlow_molps).toBeGreaterThan(0);
    expect(result.freeWaterSeparated_molps).toBeGreaterThan(30);
    expect(result.aqueousPhaseIndex).toBe(2);
    expect(result.waterSaturationFraction).toBeGreaterThan(1);
    expect(result.liquidI.moleFractions?.[1] ?? 0).toBeGreaterThan(0.9);
    expect(result.liquidII.moleFractions?.[0] ?? 0).toBeGreaterThan(0.9);
    expect((result.liquidI.totalFlow_molps ?? 0) + (result.liquidII.totalFlow_molps ?? 0)).toBeCloseTo(100, 8);
  });
});
