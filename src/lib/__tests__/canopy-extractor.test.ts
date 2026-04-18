import { describe, expect, it } from 'vitest';
import type { CanopyCompound, MaterialStream } from '../canopy/types';
import type { NRTLParams } from '../canopy/thermo';
import { solveExtractor } from '../canopy/solver';

const SOLUTE: CanopyCompound = {
  name: 'SOLUTE',
  displayName: 'Solute',
  molecularWeight: 100,
  Tc_K: 700,
  Pc_Pa: 4.5e6,
  omega: 0.2,
  dipprCoeffs: {},
  antoine: null,
};

const FEED_CARRIER: CanopyCompound = {
  name: 'FEED-CARRIER',
  displayName: 'Feed Carrier',
  molecularWeight: 90,
  Tc_K: 650,
  Pc_Pa: 4.2e6,
  omega: 0.18,
  dipprCoeffs: {},
  antoine: null,
};

const SOLVENT: CanopyCompound = {
  name: 'SOLVENT',
  displayName: 'Solvent',
  molecularWeight: 110,
  Tc_K: 760,
  Pc_Pa: 4.8e6,
  omega: 0.24,
  dipprCoeffs: {},
  antoine: null,
};

const compounds = [SOLUTE, FEED_CARRIER, SOLVENT];

const nrtl: NRTLParams = {
  A: [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
  ],
  B: [
    [0, 0.9, 0.05],
    [0.9, 0, 1.4],
    [0.05, 1.4, 0],
  ],
  alpha: [
    [0, 0.3, 0.3],
    [0.3, 0, 0.3],
    [0.3, 0.3, 0],
  ],
};

function makeStream(name: string, totalFlow_molps: number, moleFractions: number[]): MaterialStream {
  return {
    id: name,
    name,
    T_K: 298.15,
    P_Pa: 101325,
    totalFlow_molps,
    moleFractions,
    phase: 'Liquid',
    vaporFraction: 0,
    H_Jpmol: 0,
    x_liquid: moleFractions,
    y_vapor: moleFractions,
    sourceUnitId: null,
    sourcePortId: null,
    targetUnitId: null,
    targetPortId: null,
    solved: true,
  };
}

describe('solveExtractor', () => {
  it('conserves total flow and extracts more solute with more stages', () => {
    const feed = makeStream('feed', 10, [0.1, 0.9, 0]);
    const solvent = makeStream('solvent', 8, [0, 0, 1]);

    const oneStage = solveExtractor(compounds, feed, solvent, 1, 101325, {
      fluidPackage: 'NRTL',
      interactionParams: { nrtl },
    }, 2);

    const sixStage = solveExtractor(compounds, feed, solvent, 6, 101325, {
      fluidPackage: 'NRTL',
      interactionParams: { nrtl },
    }, 2);

    expect(oneStage.extractOutlet.totalFlow_molps + oneStage.raffinateOutlet.totalFlow_molps).toBeCloseTo(18, 8);
    expect(sixStage.extractOutlet.totalFlow_molps + sixStage.raffinateOutlet.totalFlow_molps).toBeCloseTo(18, 8);

    const soluteExtractOne = oneStage.extractOutlet.totalFlow_molps * oneStage.extractOutlet.moleFractions[0];
    const soluteExtractSix = sixStage.extractOutlet.totalFlow_molps * sixStage.extractOutlet.moleFractions[0];
    const soluteRaffSix = sixStage.raffinateOutlet.totalFlow_molps * sixStage.raffinateOutlet.moleFractions[0];

    expect(soluteExtractSix).toBeGreaterThan(soluteExtractOne);
    expect(soluteRaffSix).toBeLessThan(0.5);
    expect(sixStage.extractOutlet.moleFractions[2]).toBeGreaterThan(sixStage.raffinateOutlet.moleFractions[2]);
  });

  it('reduces to zero transfer under Ideal thermodynamics', () => {
    const feed = makeStream('feed', 10, [0.1, 0.9, 0]);
    const solvent = makeStream('solvent', 8, [0, 0, 1]);

    const ideal = solveExtractor(compounds, feed, solvent, 6, 101325, { fluidPackage: 'Ideal' }, 2);

    const soluteExtract = ideal.extractOutlet.totalFlow_molps * ideal.extractOutlet.moleFractions[0];
    const soluteRaffinate = ideal.raffinateOutlet.totalFlow_molps * ideal.raffinateOutlet.moleFractions[0];

    expect(soluteExtract).toBeCloseTo(0, 8);
    expect(soluteRaffinate).toBeCloseTo(1, 8);
  });
});
