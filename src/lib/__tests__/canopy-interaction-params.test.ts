import { describe, expect, it } from 'vitest';

import {
  buildEosInteractionParamsFromAspenRows,
  buildInteractionParamsForPackage,
  buildNrtlParamsFromAspenRows,
  buildUniquacParamsFromAspenRows,
  buildWilsonParamsFromAspenRows,
  type AspenBinaryParamRow,
} from '../canopy/interaction-params';
import { liquidMolarVolume_m3pmol } from '../canopy/thermo';
import type { CanopyCompound } from '../canopy/types';

function makeCompound(
  name: string,
  overrides: Partial<CanopyCompound> = {},
): CanopyCompound {
  return {
    name,
    displayName: name,
    molecularWeight: 50,
    Tc_K: 500,
    Pc_Pa: 5_000_000,
    omega: 0.2,
    dipprCoeffs: {},
    antoine: null,
    ...overrides,
  };
}

const COMPOUNDS: CanopyCompound[] = [
  makeCompound('ACETONE', { uniquac_r: 2.57, uniquac_q: 2.34, Tc_K: 508.1, Pc_Pa: 4_700_000, omega: 0.307 }),
  makeCompound('WATER', { uniquac_r: 0.92, uniquac_q: 1.4, Tc_K: 647.1, Pc_Pa: 22_064_000, omega: 0.344 }),
];

describe('canopy Aspen row parsers', () => {
  it('maps NRTL Aspen rows into canopy raw extended slots', () => {
    const rows: AspenBinaryParamRow[] = [
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'aij', Value: 1.1 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'aji', Value: 2.2 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'bij', Value: 333.0 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'bji', Value: 444.0 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'cij', Value: 0.27 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'dij', Value: 0.015 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'eij', Value: 5.5 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'eji', Value: 6.6 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'fij', Value: 0.007 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'fji', Value: 0.008 },
    ];

    const params = buildNrtlParamsFromAspenRows(COMPOUNDS, rows);

    expect(params.A).toEqual([[0, 0], [0, 0]]);
    expect(params.B).toEqual([[0, 0], [0, 0]]);
    expect(params.a_ext).toEqual([[0, 1.1], [2.2, 0]]);
    expect(params.b_ext).toEqual([[0, 333.0], [444.0, 0]]);
    expect(params.e_ext).toEqual([[0, 5.5], [6.6, 0]]);
    expect(params.f_ext).toEqual([[0, 0.007], [0.008, 0]]);
    expect(params.alpha).toEqual([[0, 0.27], [0.27, 0]]);
    expect(params.d_alpha).toEqual([[0, 0.015], [0.015, 0]]);
  });

  it('maps UNIQUAC Aspen rows into canopy extended slots and preserves structural parameters', () => {
    const rows: AspenBinaryParamRow[] = [
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'aij', Value: -12.5 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'aji', Value: 8.25 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'bij', Value: 950.0 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'bji', Value: 725.0 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'cij', Value: 0.9 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'cji', Value: 1.1 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'dij', Value: 0.002 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'dji', Value: 0.003 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'eij', Value: 15.0 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'eji', Value: 18.0 },
    ];

    const params = buildUniquacParamsFromAspenRows(COMPOUNDS, rows);

    expect(params.a).toEqual([[0, 0], [0, 0]]);
    expect(params.b).toEqual([[0, 0], [0, 0]]);
    expect(params.r).toEqual([2.57, 0.92]);
    expect(params.q).toEqual([2.34, 1.4]);
    expect(params.a_ext).toEqual([[0, -12.5], [8.25, 0]]);
    expect(params.b_ext).toEqual([[0, 950.0], [725.0, 0]]);
    expect(params.c_ext).toEqual([[0, 0.9], [1.1, 0]]);
    expect(params.d_ext).toEqual([[0, 0.002], [0.003, 0]]);
    expect(params.e_ext).toEqual([[0, 15.0], [18.0, 0]]);
  });

  it('maps Wilson Aspen rows into canopy extended slots and derives liquid molar volumes locally', () => {
    const rows: AspenBinaryParamRow[] = [
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'aij', Value: -45.0 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'aji', Value: 22.0 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'bij', Value: 1100.0 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'bji', Value: 1400.0 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'cij', Value: 0.5 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'cji', Value: 0.8 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'dij', Value: 0.0012 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'dji', Value: 0.0015 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'eij', Value: 12.0 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'eji', Value: 14.0 },
    ];

    const params = buildWilsonParamsFromAspenRows(COMPOUNDS, rows);

    expect(params.A).toEqual([[0, 0], [0, 0]]);
    expect(params.a_ext).toEqual([[0, -45.0], [22.0, 0]]);
    expect(params.b_ext).toEqual([[0, 1100.0], [1400.0, 0]]);
    expect(params.c_ext).toEqual([[0, 0.5], [0.8, 0]]);
    expect(params.d_ext).toEqual([[0, 0.0012], [0.0015, 0]]);
    expect(params.e_ext).toEqual([[0, 12.0], [14.0, 0]]);
    expect(params.molarVolumes[0]).toBeCloseTo(liquidMolarVolume_m3pmol(COMPOUNDS[0], 298.15) * 1e6, 8);
    expect(params.molarVolumes[1]).toBeCloseTo(liquidMolarVolume_m3pmol(COMPOUNDS[1], 298.15) * 1e6, 8);
  });

  it('maps EOS rows into symmetric kij matrices with optional temperature terms', () => {
    const rows: AspenBinaryParamRow[] = [
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'kij', Value: 0.125 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'aij', Value: 1.2e-3 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'bij', Value: -2.5e-5 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'cij', Value: 7.5e-8 },
    ];

    const params = buildEosInteractionParamsFromAspenRows(COMPOUNDS, rows);

    expect(params.kij).toEqual([[0, 0.125], [0.125, 0]]);
    expect(params.kij_a).toEqual([[0, 1.2e-3], [1.2e-3, 0]]);
    expect(params.kij_b).toEqual([[0, -2.5e-5], [-2.5e-5, 0]]);
    expect(params.kij_c).toEqual([[0, 7.5e-8], [7.5e-8, 0]]);
  });

  it('dispatches package-specific parsing through a single canopy helper', () => {
    const rows: AspenBinaryParamRow[] = [
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'aij', Value: 3.2 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'bij', Value: 250.0 },
      { Compound1Name: 'ACETONE', Compound2Name: 'WATER', ElementLabel: 'cij', Value: 0.33 },
    ];

    const nrtl = buildInteractionParamsForPackage('NRTL', COMPOUNDS, rows);
    const eos = buildInteractionParamsForPackage('Peng-Robinson', COMPOUNDS, rows);
    const ideal = buildInteractionParamsForPackage('Ideal', COMPOUNDS, rows);

    expect(nrtl?.nrtl?.a_ext).toEqual([[0, 3.2], [0, 0]]);
    expect(nrtl?.nrtl?.b_ext).toEqual([[0, 250.0], [0, 0]]);
    expect(nrtl?.nrtl?.alpha).toEqual([[0, 0.33], [0.33, 0]]);
    expect(eos?.kij_a).toEqual([[0, 3.2], [3.2, 0]]);
    expect(eos?.kij_b).toEqual([[0, 250.0], [250.0, 0]]);
    expect(eos?.kij_c).toEqual([[0, 0.33], [0.33, 0]]);
    expect(ideal).toBeNull();
  });
});
