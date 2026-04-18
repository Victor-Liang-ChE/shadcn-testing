import { describe, expect, it } from 'vitest';

import {
  getEditableInteractionMatrices,
  getFluidPackageDiagnostics,
} from '../canopy/property-package';
import type { InteractionParams } from '../canopy/thermo';
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
    antoine: {
      C1: 10,
      C2: -1000,
      C3: 0,
      C4: 0,
      C5: 0,
      C6: 0,
      C7: 1,
      Tmin_K: 200,
      Tmax_K: 500,
    },
    ...overrides,
  };
}

const COMPOUNDS: CanopyCompound[] = [
  makeCompound('ETHANOL', { uniquac_r: 2.1, uniquac_q: 1.97, unifacGroups: { 1: 1, 2: 1, 15: 1 } }),
  makeCompound('WATER', { uniquac_r: 0.92, uniquac_q: 1.4, unifacGroups: { 16: 1 } }),
];

describe('canopy property-package helpers', () => {
  it('prefers Aspen-form NRTL matrices when extended coefficients are present', () => {
    const params: InteractionParams = {
      nrtl: {
        A: [[0, 0], [0, 0]],
        B: [[0, 0], [0, 0]],
        alpha: [[0, 0.3], [0.3, 0]],
        a_ext: [[0, 1], [2, 0]],
        b_ext: [[0, 100], [200, 0]],
        d_alpha: [[0, 0.01], [0.01, 0]],
        e_ext: [[0, 3], [4, 0]],
        f_ext: [[0, 5], [6, 0]],
      },
    };

    const matrices = getEditableInteractionMatrices('NRTL', params, 2);

    expect(matrices.map(m => m.label)).toEqual(['aij', 'bij', 'alpha (cij)', 'dij', 'eij', 'fij']);
    expect(matrices[0].field).toBe('nrtl.a_ext');
    expect(matrices[1].matrix[0][1]).toBe(100);
  });

  it('returns EOS temperature-dependent matrices when available', () => {
    const params: InteractionParams = {
      kij: [[0, 0.1], [0.1, 0]],
      kij_a: [[0, 1], [1, 0]],
      kij_b: [[0, 2], [2, 0]],
      kij_c: [[0, 3], [3, 0]],
    };

    const matrices = getEditableInteractionMatrices('Peng-Robinson', params, 2);

    expect(matrices.map(m => m.field)).toEqual(['kij', 'kij_a', 'kij_b', 'kij_c']);
  });

  it('flags missing UNIQUAC structural parameters as an error', () => {
    const compounds = [
      makeCompound('ETHANOL', { uniquac_r: 2.1, uniquac_q: 1.97 }),
      makeCompound('WATER', { uniquac_r: undefined, uniquac_q: undefined }),
    ];
    const diagnostics = getFluidPackageDiagnostics(compounds, 'UNIQUAC', {
      uniquac: {
        a: [[0, 0], [0, 0]],
        b: [[0, 0], [0, 0]],
        r: [2.1, 1],
        q: [1.97, 1],
      },
    });

    expect(diagnostics.some(d => d.level === 'error' && d.message.includes('r/q'))).toBe(true);
  });

  it('marks UNIFAC-DMD as ready when groups and interaction data are present', () => {
    const diagnostics = getFluidPackageDiagnostics(COMPOUNDS, 'UNIFAC-DMD', {
      unifacDmd: {
        compGroups: COMPOUNDS.map(c => c.unifacGroups ?? {}),
        data: {
          subgroups: [
            { subgroupId: 1, mainGroupId: 1, R: 0.9011, Q: 0.848 },
            { subgroupId: 16, mainGroupId: 7, R: 1.4, Q: 1.2 },
          ],
          interactions: {
            '1-7': { a: 10, b: 0.1, c: 0 },
            '7-1': { a: -5, b: 0.05, c: 0 },
          },
        },
      },
    });

    expect(diagnostics.some(d => d.level === 'ok' && d.message.includes('loaded'))).toBe(true);
  });

  it('warns when EOS is running with zero binary interaction parameters', () => {
    const diagnostics = getFluidPackageDiagnostics(COMPOUNDS, 'SRK', {});
    expect(diagnostics.some(d => d.level === 'warn' && d.message.includes('zero binary interaction'))).toBe(true);
  });

  it('diagnoses Electrolyte-NRTL readiness from solvent and charge metadata', () => {
    const diagnostics = getFluidPackageDiagnostics(
      [
        { ...COMPOUNDS[0], chargeNumber: 0, isElectrolyteSolvent: true, dielectricConst: { e1: 78.36, e2: 0, e3_K: 298.15 } },
        { ...COMPOUNDS[1], chargeNumber: 1 },
      ],
      'Electrolyte-NRTL',
      {
        elecNrtl: {
          nrtl: {
            A: [[0, 0], [0, 0]],
            B: [[0, 0], [0, 0]],
            alpha: [[0, 0.2], [0.2, 0]],
          },
          charges: [0, 1],
          epsilon_r: 78.36,
          solventMW: 46.07,
        },
      },
    );

    expect(diagnostics.some(d => d.level === 'ok' && d.message.includes('Electrolyte-NRTL'))).toBe(true);
    expect(diagnostics.some(d => d.level === 'error')).toBe(false);
  });
});
