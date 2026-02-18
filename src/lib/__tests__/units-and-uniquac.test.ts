import { describe, it, expect } from 'vitest';

import { calcPurePR, calcPureSRK, gammaUni } from '../residue-curves-ode';
import { calculateUniquacGammaMulticomponent } from '../vle-calculations-multicomponent';

// Minimal mock CompoundData shapes for tests
const makeMockComp = (overrides: any = {}) => ({ name: 'X', antoine: null, uniquacParams: null, prParams: undefined, srkParams: undefined, ...overrides } as any);

describe('Unit conversions & UNIQUAC parity', () => {
  it('calcPurePR handles Pc (kPa) + Tc (°C) inputs equivalently to Pa+K', () => {
    const pr_kpa = { Tc_K: 234.95, Pc_Pa: 4700 /* kPa-like value */, omega: 0.304 } as any; // source in °C/kPa
    const pr_pa = { Tc_K: 234.95 + 273.15, Pc_Pa: 4700 * 1000, omega: 0.304 } as any; // converted

    const r = 298.15;
    const a_kpa = calcPurePR(r, pr_kpa).a;
    const b_kpa = calcPurePR(r, pr_kpa).b;

    const a_pa = calcPurePR(r, pr_pa).a;
    const b_pa = calcPurePR(r, pr_pa).b;

    expect(a_kpa).toBeCloseTo(a_pa, 12);
    expect(b_kpa).toBeCloseTo(b_pa, 12);
  });

  it('calcPureSRK handles Pc (kPa) + Tc (°C) inputs equivalently to Pa+K', () => {
    const srk_kpa = { Tc_K: 507.6, Pc_Pa: 3000 /* kPa-like value */, omega: 0.344 } as any; // source in °C/kPa
    const srk_pa = { Tc_K: 507.6 + 273.15, Pc_Pa: 3000 * 1000, omega: 0.344 } as any; // converted

    const T = 400;
    const a_kpa = calcPureSRK(T, srk_kpa).a;
    const b_kpa = calcPureSRK(T, srk_kpa).b;

    const a_pa = calcPureSRK(T, srk_pa).a;
    const b_pa = calcPureSRK(T, srk_pa).b;

    expect(a_kpa).toBeCloseTo(a_pa, 12);
    expect(b_kpa).toBeCloseTo(b_pa, 12);
  });

  it('gammaUni (residue ODE) matches calculateUniquacGammaMulticomponent (library) for the same params', () => {
    const comps = [
      makeMockComp({ uniquacParams: { r: 2.0, q: 1.2 } }),
      makeMockComp({ uniquacParams: { r: 1.5, q: 1.0 } }),
      makeMockComp({ uniquacParams: { r: 1.0, q: 0.9 } }),
    ];

    const x = [0.3, 0.5, 0.2];
    const T = 350;

    // Ternary-style A/B parameters (simple symmetric example)
    const ternaryParams = {
      A01: 0.5, B01: 0,
      A10: 0.5, B10: 0,
      A02: 0.2, B02: 0,
      A20: 0.2, B20: 0,
      A12: 0.1, B12: 0,
      A21: 0.1, B21: 0,
    } as any;

    // Build the UniquacParameterMatrix expected by calculateUniquacGammaMulticomponent
    const paramsMap = new Map<string, any>();
    const setPair = (i: number, j: number, Aij: number, Bij: number) => paramsMap.set(`${i}-${j}`, { Aij, Aji: 0, Bij, Bji: 0 });
    setPair(0, 1, ternaryParams.A01, ternaryParams.B01);
    setPair(1, 0, ternaryParams.A10, ternaryParams.B10);
    setPair(0, 2, ternaryParams.A02, ternaryParams.B02);
    setPair(2, 0, ternaryParams.A20, ternaryParams.B20);
    setPair(1, 2, ternaryParams.A12, ternaryParams.B12);
    setPair(2, 1, ternaryParams.A21, ternaryParams.B21);

    const gamma_residue = gammaUni(x, T, comps as any, ternaryParams as any);
    const gamma_lib = calculateUniquacGammaMulticomponent(comps as any, x, T, paramsMap as any);

    expect(gamma_residue).not.toBeNull();
    expect(gamma_lib).not.toBeNull();

    for (let i = 0; i < 3; i++) {
      // allow small numerical differences between the two implementations
      expect((gamma_residue as number[])[i]).toBeCloseTo((gamma_lib as number[])[i], 2);
    }
  });
});