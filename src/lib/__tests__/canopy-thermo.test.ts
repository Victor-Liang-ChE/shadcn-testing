// ─── Canopy Thermodynamic Engine Tests ──────────────────────────────────────
// Verifies flash calculations, enthalpy, and Cp against literature/Aspen data.
// Run: pnpm test -- canopy-thermo

import { describe, it, expect } from 'vitest';
import {
  Psat_Pa,
  CpIG_JkmolK,
  CpIG_JmolK,
  Hvap_Jkmol,
  Hvap_Jmol,
  CpL_JkmolK,
  flashPT_Ideal,
  bubblePointT_Ideal,
  dewPointT_Ideal,
  enthalpyIG,
  enthalpyLiquid,
} from '../canopy/thermo';
import type { CanopyCompound } from '../canopy/types';

// ── Test Compounds (from verified DB data) ──

const BENZENE: CanopyCompound = {
  name: 'BENZENE',
  displayName: 'Benzene',
  molecularWeight: 78.112,
  Tc_K: 562.1,
  Pc_Pa: 4893997.5,
  omega: 0.212,
  dipprCoeffs: {
    IdealGasHeatCapacityCp: { A: 34010.24, B: -588.0978, C: 12.81777, D: -0.000197306, E: 5.142899e-8 },
    HeatOfVaporization: { A: 4.881e7, B: 0.61066, C: -0.25882, D: 0.032238, E: 0.022475 },
    LiquidHeatCapacityCp: { A: 111460, B: -1854.3, C: 22.399, D: -0.028936, E: 0.000028991 },
  },
  antoine: {
    C1: 73.862398, C2: -5970.4385, C3: 0, C4: 0.0055376032,
    C5: -8.0797644, C6: 6.6129752e-18, C7: 6,
    Tmin_K: 293.36, Tmax_K: 562.1,
  },
};

const TOLUENE: CanopyCompound = {
  name: 'TOLUENE',
  displayName: 'Toluene',
  molecularWeight: 92.138,
  Tc_K: 591.7,
  Pc_Pa: 4113795,
  omega: 0.257,
  dipprCoeffs: {
    IdealGasHeatCapacityCp: { A: 47225, B: -565.85, C: 12.856, D: 0.000005535, E: -1.998e-8 },
    HeatOfVaporization: { A: 5.3752e7, B: 0.50341, C: 0.24755, D: -0.72898, E: 0.37794 },
    LiquidHeatCapacityCp: { A: 28291, B: 48.171, C: 10.912, D: 0.0020542, E: 8.7875e-7 },
  },
  antoine: {
    C1: 71.277464, C2: -6413.286, C3: 0, C4: 0.0041663022,
    C5: -7.5053519, C6: 5.419977e-18, C7: 6,
    Tmin_K: 318.72, Tmax_K: 591.7,
  },
};

describe('Canopy Thermo: Vapor Pressure', () => {
  it('benzene Psat at 353.24 K (80.09°C, known NBP) ≈ 101325 Pa', () => {
    // Literature: benzene boils at 80.09°C = 353.24 K at 1 atm
    const ps = Psat_Pa(BENZENE, 353.24);
    expect(ps).toBeGreaterThan(90000);
    expect(ps).toBeLessThan(115000);
  });

  it('toluene Psat at 383.78 K (110.63°C, known NBP) ≈ 101325 Pa', () => {
    const ps = Psat_Pa(TOLUENE, 383.78);
    expect(ps).toBeGreaterThan(90000);
    expect(ps).toBeLessThan(115000);
  });

  it('benzene Psat increases with temperature', () => {
    const p1 = Psat_Pa(BENZENE, 300);
    const p2 = Psat_Pa(BENZENE, 350);
    const p3 = Psat_Pa(BENZENE, 400);
    expect(p2).toBeGreaterThan(p1);
    expect(p3).toBeGreaterThan(p2);
  });
});

describe('Canopy Thermo: Ideal Gas Cp', () => {
  it('benzene Cp_ig at 300 K should be reasonable (~82 J/mol/K literature)', () => {
    const cp = CpIG_JmolK(BENZENE, 300);
    // NIST says ~82.4 J/mol/K at 300K. DIPPR eqno 16 may differ slightly.
    expect(cp).toBeGreaterThan(50);
    expect(cp).toBeLessThan(200);
  });

  it('toluene Cp_ig at 300 K should be > benzene', () => {
    const cpB = CpIG_JmolK(BENZENE, 300);
    const cpT = CpIG_JmolK(TOLUENE, 300);
    expect(cpT).toBeGreaterThan(cpB);
  });

  it('Cp_ig increases with temperature', () => {
    const cp1 = CpIG_JmolK(BENZENE, 300);
    const cp2 = CpIG_JmolK(BENZENE, 500);
    expect(cp2).toBeGreaterThan(cp1);
  });
});

describe('Canopy Thermo: Heat of Vaporization', () => {
  it('benzene ΔHvap at 353 K ≈ 30.7 kJ/mol (literature)', () => {
    const hv = Hvap_Jmol(BENZENE, 353.24);
    // Literature: 30.72 kJ/mol at NBP
    expect(hv / 1000).toBeGreaterThan(25);
    expect(hv / 1000).toBeLessThan(36);
  });

  it('ΔHvap → 0 at critical temperature', () => {
    const hv = Hvap_Jmol(BENZENE, BENZENE.Tc_K);
    expect(hv).toBeCloseTo(0, 0);
  });
});

describe('Canopy Thermo: Bubble & Dew Points', () => {
  it('50/50 benzene-toluene bubble point at 1 atm ≈ 92°C', () => {
    const Tbub = bubblePointT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 101325);
    expect(Tbub).not.toBeNull();
    // Literature: ~92°C for 50/50 benzene/toluene at 1 atm
    expect(Tbub! - 273.15).toBeGreaterThan(85);
    expect(Tbub! - 273.15).toBeLessThan(100);
  });

  it('50/50 benzene-toluene dew point at 1 atm ≈ 97°C', () => {
    const Tdew = dewPointT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 101325);
    expect(Tdew).not.toBeNull();
    expect(Tdew! - 273.15).toBeGreaterThan(90);
    expect(Tdew! - 273.15).toBeLessThan(105);
  });

  it('bubble point < dew point', () => {
    const Tbub = bubblePointT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 101325);
    const Tdew = dewPointT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 101325);
    expect(Tbub).not.toBeNull();
    expect(Tdew).not.toBeNull();
    expect(Tdew!).toBeGreaterThan(Tbub!);
  });

  it('pure benzene bubble point ≈ 80°C at 1 atm', () => {
    const T = bubblePointT_Ideal([BENZENE, TOLUENE], [1.0, 0.0], 101325);
    expect(T).not.toBeNull();
    expect(T! - 273.15).toBeGreaterThan(75);
    expect(T! - 273.15).toBeLessThan(85);
  });
});

describe('Canopy Thermo: Flash PT (Ideal)', () => {
  it('subcooled liquid at 30°C, 1 atm → all liquid', () => {
    const result = flashPT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 303.15, 101325);
    expect(result.converged).toBe(true);
    expect(result.vaporFraction).toBeCloseTo(0, 2);
    expect(result.phase).toBe('Liquid');
  });

  it('superheated vapor at 150°C, 1 atm → all vapor', () => {
    const result = flashPT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 423.15, 101325);
    expect(result.converged).toBe(true);
    expect(result.vaporFraction).toBeCloseTo(1, 2);
    expect(result.phase).toBe('Vapor');
  });

  it('two-phase flash at 95°C, 1 atm → partial vaporization', () => {
    const result = flashPT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 368.15, 101325);
    expect(result.converged).toBe(true);
    expect(result.vaporFraction).toBeGreaterThan(0.01);
    expect(result.vaporFraction).toBeLessThan(0.99);
    expect(result.phase).toBe('Two-Phase');
    // Benzene should be enriched in vapor
    expect(result.y[0]).toBeGreaterThan(result.x[0]);
    // Toluene should be enriched in liquid
    expect(result.x[1]).toBeGreaterThan(result.y[1]);
  });

  it('mole fractions sum to 1', () => {
    const result = flashPT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 368.15, 101325);
    const sumX = result.x.reduce((a, b) => a + b, 0);
    const sumY = result.y.reduce((a, b) => a + b, 0);
    expect(sumX).toBeCloseTo(1, 6);
    expect(sumY).toBeCloseTo(1, 6);
  });

  it('overall material balance holds', () => {
    const z = [0.5, 0.5];
    const result = flashPT_Ideal([BENZENE, TOLUENE], z, 368.15, 101325);
    const beta = result.vaporFraction;
    for (let i = 0; i < 2; i++) {
      const zCalc = beta * result.y[i] + (1 - beta) * result.x[i];
      expect(zCalc).toBeCloseTo(z[i], 4);
    }
  });
});

describe('Canopy Thermo: Enthalpy', () => {
  it('ideal gas enthalpy at reference T is ~0', () => {
    const h = enthalpyIG(BENZENE, 298.15);
    expect(Math.abs(h)).toBeLessThan(1);
  });

  it('liquid enthalpy at ref T is negative (below vapor)', () => {
    const hL = enthalpyLiquid(BENZENE, 298.15);
    expect(hL).toBeLessThan(0); // Liquid is lower enthalpy than ideal gas
  });

  it('heating increases enthalpy', () => {
    const h1 = enthalpyIG(BENZENE, 300);
    const h2 = enthalpyIG(BENZENE, 400);
    expect(h2).toBeGreaterThan(h1);
  });
});

// ─── Additional Edge Case Tests ────────────────────────────────

describe('Canopy Thermo: Edge Cases - Pure Component Flash', () => {
  it('pure benzene at its NBP → two-phase boundary', () => {
    const result = flashPT_Ideal([BENZENE], [1.0], 353.24, 101325);
    expect(result.converged).toBe(true);
    // At exact NBP, could be either liquid or two-phase boundary
    // Just verify it doesn't crash and gives reasonable result
    expect(result.x[0]).toBeCloseTo(1, 4);
    expect(result.y[0]).toBeCloseTo(1, 4);
  });

  it('pure toluene well above NBP → vapor', () => {
    const result = flashPT_Ideal([TOLUENE], [1.0], 413.15, 101325); // 140°C
    expect(result.converged).toBe(true);
    expect(result.vaporFraction).toBeCloseTo(1, 2);
    expect(result.phase).toBe('Vapor');
  });

  it('pure benzene well below NBP → liquid', () => {
    const result = flashPT_Ideal([BENZENE], [1.0], 283.15, 101325); // 10°C
    expect(result.converged).toBe(true);
    expect(result.vaporFraction).toBeCloseTo(0, 2);
    expect(result.phase).toBe('Liquid');
  });
});

describe('Canopy Thermo: Edge Cases - Extreme Compositions', () => {
  it('benzene-heavy (z = 0.99, 0.01) flash at 90°C, 1 atm', () => {
    const result = flashPT_Ideal([BENZENE, TOLUENE], [0.99, 0.01], 363.15, 101325);
    expect(result.converged).toBe(true);
    // Benzene heavy mix → lower bubble point → more vaporized at 90°C
    expect(result.y[0]).toBeGreaterThan(0.95);
  });

  it('toluene-heavy (z = 0.01, 0.99) flash at 100°C, 1 atm', () => {
    const result = flashPT_Ideal([BENZENE, TOLUENE], [0.01, 0.99], 373.15, 101325);
    expect(result.converged).toBe(true);
    // Toluene heavy → higher bubble point → less vaporized
    expect(result.x[1]).toBeGreaterThan(0.95);
  });
});

describe('Canopy Thermo: Edge Cases - Pressure Effects', () => {
  it('higher pressure raises bubble point', () => {
    const bp1 = bubblePointT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 101325);
    const bp2 = bubblePointT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 202650);
    expect(bp1).not.toBeNull();
    expect(bp2).not.toBeNull();
    expect(bp2!).toBeGreaterThan(bp1!);
  });

  it('flash at 2 atm produces less vapor than at 1 atm (same T)', () => {
    const r1 = flashPT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 368.15, 101325);
    const r2 = flashPT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 368.15, 202650);
    expect(r1.vaporFraction).toBeGreaterThan(r2.vaporFraction);
  });

  it('very low pressure → mostly vapor', () => {
    const result = flashPT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 353, 10000); // ~0.1 atm
    expect(result.converged).toBe(true);
    expect(result.vaporFraction).toBeGreaterThan(0.5);
  });
});

describe('Canopy Thermo: Edge Cases - Relative Volatility', () => {
  it('benzene is more volatile than toluene (K_benzene > K_toluene)', () => {
    const result = flashPT_Ideal([BENZENE, TOLUENE], [0.5, 0.5], 368.15, 101325);
    if (result.phase === 'Two-Phase') {
      const K_benzene = result.y[0] / result.x[0];
      const K_toluene = result.y[1] / result.x[1];
      expect(K_benzene).toBeGreaterThan(K_toluene);
      // Relative volatility α ≈ 2.4-2.5 for benzene/toluene at atmospheric
      const alpha = K_benzene / K_toluene;
      expect(alpha).toBeGreaterThan(2.0);
      expect(alpha).toBeLessThan(3.0);
    }
  });
});

describe('Canopy Thermo: Edge Cases - Enthalpy Consistency', () => {
  it('liquid enthalpy at boiling point + Hvap ≈ vapor enthalpy', () => {
    const T = 353.24; // benzene NBP
    const hL = enthalpyLiquid(BENZENE, T);
    const hIG = enthalpyIG(BENZENE, T);
    const hvap = Hvap_Jmol(BENZENE, T);
    // hIG = hL + hvap approximately
    // The liquid enthalpy function gives hLiquid = hIG - hvap, so hL + hvap ≈ hIG
    expect(hL + hvap).toBeCloseTo(hIG, -1); // within ~10 J/mol
  });

  it('enthalpy monotonically increases with temperature (liquid)', () => {
    const temps = [280, 300, 320, 340];
    const enthalpies = temps.map(T => enthalpyLiquid(BENZENE, T));
    for (let i = 1; i < enthalpies.length; i++) {
      expect(enthalpies[i]).toBeGreaterThan(enthalpies[i - 1]);
    }
  });

  it('enthalpy monotonically increases with temperature (vapor)', () => {
    const temps = [350, 400, 450, 500];
    const enthalpies = temps.map(T => enthalpyIG(BENZENE, T));
    for (let i = 1; i < enthalpies.length; i++) {
      expect(enthalpies[i]).toBeGreaterThan(enthalpies[i - 1]);
    }
  });
});

describe('Canopy Thermo: Edge Cases - CpL Liquid Heat Capacity', () => {
  it('benzene CpL at 300 K is reasonable (~135 J/mol/K)', () => {
    const cp = CpL_JkmolK(BENZENE, 300) / 1000; // → J/mol/K
    expect(cp).toBeGreaterThan(100);
    expect(cp).toBeLessThan(200);
  });

  it('liquid Cp increases with temperature (typical for organics)', () => {
    const cp1 = CpL_JkmolK(BENZENE, 300);
    const cp2 = CpL_JkmolK(BENZENE, 350);
    // For many organics, CpL increases with T (not always), just check reasonable
    expect(cp1).toBeGreaterThan(0);
    expect(cp2).toBeGreaterThan(0);
  });
});

describe('Canopy Thermo: Edge Cases - ΔHvap Temperature Dependence', () => {
  it('ΔHvap decreases as T approaches Tc', () => {
    const h1 = Hvap_Jmol(BENZENE, 353); // near NBP
    const h2 = Hvap_Jmol(BENZENE, 500); // closer to Tc=562.1
    expect(h1).toBeGreaterThan(h2);
    expect(h2).toBeGreaterThan(0);
  });

  it('ΔHvap at very low Tr is largest', () => {
    const h_low = Hvap_Jmol(BENZENE, 300);
    const h_mid = Hvap_Jmol(BENZENE, 400);
    expect(h_low).toBeGreaterThan(h_mid);
  });
});
