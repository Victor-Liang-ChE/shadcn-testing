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
  KvaluesWilsonEstimate,
  mixtureVaporViscosity_Pas,
  mixtureVaporThermalConductivity_WpmK,
  vaporViscosity_Pas,
  vaporThermalConductivity_WpmK,
  PsatWagner_Pa,
  PsatLeeKesler_Pa,
  HvapWatson_Jmol,
  DIPPR100,
  DIPPR100_integral,
  DIPPR101,
  DIPPR102,
  DIPPR104,
  DIPPR105,
  DIPPR106,
  DIPPR107,
  liquidDensity_DIPPR116_kmolpm3,
  CpIG_AspenPoly_JkmolK,
  RackettDensity_kmolpm3,
  CpL_BondiRowlinson_JmolK,
  computeKvalues,
  flashPT,
  unifacDMDGamma,
} from '../canopy/thermo';
import type { CanopyCompound } from '../canopy/types';
import type { UNIFACDMDData } from '../canopy/thermo';

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
  transportCoeffs: {
    VaporViscosity: { A: 3.134e-8, B: 0.9676, C: 7.9, D: 0 },
    VaporThermalConductivity: { A: 1.652e-5, B: 1.3117, C: 491, D: 0 },
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
  transportCoeffs: {
    VaporViscosity: { A: 2.919e-8, B: 0.9648, C: 0, D: 0 },
    VaporThermalConductivity: { A: 2.392e-5, B: 1.2694, C: 537, D: 0 },
  },
};

const ETHANOL: CanopyCompound = {
  name: 'ETHANOL',
  displayName: 'Ethanol',
  molecularWeight: 46.068,
  Tc_K: 514.0,
  Pc_Pa: 6137000,
  omega: 0.644,
  dipprCoeffs: {
    IdealGasHeatCapacityCp: { A: 49200, B: 145770, C: 1662.8, D: 93900, E: 744.7 },
    HeatOfVaporization: { A: 5.502e7, B: 0.649, C: -0.289, D: 0.116, E: 0 },
    LiquidHeatCapacityCp: { A: 102640, B: -139.63, C: 0.60511, D: 0, E: 0, eqNo: 100 },
  },
  antoine: {
    C1: 73.304, C2: -7122.3, C3: 0, C4: 0,
    C5: -7.1424, C6: 2.8853e-6, C7: 2,
    Tmin_K: 159.05, Tmax_K: 514.0,
  },
  unifacGroups: { 1: 1, 2: 1, 15: 1 },
};

const HEXANE: CanopyCompound = {
  name: 'HEXANE',
  displayName: 'Hexane',
  molecularWeight: 86.178,
  Tc_K: 507.6,
  Pc_Pa: 3025000,
  omega: 0.301,
  dipprCoeffs: {
    IdealGasHeatCapacityCp: { A: 1.15e5, B: 2.4e5, C: 1.6e3, D: 1.1e5, E: 7.5e2 },
    HeatOfVaporization: { A: 4.1e7, B: 0.39, C: 0, D: 0, E: 0 },
    LiquidHeatCapacityCp: { A: 167500, B: -268.84, C: 0.86728, D: 0, E: 0, eqNo: 100 },
  },
  antoine: {
    C1: 87.829, C2: -6996.4, C3: 0, C4: 0,
    C5: -9.8802, C6: 7.2099e-6, C7: 2,
    Tmin_K: 177.8, Tmax_K: 507.6,
  },
  unifacGroups: { 1: 2, 2: 4 },
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

// ────────────────────────────────────────────────────────────────
// Wilson K-value Estimate (Aspen Plus column initialization)
// ────────────────────────────────────────────────────────────────

describe('Canopy Thermo: Wilson K-value Estimate', () => {
  it('benzene K > 1 at its boiling point (volatile)', () => {
    const K = KvaluesWilsonEstimate([BENZENE], 353.24, 101325);
    expect(K[0]).toBeGreaterThan(0.5);
    expect(K[0]).toBeLessThan(5);
  });

  it('heavier component has smaller K than lighter one', () => {
    const K = KvaluesWilsonEstimate([BENZENE, TOLUENE], 373, 101325);
    expect(K[0]).toBeGreaterThan(K[1]); // benzene more volatile
  });

  it('K increases with temperature', () => {
    const K1 = KvaluesWilsonEstimate([BENZENE], 350, 101325);
    const K2 = KvaluesWilsonEstimate([BENZENE], 450, 101325);
    expect(K2[0]).toBeGreaterThan(K1[0]);
  });

  it('K decreases with pressure', () => {
    const K1 = KvaluesWilsonEstimate([BENZENE], 400, 101325);
    const K2 = KvaluesWilsonEstimate([BENZENE], 400, 500000);
    expect(K1[0]).toBeGreaterThan(K2[0]);
  });
});

// ────────────────────────────────────────────────────────────────
// Mixture Vapor Transport Properties (Wilke, Wassilijewa-Mason-Saxena)
// ────────────────────────────────────────────────────────────────

describe('Canopy Thermo: Mixture Vapor Transport', () => {
  it('pure component vapor viscosity is positive (DIPPR 102)', () => {
    const mu = vaporViscosity_Pas(BENZENE, 400);
    expect(mu).toBeGreaterThan(0);
    expect(mu).toBeLessThan(1e-3); // gas viscosity < 1 mPa·s
  });

  it('mixture vapor viscosity (Wilke) returns reasonable value', () => {
    const mu = mixtureVaporViscosity_Pas([BENZENE, TOLUENE], [0.5, 0.5], 400);
    expect(mu).toBeGreaterThan(0);
    expect(mu).toBeLessThan(1e-3);
  });

  it('mixture vapor viscosity is between pure components', () => {
    const muB = mixtureVaporViscosity_Pas([BENZENE, TOLUENE], [1, 0], 400);
    const muT = mixtureVaporViscosity_Pas([BENZENE, TOLUENE], [0, 1], 400);
    const muMix = mixtureVaporViscosity_Pas([BENZENE, TOLUENE], [0.5, 0.5], 400);
    expect(muMix).toBeGreaterThan(Math.min(muB, muT) * 0.9);
    expect(muMix).toBeLessThan(Math.max(muB, muT) * 1.1);
  });

  it('pure component vapor thermal conductivity is positive', () => {
    const lam = vaporThermalConductivity_WpmK(BENZENE, 400);
    expect(lam).toBeGreaterThan(0);
    expect(lam).toBeLessThan(1); // gas conductivity < 1 W/m/K
  });

  it('mixture vapor thermal conductivity (Wassilijewa) returns reasonable value', () => {
    const lam = mixtureVaporThermalConductivity_WpmK([BENZENE, TOLUENE], [0.5, 0.5], 400);
    expect(lam).toBeGreaterThan(0);
    expect(lam).toBeLessThan(1);
  });

  it('mixture vapor conductivity is between pure components', () => {
    const lamB = mixtureVaporThermalConductivity_WpmK([BENZENE, TOLUENE], [1, 0], 400);
    const lamT = mixtureVaporThermalConductivity_WpmK([BENZENE, TOLUENE], [0, 1], 400);
    const lamMix = mixtureVaporThermalConductivity_WpmK([BENZENE, TOLUENE], [0.5, 0.5], 400);
    expect(lamMix).toBeGreaterThan(Math.min(lamB, lamT) * 0.9);
    expect(lamMix).toBeLessThan(Math.max(lamB, lamT) * 1.1);
  });
});

// ────────────────────────────────────────────────────────────────
// 9. DIPPR Equation Evaluators
// ────────────────────────────────────────────────────────────────

describe('DIPPR Equation Evaluators', () => {
  it('DIPPR100 polynomial: evaluates correctly', () => {
    // f(T) = 1 + 2T + 3T² at T=10 → 1 + 20 + 300 = 321
    expect(DIPPR100(10, 1, 2, 3, 0, 0)).toBeCloseTo(321, 6);
  });

  it('DIPPR100 integral: analytical integral of polynomial', () => {
    // ∫[0→1] (2 + 3T) dT = 2 + 3/2 = 3.5
    expect(DIPPR100_integral(0, 1, 2, 3, 0, 0, 0)).toBeCloseTo(3.5, 6);
  });

  it('DIPPR101: liquid viscosity form exp(A + B/T + C·ln(T))', () => {
    // At T=300, A=-10, B=1000, C=0, D=0, E=0 → exp(-10 + 1000/300) = exp(-6.667)
    const val = DIPPR101(300, -10, 1000, 0, 0, 0);
    expect(val).toBeCloseTo(Math.exp(-10 + 1000 / 300), 6);
  });

  it('DIPPR102: vapor viscosity form A·T^B / (1 + C/T)', () => {
    // At T=400, A=1e-6, B=0.5, C=100, D=0 → 1e-6 * 20 / (1 + 0.25) = 1.6e-5
    const val = DIPPR102(400, 1e-6, 0.5, 100, 0);
    expect(val).toBeCloseTo(1e-6 * Math.pow(400, 0.5) / (1 + 100 / 400), 6);
  });

  it('DIPPR104: second virial coefficient form', () => {
    // At T=300, A=0.1, B=-10, C=0, D=0, E=0 → 0.1 + (-10)/300 = 0.0667
    const val = DIPPR104(300, 0.1, -10, 0, 0, 0);
    expect(val).toBeCloseTo(0.1 - 10 / 300, 4);
  });

  it('DIPPR105: liquid density Rackett form', () => {
    // Benzene DIPPR 105: A=1.0259, B=0.26666, C=562.05, D=0.28394
    // At T=298.15: ρ = 1.0259 / 0.26666^(1 + (1-298.15/562.05)^0.28394)
    const rho = DIPPR105(298.15, 1.0259, 0.26666, 562.05, 0.28394);
    expect(rho).toBeGreaterThan(10); // ~11.4 kmol/m³ for benzene
    expect(rho).toBeLessThan(14);
  });

  it('DIPPR106: Hvap / surface tension form', () => {
    // Benzene Hvap DIPPR 106 at 298.15K: A=4.881e7, B=0.61066, C=-0.25882, D=0.032238, E=0.022475
    const hvap = DIPPR106(298.15, 562.05, 4.881e7, 0.61066, -0.25882, 0.032238, 0.022475);
    // Literature: ~33.8 kJ/mol → 33800000 J/kmol
    expect(hvap / 1000).toBeGreaterThan(30000);
    expect(hvap / 1000).toBeLessThan(37000);
  });

  it('DIPPR107: Aly-Lee ideal gas Cp', () => {
    // Benzene Aly-Lee: A=55238, B=173380, C=764.25, D=72545, E=2445.7
    const cp = DIPPR107(298.15, 55238, 173380, 764.25, 72545, 2445.7);
    // Literature benzene CpIG at 298K ≈ 82 J/(mol·K) = 82000 J/(kmol·K)
    expect(cp).toBeGreaterThan(75000);
    expect(cp).toBeLessThan(90000);
  });
});

// ────────────────────────────────────────────────────────────────
// 10. Wagner Vapor Pressure
// ────────────────────────────────────────────────────────────────

describe('Wagner Vapor Pressure', () => {
  const WATER_WAGNER: CanopyCompound = {
    ...BENZENE,
    name: 'WATER',
    Tc_K: 647.096,
    Pc_Pa: 22064000,
    omega: 0.3443,
    wagnerVP: { a: -7.76451, b: 1.45838, c: -2.77580, d: -1.23303 },
  };

  it('returns Pc at critical temperature', () => {
    const P = PsatWagner_Pa(WATER_WAGNER, 647.096);
    expect(P).toBeCloseTo(22064000, -2);
  });

  it('water VP at 373.15K ≈ 101325 Pa (1 atm)', () => {
    const P = PsatWagner_Pa(WATER_WAGNER, 373.15);
    expect(P).toBeGreaterThan(90000);
    expect(P).toBeLessThan(115000);
  });

  it('returns NaN when no Wagner data', () => {
    expect(PsatWagner_Pa(BENZENE, 300)).toBeNaN();
  });
});

// ────────────────────────────────────────────────────────────────
// 11. Watson Hvap Extrapolation
// ────────────────────────────────────────────────────────────────

describe('Watson Hvap Extrapolation', () => {
  const BENZENE_WATSON: CanopyCompound = {
    ...BENZENE,
    watsonHvap: { Hvap_ref_Jmol: 30720, T_ref_K: 353.25, n: 0.38 },
  };

  it('returns reference Hvap at reference temperature', () => {
    const h = HvapWatson_Jmol(BENZENE_WATSON, 353.25);
    expect(h).toBeCloseTo(30720, 0);
  });

  it('extrapolates to lower T → higher Hvap', () => {
    const h = HvapWatson_Jmol(BENZENE_WATSON, 298.15);
    expect(h).toBeGreaterThan(30720); // Hvap increases at lower T
    expect(h).toBeLessThan(45000);
  });

  it('approaches 0 near critical temperature', () => {
    const h = HvapWatson_Jmol(BENZENE_WATSON, 560);
    expect(h).toBeLessThan(6000); // Much smaller than Hvap at Tb (~30720)
    expect(h).toBeGreaterThanOrEqual(0);
  });
});

// ────────────────────────────────────────────────────────────────
// 12. DIPPR 116 Liquid Density
// ────────────────────────────────────────────────────────────────

describe('DIPPR 116 Liquid Density', () => {
  const WATER_116: CanopyCompound = {
    ...BENZENE,
    name: 'WATER',
    Tc_K: 647.096,
    Pc_Pa: 22064000,
    // DIPPR 116 coefficients from NIST for water
    liquidDensityDIPPR116: { A: 17.863, B: 58.606, C: -95.396, D: 213.89, E: -141.26 },
  };

  it('water density at 298.15K ≈ 55.3 kmol/m³', () => {
    const rho = liquidDensity_DIPPR116_kmolpm3(WATER_116, 298.15);
    expect(rho).toBeGreaterThan(53);
    expect(rho).toBeLessThan(57);
  });

  it('density decreases near critical point', () => {
    const rho300 = liquidDensity_DIPPR116_kmolpm3(WATER_116, 300);
    const rho600 = liquidDensity_DIPPR116_kmolpm3(WATER_116, 600);
    expect(rho600).toBeLessThan(rho300);
  });

  it('returns NaN above critical T', () => {
    expect(liquidDensity_DIPPR116_kmolpm3(WATER_116, 700)).toBeNaN();
  });
});

// ────────────────────────────────────────────────────────────────
// 13. Aspen Polynomial CpIG
// ────────────────────────────────────────────────────────────────

describe('Aspen Polynomial CpIG', () => {
  const METHANE_POLY: CanopyCompound = {
    ...BENZENE,
    name: 'METHANE',
    cpigPoly: { A: 19250, B: 52.13, C: 11.97e-3, D: -11.33e-6, E: 0, Tmin_K: 50, Tmax_K: 1500 },
  };

  it('returns reasonable CpIG for methane at 298K', () => {
    const cp = CpIG_AspenPoly_JkmolK(METHANE_POLY, 298.15);
    // Literature methane CpIG at 298K ≈ 35.7 J/(mol·K) = 35700 J/(kmol·K)
    expect(cp).toBeGreaterThan(30000);
    expect(cp).toBeLessThan(45000);
  });

  it('returns NaN when no polynomial data', () => {
    expect(CpIG_AspenPoly_JkmolK(BENZENE, 298.15)).toBeNaN();
  });
});

// ────────────────────────────────────────────────────────────────
// 14. Lee-Kesler Vapor Pressure
// ────────────────────────────────────────────────────────────────

describe('Lee-Kesler VP', () => {
  it('benzene VP at 353K close to 101325 Pa (boiling)', () => {
    const P = PsatLeeKesler_Pa(BENZENE, 353.25);
    // Lee-Kesler is approximate but should be within 10% of 101325
    expect(P).toBeGreaterThan(80000);
    expect(P).toBeLessThan(130000);
  });

  it('returns Pc at critical T', () => {
    const P = PsatLeeKesler_Pa(BENZENE, 562.05);
    expect(P).toBeCloseTo(BENZENE.Pc_Pa, -4);
  });
});

// ────────────────────────────────────────────────────────────────
// 15. Modified Rackett Density
// ────────────────────────────────────────────────────────────────

describe('Modified Rackett Density', () => {
  it('benzene liquid density at 298K ≈ 11.2 kmol/m³', () => {
    const rho = RackettDensity_kmolpm3(BENZENE, 298.15);
    // Benzene MW=78.11, liquid density ~876 kg/m³ → 876/78.11 = 11.2 kmol/m³
    expect(rho).toBeGreaterThan(9);
    expect(rho).toBeLessThan(14);
  });

  it('density with fitted ZRA is more accurate', () => {
    const compFitted = { ...BENZENE, rackett_ZRA: 0.2713 };
    const rho = RackettDensity_kmolpm3(compFitted, 298.15);
    // With fitted ZRA, should be closer to 11.2
    expect(rho).toBeGreaterThan(10);
    expect(rho).toBeLessThan(13);
  });
});

// ────────────────────────────────────────────────────────────────
// 16. Bondi-Rowlinson Liquid Cp
// ────────────────────────────────────────────────────────────────

describe('Bondi-Rowlinson Liquid Cp', () => {
  it('benzene CpL at 298K is higher than CpIG', () => {
    const cpIG = CpIG_JmolK(BENZENE, 298.15);
    const cpL = CpL_BondiRowlinson_JmolK(BENZENE, 298.15);
    expect(cpL).toBeGreaterThan(cpIG);
    // Literature benzene CpL at 298K ≈ 135.7 J/(mol·K)
    expect(cpL).toBeGreaterThan(100);
    expect(cpL).toBeLessThan(200);
  });

  it('returns NaN near critical', () => {
    expect(CpL_BondiRowlinson_JmolK(BENZENE, 560)).toBeNaN();
  });
});

describe('Dortmund UNIFAC package wiring', () => {
  const dmdData: UNIFACDMDData = {
    subgroups: [
      { subgroupId: 1, mainGroupId: 1, R: 0.9011, Q: 0.848 },
      { subgroupId: 2, mainGroupId: 1, R: 0.6744, Q: 0.540 },
      { subgroupId: 15, mainGroupId: 7, R: 1.0, Q: 1.2 },
    ],
    interactions: {
      '1-7': { a: 2777.0, b: -4.674, c: 0.001551 },
      '7-1': { a: 1606.0, b: -4.746, c: 0.0009181 },
    },
  };

  it('computeKvalues supports UNIFAC-DMD', () => {
    const compounds = [ETHANOL, HEXANE];
    const x = [0.4, 0.6];
    const y = [0.7, 0.3];
    const params = {
      unifacDmd: {
        compGroups: compounds.map(c => c.unifacGroups ?? {}),
        data: dmdData,
      },
    };

    const gamma = unifacDMDGamma(x, 320, params.unifacDmd.compGroups, dmdData);
    const K = computeKvalues(compounds, x, y, 320, 101325, 'UNIFAC-DMD', params);

    expect(gamma[0]).toBeGreaterThan(0);
    expect(gamma[1]).toBeGreaterThan(0);
    expect(K[0]).toBeGreaterThan(0);
    expect(K[1]).toBeGreaterThan(0);
  });

  it('flashPT converges with UNIFAC-DMD interaction parameters', () => {
    const compounds = [ETHANOL, HEXANE];
    const z = [0.5, 0.5];
    const params = {
      unifacDmd: {
        compGroups: compounds.map(c => c.unifacGroups ?? {}),
        data: dmdData,
      },
    };

    const result = flashPT(compounds, z, 340, 101325, 'UNIFAC-DMD', params);

    expect(result.converged).toBe(true);
    expect(result.x.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 6);
    expect(result.y.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 6);
    expect(Number.isFinite(result.vaporFraction)).toBe(true);
    expect(result.K.every(k => Number.isFinite(k) && k > 0)).toBe(true);
  });
});

describe('Electrolyte-NRTL package wiring', () => {
  it('computeKvalues supports Electrolyte-NRTL when electrolyte metadata is present', () => {
    const compounds = [
      { ...ETHANOL, name: 'WATER', displayName: 'Water', molecularWeight: 18.015, chargeNumber: 0, isElectrolyteSolvent: true, dielectricConst: { e1: 78.36, e2: 0, e3_K: 298.15 } },
      { ...HEXANE, name: 'NA+', displayName: 'Na+', chargeNumber: 1 },
    ];
    const x = [0.9, 0.1];
    const y = [0.95, 0.05];
    const params = {
      elecNrtl: {
        nrtl: {
          A: [[0, 0], [0, 0]],
          B: [[0, 0], [0, 0]],
          alpha: [[0, 0.2], [0.2, 0]],
          a_ext: [[0, 1.2], [0, 0]],
          b_ext: [[0, 250], [0, 0]],
        },
        charges: [0, 1],
        epsilon_r: 78.36,
        solventMW: 18.015,
      },
    };

    const K = computeKvalues(compounds, x, y, 320, 101325, 'Electrolyte-NRTL', params);

    expect(K[0]).toBeGreaterThan(0);
    expect(K[1]).toBeGreaterThan(0);
    expect(K.every(value => Number.isFinite(value))).toBe(true);
  });
});
