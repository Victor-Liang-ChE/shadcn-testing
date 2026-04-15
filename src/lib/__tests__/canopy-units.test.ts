// ─── Canopy Unit Conversion Tests ────────────────────────────────────────────
import { describe, it, expect } from 'vitest';
import {
  fromSI_T, toSI_T, fromSI_P, toSI_P, fromSI_F, toSI_F, fromSI_Q,
  tempLabel, pressLabel, flowLabel, energyLabel,
} from '../canopy/units';

describe('Canopy Units: Temperature', () => {
  it('K → °C: 373.15 K = 100 °C', () => {
    expect(fromSI_T(373.15, 'C')).toBeCloseTo(100, 2);
  });
  it('K → K: identity', () => {
    expect(fromSI_T(300, 'K')).toBeCloseTo(300, 2);
  });
  it('K → °F: 373.15 K = 212 °F', () => {
    expect(fromSI_T(373.15, 'F')).toBeCloseTo(212, 1);
  });
  it('°C → K roundtrip', () => {
    expect(toSI_T(100, 'C')).toBeCloseTo(373.15, 2);
    expect(fromSI_T(toSI_T(25, 'C'), 'C')).toBeCloseTo(25, 5);
  });
  it('°F → K roundtrip', () => {
    expect(toSI_T(32, 'F')).toBeCloseTo(273.15, 2);
    expect(fromSI_T(toSI_T(72, 'F'), 'F')).toBeCloseTo(72, 5);
  });
  it('absolute zero: 0 K = -273.15 °C', () => {
    expect(fromSI_T(0, 'C')).toBeCloseTo(-273.15, 2);
  });
});

describe('Canopy Units: Pressure', () => {
  it('Pa → bar: 1e5 Pa = 1 bar', () => {
    expect(fromSI_P(1e5, 'bar')).toBeCloseTo(1, 5);
  });
  it('Pa → atm: 101325 Pa = 1 atm', () => {
    expect(fromSI_P(101325, 'atm')).toBeCloseTo(1, 5);
  });
  it('Pa → kPa: 101325 Pa = 101.325 kPa', () => {
    expect(fromSI_P(101325, 'kPa')).toBeCloseTo(101.325, 3);
  });
  it('Pa → psi: 101325 Pa ≈ 14.696 psi', () => {
    expect(fromSI_P(101325, 'psi')).toBeCloseTo(14.696, 1);
  });
  it('Pa → Pa: identity', () => {
    expect(fromSI_P(200000, 'Pa')).toBeCloseTo(200000, 0);
  });
  it('roundtrip bar', () => {
    expect(fromSI_P(toSI_P(5, 'bar'), 'bar')).toBeCloseTo(5, 5);
  });
  it('roundtrip psi', () => {
    expect(fromSI_P(toSI_P(30, 'psi'), 'psi')).toBeCloseTo(30, 5);
  });
});

describe('Canopy Units: Flow', () => {
  it('mol/s → mol/s: identity', () => {
    expect(fromSI_F(1, 'mol/s')).toBeCloseTo(1, 5);
  });
  it('mol/s → kmol/h: 1 mol/s = 3.6 kmol/h', () => {
    expect(fromSI_F(1, 'kmol/h')).toBeCloseTo(3.6, 3);
  });
  it('roundtrip kmol/h', () => {
    expect(fromSI_F(toSI_F(10, 'kmol/h'), 'kmol/h')).toBeCloseTo(10, 5);
  });
});

describe('Canopy Units: Energy', () => {
  it('W → kW: 1000 W = 1 kW', () => {
    expect(fromSI_Q(1000, 'kW')).toBeCloseTo(1, 5);
  });
  it('W → W: identity', () => {
    expect(fromSI_Q(500, 'W')).toBeCloseTo(500, 0);
  });
  it('W → BTU/h: 1000 W ≈ 3412 BTU/h', () => {
    expect(fromSI_Q(1000, 'BTU/h')).toBeCloseTo(3412.14, 0);
  });
  it('W → kcal/h: 1000 W ≈ 860 kcal/h', () => {
    expect(fromSI_Q(1000, 'kcal/h')).toBeCloseTo(860, 0);
  });
});

describe('Canopy Units: Labels', () => {
  it('tempLabel °C', () => expect(tempLabel('C')).toBe('°C'));
  it('tempLabel K', () => expect(tempLabel('K')).toBe('K'));
  it('tempLabel °F', () => expect(tempLabel('F')).toBe('°F'));
  it('pressLabel bar', () => expect(pressLabel('bar')).toBe('bar'));
  it('flowLabel mol/s', () => expect(flowLabel('mol/s')).toBe('mol/s'));
  it('energyLabel kW', () => expect(energyLabel('kW')).toBe('kW'));
});
