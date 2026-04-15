// ─── Unit Conversions for Canopy ─────────────────────────────────────────────
// All internal values are SI. These convert for display only.

import type { UnitSystem } from './store';

// Temperature
export function fromSI_T(T_K: number, unit: UnitSystem['temperature']): number {
  switch (unit) {
    case 'K': return T_K;
    case 'C': return T_K - 273.15;
    case 'F': return (T_K - 273.15) * 9 / 5 + 32;
  }
}
export function toSI_T(val: number, unit: UnitSystem['temperature']): number {
  switch (unit) {
    case 'K': return val;
    case 'C': return val + 273.15;
    case 'F': return (val - 32) * 5 / 9 + 273.15;
  }
}

// Pressure
export function fromSI_P(P_Pa: number, unit: UnitSystem['pressure']): number {
  switch (unit) {
    case 'Pa': return P_Pa;
    case 'kPa': return P_Pa / 1e3;
    case 'bar': return P_Pa / 1e5;
    case 'atm': return P_Pa / 101325;
    case 'psi': return P_Pa / 6894.757;
  }
}
export function toSI_P(val: number, unit: UnitSystem['pressure']): number {
  switch (unit) {
    case 'Pa': return val;
    case 'kPa': return val * 1e3;
    case 'bar': return val * 1e5;
    case 'atm': return val * 101325;
    case 'psi': return val * 6894.757;
  }
}

// Molar/mass flow
export function fromSI_F(F_molps: number, unit: UnitSystem['flow'], mw?: number): number {
  switch (unit) {
    case 'mol/s': return F_molps;
    case 'kmol/h': return F_molps * 3.6;
    case 'kg/s': return F_molps * (mw ?? 1) / 1000;
    case 'kg/h': return F_molps * (mw ?? 1) / 1000 * 3600;
  }
}
export function toSI_F(val: number, unit: UnitSystem['flow'], mw?: number): number {
  switch (unit) {
    case 'mol/s': return val;
    case 'kmol/h': return val / 3.6;
    case 'kg/s': return val / ((mw ?? 1) / 1000);
    case 'kg/h': return val / ((mw ?? 1) / 1000 * 3600);
  }
}

// Energy (duty)
export function fromSI_Q(Q_W: number, unit: UnitSystem['energy']): number {
  switch (unit) {
    case 'W': return Q_W;
    case 'kW': return Q_W / 1000;
    case 'BTU/h': return Q_W * 3.412142;
    case 'kcal/h': return Q_W * 0.860421;
  }
}

// Unit labels
export function tempLabel(unit: UnitSystem['temperature']): string {
  switch (unit) {
    case 'K': return 'K';
    case 'C': return '°C';
    case 'F': return '°F';
  }
}
export function pressLabel(unit: UnitSystem['pressure']): string { return unit; }
export function flowLabel(unit: UnitSystem['flow']): string { return unit; }
export function energyLabel(unit: UnitSystem['energy']): string { return unit; }
