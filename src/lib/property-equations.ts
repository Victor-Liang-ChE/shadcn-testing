export interface PropertyDefinition {
  displayName: string;
  jsonKey: string;
  yAxisIndex: number;
  targetUnitName: string; // Base unit for calculations and storage
  color: string;
  coeffs: string[];
  requiresTc?: boolean;
  requiresMolarMass?: boolean; // Indicates if the equation itself needs molar mass
  conversionFactor?: number; // Factor to convert equation's raw output to targetUnitName
  equationTemplate?: string;
  symbol?: string; // For display in dropdown
  isValueAtPressure?: boolean; // If true, "Get Value at..." section is pressure-based
  availableUnits?: Array<{ 
    unit: string; // Display unit
    // Factor to convert from targetUnitName to this display unit.
    // Can be a number for direct scaling, or an object for dynamic conversion (e.g., involving molar mass).
    conversionFactorFromBase: number | { 
        operation: 'divide_by_mw' | 'multiply_by_mw'; // MW is in kg/kmol
        factor?: number; // Additional scaling factor applied AFTER mw operation
    };
    displayName?: string; // Optional: if unit string itself isn't descriptive enough
  }>;
}

export const baseColors = ['#5470C6', '#91CC75', '#FAC858', '#EE6666', '#73C0DE', '#3BA272', '#FC8452', '#9A60B4', '#EA7CCC'];

export const propertiesToPlotConfig: PropertyDefinition[] = [
  { 
    displayName: "Vapor Pressure", jsonKey: "VaporPressure", symbol: "P", yAxisIndex: 0, targetUnitName: "Pa", 
    availableUnits: [
        { unit: "bar", conversionFactorFromBase: 1e-5 },
        { unit: "Pa", conversionFactorFromBase: 1 },
        { unit: "kPa", conversionFactorFromBase: 1e-3 },
        { unit: "MPa", conversionFactorFromBase: 1e-6 },
        { unit: "atm", conversionFactorFromBase: 1/101325 },
        { unit: "mmHg", conversionFactorFromBase: 1/133.322 }
    ],
    color: baseColors[0], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "exp(A + B/T + C ln(T) + D T<sup>E</sup>)" 
  },
  { 
    displayName: "Boiling Point", 
    jsonKey: "VaporPressure", // Uses the same data source as Vapor Pressure
    symbol: "T_b", // Conceptual symbol for the dropdown
    yAxisIndex: 0, 
    targetUnitName: "Pa", // Base unit for the *pressure* input to the solver
    isValueAtPressure: true, // Indicates this property's "Get Value at" is pressure-based
    availableUnits: [ // These are PRESSURE units for the X-axis
        { unit: "bar", conversionFactorFromBase: 1e-5 },
        { unit: "Pa", conversionFactorFromBase: 1 },
        { unit: "kPa", conversionFactorFromBase: 1e-3 },
        { unit: "MPa", conversionFactorFromBase: 1e-6 },
        { unit: "atm", conversionFactorFromBase: 1/101325 },
        { unit: "mmHg", conversionFactorFromBase: 1/133.322 }
    ],
    color: baseColors[6] || '#FC8452', // A distinct color (e.g., Orange-Red from baseColors)
    coeffs: ['A', 'B', 'C', 'D', 'E'], 
    equationTemplate: "exp(A + B/T + C ln(T) + D T<sup>E</sup>)" // Original equation string for reference
  },
  { 
    displayName: "Liquid Density", jsonKey: "LiquidDensity", symbol: "ρ_L", yAxisIndex: 0, targetUnitName: "kmol/m³", 
    availableUnits: [
        { unit: "kmol/m³", conversionFactorFromBase: 1 },
        { unit: "mol/L", conversionFactorFromBase: 1 }, // 1 kmol/m³ = 1 mol/L
        { unit: "mol/m³", conversionFactorFromBase: 1000 },
        { unit: "kg/m³", conversionFactorFromBase: { operation: 'multiply_by_mw' } },
        { unit: "g/m³", conversionFactorFromBase: { operation: 'multiply_by_mw', factor: 1000 } },
        { unit: "g/L", conversionFactorFromBase: { operation: 'multiply_by_mw' } }, // kg/m³ is numerically equal to g/L
        { unit: "g/cm³", conversionFactorFromBase: { operation: 'multiply_by_mw', factor: 0.001 } } // kg/m³ * 0.001 = g/cm³
    ],
    color: baseColors[1], coeffs: ['A', 'B', 'C', 'D'], requiresMolarMass: false, requiresTc: true, conversionFactor: 1, equationTemplate: "A / B<sup>(1+(1-T/Tc)<sup>D</sup>)</sup>" 
  },
  { 
    displayName: "Liquid Heat Capacity", jsonKey: "LiquidHeatCapacityCp", symbol: "Cp", yAxisIndex: 0, targetUnitName: "J/kmol/K", 
    conversionFactor: 1, 
    availableUnits: [
        { unit: "J/mol/K", conversionFactorFromBase: 0.001 },
        { unit: "J/kmol/K", conversionFactorFromBase: 1 },
        { unit: "kJ/kmol/K", conversionFactorFromBase: 0.001 },
        { unit: "kJ/mol/K", conversionFactorFromBase: 1e-6 },
        { unit: "J/kg/K", conversionFactorFromBase: { operation: 'divide_by_mw' } },
        { unit: "kJ/kg/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } }, // Corrected kJ/kg/K
        { unit: "J/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } },
        { unit: "kJ/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 1e-6 } }
    ],
    color: baseColors[2], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + exp(B/T + C + D T + E T<sup>2</sup>)" 
  },
  { 
    displayName: "Liquid Viscosity", jsonKey: "LiquidViscosity", symbol: "μ", yAxisIndex: 0, targetUnitName: "Pa·s", 
    availableUnits: [
        { unit: "cP", conversionFactorFromBase: 1000 },
        { unit: "Pa·s", conversionFactorFromBase: 1 },
        { unit: "mPa·s", conversionFactorFromBase: 1000 }
    ],
    color: baseColors[3], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "exp(A + B/T + C ln(T) + D T<sup>E</sup>)" 
  },
  { 
    displayName: "Heat of Vaporization", jsonKey: "HeatOfVaporization", symbol: "ΔH_v", yAxisIndex: 0, targetUnitName: "J/kmol", 
    conversionFactor: 1, // Assuming eq output is already J/kmol, or adjust if it's J/mol from DB
    availableUnits: [
        { unit: "J/mol", conversionFactorFromBase: 0.001 },
        { unit: "J/kmol", conversionFactorFromBase: 1 },
        { unit: "kJ/kmol", conversionFactorFromBase: 0.001 },
        { unit: "kJ/mol", conversionFactorFromBase: 1e-6 },
        { unit: "J/kg", conversionFactorFromBase: { operation: 'divide_by_mw' } },
        { unit: "kJ/kg", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } },
        { unit: "J/g", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } },
        { unit: "kJ/g", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 1e-6 } }
    ],
    color: baseColors[4], coeffs: ['A', 'B', 'C', 'D', 'E'], requiresTc: true, equationTemplate: "A(1-T/Tc)<sup>(B+C(T/Tc)+D(T/Tc)<sup>2</sup>+E(T/Tc)<sup>3</sup>)</sup>" 
  },
  { 
    displayName: "Ideal Gas Heat Capacity", jsonKey: "IdealGasHeatCapacityCp", symbol: "Cp^0", yAxisIndex: 0, targetUnitName: "J/kmol/K", 
    conversionFactor: 1, 
    availableUnits: [
        { unit: "J/mol/K", conversionFactorFromBase: 0.001 },
        { unit: "J/kmol/K", conversionFactorFromBase: 1 },
        { unit: "kJ/kmol/K", conversionFactorFromBase: 0.001 },
        { unit: "kJ/mol/K", conversionFactorFromBase: 1e-6 },
        { unit: "J/kg/K", conversionFactorFromBase: { operation: 'divide_by_mw' } },
        { unit: "kJ/kg/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } }, // Corrected kJ/kg/K
        { unit: "J/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } },
        { unit: "kJ/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 1e-6 } }
    ],
    color: baseColors[5], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + exp(B/T + C + D T + E T<sup>2</sup>)"
  },
  { 
    displayName: "Liquid Thermal Conductivity", jsonKey: "LiquidThermalConductivity", symbol: "k_L", yAxisIndex: 0, targetUnitName: "W/m/K", 
    availableUnits: [
        { unit: "W/m/K", conversionFactorFromBase: 1 },
        { unit: "mW/m/K", conversionFactorFromBase: 1000 }
    ],
    color: baseColors[7], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + exp(B/T + C + D T + E T<sup>2</sup>)" 
  },
  { 
    displayName: "Second Virial Coefficient", jsonKey: "SecondVirialCoefficient", symbol: "B_v", yAxisIndex: 0, targetUnitName: "m³/kmol", 
    availableUnits: [
        { unit: "cm³/mol", conversionFactorFromBase: 1000 }, 
        { unit: "m³/kmol", conversionFactorFromBase: 1 },
        { unit: "m³/mol", conversionFactorFromBase: 0.001 },
        { unit: "L/mol", conversionFactorFromBase: 1 }, // 1 m³/kmol = 1 L/mol
        { unit: "m³/kg", conversionFactorFromBase: { operation: 'divide_by_mw' } },
        { unit: "cm³/g", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 1000 } } 
    ],
    color: baseColors[8], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + B/T + C/T<sup>2</sup> + D/T<sup>8</sup> + E/T<sup>9</sup>" 
  },
  { 
    displayName: "Solid Density", jsonKey: "SolidDensity", symbol: "ρ_S", yAxisIndex: 0, targetUnitName: "kmol/m³", 
    availableUnits: [
        { unit: "kmol/m³", conversionFactorFromBase: 1 },
        { unit: "mol/L", conversionFactorFromBase: 1 },
        { unit: "mol/m³", conversionFactorFromBase: 1000 },
        { unit: "kg/m³", conversionFactorFromBase: { operation: 'multiply_by_mw' } },
        { unit: "g/m³", conversionFactorFromBase: { operation: 'multiply_by_mw', factor: 1000 } },
        { unit: "g/L", conversionFactorFromBase: { operation: 'multiply_by_mw' } },
        { unit: "g/cm³", conversionFactorFromBase: { operation: 'multiply_by_mw', factor: 0.001 } }
    ],
    color: '#808080', coeffs: ['A', 'B'], requiresMolarMass: false, conversionFactor: 1, equationTemplate: "A + B T" 
  },
  { 
    displayName: "Solid Heat Capacity", jsonKey: "SolidHeatCapacityCp", symbol: "Cp_S", yAxisIndex: 0, targetUnitName: "J/kmol/K", 
    conversionFactor: 1, 
    availableUnits: [
        { unit: "J/mol/K", conversionFactorFromBase: 0.001 },
        { unit: "J/kmol/K", conversionFactorFromBase: 1 },
        { unit: "kJ/kmol/K", conversionFactorFromBase: 0.001 },
        { unit: "kJ/mol/K", conversionFactorFromBase: 1e-6 },
        { unit: "J/kg/K", conversionFactorFromBase: { operation: 'divide_by_mw' } },
        { unit: "kJ/kg/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } }, // Corrected kJ/kg/K
        { unit: "J/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } },
        { unit: "kJ/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 1e-6 } }
    ],
    color: '#FFD700', coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + B T + C T<sup>2</sup> + D T<sup>3</sup> + E T<sup>4</sup>" 
  },
  { 
    displayName: "Surface Tension", jsonKey: "SurfaceTension", symbol: "σ", yAxisIndex: 0, targetUnitName: "N/m", 
    availableUnits: [
        { unit: "N/m", conversionFactorFromBase: 1 },
        { unit: "mN/m", conversionFactorFromBase: 1000 },
        { unit: "dyn/cm", conversionFactorFromBase: 1000 }
    ],
    color: '#00CED1', coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + exp(B/T + C + D T + E T<sup>2</sup>)" 
  },
  { 
    displayName: "Vapor Thermal Conductivity", jsonKey: "VaporThermalConductivity", symbol: "k_V", yAxisIndex: 0, targetUnitName: "W/m/K", 
    availableUnits: [
        { unit: "W/m/K", conversionFactorFromBase: 1 },
        { unit: "mW/m/K", conversionFactorFromBase: 1000 }
    ],
    color: '#DA70D6', coeffs: ['A', 'B', 'C', 'D'], equationTemplate: "(A T<sup>B</sup>) / (1 + C/T + D/T<sup>2</sup>)" 
  },
  { 
    displayName: "Vapor Viscosity", jsonKey: "VaporViscosity", symbol: "μ_V", yAxisIndex: 0, targetUnitName: "Pa·s", 
    availableUnits: [
        { unit: "cP", conversionFactorFromBase: 1000 },
        { unit: "Pa·s", conversionFactorFromBase: 1 },
        { unit: "mPa·s", conversionFactorFromBase: 1000 }
    ],
    color: '#6A5ACD', coeffs: ['A', 'B', 'C', 'D'], equationTemplate: "(A T<sup>B</sup>) / (1 + C/T + D/T<sup>2</sup>)" 
  },
  { 
    displayName: "Relative Static Permittivity", jsonKey: "RelativeStaticPermittivity", symbol: "ε_r", yAxisIndex: 0, targetUnitName: "-", 
    // NOTE: Coefficients for "Water" for this property in the database (as of current data)
    // seem to produce highly incorrect values when used with calculatePolynomial.
    // Updated to use calculateEq121 based on eqno: 121 in the database.
    availableUnits: [
        { unit: "-", conversionFactorFromBase: 1 }
    ],
    color: '#FF4500', coeffs: ['A', 'B', 'C', 'D'], equationTemplate: "A + B/T + C ln(T) + D T" 
  },
  { 
    displayName: "Antoine Vapor Pressure", jsonKey: "AntoineVaporPressure", symbol: "P", yAxisIndex: 0, targetUnitName: "Pa", 
    availableUnits: [
        { unit: "bar", conversionFactorFromBase: 1e-5 },
        { unit: "Pa", conversionFactorFromBase: 1 },
        { unit: "kPa", conversionFactorFromBase: 1e-3 },
        { unit: "MPa", conversionFactorFromBase: 1e-6 },
        { unit: "atm", conversionFactorFromBase: 1/101325 },
        { unit: "mmHg", conversionFactorFromBase: 1/133.322 }
    ],
    color: '#FF6347', coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "exp(A + B/T + C ln(T) + D T<sup>E</sup>)" 
  },
  { 
    displayName: "Ideal Gas Heat Capacity (RPP)", jsonKey: "RPPHeatCapacityCp", symbol: "Cp^0_RPP", yAxisIndex: 0, targetUnitName: "J/kmol/K", 
    conversionFactor: 1, 
    availableUnits: [
        { unit: "J/mol/K", conversionFactorFromBase: 0.001 },
        { unit: "J/kmol/K", conversionFactorFromBase: 1 },
        { unit: "kJ/kmol/K", conversionFactorFromBase: 0.001 },
        { unit: "kJ/mol/K", conversionFactorFromBase: 1e-6 },
        { unit: "J/kg/K", conversionFactorFromBase: { operation: 'divide_by_mw' } },
        { unit: "kJ/kg/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } }, // Corrected kJ/kg/K
        { unit: "J/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } },
        { unit: "kJ/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 1e-6 } }
    ],
    color: '#4682B4', coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + B T + C T<sup>2</sup> + D T<sup>3</sup> + E T<sup>4</sup>" 
  },
  { 
    displayName: "Liquid Viscosity (RPS)", jsonKey: "LiquidViscosityRPS", symbol: "μ_RPS", yAxisIndex: 0, targetUnitName: "Pa·s", 
    availableUnits: [
        { unit: "cP", conversionFactorFromBase: 1000 },
        { unit: "Pa·s", conversionFactorFromBase: 1 },
        { unit: "mPa·s", conversionFactorFromBase: 1000 }
    ],
    color: '#8A2BE2', coeffs: ['A','B','C'], equationTemplate: "exp(A + B T + C T<sup>2</sup>)" 
  }
];

//
// All other exports (PropertyDefinition, baseColors, propertiesToPlotConfig) remain the same
//

// NEW: Master dispatcher function
/**
 * Calculates a property by dynamically calling the correct equation function based on eqno.
 * @param eqno The equation number as a string.
 * @param T Temperature in Kelvin.
 * @param coeffs An object containing the coefficients (A, B, C, D, E).
 * @param Tc Critical Temperature in Kelvin (required for some equations).
 * @returns The calculated property value, or null if the equation is not supported or calculation fails.
 */
export function calculatePropertyByEquation(
  eqno: string,
  T: number,
  coeffs: { [key: string]: number | null | undefined },
  Tc?: number | null
): number | null {
  const A = coeffs.A ?? 0;
  const B = coeffs.B ?? 0;
  const C = coeffs.C ?? 0;
  const D = coeffs.D ?? 0;
  const E = coeffs.E ?? 0;
  const eq = parseInt(eqno, 10);

  // A Tr value that is safe to calculate even if Tc is not provided for all equations.
  const Tr = (Tc && T) ? T / Tc : 0;

  switch (eq) {
    case 1: return A;
    case 2: return A + B * T;
    case 3: return A + B * T + C * Math.pow(T, 2);
    case 4: return A + B * T + C * Math.pow(T, 2) + D * Math.pow(T, 3);
    case 5: return A + B * T + C * Math.pow(T, 2) + D * Math.pow(T, 3) + E * Math.pow(T, 4);
    case 6: return A + B * T + C * Math.pow(T, 2) + D * Math.pow(T, 3) + E / Math.pow(T, 2);
    case 10: return Math.exp(A - B / (T + C));
    case 11: return Math.exp(A);
    case 12: return Math.exp(A + B * T);
    case 13: return Math.exp(A + B * T + C * Math.pow(T, 2));
    case 14: return Math.exp(A + B * T + C * Math.pow(T, 2) + D * Math.pow(T, 3));
    case 15: return Math.exp(A + B * T + C * Math.pow(T, 2) + D * Math.pow(T, 3) + E * Math.pow(T, 4));
    case 16: return A + Math.exp(B / T + C + D * T + E * Math.pow(T, 2));
    case 17: return A + Math.exp(B + C * T + D * Math.pow(T, 2) + E * Math.pow(T, 3));
    case 45: return A * T + B * Math.pow(T, 2) / 2 + C * Math.pow(T, 3) / 3 + D * Math.pow(T, 4) / 4 + E * Math.pow(T, 5) / 5;
    case 75: return B + 2 * C * T + 3 * D * Math.pow(T, 2) + 4 * E * Math.pow(T, 3);
    case 100: return A + B * T + C * Math.pow(T, 2) + D * Math.pow(T, 3) + E * Math.pow(T, 4);
    case 101: return Math.exp(A + B / T + C * Math.log(T) + D * Math.pow(T, E));
    case 102: return A * Math.pow(T, B) / (1 + C / T + D / Math.pow(T, 2));
    case 103: return A + B * Math.exp(-C / Math.pow(T, D));
    case 104: return A + B / T + C / Math.pow(T, 3) + D / Math.pow(T, 8) + E / Math.pow(T, 9);
    case 105: return A / Math.pow(B, 1 + Math.pow(1 - T / C, D));
    case 106: return (Tc) ? A * Math.pow(1 - Tr, B + C * Tr + D * Math.pow(Tr, 2) + E * Math.pow(Tr, 3)) : null; // Tc is required
    case 107: return A + B * Math.pow(C / T / Math.sinh(C / T), 2) + D * Math.pow(E / T / Math.cosh(E / T), 2);
    case 114: return (Tc) ? (Math.pow(A, 2) * (1 - Tr)) + B - (2 * A * C * (1 - Tr)) - (A * D * Math.pow(1 - Tr, 2)) - (Math.pow(C, 2) * Math.pow(1 - Tr, 3) / 3) - (C * D * Math.pow(1 - Tr, 4) / 2) - (Math.pow(D, 2) * Math.pow(1 - Tr, 5) / 5): null;
    case 115: return Math.exp(A + B / T + C * Math.log(T) + D * Math.pow(T, 2) + E / Math.pow(T, 2));
    case 116: return (Tc) ? A + B * Math.pow(1 - Tr, 0.35) + C * Math.pow(1 - Tr, 2 / 3) + D * (1 - Tr) + E * Math.pow(1 - Tr, 4 / 3) : null;
    case 117: return A * T + B * (C / T) / Math.tanh(C / T) - D * (E / T) / Math.tanh(E / T);
    case 119: return Math.exp(A / T + B + C * T + D * Math.pow(T, 2) + E * Math.log(T));
    case 120: return A - B / (T + C);
    case 121: return A + B / T + C * Math.log(T) + D * Math.pow(T, E);
    case 122: return A + B / T + C * Math.log(T) + D * Math.pow(T, 2) + E / Math.pow(T, 2);
    case 207: return Math.exp(A - B / (T + C));
    case 208: return Math.pow(10, A - B / (T + C));
    case 209: return Math.pow(10, A * (1 / T - 1 / B));
    case 210: return Math.pow(10, A + B / T + C * T + D * Math.pow(T, 2));
    case 211: return A * Math.pow((B - T) / (B - C), D);
    case 212: return (E) ? Math.exp((E / T) * (A * (1 - T / E) + B * Math.pow(1 - T / E, 1.5) + C * Math.pow(1 - T / E, 3) + D * Math.pow(1 - T / E, 6))) : null;
    case 213: return (E) ? (E / T) * (A * (1 - T / E) + B * Math.pow(1 - T / E, 1.5) + C * Math.pow(1 - T / E, 3) + D * Math.pow(1 - T / E, 6)) : null;
    case 221: return -B / Math.pow(T, 2) + C / T + D * E * Math.pow(T, E - 1);
    case 230: return -B / Math.pow(T, 2) + C / T + D - 2 * E / Math.pow(T, 3);
    case 231: return B - C / Math.pow(T - D, 2);
    default:
        console.warn(`Equation number ${eqno} is not supported.`);
        return null;
  }
}

/**
 * Helper function to calculate the D*T^E term and its derivative for boiling point calculation.
 * This is kept internal to the boiling point calculation logic.
 */
function _calculate_dte_term_and_derivative_for_boiling_point(
  T: number,
  coeffD: number,
  coeffE_exponent: number | undefined
): { value: number; derivative: number } {
  let value = 0;
  let derivative = 0;

  // Only calculate if coeffE_exponent is a number.
  // If coeffE_exponent is undefined, the D*T^E term is considered absent.
  if (typeof coeffE_exponent === 'number') {
    if (coeffE_exponent === 0) {
      value = coeffD; // D * T^0 = D
      derivative = 0; // d(D)/dT = 0
    } else {
      value = coeffD * Math.pow(T, coeffE_exponent);
      derivative = coeffD * coeffE_exponent * Math.pow(T, coeffE_exponent - 1);
    }
  }
  return { value, derivative };
}

/**
 * Calculates Boiling Point Temperature for a given pressure using Equation 101 form.
 * Solves for T in P = exp(A + B/T + C*ln(T) + D*T^E) using Newton-Raphson method.
 * @param targetPressurePa Target pressure in Pascals.
 * @param coeffA Coefficient A.
 * @param coeffB Coefficient B (for B/T term).
 * @param coeffC Coefficient C (for C*ln(T) term).
 * @param coeffD Coefficient D (for D*T^E term).
 * @param coeffE_exponent Exponent E (for D*T^E term). Can be undefined if the term is not present.
 * @param initialTempGuessK Initial guess for temperature in Kelvin.
 * @param maxIter Maximum number of iterations for Newton-Raphson.
 * @param tolerance Convergence tolerance.
 * @returns Calculated temperature in Kelvin, or null if not converged or error.
 */
export function calculateBoilingPointEq101(
  targetPressurePa: number,
  coeffA: number,
  coeffB: number,
  coeffC: number,
  coeffD: number, // This is the coefficient for the T^E term
  coeffE_exponent: number | undefined, // This is the E exponent
  initialTempGuessK: number = 373.15,
  maxIter: number = 100,
  tolerance: number = 1e-7
): number | null {
  if (targetPressurePa <= 0) {
    console.warn("calculateBoilingPointEq101: Target pressure must be positive.");
    return null;
  }
  if (initialTempGuessK <= 0) {
    console.warn("calculateBoilingPointEq101: Initial temperature guess must be positive.");
    return null; // Or use a default positive guess
  }

  const lnP_target = Math.log(targetPressurePa);
  let T_current = initialTempGuessK;

  for (let i = 0; i < maxIter; i++) {
    if (T_current <= 0) {
      // Temperature became non-positive, which is invalid for log(T) or T in denominator
      console.warn(`calculateBoilingPointEq101: Temperature became non-positive during iteration (${T_current.toFixed(2)} K).`);
      return null;
    }

    const { value: term_DTE, derivative: deriv_term_DTE } = _calculate_dte_term_and_derivative_for_boiling_point(
      T_current,
      coeffD,
      coeffE_exponent
    );

    // f(T) = A + B/T + C*ln(T) + D*T^E - ln(P_target)
    const f_T = coeffA + (coeffB / T_current) + (coeffC * Math.log(T_current)) + term_DTE - lnP_target;

    if (isNaN(f_T) || !isFinite(f_T)) {
        console.warn(`calculateBoilingPointEq101: f(T) is NaN or non-finite at T=${T_current.toFixed(2)}K.`);
        return null;
    }

    if (Math.abs(f_T) < tolerance) {
      return T_current; // Converged
    }

    // f'(T) = -B/T^2 + C/T + D*E*T^(E-1)
    const f_prime_T = (-coeffB / Math.pow(T_current, 2)) + (coeffC / T_current) + deriv_term_DTE;

    if (isNaN(f_prime_T) || !isFinite(f_prime_T)) {
        console.warn(`calculateBoilingPointEq101: f'(T) is NaN or non-finite at T=${T_current.toFixed(2)}K.`);
        return null;
    }

    if (Math.abs(f_prime_T) < 1e-10) {
      // Derivative is too small, Newton-Raphson may fail or oscillate
      console.warn(`calculateBoilingPointEq101: Derivative too small at T=${T_current.toFixed(2)}K. f'(T)=${f_prime_T}.`);
      return null;
    }

    const T_next = T_current - f_T / f_prime_T;

    // Basic step damping: if T_next is drastically different or invalid, adjust or fail.
    // For simplicity, we'll rely on T_current <= 0 check at the start of the next iteration.
    // More advanced solvers might limit the step size or use bounds.

    T_current = T_next;
  }

  console.warn(`calculateBoilingPointEq101: Did not converge after ${maxIter} iterations for P=${targetPressurePa} Pa.`);
  return null; // Did not converge
}

/**
 * Parses a coefficient that might be a direct number or an object like { value: number }
 */
export function parseCoefficient(value: any): number | null {
    if (typeof value === 'number') {
        return value;
    }
  if (typeof value === 'object' && value !== null) {
    if (typeof value.value === 'number') {
      return value.value;
    }
    if (typeof value.value === 'string') {
      const n = parseFloat(value.value);
      if (!isNaN(n)) return n;
    }
  }
    if (typeof value === 'string') {
        const num = parseFloat(value);
        if (!isNaN(num)) return num;
    }
    return null; // Changed from undefined to null
}

