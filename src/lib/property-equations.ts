export interface PropertyDefinition {
  displayName: string;
  jsonKey: string;
  equationType: 'eq101' | 'eq105' | 'polynomial' | 'eq106' | 'eq102_cv' 
    | 'eq104_virial' | 'eq16_complex' | 'eq105_molar' | 'eq121' | 'eq13';
  yAxisIndex: number;
  targetUnitName: string; // Base unit for calculations and storage
  color: string;
  coeffs: string[];
  requiresTc?: boolean;
  requiresMolarMass?: boolean; // Indicates if the equation itself needs molar mass
  conversionFactor?: number; // Factor to convert equation's raw output to targetUnitName
  equationTemplate?: string;
  symbol?: string; // For display in dropdown
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
    displayName: "Vapor Pressure", jsonKey: "Vapour pressure", symbol: "P", equationType: "eq101", yAxisIndex: 0, targetUnitName: "Pa", 
    availableUnits: [
        { unit: "bar", conversionFactorFromBase: 1e-5 },
        { unit: "Pa", conversionFactorFromBase: 1 },
        { unit: "kPa", conversionFactorFromBase: 1e-3 },
        { unit: "atm", conversionFactorFromBase: 1/101325 },
        { unit: "mmHg", conversionFactorFromBase: 1/133.322 }
    ],
    color: baseColors[0], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "exp(A + B/T + C ln(T) + D T<sup>E</sup>)" 
  },
  { 
    displayName: "Liquid Density", jsonKey: "Liquid density", symbol: "ρ_L", equationType: "eq105_molar", yAxisIndex: 0, targetUnitName: "kmol/m³", 
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
    displayName: "Liquid Heat Capacity", jsonKey: "Liquid heat capacity", symbol: "Cp", equationType: "eq16_complex", yAxisIndex: 0, targetUnitName: "J/kmol/K", 
    conversionFactor: 1, // Raw eq output (with DB coeffs) assumed to be in J/kmol/K (targetUnitName)
    availableUnits: [
        { unit: "J/mol/K", conversionFactorFromBase: 0.001 },
        { unit: "J/kmol/K", conversionFactorFromBase: 1 },
        { unit: "kJ/kmol/K", conversionFactorFromBase: 0.001 },
        { unit: "kJ/mol/K", conversionFactorFromBase: 1e-6 },
        { unit: "J/kg/K", conversionFactorFromBase: { operation: 'divide_by_mw' } },
        { unit: "kJ/kg", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } },
        { unit: "J/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } },
        { unit: "kJ/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 1e-6 } }
    ],
    color: baseColors[2], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + exp(B/T + C + D T + E T<sup>2</sup>)" 
  },
  { 
    displayName: "Liquid Viscosity", jsonKey: "Liquid viscosity", symbol: "μ", equationType: "eq101", yAxisIndex: 0, targetUnitName: "Pa·s", 
    availableUnits: [
        { unit: "cP", conversionFactorFromBase: 1000 },
        { unit: "Pa·s", conversionFactorFromBase: 1 },
        { unit: "mPa·s", conversionFactorFromBase: 1000 }
    ],
    color: baseColors[3], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "exp(A + B/T + C ln(T) + D T<sup>E</sup>)" 
  },
  { 
    displayName: "Heat of Vaporization", jsonKey: "Heat of vaporization", symbol: "ΔH_v", equationType: "eq106", yAxisIndex: 0, targetUnitName: "J/kmol", 
    conversionFactor: 1000, // Raw eq output likely J/mol, convert to J/kmol
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
    displayName: "Ideal Gas Heat Capacity", jsonKey: "Ideal gas heat capacity", symbol: "Cp^0", equationType: "eq16_complex", yAxisIndex: 0, targetUnitName: "J/kmol/K", 
    conversionFactor: 1, // Raw eq output (with DB coeffs) assumed to be in J/kmol/K (targetUnitName)
    availableUnits: [
        { unit: "J/mol/K", conversionFactorFromBase: 0.001 },
        { unit: "J/kmol/K", conversionFactorFromBase: 1 },
        { unit: "kJ/kmol/K", conversionFactorFromBase: 0.001 },
        { unit: "kJ/mol/K", conversionFactorFromBase: 1e-6 },
        { unit: "J/kg/K", conversionFactorFromBase: { operation: 'divide_by_mw' } },
        { unit: "kJ/kg", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } }, 
        { unit: "J/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } },
        { unit: "kJ/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 1e-6 } }
    ],
    color: baseColors[5], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + exp(B/T + C + D T + E T<sup>2</sup>)"
  },
  { 
    displayName: "Liquid Thermal Conductivity", jsonKey: "Liquid thermal conductivity", symbol: "k_L", equationType: "eq16_complex", yAxisIndex: 0, targetUnitName: "W/m/K", 
    availableUnits: [
        { unit: "W/m/K", conversionFactorFromBase: 1 },
        { unit: "mW/m/K", conversionFactorFromBase: 1000 }
    ],
    color: baseColors[7], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + exp(B/T + C + D T + E T<sup>2</sup>)" 
  },
  { 
    displayName: "Second Virial Coefficient", jsonKey: "Second virial coefficient", symbol: "B_v", equationType: "eq104_virial", yAxisIndex: 0, targetUnitName: "m³/kmol", 
    availableUnits: [
        { unit: "cm³/mol", conversionFactorFromBase: 1000 }, // 1 m³/kmol * 10^6 cm³/m³ / 1000 mol/kmol = 1000 cm³/mol
        { unit: "m³/kmol", conversionFactorFromBase: 1 },
        { unit: "m³/kg", conversionFactorFromBase: { operation: 'divide_by_mw' } },
        { unit: "cm³/g", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 1000 } } // m³/kg * 10^6 cm³/m³ / 1000 g/kg = cm³/g * 1000
    ],
    color: baseColors[8], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + B/T + C/T<sup>2</sup> + D/T<sup>8</sup> + E/T<sup>9</sup>" 
  },
  { 
    displayName: "Solid Density", jsonKey: "Solid density", symbol: "ρ_S", equationType: "polynomial", yAxisIndex: 0, targetUnitName: "kmol/m³", 
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
    displayName: "Solid Heat Capacity", jsonKey: "Solid heat capacity", symbol: "Cp_S", equationType: "polynomial", yAxisIndex: 0, targetUnitName: "J/kmol/K", 
    conversionFactor: 1, // Raw eq output (with DB coeffs) assumed to be in J/kmol/K (targetUnitName)
    availableUnits: [
        { unit: "J/mol/K", conversionFactorFromBase: 0.001 },
        { unit: "J/kmol/K", conversionFactorFromBase: 1 },
        { unit: "kJ/kmol/K", conversionFactorFromBase: 0.001 },
        { unit: "kJ/mol/K", conversionFactorFromBase: 1e-6 },
        { unit: "J/kg/K", conversionFactorFromBase: { operation: 'divide_by_mw' } },
        { unit: "kJ/kg", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } }, 
        { unit: "J/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } },
        { unit: "kJ/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 1e-6 } }
    ],
    color: '#FFD700', coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + B T + C T<sup>2</sup> + D T<sup>3</sup> + E T<sup>4</sup>" 
  },
  { 
    displayName: "Surface Tension", jsonKey: "Surface tension", symbol: "σ", equationType: "eq16_complex", yAxisIndex: 0, targetUnitName: "N/m", 
    availableUnits: [
        { unit: "N/m", conversionFactorFromBase: 1 },
        { unit: "mN/m", conversionFactorFromBase: 1000 },
        { unit: "dyn/cm", conversionFactorFromBase: 1000 }
    ],
    color: '#00CED1', coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + exp(B/T + C + D T + E T<sup>2</sup>)" 
  },
  { 
    displayName: "Vapour Thermal Conductivity", jsonKey: "Vapour thermal conductivity", symbol: "k_V", equationType: "eq102_cv", yAxisIndex: 0, targetUnitName: "W/m/K", 
    availableUnits: [
        { unit: "W/m/K", conversionFactorFromBase: 1 },
        { unit: "mW/m/K", conversionFactorFromBase: 1000 }
    ],
    color: '#DA70D6', coeffs: ['A', 'B', 'C', 'D'], equationTemplate: "(A T<sup>B</sup>) / (1 + C/T + D/T<sup>2</sup>)" 
  },
  { 
    displayName: "Vapour Viscosity", jsonKey: "Vapour viscosity", symbol: "μ_V", equationType: "eq102_cv", yAxisIndex: 0, targetUnitName: "Pa·s", 
    availableUnits: [
        { unit: "cP", conversionFactorFromBase: 1000 },
        { unit: "Pa·s", conversionFactorFromBase: 1 },
        { unit: "mPa·s", conversionFactorFromBase: 1000 }
    ],
    color: '#6A5ACD', coeffs: ['A', 'B', 'C', 'D'], equationTemplate: "(A T<sup>B</sup>) / (1 + C/T + D/T<sup>2</sup>)" 
  },
  { 
    displayName: "Relative Static Permittivity", jsonKey: "Relative static permittivity", symbol: "ε_r", equationType: "eq121", yAxisIndex: 0, targetUnitName: "-", 
    // NOTE: Coefficients for "Water" for this property in the database (as of current data)
    // seem to produce highly incorrect values when used with calculatePolynomial.
    // Updated to use calculateEq121 based on eqno: 121 in the database.
    availableUnits: [
        { unit: "-", conversionFactorFromBase: 1 }
    ],
    color: '#FF4500', coeffs: ['A', 'B', 'C', 'D'], equationTemplate: "A + B/T + C ln(T) + D T" 
  },
  { 
    displayName: "Antoine Vapor Pressure", jsonKey: "Antoine vapor pressure", symbol: "P", equationType: "eq101", yAxisIndex: 0, targetUnitName: "Pa", 
    availableUnits: [
        { unit: "Pa", conversionFactorFromBase: 1 },
        { unit: "kPa", conversionFactorFromBase: 1e-3 },
        { unit: "bar", conversionFactorFromBase: 1e-5 },
        { unit: "atm", conversionFactorFromBase: 1/101325 },
        { unit: "mmHg", conversionFactorFromBase: 1/133.322 }
    ],
    color: '#FF6347', coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "exp(A + B/T + C ln(T) + D T<sup>E</sup>)" 
  },
  { 
    displayName: "Ideal Gas Heat Capacity (RPP)", jsonKey: "Ideal gas heat capacity (RPP)", symbol: "Cp^0_RPP", equationType: "polynomial", yAxisIndex: 0, targetUnitName: "J/kmol/K", 
    conversionFactor: 1, // Raw eq output (with DB coeffs) assumed to be in J/kmol/K (targetUnitName)
    availableUnits: [
        { unit: "J/mol/K", conversionFactorFromBase: 0.001 },
        { unit: "J/kmol/K", conversionFactorFromBase: 1 },
        { unit: "kJ/kmol/K", conversionFactorFromBase: 0.001 },
        { unit: "kJ/mol/K", conversionFactorFromBase: 1e-6 },
        { unit: "J/kg/K", conversionFactorFromBase: { operation: 'divide_by_mw' } },
        { unit: "kJ/kg", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } }, 
        { unit: "J/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 0.001 } },
        { unit: "kJ/g/K", conversionFactorFromBase: { operation: 'divide_by_mw', factor: 1e-6 } }
    ],
    color: '#4682B4', coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "A + B T + C T<sup>2</sup> + D T<sup>3</sup> + E T<sup>4</sup>" 
  },
  { 
    displayName: "Liquid Viscosity (RPS)", jsonKey: "Liquid viscosity (RPS)", symbol: "μ_RPS", equationType: "eq13", yAxisIndex: 0, targetUnitName: "Pa·s", 
    availableUnits: [
        { unit: "cP", conversionFactorFromBase: 1000 },
        { unit: "Pa·s", conversionFactorFromBase: 1 },
        { unit: "mPa·s", conversionFactorFromBase: 1000 }
    ],
    color: '#8A2BE2', coeffs: ['A','B','C'], equationTemplate: "exp(A + B T + C T<sup>2</sup>)" 
  }
];

/**
 * Calculates property using equation 101: P = exp(A + B/T + C*ln(T) + D*T^E)
 * Commonly used for Vapor Pressure.
 * Coefficients A, B, C, D, E.
 * T in Kelvin.
 * Result needs to be handled based on original units (e.g. Pa for pressure).
 */
export function calculateEq101(T: number, A: number, B: number, C: number, D: number, E: number | undefined): number | null {
  if (T <= 0) return null;
  try {
    // The form in ChemSep docs for Eq 101 is often ln(P) = A + B/T + C*ln(T) + D*T^E
    // So P = exp(A + B/T + C*ln(T) + D*T^E)
    // If E is not provided or is for T^E where E=1, D*T is used.
    // The example JSON for Vapour Pressure has E=2, so it's D*T^2.
    // The general form is D*T^E. If E is not present in JSON, it might imply E=1 or the term is not D*T^E.
    // For now, assuming E is provided if the term is T^E. If D is the last coeff and E is not given, it might be D*T.
    // The example JSON for Vapour Pressure has A, B, C, D, E.
    let termE = 0;
    if (typeof D === 'number' && typeof E === 'number') { // Check if D and E are numbers
        termE = D * Math.pow(T, E);
    } else if (typeof D === 'number' && E === undefined) { // If E is not provided, assume D is the coefficient for T (E=1)
        // This case needs to be confirmed based on how coefficients are stored for eq 101
        // For now, if E is undefined, we assume the term is D*T or the D parameter is for a different term.
        // The example has E=2 for Vapour Pressure, so this branch might not be hit for that specific property.
        // If it's just 4 coefficients A,B,C,D, then D might be for T.
        // Let's assume if E is undefined, the D term is D*T.
        // However, the example for Liquid Viscosity (eq 101) also has 5 params A,B,C,D,E.
        // So, it's safer to assume E is always present if that term exists.
        // If D is present but E is not, it's ambiguous. For safety, we'll only calculate if E is present for the T^E term.
        // Update: The prompt's example JSON for Vapour Pressure and Liquid Viscosity (eq 101) both have A,B,C,D,E.
        // So, E should generally be present.
         console.warn("calculateEq101: D is present but E is undefined. Term D*T^E cannot be calculated accurately.");
    }


    const lnP = A + (B / T) + (C * Math.log(T)) + termE;
    return Math.exp(lnP);
  } catch (e) {
    console.error("Error in calculateEq101:", e);
    return null;
  }
}

/**
 * Calculates Liquid Density using Rackett equation (modified form, eq 105): rho_L = (A / (B^(1 + (1-T/Tc)^D)))
 * A, B, D are coefficients. Tc (Critical Temperature) is C from JSON.
 * T in Kelvin.
 * Result in kmol/m^3. Multiply by Molar Mass (kg/kmol) for kg/m^3.
 */
export function calculateEq105(T: number, A: number, B: number, Tc: number, D: number, molarMass_kg_kmol: number): number | null {
  if (T <= 0 || T > Tc || B <= 0) return null; // Density is not typically defined above Tc with this form.
  try {
    const Tr = T / Tc;
    if (Tr > 1 && D > 0) { // Avoid issues with (1-Tr)^D if Tr > 1
        // This equation is generally for T < Tc
        // return null;
    }
    const exponent = 1 + Math.pow(1 - Tr, D);
    const molarDensity_kmol_m3 = A / Math.pow(B, exponent);
    return molarDensity_kmol_m3 * molarMass_kg_kmol; // Convert to kg/m^3
  } catch (e) {
    console.error("Error in calculateEq105:", e);
    return null;
  }
}

// For calculating molar density directly using Eq105-like coefficients
// ρ_molar = A / B^(1+(1-T/Tc)^D)  (output in kmol/m³)
export function calculateEq105_molar(T: number, A: number, B: number, Tc: number, D_coeff: number): number | null {
    const debugInfo = `(T=${T.toFixed(2)}, A=${A}, B=${B}, Tc=${Tc.toFixed(2)}, D=${D_coeff})`; // For logging
    if (T <= 0 || Tc <= 0 || B <= 0) {
        console.warn(`calculateEq105_molar: Returning null due to invalid initial params (T, Tc, or B <= 0). ${debugInfo}`);
        return null;
    }

    const oneMinusTr = 1 - T / Tc;
    let powerTerm;

    if (oneMinusTr < 0) {
        // If D_coeff is, for example, 0.35, then (-ve)^0.35 is complex.
        // For now, let's allow calculation but be aware it might yield NaN if not handled carefully by specific D values.
        powerTerm = Math.pow(oneMinusTr, D_coeff);
        if (isNaN(powerTerm) && oneMinusTr < 0) { // If it resulted in NaN, it might be an issue with negative base to fractional power
            console.warn(`calculateEq105_molar: Returning null due to NaN from (1-T/Tc)^D for T > Tc. ${debugInfo}, 1-Tr=${oneMinusTr.toFixed(4)}, powerTerm=${powerTerm}`);
            return null;
        }
    } else {
        powerTerm = Math.pow(oneMinusTr, D_coeff);
    }

    if (isNaN(powerTerm)) {
        console.warn(`calculateEq105_molar: Returning null due to powerTerm being NaN. ${debugInfo}, 1-Tr=${oneMinusTr.toFixed(4)}, powerTerm=${powerTerm}`);
        return null;
    }

    const term_B_exponent = 1 + powerTerm;
    const term_B = Math.pow(B, term_B_exponent);

    if (term_B === 0 || isNaN(term_B) || !isFinite(term_B)) {
        console.warn(`calculateEq105_molar: Returning null due to term_B being zero, NaN, or non-finite. ${debugInfo}, term_B_exponent=${term_B_exponent.toFixed(4)}, term_B=${term_B}`);
        return null;
    }
    const result = A / term_B;
    if (isNaN(result) || !isFinite(result)) {
        console.warn(`calculateEq105_molar: Returning null due to final result being NaN or non-finite. ${debugInfo}, result=${result}`);
        return null;
    }
    return result;
}

/**
 * Calculates property using polynomial equation (eq 16 or 100): Val = A + BT + CT^2 + DT^3 + ET^4
 * Commonly used for Heat Capacity, Surface Tension, Liquid Thermal Conductivity.
 * Coefficients A, B, C, D, E.
 * T in Kelvin.
 */
export function calculatePolynomial(T: number, A: number, B?: number | null, C?: number | null, D?: number | null, E?: number | null): number | null {
  try {
    let value = A;
    if (typeof B === 'number') value += B * T;
    if (typeof C === 'number') value += C * Math.pow(T, 2);
    if (typeof D === 'number') value += D * Math.pow(T, 3);
    if (typeof E === 'number') value += E * Math.pow(T, 4);
    return value;
  } catch (e) {
    console.error("Error in calculatePolynomial:", e);
    return null;
  }
}


/**
 * Calculates property using DWSIM's ChemSep eqno 16: Val = A + Exp(B/T + C + D*T + E*T^2)
 * Used for Liquid Heat Capacity, Liquid Thermal Conductivity, Surface Tension.
 * Coefficients A, B, C, D, E.
 * T in Kelvin.
 */
export function calculateEq16Complex(T: number, A: number, B: number, C: number, D: number, E: number): number | null {
    if (T <= 0) return null;
    // Eq 16: A + exp(B/T + C + DT + ET^2)
    return A + Math.exp(B / T + C + D * T + E * Math.pow(T, 2));
}

/**
 * Calculates Heat of Vaporization using Watson-like equation (eq 106): Hv = A * (1-Tr)^(B + C*Tr + D*Tr^2 + E*Tr^3)
 * A, B, C, D, E are coefficients. Tc is Critical Temperature.
 * T in Kelvin. Tr = T/Tc.
 * Result in J/kmol. Divide by 1000 for J/mol.
 */
export function calculateEq106(T: number, A: number, B: number, C: number, D: number, E: number, Tc: number): number | null {
  if (T <= 0 || T >= Tc) return null; // Hv is not defined at or above Tc.
  try {
    const Tr = T / Tc;
    const exponent = B + C * Tr + D * Math.pow(Tr, 2) + E * Math.pow(Tr, 3);
    const Hv_J_kmol = A * Math.pow(1 - Tr, exponent);
    return Hv_J_kmol / 1000; // Convert to J/mol
  } catch (e) {
    console.error("Error in calculateEq106:", e);
    return null;
  }
}

/**
 * Calculates property using equation 102: Val = (A * T^B) / (1 + C/T + D/T^2)
 * Commonly used for Vapour Thermal Conductivity and Vapour Viscosity.
 * Coefficients A, B, C, D.
 * T in Kelvin.
 */
export function calculateEq102_conductivity_viscosity(T: number, A: number, B: number, C: number, D: number): number | null {
  if (T <= 0) return null;
  try {
    const numerator = A * Math.pow(T, B);
    const denominator = 1 + (C / T) + (D / Math.pow(T, 2));
    if (Math.abs(denominator) < 1e-9) return null; // Avoid division by zero
    return numerator / denominator;
  } catch (e) {
    console.error("Error in calculateEq102_conductivity_viscosity:", e);
    return null;
  }
}

/**
 * Calculates Second Virial Coefficient using equation 104: B = A + B_coeff/T + C_coeff/T^2 + D_coeff/T^8 + E_coeff/T^9
 * Note: The 'B' in "B_virial = ..." is the result, A,B_coeff,C_coeff,D_coeff,E_coeff are coefficients from JSON.
 * Coefficients A, B_coeff, C_coeff, D_coeff, E_coeff.
 * T in Kelvin.
 * Result in m^3/kmol.
 */
export function calculateEq104_virial(T: number, A: number, B_coeff: number, C_coeff: number, D_coeff: number, E_coeff: number): number | null {
  if (T <= 0) return null;
  try {
    return A + (B_coeff / T) + (C_coeff / Math.pow(T, 2)) + (D_coeff / Math.pow(T, 8)) + (E_coeff / Math.pow(T, 9));
  } catch (e) {
    console.error("Error in calculateEq104_virial:", e);
    return null;
  }
}

/**
 * Parses a coefficient that might be a direct number or an object like { value: number }
 */
export function parseCoefficient(value: any): number | null {
    if (typeof value === 'number') {
        return value;
    }
    if (typeof value === 'object' && value !== null && typeof value.value === 'number') {
        return value.value;
    }
    if (typeof value === 'string') {
        const num = parseFloat(value);
        if (!isNaN(num)) return num;
    }
    return null; // Changed from undefined to null
}

/**
 * Calculates property using equation 121: Val = A + B/T + C*ln(T) + D*T
 * Commonly used for properties like Relative Static Permittivity.
 * Coefficients A, B, C, D.
 * T in Kelvin.
 */
export function calculateEq121(
  T: number,
  A: number,
  B: number,
  C: number,
  D: number
): number | null {
  if (T <= 0) {
    console.warn("calculateEq121: Temperature must be positive Kelvin.");
    return null;
  }
  try {
    const result = A + (B / T) + (C * Math.log(T)) + (D * T);

    if (isNaN(result) || !isFinite(result)) {
      console.warn(
        `calculateEq121: Result is NaN or Infinite for T=${T.toFixed(
          2
        )}K. Inputs: A=${A}, B=${B}, C=${C}, D=${D}.`
      );
      return null;
    }
    return result;
  } catch (e) {
    console.error("Error in calculateEq121:", e);
    return null;
  }
}

/**
 * Calculates property using equation 13: Val = exp(A + B T + C T^2)
 * Commonly used for RPS liquid viscosity.
 * Coefficients A, B, C.
 * T in Kelvin.
 */
export function calculateEq13(
  T: number,
  A: number,
  B: number,
  C: number
): number | null {
  if (T <= 0) {
    console.warn("calculateEq13: Temperature must be positive Kelvin.");
    return null;
  }
  try {
    const lnProperty = A + (B * T) + (C * Math.pow(T, 2));
    const result = Math.exp(lnProperty);
    if (!isFinite(result)) {
      console.warn(
        `calculateEq13: Non-finite result at T=${T}. Inputs: A=${A}, B=${B}, C=${C}.`
      );
      return null;
    }
    return result;
  } catch (e) {
    console.error("Error in calculateEq13:", e);
    return null;
  }
}

