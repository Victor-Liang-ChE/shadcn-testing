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
  try {
    const expTerm = Math.exp((B / T) + C + (D * T) + (E * Math.pow(T, 2)));
    return A + expTerm;
  } catch (e) {
    console.error("Error in calculateEq16Complex:", e);
    return null;
  }
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
export function parseCoefficient(coeff: any): number | null { // Changed return type to number | null
    if (typeof coeff === 'number') {
        return coeff;
    }
    if (typeof coeff === 'object' && coeff !== null && typeof coeff.value === 'number') {
        return coeff.value;
    }
    if (typeof coeff === 'string') {
        const num = parseFloat(coeff);
        if (!isNaN(num)) return num;
    }
    return null; // Changed from undefined to null
}

