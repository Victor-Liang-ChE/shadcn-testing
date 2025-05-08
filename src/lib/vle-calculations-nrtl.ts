import type { SupabaseClient } from '@supabase/supabase-js';
import type { AntoineParams, CompoundData, BubbleDewResult } from './vle-types'; // Import shared types
import { calculatePsat_Pa } from './vle-calculations-unifac'; // Psat can still be imported from unifac or moved to a common util

// const R_J_molK = 8.31446261815324; // Gas constant in J/mol·K - REMOVED
const R_cal_molK = 1.98720425864083; // Gas constant in cal/mol·K - ADDED

export interface NrtlInteractionParams {
    A12: number;    // NRTL parameter A12 (e.g., J/mol)
    A21: number;    // NRTL parameter A21 (e.g., J/mol)
    alpha: number; // NRTL non-randomness factor (alpha12 = alpha21)
    casn1?: string; // For reference
    casn2?: string; // For reference
}

/**
 * Fetches NRTL interaction parameters for a binary pair from Supabase.
 */
export async function fetchNrtlParameters(
    supabase: SupabaseClient,
    casn1: string,
    casn2: string
): Promise<NrtlInteractionParams> {
    if (!casn1 || !casn2) {
        throw new Error("CAS numbers for both compounds are required to fetch NRTL parameters.");
    }
    if (casn1 === casn2) { 
        console.warn("fetchNrtlParameters called for identical CAS numbers. Returning placeholder params.");
        return { A12: 0, A21: 0, alpha: 0.1, casn1, casn2 }; 
    }

    console.log(`Fetching NRTL parameters for pair: casn1='${casn1}' (type: ${typeof casn1}, len: ${casn1?.length}), casn2='${casn2}' (type: ${typeof casn2}, len: ${casn2?.length})`);

    const query = supabase
        .from('nrtl parameters')
        .select('"A12", "A21", "alpha12", "CASN1", "CASN2"')
        .or(`and(CASN1.eq.${casn1},CASN2.eq.${casn2}),and(CASN1.eq.${casn2},CASN2.eq.${casn1})`)
        .limit(1);

    // For debugging, you could log the query string if the library supports it, e.g., query.toString() or similar.
    // console.log("Supabase NRTL query:", query.toString()); // This might not be a standard Supabase JS client feature.

    const { data, error } = await query;

    // ADDED: Log raw Supabase response
    console.log("Supabase NRTL query response - data:", JSON.stringify(data, null, 2));
    console.log("Supabase NRTL query response - error:", JSON.stringify(error, null, 2));


    if (error) {
        // Log the detailed error before throwing a generic one
        console.error(`Supabase NRTL query raw error object:`, error);
        throw new Error(`Supabase NRTL query error: ${error.message}`);
    }

    if (!data || data.length === 0) {
        throw new Error(`NRTL parameters not found for pair ${casn1}/${casn2}. Check 'nrtl parameters' table.`);
    }

    const params = data[0];
    if (params.CASN1 === casn1) {
        // Order matches: casn1 is CASN1, casn2 is CASN2
        return {
            A12: params.A12,
            A21: params.A21,
            alpha: params.alpha12, // Assuming alpha12 in table is the symmetric alpha
            casn1: params.CASN1,
            casn2: params.CASN2,
        };
    } else {
        // Order is reversed: casn1 is CASN2, casn2 is CASN1
        // So, the table's A12 is our A21, and table's A21 is our A12
        return {
            A12: params.A21,
            A21: params.A12,
            alpha: params.alpha12, // Assuming alpha12 in table is the symmetric alpha
            casn1: params.CASN2, // original casn1 from table
            casn2: params.CASN1, // original casn2 from table
        };
    }
}

/**
 * Calculates NRTL activity coefficients gamma_i for a binary mixture.
 * @param components Array of CompoundData objects (currently unused in this NRTL gamma calculation but kept for signature consistency).
 * @param x Array of mole fractions [x1, x2].
 * @param T_K Temperature in Kelvin.
 * @param params NRTL interaction parameters.
 * @returns A tuple [gamma1, gamma2] or null if calculation fails.
 */
export function calculateNrtlGamma(
    components: CompoundData[], // Added for signature consistency, though not directly used in this specific gamma calculation
    x: number[],
    T_K: number,
    params: NrtlInteractionParams
): [number, number] | null {
    if (x.length !== 2) return null;
    const x1 = x[0];
    const x2 = x[1];

    if (Math.abs(x1 + x2 - 1.0) > 1e-6) {
        console.error("Mole fractions do not sum to 1 in calculateNrtlGamma:", x1, x2);
        return null;
    }
     // Handle pure components: gamma should be 1
    if (Math.abs(x1 - 1.0) < 1e-9 || Math.abs(x2 - 1.0) < 1e-9) {
        return [1.0, 1.0];
    }
    if (!params) {
        console.error("NRTL parameters are missing in calculateNrtlGamma");
        return null;
    }


    const { A12, A21, alpha } = params; // A12 and A21 are now assumed to be in cal/mol

    const tau12 = A12 / (R_cal_molK * T_K); // USE R_cal_molK
    const tau21 = A21 / (R_cal_molK * T_K); // USE R_cal_molK

    const G12 = Math.exp(-alpha * tau12);
    const G21 = Math.exp(-alpha * tau21);

    const term1_denom = x1 + x2 * G21;
    const term2_denom = x2 + x1 * G12;

    if (term1_denom === 0 || term2_denom === 0) {
        console.warn("NRTL calculation encountered division by zero. Check parameters/composition.");
        return null; // Avoid division by zero
    }

    const lnGamma1 = x2 * x2 * (tau21 * (G21 / term1_denom) * (G21 / term1_denom) + tau12 * G12 / (term2_denom * term2_denom));
    const lnGamma2 = x1 * x1 * (tau12 * (G12 / term2_denom) * (G12 / term2_denom) + tau21 * G21 / (term1_denom * term1_denom));
    
    const gamma1 = Math.exp(lnGamma1);
    const gamma2 = Math.exp(lnGamma2);

    if (isNaN(gamma1) || isNaN(gamma2) || !isFinite(gamma1) || !isFinite(gamma2)) {
        console.warn(`NRTL gamma calculation resulted in NaN/Infinity. x1=${x1.toFixed(4)}, T=${T_K.toFixed(2)}K, A12=${A12}, A21=${A21}, alpha=${alpha}, tau12=${tau12.toFixed(3)}, tau21=${tau21.toFixed(3)}, G12=${G12.toExponential(3)}, G21=${G21.toExponential(3)}`);
        return null;
    }

    return [gamma1, gamma2];
}


/**
 * NRTL-based Bubble Point Temperature Calculation.
 * Iteratively solves for T given P_total and x1.
 */
export function calculateBubbleTemperatureNrtl(
    components: CompoundData[], // Uses imported CompoundData
    x1_feed: number,
    P_system_Pa: number,
    nrtlParams: NrtlInteractionParams,
    initialTempGuess_K: number = 350,
    maxIter: number = 50,
    tolerance: number = 1e-4
): BubbleDewResult | null { // Uses imported BubbleDewResult
    console.log(`NRTL Bubble T: x1=${x1_feed.toFixed(3)}, P=${(P_system_Pa/1000).toFixed(1)}kPa, Init_T=${initialTempGuess_K.toFixed(1)}K`);
    let T = initialTempGuess_K;
    const x = [x1_feed, 1 - x1_feed];
    let y = [...x]; // Initial guess for vapor composition

    if (!components[0]?.antoine || !components[1]?.antoine) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: NaN, P_Pa: P_system_Pa, error: "NRTL: Antoine parameters missing for one or both components.", calculationType: 'bubbleT' };
    }

    for (let iter = 0; iter < maxIter; iter++) {
        const Psat1 = calculatePsat_Pa(components[0].antoine, T, components[0].name);
        const Psat2 = calculatePsat_Pa(components[1].antoine, T, components[1].name);
        if (isNaN(Psat1) || isNaN(Psat2)) {
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T, P_Pa: P_system_Pa, error: `NRTL: Psat calculation failed at T=${T.toFixed(1)}K`, calculationType: 'bubbleT', iterations: iter };
        }

        const gamma = calculateNrtlGamma(components, x, T, nrtlParams);
        if (!gamma) {
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T, P_Pa: P_system_Pa, error: `NRTL: Gamma calculation failed at T=${T.toFixed(1)}K`, calculationType: 'bubbleT', iterations: iter };
        }

        // P_calc = y1*P + y2*P = (gamma1*x1*Psat1) + (gamma2*x2*Psat2)
        // Target: P_calc = P_system_Pa
        // f(T) = gamma1*x1*Psat1 + gamma2*x2*Psat2 - P_system_Pa = 0
        const P_calc = gamma[0] * x[0] * Psat1 + gamma[1] * x[1] * Psat2;
        const f_T = P_calc - P_system_Pa;

        if (Math.abs(f_T) < tolerance * P_system_Pa || Math.abs(f_T) < 1) { // Relative or absolute tolerance
            y = [ (gamma[0] * x[0] * Psat1) / P_calc, (gamma[1] * x[1] * Psat2) / P_calc ];
            // Normalize y just in case of small numerical errors, though P_calc should be P_system_Pa here
            const sumY = y[0] + y[1];
            if (sumY !== 0) { y[0] /= sumY; y[1] /= sumY;}


            return {
                comp1_feed: x1_feed,
                comp1_equilibrium: y[0],
                T_K: T,
                P_Pa: P_system_Pa, // P_calc should be very close to P_system_Pa
                iterations: iter + 1,
                calculationType: 'bubbleT'
            };
        }

        // Temperature adjustment (simplified - Newton or Brent would be better)
        // Estimate d(f_T)/dT. d(gamma)/dT is complex. Approx d(Psat)/dT.
        // For simplicity, use a bisection-like or scaled step.
        // If f_T > 0 (P_calc > P_system), T is too high, decrease T.
        // If f_T < 0 (P_calc < P_system), T is too low, increase T.
        const T_step_factor = 0.02; // Smaller step factor
        let T_change = -f_T / P_system_Pa * T * T_step_factor; // Heuristic step, scaled by error
        
        // Limit step size
        const max_T_change = T * 0.1; // Max 10% change per iteration
        if (Math.abs(T_change) > max_T_change) {
            T_change = Math.sign(T_change) * max_T_change;
        }
        if (Math.abs(T_change) < 0.01 && f_T !==0) { // Min step if not converged
             T_change = Math.sign(-f_T) * 0.01;
        }


        let T_new = T + T_change;
        if (T_new <= 0) T_new = T * 0.5; // Prevent non-physical temperatures

        T = T_new;
    }

    return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T, P_Pa: P_system_Pa, error: "NRTL: Bubble T max iterations reached", calculationType: 'bubbleT', iterations: maxIter };
}


/**
 * NRTL-based Bubble Point Pressure Calculation.
 * Calculates P_total and y1 given T and x1.
 */
export function calculateBubblePressureNrtl(
    components: CompoundData[], // Uses imported CompoundData
    x1_feed: number,
    T_system_K: number,
    nrtlParams: NrtlInteractionParams,
    initialPressureGuess_Pa: number = 101325, // Typical atmospheric pressure
    maxIter: number = 20,
    tolerance: number = 1e-5
): BubbleDewResult | null { // Uses imported BubbleDewResult
    console.log(`NRTL Bubble P: x1=${x1_feed.toFixed(3)}, T=${T_system_K.toFixed(1)}K, Init_P=${(initialPressureGuess_Pa/1000).toFixed(1)}kPa`);
    // For Bubble P, P is the unknown. Gamma is often assumed independent of P for liquids.
    // P = sum(gamma_i * x_i * Psat_i)

    if (!components[0]?.antoine || !components[1]?.antoine) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "NRTL: Antoine parameters missing for one or both components.", calculationType: 'bubbleP' };
    }

    const Psat1 = calculatePsat_Pa(components[0].antoine, T_system_K, components[0].name);
    const Psat2 = calculatePsat_Pa(components[1].antoine, T_system_K, components[1].name);

    if (isNaN(Psat1) || isNaN(Psat2)) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: `NRTL: Psat calculation failed at T=${T_system_K.toFixed(1)}K`, calculationType: 'bubbleP' };
    }

    const x = [x1_feed, 1 - x1_feed];
    const gamma = calculateNrtlGamma(components, x, T_system_K, nrtlParams);

    if (!gamma) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: `NRTL: Gamma calculation failed at T=${T_system_K.toFixed(1)}K`, calculationType: 'bubbleP' };
    }

    const P_calculated = gamma[0] * x[0] * Psat1 + gamma[1] * x[1] * Psat2;

    if (isNaN(P_calculated) || P_calculated < 0) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: P_calculated, error: `NRTL: Calculated P is NaN or negative (${P_calculated.toFixed(1)} Pa)`, calculationType: 'bubbleP' };
    }
    
    const y1_calculated = (gamma[0] * x[0] * Psat1) / P_calculated;
    // const y2_calculated = (gamma[1] * x[1] * Psat2) / P_calculated; // y2 = 1 - y1

    if (isNaN(y1_calculated)) {
         return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: P_calculated, error: `NRTL: Calculated y1 is NaN`, calculationType: 'bubbleP' };
    }


    return {
        comp1_feed: x1_feed,
        comp1_equilibrium: y1_calculated,
        T_K: T_system_K,
        P_Pa: P_calculated,
        iterations: 1, // Direct calculation for ideal gas phase assumption
        calculationType: 'bubbleP'
    };
}

// Placeholder for Dew Point calculations with NRTL - these are more complex to implement
export function calculateDewTemperatureNrtl(
    _components: CompoundData[], y1_feed: number, _P_total_Pa: number, _nrtlParams: NrtlInteractionParams, _initialTempGuess_K?: number
): BubbleDewResult {
    return { comp1_feed: y1_feed, comp1_equilibrium: NaN, error: "NRTL Dew-T not yet implemented.", calculationType: 'dewT' };
}
export function calculateDewPressureNrtl(
    _components: CompoundData[], y1_feed: number, _T_K: number, _nrtlParams: NrtlInteractionParams, _initialPressureGuess_Pa?: number
): BubbleDewResult {
    return { comp1_feed: y1_feed, comp1_equilibrium: NaN, error: "NRTL Dew-P not yet implemented.", calculationType: 'dewP' };
}
