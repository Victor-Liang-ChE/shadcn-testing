import type { SupabaseClient } from '@supabase/supabase-js';
import type { CompoundData, WilsonPureComponentParams, BubbleDewResult, AntoineParams } from './vle-types';
import { calculatePsat_Pa } from './vle-calculations-unifac'; // Re-use Psat for initial guesses

// --- Constants ---
const R_gas_const_J_molK = 8.31446261815324; // J/mol·K

// --- Interfaces ---
export interface WilsonInteractionParams {
    // a12 and a21 are (lambda_12 - lambda_11) and (lambda_21 - lambda_22) respectively.
    // Units depend on database, typically J/mol or cal/mol. Assuming J/mol for consistency with R_gas_const_J_molK.
    // The table uses A12, A21. We map them to a12, a21 here.
    a12_J_mol: number; 
    a21_J_mol: number; 
    casn1?: string;
    casn2?: string;
}

export async function fetchWilsonInteractionParams(
    supabase: SupabaseClient,
    casn1: string,
    casn2: string
): Promise<WilsonInteractionParams> {
    if (!casn1 || !casn2) {
        console.warn("Wilson a_ij fetch: CAS numbers for both compounds are required. Defaulting a_ij = 0.");
        return { a12_J_mol: 0, a21_J_mol: 0 };
    }
    if (casn1 === casn2) { // For pure component, a_ii related terms effectively cancel or are zero.
        return { a12_J_mol: 0, a21_J_mol: 0 };
    }

    console.log(`Fetching Wilson interaction parameters for pair: casn1='${casn1}', casn2='${casn2}'`);
    // The table is "wilson parameters" and columns are "A12", "A21"
    const query = supabase
        .from('wilson parameters') 
        .select('"A12", "A21", "CASN1", "CASN2"') // Use quoted column names if they contain capitals
        .or(`and(CASN1.eq.${casn1},CASN2.eq.${casn2}),and(CASN1.eq.${casn2},CASN2.eq.${casn1})`)
        .limit(1);

    const { data, error } = await query;

    if (error) {
        console.error(`Supabase Wilson interaction query error: ${error.message}`, error);
        throw new Error(`Supabase Wilson interaction query error: ${error.message}`);
    }
    if (!data || data.length === 0) {
        console.warn(`Wilson interaction parameters (A12, A21) not found for pair ${casn1}/${casn2}. Defaulting to a12=0, a21=0.`);
        return { a12_J_mol: 0, a21_J_mol: 0, casn1, casn2 };
    }
    const params = data[0];
    // Assuming params.A12 and params.A21 are the energy parameters a12 and a21
    return {
        a12_J_mol: typeof params.A12 === 'number' ? params.A12 : 0,
        a21_J_mol: typeof params.A21 === 'number' ? params.A21 : 0,
        casn1: params.CASN1,
        casn2: params.CASN2,
    };
}

export function calculateWilsonGamma(
    components: CompoundData[],
    x: number[], // mole fractions [x1, x2]
    T_K: number,
    interactionParams: WilsonInteractionParams
): [number, number] | null { // Returns [gamma1, gamma2] or null on error
    if (components.length !== 2 || x.length !== 2) {
        console.error("Wilson Gamma calculation currently supports only binary mixtures.");
        return null;
    }
    if (!components[0].wilsonParams || !components[1].wilsonParams) {
        console.error("Wilson pure component parameters (V_L) missing for one or both components.");
        return null;
    }
    if (T_K <= 0) {
        console.error("Wilson Gamma: Temperature must be positive Kelvin.");
        return null;
    }

    const V1_m3mol = components[0].wilsonParams.V_L_m3mol;
    const V2_m3mol = components[1].wilsonParams.V_L_m3mol;

    if (V1_m3mol <= 0 || V2_m3mol <= 0) {
        console.error("Wilson Gamma: Molar volumes must be positive.");
        return null;
    }

    const x1 = x[0];
    const x2 = x[1];

    const { a12_J_mol, a21_J_mol } = interactionParams;

    // Lambda parameters
    // Assuming a12_J_mol and a21_J_mol are (lambda_12 - lambda_11) and (lambda_21 - lambda_22)
    // If they are g_ij - g_ii type parameters (interaction energies)
    const Lambda12 = (V2_m3mol / V1_m3mol) * Math.exp(-a12_J_mol / (R_gas_const_J_molK * T_K));
    const Lambda21 = (V1_m3mol / V2_m3mol) * Math.exp(-a21_J_mol / (R_gas_const_J_molK * T_K));

    const term1_den = x1 + Lambda12 * x2;
    const term2_den = Lambda21 * x1 + x2;

    if (term1_den === 0 || term2_den === 0) {
        console.warn(`Wilson: Division by zero in gamma calculation. T=${T_K.toFixed(1)}, x1=${x1.toFixed(3)}`, {Lambda12, Lambda21, term1_den, term2_den});
        return null; // Avoid division by zero
    }

    const lngamma1 = -Math.log(term1_den) + x2 * (Lambda12 / term1_den - Lambda21 / term2_den);
    const lngamma2 = -Math.log(term2_den) - x1 * (Lambda12 / term1_den - Lambda21 / term2_den); // Note: it's -x1 * (...) or +x1 * (Lambda21/term2_den - Lambda12/term1_den)
    
    const gamma1 = Math.exp(lngamma1);
    const gamma2 = Math.exp(lngamma2);

    if (isNaN(gamma1) || isNaN(gamma2) || !isFinite(gamma1) || !isFinite(gamma2)) {
        console.warn(`Wilson: NaN/Inf gamma. T=${T_K.toFixed(1)}, x1=${x1.toFixed(3)}`, {gamma1, gamma2, Lambda12, Lambda21});
        return null;
    }
    return [gamma1, gamma2];
}


export function calculateBubbleTemperatureWilson(
    components: CompoundData[],
    x1_feed: number,
    P_system_Pa: number,
    wilsonInteractionParams: WilsonInteractionParams,
    initialTempGuess_K: number,
    maxIter: number = 50,
    tolerance: number = 1e-5
): BubbleDewResult | null {
    console.log(`Wilson Bubble T: x1=${x1_feed.toFixed(3)}, P=${(P_system_Pa/1000).toFixed(1)}kPa, Init_T=${initialTempGuess_K.toFixed(1)}K`);
    let T_K = initialTempGuess_K;
    const x = [x1_feed, 1 - x1_feed];

    if (!components[0]?.wilsonParams || !components[1]?.wilsonParams) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: NaN, P_Pa: P_system_Pa, error: "Wilson: Pure component V_L parameters missing.", calculationType: 'bubbleT' };
    }
    if (!components[0]?.antoine || !components[1]?.antoine) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: NaN, P_Pa: P_system_Pa, error: "Wilson: Antoine parameters missing for Psat.", calculationType: 'bubbleT' };
    }

    for (let iter = 0; iter < maxIter; iter++) {
        const gammas = calculateWilsonGamma(components, x, T_K, wilsonInteractionParams);
        if (!gammas) {
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: `Wilson: Gamma calculation failed at T=${T_K.toFixed(1)}K.`, calculationType: 'bubbleT', iterations: iter };
        }
        const Psat1 = calculatePsat_Pa(components[0].antoine, T_K);
        const Psat2 = calculatePsat_Pa(components[1].antoine, T_K);
        if (isNaN(Psat1) || isNaN(Psat2)) {
             return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: `Wilson: Psat calculation failed at T=${T_K.toFixed(1)}K.`, calculationType: 'bubbleT', iterations: iter };
        }

        const P_calc = x[0] * gammas[0] * Psat1 + x[1] * gammas[1] * Psat2;
        const error_func = P_calc - P_system_Pa;

        if (Math.abs(error_func) < tolerance * P_system_Pa || Math.abs(error_func) < 0.1 /*Pa abs tolerance*/) {
            const y1 = (x[0] * gammas[0] * Psat1) / P_calc;
            return { comp1_feed: x1_feed, comp1_equilibrium: y1, T_K: T_K, P_Pa: P_system_Pa, iterations: iter + 1, calculationType: 'bubbleT' };
        }
        
        let T_change = -error_func / P_system_Pa * (T_K * 0.05); 
        const max_T_step = 5.0; 
        if (Math.abs(T_change) > max_T_step) T_change = Math.sign(T_change) * max_T_step;
        if (iter < 3) T_change *= 0.5; 

        let T_K_new = T_K + T_change;
        if (T_K_new <= 0) T_K_new = T_K * 0.9; 
        T_K = T_K_new;
    }
    return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: "Wilson: Bubble T max iterations reached", calculationType: 'bubbleT', iterations: maxIter };
}

export function calculateBubblePressureWilson(
    components: CompoundData[],
    x1_feed: number,
    T_system_K: number,
    wilsonInteractionParams: WilsonInteractionParams,
    initialPressureGuess_Pa: number, 
    maxIter: number = 10, 
    tolerance: number = 1e-5
): BubbleDewResult | null {
    console.log(`Wilson Bubble P: x1=${x1_feed.toFixed(3)}, T=${T_system_K.toFixed(1)}K`);
    const x = [x1_feed, 1 - x1_feed];

    if (!components[0]?.wilsonParams || !components[1]?.wilsonParams) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "Wilson: Pure component V_L parameters missing.", calculationType: 'bubbleP' };
    }
    if (!components[0]?.antoine || !components[1]?.antoine) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "Wilson: Antoine parameters missing for Psat.", calculationType: 'bubbleP' };
    }
    
    const gammas = calculateWilsonGamma(components, x, T_system_K, wilsonInteractionParams);
    if (!gammas) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "Wilson: Gamma calculation failed.", calculationType: 'bubbleP', iterations: 0 };
    }
    const Psat1 = calculatePsat_Pa(components[0].antoine, T_system_K);
    const Psat2 = calculatePsat_Pa(components[1].antoine, T_system_K);
     if (isNaN(Psat1) || isNaN(Psat2)) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "Wilson: Psat calculation failed.", calculationType: 'bubbleP', iterations: 0 };
    }

    const P_calc = x[0] * gammas[0] * Psat1 + x[1] * gammas[1] * Psat2;
    if (isNaN(P_calc) || P_calc <=0) {
         return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: P_calc, error: "Wilson: Calculated pressure is invalid.", calculationType: 'bubbleP', iterations: 0 };
    }
    const y1 = (x[0] * gammas[0] * Psat1) / P_calc;

    return {
        comp1_feed: x1_feed,
        comp1_equilibrium: y1,
        T_K: T_system_K,
        P_Pa: P_calc,
        iterations: 1, 
        calculationType: 'bubbleP'
    };
}
