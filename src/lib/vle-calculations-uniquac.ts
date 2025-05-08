import type { SupabaseClient } from '@supabase/supabase-js';
import type { CompoundData, UniquacPureComponentParams, BubbleDewResult, AntoineParams } from './vle-types';
import { calculatePsat_Pa } from './vle-calculations-unifac'; // Re-use Psat for initial guesses

// --- Constants ---
const R_gas_const = 8.31446261815324; // J/mol·K (used if Aij were energy, but here Aij are in K)
const Z_coord_number = 10; // UNIQUAC coordination number

// --- Interfaces ---
export interface UniquacInteractionParams {
    A12: number; // Interaction parameter a_12 (units of K)
    A21: number; // Interaction parameter a_21 (units of K)
    casn1?: string;
    casn2?: string;
}

export async function fetchUniquacInteractionParams(
    supabase: SupabaseClient,
    casn1: string,
    casn2: string
): Promise<UniquacInteractionParams> {
    if (!casn1 || !casn2) {
        console.warn("UNIQUAC Aij fetch: CAS numbers for both compounds are required. Defaulting Aij = 0.");
        return { A12: 0, A21: 0 };
    }
    if (casn1 === casn2) { // For pure component, Aii = 0
        return { A12: 0, A21: 0 };
    }

    console.log(`Fetching UNIQUAC interaction parameters for pair: casn1='${casn1}', casn2='${casn2}'`);
    const query = supabase
        .from('uniquac parameters') // Table name from user
        .select('"A12", "A21", "CASN1", "CASN2"')
        .or(`and(CASN1.eq.${casn1},CASN2.eq.${casn2}),and(CASN1.eq.${casn2},CASN2.eq.${casn1})`)
        .limit(1);

    const { data, error } = await query;

    if (error) {
        console.error(`Supabase UNIQUAC interaction query error: ${error.message}`, error);
        throw new Error(`Supabase UNIQUAC interaction query error: ${error.message}`);
    }
    if (!data || data.length === 0) {
        console.warn(`UNIQUAC interaction parameters (A12, A21) not found for pair ${casn1}/${casn2}. Defaulting to A12=0, A21=0.`);
        return { A12: 0, A21: 0, casn1, casn2 };
    }
    const params = data[0];
    return {
        A12: typeof params.A12 === 'number' ? params.A12 : 0,
        A21: typeof params.A21 === 'number' ? params.A21 : 0,
        casn1: params.CASN1,
        casn2: params.CASN2,
    };
}

export function calculateUniquacGamma(
    components: CompoundData[],
    x: number[], // mole fractions [x1, x2]
    T_K: number,
    interactionParams: UniquacInteractionParams
): [number, number] | null { // Returns [gamma1, gamma2] or null on error
    if (components.length !== 2 || x.length !== 2) {
        console.error("UNIQUAC Gamma calculation currently supports only binary mixtures.");
        return null;
    }
    if (!components[0].uniquacParams || !components[1].uniquacParams) {
        console.error("UNIQUAC pure component parameters (r, q) missing for one or both components.");
        return null;
    }

    const r1 = components[0].uniquacParams.r;
    const q1 = components[0].uniquacParams.q;
    const r2 = components[1].uniquacParams.r;
    const q2 = components[1].uniquacParams.q;

    const x1 = x[0];
    const x2 = x[1];

    // Combinatorial part
    const Phi1 = (x1 * r1) / (x1 * r1 + x2 * r2);
    const Phi2 = (x2 * r2) / (x1 * r1 + x2 * r2);
    const Theta1 = (x1 * q1) / (x1 * q1 + x2 * q2);
    const Theta2 = (x2 * q2) / (x1 * q1 + x2 * q2);

    const l1 = (Z_coord_number / 2) * (r1 - q1) - (r1 - 1);
    const l2 = (Z_coord_number / 2) * (r2 - q2) - (r2 - 1);

    const lngamma1_C = Math.log(Phi1 / x1) + (Z_coord_number / 2) * q1 * Math.log(Theta1 / Phi1) + l1 - (Phi1 / x1) * (x1 * l1 + x2 * l2);
    const lngamma2_C = Math.log(Phi2 / x2) + (Z_coord_number / 2) * q2 * Math.log(Theta2 / Phi2) + l2 - (Phi2 / x2) * (x1 * l1 + x2 * l2);

    // Residual part
    const { A12, A21 } = interactionParams; // These are a_ij in K
    const tau12 = Math.exp(-A12 / T_K); // A12 is a_12
    const tau21 = Math.exp(-A21 / T_K); // A21 is a_21
    // tau_ii = exp(-Aii/T_K) = exp(0) = 1

    const sum_Theta_tau_1 = Theta1 + Theta2 * tau21; // sum_j Theta_j * tau_j1 (tau_11 = 1)
    const sum_Theta_tau_2 = Theta1 * tau12 + Theta2; // sum_j Theta_j * tau_j2 (tau_22 = 1)

    const lngamma1_R = q1 * (1 - Math.log(sum_Theta_tau_1) - (Theta1 / sum_Theta_tau_1) - (Theta2 * tau12 / sum_Theta_tau_2));
    const lngamma2_R = q2 * (1 - Math.log(sum_Theta_tau_2) - (Theta1 * tau21 / sum_Theta_tau_1) - (Theta2 / sum_Theta_tau_2));
    
    const gamma1 = Math.exp(lngamma1_C + lngamma1_R);
    const gamma2 = Math.exp(lngamma2_C + lngamma2_R);

    if (isNaN(gamma1) || isNaN(gamma2) || !isFinite(gamma1) || !isFinite(gamma2)) {
        console.warn(`UNIQUAC: NaN/Inf gamma. T=${T_K.toFixed(1)}, x1=${x1.toFixed(3)}`, {gamma1, gamma2, tau12, tau21});
        return null;
    }
    return [gamma1, gamma2];
}


export function calculateBubbleTemperatureUniquac(
    components: CompoundData[],
    x1_feed: number,
    P_system_Pa: number,
    uniquacInteractionParams: UniquacInteractionParams,
    initialTempGuess_K: number,
    maxIter: number = 50,
    tolerance: number = 1e-5
): BubbleDewResult | null {
    console.log(`UNIQUAC Bubble T: x1=${x1_feed.toFixed(3)}, P=${(P_system_Pa/1000).toFixed(1)}kPa, Init_T=${initialTempGuess_K.toFixed(1)}K`);
    let T_K = initialTempGuess_K;
    const x = [x1_feed, 1 - x1_feed];

    if (!components[0]?.uniquacParams || !components[1]?.uniquacParams) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: NaN, P_Pa: P_system_Pa, error: "UNIQUAC: Pure component r,q parameters missing.", calculationType: 'bubbleT' };
    }
    if (!components[0]?.antoine || !components[1]?.antoine) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: NaN, P_Pa: P_system_Pa, error: "UNIQUAC: Antoine parameters missing for Psat.", calculationType: 'bubbleT' };
    }

    for (let iter = 0; iter < maxIter; iter++) {
        const gammas = calculateUniquacGamma(components, x, T_K, uniquacInteractionParams);
        if (!gammas) {
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: `UNIQUAC: Gamma calculation failed at T=${T_K.toFixed(1)}K.`, calculationType: 'bubbleT', iterations: iter };
        }
        const Psat1 = calculatePsat_Pa(components[0].antoine, T_K);
        const Psat2 = calculatePsat_Pa(components[1].antoine, T_K);
        if (isNaN(Psat1) || isNaN(Psat2)) {
             return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: `UNIQUAC: Psat calculation failed at T=${T_K.toFixed(1)}K.`, calculationType: 'bubbleT', iterations: iter };
        }

        const P_calc = x[0] * gammas[0] * Psat1 + x[1] * gammas[1] * Psat2;
        const error_func = P_calc - P_system_Pa;

        if (Math.abs(error_func) < tolerance * P_system_Pa || Math.abs(error_func) < 0.1 /*Pa abs tolerance*/) {
            const y1 = (x[0] * gammas[0] * Psat1) / P_calc;
            return { comp1_feed: x1_feed, comp1_equilibrium: y1, T_K: T_K, P_Pa: P_system_Pa, iterations: iter + 1, calculationType: 'bubbleT' };
        }

        // Simple temperature adjustment - can be improved (e.g., Newton or Brent)
        // If P_calc > P_system, T is too high, decrease T.
        // If P_calc < P_system, T is too low, increase T.
        let T_change = -error_func / P_system_Pa * (T_K * 0.05); // Proportional step, scaled by current T
        const max_T_step = 5.0; // Max 5K change
        if (Math.abs(T_change) > max_T_step) T_change = Math.sign(T_change) * max_T_step;
        if (iter < 3) T_change *= 0.5; // Damp early steps

        let T_K_new = T_K + T_change;
        if (T_K_new <= 0) T_K_new = T_K * 0.9; // Prevent non-physical T
        T_K = T_K_new;
    }
    return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: "UNIQUAC: Bubble T max iterations reached", calculationType: 'bubbleT', iterations: maxIter };
}

export function calculateBubblePressureUniquac(
    components: CompoundData[],
    x1_feed: number,
    T_system_K: number,
    uniquacInteractionParams: UniquacInteractionParams,
    initialPressureGuess_Pa: number, // Not directly used in this formulation but good for consistency
    maxIter: number = 10, // Usually converges fast
    tolerance: number = 1e-5 // Relative tolerance for P
): BubbleDewResult | null {
    console.log(`UNIQUAC Bubble P: x1=${x1_feed.toFixed(3)}, T=${T_system_K.toFixed(1)}K`);
    const x = [x1_feed, 1 - x1_feed];

    if (!components[0]?.uniquacParams || !components[1]?.uniquacParams) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "UNIQUAC: Pure component r,q parameters missing.", calculationType: 'bubbleP' };
    }
    if (!components[0]?.antoine || !components[1]?.antoine) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "UNIQUAC: Antoine parameters missing for Psat.", calculationType: 'bubbleP' };
    }
    
    // Bubble pressure is a direct calculation for activity coefficient models if T is fixed
    const gammas = calculateUniquacGamma(components, x, T_system_K, uniquacInteractionParams);
    if (!gammas) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "UNIQUAC: Gamma calculation failed.", calculationType: 'bubbleP', iterations: 0 };
    }
    const Psat1 = calculatePsat_Pa(components[0].antoine, T_system_K);
    const Psat2 = calculatePsat_Pa(components[1].antoine, T_system_K);
     if (isNaN(Psat1) || isNaN(Psat2)) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "UNIQUAC: Psat calculation failed.", calculationType: 'bubbleP', iterations: 0 };
    }

    const P_calc = x[0] * gammas[0] * Psat1 + x[1] * gammas[1] * Psat2;
    if (isNaN(P_calc) || P_calc <=0) {
         return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: P_calc, error: "UNIQUAC: Calculated pressure is invalid.", calculationType: 'bubbleP', iterations: 0 };
    }
    const y1 = (x[0] * gammas[0] * Psat1) / P_calc;

    return {
        comp1_feed: x1_feed,
        comp1_equilibrium: y1,
        T_K: T_system_K,
        P_Pa: P_calc,
        iterations: 1, // Direct calculation
        calculationType: 'bubbleP'
    };
}
