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
    casn1_query: string, // CAS number of your "component 0" or "component 1"
    casn2_query: string  // CAS number of your "component 1" or "component 2"
): Promise<UniquacInteractionParams> {
    if (!casn1_query || !casn2_query) {
        console.warn("UNIQUAC Aij fetch: CAS numbers for both compounds are required. Defaulting Aij = 0.");
        return { A12: 0, A21: 0, casn1: casn1_query, casn2: casn2_query };
    }
    if (casn1_query === casn2_query) { // For pure component, Aii = 0
        return { A12: 0, A21: 0, casn1: casn1_query, casn2: casn2_query };
    }

    console.log(`Fetching UNIQUAC interaction parameters for query pair: comp1_cas='${casn1_query}', comp2_cas='${casn2_query}'`);
    const query = supabase
        .from('uniquac parameters') // Table name from user
        .select('"A12", "A21", "CASN1", "CASN2"') // Ensure these column names are exact
        .or(`and(CASN1.eq.${casn1_query},CASN2.eq.${casn2_query}),and(CASN1.eq.${casn2_query},CASN2.eq.${casn1_query})`)
        .limit(1);

    const { data, error } = await query;

    if (error) {
        console.error(`Supabase UNIQUAC interaction query error: ${error.message}`, error);
        throw new Error(`Supabase UNIQUAC interaction query error: ${error.message}`);
    }

    if (!data || data.length === 0) {
        console.warn(`UNIQUAC interaction parameters (A12, A21) not found for query pair ${casn1_query}/${casn2_query}. Defaulting to A12=0, A21=0.`);
        return { A12: 0, A21: 0, casn1: casn1_query, casn2: casn2_query };
    }

    const dbRow = data[0];
    let fetched_A12: number; // This will be A_param_for_casn1_query_TO_casn2_query
    let fetched_A21: number; // This will be A_param_for_casn2_query_TO_casn1_query

    // Check if the first CAS number in the query matches the first CAS number in the DB row
    if (dbRow.CASN1 === casn1_query) {
        // Order matches: casn1_query is DB's CASN1, casn2_query is DB's CASN2
        // So, DB's A12 is for casn1_query -> casn2_query
        // And DB's A21 is for casn2_query -> casn1_query
        fetched_A12 = typeof dbRow.A12 === 'number' ? dbRow.A12 : 0;
        fetched_A21 = typeof dbRow.A21 === 'number' ? dbRow.A21 : 0;
    } else {
        // Order is swapped: casn1_query is DB's CASN2, casn2_query is DB's CASN1
        // So, DB's A12 is for casn2_query -> casn1_query (this should be our A21)
        // And DB's A21 is for casn1_query -> casn2_query (this should be our A12)
        fetched_A12 = typeof dbRow.A21 === 'number' ? dbRow.A21 : 0; // Swapped
        fetched_A21 = typeof dbRow.A12 === 'number' ? dbRow.A12 : 0; // Swapped
    }
    
    console.log(`UNIQUAC_PARAMS_FETCHED: For query ${casn1_query}/${casn2_query}, DB row was ${dbRow.CASN1}/${dbRow.CASN2}. Assigned A12(K)=${fetched_A12}, A21(K)=${fetched_A21}`);

    return {
        A12: fetched_A12,
        A21: fetched_A21,
        casn1: casn1_query, // Return based on query order
        casn2: casn2_query
    };
}

export function calculateUniquacGamma(
    components: CompoundData[],
    x: number[], // mole fractions [x1, x2]
    T_K: number,
    interactionParams: UniquacInteractionParams
): [number, number] | null { // Returns [gamma1, gamma2] or null on error
    console.log("!!!!!!!!!! EXECUTING calculateUniquacGamma (UNIQUAC V4) !!!!!!!!!!");

    if (components.length !== 2 || x.length !== 2) {
        console.error("UNIQUAC Gamma: Supports only binary mixtures.");
        return null;
    }
    if (!components[0].uniquacParams || !components[1].uniquacParams) {
        console.error("UNIQUAC Gamma: Pure component r,q missing.");
        return null;
    }

    const r_params = [components[0].uniquacParams.r, components[1].uniquacParams.r];
    const q_params = [components[0].uniquacParams.q, components[1].uniquacParams.q];
    const x0 = x[0]; // mole fraction of component 0
    const x1 = x[1]; // mole fraction of component 1

    const Z_coord_number_const = 10.0; // UNIQUAC coordination number
    const Z_by_2 = Z_coord_number_const / 2.0;

    const l_val = [
        Z_by_2 * (r_params[0] - q_params[0]) - (r_params[0] - 1.0),
        Z_by_2 * (r_params[1] - q_params[1]) - (r_params[1] - 1.0)
    ];

    const lngamma_C: number[] = [0.0, 0.0];
    const lngamma_R: number[] = [0.0, 0.0]; // Initialize here

    const epsilon = 1e-12;

    // --- Combinatorial Part ---
    // Using the formulation: ln γᵢᶜ = (ln(Φᵢ/xᵢ) + 1 - Φᵢ/xᵢ) - (Z/2)qᵢ * (ln(Φᵢ/Θᵢ) + 1 - Φᵢ/Θᵢ)
    // which DWSIM uses as: (1 - J_i) + ln(J_i) - (Z/2)q_i * (1 - J_i/L_i' + ln(J_i/L_i'))
    // where J_i = r_i / sum(x_k * r_k) = Φ_i / x_i (if x_i !=0)
    // L_i' = q_i / sum(x_k * q_k) = Θ_i / x_i (if x_i !=0)
    // J_i/L_i' = (r_i * sum(x_k * q_k)) / (q_i * sum(x_k * r_k)) = Φ_i / Θ_i

    if (Math.abs(x0) < 1e-9) { // Component 0 (index 0) is infinitely dilute in component 1 (index 1)
        // i=0 (solute), solvent=1
        if (r_params[1] < epsilon || q_params[1] < epsilon) {
            console.warn("UNIQUAC InfDil C0: Solvent r1 or q1 is zero/small.");
            lngamma_C[0] = l_val[0]; // A fallback, though likely not ideal
        } else {
            const J0_limit = r_params[0] / r_params[1]; // Limit of Φ₀/x₀
            const L0_prime_limit = (q_params[1] > epsilon) ? (q_params[0] / q_params[1]) : Infinity; // Limit of Θ₀/x₀ (if q_params[0] is also > epsilon)

            const log_J0_limit = (J0_limit > epsilon) ? Math.log(J0_limit) : -Infinity; // Avoid log(0)

            let term_q_part = 0;
            if (q_params[0] > epsilon) { // If solute q is zero, this whole term is zero
                const J_div_L0_limit = (L0_prime_limit > epsilon && isFinite(L0_prime_limit)) ? (J0_limit / L0_prime_limit) : ( (J0_limit > epsilon && L0_prime_limit === 0) ? Infinity : 0); // Limit of Φ₀/Θ₀
                const log_J_div_L0_limit = (J_div_L0_limit > epsilon && isFinite(J_div_L0_limit)) ? Math.log(J_div_L0_limit) : (J_div_L0_limit === Infinity ? Infinity : -Infinity);
                
                let val_J_div_L0 = J_div_L0_limit;
                if (J_div_L0_limit === Infinity) val_J_div_L0 = Number.MAX_VALUE; // Approximation for (1 - Infinity)
                else if (!isFinite(J_div_L0_limit) && J_div_L0_limit < 0) val_J_div_L0 = -Number.MAX_VALUE; // Should not happen if r,q positive

                term_q_part = Z_by_2 * q_params[0] * (1 - val_J_div_L0 + (isFinite(log_J_div_L0_limit) ? log_J_div_L0_limit : 0) );
                if (log_J_div_L0_limit === Infinity && q_params[0] > epsilon) term_q_part = -Infinity; 
                else if (log_J_div_L0_limit === -Infinity && q_params[0] > epsilon && Math.abs(1 - val_J_div_L0 - 1) < epsilon ) term_q_part = Infinity;
            }
            lngamma_C[0] = (1 - J0_limit) + (isFinite(log_J0_limit) ? log_J0_limit : 0) - term_q_part;
        }
        lngamma_C[1] = 0.0; // Pure component 1
    } else if (Math.abs(x1) < 1e-9) {  // Component 1 (index 1) is infinitely dilute in component 0 (index 0)
        // i=1 (solute), solvent=0
        if (r_params[0] < epsilon || q_params[0] < epsilon) {
            console.warn("UNIQUAC InfDil C1: Solvent r0 or q0 is zero/small.");
            lngamma_C[1] = l_val[1];
        } else {
            const J1_limit = r_params[1] / r_params[0];
            const L1_prime_limit = (q_params[0] > epsilon) ? (q_params[1] / q_params[0]) : Infinity;

            const log_J1_limit = (J1_limit > epsilon) ? Math.log(J1_limit) : -Infinity;
            
            let term_q_part = 0;
            if (q_params[1] > epsilon) {
                const J_div_L1_limit = (L1_prime_limit > epsilon && isFinite(L1_prime_limit)) ? (J1_limit / L1_prime_limit) : ( (J1_limit > epsilon && L1_prime_limit === 0) ? Infinity : 0);
                const log_J_div_L1_limit = (J_div_L1_limit > epsilon && isFinite(J_div_L1_limit)) ? Math.log(J_div_L1_limit) : (J_div_L1_limit === Infinity ? Infinity : -Infinity);

                let val_J_div_L1 = J_div_L1_limit;
                if (J_div_L1_limit === Infinity) val_J_div_L1 = Number.MAX_VALUE;
                else if (!isFinite(J_div_L1_limit) && J_div_L1_limit < 0) val_J_div_L1 = -Number.MAX_VALUE;
                
                term_q_part = Z_by_2 * q_params[1] * (1 - val_J_div_L1 + (isFinite(log_J_div_L1_limit) ? log_J_div_L1_limit : 0) );
                if (log_J_div_L1_limit === Infinity && q_params[1] > epsilon) term_q_part = -Infinity;
                else if (log_J_div_L1_limit === -Infinity && q_params[1] > epsilon && Math.abs(1 - val_J_div_L1 - 1) < epsilon) term_q_part = Infinity;
            }
            lngamma_C[1] = (1 - J1_limit) + (isFinite(log_J1_limit) ? log_J1_limit : 0) - term_q_part;
        }
        lngamma_C[0] = 0.0;  // Pure component 0
    } else {   // Mixture
        const sum_xr = x0 * r_params[0] + x1 * r_params[1];
        const sum_xq = x0 * q_params[0] + x1 * q_params[1];

        if (sum_xr < epsilon || sum_xq < epsilon) {
            console.warn("UNIQUAC Combinatorial Mixture: sum_xr or sum_xq is near zero.");
            return null;
        }

        const Phi0_mix = (x0 * r_params[0]) / sum_xr;
        const Phi1_mix = (x1 * r_params[1]) / sum_xr;
        const Theta0_mix = (x0 * q_params[0]) / sum_xq;
        const Theta1_mix = (x1 * q_params[1]) / sum_xq;
        
        const sum_xl_mix = x0 * l_val[0] + x1 * l_val[1];

        // Component 0
        if (Phi0_mix <= epsilon || x0 <= epsilon || Theta0_mix <= epsilon || Phi0_mix/x0 <= epsilon || Theta0_mix/Phi0_mix <= epsilon ) {
            console.warn(`UNIQUAC Mixture C0: Non-positive or very small arg for log/div. x0=${x0.toExponential(2)}, Phi0=${Phi0_mix.toExponential(2)}, Theta0=${Theta0_mix.toExponential(2)}`);
            lngamma_C[0] = 0; // Fallback
        } else {
            lngamma_C[0] = Math.log(Phi0_mix / x0) + 
                           Z_by_2 * q_params[0] * Math.log(Theta0_mix / Phi0_mix) + 
                           l_val[0] - (Phi0_mix / x0) * sum_xl_mix;
        }

        // Component 1
        if (Phi1_mix <= epsilon || x1 <= epsilon || Theta1_mix <= epsilon || Phi1_mix/x1 <= epsilon || Theta1_mix/Phi1_mix <= epsilon) {
            console.warn(`UNIQUAC Mixture C1: Non-positive or very small arg for log/div. x1=${x1.toExponential(2)}, Phi1=${Phi1_mix.toExponential(2)}, Theta1=${Theta1_mix.toExponential(2)}`);
            lngamma_C[1] = 0; // Fallback
        } else {
            lngamma_C[1] = Math.log(Phi1_mix / x1) + 
                           Z_by_2 * q_params[1] * Math.log(Theta1_mix / Phi1_mix) + 
                           l_val[1] - (Phi1_mix / x1) * sum_xl_mix;
        }
    }
    
    console.log(`UNIQUAC_DETAILS: ln_gamma_C = [${lngamma_C[0]?.toFixed(4)}, ${lngamma_C[1]?.toFixed(4)}]`);

    // --- Residual Part ---
    const { A12, A21 } = interactionParams;
    const tau12 = Math.exp(-A12 / T_K); 
    const tau21 = Math.exp(-A21 / T_K);

    const sum_xq_res_eff = (x0 * q_params[0] + x1 * q_params[1]);
    let Theta0_res = 0.0, Theta1_res = 0.0;

    if (sum_xq_res_eff > epsilon) {
        Theta0_res = (x0 * q_params[0]) / sum_xq_res_eff;
        Theta1_res = (x1 * q_params[1]) / sum_xq_res_eff;
    } else if (x0 > 0.5) { // Fallback if sum_xq_res_eff is zero (e.g. all q are zero)
        Theta0_res = 1.0; Theta1_res = 0.0;
    } else if (x1 > 0.5) {
        Theta0_res = 0.0; Theta1_res = 1.0;
    }

    const s1_res = Theta0_res * 1.0 + Theta1_res * tau21; // Sum_j Theta_j * tau_j0 (tau00=1)
    const s2_res = Theta0_res * tau12 + Theta1_res * 1.0; // Sum_j Theta_j * tau_j1 (tau11=1)
    
    let term_s1_log_c0 = 0;
    if (s1_res > epsilon) term_s1_log_c0 = Math.log(s1_res);
    else if (s1_res < -epsilon) { console.warn("UNIQUAC Residual C0: s1_res is negative for log."); return null; }
    // if s1_res is near zero and positive, log will be large negative.
    
    let term0_r_frac1 = 0;
    if (s1_res > epsilon) term0_r_frac1 = Theta0_res / s1_res;
    else if (Theta0_res > epsilon && Math.abs(s1_res) < epsilon) { console.warn("UNIQUAC Residual C0: Division by near-zero s1_res for term0_r_frac1."); return null; }
    
    let term0_r_frac2 = 0;
    if (s2_res > epsilon) term0_r_frac2 = (Theta1_res * tau12) / s2_res;
    else if (Theta1_res * tau12 > epsilon && Math.abs(s2_res) < epsilon) { console.warn("UNIQUAC Residual C0: Division by near-zero s2_res for term0_r_frac2."); return null; }
    
    lngamma_R[0] = q_params[0] * (1.0 - term_s1_log_c0 - term0_r_frac1 - term0_r_frac2);

    let term_s2_log_c1 = 0;
    if (s2_res > epsilon) term_s2_log_c1 = Math.log(s2_res);
    else if (s2_res < -epsilon) { console.warn("UNIQUAC Residual C1: s2_res is negative for log."); return null; }

    let term1_r1_frac1 = 0;
    if (s1_res > epsilon) term1_r1_frac1 = (Theta0_res * tau21) / s1_res;
    else if (Theta0_res * tau21 > epsilon && Math.abs(s1_res) < epsilon) { console.warn("UNIQUAC Residual C1: Division by near-zero s1_res for term1_r1_frac1."); return null; }
    
    let term1_r1_frac2 = 0;
    if (s2_res > epsilon) term1_r1_frac2 = Theta1_res / s2_res;
    else if (Theta1_res > epsilon && Math.abs(s2_res) < epsilon) { console.warn("UNIQUAC Residual C1: Division by near-zero s2_res for term1_r1_frac2."); return null; }
    
    lngamma_R[1] = q_params[1] * (1.0 - term_s2_log_c1 - term1_r1_frac1 - term1_r1_frac2);
    
    const gamma1_val = Math.exp(lngamma_C[0] + lngamma_R[0]);
    const gamma2_val = Math.exp(lngamma_C[1] + lngamma_R[1]);

    if (isNaN(gamma1_val) || isNaN(gamma2_val) || !isFinite(gamma1_val) || !isFinite(gamma2_val)) {
        console.warn(`UNIQUAC: NaN/Inf gamma produced. T=${T_K.toFixed(1)}, x0=${x0.toFixed(3)}`, {x_input: x, lngamma_C, lngamma_R, tau12, tau21, s1_res, s2_res, Theta0_res, Theta1_res});
        return null;
    }
    return [gamma1_val, gamma2_val];
}

// Re-use or define numericalDerivative if not already globally available
function numericalDerivative(
    func: (t: number) => number,
    t: number,
    h_relative_step: number = 1e-4 // Relative step for h
): number {
    const h = Math.max(1e-5, Math.abs(t * h_relative_step)); // Absolute minimum h, else relative
    const f_plus_h = func(t + h);
    const f_minus_h = func(t - h);

    if (isNaN(f_plus_h) || isNaN(f_minus_h)) {
        return NaN;
    }
    if (Math.abs(2 * h) < 1e-12) { // Avoid division by zero
        return NaN;
    }
    return (f_plus_h - f_minus_h) / (2 * h);
}


export function calculateBubbleTemperatureUniquac(
    components: CompoundData[],
    x1_feed: number,
    P_system_Pa: number,
    uniquacInteractionParams: UniquacInteractionParams,
    initialTempGuess_K: number,
    maxIter: number = 100, // Increased maxIter from 50
    tolerance_f: number = 1e-6, // Tolerance for the objective function |f(T)|
    tolerance_T_step: number = 1e-4 // Tolerance for |T_new - T_old|
): BubbleDewResult | null {
    console.log(`UNIQUAC Bubble T (Newton): x1=${x1_feed.toFixed(3)}, P=${(P_system_Pa / 1000).toFixed(1)}kPa, Init_T=${initialTempGuess_K.toFixed(1)}K`);
    let T_K = initialTempGuess_K;
    const x = [x1_feed, 1 - x1_feed];

    if (!components[0]?.uniquacParams || !components[1]?.uniquacParams) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: NaN, P_Pa: P_system_Pa, error: "UNIQUAC: Pure component r,q parameters missing.", calculationType: 'bubbleT' };
    }
    if (!components[0]?.antoine || !components[1]?.antoine) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: NaN, P_Pa: P_system_Pa, error: "UNIQUAC: Antoine parameters missing for Psat.", calculationType: 'bubbleT' };
    }

    // Objective function: f(T) = P_calc - P_system_Pa
    const objectiveFunction = (T_current_K: number): number => {
        const gammas = calculateUniquacGamma(components, x, T_current_K, uniquacInteractionParams);
        if (!gammas) return NaN;

        const Psat1 = calculatePsat_Pa(components[0].antoine!, T_current_K, components[0].name);
        const Psat2 = calculatePsat_Pa(components[1].antoine!, T_current_K, components[1].name);
        if (isNaN(Psat1) || isNaN(Psat2) || Psat1 <=0 || Psat2 <=0) return NaN;

        const P_calc = x[0] * gammas[0] * Psat1 + x[1] * gammas[1] * Psat2;
        return P_calc - P_system_Pa;
    };

    for (let iter = 0; iter < maxIter; iter++) {
        const fT = objectiveFunction(T_K);

        if (isNaN(fT)) {
            // console.error(`UNIQUAC BubbleT: Objective function returned NaN at T=${T_K.toFixed(2)}K, iter=${iter}`);
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, iterations: iter, error: "UNIQUAC BubbleT: Objective function NaN", calculationType: 'bubbleT' };
        }

        if (Math.abs(fT) < tolerance_f * P_system_Pa || Math.abs(fT) < 0.1 /*Pa abs tol*/) { // Relative or small absolute
            const gammas_conv = calculateUniquacGamma(components, x, T_K, uniquacInteractionParams);
            const Psat1_conv = calculatePsat_Pa(components[0].antoine!, T_K, components[0].name);
            const Psat2_conv = calculatePsat_Pa(components[1].antoine!, T_K, components[1].name);
            if (!gammas_conv || isNaN(Psat1_conv) || isNaN(Psat2_conv) || Psat1_conv <=0 || Psat2_conv <=0) {
                 return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, iterations: iter, error: "UNIQUAC BubbleT: Final prop calc failed", calculationType: 'bubbleT' };
            }
            const P_final = x[0] * gammas_conv[0] * Psat1_conv + x[1] * gammas_conv[1] * Psat2_conv;
            if (P_final <= 0) { // Avoid division by zero or negative pressure
                return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, iterations: iter, error: "UNIQUAC BubbleT: Final P_calc non-positive", calculationType: 'bubbleT' };
            }
            const y1 = (x[0] * gammas_conv[0] * Psat1_conv) / P_final;
            return { comp1_feed: x1_feed, comp1_equilibrium: Math.min(Math.max(y1,0),1), T_K: T_K, P_Pa: P_system_Pa, iterations: iter + 1, calculationType: 'bubbleT' };
        }

        const dfdT = numericalDerivative(objectiveFunction, T_K);

        if (isNaN(dfdT) || Math.abs(dfdT) < 1e-9) {
            // console.error(`UNIQUAC BubbleT: Derivative is NaN or too small at T=${T_K.toFixed(2)}K, dfdT=${dfdT?.toExponential(2)}, iter=${iter}`);
             // Try a small fixed step if derivative is bad
            T_K += (fT > 0 ? -0.5 : 0.5) * (iter > 10 ? 0.1 : 1.0) ; // Small kick, reduce kick after some iterations
            if (T_K <= 0) T_K = initialTempGuess_K * (0.9 + Math.random()*0.2); // Reset if T goes wild
            continue;
        }

        let deltaT = -fT / dfdT;
        const max_deltaT_abs = T_K * 0.15; // Max 15% change of current T
        if (Math.abs(deltaT) > max_deltaT_abs) {
            deltaT = Math.sign(deltaT) * max_deltaT_abs;
        }
         if (iter < 5) { // Stronger damping for initial iterations
            deltaT *= 0.6;
        }

        const T_K_new = T_K + deltaT;
        // console.debug(`UNIQUAC Bubble T iter ${iter}: T_old=${T_K.toFixed(2)}, f(T)=${fT.toExponential(3)}, dfdT=${dfdT.toExponential(3)}, deltaT=${deltaT.toFixed(3)}, T_new=${T_K_new.toFixed(2)}`);

        if (T_K_new <= 0) {
            // console.warn(`UNIQUAC Bubble T: Newton step resulted in T_new <= 0 (${T_K_new.toFixed(2)}). Halving current T_change.`);
            T_K = T_K + deltaT / 2; // Try smaller step from T_K, not T_K_new
            if (T_K <= 0) T_K = initialTempGuess_K * (0.8 + Math.random()*0.4); // Fallback
            continue;
        }
        
        if (Math.abs(deltaT) < tolerance_T_step && Math.abs(deltaT) < tolerance_T_step * Math.abs(T_K_new) && iter > 0 ) {
            const gammas_conv = calculateUniquacGamma(components, x, T_K_new, uniquacInteractionParams);
            const Psat1_conv = calculatePsat_Pa(components[0].antoine!, T_K_new, components[0].name);
            const Psat2_conv = calculatePsat_Pa(components[1].antoine!, T_K_new, components[1].name);
             if (!gammas_conv || isNaN(Psat1_conv) || isNaN(Psat2_conv) || Psat1_conv <=0 || Psat2_conv <=0) {
                 return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K_new, P_Pa: P_system_Pa, iterations: iter + 1, error: "UNIQUAC BubbleT: Final prop calc failed on T_step convergence", calculationType: 'bubbleT' };
             }
            const P_final = x[0] * gammas_conv[0] * Psat1_conv + x[1] * gammas_conv[1] * Psat2_conv;
            if (P_final <= 0) {
                return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K_new, P_Pa: P_system_Pa, iterations: iter + 1, error: "UNIQUAC BubbleT: Final P_calc non-positive on T_step convergence", calculationType: 'bubbleT' };
            }
            const y1 = (x[0] * gammas_conv[0] * Psat1_conv) / P_final;
            return { comp1_feed: x1_feed, comp1_equilibrium: Math.min(Math.max(y1,0),1), T_K: T_K_new, P_Pa: P_system_Pa, iterations: iter + 1, calculationType: 'bubbleT' };
        }
        T_K = T_K_new;
    }

    // console.warn(`UNIQUAC Bubble T: Failed to converge for x1=${x1_feed.toFixed(4)} after ${maxIter} iterations. Last T=${T_K.toFixed(2)}K, f(T)=${objectiveFunction(T_K)?.toExponential(2)}`);
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
