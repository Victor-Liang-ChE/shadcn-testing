import type { SupabaseClient } from '@supabase/supabase-js';
import type { CompoundData, SrkPureComponentParams, BubbleDewResult, AntoineParams } from './vle-types';
import { calculatePsat_Pa } from './vle-calculations-unifac'; // For initial guesses

// --- Constants ---
const R_J_molK = 8.31446261815324; // Gas constant in J/mol·K

// --- Interfaces ---
export interface SrkInteractionParams {
    k12: number;
    casn1?: string;
    casn2?: string;
}

// Re-use or adapt the cubic solver from PR calculations.
// For SRK EOS: Z^3 - Z^2 + (A - B - B^2)Z - AB = 0
// Mapping to Z^3 + pZ^2 + qZ + r = 0:
// p_coeff = -1.0
// q_coeff = A_mix - B_mix - B_mix * B_mix
// r_coeff = -A_mix * B_mix
function solveCubicRealRoots(p: number, q: number, r: number): number[] {
    // This is the same robust cubic solver from vle-calculations-pr.ts
    // Convert to depressed cubic: y^3 + Ay + B = 0
    // Z = y - p/3
    const A_dep = q - (p * p) / 3.0;
    const B_dep = (2.0 * p * p * p) / 27.0 - (p * q) / 3.0 + r;

    let y_roots: number[] = [];

    const discriminant_term1 = (B_dep / 2.0) * (B_dep / 2.0);
    const discriminant_term2 = (A_dep / 3.0) * (A_dep / 3.0) * (A_dep / 3.0);
    let Delta = discriminant_term1 + discriminant_term2;

    if (Math.abs(Delta) < 1e-12) Delta = 0;

    if (Delta > 0) {
        const sqrtDelta = Math.sqrt(Delta);
        const term1_val = -B_dep / 2.0 + sqrtDelta;
        const term2_val = -B_dep / 2.0 - sqrtDelta;
        const y1 = Math.cbrt(term1_val) + Math.cbrt(term2_val);
        y_roots.push(y1);
    } else if (Delta === 0) {
        if (Math.abs(A_dep) < 1e-9 && Math.abs(B_dep) < 1e-9) {
            y_roots.push(0);
        } else {
            const y1 = -2.0 * Math.cbrt(B_dep / 2.0);
            const y2 = Math.cbrt(B_dep / 2.0);
            y_roots.push(y1);
            if (Math.abs(y1 - y2) > 1e-7) {
                 y_roots.push(y2);
            }
        }
    } else {
        const term_cos_arg_num = (-B_dep / 2.0);
        const term_cos_arg_den = Math.sqrt(-discriminant_term2);
        let cos_arg = term_cos_arg_num / term_cos_arg_den;
        if (cos_arg > 1.0) cos_arg = 1.0;
        if (cos_arg < -1.0) cos_arg = -1.0;
        const angle_third = Math.acos(cos_arg) / 3.0;
        const factor = 2.0 * Math.sqrt(-A_dep / 3.0);
        y_roots.push(factor * Math.cos(angle_third));
        y_roots.push(factor * Math.cos(angle_third - (2.0 * Math.PI) / 3.0));
        y_roots.push(factor * Math.cos(angle_third + (2.0 * Math.PI) / 3.0));
    }

    let Z_roots = y_roots.map(y_root => y_root - p / 3.0);
    Z_roots = Z_roots.filter(z => z > 1e-9);
    if (Z_roots.length > 1) {
        Z_roots.sort((a, b) => a - b);
        Z_roots = Z_roots.filter((z, index, self) => index === 0 || Math.abs(z - self[index-1]) > 1e-7);
    }
    Z_roots.sort((a, b) => a - b);
    return Z_roots;
}

export async function fetchSrkInteractionParams(
    supabase: SupabaseClient,
    casn1: string,
    casn2: string
): Promise<SrkInteractionParams> {
    if (!casn1 || !casn2) {
        console.warn("SRK k_ij fetch: CAS numbers for both compounds are required. Defaulting k_ij = 0.");
        return { k12: 0 };
    }
    if (casn1 === casn2) {
        return { k12: 0 };
    }

    console.log(`Fetching SRK interaction parameters for pair: casn1='${casn1}', casn2='${casn2}'`);
    const query = supabase
        .from('srk paramaters') // Note the space in "srk paramaters"
        .select('"k12", "CASN1", "CASN2"')
        .or(`and(CASN1.eq.${casn1},CASN2.eq.${casn2}),and(CASN1.eq.${casn2},CASN2.eq.${casn1})`)
        .limit(1);

    const { data, error } = await query;
    console.log("Supabase SRK interaction query response - data:", JSON.stringify(data, null, 2));
    console.log("Supabase SRK interaction query response - error:", JSON.stringify(error, null, 2));

    if (error) {
        console.error(`Supabase SRK interaction query raw error object:`, error);
        throw new Error(`Supabase SRK interaction query error: ${error.message}`);
    }
    if (!data || data.length === 0) {
        console.warn(`SRK interaction parameter k12 not found for pair ${casn1}/${casn2}. Defaulting to k12 = 0.`);
        return { k12: 0, casn1, casn2 };
    }
    const params = data[0];
    return {
        k12: typeof params.k12 === 'number' ? params.k12 : 0,
        casn1: params.CASN1,
        casn2: params.CASN2,
    };
}

export function calculateSrkFugacityCoefficients(
    components: CompoundData[],
    x_or_y: number[],
    T_K: number,
    P_Pa: number,
    interactionParams: SrkInteractionParams,
    phase: 'liquid' | 'vapor'
): [number, number] | null {
    if (components.length !== 2 || x_or_y.length !== 2) return null;
    if (!components[0].srkParams || !components[1].srkParams) return null;

    const x1 = x_or_y[0];
    const x2 = x_or_y[1];
    const { k12 } = interactionParams;

    const srkPure = components.map(c => {
        const { Tc_K, Pc_Pa, omega } = c.srkParams!;
        const Tr = T_K / Tc_K;
        const m = 0.480 + 1.574 * omega - 0.176 * omega * omega;
        const alpha = Math.pow(1 + m * (1 - Math.sqrt(Tr)), 2);
        const ac_i = 0.42748 * Math.pow(R_J_molK * Tc_K, 2) / Pc_Pa * alpha;
        const b_i = 0.08664 * R_J_molK * Tc_K / Pc_Pa;
        return { ac_i, b_i, sqrt_ac_i: Math.sqrt(ac_i) }; // Store sqrt_ac_i for mixing rule
    });

    const { ac_i: ac1, b_i: b1, sqrt_ac_i: sqrt_ac1 } = srkPure[0];
    const { ac_i: ac2, b_i: b2, sqrt_ac_i: sqrt_ac2 } = srkPure[1];

    const b_mix = x1 * b1 + x2 * b2;
    const a12_mix_term = sqrt_ac1 * sqrt_ac2 * (1 - k12); // This is (a_ij) for i!=j
    const a_mix = x1 * x1 * ac1 + x2 * x2 * ac2 + 2 * x1 * x2 * a12_mix_term;

    const A_mix_eos = a_mix * P_Pa / Math.pow(R_J_molK * T_K, 2);
    const B_mix_eos = b_mix * P_Pa / (R_J_molK * T_K);

    const p_coeff = -1.0;
    const q_coeff = A_mix_eos - B_mix_eos - B_mix_eos * B_mix_eos;
    const r_coeff = -A_mix_eos * B_mix_eos;

    const Z_roots = solveCubicRealRoots(p_coeff, q_coeff, r_coeff);
    let Z: number;

    if (Z_roots.length === 0) {
        console.warn(`SRK EOS: No positive real roots for Z at T=${T_K.toFixed(1)}, P=${(P_Pa/1000).toFixed(1)}kPa. Phase: ${phase}`, {A_mix_eos, B_mix_eos, x1});
        return null;
    } else if (Z_roots.length === 1) {
        Z = Z_roots[0];
        if (phase === 'liquid' && Z <= B_mix_eos) return null;
    } else {
        if (phase === 'liquid') {
            const liquidCandidates = Z_roots.filter(r => r > B_mix_eos);
            if (liquidCandidates.length > 0) Z = liquidCandidates[0];
            else return null;
        } else {
            Z = Z_roots[Z_roots.length - 1];
        }
    }
    if (Z <= B_mix_eos) return null; // Final check for both phases

    // Fugacity coefficients for SRK
    // ln(phi_i) = (b_i/b_mix)*(Z-1) - ln(Z-B_mix) - (A_mix_eos/B_mix_eos) * ( (2 * sum_j(x_j * sqrt(ac_i*ac_j)*(1-k_ij)))/a_mix - b_i/b_mix ) * ln(1 + B_mix_eos/Z)
    // For component 1: term_sum_j = x1*ac1 + x2*a12_mix_term
    // For component 2: term_sum_j = x2*ac2 + x1*a12_mix_term
    
    const term_A_div_B = (B_mix_eos === 0) ? 0 : A_mix_eos / B_mix_eos; // Avoid division by zero if B_mix_eos is zero
    const term_ln_Z_B = Math.log(1 + B_mix_eos / Z);

    const two_sum_xj_aij_div_amix_1 = (a_mix === 0) ? 0 : (2 * (x1 * ac1 + x2 * a12_mix_term)) / a_mix;
    const two_sum_xj_aij_div_amix_2 = (a_mix === 0) ? 0 : (2 * (x2 * ac2 + x1 * a12_mix_term)) / a_mix;
    
    const ln_phi1 = (b1 / b_mix) * (Z - 1) - Math.log(Z - B_mix_eos) - term_A_div_B * (two_sum_xj_aij_div_amix_1 - (b1 / b_mix)) * term_ln_Z_B;
    const ln_phi2 = (b2 / b_mix) * (Z - 1) - Math.log(Z - B_mix_eos) - term_A_div_B * (two_sum_xj_aij_div_amix_2 - (b2 / b_mix)) * term_ln_Z_B;

    if (isNaN(ln_phi1) || isNaN(ln_phi2) || !isFinite(ln_phi1) || !isFinite(ln_phi2)) {
        console.warn(`SRK EOS: NaN/Inf fugacity coefficient. Z=${Z.toFixed(4)}, T=${T_K.toFixed(1)}, P=${(P_Pa/1000).toFixed(1)}kPa, x1=${x1.toFixed(3)}`, {ln_phi1, ln_phi2, A_mix_eos, B_mix_eos, phase});
        return null;
    }
    return [Math.exp(ln_phi1), Math.exp(ln_phi2)];
}

export function calculateBubbleTemperatureSrk(
    components: CompoundData[],
    x1_feed: number,
    P_system_Pa: number,
    srkInteractionParams: SrkInteractionParams,
    initialTempGuess_K: number,
    maxIter: number = 100,
    tolerance: number = 1e-5
): BubbleDewResult | null {
    console.log(`SRK Bubble T: x1=${x1_feed.toFixed(3)}, P=${(P_system_Pa/1000).toFixed(1)}kPa, Init_T=${initialTempGuess_K.toFixed(1)}K`);
    let T_K = initialTempGuess_K;
    const x = [x1_feed, 1 - x1_feed];
    let y = [...x];

    if (!components[0]?.srkParams || !components[1]?.srkParams) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: NaN, P_Pa: P_system_Pa, error: "SRK: Pure component SRK parameters missing.", calculationType: 'bubbleT' };
    }

    let K_values = [1.0, 1.0];

    for (let iter = 0; iter < maxIter; iter++) {
        const phi_L = calculateSrkFugacityCoefficients(components, x, T_K, P_system_Pa, srkInteractionParams, 'liquid');
        if (!phi_L) {
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: `SRK: Failed to calculate liquid fugacities at T=${T_K.toFixed(1)}K.`, calculationType: 'bubbleT', iterations: iter };
        }

        if (iter === 0) {
            const Psat1_est = components[0].antoine ? calculatePsat_Pa(components[0].antoine, T_K) : P_system_Pa;
            const Psat2_est = components[1].antoine ? calculatePsat_Pa(components[1].antoine, T_K) : P_system_Pa;
            const K1_est = Psat1_est / P_system_Pa;
            const K2_est = Psat2_est / P_system_Pa;
            y = [K1_est * x[0], K2_est * x[1]];
            const sumY_est = y[0] + y[1];
            if (sumY_est > 0) { y[0] /= sumY_est; y[1] /= sumY_est; } else { y = [...x]; }
        } else {
             y = [K_values[0] * x[0], K_values[1] * x[1]];
             const sumY_iter = y[0] + y[1];
             if (sumY_iter > 0) { y[0] /= sumY_iter; y[1] /= sumY_iter; }
        }
        
        const phi_V = calculateSrkFugacityCoefficients(components, y, T_K, P_system_Pa, srkInteractionParams, 'vapor');
        if (!phi_V) {
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: `SRK: Failed to calculate vapor fugacities at T=${T_K.toFixed(1)}K with y=[${y[0].toFixed(3)},${y[1].toFixed(3)}]`, calculationType: 'bubbleT', iterations: iter };
        }

        K_values = [phi_L[0] / phi_V[0], phi_L[1] / phi_V[1]];
        if (isNaN(K_values[0]) || isNaN(K_values[1]) || !isFinite(K_values[0]) || !isFinite(K_values[1])) {
             return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: `SRK: K-value calculation failed (NaN/Inf) at T=${T_K.toFixed(1)}K.`, calculationType: 'bubbleT', iterations: iter };
        }

        const sum_Ki_xi = K_values[0] * x[0] + K_values[1] * x[1];
        const error_func = sum_Ki_xi - 1.0;

        if (Math.abs(error_func) < tolerance) {
            y = [K_values[0] * x[0], K_values[1] * x[1]];
            const sumY_final = y[0] + y[1];
            if (sumY_final > 0) { y[0] /= sumY_final; y[1] /= sumY_final; }
            return { comp1_feed: x1_feed, comp1_equilibrium: y[0], T_K: T_K, P_Pa: P_system_Pa, iterations: iter + 1, calculationType: 'bubbleT' };
        }

        let T_change = error_func * (T_K * 0.1);
        const max_T_step = 10.0;
        if (Math.abs(T_change) > max_T_step) T_change = Math.sign(T_change) * max_T_step;
        if (iter < 5 || Math.abs(error_func) > 0.1) T_change *= 0.5;
        
        let T_K_new = T_K + T_change;
        if (T_K_new <= 0) T_K_new = T_K * 0.9;
        T_K = T_K_new;
    }
    return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: "SRK: Bubble T max iterations reached", calculationType: 'bubbleT', iterations: maxIter };
}

export function calculateBubblePressureSrk(
    components: CompoundData[],
    x1_feed: number,
    T_system_K: number,
    srkInteractionParams: SrkInteractionParams,
    initialPressureGuess_Pa: number,
    maxIter: number = 100,
    tolerance: number = 1e-5
): BubbleDewResult | null {
    console.log(`SRK Bubble P: x1=${x1_feed.toFixed(3)}, T=${T_system_K.toFixed(1)}K, Init_P=${(initialPressureGuess_Pa/1000).toFixed(1)}kPa`);
    let P_Pa = initialPressureGuess_Pa;
    const x = [x1_feed, 1 - x1_feed];
    let y = [...x];

    if (!components[0]?.srkParams || !components[1]?.srkParams) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "SRK: Pure component SRK parameters missing.", calculationType: 'bubbleP' };
    }
    
    let K_values = [1.0, 1.0];

    for (let iter = 0; iter < maxIter; iter++) {
        const phi_L = calculateSrkFugacityCoefficients(components, x, T_system_K, P_Pa, srkInteractionParams, 'liquid');
        if (!phi_L) {
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: P_Pa, error: `SRK: Failed to calculate liquid fugacities at P=${(P_Pa/1000).toFixed(1)}kPa.`, calculationType: 'bubbleP', iterations: iter };
        }

        if (iter === 0) {
             const Psat1_est = components[0].antoine ? calculatePsat_Pa(components[0].antoine, T_system_K) : P_Pa;
             const Psat2_est = components[1].antoine ? calculatePsat_Pa(components[1].antoine, T_system_K) : P_Pa;
             const K1_est = Psat1_est / P_Pa; // Approximation
             const K2_est = Psat2_est / P_Pa; // Approximation
             y = [K1_est * x[0], K2_est * x[1]];
             const sumY_est = y[0] + y[1];
             if (sumY_est > 0) { y[0] /= sumY_est; y[1] /= sumY_est; } else { y = [...x]; }
        } else {
             y = [K_values[0] * x[0], K_values[1] * x[1]];
             const sumY_iter = y[0] + y[1];
             if (sumY_iter > 0) { y[0] /= sumY_iter; y[1] /= sumY_iter; }
        }

        const phi_V = calculateSrkFugacityCoefficients(components, y, T_system_K, P_Pa, srkInteractionParams, 'vapor');
        if (!phi_V) {
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: P_Pa, error: `SRK: Failed to calculate vapor fugacities at P=${(P_Pa/1000).toFixed(1)}kPa with y=[${y[0].toFixed(3)},${y[1].toFixed(3)}]`, calculationType: 'bubbleP', iterations: iter };
        }
        
        K_values = [phi_L[0] / phi_V[0], phi_L[1] / phi_V[1]];
        if (isNaN(K_values[0]) || isNaN(K_values[1]) || !isFinite(K_values[0]) || !isFinite(K_values[1])) {
             return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: P_Pa, error: `SRK: K-value calculation failed (NaN/Inf) at P=${(P_Pa/1000).toFixed(1)}kPa.`, calculationType: 'bubbleP', iterations: iter };
        }

        // Objective: P_calc = P_system. Here, P_system is P_Pa (iterated).
        // P_calc = sum (y_i * P_system) = P_system * sum(y_i)
        // We need sum(y_i) = 1, which means sum(K_i * x_i) = 1.
        // The iteration is on P_Pa.
        // P_new = P_old * sum(K_i * x_i)
        const sum_Ki_xi = K_values[0] * x[0] + K_values[1] * x[1];
        const P_Pa_new = P_Pa * sum_Ki_xi;
        const error_func = P_Pa_new - P_Pa;

        if (Math.abs(error_func) < tolerance * P_Pa || Math.abs(error_func) < 1.0 /*Pa absolute tolerance*/) {
            y = [K_values[0] * x[0], K_values[1] * x[1]]; // Sum should be P_Pa_new / P_Pa
            const sumY_final = y[0] + y[1]; // This is sum_Ki_xi
            // Normalize y based on P_Pa_new
            if (sumY_final > 0) { y[0] /= sumY_final; y[1] /= sumY_final; }

            return { comp1_feed: x1_feed, comp1_equilibrium: y[0], T_K: T_system_K, P_Pa: P_Pa_new, iterations: iter + 1, calculationType: 'bubbleP' };
        }
        
        let P_change_factor = 0.5; // Damping factor
        P_Pa = P_Pa + P_change_factor * error_func;

        if (P_Pa <= 0) {
            P_Pa = initialPressureGuess_Pa * 0.1; // Reset if P goes non-physical
            console.warn(`SRK Bubble P: New P was <=0. Resetting P to ${(P_Pa/1000).toFixed(1)}kPa`);
        }
    }
    return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: P_Pa, error: "SRK: Bubble P max iterations reached", calculationType: 'bubbleP', iterations: maxIter };
}
