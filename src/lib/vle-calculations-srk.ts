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

// Helper function for numerical derivative (central difference)
function numericalDerivative(
    func: (t: number) => number,
    t: number,
    h: number = 1e-3 // for temperature, can be tuned
): number {
    const f_plus_h = func(t + h);
    const f_minus_h = func(t - h);

    if (isNaN(f_plus_h) || isNaN(f_minus_h)) {
        console.warn(`Numerical derivative: NaN encountered at T=${t.toFixed(2)}`);
        return NaN;
    }
    if (Math.abs(2 * h) < 1e-9) { // Avoid division by zero if h is too small
        console.warn(`Numerical derivative: Step size h is too small or zero.`);
        return NaN;
    }
    return (f_plus_h - f_minus_h) / (2 * h);
}

export function calculateBubbleTemperatureSrk(
    components: CompoundData[],
    x1_feed: number,
    P_system_Pa: number,
    srkInteractionParams: SrkInteractionParams,
    initialTempGuess_K: number,
    maxIter: number = 100,
    tolerance_f: number = 1e-6, // for the objective function |f(T)|
    tolerance_T_step: number = 1e-4 // for |T_new - T_old|
): BubbleDewResult | null {
    console.log(`SRK Bubble T (Newton): x1=${x1_feed.toFixed(3)}, P=${(P_system_Pa / 1000).toFixed(1)}kPa, Init_T=${initialTempGuess_K.toFixed(1)}K`);
    let T_K = initialTempGuess_K;
    const x = [x1_feed, 1 - x1_feed]; // Liquid phase composition is fixed for bubble point

    if (!components[0]?.srkParams || !components[1]?.srkParams) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: NaN, P_Pa: P_system_Pa, error: "SRK: Pure component SRK parameters missing.", calculationType: 'bubbleT' };
    }

    let K_values_iter: number[] = [1.0, 1.0]; // K-values from current iteration
    let y_iter: number[] = [...x];       // y-values from current iteration

    // Objective function: f(T) = sum(K_i * x_i) - 1
    const objectiveFunction = (T_current_K: number): number => {
        // Calculate K-values at T_current_K, P_system_Pa, and for liquid phase x
        const phi_L = calculateSrkFugacityCoefficients(components, x, T_current_K, P_system_Pa, srkInteractionParams, 'liquid');
        if (!phi_L) {
            console.warn(`SRK Bubble T ObjFunc: Failed to get phi_L at T=${T_current_K.toFixed(2)}K`);
            return NaN; // Indicate failure
        }

        // Estimate y based on current K_values_iter
        let y_for_phi_V = [K_values_iter[0] * x[0], K_values_iter[1] * x[1]];
        const sumY_phi_V = y_for_phi_V[0] + y_for_phi_V[1];
        if (sumY_phi_V > 1e-9) {
            y_for_phi_V = [y_for_phi_V[0] / sumY_phi_V, y_for_phi_V[1] / sumY_phi_V];
        } else {
            y_for_phi_V = [...x]; // Fallback if K-values are problematic
        }
        
        const phi_V = calculateSrkFugacityCoefficients(components, y_for_phi_V, T_current_K, P_system_Pa, srkInteractionParams, 'vapor');
        if (!phi_V) {
            console.warn(`SRK Bubble T ObjFunc: Failed to get phi_V at T=${T_current_K.toFixed(2)}K for y=[${y_for_phi_V[0].toFixed(3)},${y_for_phi_V[1].toFixed(3)}]`);
            return NaN; // Indicate failure
        }

        const K_at_T = [phi_L[0] / phi_V[0], phi_L[1] / phi_V[1]];
        if (K_at_T.some(k => isNaN(k) || !isFinite(k))) {
            console.warn(`SRK Bubble T ObjFunc: Invalid K-values at T=${T_current_K.toFixed(2)}K`);
            return NaN;
        }

        (objectiveFunction as any)._last_K_values = K_at_T;
        const sum_Kx_for_y = K_at_T[0]*x[0] + K_at_T[1]*x[1];
        if (sum_Kx_for_y > 1e-9) {
            (objectiveFunction as any)._last_y_values = [K_at_T[0]*x[0]/sum_Kx_for_y, K_at_T[1]*x[1]/sum_Kx_for_y];
        } else {
             (objectiveFunction as any)._last_y_values = [...x]; // Fallback
        }

        return K_at_T[0] * x[0] + K_at_T[1] * x[1] - 1.0;
    };

    // Initialize K_values_iter and y_iter for the first call to objectiveFunction
    if (components[0].antoine && components[1].antoine) {
        const Psat1_init = calculatePsat_Pa(components[0].antoine, T_K); // Removed name argument
        const Psat2_init = calculatePsat_Pa(components[1].antoine, T_K); // Removed name argument
        if (!isNaN(Psat1_init) && !isNaN(Psat2_init) && Psat1_init > 0 && Psat2_init > 0) {
            K_values_iter = [Psat1_init / P_system_Pa, Psat2_init / P_system_Pa];
            const sumY_init = K_values_iter[0]*x[0] + K_values_iter[1]*x[1];
            if(sumY_init > 1e-9) {
                y_iter = [K_values_iter[0]*x[0]/sumY_init, K_values_iter[1]*x[1]/sumY_init];
            }
        }
    }

    for (let iter = 0; iter < maxIter; iter++) {
        const fT = objectiveFunction(T_K);

        if ((objectiveFunction as any)._last_K_values) {
            K_values_iter = (objectiveFunction as any)._last_K_values;
        }
        if ((objectiveFunction as any)._last_y_values) {
            y_iter = (objectiveFunction as any)._last_y_values;
        }

        if (isNaN(fT)) {
            // console.error(`SRK Bubble T: Objective function returned NaN at T=${T_K.toFixed(2)}K, iter=${iter}`);
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, iterations: iter, error: "SRK BubbleT: Objective function NaN", calculationType: 'bubbleT' };
        }

        if (Math.abs(fT) < tolerance_f) {
            // console.log(`SRK Bubble T: Converged on f(T) at iter=${iter}, T=${T_K.toFixed(2)}, f(T)=${fT.toExponential(3)}`);
            return {
                comp1_feed: x1_feed,
                comp1_equilibrium: Math.min(Math.max(y_iter[0], 0), 1),
                T_K: T_K,
                P_Pa: P_system_Pa,
                iterations: iter + 1,
                calculationType: 'bubbleT'
            };
        }

        const dfdT = numericalDerivative(objectiveFunction, T_K, T_K * 1e-4); // Relative step for h

        if (isNaN(dfdT) || Math.abs(dfdT) < 1e-9) {
            // console.error(`SRK Bubble T: Derivative is NaN or too small at T=${T_K.toFixed(2)}K, dfdT=${dfdT?.toExponential(2)}, iter=${iter}`);
            T_K += (fT > 0 ? -0.1 : 0.1); // Heuristic kick, same logic as PR
            if (T_K <= 0) T_K = initialTempGuess_K * (Math.random() * 0.2 + 0.9); // Reset
            // console.warn(`SRK Bubble T: Derivative issue. Perturbed T to ${T_K.toFixed(2)}K`);
            continue; // Skip update and retry
        }

        let deltaT = -fT / dfdT;
        const max_deltaT_abs = T_K * 0.20; // Max 20% change of current T
        if (Math.abs(deltaT) > max_deltaT_abs) {
            deltaT = Math.sign(deltaT) * max_deltaT_abs;
        }
        if (iter < 3) { // Stronger damping for initial iterations
            deltaT *= 0.5;
        }
        
        const T_K_new = T_K + deltaT;

        // console.debug(`SRK Bubble T iter ${iter}: T_old=${T_K.toFixed(2)}, f(T)=${fT.toExponential(3)}, dfdT=${dfdT.toExponential(3)}, deltaT=${deltaT.toFixed(3)}, T_new=${T_K_new.toFixed(2)}`);

        if (T_K_new <= 0) {
            // console.warn(`SRK Bubble T: Newton step resulted in T_new <= 0 (${T_K_new.toFixed(2)}). Halving step or resetting.`);
            deltaT /= 2;
            T_K = T_K + deltaT; // Re-calculate T_K with halved deltaT
            if (T_K <= 0) T_K = (T_K + initialTempGuess_K) / 2; // Fallback
            // console.warn(`SRK Bubble T: Adjusted T to ${T_K.toFixed(2)}`);
            continue;
        }
        
        if (Math.abs(deltaT) < tolerance_T_step && Math.abs(deltaT) < tolerance_T_step * Math.abs(T_K_new) ) {
            // console.log(`SRK Bubble T: Converged on T_step at iter=${iter}, T=${T_K_new.toFixed(2)}, deltaT=${deltaT.toExponential(3)}`);
            objectiveFunction(T_K_new); // Final update for K and y
            K_values_iter = (objectiveFunction as any)._last_K_values || K_values_iter;
            y_iter        = (objectiveFunction as any)._last_y_values || y_iter;
            return {
                comp1_feed: x1_feed,
                comp1_equilibrium: Math.min(Math.max(y_iter[0], 0), 1),
                T_K: T_K_new,
                P_Pa: P_system_Pa,
                iterations: iter + 1,
                calculationType: 'bubbleT'
            };
        }
        T_K = T_K_new;
    }

    // console.warn(`SRK Bubble T: Failed to converge for x1=${x1_feed.toFixed(4)} after ${maxIter} iterations. Last T=${T_K.toFixed(2)}K, f(T)=${objectiveFunction(T_K)?.toExponential(2)}`);
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
