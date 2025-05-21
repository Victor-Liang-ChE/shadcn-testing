import type { CompoundData, PrPureComponentParams, AntoineParams } from '@/lib/vle-types';
import { R_gas_const_J_molK, solveCubicEOS, type PrInteractionParams as BinaryPrInteractionParams } from './vle-calculations-pr'; // Import R and solver
import { calculatePsat_Pa } from './vle-calculations-unifac';

// --- Constants ---
// R_gas_const_J_molK is imported
const MAX_BUBBLE_T_ITERATIONS = 50; // Max iterations for outer temperature loop
const MAX_K_VALUE_ITERATIONS = 30;  // Max iterations for inner K-value/flash loop
const BUBBLE_T_TOLERANCE_K = 0.001; // Tolerance for temperature convergence
const K_VALUE_TOLERANCE = 1e-6;     // Tolerance for K-value convergence (sqrt(sum(Knew-Kold)^2))
const MIN_MOLE_FRACTION = 1e-9;
const BUBBLE_T_CONVERGENCE_TOLERANCE = 1e-5; // Tolerance for sum Ki*xi - 1
// const azeotrope_dxdt_threshold = 1e-7; // Azeotrope detection threshold // REMOVED

// --- Interfaces ---
export interface ResidueCurvePointODE {
    x: number[]; // Mole fractions [x0, x1, x2]
    T_K: number; // Temperature at this point
    step?: number; // Changed from xi to step
    // dxdt_norm?: number; // REMOVED
    // isPotentialAzeotrope?: boolean; // REMOVED
}
export type ResidueCurveODE = ResidueCurvePointODE[];

// This interface should match what page.tsx provides as TernaryPrParams
export interface TernaryPrParams {
    k01: number; k10?: number; 
    k02: number; k20?: number; 
    k12: number; k21?: number; 
}

// Add AzeotropeResult interface
export interface AzeotropeResult {
    x: number[];
    T_K: number;
    iterations: number;
    converged: boolean;
    errorNorm?: number;
}

function normalizeMoleFractions(x: number[]): number[] {
    const sum = x.reduce((acc, val) => acc + val, 0);
    if (sum === 0) return x.map(() => 1 / x.length); 
    const upper_bound_val = 1.0 - MIN_MOLE_FRACTION * (x.length > 1 ? (x.length - 1) : 0);
    return x.map(val => Math.max(MIN_MOLE_FRACTION, Math.min(upper_bound_val, val / sum)));
}

// --- Helper functions ---
function getKij(idx1: number, idx2: number, params: TernaryPrParams): number {
    const key = `${Math.min(idx1, idx2)}${Math.max(idx1, idx2)}`;
    if (idx1 === idx2) return 0; // k_ii = 0
    if (key === '01') return params.k01;
    if (key === '02') return params.k02;
    if (key === '12') return params.k12;
    // Support for k_ji if defined and different, though typically k_ij = k_ji
    if (idx1 === 1 && idx2 === 0 && params.k10 !== undefined) return params.k10;
    if (idx1 === 2 && idx2 === 0 && params.k20 !== undefined) return params.k20;
    if (idx1 === 2 && idx2 === 1 && params.k21 !== undefined) return params.k21;
    // Fallback to symmetric if specific k_ji not provided
    if (key === '01') return params.k01; // k10 defaults to k01
    if (key === '02') return params.k02; // k20 defaults to k02
    if (key === '12') return params.k12; // k21 defaults to k12
    return 0; 
}

function calculatePrPureComponentParams(
    T_K: number,
    prParams: PrPureComponentParams
): { ac_i: number, b_i: number, m_i: number, alpha_i: number, a_i: number } {
    const Tr_i = T_K / prParams.Tc_K;
    const m_i = 0.37464 + 1.54226 * prParams.omega - 0.26992 * prParams.omega * prParams.omega;
    const alpha_i = Math.pow(1 + m_i * (1 - Math.sqrt(Tr_i)), 2);
    const ac_i = 0.45724 * R_gas_const_J_molK * R_gas_const_J_molK * prParams.Tc_K * prParams.Tc_K / prParams.Pc_Pa;
    const a_i = ac_i * alpha_i;
    const b_i = 0.07780 * R_gas_const_J_molK * prParams.Tc_K / prParams.Pc_Pa;
    return { ac_i, b_i, m_i, alpha_i, a_i };
}

function calculatePrMixtureParamsTernary(
    phaseCompositions: number[], 
    T_K: number,
    componentsData: CompoundData[],
    ternaryKij: TernaryPrParams
): { mixture_a: number, mixture_b: number, component_a_values: number[], component_b_values: number[] } {
    const numComponents = componentsData.length;
    const component_a_values = componentsData.map(c => calculatePrPureComponentParams(T_K, c.prParams!).a_i);
    const component_b_values = componentsData.map(c => calculatePrPureComponentParams(T_K, c.prParams!).b_i);

    let mixture_a = 0;
    for (let i = 0; i < numComponents; i++) {
        for (let j = 0; j < numComponents; j++) {
            const kij = getKij(i, j, ternaryKij);
            mixture_a += phaseCompositions[i] * phaseCompositions[j] * Math.sqrt(component_a_values[i] * component_a_values[j]) * (1 - kij);
        }
    }
    const mixture_b = phaseCompositions.reduce((sum, x_i, i) => sum + x_i * component_b_values[i], 0);
    return { mixture_a, mixture_b, component_a_values, component_b_values };
}

function calculatePrCoefficients(mixture_a: number, mixture_b: number, T_K: number, P_Pa: number): { p: number, q: number, r: number, A_mix_eos: number, B_mix_eos: number } {
    const A_mix_eos = mixture_a * P_Pa / Math.pow(R_gas_const_J_molK * T_K, 2);
    const B_mix_eos = mixture_b * P_Pa / (R_gas_const_J_molK * T_K);

    // For PR EOS: Z^3 - (1-B)Z^2 + (A-3B^2-2B)Z - (AB-B^2-B^3) = 0
    // Z^3 + pZ^2 + qZ + r = 0
    const p_coeff = -(1.0 - B_mix_eos);
    const q_coeff = A_mix_eos - 3.0 * B_mix_eos * B_mix_eos - 2.0 * B_mix_eos;
    const r_coeff = -(A_mix_eos * B_mix_eos - B_mix_eos * B_mix_eos - B_mix_eos * B_mix_eos * B_mix_eos);

    return { p: p_coeff, q: q_coeff, r: r_coeff, A_mix_eos, B_mix_eos };
}

// --- TERNARY FUGACITY CALCULATION ---
function calculatePrFugacityCoefficientsTernary(
    phaseCompositions: number[], // x_i or y_i (should be 3 components)
    T_K: number, P_Pa: number, Z_phase: number,
    mixture_a: number, mixture_b: number,
    component_a_values: number[], component_b_values: number[], // pure component a_i and b_i values at T_K
    ternaryKij: TernaryPrParams // Use the ternary k_ij params
): number[] { // Returns array of 3 phi_i
    const numComponents = phaseCompositions.length;
    if (numComponents !== 3) {
        console.error("calculatePrFugacityCoefficientsTernary called with non-ternary composition");
        return [NaN, NaN, NaN];
    }
    const A_mix_eos = mixture_a * P_Pa / Math.pow(R_gas_const_J_molK * T_K, 2);
    const B_mix_eos = mixture_b * P_Pa / (R_gas_const_J_molK * T_K);
    const fugacity_coeffs: number[] = [];

    for (let k = 0; k < numComponents; k++) {
        let sum_term_a_ki = 0;
        for (let i = 0; i < numComponents; i++) {
            const kij_ki = getKij(k, i, ternaryKij); // Correctly gets k_ki
            sum_term_a_ki += phaseCompositions[i] * Math.sqrt(component_a_values[k] * component_a_values[i]) * (1 - kij_ki);
        }

        const term1 = component_b_values[k] / mixture_b * (Z_phase - 1);

        const Z_minus_B = Z_phase - B_mix_eos;
        if (Z_minus_B <= MIN_MOLE_FRACTION) {
            fugacity_coeffs.push(NaN);
            continue;
        }
        const term2 = -Math.log(Z_minus_B);

        let term3_factor = 0;
        const denominator_term3_factor = 2 * Math.sqrt(2) * B_mix_eos;
        if (Math.abs(denominator_term3_factor) > 1e-12) { // Normal case: B_mix_eos is not zero
            term3_factor = A_mix_eos / denominator_term3_factor;
        } else if (Math.abs(A_mix_eos) < 1e-12) { // B_mix_eos is zero (or very small) AND A_mix_eos is also zero
            term3_factor = 0; // Both A and B terms are zero, so factor is zero
        } else { // B_mix_eos is zero (or very small) but A_mix_eos is not.
            // This will lead to a very large (or Inf) term3_factor.
            // Let it proceed; subsequent checks for NaN/Inf on ln_phi_k or K_values should catch it.
            term3_factor = A_mix_eos / denominator_term3_factor; // This will likely be Infinity or -Infinity
        }

        const term3_main_num = 2 * sum_term_a_ki; // This is sum_i (x_i * a_ki_sqrt) where a_ki_sqrt = sqrt(a_k a_i)(1-k_ki)
        let term3_main_frac1 = 0;
        if (Math.abs(mixture_a) > 1e-12) term3_main_frac1 = term3_main_num / mixture_a;

        let term3_main_frac2 = 0;
        if (Math.abs(mixture_b) > 1e-12) term3_main_frac2 = component_b_values[k] / mixture_b;
        
        const term3_main = term3_main_frac1 - term3_main_frac2;

        const log_arg_num = Z_phase + (1 + Math.sqrt(2)) * B_mix_eos;
        const log_arg_den = Z_phase + (1 - Math.sqrt(2)) * B_mix_eos;
        let term3_log = 0;
        if (log_arg_den > MIN_MOLE_FRACTION && log_arg_num / log_arg_den > MIN_MOLE_FRACTION) { // Ensure arg of log is positive
            term3_log = Math.log(log_arg_num / log_arg_den);
        } else {
            fugacity_coeffs.push(NaN);
            continue;
        }

        const ln_phi_k = term1 + term2 - term3_factor * term3_main * term3_log;
        if (isNaN(ln_phi_k) || !isFinite(ln_phi_k)) {
            fugacity_coeffs.push(NaN);
        } else {
            fugacity_coeffs.push(Math.exp(ln_phi_k));
        }
    }
    return fugacity_coeffs;
}

export function findBubblePointTemperatureTernaryPr(
    x_liq: number[],
    P_system_Pa: number,
    components: CompoundData[],
    ternaryPrKij: TernaryPrParams, 
    initialTempGuess_K: number
): { T_K: number, y_vap: number[], K_values: number[], ZL: number, ZV: number } | null {
    // console.log(`--- findBubblePointTemperatureTernaryPr START ---`);
    // console.log(`Initial x_liq: [${x_liq.join(', ')}], P_sys: ${P_system_Pa}, Init_T_guess: ${initialTempGuess_K}`);

    let T_K = initialTempGuess_K;
    const x = normalizeMoleFractions(x_liq);
    // console.log(`Normalized x: [${x.join(', ')}]`);

    if (components.some(c => !c.antoine || !c.prParams) || components.length !== 3) {
        // console.error("PR BubbleT ERR: Invalid component data or not 3 components.");
        return null;
    }

    // let y_vap = [...x]; // y_vap will be calculated inside K-loop as y_vap_final
    let K_values = components.map(c => {
        const psat = calculatePsat_Pa(c.antoine!, T_K);
        return isNaN(psat) || P_system_Pa === 0 ? 1.0 : psat / P_system_Pa; 
    });
    if (K_values.some(isNaN)) K_values = [1, 1, 1]; 
    // console.log(`Initial K_values (Raoult): [${K_values.map(k => k.toFixed(3)).join(', ')}]`);

    let sum_Kx_from_K_loop = 0; 
    let ZL_K_loop: number | null = null; 
    let ZV_K_loop: number | null = null; 

    for (let T_iter = 0; T_iter < MAX_BUBBLE_T_ITERATIONS; T_iter++) {
        // console.log(`\nT_iter: ${T_iter}, Current T_K: ${T_K.toFixed(3)}`);
        
        let K_converged_inner_loop = false;
        // K_values holds the result from the previous T_iter or the initial Raoult's guess.

        for (let K_iter = 0; K_iter < MAX_K_VALUE_ITERATIONS; K_iter++) {
            const K_values_at_K_iter_start = [...K_values]; 

            // 1. Calculate y_vap_final based on K_values_at_K_iter_start
            const y_unnormalized_for_iter = K_values_at_K_iter_start.map((Ki, i) => Ki * x[i]);
            const sum_Kx_for_y_norm = y_unnormalized_for_iter.reduce((sum, val) => sum + val, 0);

            if (Math.abs(sum_Kx_for_y_norm) < MIN_MOLE_FRACTION * 1e-3) {
                ZL_K_loop = null; ZV_K_loop = null; // Ensure failure state for T-loop
                break; // K-loop fails
            }
            const y_vap_final = normalizeMoleFractions(y_unnormalized_for_iter);
            // console.log(`  K_iter: ${K_iter}, y_vap_final guess: [${y_vap_final.map(yi => yi.toFixed(3)).join(', ')}]`);
            if (y_vap_final.some(isNaN)) { ZL_K_loop = null; ZV_K_loop = null; break; }


            // 2. Calculate mixture params, ZL, ZV based on y_vap_final
            const liqMix = calculatePrMixtureParamsTernary(x, T_K, components, ternaryPrKij);
            const vapMix = calculatePrMixtureParamsTernary(y_vap_final, T_K, components, ternaryPrKij);
            const liqEosCoeffs = calculatePrCoefficients(liqMix.mixture_a, liqMix.mixture_b, T_K, P_system_Pa);
            const vapEosCoeffs = calculatePrCoefficients(vapMix.mixture_a, vapMix.mixture_b, T_K, P_system_Pa);
            const ZL_all_roots = solveCubicEOS(1, liqEosCoeffs.p, liqEosCoeffs.q, liqEosCoeffs.r);
            const ZV_all_roots = solveCubicEOS(1, vapEosCoeffs.p, vapEosCoeffs.q, vapEosCoeffs.r);
            
            ZL_K_loop = null; ZV_K_loop = null; // Reset for this K-iteration's Z calculation
            if (ZL_all_roots?.length) {
                const valid_ZL = ZL_all_roots.filter(z => z > liqEosCoeffs.B_mix_eos + MIN_MOLE_FRACTION * 1e-3);
                if (valid_ZL.length > 0) ZL_K_loop = Math.min(...valid_ZL);
            }
            if (ZV_all_roots?.length) {
                const valid_ZV = ZV_all_roots.filter(z => z > vapEosCoeffs.B_mix_eos + MIN_MOLE_FRACTION * 1e-3);
                if (valid_ZV.length > 0) ZV_K_loop = Math.max(...valid_ZV);
            }
            // console.log(`    Selected ZL=${ZL_K_loop?.toFixed(4)}, ZV=${ZV_K_loop?.toFixed(4)} (B_liq=${liqEosCoeffs.B_mix_eos.toFixed(4)}, B_vap=${vapEosCoeffs.B_mix_eos.toFixed(4)})`);
            if (ZL_K_loop === null || ZV_K_loop === null) {
                // console.warn(`  PR Bubble T (K-loop) ERR: No valid ZL/ZV. T_K=${T_K.toFixed(2)}, K_iter=${K_iter}`);
                break; // K-loop fails, ZL_K_loop/ZV_K_loop are null for T-loop check
            }

            // 3. Calculate phi_L and phi_V
            const phi_L = calculatePrFugacityCoefficientsTernary(x, T_K, P_system_Pa, ZL_K_loop, liqMix.mixture_a, liqMix.mixture_b, liqMix.component_a_values, liqMix.component_b_values, ternaryPrKij);
            const phi_V = calculatePrFugacityCoefficientsTernary(y_vap_final, T_K, P_system_Pa, ZV_K_loop, vapMix.mixture_a, vapMix.mixture_b, vapMix.component_a_values, vapMix.component_b_values, ternaryPrKij);
            // console.log(`    phi_L: [${phi_L.map(p => p?.toFixed(3) ?? 'NaN').join(', ')}], phi_V: [${phi_V.map(p => p?.toFixed(3) ?? 'NaN').join(', ')}]`);
            if (phi_L.some(isNaN) || phi_V.some(isNaN) || phi_V.some(p => p < MIN_MOLE_FRACTION * 1e-3)) {
                // console.error(`  PR Bubble T (K-loop) ERR: NaN or too small phi_V. T_K=${T_K.toFixed(2)}`);
                ZL_K_loop = null; ZV_K_loop = null; // Signal failure
                break; 
            }
            
            // 4. Calculate new K-values for this K-iteration
            const K_values_updated_this_K_iter = phi_L.map((phiL_i, i) => Math.max(MIN_MOLE_FRACTION, phiL_i / phi_V[i]));
            if (K_values_updated_this_K_iter.some(k => isNaN(k) || !isFinite(k))) {
                // console.error(`  PR Bubble T (K-loop) ERR: NaN/Inf in K_values_updated_this_K_iter. T_K=${T_K.toFixed(2)}`);
                ZL_K_loop = null; ZV_K_loop = null; // Signal failure
                break;
            }
            // console.log(`    K_values_updated_this_K_iter: [${K_values_updated_this_K_iter.map(k => k.toFixed(3)).join(', ')}]`);

            // 5. Check convergence: Compare K_values_updated_this_K_iter with K_values_at_K_iter_start
            const K_diff_sq_sum = K_values_at_K_iter_start.reduce((sum, oldK_for_iter, i) => sum + (K_values_updated_this_K_iter[i] - oldK_for_iter)**2, 0);
            
            // Update K_values for the *next* K-iteration or for the T-loop using damping
            const K_damping = (K_iter < 3 ? 0.3 : 0.7); // More damping initially
            K_values = K_values_at_K_iter_start.map((oldK, idx) => oldK + K_damping * (K_values_updated_this_K_iter[idx] - oldK));
            
            // This sum_Kx_from_K_loop is what the T-loop tries to drive to 1.
            sum_Kx_from_K_loop = K_values.reduce((sum, Ki, i) => sum + Ki * x[i], 0);

            // console.log(`    Damped K_values for next iter/T-loop: [${K_values.map(k => k.toFixed(3)).join(', ')}]`);
            // console.log(`    Sum Kx (for T-loop) = ${sum_Kx_from_K_loop.toFixed(5)}, K_diff_sqrt (updated vs start_of_K_iter) = ${Math.sqrt(K_diff_sq_sum).toExponential(3)}`);

            if (Math.sqrt(K_diff_sq_sum) < K_VALUE_TOLERANCE && K_iter > 0) {
                K_converged_inner_loop = true;
                // console.log(`  K-loop converged at K_iter=${K_iter}.`);
                break; 
            }
             if (K_iter === MAX_K_VALUE_ITERATIONS - 1) {
                // console.warn(`  PR Bubble T: K-value loop MAX iter at T=${T_K.toFixed(2)}K.`);
             }
        } // End K-value loop

        if (ZL_K_loop === null || ZV_K_loop === null || !K_converged_inner_loop) { 
            // console.warn(`  Adjusting T due to K-loop failure/non-convergence. Old T_K=${T_K.toFixed(2)}. ZL/ZV null: ${ZL_K_loop === null || ZV_K_loop === null}, K_converged: ${K_converged_inner_loop}`);
            T_K *= (1 + (Math.random() - 0.5) * 0.04); 
            if (T_K <= 100 || isNaN(T_K)) T_K = initialTempGuess_K * (0.8 + Math.random() * 0.4);
            K_values = components.map(c => {
                const psat = calculatePsat_Pa(c.antoine!, T_K); 
                return isNaN(psat) || P_system_Pa === 0 ? 1.0 : psat / P_system_Pa;
            });
            if (K_values.some(isNaN)) K_values = [1, 1, 1];
            // console.warn(`  New T_K for next T_iter = ${T_K.toFixed(2)}, K-values reset to Raoult's.`);
            continue; 
        }

        const error_T_loop = sum_Kx_from_K_loop - 1.0;

        if (Math.abs(error_T_loop) < BUBBLE_T_CONVERGENCE_TOLERANCE) {
            const final_y_vap_for_return = normalizeMoleFractions(K_values.map((Ki, i) => Ki * x[i]));
            return { T_K, y_vap: final_y_vap_for_return, K_values, ZL: ZL_K_loop!, ZV: ZV_K_loop! };
        }

        if (Math.abs(error_T_loop) > 0.2 && T_iter > 1) { 
            // console.warn(`  Sum Kx = ${sum_Kx_from_K_loop.toFixed(3)} is far from 1. Resetting K-values to Raoult's.`);
            K_values = components.map(c => {
                const psat = calculatePsat_Pa(c.antoine!, T_K); 
                return isNaN(psat) || P_system_Pa === 0 ? 1.0 : psat / P_system_Pa;
            });
            if (K_values.some(isNaN)) K_values = [1, 1, 1];
        }
        
        const target_T_ratio = 1.0 / sum_Kx_from_K_loop;
        let damped_T_ratio_change = (target_T_ratio - 1.0) * (T_iter < 3 ? 0.05 : (T_iter < 10 ? 0.1 : 0.15));
        
        const max_T_change_factor = 0.05; 
        damped_T_ratio_change = Math.max(-max_T_change_factor, Math.min(max_T_change_factor, damped_T_ratio_change));

        let T_K_new_ratio = 1.0 + damped_T_ratio_change;
        const T_K_new = T_K * T_K_new_ratio;

        if (Math.abs(T_K_new - T_K) < BUBBLE_T_TOLERANCE_K && T_iter > 2) {
            if (Math.abs(error_T_loop) < BUBBLE_T_CONVERGENCE_TOLERANCE * 20) { // Looser tolerance if T is stable
                const final_y_vap_for_return = normalizeMoleFractions(K_values.map((Ki, i) => Ki * x[i]));
                return { T_K, y_vap: final_y_vap_for_return, K_values, ZL: ZL_K_loop!, ZV: ZV_K_loop! };
            }
        }
        T_K = Math.max(150, Math.min(800, T_K_new)); // Ensure T_K stays within reasonable bounds
        if (isNaN(T_K) || T_K <= 0) {
            T_K = initialTempGuess_K * (0.7 + Math.random() * 0.6); // More aggressive reset if T becomes invalid
            if (T_K <= 0) T_K = 273.15; // Absolute fallback
        }
    }
    // console.warn(`PR Bubble T ERR: Max T-iterations reached for x=[${x.join(', ')}]. Last T_K=${T_K.toFixed(2)}, Last Sum Kx=${sum_Kx_from_K_loop.toFixed(4)}`);
    return null;
}

function calculateVaporCompositionEos(
    x_liq: number[], K_values: number[]
): number[] | null {
    const y_unnormalized = x_liq.map((xi: number, i: number) => K_values[i] * xi);
    const sum_y_unnormalized = y_unnormalized.reduce((s: number, val: number) => s + val, 0);
    if (sum_y_unnormalized === 0 || isNaN(sum_y_unnormalized)) return null;
    const y_vap = y_unnormalized.map((val: number) => val / sum_y_unnormalized);
    return normalizeMoleFractions(y_vap);
}

function residue_curve_ode_dxdxi_pr(
    x_liq_current: number[],
    P_system_Pa: number,
    components: CompoundData[],
    ternaryPrParams: TernaryPrParams,
    initialTempGuess_K: number
): { dxdxi: number[], T_bubble_K: number, y_vap: number[] } | null { // Added y_vap to return
    const x_norm = normalizeMoleFractions(x_liq_current);

    const bubbleResult = findBubblePointTemperatureTernaryPr(x_norm, P_system_Pa, components, ternaryPrParams, initialTempGuess_K);
    if (!bubbleResult) return null;

    const { T_K: T_bubble_K, K_values } = bubbleResult;

    const y_vap = calculateVaporCompositionEos(x_norm, K_values);
    if (!y_vap) return null;

    const dxdxi = [x_norm[0] - y_vap[0], x_norm[1] - y_vap[1], x_norm[2] - y_vap[2]];
    return { dxdxi, T_bubble_K, y_vap }; // Return y_vap
}

// --- Constants for simulateSingleResidueCurveODE_Pr ---
const AZEOTROPE_NORM_THRESHOLD = 1e-5; // If ||x-y|| is below this, potential azeotrope
const STAGNATION_DXDXI_NORM_THRESHOLD = 1e-7; // If ||x-y|| is very small, strong stagnation/azeotrope
const MIN_DIST_SQ_TO_ADD_POINT_SIM = 1e-10 * 1e-10; // sqrt is 1e-10; for adding points to curve
const PURE_COMPONENT_THRESHOLD = 0.9999; // Threshold for a component to be considered "pure"
const OTHER_COMPONENTS_TINY_THRESHOLD = 1e-4; // Threshold for other components when one is pure

export async function simulateSingleResidueCurveODE_Pr(
    initial_x: number[], 
    P_system_Pa: number,
    componentsData: CompoundData[], 
    ternaryPrParams: TernaryPrParams, 
    fixed_d_xi_step: number, // Renamed from d_xi_step
    max_steps_per_direction: number,
    initialTemp_K: number
): Promise<ResidueCurveODE | null> {

    const run_direction = async (start_x: number[], forward: boolean, current_T_guess: number): Promise<ResidueCurveODE> => {
        const curve: ResidueCurveODE = [];
        let current_x = [...start_x];
        let last_T_K = current_T_guess;

        for (let loop_step = 0; loop_step < max_steps_per_direction; loop_step++) { 
            const ode_result = residue_curve_ode_dxdxi_pr(current_x, P_system_Pa, componentsData, ternaryPrParams, last_T_K);
            if (!ode_result) break;

            const { dxdxi, T_bubble_K, y_vap } = ode_result; 
            last_T_K = T_bubble_K;

            const dxdxi_norm = Math.sqrt(dxdxi[0]**2 + dxdxi[1]**2 + dxdxi[2]**2);

            const currentPointData: ResidueCurvePointODE = {
                x: normalizeMoleFractions(current_x), // Ensure current_x is normalized before storing
                T_K: T_bubble_K,
                step: (forward ? 1 : -1) * loop_step * fixed_d_xi_step,
            };
            
            if (curve.length > 0) {
                const lastPoint = curve[curve.length - 1];
                const x_diff_sq = lastPoint.x.reduce((sum, val, i) => sum + (val - currentPointData.x[i])**2, 0);
                if (x_diff_sq < MIN_DIST_SQ_TO_ADD_POINT_SIM && Math.abs(lastPoint.T_K - currentPointData.T_K) < BUBBLE_T_TOLERANCE_K * 10) {
                    // Point is virtually identical to the last one, likely stagnation
                    if (dxdxi_norm < STAGNATION_DXDXI_NORM_THRESHOLD * 10) break; // More aggressive break if dxdxi is also tiny
                    // otherwise, might just be a very small step, continue but don't add
                } else {
                    curve.push(currentPointData);
                }
            } else { // First point of this direction
                curve.push(currentPointData);
            }

            if (dxdxi_norm < STAGNATION_DXDXI_NORM_THRESHOLD) {
                console.log(`Stagnation detected (dxdxi_norm < ${STAGNATION_DXDXI_NORM_THRESHOLD.toExponential()}) at step ${loop_step}`);
                break; 
            }
            if (loop_step > 5 && dxdxi_norm < AZEOTROPE_NORM_THRESHOLD) { // Check after a few steps
                console.log(`Potential azeotrope detected (dxdxi_norm < ${AZEOTROPE_NORM_THRESHOLD.toExponential()}) at step ${loop_step}`);
                break;
            }

            const next_x_unnormalized = [
                current_x[0] + (forward ? 1 : -1) * dxdxi[0] * fixed_d_xi_step,
                current_x[1] + (forward ? 1 : -1) * dxdxi[1] * fixed_d_xi_step,
                current_x[2] + (forward ? 1 : -1) * dxdxi[2] * fixed_d_xi_step,
            ];
            current_x = normalizeMoleFractions(next_x_unnormalized);

            // Pure component check
            const isPureComponentTermination = current_x.some((val, idx, arr) => {
                if (val > PURE_COMPONENT_THRESHOLD) {
                    return arr.every((otherVal, otherIdx) => otherIdx === idx || otherVal < OTHER_COMPONENTS_TINY_THRESHOLD);
                }
                if (val < (1.0 - PURE_COMPONENT_THRESHOLD)) { // e.g. val < 0.0001
                    // Check if sum of other two is high
                    let sumOfOthers = 0;
                    for(let k=0; k<arr.length; k++) if(k !== idx) sumOfOthers += arr[k];
                    if (sumOfOthers > PURE_COMPONENT_THRESHOLD) { // This implies one of the others is high
                        // This condition is a bit complex, let's simplify:
                        // If one component is very high OR very low (and others make up the bulk)
                    }
                }
                return false;
            });
            
            // Simpler pure component check: if any component is > PURE_COMPONENT_THRESHOLD or < (1-PURE_COMPONENT_THRESHOLD)
            // and the others are small enough.
            let pureFound = false;
            for(let i=0; i<current_x.length; i++) {
                if (current_x[i] > PURE_COMPONENT_THRESHOLD) {
                    let othersSmall = true;
                    for(let j=0; j<current_x.length; j++) {
                        if (i === j) continue;
                        if (current_x[j] > OTHER_COMPONENTS_TINY_THRESHOLD) {
                            othersSmall = false;
                            break;
                        }
                    }
                    if (othersSmall) { pureFound = true; break; }
                }
            }


            if (pureFound) { 
                const final_x_normalized = normalizeMoleFractions(current_x); // Ensure it's fully normalized to the vertex
                 // Add this terminal point if it's different enough from the last added one
                if (curve.length === 0 || curve[curve.length-1].x.reduce((s,v,i) => s + (v-final_x_normalized[i])**2, 0) > MIN_DIST_SQ_TO_ADD_POINT_SIM) {
                    curve.push({ 
                        x: final_x_normalized, 
                        T_K: last_T_K, 
                        step: (forward ? 1 : -1) * (loop_step + 1) * fixed_d_xi_step, 
                    });
                }
                console.log(`Pure component vertex reached at step ${loop_step}`);
                break;
            }
            
            // Stagnation check based on change in x per step
            if (loop_step > 0 && curve.length > 0) { 
                const lastAddedPoint = curve[curve.length - 1];
                // x_before_euler_step was the x that *produced* lastAddedPoint via ODE
                // So, compare lastAddedPoint.x with current_x (which is x_after_euler_step)
                const diff_sq_from_last_added_x = lastAddedPoint.x.reduce((sum, px, i) => sum + (px - current_x[i]) ** 2, 0);
                
                // If the change in x from the last *added* point to the current *candidate* point (after Euler step)
                // is much smaller than the step size itself, it implies stagnation.
                if (Math.sqrt(diff_sq_from_last_added_x) < fixed_d_xi_step * 1e-4) { 
                    console.log(`Stagnation (small dx) detected at step ${loop_step}`);
                    break; 
                }
            }
        }
        return curve;
    };

    const initial_x_norm = normalizeMoleFractions(initial_x);
    const initialBubbleResult = findBubblePointTemperatureTernaryPr(initial_x_norm, P_system_Pa, componentsData, ternaryPrParams, initialTemp_K);
    if (!initialBubbleResult) {
        console.warn("PR simulateSingle: Initial bubble point failed. Returning null.");
        return null;
    }

    const T_start_K = initialBubbleResult.T_K;

    const curve_forward = await run_direction(initial_x_norm, true, T_start_K);
    const curve_backward = await run_direction(initial_x_norm, false, T_start_K);

    let combined_curve: ResidueCurveODE = [];

    if (curve_backward.length > 0) {
        combined_curve = curve_backward.reverse(); // Starts far, ends at P0
    }

    if (curve_forward.length > 0) {
        if (combined_curve.length > 0) { 
            // curve_forward starts with P0. combined_curve currently ends with P0.
            // Pop P0 from combined_curve before appending curve_forward (which starts with P0)
            combined_curve.pop(); 
            combined_curve.push(...curve_forward);
        } else {
            // Only forward ran
            combined_curve = curve_forward;
        }
    }
    // If only backward ran, combined_curve is already correctly set.
    // If neither ran (e.g. run_direction returned empty for both), combined_curve is [].

    // If the combined curve is empty after attempts, but we had a valid start point, add it.
    // This ensures at least the starting point is plotted if simulations fail immediately.
    if (combined_curve.length === 0 && initial_x_norm) {
         console.log("PR simulateSingle: Both directions yielded no points, adding initial point.");
         combined_curve.push({ x: initial_x_norm, T_K: T_start_K, step: 0 });
    }
    
    console.log(`PR simulateSingle: Combined curve length: ${combined_curve.length}. Forward: ${curve_forward.length}, Backward: ${curve_backward.length}`);
    return combined_curve.length > 0 ? combined_curve : null;
}

// Assuming normalizeMoleFractions and findBubblePointTemperatureTernaryPr are defined in this file or imported.
// Ensure CompoundData and TernaryPrParams types are correctly defined/imported.

export function findAzeotropePr(
    P_system_Pa: number,
    components: CompoundData[], // Each component must have prParams
    ternaryPrKij: TernaryPrParams,
    initialGuess_x: number[], // Initial guess for mole fractions [x0, x1, x2]
    initialGuess_T_K: number, // Initial guess for temperature
    maxIterations: number = 150, // Max iterations for the azeotrope search (increased slightly)
    tolerance: number = 1e-7,    // Tolerance for ||y-x|| norm for azeotrope convergence
    dampingFactor: number = 0.1 // Damping factor for x updates; Reduced from 0.25
): AzeotropeResult | null {

    if (components.length !== 3 || components.some(c => !c.prParams)) {
        console.error("Azeotrope search (PR): Requires 3 components with PR parameters.");
        return null;
    }
    if (initialGuess_x.length !== 3) {
        console.error("Azeotrope search (PR): Initial guess for x must have 3 components.");
        return null;
    }

    let current_x = normalizeMoleFractions([...initialGuess_x]);
    let current_T_K = initialGuess_T_K;
    let last_errorNorm = Infinity;

    console.log(`--- findAzeotropePr START ---`);
    console.log(`Initial x_guess: [${current_x.join(', ')}], Init_T_guess: ${current_T_K.toFixed(2)}`);


    for (let iter = 0; iter < maxIterations; iter++) {
        // 1. Calculate bubble point T and corresponding vapor composition y_vap
        //    for the current liquid composition estimate (current_x) at P_system_Pa.
        const bubbleResult = findBubblePointTemperatureTernaryPr(
            current_x,
            P_system_Pa,
            components,
            ternaryPrKij,
            current_T_K // Use current T_K as the initial guess for the bubble T solver
        );

        if (!bubbleResult || !bubbleResult.y_vap) {
            console.warn(`Azeotrope search (PR, iter ${iter}): Bubble point calculation failed for x = [${current_x.map(v => v.toFixed(4)).join(', ')}], T_guess=${current_T_K.toFixed(2)}`);
            // Attempt a slight perturbation of T_K if bubble point fails, or return current state.
            // This might indicate the guess is too far off or in a non-physical region.
            current_T_K *= (iter % 2 === 0 ? 1.02 : 0.98); // Small nudge
            if (iter > 5 && last_errorNorm === Infinity) { // If consistently failing early
                 return { x: current_x, T_K: initialGuess_T_K, iterations: iter, converged: false, errorNorm: last_errorNorm };
            }
            continue; // Try next iteration with potentially adjusted T_K if bubbleResult was null
        }

        current_T_K = bubbleResult.T_K; // Update T_K to the converged bubble point temperature
        const y_vap_at_bubble_T = normalizeMoleFractions(bubbleResult.y_vap);

        // 2. Check for azeotropic condition: ||y - x|| < tolerance
        const errorVector = [
            y_vap_at_bubble_T[0] - current_x[0],
            y_vap_at_bubble_T[1] - current_x[1],
            y_vap_at_bubble_T[2] - current_x[2]
        ];
        const errorNorm = Math.sqrt(
            errorVector[0]**2 + errorVector[1]**2 + errorVector[2]**2
        );
        last_errorNorm = errorNorm;
        
        console.log(`  AzeoIter ${iter}: T_K=${current_T_K.toFixed(2)}, x=[${current_x.map(v=>v.toFixed(4)).join(', ')}], y=[${y_vap_at_bubble_T.map(v=>v.toFixed(4)).join(', ')}], ||y-x||=${errorNorm.toExponential(3)}`);


        if (errorNorm < tolerance) {
            console.log(`Azeotrope (PR) converged at iter ${iter}.`);
            return { x: normalizeMoleFractions(current_x), T_K: current_T_K, iterations: iter, converged: true, errorNorm };
        }

        // 3. Update liquid composition estimate (current_x) for the next iteration.
        //    Using damped fixed-point iteration: x_new = x_old + lambda * (y_calc - x_old)
        let next_x_unnormalized = [
            current_x[0] + dampingFactor * errorVector[0],
            current_x[1] + dampingFactor * errorVector[1],
            current_x[2] + dampingFactor * errorVector[2],
        ];
        current_x = normalizeMoleFractions(next_x_unnormalized);
        
        // Optional: Check for stagnation if x isn't changing much but error is still high
        if (iter > 10 && errorNorm > last_errorNorm * 0.999 && errorNorm > tolerance * 100) {
             // console.warn(`Azeotrope (PR) search may be stagnating or diverging. errorNorm: ${errorNorm.toExponential(3)}`);
        }

    }

    console.warn(`Azeotrope search (PR): Did not converge after ${maxIterations} iterations. Final ||y-x||: ${last_errorNorm.toExponential(3)} for x=[${current_x.map(v=>v.toFixed(4)).join(',')}] T=${current_T_K.toFixed(2)}`);
    return { x: current_x, T_K: current_T_K, iterations: maxIterations, converged: false, errorNorm: last_errorNorm };
}

