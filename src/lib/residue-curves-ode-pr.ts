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

// --- Interfaces ---
export interface ResidueCurvePointODE {
    x: number[]; // Mole fractions [x0, x1, x2]
    T_K: number; // Temperature at this point
    step?: number; // Changed from xi to step
}
export type ResidueCurveODE = ResidueCurvePointODE[];

// This interface should match what page.tsx provides as TernaryPrParams
export interface TernaryPrParams {
    k01: number; k10?: number; 
    k02: number; k20?: number; 
    k12: number; k21?: number; 
}

function normalizeMoleFractions(x: number[]): number[] {
    const sum = x.reduce((acc, val) => acc + val, 0);
    if (sum === 0) return x.map(() => 1 / x.length); // Should not happen with valid inputs
    return x.map(val => Math.max(MIN_MOLE_FRACTION, Math.min(1.0 - MIN_MOLE_FRACTION * (x.length - 1), val / sum)));
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
    ternaryPrKij: TernaryPrParams, // Changed variable name and ensured type is TernaryPrParams
    initialTempGuess_K: number
): { T_K: number, y_vap: number[], K_values: number[], ZL: number, ZV: number } | null {
    console.log(`--- findBubblePointTemperatureTernaryPr START ---`);
    console.log(`Initial x_liq: [${x_liq.join(', ')}], P_sys: ${P_system_Pa}, Init_T_guess: ${initialTempGuess_K}`);

    let T_K = initialTempGuess_K;
    const x = normalizeMoleFractions(x_liq);
    console.log(`Normalized x: [${x.join(', ')}]`);

    if (components.some(c => !c.antoine || !c.prParams) || components.length !== 3) {
        console.error("PR BubbleT ERR: Invalid component data or not 3 components.");
        return null;
    }

    let y_vap = [...x]; // Initialize y_vap
    let K_values = components.map(c => {
        const psat = calculatePsat_Pa(c.antoine!, T_K);
        return isNaN(psat) || P_system_Pa === 0 ? 1.0 : psat / P_system_Pa; // Corrected typo
    });
    if (K_values.some(isNaN)) K_values = [1, 1, 1]; // Fallback if any K is NaN
    console.log(`Initial K_values (Raoult): [${K_values.map(k => k.toFixed(3)).join(', ')}]`);

    let sum_Kx_from_K_loop = 0; // Will hold the sum Kx from the K-loop
    let ZL_K_loop: number | null = null; // To store ZL from K-loop for final check
    let ZV_K_loop: number | null = null; // To store ZV from K-loop for final check

    for (let T_iter = 0; T_iter < MAX_BUBBLE_T_ITERATIONS; T_iter++) {
        console.log(`\nT_iter: ${T_iter}, Current T_K: ${T_K.toFixed(3)}`);
        
        for (let K_iter = 0; K_iter < MAX_K_VALUE_ITERATIONS; K_iter++) {
            y_vap = normalizeMoleFractions(K_values.map((Ki, i) => Ki * x[i]));
            console.log(`  K_iter: ${K_iter}, y_vap guess: [${y_vap.map(yi => yi.toFixed(3)).join(', ')}]`);

            const liqMix = calculatePrMixtureParamsTernary(x, T_K, components, ternaryPrKij);
            const vapMix = calculatePrMixtureParamsTernary(y_vap, T_K, components, ternaryPrKij);

            const liqEosCoeffs = calculatePrCoefficients(liqMix.mixture_a, liqMix.mixture_b, T_K, P_system_Pa);
            const vapEosCoeffs = calculatePrCoefficients(vapMix.mixture_a, vapMix.mixture_b, T_K, P_system_Pa);

            const ZL_all_roots = solveCubicEOS(1, liqEosCoeffs.p, liqEosCoeffs.q, liqEosCoeffs.r);
            const ZV_all_roots = solveCubicEOS(1, vapEosCoeffs.p, vapEosCoeffs.q, vapEosCoeffs.r);
            
            console.log(`    Liq EOS: A=${liqEosCoeffs.A_mix_eos.toExponential(2)}, B=${liqEosCoeffs.B_mix_eos.toExponential(2)}, p=${liqEosCoeffs.p.toFixed(3)}, q=${liqEosCoeffs.q.toFixed(3)}, r=${liqEosCoeffs.r.toFixed(3)}, Roots=${ZL_all_roots?.map(r => r.toFixed(3))}`);
            console.log(`    Vap EOS: A=${vapEosCoeffs.A_mix_eos.toExponential(2)}, B=${vapEosCoeffs.B_mix_eos.toExponential(2)}, p=${vapEosCoeffs.p.toFixed(3)}, q=${vapEosCoeffs.q.toFixed(3)}, r=${vapEosCoeffs.r.toFixed(3)}, Roots=${ZV_all_roots?.map(r => r.toFixed(3))}`);

            ZL_K_loop = null; // Reset for this K-iteration
            if (ZL_all_roots && ZL_all_roots.length > 0) {
                const valid_ZL = ZL_all_roots.filter(z => z > liqEosCoeffs.B_mix_eos + MIN_MOLE_FRACTION); // Ensure Z > B_mix
                if (valid_ZL.length > 0) ZL_K_loop = valid_ZL[0]; // Smallest positive root > B_mix
            }

            ZV_K_loop = null; // Reset for this K-iteration
            if (ZV_all_roots && ZV_all_roots.length > 0) {
                const valid_ZV = ZV_all_roots.filter(z => z > vapEosCoeffs.B_mix_eos + MIN_MOLE_FRACTION); // Ensure Z > B_mix
                if (valid_ZV.length > 0) ZV_K_loop = valid_ZV[valid_ZV.length - 1]; // Largest positive root > B_mix
            }
            console.log(`    Selected ZL=${ZL_K_loop?.toFixed(4)}, ZV=${ZV_K_loop?.toFixed(4)} (B_liq=${liqEosCoeffs.B_mix_eos.toFixed(4)}, B_vap=${vapEosCoeffs.B_mix_eos.toFixed(4)})`);

            if (ZL_K_loop === null || ZV_K_loop === null) {
                console.warn(`  PR Bubble T (K-loop) ERR: No valid ZL/ZV. T_K=${T_K.toFixed(2)}, K_iter=${K_iter}`);
                break; 
            }

            const phi_L = calculatePrFugacityCoefficientsTernary(x, T_K, P_system_Pa, ZL_K_loop, liqMix.mixture_a, liqMix.mixture_b, liqMix.component_a_values, liqMix.component_b_values, ternaryPrKij);
            const phi_V = calculatePrFugacityCoefficientsTernary(y_vap, T_K, P_system_Pa, ZV_K_loop, vapMix.mixture_a, vapMix.mixture_b, vapMix.component_a_values, vapMix.component_b_values, ternaryPrKij);
            console.log(`    phi_L: [${phi_L.map(p => p?.toFixed(3) ?? 'NaN').join(', ')}], phi_V: [${phi_V.map(p => p?.toFixed(3) ?? 'NaN').join(', ')}]`);

            if (phi_L.some(isNaN) || phi_V.some(isNaN)) {
                console.error(`  PR Bubble T (K-loop) ERR: NaN in phi_L or phi_V. T_K=${T_K.toFixed(2)}`);
                ZL_K_loop = null; ZV_K_loop = null; // Mark as failed for outer loop check
                break; 
            }

            const new_K_values = phi_L.map((phiL_i, i) => phi_V[i] === 0 ? Infinity : phiL_i / phi_V[i]); // Corrected typo
             if (new_K_values.some(k => isNaN(k) || !isFinite(k))) {
                console.error(`  PR Bubble T (K-loop) ERR: NaN/Inf in new K_values. T_K=${T_K.toFixed(2)} K_values were: [${new_K_values.map(kv => kv?.toString()).join(',')}]`);
                 ZL_K_loop = null; ZV_K_loop = null; // Mark as failed
                 break;
            }
            console.log(`    New K_values: [${new_K_values.map(k => k.toFixed(3)).join(', ')}]`);

            const K_diff_sq_sum = K_values.reduce((sum, oldK, i) => sum + (new_K_values[i] - oldK)**2, 0);
            K_values = new_K_values;
            sum_Kx_from_K_loop = K_values.reduce((sum, Ki, i) => sum + Ki * x[i], 0);
            console.log(`    Sum Kx = ${sum_Kx_from_K_loop.toFixed(5)}, K_diff_sqrt = ${Math.sqrt(K_diff_sq_sum).toExponential(3)}`);

            if (Math.sqrt(K_diff_sq_sum) < K_VALUE_TOLERANCE * K_values.length && K_iter > 0) {
                console.log(`  K-loop converged at K_iter=${K_iter}.`);
                break; 
            }
            if (K_iter === MAX_K_VALUE_ITERATIONS - 1) {
                console.warn(`  PR Bubble T: K-value loop MAX iter at T=${T_K.toFixed(2)}K.`);
            }
        } // End K-value loop

        if (ZL_K_loop === null || ZV_K_loop === null) { // If K-loop broke due to ZL/ZV failure
            let T_K_change_factor = (T_iter % 2 === 0) ? 1.05 : 0.95; // Try nudging T
            console.warn(`  Adjusting T due to ZL/ZV failure in K-loop. Old T_K=${T_K.toFixed(2)}, Factor=${T_K_change_factor}`);
            T_K *= T_K_change_factor;
            if (T_K <= 100 || isNaN(T_K)) T_K = initialTempGuess_K * (0.8 + Math.random() * 0.4); // Ensure T_K is reasonable
             console.warn(`  New T_K for next T_iter = ${T_K.toFixed(2)}`);
            continue; 
        }

        const error_T_loop = sum_Kx_from_K_loop - 1.0;

        if (Math.abs(error_T_loop) < BUBBLE_T_CONVERGENCE_TOLERANCE) {
            const final_y_vap = normalizeMoleFractions(K_values.map((Ki, i) => Ki * x[i]));
            return { T_K, y_vap: final_y_vap, K_values, ZL: ZL_K_loop!, ZV: ZV_K_loop! };
        }

        // conservative damped temperature update
        const target_T_ratio = 1.0 / sum_Kx_from_K_loop;
        let damped_T_ratio_change = (target_T_ratio - 1.0) * (T_iter < 5 ? 0.1 : 0.25);
        let T_K_new_ratio = 1.0 + damped_T_ratio_change;
        T_K_new_ratio = Math.max(0.90, Math.min(1.10, T_K_new_ratio));
        const T_K_new = T_K * T_K_new_ratio;

        if (Math.abs(T_K_new - T_K) < BUBBLE_T_TOLERANCE_K && T_iter > 2) {
            if (Math.abs(error_T_loop) < BUBBLE_T_CONVERGENCE_TOLERANCE * 20) {
                const final_y_vap = normalizeMoleFractions(K_values.map((Ki, i) => Ki * x[i]));
                return { T_K, y_vap: final_y_vap, K_values, ZL: ZL_K_loop!, ZV: ZV_K_loop! };
            }
        }
        T_K = Math.max(150, Math.min(800, T_K_new));
        if (isNaN(T_K) || T_K <= 0) {
            T_K = initialTempGuess_K * (0.7 + Math.random() * 0.6);
            if (T_K <= 0) T_K = 273.15;
        }
    }
    console.warn(`PR Bubble T ERR: Max T-iterations reached for x=[${x.join(', ')}]. Last T_K=${T_K.toFixed(2)}, Last Sum Kx=${sum_Kx_from_K_loop.toFixed(4)}`);
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
): { dxdxi: number[], T_bubble_K: number } | null {
    const x_norm = normalizeMoleFractions(x_liq_current);

    const bubbleResult = findBubblePointTemperatureTernaryPr(x_norm, P_system_Pa, components, ternaryPrParams, initialTempGuess_K);
    if (!bubbleResult) return null;

    const { T_K: T_bubble_K, K_values } = bubbleResult;

    const y_vap = calculateVaporCompositionEos(x_norm, K_values);
    if (!y_vap) return null;

    const dxdxi = [x_norm[0] - y_vap[0], x_norm[1] - y_vap[1], x_norm[2] - y_vap[2]];
    return { dxdxi, T_bubble_K };
}

export async function simulateSingleResidueCurveODE_Pr(
    initial_x: number[], 
    P_system_Pa: number,
    componentsData: CompoundData[], 
    ternaryPrParams: TernaryPrParams, // Changed type to match internal usage
    d_xi_step: number,
    max_steps_per_direction: number,
    initialTemp_K: number
): Promise<ResidueCurveODE | null> {

    const run_direction = async (start_x: number[], forward: boolean, current_T_guess: number): Promise<ResidueCurveODE> => {
        const curve: ResidueCurveODE = [];
        let current_x = [...start_x];
        let last_T_K = current_T_guess;

        for (let loop_step = 0; loop_step < max_steps_per_direction; loop_step++) { // Renamed step to loop_step
            const ode_result = residue_curve_ode_dxdxi_pr(current_x, P_system_Pa, componentsData, ternaryPrParams, last_T_K);
            if (!ode_result) break;

            const { dxdxi, T_bubble_K } = ode_result;
            last_T_K = T_bubble_K;

            curve.push({ x: normalizeMoleFractions(current_x), T_K: T_bubble_K, step: (forward ? 1 : -1) * loop_step * d_xi_step }); // Use loop_step for 'step'

            const next_x_unnormalized = [
                current_x[0] + (forward ? 1 : -1) * dxdxi[0] * d_xi_step,
                current_x[1] + (forward ? 1 : -1) * dxdxi[1] * d_xi_step,
                current_x[2] + (forward ? 1 : -1) * dxdxi[2] * d_xi_step,
            ];
            current_x = normalizeMoleFractions(next_x_unnormalized);

            if (current_x.some(val => val >= 1.0 - MIN_MOLE_FRACTION * 2) || current_x.some(val => val <= MIN_MOLE_FRACTION * 2)) { // Use val
                curve.push({ x: normalizeMoleFractions(current_x), T_K: T_bubble_K, step: (forward ? 1 : -1) * (loop_step + 1) * d_xi_step }); // Use loop_step
                break;
            }
            if (loop_step > 0) { // prevent stagnation
                const prev_x_in_curve = curve[curve.length - 1].x;
                const diff_sq = prev_x_in_curve.reduce((sum, px, i) => sum + (px - current_x[i]) ** 2, 0);
                if (Math.sqrt(diff_sq) < d_xi_step * 1e-4) { // was 1e-3
                    break;
                }
            }
        }
        return curve;
    };

    const initial_x_norm = normalizeMoleFractions(initial_x);
    const initialBubbleResult = findBubblePointTemperatureTernaryPr(initial_x_norm, P_system_Pa, componentsData, ternaryPrParams, initialTemp_K);
    if (!initialBubbleResult) return null;

    const T_start_K = initialBubbleResult.T_K;

    const curve_forward = await run_direction(initial_x_norm, true, T_start_K);
    const curve_backward = await run_direction(initial_x_norm, false, T_start_K);

    const combined_curve: ResidueCurveODE = [...curve_backward.reverse().slice(1), ...curve_forward];

    return combined_curve.length < 2 ? (curve_forward.length > 0 ? curve_forward : (curve_backward.length > 0 ? curve_backward.reverse() : null)) : combined_curve;
}

