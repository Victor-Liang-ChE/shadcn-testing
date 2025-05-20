import type { CompoundData, AntoineParams } from './vle-types';
// calculatePsat_Pa is needed. NrtlInteractionParams and calculateNrtlGamma (binary) are not directly used by the new ternary logic.
import { calculatePsat_Pa } from './vle-calculations-nrtl';

// --- Constants ---
const R_gas_const_J_molK = 8.31446261815324; // J/mol·K
const MAX_BUBBLE_T_ITERATIONS = 50;
const BUBBLE_T_TOLERANCE_K = 0.01;
const MIN_MOLE_FRACTION = 1e-6; // Adjusted from 1e-9 to stabilize calculations
const azeotrope_dxdt_threshold = 1e-7; // Azeotrope detection threshold

// --- Interfaces ---
export interface ResidueCurvePointODE {
    x: number[]; // Mole fractions [x0, x1, x2]
    T_K: number; // Temperature at this point
    step?: number; // Optional integration variable value, changed from xi
    dxdt_norm?: number;
    isPotentialAzeotrope?: boolean;
}
export type ResidueCurveODE = ResidueCurvePointODE[];

export interface AzeotropeResult {
    x: number[]; // mole fractions [x0, x1, x2]
    T_K: number; // temperature
    iterations: number;
    converged: boolean;
    errorNorm?: number; // final norm of ||y-x||
}

export interface TernaryNrtlParams {
    // Assuming g_ij_J_mol and alpha_ij for pairs 0-1, 0-2, 1-2
    // g_ij means interaction from i to j.
    g01_J_mol: number; g10_J_mol: number; alpha01: number; // Pair (0,1)
    g02_J_mol: number; g20_J_mol: number; alpha02: number; // Pair (0,2)
    g12_J_mol: number; g21_J_mol: number; alpha12: number; // Pair (1,2)
}

function normalizeMoleFractions(x: number[]): number[] {
    const sum = x.reduce((acc, val) => acc + val, 0);
    if (sum === 0) return x.map(() => 1 / x.length);
    // For a 3-component system, (x.length - 1) is 2.
    const minUpperClamp = 1.0 - MIN_MOLE_FRACTION * (x.length > 1 ? (x.length - 1) : 0);
    return x.map(val => Math.max(MIN_MOLE_FRACTION, Math.min(minUpperClamp, val / sum)));
}

export function findBubblePointTemperatureTernaryNrtl(
    x_liq: number[], // mole fractions [x0, x1, x2]
    P_system_Pa: number,
    components: CompoundData[], // Antoine params needed
    ternaryNrtlParams: TernaryNrtlParams,
    initialTempGuess_K: number
): { T_K: number, gammas: number[], Psats: number[] } | null {
    let T_K = initialTempGuess_K;
    const x = normalizeMoleFractions(x_liq); // x is an array [x0, x1, x2]

    if (components.some(c => !c.antoine)) return null;
    if (components.length !== 3 || x.length !== 3) return null;

    // Prepare g and alpha matrices from ternaryNrtlParams
    // g[i][j] is g_ij (interaction energy parameter for i-j pair, g_ii = 0)
    // alpha[i][j] is alpha_ij (non-randomness factor for i-j pair, alpha_ii = 0, alpha_ij = alpha_ji)
    const g_J_mol = [
        [0, ternaryNrtlParams.g01_J_mol, ternaryNrtlParams.g02_J_mol],
        [ternaryNrtlParams.g10_J_mol, 0, ternaryNrtlParams.g12_J_mol],
        [ternaryNrtlParams.g20_J_mol, ternaryNrtlParams.g21_J_mol, 0]
    ];
    const alpha = [
        [0, ternaryNrtlParams.alpha01, ternaryNrtlParams.alpha02],
        [ternaryNrtlParams.alpha01, 0, ternaryNrtlParams.alpha12],
        [ternaryNrtlParams.alpha02, ternaryNrtlParams.alpha12, 0]
    ];

    for (let iter = 0; iter < MAX_BUBBLE_T_ITERATIONS; iter++) {
        const Psats = components.map(c => calculatePsat_Pa(c.antoine!, T_K));
        if (Psats.some(isNaN)) return null;

        // Calculate tau and G matrices based on current T_K
        const tau = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        const G = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        const RT = R_gas_const_J_molK * T_K;

        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                if (i === j) continue;
                tau[i][j] = g_J_mol[i][j] / RT;
                G[i][j] = Math.exp(-alpha[i][j] * tau[i][j]);
            }
        }

        const gammas = [0, 0, 0];
        for (let i = 0; i < 3; i++) { // For each component i
            let term1_num = 0;
            let term1_den = 0;
            for (let j = 0; j < 3; j++) { // Sum over j
                term1_num += x[j] * tau[j][i] * G[j][i]; // tau_ji * G_ji
                term1_den += x[j] * G[j][i];             // G_ji
            }
            const term1 = term1_den === 0 ? 0 : term1_num / term1_den;

            let term2_sum = 0;
            for (let j = 0; j < 3; j++) { // Sum over j
                if (term1_den === 0) continue; // Avoid issues if G_ki are all zero for some k,i pair

                let factor1_num = x[j] * G[i][j]; // x_j * G_ij
                let factor1_den = 0;
                for (let k = 0; k < 3; k++) { // Sum over k for denominator of factor1
                    factor1_den += x[k] * G[k][j]; // x_k * G_kj
                }
                const factor1 = factor1_den === 0 ? 0 : factor1_num / factor1_den;

                let factor2_num = 0;
                for (let m = 0; m < 3; m++) { // Sum over m for numerator of factor2
                    factor2_num += x[m] * tau[m][j] * G[m][j]; // x_m * tau_mj * G_mj
                }
                // Denominator for factor2 is the same as factor1_den
                const factor2 = factor1_den === 0 ? 0 : factor2_num / factor1_den;

                term2_sum += factor1 * (tau[i][j] - factor2); // tau_ij
            }
            const ln_gamma_i = term1 + term2_sum;
            gammas[i] = Math.exp(ln_gamma_i);
            if (isNaN(gammas[i]) || !isFinite(gammas[i])) {
                return null;
            }
        }

        const P_calc = x[0] * gammas[0] * Psats[0] + x[1] * gammas[1] * Psats[1] + x[2] * gammas[2] * Psats[2];
        const error_P = P_calc - P_system_Pa;

        if (Math.abs(error_P / P_system_Pa) < 1e-5 || Math.abs(error_P) < 0.1) {
            return { T_K, gammas, Psats };
        }

        let T_change = -error_P / P_system_Pa * (T_K * 0.05);
        const max_T_step = 5.0;
        if (Math.abs(T_change) > max_T_step) T_change = Math.sign(T_change) * max_T_step;
        if (iter < 5) T_change *= 0.5;

        let T_K_new = T_K + T_change;
        if (T_K_new <= 0) T_K_new = T_K * 0.9;
        if (Math.abs(T_K - T_K_new) < BUBBLE_T_TOLERANCE_K && iter > 3) {
            if (Math.abs(error_P / P_system_Pa) < 5e-5) return { T_K, gammas, Psats };
        }
        T_K = T_K_new;
    }
    return null;
}

function calculateVaporCompositionTernary(
    x_liq: number[], P_system_Pa: number, gammas: number[], Psats: number[]
): number[] | null {
    const y_unnormalized = [
        x_liq[0] * gammas[0] * Psats[0] / P_system_Pa,
        x_liq[1] * gammas[1] * Psats[1] / P_system_Pa,
        x_liq[2] * gammas[2] * Psats[2] / P_system_Pa,
    ];
    return normalizeMoleFractions(y_unnormalized);
}

function residue_curve_ode_dxdxi_nrtl(
    x_liq_current: number[],
    P_system_Pa: number,
    components: CompoundData[],
    ternaryNrtlParams: TernaryNrtlParams,
    initialTempGuess_K: number
): { dxdxi: number[], T_bubble_K: number } | null {
    const x_norm = normalizeMoleFractions(x_liq_current);

    const bubbleResult = findBubblePointTemperatureTernaryNrtl(x_norm, P_system_Pa, components, ternaryNrtlParams, initialTempGuess_K);
    if (!bubbleResult) return null;

    const { T_K: T_bubble_K, gammas, Psats } = bubbleResult;

    const y_vap = calculateVaporCompositionTernary(x_norm, P_system_Pa, gammas, Psats);
    if (!y_vap) return null;

    const dxdxi = [x_norm[0] - y_vap[0], x_norm[1] - y_vap[1], x_norm[2] - y_vap[2]];
    return { dxdxi, T_bubble_K };
}

export async function simulateSingleResidueCurveODE_Nrtl(
    initial_x: number[],
    P_system_Pa: number,
    componentsData: CompoundData[],
    ternaryNrtlParams: TernaryNrtlParams,
    d_xi_step: number,
    max_steps_per_direction: number,
    initialTemp_K: number
): Promise<ResidueCurveODE | null> {

    const run_direction = async (start_x: number[], forward: boolean, current_T_guess: number): Promise<ResidueCurveODE> => {
        const curve: ResidueCurveODE = [];
        let current_x = [...start_x];
        let last_T_K = current_T_guess;

        for (let loop_step = 0; loop_step < max_steps_per_direction; loop_step++) { 
            const ode_result = residue_curve_ode_dxdxi_nrtl(current_x, P_system_Pa, componentsData, ternaryNrtlParams, last_T_K);
            if (!ode_result) break;

            const { dxdxi, T_bubble_K } = ode_result;
            last_T_K = T_bubble_K;

            const y_vap = [current_x[0] - dxdxi[0], current_x[1] - dxdxi[1], current_x[2] - dxdxi[2]];
            const dxdt = dxdxi;
            const dxdt_norm = Math.sqrt(dxdt[0]**2 + dxdt[1]**2 + dxdt[2]**2);
            
            const currentPointData: ResidueCurvePointODE = {
                x: normalizeMoleFractions(current_x),
                T_K: T_bubble_K,
                step: (forward ? 1 : -1) * loop_step * d_xi_step, 
                dxdt_norm,
                isPotentialAzeotrope: false // Initialize to false
            };

            if (curve.length > 0) {
                const lastPoint = curve[curve.length - 1];
                const x_diff_sq = lastPoint.x.reduce((sum, val, i) => sum + (val - currentPointData.x[i])**2, 0);
                if (Math.sqrt(x_diff_sq) < 1e-7 && Math.abs(lastPoint.T_K - currentPointData.T_K) < 1e-4) {
                    if (dxdt_norm < azeotrope_dxdt_threshold) { 
                        lastPoint.isPotentialAzeotrope = true; 
                        lastPoint.dxdt_norm = Math.min(lastPoint.dxdt_norm ?? Infinity, dxdt_norm);
                        break; 
                    }
                } else {
                    curve.push(currentPointData);
                }
            } else {
                curve.push(currentPointData);
            }

            if (dxdt_norm < azeotrope_dxdt_threshold) {
                if (curve.length > 0) {
                    curve[curve.length-1].isPotentialAzeotrope = true;
                }
                break; 
            }

            const next_x_unnormalized = [
                current_x[0] + (forward ? 1 : -1) * dxdxi[0] * d_xi_step,
                current_x[1] + (forward ? 1 : -1) * dxdxi[1] * d_xi_step,
                current_x[2] + (forward ? 1 : -1) * dxdxi[2] * d_xi_step,
            ];
            const prev_current_x_for_stagnation_check = [...current_x];
            current_x = normalizeMoleFractions(next_x_unnormalized);

            const upper_bound_limit = 1.0 - MIN_MOLE_FRACTION * (componentsData.length > 1 ? (componentsData.length - 1) : 0);
            if (current_x.some(val => val >= upper_bound_limit ) || current_x.some(val => val <= MIN_MOLE_FRACTION)) { 
                const final_x_normalized = normalizeMoleFractions(current_x);
                let lastPushedPoint = curve.length > 0 ? curve[curve.length - 1] : null;
                if (lastPushedPoint && final_x_normalized.every((val, i) => Math.abs(val - lastPushedPoint!.x[i]) < 1e-7)) {
                    lastPushedPoint.isPotentialAzeotrope = true; 
                    lastPushedPoint.dxdt_norm = 0;
                } else {
                    curve.push({ 
                        x: final_x_normalized, 
                        T_K: last_T_K, 
                        step: (forward ? 1 : -1) * (loop_step + 1) * d_xi_step, 
                        dxdt_norm: 0,
                        isPotentialAzeotrope: true
                    });
                }
                break;
            }
            if (loop_step > 0 && curve.length > 0) { 
                const lastAddedPoint = curve[curve.length - 1];
                const x_before_euler_step = curve.length > 1 ? curve[curve.length - 2].x : start_x;
                
                const diff_sq_from_prev_step = lastAddedPoint.x.reduce((sum, px, i) => sum + (px - x_before_euler_step[i]) ** 2, 0);
                if (Math.sqrt(diff_sq_from_prev_step) < d_xi_step * 1e-5) {
                     if (lastAddedPoint.dxdt_norm !== undefined && lastAddedPoint.dxdt_norm < azeotrope_dxdt_threshold) {
                        lastAddedPoint.isPotentialAzeotrope = true;
                    }
                    break; 
                }
            }
        }
        if (curve.length > 0) {
            const lastCurvePoint = curve[curve.length - 1];
            if (lastCurvePoint.dxdt_norm !== undefined && lastCurvePoint.dxdt_norm < azeotrope_dxdt_threshold && !lastCurvePoint.isPotentialAzeotrope) {
                lastCurvePoint.isPotentialAzeotrope = true;
            }
        }
        return curve;
    };

    const initial_x_norm = normalizeMoleFractions(initial_x);
    const initialBubbleResult = findBubblePointTemperatureTernaryNrtl(initial_x_norm, P_system_Pa, componentsData, ternaryNrtlParams, initialTemp_K);
    if (!initialBubbleResult) return null;

    const T_start_K = initialBubbleResult.T_K;

    const curve_forward = await run_direction(initial_x_norm, true, T_start_K);
    const curve_backward = await run_direction(initial_x_norm, false, T_start_K);

    const combined_curve: ResidueCurveODE = [...curve_backward.reverse().slice(1), ...curve_forward];

    return combined_curve.length < 2 ? (curve_forward.length > 0 ? curve_forward : (curve_backward.length > 0 ? curve_backward.reverse() : null)) : combined_curve;
}

export function findAzeotropeNRTL(
    P_system_Pa: number,
    components: CompoundData[],
    ternaryNrtlParams: TernaryNrtlParams,
    initialGuess_x: number[], // Initial guess for mole fractions [x0, x1, x2]
    initialGuess_T_K: number, // Initial guess for temperature
    maxIterations: number = 100,
    tolerance: number = 1e-7, // Tolerance for ||y-x|| norm
    dampingFactor: number = 0.5 // Damping factor for x updates (0 < dampingFactor <= 1)
): AzeotropeResult | null {

    if (components.length !== 3) {
        console.error("Azeotrope search: Requires 3 components.");
        return null;
    }
    if (initialGuess_x.length !== 3) {
        console.error("Azeotrope search: Initial guess for x must have 3 components.");
        return null;
    }

    let current_x = normalizeMoleFractions([...initialGuess_x]);
    let current_T_K = initialGuess_T_K;
    let last_errorNorm = Infinity;

    for (let iter = 0; iter < maxIterations; iter++) {
        // 1. Calculate bubble point T for current_x at P_system_Pa
        // This function (findBubblePointTemperatureTernaryNrtl) needs to be robust.
        // It internally uses normalizeMoleFractions for its input x_liq.
        const bubbleResult = findBubblePointTemperatureTernaryNrtl(
            current_x,
            P_system_Pa,
            components,
            ternaryNrtlParams,
            current_T_K // Use current T as guess for next bubble T calc
        );

        if (!bubbleResult) {
            console.warn(`Azeotrope search (iter ${iter}): Bubble point calculation failed for x = [${current_x.map(v => v.toFixed(4)).join(', ')}]`);
            return { x: current_x, T_K: current_T_K, iterations: iter, converged: false, errorNorm: last_errorNorm };
        }

        current_T_K = bubbleResult.T_K;
        const gammas = bubbleResult.gammas;
        const Psats = bubbleResult.Psats;

        // 2. Calculate vapor composition y for this (current_x, current_T_K)
        const y_vap = calculateVaporCompositionTernary(current_x, P_system_Pa, gammas, Psats);

        if (!y_vap) {
            console.warn(`Azeotrope search (iter ${iter}): Vapor composition calculation failed.`);
            return { x: current_x, T_K: current_T_K, iterations: iter, converged: false, errorNorm: last_errorNorm };
        }

        // 3. Check for convergence: ||y - x|| < tolerance
        const errorVector = [
            y_vap[0] - current_x[0],
            y_vap[1] - current_x[1],
            y_vap[2] - current_x[2]
        ];
        const errorNorm = Math.sqrt(
            errorVector[0]**2 + errorVector[1]**2 + errorVector[2]**2
        );

        last_errorNorm = errorNorm;

        if (errorNorm < tolerance) {
            return { x: normalizeMoleFractions(current_x), T_K: current_T_K, iterations: iter, converged: true, errorNorm };
        }

        // 4. Update x for next iteration (damped fixed-point iteration: x_new = x_old + lambda * (y_old - x_old))
        let next_x_unnormalized = [
            current_x[0] + dampingFactor * errorVector[0],
            current_x[1] + dampingFactor * errorVector[1],
            current_x[2] + dampingFactor * errorVector[2],
        ];
        
        const previous_x = [...current_x];
        current_x = normalizeMoleFractions(next_x_unnormalized);

        // Optional: Check for stagnation or very small change in x
        const x_change_norm = Math.sqrt(
            (current_x[0] - previous_x[0])**2 +
            (current_x[1] - previous_x[1])**2 +
            (current_x[2] - previous_x[2])**2
        );
        if (x_change_norm < tolerance * 1e-1 && iter > 5) { // Stricter tolerance for x_change
            // console.log(`Azeotrope search: x update step too small at iter ${iter}. Current errorNorm: ${errorNorm.toExponential(3)}`);
            // If errorNorm is also small but not quite met, might be close enough.
            // Otherwise, it might be stuck.
            if (errorNorm < tolerance * 10) { // Looser condition if x is stuck but error is low
                return { x: current_x, T_K: current_T_K, iterations: iter, converged: true, errorNorm };
            }
            break; // Or break if truly stuck with high error
        }
    }

    console.warn(`Azeotrope search: Did not converge after ${maxIterations} iterations. Final ||y-x||: ${last_errorNorm.toExponential(3)}`);
    return { x: current_x, T_K: current_T_K, iterations: maxIterations, converged: false, errorNorm: last_errorNorm };
}
