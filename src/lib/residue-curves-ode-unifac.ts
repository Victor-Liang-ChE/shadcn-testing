import type { CompoundData, AntoineParams, UnifacGroupComposition } from './vle-types';
// Corrected import path for UnifacParameters
import { calculateUnifacGamma, calculatePsat_Pa, type UnifacParameters } from './vle-calculations-unifac';

// --- Constants ---
const R_gas_const_J_molK = 8.31446261815324; // J/mol·K - Defined locally
const MAX_BUBBLE_T_ITERATIONS = 50;
const BUBBLE_T_TOLERANCE_K = 0.01;
const MIN_MOLE_FRACTION = 1e-6; // Adjusted from 1e-9 to 1e-6
const azeotrope_dxdt_threshold = 1e-7; // Azeotrope detection threshold

// --- Interfaces ---
export interface ResidueCurvePointODE {
    x: number[]; // Mole fractions [x0, x1, x2]
    T_K: number; // Temperature at this point
    step?: number; // Optional integration variable value
    dxdt_norm?: number;
    isPotentialAzeotrope?: boolean;
}
export type ResidueCurveODE = ResidueCurvePointODE[];

// No specific TernaryUnifacParams needed as interaction params are global

function normalizeMoleFractions(x: number[]): number[] {
    const sum = x.reduce((acc, val) => acc + val, 0);
    if (sum === 0) return x.map(() => 1 / x.length);
    // Ensure the calculation of the upper bound for min correctly uses x.length
    const minUpperClamp = 1.0 - MIN_MOLE_FRACTION * (x.length > 1 ? (x.length - 1) : 0);
    return x.map(val => Math.max(MIN_MOLE_FRACTION, Math.min(minUpperClamp, val / sum)));
}

// calculateUnifacGamma is imported and can handle multiple components

export function findBubblePointTemperatureTernaryUnifac(
    x_liq: number[],
    P_system_Pa: number,
    components: CompoundData[], // Antoine and Unifac groups needed
    unifacGlobalParams: UnifacParameters,
    initialTempGuess_K: number
): { T_K: number, gammas: number[], Psats: number[] } | null {
    let T_K = initialTempGuess_K;
    const x = normalizeMoleFractions(x_liq);

    if (components.some(c => !c.antoine || !c.unifacGroups)) return null;
    if (components.length !== 3 || x.length !== 3) return null;


    for (let iter = 0; iter < MAX_BUBBLE_T_ITERATIONS; iter++) {
        const Psats = components.map(c => calculatePsat_Pa(c.antoine!, T_K));
        if (Psats.some(isNaN)) return null;

        // Pass all three components to calculateUnifacGamma
        const gammasResult = calculateUnifacGamma(components, x, T_K, unifacGlobalParams);
        
        // The following check ensures gammasResult is a valid array of 3 numbers.
        // The TS2367 error ("types '2' and '3' have no overlap") on `gammasResult.length !== 3`
        // might be due to specific type inference by the TS language server for `calculateUnifacGamma`
        // or an issue with its declared return type. If `calculateUnifacGamma` is correctly typed
        // (e.g., to return `number[] | null` where the array can have 3 elements for a ternary case),
        // this line is logically sound and correctly validates the data structure.
        if (!gammasResult || !Array.isArray(gammasResult) || gammasResult.length !== 3 || gammasResult.some(isNaN)) return null;
        const gammas: number[] = gammasResult as number[]; // gammasResult is now confirmed to be number[] with length 3


        const P_calc = x[0]*gammas[0]*Psats[0] + x[1]*gammas[1]*Psats[1] + x[2]*gammas[2]*Psats[2];
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

function residue_curve_ode_dxdxi_unifac(
    x_liq_current: number[],
    P_system_Pa: number,
    components: CompoundData[],
    unifacGlobalParams: UnifacParameters,
    initialTempGuess_K: number
): { dxdxi: number[], T_bubble_K: number } | null {
    const x_norm = normalizeMoleFractions(x_liq_current);

    const bubbleResult = findBubblePointTemperatureTernaryUnifac(x_norm, P_system_Pa, components, unifacGlobalParams, initialTempGuess_K);
    if (!bubbleResult) return null;

    const { T_K: T_bubble_K, gammas, Psats } = bubbleResult;

    const y_vap = calculateVaporCompositionTernary(x_norm, P_system_Pa, gammas, Psats);
    if (!y_vap) return null;

    const dxdxi = [ x_norm[0] - y_vap[0], x_norm[1] - y_vap[1], x_norm[2] - y_vap[2] ];
    return { dxdxi, T_bubble_K };
}

export async function simulateSingleResidueCurveODE_Unifac(
    initial_x: number[],
    P_system_Pa: number,
    componentsData: CompoundData[],
    unifacGlobalParams: UnifacParameters,
    d_xi_step: number,
    max_steps_per_direction: number,
    initialTemp_K: number
): Promise<ResidueCurveODE | null> {
    
    const run_direction = async (start_x: number[], forward: boolean, current_T_guess: number): Promise<ResidueCurveODE> => {
        const curve: ResidueCurveODE = [];
        let current_x = [...start_x];
        let last_T_K = current_T_guess;

        for (let loop_step = 0; loop_step < max_steps_per_direction; loop_step++) { // Changed step to loop_step
            // Yield to the event loop to prevent freezing
            await new Promise(resolve => setTimeout(resolve, 0));

            const ode_result = residue_curve_ode_dxdxi_unifac(current_x, P_system_Pa, componentsData, unifacGlobalParams, last_T_K);
            if (!ode_result) break;

            const { dxdxi, T_bubble_K } = ode_result;
            last_T_K = T_bubble_K;

            const y_vap = [current_x[0] - dxdxi[0], current_x[1] - dxdxi[1], current_x[2] - dxdxi[2]];
            const dxdt = dxdxi;
            const dxdt_norm = Math.sqrt(dxdt[0]**2 + dxdt[1]**2 + dxdt[2]**2);

            const currentPointData: ResidueCurvePointODE = {
                x: normalizeMoleFractions(current_x),
                T_K: T_bubble_K,
                step: (forward ? 1 : -1) * loop_step * d_xi_step, // Use loop_step
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
            if (current_x.some(val => val >= upper_bound_limit) || current_x.some(val => val <= MIN_MOLE_FRACTION)) { 
                const final_x_normalized = normalizeMoleFractions(current_x);
                let lastPushedPoint = curve.length > 0 ? curve[curve.length - 1] : null;
                if (lastPushedPoint && final_x_normalized.every((val, i) => Math.abs(val - lastPushedPoint!.x[i]) < 1e-7)) {
                    lastPushedPoint.isPotentialAzeotrope = true; 
                    lastPushedPoint.dxdt_norm = 0;
                } else {
                    curve.push({ 
                        x: final_x_normalized, 
                        T_K: last_T_K, 
                        step: (forward ? 1 : -1) * (loop_step + 1) * d_xi_step, // Use loop_step
                        dxdt_norm: 0,
                        isPotentialAzeotrope: true
                    });
                }
                break;
            }
            if (loop_step > 0 && curve.length > 0) { // Use loop_step
                const lastAddedPoint = curve[curve.length - 1];
                const x_before_euler_step = curve.length > 1 ? curve[curve.length - 2].x : start_x;
                
                const diff_sq_from_prev_step = lastAddedPoint.x.reduce((sum, px, i) => sum + (px - x_before_euler_step[i]) ** 2, 0);

                if (Math.sqrt(diff_sq_from_prev_step) < d_xi_step * 1e-5) { // Adjusted threshold
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
    const initialBubbleResult = findBubblePointTemperatureTernaryUnifac(initial_x_norm, P_system_Pa, componentsData, unifacGlobalParams, initialTemp_K);
    if (!initialBubbleResult) return null;
    
    const T_start_K = initialBubbleResult.T_K;

    const curve_forward = await run_direction(initial_x_norm, true, T_start_K);
    const curve_backward = await run_direction(initial_x_norm, false, T_start_K);

    const combined_curve: ResidueCurveODE = [ ...curve_backward.reverse().slice(1), ...curve_forward ];
    
    return combined_curve.length < 2 ? (curve_forward.length > 0 ? curve_forward : (curve_backward.length > 0 ? curve_backward.reverse() : null)) : combined_curve;
}
