import type { CompoundData, AntoineParams, UnifacGroupComposition } from './vle-types';
// Corrected import path for UnifacParameters
import { calculateUnifacGamma, calculatePsat_Pa, type UnifacParameters } from './vle-calculations-unifac';

// --- Constants ---
// R_gas_const_J_molK = 8.31446261815324; // J/mol·K - Defined locally // Not used locally, assuming imported if needed elsewhere
const MAX_BUBBLE_T_ITERATIONS_UNIFAC = 35; // Reduced
const BUBBLE_T_TOLERANCE_K_UNIFAC = 0.02;  // Slightly looser
const MIN_MOLE_FRACTION_UNIFAC = 1e-7; // Adjusted slightly
const AZEOTROPE_DXDT_THRESHOLD_UNIFAC = 1e-6; // More sensitive for azeotrope detection

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

function normalizeMoleFractions(x: number[]): number[] { // Using a local normalize if specific MIN_MOLE_FRACTION is needed
    const sum = x.reduce((acc, val) => acc + val, 0);
    if (sum === 0) return x.map(() => 1 / x.length);
    // Ensure the calculation of the upper bound for min correctly uses x.length
    const minUpperClamp = 1.0 - MIN_MOLE_FRACTION_UNIFAC * (x.length > 1 ? (x.length - 1) : 0);
    return x.map(val => Math.max(MIN_MOLE_FRACTION_UNIFAC, Math.min(minUpperClamp, val / sum)));
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

    if (components.some(c => !c.antoine || !c.unifacGroups) || components.length !== 3 || x.length !== 3) return null;


    for (let iter = 0; iter < MAX_BUBBLE_T_ITERATIONS_UNIFAC; iter++) {
        const Psats = components.map(c => calculatePsat_Pa(c.antoine!, T_K));
        if (Psats.some(p => isNaN(p) || p <= 0)) { // Added check for p <= 0
            // Try a small T adjustment if Psat fails, might be out of Antoine range
            T_K *= (iter % 2 === 0 ? 1.005 : 0.995); 
            T_K = Math.max(150, Math.min(700, T_K)); // Bound T
            if (iter > 5) return null; // Give up if consistently failing
            continue;
        }

        // Pass all three components to calculateUnifacGamma
        const gammasResult = calculateUnifacGamma(components, x, T_K, unifacGlobalParams);
        
        if (!gammasResult || !Array.isArray(gammasResult) || gammasResult.length !== 3 || gammasResult.some(g => isNaN(g) || g <=0)) { // Added check for g <= 0
            // Similar T adjustment if gamma fails
            T_K *= (iter % 2 === 0 ? 1.005 : 0.995);
            T_K = Math.max(150, Math.min(700, T_K));
            if (iter > 5) return null;
            continue;
        }
        const gammas: number[] = gammasResult as number[]; // gammasResult is now confirmed to be number[] with length 3


        const P_calc = x[0]*gammas[0]*Psats[0] + x[1]*gammas[1]*Psats[1] + x[2]*gammas[2]*Psats[2];
        const error_P = P_calc - P_system_Pa;

        // Relative tolerance for P_calc vs P_system_Pa
        if (Math.abs(error_P / P_system_Pa) < 1e-4 || Math.abs(error_P) < 0.5 /* Pa absolute*/) { // Looser P tolerance
            return { T_K, gammas, Psats };
        }

        // More conservative T_change
        let T_change_factor = -error_P / P_system_Pa; // Proportional to relative error
        // Bound the factor to prevent excessive jumps
        const max_factor = 0.1; // Max 10% relative change in P implies proportional change in T
        T_change_factor = Math.max(-max_factor, Math.min(max_factor, T_change_factor));
        
        let T_change = T_change_factor * (T_K * 0.1); // Change T by up to 1% of current T, scaled by error factor
                                                     // Reduced from 0.05 to 0.1 for sensitivity scaling
        
        const max_abs_T_step = 2.0; // Max 2K change per step
        if (Math.abs(T_change) > max_abs_T_step) T_change = Math.sign(T_change) * max_abs_T_step;
        if (iter < 3) T_change *= 0.3; // Heavy damping for first few steps

        let T_K_new = T_K + T_change;
        if (T_K_new <= 0) T_K_new = T_K * 0.95; // Ensure positive
        T_K_new = Math.max(150, Math.min(700, T_K_new)); // Bound T

         if (Math.abs(T_K - T_K_new) < BUBBLE_T_TOLERANCE_K_UNIFAC && iter > 2) { // iter > 2 for stability
             if (Math.abs(error_P / P_system_Pa) < 5e-4) return { T_K, gammas, Psats }; // Looser final convergence
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
    if (y_unnormalized.some(isNaN)) return null;
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
    d_xi_step: number, // This is the step size passed from page.tsx
    max_steps_per_direction: number,
    initialTemp_K: number
): Promise<ResidueCurveODE | null> {
    
    const run_direction = async (start_x: number[], forward: boolean, current_T_guess: number): Promise<ResidueCurveODE> => {
        const curve: ResidueCurveODE = [];
        let current_x = normalizeMoleFractions(start_x); // Normalize at start of direction
        let last_T_K = current_T_guess;

        for (let loop_step = 0; loop_step < max_steps_per_direction; loop_step++) {
            // Conditional yield to event loop
            if (loop_step > 0 && loop_step % 20 === 0) { // Yield every 20 steps
                await new Promise(resolve => setTimeout(resolve, 0));
            }

            const ode_result = residue_curve_ode_dxdxi_unifac(current_x, P_system_Pa, componentsData, unifacGlobalParams, last_T_K);
            if (!ode_result) break;

            const { dxdxi, T_bubble_K } = ode_result;
            last_T_K = T_bubble_K;

            // const y_vap = [current_x[0] - dxdxi[0], current_x[1] - dxdxi[1], current_x[2] - dxdxi[2]]; // Not directly used later
            // const dxdt = dxdxi; // Redundant
            const dxdt_norm = Math.sqrt(dxdxi[0]**2 + dxdxi[1]**2 + dxdxi[2]**2);

            const currentPointData: ResidueCurvePointODE = {
                x: [...current_x], // Store a copy
                T_K: T_bubble_K,
                step: (forward ? 1 : -1) * loop_step, // Simpler step count
                dxdt_norm,
                isPotentialAzeotrope: false // Initialize to false
            };
            
            if (curve.length > 0) {
                const lastPoint = curve[curve.length - 1];
                const x_diff_sq = lastPoint.x.reduce((sum, val, i) => sum + (val - currentPointData.x[i])**2, 0);
                if (Math.sqrt(x_diff_sq) < 1e-8 && Math.abs(lastPoint.T_K - currentPointData.T_K) < 1e-5) { // Stricter check for adding point
                    if (dxdt_norm < AZEOTROPE_DXDT_THRESHOLD_UNIFAC) { 
                        lastPoint.isPotentialAzeotrope = true; 
                        lastPoint.dxdt_norm = Math.min(lastPoint.dxdt_norm ?? Infinity, dxdt_norm);
                    }
                    break; 
                } else {
                    curve.push(currentPointData);
                }
            } else {
                curve.push(currentPointData);
            }

            if (dxdt_norm < AZEOTROPE_DXDT_THRESHOLD_UNIFAC && loop_step > 2) { // Allow few steps before azeo check
                if (curve.length > 0) {
                    curve[curve.length-1].isPotentialAzeotrope = true;
                }
                break; 
            }
            
            const next_x_unnormalized = current_x.map((xi, idx) => xi + (forward ? 1 : -1) * dxdxi[idx] * d_xi_step);
            // const prev_current_x_for_stagnation_check = [...current_x]; // Not used
            current_x = normalizeMoleFractions(next_x_unnormalized);

            const numComponents = componentsData.length;
            const PURE_COMP_HIGH_THRESHOLD_UNIFAC = 1.0 - (numComponents - 1) * (MIN_MOLE_FRACTION_UNIFAC * 2.0);
            const PURE_COMP_LOW_THRESHOLD_UNIFAC = MIN_MOLE_FRACTION_UNIFAC * 2.0;

            for (let i = 0; i < numComponents; i++) {
                if (current_x[i] >= PURE_COMP_HIGH_THRESHOLD_UNIFAC) {
                    let others_are_tiny = true;
                    for (let j = 0; j < numComponents; j++) {
                        if (i !== j && current_x[j] > PURE_COMP_LOW_THRESHOLD_UNIFAC * 100) { // Increased multiplier for stricter "tiny" check
                            others_are_tiny = false;
                            break;
                        }
                    }
                    if (others_are_tiny) {
                        const final_pure_x = Array(numComponents).fill(0.0).map((_,idx)=> idx === i ? 1.0 : 0.0);
                        const normalized_final_pure_x = normalizeMoleFractions(final_pure_x); // Ensure it's normalized
                         if (curve.length > 0) {
                             curve[curve.length - 1].x = normalized_final_pure_x;
                             curve[curve.length - 1].isPotentialAzeotrope = true; // Reaching pure is a terminal point
                             curve[curve.length - 1].dxdt_norm = 0;
                         } else {
                             curve.push({x: normalized_final_pure_x, T_K: last_T_K, step: (forward ? 1 : -1) * (loop_step +1), dxdt_norm: 0, isPotentialAzeotrope: true});
                         }
                        return curve; // Exit run_direction
                    }
                }
            }
            // Stagnation if x doesn't change much
            if (loop_step > 0 && curve.length > 0) { // Use loop_step
                const lastAddedPoint = curve[curve.length - 1];
                const prev_x = loop_step > 1 && curve.length > 1 ? curve[curve.length - 2].x : start_x; // Use start_x if only one point in curve
                
                const actual_x_change_sq = lastAddedPoint.x.reduce((s,v,i)=> s + (v-prev_x[i])**2, 0);

                if (Math.sqrt(actual_x_change_sq) < d_xi_step * 1e-5 && dxdt_norm < AZEOTROPE_DXDT_THRESHOLD_UNIFAC * 10) { // Adjusted threshold
                    if (lastAddedPoint.dxdt_norm !== undefined && lastAddedPoint.dxdt_norm < AZEOTROPE_DXDT_THRESHOLD_UNIFAC) {
                        lastAddedPoint.isPotentialAzeotrope = true;
                    }
                    break; 
                }
            }
        }
        // Ensure last point's azeotrope status is set if norm is low
        if (curve.length > 0) {
            const lastCurvePoint = curve[curve.length - 1];
            if (lastCurvePoint.dxdt_norm !== undefined && lastCurvePoint.dxdt_norm < AZEOTROPE_DXDT_THRESHOLD_UNIFAC && !lastCurvePoint.isPotentialAzeotrope) {
                lastCurvePoint.isPotentialAzeotrope = true;
            }
        }
        return curve;
    };

    const initial_x_norm = normalizeMoleFractions(initial_x);
    const initialBubbleResult = findBubblePointTemperatureTernaryUnifac(initial_x_norm, P_system_Pa, componentsData, unifacGlobalParams, initialTemp_K);
    if (!initialBubbleResult) return null;
    
    const T_start_K = initialBubbleResult.T_K;

    // Check if initial point is critical
    let is_initial_point_critical = false;
    const initial_ode_check_val = residue_curve_ode_dxdxi_unifac(initial_x_norm, P_system_Pa, componentsData, unifacGlobalParams, T_start_K);
    if (initial_ode_check_val) {
        const dxdxi_norm_initial = Math.sqrt(initial_ode_check_val.dxdxi.reduce((s, v) => s + v*v, 0));
        if (dxdxi_norm_initial < AZEOTROPE_DXDT_THRESHOLD_UNIFAC) is_initial_point_critical = true;
    }
    const PURE_COMP_THRESH_INIT_CHECK_UNIFAC = 1.0 - (componentsData.length - 1) * (MIN_MOLE_FRACTION_UNIFAC * 10.0);
    for (let i = 0; i < componentsData.length; i++) {
        if (initial_x_norm[i] >= PURE_COMP_THRESH_INIT_CHECK_UNIFAC) {
            let others_tiny = true;
            for(let j=0; j<componentsData.length; ++j) {
                if (i!==j && initial_x_norm[j] > MIN_MOLE_FRACTION_UNIFAC*100) {
                    others_tiny = false;
                    break;
                }
            }
            if(others_tiny) {
                is_initial_point_critical = true;
                break;
            }
        }
    }


    const curve_forward = await run_direction(initial_x_norm, true, T_start_K);
    let curve_backward: ResidueCurveODE = [];

    if (!is_initial_point_critical || curve_forward.length <= 1) { // If initial is critical and forward found something, backward might be redundant
        curve_backward = await run_direction(initial_x_norm, false, T_start_K);
    }
    
    // Curve combination logic (similar to robust SRK/PR)
    const first_point_data: ResidueCurvePointODE = { 
        x: initial_x_norm, 
        T_K: T_start_K, 
        step: 0, 
        dxdt_norm: initial_ode_check_val ? Math.sqrt(initial_ode_check_val.dxdxi.reduce((s,v)=>s+v*v,0)) : undefined, 
        isPotentialAzeotrope: is_initial_point_critical
    };
    let combined_curve: ResidueCurveODE = [];
    const MIN_DIST_SQ_ADD_COMBINED = (1e-10)**2; // Very small distance to consider points identical

    // Process backward curve: reverse it, and remove the first point if it's identical to first_point_data (which is the starting point)
    const backward_processed = curve_backward.length > 0 ? curve_backward.reverse() : [];
    if (backward_processed.length > 0) {
        const lastBackwardX = backward_processed[backward_processed.length -1].x; // This is the point closest to start_x from backward run
        const firstPointX = first_point_data.x;
        const distSqToStart = lastBackwardX.reduce((s,v,idx) => s + (v-firstPointX[idx])**2, 0);
        if (distSqToStart < MIN_DIST_SQ_ADD_COMBINED) {
            combined_curve.push(...backward_processed.slice(0, -1)); // Exclude the last point (which is like start_x)
        } else {
            combined_curve.push(...backward_processed);
        }
    }
    
    // Add the initial starting point (first_point_data)
    let addFirstPointToCombined = true;
    if (combined_curve.length > 0) { // If backward curve added points
        const lastCombinedPtX = combined_curve[combined_curve.length - 1].x;
        const distSq = lastCombinedPtX.reduce((s,v,idx) => s + (v-first_point_data.x[idx])**2, 0);
        if (distSq < MIN_DIST_SQ_ADD_COMBINED) addFirstPointToCombined = false;
    }
    if (addFirstPointToCombined) {
        combined_curve.push(first_point_data);
    }
    
    // Process forward curve: remove the first point if it's identical to first_point_data
    const forward_processed = curve_forward.length > 0 ? curve_forward : [];
    if (forward_processed.length > 0) {
        const firstForwardX = forward_processed[0].x;
        const firstPointX = first_point_data.x;
        const distSqToStart = firstForwardX.reduce((s,v,idx) => s + (v-firstPointX[idx])**2, 0);
        if (distSqToStart < MIN_DIST_SQ_ADD_COMBINED) {
            combined_curve.push(...forward_processed.slice(1)); // Exclude the first point
        } else {
            combined_curve.push(...forward_processed);
        }
    }

    // If combined_curve is still empty (e.g. initial point was critical and both directions returned nothing/single point)
    // but we had a valid start point, add it.
    if (combined_curve.length === 0 && initialBubbleResult) {
        combined_curve.push({x: initial_x_norm, T_K: T_start_K, step: 0, dxdt_norm: first_point_data.dxdt_norm, isPotentialAzeotrope: first_point_data.isPotentialAzeotrope});
    }

    if (combined_curve.length === 0) return null;
    return combined_curve;
}
