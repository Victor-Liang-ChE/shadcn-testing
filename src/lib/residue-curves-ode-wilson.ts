import type { CompoundData, WilsonPureComponentParams, AntoineParams } from './vle-types'; // WilsonInteractionParams removed from here
import { type WilsonInteractionParams as BinaryWilsonInteractionParams, R_gas_const_J_molK as R } from './vle-calculations-wilson'; // Corrected import
import { calculatePsat_Pa } from './vle-calculations-unifac'; 

// --- Constants ---
const MAX_BUBBLE_T_ITERATIONS = 50;
const BUBBLE_T_TOLERANCE_K = 0.01; // Tolerance for bubble T convergence
const MIN_MOLE_FRACTION = 1e-9; // Clip mole fractions to avoid log(0) or division by zero
const azeotrope_dxdt_threshold = 1e-7; // Azeotrope detection threshold

// --- Interfaces ---
export interface ResidueCurvePointODE {
    x: number[]; // Mole fractions [x0, x1, x2]
    T_K: number; // Temperature in Kelvin
    P_Pa?: number; // Pressure in Pascals (optional, as it's fixed for the simulation but good to store)
    step?: number; // Step number in the ODE integration
    dxdt_norm?: number; // Norm of (x - y)
    isPotentialAzeotrope?: boolean; // Flag for potential azeotropes
}
export type ResidueCurveODE = ResidueCurvePointODE[];

export interface TernaryWilsonParams {
    // For pair comp0-comp1 (components[0] and components[1])
    a01_J_mol: number; // Interaction from 0 towards 1 (e.g., lambda_01 - lambda_00)
    a10_J_mol: number; // Interaction from 1 towards 0 (e.g., lambda_10 - lambda_11)
    // For pair comp0-comp2 (components[0] and components[2])
    a02_J_mol: number;
    a20_J_mol: number;
    // For pair comp1-comp2 (components[1] and components[2])
    a12_J_mol: number;
    a21_J_mol: number;
}

function normalizeMoleFractions(x: number[]): number[] {
    const sum = x.reduce((acc, val) => acc + val, 0);
    if (sum === 0) return x.map(() => 1 / x.length); // Avoid division by zero, distribute equally
    // Consistent upper bound calculation
    const upper_bound_val = 1.0 - MIN_MOLE_FRACTION * (x.length > 1 ? (x.length - 1) : 0);
    return x.map(val => Math.max(MIN_MOLE_FRACTION, Math.min(upper_bound_val, val / sum)));
}


function calculateWilsonLambdas(
    V_L_m3mol: number[], // Molar volumes [V0, V1, V2]
    T_K: number,
    params: TernaryWilsonParams
): { L01: number, L10: number, L02: number, L20: number, L12: number, L21: number } {
    const RT = R * T_K;
    // Lambda_ij = (Vj/Vi) * exp(-a_ij / RT)
    // Where a_ij is (lambda_ij - lambda_ii) type param
    return {
        L01: (V_L_m3mol[1] / V_L_m3mol[0]) * Math.exp(-params.a01_J_mol / RT),
        L10: (V_L_m3mol[0] / V_L_m3mol[1]) * Math.exp(-params.a10_J_mol / RT),
        L02: (V_L_m3mol[2] / V_L_m3mol[0]) * Math.exp(-params.a02_J_mol / RT),
        L20: (V_L_m3mol[0] / V_L_m3mol[2]) * Math.exp(-params.a20_J_mol / RT),
        L12: (V_L_m3mol[2] / V_L_m3mol[1]) * Math.exp(-params.a12_J_mol / RT),
        L21: (V_L_m3mol[1] / V_L_m3mol[2]) * Math.exp(-params.a21_J_mol / RT),
    };
}

export function calculateWilsonGammaTernary(
    x: number[], // mole fractions [x0, x1, x2]
    T_K: number,
    components: CompoundData[], // Must have wilsonParams.V_L_m3mol
    ternaryParams: TernaryWilsonParams
): number[] | null { // Returns [gamma0, gamma1, gamma2] or null on error
    if (x.length !== 3 || components.length !== 3) return null;
    if (components.some(c => !c.wilsonParams?.V_L_m3mol || c.wilsonParams.V_L_m3mol <= 0)) return null;

    const V_L = components.map(c => c.wilsonParams!.V_L_m3mol);
    const L = calculateWilsonLambdas(V_L, T_K, ternaryParams);

    const x0 = x[0], x1 = x[1], x2 = x[2];

    // Sum_j (x_j * Lambda_ij) terms
    const S0 = x0 + x1 * L.L01 + x2 * L.L02;
    const S1 = x0 * L.L10 + x1 + x2 * L.L12;
    const S2 = x0 * L.L20 + x1 * L.L21 + x2;

    if (S0 === 0 || S1 === 0 || S2 === 0) return null; // Avoid division by zero

    // ln gamma_i = 1 - ln(Si) - sum_k ( (xk * Lambda_ki) / Sk )
    const lngamma0 = 1 - Math.log(S0) - ( (x0 * 1    / S0) + (x1 * L.L10 / S1) + (x2 * L.L20 / S2) );
    const lngamma1 = 1 - Math.log(S1) - ( (x0 * L.L01 / S0) + (x1 * 1    / S1) + (x2 * L.L21 / S2) );
    const lngamma2 = 1 - Math.log(S2) - ( (x0 * L.L02 / S0) + (x1 * L.L12 / S1) + (x2 * 1    / S2) );
    
    return [Math.exp(lngamma0), Math.exp(lngamma1), Math.exp(lngamma2)];
}


export function findBubblePointTemperatureTernary(
    x_liq: number[], // [x0, x1, x2]
    P_system_Pa: number,
    components: CompoundData[], // Antoine and Wilson V_L params needed
    ternaryWilsonParams: TernaryWilsonParams,
    initialTempGuess_K: number
): { T_K: number, gammas: number[], Psats: number[] } | null {
    let T_K = initialTempGuess_K;
    const x = normalizeMoleFractions(x_liq);

    for (let iter = 0; iter < MAX_BUBBLE_T_ITERATIONS; iter++) {
        const Psats = components.map(c => calculatePsat_Pa(c.antoine!, T_K));
        if (Psats.some(isNaN)) return null;

        const gammas = calculateWilsonGammaTernary(x, T_K, components, ternaryWilsonParams);
        if (!gammas || gammas.some(isNaN)) return null;

        const P_calc = x[0]*gammas[0]*Psats[0] + x[1]*gammas[1]*Psats[1] + x[2]*gammas[2]*Psats[2];
        const error_P = P_calc - P_system_Pa;

        if (Math.abs(error_P / P_system_Pa) < 1e-5 || Math.abs(error_P) < 0.1 /*Pa abs tolerance*/) {
            return { T_K, gammas, Psats };
        }

        // Simple temperature adjustment - can be improved with a derivative
        let T_change = -error_P / P_system_Pa * (T_K * 0.05); // Proportional step
        const max_T_step = 5.0; // K
        if (Math.abs(T_change) > max_T_step) T_change = Math.sign(T_change) * max_T_step;
        if (iter < 5) T_change *= 0.5; // Smaller steps initially

        let T_K_new = T_K + T_change;
        if (T_K_new <= 0) T_K_new = T_K * 0.9; // Avoid non-physical temps
        if (Math.abs(T_K - T_K_new) < BUBBLE_T_TOLERANCE_K && iter > 3) {
             // If T change is too small, might be stuck or converged
             if (Math.abs(error_P / P_system_Pa) < 5e-5) return { T_K, gammas, Psats }; // Looser tolerance if T stuck
             // Could try a more aggressive jump or stop
        }
        T_K = T_K_new;
    }
    return null; // Failed to converge
}

function calculateVaporCompositionTernary(
    x_liq: number[], // [x0, x1, x2]
    P_system_Pa: number,
    gammas: number[], // [g0, g1, g2]
    Psats: number[]   // [Ps0, Ps1, Ps2]
): number[] | null {
    const P_total_check = x_liq[0]*gammas[0]*Psats[0] + x_liq[1]*gammas[1]*Psats[1] + x_liq[2]*gammas[2]*Psats[2];
    // if (Math.abs(P_total_check - P_system_Pa) / P_system_Pa > 0.01) { // Check if P_total is close to P_system
    //     console.warn("Vapor comp calc: P_total deviates significantly from P_system");
    //     // This indicates that T_bubble might not be accurate
    // }

    const y_unnormalized = [
        x_liq[0] * gammas[0] * Psats[0] / P_system_Pa,
        x_liq[1] * gammas[1] * Psats[1] / P_system_Pa,
        x_liq[2] * gammas[2] * Psats[2] / P_system_Pa,
    ];
    return normalizeMoleFractions(y_unnormalized);
}

export function residue_curve_ode_dxdxi(
    x_liq_current: number[], // [x0, x1, x2]
    P_system_Pa: number,
    components: CompoundData[],
    ternaryWilsonParams: TernaryWilsonParams,
    initialTempGuess_K: number // For speeding up bubble T calc
): { dxdxi: number[], T_bubble_K: number } | null {
    const x_norm = normalizeMoleFractions(x_liq_current);

    const bubbleResult = findBubblePointTemperatureTernary(x_norm, P_system_Pa, components, ternaryWilsonParams, initialTempGuess_K);
    if (!bubbleResult) return null;

    const { T_K: T_bubble_K, gammas, Psats } = bubbleResult;

    const y_vap = calculateVaporCompositionTernary(x_norm, P_system_Pa, gammas, Psats);
    if (!y_vap) return null;

    const dxdxi = [
        x_norm[0] - y_vap[0],
        x_norm[1] - y_vap[1],
        x_norm[2] - y_vap[2],
    ];
    return { dxdxi, T_bubble_K };
}

export async function simulateSingleResidueCurveODE(
    initial_x: number[], // [x0, x1, x2]
    P_system_Pa: number,
    componentsData: CompoundData[],
    ternaryWilsonParams: TernaryWilsonParams,
    d_xi_step: number,
    max_steps_per_direction: number,
    initialTemp_K: number // A reasonable starting guess for temperature
): Promise<ResidueCurveODE | null> {
    
    const run_direction = async (start_x: number[], forward: boolean, current_T_guess: number): Promise<ResidueCurveODE> => {
        const curve: ResidueCurveODE = [];
        let current_x = [...start_x];
        let last_T_K = current_T_guess;

        for (let loop_idx = 0; loop_idx < max_steps_per_direction; loop_idx++) {
            const ode_result = residue_curve_ode_dxdxi(current_x, P_system_Pa, componentsData, ternaryWilsonParams, last_T_K);
            if (!ode_result) break; 

            const { dxdxi, T_bubble_K } = ode_result;
            last_T_K = T_bubble_K; 

            const y_vap = [current_x[0] - dxdxi[0], current_x[1] - dxdxi[1], current_x[2] - dxdxi[2]];
            const dxdt = dxdxi; 
            const dxdt_norm = Math.sqrt(dxdt[0]**2 + dxdt[1]**2 + dxdt[2]**2);

            const currentPointData: ResidueCurvePointODE = {
                x: normalizeMoleFractions(current_x),
                T_K: T_bubble_K,
                step: (forward ? 1 : -1) * loop_idx * d_xi_step,
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
            
            // Euler step
            const next_x_unnormalized = [
                current_x[0] + (forward ? 1 : -1) * dxdxi[0] * d_xi_step,
                current_x[1] + (forward ? 1 : -1) * dxdxi[1] * d_xi_step,
                current_x[2] + (forward ? 1 : -1) * dxdxi[2] * d_xi_step,
            ];
            
            const prev_current_x_for_stagnation_check = [...current_x]; // Save before normalization for stagnation check
            current_x = normalizeMoleFractions(next_x_unnormalized);

            // Termination conditions
            const upper_bound_limit = 1.0 - MIN_MOLE_FRACTION * (componentsData.length > 1 ? (componentsData.length - 1) : 0);
            if (current_x.some(val => val >= upper_bound_limit) || current_x.some(val => val <= MIN_MOLE_FRACTION)) {
                const final_x_normalized = normalizeMoleFractions(current_x); // Ensure it's normalized
                let lastPushedPoint = curve.length > 0 ? curve[curve.length - 1] : null;
                if (lastPushedPoint && final_x_normalized.every((val, i) => Math.abs(val - lastPushedPoint!.x[i]) < 1e-7)) {
                    lastPushedPoint.isPotentialAzeotrope = true; 
                    lastPushedPoint.dxdt_norm = 0;
                } else {
                    curve.push({ 
                        x: final_x_normalized, 
                        T_K: last_T_K, // Use T from before this step
                        step: (forward ? 1 : -1) * (loop_idx + 1) * d_xi_step,
                        dxdt_norm: 0, // At pure component vertex
                        isPotentialAzeotrope: true // Pure components are nodes
                    });
                }
                break;
            }

            if (loop_idx > 0 && curve.length > 0) {
                const lastAddedPoint = curve[curve.length - 1];
                // Point before the one just added. If only one point in curve, compare to start_x.
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
        // After loop, check last point if not already marked by break conditions
        if (curve.length > 0) {
            const lastCurvePoint = curve[curve.length - 1];
            if (lastCurvePoint.dxdt_norm !== undefined && lastCurvePoint.dxdt_norm < azeotrope_dxdt_threshold && !lastCurvePoint.isPotentialAzeotrope) {
                lastCurvePoint.isPotentialAzeotrope = true;
            }
        }
        return curve;
    };

    const initial_x_norm = normalizeMoleFractions(initial_x);
    
    // Determine initial temperature for the starting point
    const initialBubbleResult = findBubblePointTemperatureTernary(initial_x_norm, P_system_Pa, componentsData, ternaryWilsonParams, initialTemp_K);
    if (!initialBubbleResult) {
        console.error("Failed to calculate initial bubble point for residue curve start.");
        return null;
    }
    const T_start_K = initialBubbleResult.T_K;

    const curve_forward = await run_direction(initial_x_norm, true, T_start_K);
    const curve_backward = await run_direction(initial_x_norm, false, T_start_K);

    // Combine: backward (reversed) + initial point (from forward[0] or backward[0]) + forward (skipping its first point)
    const combined_curve: ResidueCurveODE = [
        ...curve_backward.reverse().slice(1), // All but last point of backward (which is start_x)
        ...curve_forward
    ];
    
    if (combined_curve.length < 2) {
        console.warn("Residue curve too short, possibly failed or started at an extremum.", initial_x);
        // Return at least the starting point if available
        return curve_forward.length > 0 ? curve_forward : (curve_backward.length > 0 ? curve_backward.reverse() : null);
    }

    return combined_curve;
}
