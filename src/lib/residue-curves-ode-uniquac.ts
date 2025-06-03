import type { CompoundData, AntoineParams, UniquacPureComponentParams } from './vle-types';
// UniquacInteractionParams might be binary, we use TernaryUniquacParams for ternary calculations
// import { type UniquacInteractionParams as BinaryUniquacInteractionParams } from './vle-calculations-uniquac'; // Not strictly needed if only using ternary
import { calculatePsat_Pa } from './vle-calculations-unifac'; // Shared Psat utility

// --- Constants ---
const R_gas_const_J_molK = 8.31446261815324; 
const MAX_BUBBLE_T_ITERATIONS = 50;
const BUBBLE_T_TOLERANCE_K = 0.01;
const MIN_MOLE_FRACTION = 1e-9; 
const Z_COORD_NUMBER = 10; 
// const azeotrope_dxdt_threshold = 1e-7; // REMOVED
const NUMERICAL_JACOBIAN_EPS = 1e-7; // For embedded Newton-Raphson

// --- Interfaces ---
export interface ResidueCurvePointODE {
    x: number[]; 
    T_K: number; 
    step?: number; 
    // dxdt_norm?: number; // REMOVED
    // isPotentialAzeotrope?: boolean; // REMOVED
}
export type ResidueCurveODE = ResidueCurvePointODE[];

export interface TernaryUniquacParams {
    a01_J_mol: number; a10_J_mol: number; 
    a02_J_mol: number; a20_J_mol: number; 
    a12_J_mol: number; a21_J_mol: number; 
}

export interface AzeotropeResult { // Consistent result type
    x: number[];
    T_K: number;
    P_Pa: number;
    converged: boolean;
    errorNorm?: number;
    iterations?: number;
    message?: string;
}

// --- START: Embedded Newton-Raphson Utilities ---
function calculateNumericalJacobian(
    func: (vars: number[]) => number[] | null,
    x: number[]
): number[][] | null {
    const n = x.length;
    const jacobian: number[][] = Array(n).fill(null).map(() => Array(n).fill(0));
    const fx = func(x);

    if (!fx) { return null; }

    for (let j = 0; j < n; j++) {
        const x_plus_h = [...x];
        x_plus_h[j] += NUMERICAL_JACOBIAN_EPS;
        const fx_plus_h = func(x_plus_h);

        if (!fx_plus_h) { return null; }

        for (let i = 0; i < n; i++) {
            jacobian[i][j] = (fx_plus_h[i] - fx[i]) / NUMERICAL_JACOBIAN_EPS;
            if (isNaN(jacobian[i][j]) || !isFinite(jacobian[i][j])) {
                if (x[j] > NUMERICAL_JACOBIAN_EPS * 2) {
                    const x_minus_h = [...x];
                    x_minus_h[j] -= NUMERICAL_JACOBIAN_EPS;
                    const fx_minus_h = func(x_minus_h);
                    if (fx_minus_h) {
                        const backward_deriv = (fx[i] - fx_minus_h[i]) / NUMERICAL_JACOBIAN_EPS;
                        if (!isNaN(backward_deriv) && isFinite(backward_deriv)) {
                            jacobian[i][j] = backward_deriv;
                        } else { return null; }
                    } else { return null; }
                } else { return null; }
            }
        }
    }
    return jacobian;
}

function solveLinearSystem(A_orig: number[][], b_orig: number[]): number[] | null {
    const n = A_orig.length;
    const A = A_orig.map(row => [...row]); 
    const b = [...b_orig]; 

    for (let i = 0; i < n; i++) { A[i].push(b[i]); }

    for (let i = 0; i < n; i++) {
        let maxRow = i;
        for (let k = i + 1; k < n; k++) {
            if (Math.abs(A[k][i]) > Math.abs(A[maxRow][i])) maxRow = k;
        }
        [A[i], A[maxRow]] = [A[maxRow], A[i]];
        if (Math.abs(A[i][i]) < 1e-12) return null;
        for (let k = i + 1; k < n + 1; k++) A[i][k] /= A[i][i];
        A[i][i] = 1;
        for (let k = 0; k < n; k++) {
            if (k !== i) {
                const factor = A[k][i];
                for (let j = i; j < n + 1; j++) A[k][j] -= factor * A[i][j];
            }
        }
    }
    const x_sol: number[] = Array(n);
    for (let i = 0; i < n; i++) {
        x_sol[i] = A[i][n];
        if (isNaN(x_sol[i]) || !isFinite(x_sol[i])) return null;
    }
    return x_sol;
}

function solveNewtonRaphsonSystem(
    systemFunction: (vars: number[]) => number[] | null,
    initialVars: number[], maxIterations: number, tolerance: number
): { vars: number[], converged: boolean, iterations?: number, errorNorm?: number, message?: string } {
    let x_k = [...initialVars];
    let fx_k = systemFunction(x_k);
    if (!fx_k) return { vars: x_k, converged: false, message: "NR: Initial function evaluation failed." };

    for (let iter = 0; iter < maxIterations; iter++) {
        const errorNorm = Math.sqrt(fx_k.reduce((sum, val) => sum + val * val, 0));
        if (errorNorm < tolerance) return { vars: x_k, converged: true, iterations: iter, errorNorm, message: "NR: Converged." };
        const J_k = calculateNumericalJacobian(systemFunction, x_k);
        if (!J_k) return { vars: x_k, converged: false, iterations: iter, errorNorm, message: "NR: Jacobian calculation failed." };
        const neg_fx_k = fx_k.map(val => -val);
        const delta_x = solveLinearSystem(J_k, neg_fx_k);
        if (!delta_x) return { vars: x_k, converged: false, iterations: iter, errorNorm, message: "NR: Linear system solution failed." };
        let lambda = 1.0;
        let x_k_plus_1 = x_k.map((val, i) => val + lambda * delta_x[i]);
        let fx_k_plus_1 = systemFunction(x_k_plus_1);
        if (fx_k_plus_1) {
            let newErrorNorm = Math.sqrt(fx_k_plus_1.reduce((sum, val) => sum + val * val, 0));
            let attempts = 0;
            while (newErrorNorm >= errorNorm && attempts < 5 && lambda > 1e-6) { 
                lambda /= 2;
                x_k_plus_1 = x_k.map((val, i) => val + lambda * delta_x[i]);
                const temp_fx = systemFunction(x_k_plus_1);
                if (!temp_fx) { fx_k_plus_1 = null; break; } 
                fx_k_plus_1 = temp_fx;
                newErrorNorm = Math.sqrt(fx_k_plus_1.reduce((sum, val) => sum + val * val, 0));
                attempts++;
            }
        }
        if (!fx_k_plus_1) return { vars: x_k, converged: false, iterations: iter, errorNorm, message: "NR: Func eval failed after update." };
        x_k = x_k_plus_1;
        fx_k = fx_k_plus_1;
    }
    const finalErrorNorm = fx_k ? Math.sqrt(fx_k.reduce((sum, val) => sum + val * val, 0)) : Infinity;
    return { vars: x_k, converged: finalErrorNorm < tolerance, iterations: maxIterations, errorNorm: finalErrorNorm, message: "NR: Max iterations." };
}
// --- END: Embedded Newton-Raphson Utilities ---


function normalizeMoleFractions(x: number[]): number[] { 
    const sum = x.reduce((acc, val) => acc + val, 0);
    if (sum === 0 || isNaN(sum) || !isFinite(sum)) { return x.map(() => 1 / x.length); }
    const upper_bound_val = 1.0 - MIN_MOLE_FRACTION * (x.length > 1 ? (x.length - 1) : 0);
    return x.map(val => Math.max(MIN_MOLE_FRACTION, Math.min(upper_bound_val, val / sum)));
}

function calculateUniquacTaus(
    T_K: number,
    params: TernaryUniquacParams
): { tau01: number, tau10: number, tau02: number, tau20: number, tau12: number, tau21: number } {
    const RT = R_gas_const_J_molK * T_K;
    return {
        tau01: Math.exp(-params.a01_J_mol / RT), 
        tau10: Math.exp(-params.a10_J_mol / RT), 
        tau02: Math.exp(-params.a02_J_mol / RT),
        tau20: Math.exp(-params.a20_J_mol / RT),
        tau12: Math.exp(-params.a12_J_mol / RT),
        tau21: Math.exp(-params.a21_J_mol / RT),
    };
}

export function calculateUniquacGammaTernary(
    x: number[], 
    T_K: number,
    components: CompoundData[], 
    ternaryParams: TernaryUniquacParams
): number[] | null { 
    if (x.length !== 3 || components.length !== 3) return null;
    if (components.some(c => !c.uniquacParams)) return null;

    const r = components.map(c => c.uniquacParams!.r);
    const q = components.map(c => c.uniquacParams!.q);
    const taus = calculateUniquacTaus(T_K, ternaryParams);
    const x0 = x[0], x1 = x[1], x2 = x[2];

    const sum_rx = r[0]*x0 + r[1]*x1 + r[2]*x2;
    if (sum_rx === 0) return null;
    const phi = [r[0]*x0/sum_rx, r[1]*x1/sum_rx, r[2]*x2/sum_rx];

    const sum_qx = q[0]*x0 + q[1]*x1 + q[2]*x2;
    if (sum_qx === 0) return null;
    const theta = [q[0]*x0/sum_qx, q[1]*x1/sum_qx, q[2]*x2/sum_qx];
    
    const l = r.map((ri, i) => (Z_COORD_NUMBER/2)*(ri - q[i]) - (ri - 1));

    const sum_xl = x0*l[0] + x1*l[1] + x2*l[2];
    const lngamma_C = phi.map((phi_i, i) => 
        Math.log(phi_i/x[i]) + (Z_COORD_NUMBER/2)*q[i]*Math.log(theta[i]/phi_i) + l[i] - (phi_i/x[i])*sum_xl
    );

    const sum_theta_tau_k0 = theta[0] + theta[1]*taus.tau10 + theta[2]*taus.tau20;
    const sum_theta_tau_k1 = theta[0]*taus.tau01 + theta[1] + theta[2]*taus.tau21;
    const sum_theta_tau_k2 = theta[0]*taus.tau02 + theta[1]*taus.tau12 + theta[2];

    if (sum_theta_tau_k0 === 0 || sum_theta_tau_k1 === 0 || sum_theta_tau_k2 === 0) return null;

    const lngamma_R = [
        q[0]*(1 - Math.log(sum_theta_tau_k0) - (theta[0]/sum_theta_tau_k0 + theta[1]*taus.tau01/sum_theta_tau_k1 + theta[2]*taus.tau02/sum_theta_tau_k2) ),
        q[1]*(1 - Math.log(sum_theta_tau_k1) - (theta[0]*taus.tau10/sum_theta_tau_k0 + theta[1]/sum_theta_tau_k1 + theta[2]*taus.tau12/sum_theta_tau_k2) ),
        q[2]*(1 - Math.log(sum_theta_tau_k2) - (theta[0]*taus.tau20/sum_theta_tau_k0 + theta[1]*taus.tau21/sum_theta_tau_k1 + theta[2]/sum_theta_tau_k2) )
    ];
    
    return lngamma_C.map((lngC, i) => Math.exp(lngC + lngamma_R[i]));
}

export function findBubblePointTemperatureTernaryUniquac(
    x_liq: number[],
    P_system_Pa: number,
    components: CompoundData[], 
    ternaryUniquacParams: TernaryUniquacParams,
    initialTempGuess_K: number
): { T_K: number, gammas: number[], Psats: number[] } | null {
    let T_K = initialTempGuess_K;
    if (!x_liq || x_liq.length !== components.length || x_liq.some(isNaN)) {
        return null;
    }
    const x = normalizeMoleFractions(x_liq); 

    if (components.some(c => !c.antoine || !c.uniquacParams)) return null;

    for (let iter = 0; iter < MAX_BUBBLE_T_ITERATIONS; iter++) {
        const Psats = components.map(c => calculatePsat_Pa(c.antoine!, T_K));
        if (Psats.some(isNaN)) return null;

        const gammas = calculateUniquacGammaTernary(x, T_K, components, ternaryUniquacParams);
        if (!gammas || gammas.some(isNaN)) return null;

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

export function residue_curve_ode_dxdxi_uniquac(
    x_liq_current: number[],
    P_system_Pa: number,
    components: CompoundData[],
    ternaryUniquacParams: TernaryUniquacParams,
    initialTempGuess_K: number
): { dxdxi: number[], T_bubble_K: number } | null {
    if (!x_liq_current || x_liq_current.length !== components.length || x_liq_current.some(isNaN)) {
        return null;
    }
    const x_norm = normalizeMoleFractions(x_liq_current);

    const bubbleResult = findBubblePointTemperatureTernaryUniquac(x_norm, P_system_Pa, components, ternaryUniquacParams, initialTempGuess_K);
    if (!bubbleResult) return null;

    const { T_K: T_bubble_K, gammas, Psats } = bubbleResult;

    const y_vap = calculateVaporCompositionTernary(x_norm, P_system_Pa, gammas, Psats);
    if (!y_vap) return null;

    const dxdxi = [ x_norm[0] - y_vap[0], x_norm[1] - y_vap[1], x_norm[2] - y_vap[2] ];
    return { dxdxi, T_bubble_K };
}


// --- Direct Azeotrope Solver for UNIQUAC ---
export function findAzeotropeUniquac(
    P_system_Pa: number,
    components: CompoundData[],
    ternaryUniquacParams: TernaryUniquacParams,
    initialGuess_x_independent: number[], // Initial guess for [x0, x1]
    initialGuess_T_K: number,
    maxIterations: number = 100,
    tolerance: number = 1e-8
): AzeotropeResult | null {
    if (components.length !== 3) { return null; }
    if (components.some(c => !c.antoine || !c.uniquacParams)) { return null; }

    const systemFunction = (vars: number[]): number[] | null => {
        let x0_g = vars[0];
        let x1_g = vars[1];
        const T_K_g = vars[2];

        x0_g = Math.max(MIN_MOLE_FRACTION, Math.min(1.0 - 2 * MIN_MOLE_FRACTION, x0_g));
        x1_g = Math.max(MIN_MOLE_FRACTION, Math.min(1.0 - x0_g - MIN_MOLE_FRACTION, x1_g));
        const x_guess_arr = [x0_g, x1_g, 1.0 - x0_g - x1_g];
        const x_norm = normalizeMoleFractions(x_guess_arr); 

        const Psats = components.map(c => calculatePsat_Pa(c.antoine!, T_K_g));
        if (Psats.some(isNaN)) return null;

        const gammas = calculateUniquacGammaTernary(x_norm, T_K_g, components, ternaryUniquacParams);
        if (!gammas || gammas.some(isNaN)) return null;

        const y_vap_unnorm = components.map((c, i) => x_norm[i] * gammas[i] * Psats[i]);
        const P_calc = y_vap_unnorm.reduce((s, val) => s + val, 0);
        if (P_calc === 0 || isNaN(P_calc)) return null;

        const y_norm = y_vap_unnorm.map(yiP => yiP / P_calc);

        const F = [
            y_norm[0] - x_norm[0],
            y_norm[1] - x_norm[1],
            P_calc - P_system_Pa
        ];
        return F;
    };

    const initialVars = [
        Math.max(MIN_MOLE_FRACTION, Math.min(1.0 - MIN_MOLE_FRACTION, initialGuess_x_independent[0])),
        Math.max(MIN_MOLE_FRACTION, Math.min(1.0 - initialGuess_x_independent[0] - MIN_MOLE_FRACTION, initialGuess_x_independent[1])),
        Math.max(150, initialGuess_T_K)
    ];

    try {
        const solution = solveNewtonRaphsonSystem(systemFunction, initialVars, maxIterations, tolerance);
        if (solution && solution.converged) {
            const x0_s = solution.vars[0];
            const x1_s = solution.vars[1];
            const T_K_s = solution.vars[2];
            const x_s_final = normalizeMoleFractions([x0_s, x1_s, 1.0 - x0_s - x1_s]);
            const errorsAtSol = systemFunction(solution.vars);
            const errorNormAtSol = errorsAtSol ? Math.sqrt(errorsAtSol.reduce((s, e) => s + e * e, 0)) : Infinity;
            return {
                x: x_s_final, T_K: T_K_s, P_Pa: P_system_Pa, converged: true,
                errorNorm: errorNormAtSol, iterations: solution.iterations, message: "Converged (UNIQUAC)."
            };
        } else {
            return {
                x: [initialVars[0], initialVars[1], 1.0 - initialVars[0] - initialVars[1]], T_K: initialVars[2],
                P_Pa: P_system_Pa, converged: false, message: solution?.message || "NR (UNIQUAC) did not converge.",
                iterations: solution?.iterations
            };
        }
    } catch (e: any) {
        return {
            x: [initialVars[0], initialVars[1], 1.0 - initialVars[0] - initialVars[1]], T_K: initialVars[2],
            P_Pa: P_system_Pa, converged: false, message: `Solver exception (UNIQUAC): ${e.message}`
        };
    }
}


export async function simulateSingleResidueCurveODE_Uniquac(
    initial_x: number[], P_system_Pa: number, componentsData: CompoundData[],
    ternaryUniquacParams: TernaryUniquacParams, d_xi_step: number,
    max_steps_per_direction: number, initialTemp_K: number
): Promise<ResidueCurveODE | null> {

    const run_direction = async (start_x: number[], forward: boolean, current_T_guess: number): Promise<ResidueCurveODE> => {
        const curve: ResidueCurveODE = [];
        let current_x = [...start_x];
        let last_T_K = current_T_guess;

        for (let loop_step = 0; loop_step < max_steps_per_direction; loop_step++) {
            const ode_result = residue_curve_ode_dxdxi_uniquac(current_x, P_system_Pa, componentsData, ternaryUniquacParams, last_T_K);
            if (!ode_result) break;

            const { dxdxi, T_bubble_K } = ode_result;
            last_T_K = T_bubble_K;
            // const dxdt_norm = Math.sqrt(dxdxi[0]**2 + dxdxi[1]**2 + dxdxi[2]**2); // REMOVED

            const currentPointData: ResidueCurvePointODE = {
                x: normalizeMoleFractions(current_x),
                T_K: T_bubble_K,
                step: (forward ? 1 : -1) * loop_step * d_xi_step,
                // dxdt_norm, // REMOVED
                // isPotentialAzeotrope: false // REMOVED
            };
            
            if (curve.length > 0) {
                const lastPoint = curve[curve.length - 1];
                const x_diff_sq = lastPoint.x.reduce((sum, val, i) => sum + (val - currentPointData.x[i])**2, 0);
                if (Math.sqrt(x_diff_sq) < 1e-7 && Math.abs(lastPoint.T_K - currentPointData.T_K) < 1e-4) {
                    // Stagnation detected
                    // if (dxdt_norm < azeotrope_dxdt_threshold) { // REMOVED
                        // lastPoint.isPotentialAzeotrope = true; // REMOVED
                        // lastPoint.dxdt_norm = Math.min(lastPoint.dxdt_norm ?? Infinity, dxdt_norm); // REMOVED
                        break; 
                    // } // REMOVED
                } else { curve.push(currentPointData); }
            } else { curve.push(currentPointData); }

            // const pointToCheck = curve[curve.length - 1]; // REMOVED BLOCK
            // if (pointToCheck.dxdt_norm !== undefined && pointToCheck.dxdt_norm < azeotrope_dxdt_threshold) { // REMOVED
                // pointToCheck.isPotentialAzeotrope = true; // REMOVED
                // break;  // REMOVED
            // } // REMOVED
            
            const next_x_unnormalized: number[] = [
                current_x[0] + (forward ? 1 : -1) * dxdxi[0] * d_xi_step,
                current_x[1] + (forward ? 1 : -1) * dxdxi[1] * d_xi_step,
                current_x[2] + (forward ? 1 : -1) * dxdxi[2] * d_xi_step,
            ];
            current_x = normalizeMoleFractions(next_x_unnormalized);

            const upper_bound_limit = 1.0 - MIN_MOLE_FRACTION * (componentsData.length > 1 ? (componentsData.length - 1) : 0);
            if (current_x.some(val => val >= upper_bound_limit) || current_x.some(val => val <= MIN_MOLE_FRACTION)) {
                const final_x_normalized = normalizeMoleFractions(current_x);
                let lastPushedPoint = curve.length > 0 ? curve[curve.length - 1] : null;
                if (lastPushedPoint && final_x_normalized.every((val, i) => Math.abs(val - lastPushedPoint!.x[i]) < 1e-7)) {
                    // lastPushedPoint.isPotentialAzeotrope = true; // REMOVED
                    // lastPushedPoint.dxdt_norm = 0; // REMOVED
                } else {
                    curve.push({ 
                        x: final_x_normalized, 
                        T_K: last_T_K, 
                        step: (forward ? 1 : -1) * (loop_step + 1) * d_xi_step,
                        // dxdt_norm: 0, // REMOVED
                        // isPotentialAzeotrope: true // REMOVED
                    });
                }
                break;
            }
            if (loop_step > 0 && curve.length > 0) {
                const lastAddedPoint = curve[curve.length - 1];
                const x_before_euler_step = curve.length > 1 ? curve[curve.length - 2].x : start_x;
                
                const diff_sq_from_prev_step = lastAddedPoint.x.reduce((sum, px, i) => sum + (px - x_before_euler_step[i]) ** 2, 0);

                if (Math.sqrt(diff_sq_from_prev_step) < d_xi_step * 1e-5) { 
                    // if (lastAddedPoint.dxdt_norm !== undefined && lastAddedPoint.dxdt_norm < azeotrope_dxdt_threshold) { // REMOVED
                        // lastAddedPoint.isPotentialAzeotrope = true; // REMOVED
                    // } // REMOVED
                    break; 
                }
            }
        }
        // if (curve.length > 0) { // REMOVED BLOCK
            // const lastCurvePoint = curve[curve.length - 1]; // REMOVED
            // if (lastCurvePoint.dxdt_norm !== undefined && lastCurvePoint.dxdt_norm < azeotrope_dxdt_threshold && !lastCurvePoint.isPotentialAzeotrope) { // REMOVED
                // lastCurvePoint.isPotentialAzeotrope = true; // REMOVED
            // } // REMOVED
        // } // REMOVED
        return curve;
    };

    if (!initial_x || initial_x.length !== componentsData.length || initial_x.some(isNaN)) {
        return null;
    }
    const initial_x_norm = normalizeMoleFractions(initial_x);
    const initialBubbleResult = findBubblePointTemperatureTernaryUniquac(initial_x_norm, P_system_Pa, componentsData, ternaryUniquacParams, initialTemp_K);
    if (!initialBubbleResult) return null;
    const T_start_K = initialBubbleResult.T_K;

    const curve_forward = await run_direction(initial_x_norm, true, T_start_K);
    const curve_backward = await run_direction(initial_x_norm, false, T_start_K);

    const combined_curve: ResidueCurveODE = [ ...curve_backward.reverse().slice(1), ...curve_forward ];
    // if (combined_curve.length < 2 && curve_forward.length === 1 && curve_forward[0].isPotentialAzeotrope) return curve_forward; // REMOVED
    if (combined_curve.length < 2) return curve_forward.length > 0 ? curve_forward : (curve_backward.length > 0 ? curve_backward.reverse() : null);
    return combined_curve;
}

// --- END OF FILE residue-curves-ode-uniquac.ts ---
