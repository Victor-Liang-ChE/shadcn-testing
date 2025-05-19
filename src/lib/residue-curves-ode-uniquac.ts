import type { CompoundData, AntoineParams, UniquacPureComponentParams } from './vle-types';
import { type UniquacInteractionParams as BinaryUniquacInteractionParams } from './vle-calculations-uniquac';
import { calculatePsat_Pa } from './vle-calculations-unifac';

// --- Constants ---
const R_gas_const_J_molK = 8.31446261815324; // J/mol·K - Defined locally
const MAX_BUBBLE_T_ITERATIONS = 50;
const BUBBLE_T_TOLERANCE_K = 0.01;
const MIN_MOLE_FRACTION = 1e-9;
const Z_COORD_NUMBER = 10; // UNIQUAC coordination number

// --- Interfaces ---
export interface ResidueCurvePointODE {
    x: number[]; // Mole fractions [x0, x1, x2]
    T_K: number; // Temperature at this point
    step?: number; // Optional integration variable value - CHANGED FROM xi
}
export type ResidueCurveODE = ResidueCurvePointODE[];

export interface TernaryUniquacParams {
    // Interaction energy parameters u_ij - u_jj (or similar, depends on DB source for a_ij)
    // Units: J/mol
    a01_J_mol: number; a10_J_mol: number; // For pair 0-1
    a02_J_mol: number; a20_J_mol: number; // For pair 0-2
    a12_J_mol: number; a21_J_mol: number; // For pair 1-2
}

function normalizeMoleFractions(x: number[]): number[] {
    const sum = x.reduce((acc, val) => acc + val, 0);
    if (sum === 0) return x.map(() => 1 / x.length);
    return x.map(val => Math.max(MIN_MOLE_FRACTION, Math.min(1.0 - MIN_MOLE_FRACTION * (x.length -1), val / sum)));
}

function calculateUniquacTaus(
    T_K: number,
    params: TernaryUniquacParams
): { tau01: number, tau10: number, tau02: number, tau20: number, tau12: number, tau21: number } {
    const RT = R_gas_const_J_molK * T_K;
    return {
        tau01: Math.exp(-params.a01_J_mol / RT), // exp(-(u_01 - u_11)/RT) or exp(-a_01/RT)
        tau10: Math.exp(-params.a10_J_mol / RT), // exp(-(u_10 - u_00)/RT) or exp(-a_10/RT)
        tau02: Math.exp(-params.a02_J_mol / RT),
        tau20: Math.exp(-params.a20_J_mol / RT),
        tau12: Math.exp(-params.a12_J_mol / RT),
        tau21: Math.exp(-params.a21_J_mol / RT),
    };
}

export function calculateUniquacGammaTernary(
    x: number[], // mole fractions [x0, x1, x2]
    T_K: number,
    components: CompoundData[], // Must have uniquacParams (r, q)
    ternaryParams: TernaryUniquacParams
): number[] | null { // Returns [gamma0, gamma1, gamma2] or null on error
    if (x.length !== 3 || components.length !== 3) return null;
    if (components.some(c => !c.uniquacParams)) return null;

    const r = components.map(c => c.uniquacParams!.r);
    const q = components.map(c => c.uniquacParams!.q);
    const taus = calculateUniquacTaus(T_K, ternaryParams);
    const x0 = x[0], x1 = x[1], x2 = x[2];

    // Segment fractions phi_i = r_i * x_i / sum(r_j * x_j)
    const sum_rx = r[0]*x0 + r[1]*x1 + r[2]*x2;
    if (sum_rx === 0) return null;
    const phi = [r[0]*x0/sum_rx, r[1]*x1/sum_rx, r[2]*x2/sum_rx];

    // Surface area fractions theta_i = q_i * x_i / sum(q_j * x_j)
    const sum_qx = q[0]*x0 + q[1]*x1 + q[2]*x2;
    if (sum_qx === 0) return null;
    const theta = [q[0]*x0/sum_qx, q[1]*x1/sum_qx, q[2]*x2/sum_qx];
    
    // l_i = (Z/2)*(r_i - q_i) - (r_i - 1)
    const l = r.map((ri, i) => (Z_COORD_NUMBER/2)*(ri - q[i]) - (ri - 1));

    // Combinatorial part: ln(gamma_i^C) = ln(phi_i/x_i) + (Z/2)*q_i*ln(theta_i/phi_i) + l_i - (phi_i/x_i)*sum(x_j*l_j)
    const sum_xl = x0*l[0] + x1*l[1] + x2*l[2];
    const lngamma_C = phi.map((phi_i, i) => 
        Math.log(phi_i/x[i]) + (Z_COORD_NUMBER/2)*q[i]*Math.log(theta[i]/phi_i) + l[i] - (phi_i/x[i])*sum_xl
    );

    // Residual part: ln(gamma_i^R) = q_i * (1 - ln(sum_j theta_j * tau_ji) - sum_j (theta_j * tau_ij / sum_k theta_k * tau_kj) )
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
    components: CompoundData[], // Antoine and Uniquac r,q params needed
    ternaryUniquacParams: TernaryUniquacParams,
    initialTempGuess_K: number
): { T_K: number, gammas: number[], Psats: number[] } | null {
    let T_K = initialTempGuess_K;
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

// Export this function for azeotrope detection
export function residue_curve_ode_dxdxi_uniquac(
    x_liq_current: number[],
    P_system_Pa: number,
    components: CompoundData[],
    ternaryUniquacParams: TernaryUniquacParams,
    initialTempGuess_K: number
): { dxdxi: number[], T_bubble_K: number } | null {
    const x_norm = normalizeMoleFractions(x_liq_current);

    const bubbleResult = findBubblePointTemperatureTernaryUniquac(x_norm, P_system_Pa, components, ternaryUniquacParams, initialTempGuess_K);
    if (!bubbleResult) return null;

    const { T_K: T_bubble_K, gammas, Psats } = bubbleResult;

    const y_vap = calculateVaporCompositionTernary(x_norm, P_system_Pa, gammas, Psats);
    if (!y_vap) return null;

    const dxdxi = [ x_norm[0] - y_vap[0], x_norm[1] - y_vap[1], x_norm[2] - y_vap[2] ];
    return { dxdxi, T_bubble_K };
}

export async function simulateSingleResidueCurveODE_Uniquac(
    initial_x: number[],
    P_system_Pa: number,
    componentsData: CompoundData[],
    ternaryUniquacParams: TernaryUniquacParams,
    d_xi_step: number,
    max_steps_per_direction: number,
    initialTemp_K: number
): Promise<ResidueCurveODE | null> {
    
    const run_direction = async (start_x: number[], forward: boolean, current_T_guess: number): Promise<ResidueCurveODE> => {
        const curve: ResidueCurveODE = [];
        let current_x = [...start_x];
        let last_T_K = current_T_guess;

        for (let loop_step = 0; loop_step < max_steps_per_direction; loop_step++) { // Renamed loop variable for clarity
            const ode_result = residue_curve_ode_dxdxi_uniquac(current_x, P_system_Pa, componentsData, ternaryUniquacParams, last_T_K);
            if (!ode_result) break;

            const { dxdxi, T_bubble_K } = ode_result;
            last_T_K = T_bubble_K;

            curve.push({ x: normalizeMoleFractions(current_x), T_K: T_bubble_K, step: (forward ? 1 : -1) * loop_step * d_xi_step }); // CHANGED xi to step, used loop_step
            
            const next_x_unnormalized = [
                current_x[0] + (forward ? 1 : -1) * dxdxi[0] * d_xi_step,
                current_x[1] + (forward ? 1 : -1) * dxdxi[1] * d_xi_step,
                current_x[2] + (forward ? 1 : -1) * dxdxi[2] * d_xi_step,
            ];
            current_x = normalizeMoleFractions(next_x_unnormalized);

            if (current_x.some(val => val >= 1.0 - MIN_MOLE_FRACTION * 2) || current_x.some(val => val <= MIN_MOLE_FRACTION * 2)) { // Used val instead of xi
                 curve.push({ x: normalizeMoleFractions(current_x), T_K: T_bubble_K, step: (forward ? 1 : -1) * (loop_step +1) * d_xi_step }); // CHANGED xi to step, used loop_step
                break;
            }
            if (loop_step > 0) { // Used loop_step
                const prev_x = curve[curve.length-1].x;
                const diff_sq = prev_x.reduce((sum, px, i) => sum + (px - current_x[i])**2, 0);
                if (Math.sqrt(diff_sq) < d_xi_step * 1e-3) break;
            }
        }
        return curve;
    };

    const initial_x_norm = normalizeMoleFractions(initial_x);
    const initialBubbleResult = findBubblePointTemperatureTernaryUniquac(initial_x_norm, P_system_Pa, componentsData, ternaryUniquacParams, initialTemp_K);
    if (!initialBubbleResult) return null;
    
    const T_start_K = initialBubbleResult.T_K;

    const curve_forward = await run_direction(initial_x_norm, true, T_start_K);
    const curve_backward = await run_direction(initial_x_norm, false, T_start_K);

    const combined_curve: ResidueCurveODE = [ ...curve_backward.reverse().slice(1), ...curve_forward ];
    
    return combined_curve.length < 2 ? (curve_forward.length > 0 ? curve_forward : (curve_backward.length > 0 ? curve_backward.reverse() : null)) : combined_curve;
}
