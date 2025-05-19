import type { CompoundData, AntoineParams } from './vle-types';
// calculatePsat_Pa is needed. NrtlInteractionParams and calculateNrtlGamma (binary) are not directly used by the new ternary logic.
import { calculatePsat_Pa } from './vle-calculations-nrtl';

// --- Constants ---
const R_gas_const_J_molK = 8.31446261815324; // J/mol·K
const MAX_BUBBLE_T_ITERATIONS = 50;
const BUBBLE_T_TOLERANCE_K = 0.01;
const MIN_MOLE_FRACTION = 1e-6; // Adjusted from 1e-9 to stabilize calculations

// --- Interfaces ---
export interface ResidueCurvePointODE {
    x: number[]; // Mole fractions [x0, x1, x2]
    T_K: number; // Temperature at this point
    xi?: number; // Optional integration variable value
}
export type ResidueCurveODE = ResidueCurvePointODE[];

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

        for (let step = 0; step < max_steps_per_direction; step++) {
            const ode_result = residue_curve_ode_dxdxi_nrtl(current_x, P_system_Pa, componentsData, ternaryNrtlParams, last_T_K);
            if (!ode_result) break;

            const { dxdxi, T_bubble_K } = ode_result;
            last_T_K = T_bubble_K;

            curve.push({ x: normalizeMoleFractions(current_x), T_K: T_bubble_K, xi: (forward ? 1 : -1) * step * d_xi_step });

            const next_x_unnormalized = [
                current_x[0] + (forward ? 1 : -1) * dxdxi[0] * d_xi_step,
                current_x[1] + (forward ? 1 : -1) * dxdxi[1] * d_xi_step,
                current_x[2] + (forward ? 1 : -1) * dxdxi[2] * d_xi_step,
            ];
            current_x = normalizeMoleFractions(next_x_unnormalized);

            // Consistent stopping condition:
            // Stop if any component reaches MIN_MOLE_FRACTION or its corresponding upper bound.
            // For a 3-component system, (componentsData.length - 1) is 2.
            const upper_bound_limit = 1.0 - MIN_MOLE_FRACTION * (componentsData.length > 1 ? (componentsData.length - 1) : 0);
            if (current_x.some(xi => xi >= upper_bound_limit ) || current_x.some(xi => xi <= MIN_MOLE_FRACTION)) {
                curve.push({ x: normalizeMoleFractions(current_x), T_K: T_bubble_K, xi: (forward ? 1 : -1) * (step + 1) * d_xi_step });
                break;
            }
            if (step > 0) {
                const prev_x = curve[curve.length - 1].x;
                const diff_sq = prev_x.reduce((sum, px, i) => sum + (px - current_x[i]) ** 2, 0);
                if (Math.sqrt(diff_sq) < d_xi_step * 1e-4) break;
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
