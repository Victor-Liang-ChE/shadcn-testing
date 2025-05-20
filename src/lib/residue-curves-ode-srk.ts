import type { CompoundData, SrkPureComponentParams, AntoineParams } from '@/lib/vle-types';
import { R_gas_const_J_molK, solveCubicEOS, type SrkInteractionParams as BinarySrkInteractionParams } from './vle-calculations-srk';
import { calculatePsat_Pa } from './vle-calculations-unifac';

// --- Constants ---
const MAX_BUBBLE_T_ITERATIONS = 50;
const MAX_K_VALUE_ITERATIONS = 30;
const BUBBLE_T_TOLERANCE_K = 0.001;
const K_VALUE_TOLERANCE = 1e-6;
const MIN_MOLE_FRACTION = 1e-9;
const BUBBLE_T_CONVERGENCE_TOLERANCE = 1e-5;
const NUMERICAL_JACOBIAN_EPS = 1e-7; // For embedded Newton-Raphson

// --- Interfaces ---
export interface ResidueCurvePointODE {
    x: number[]; 
    T_K: number; 
    step?: number; 
}
export type ResidueCurveODE = ResidueCurvePointODE[];

export interface TernarySrkParams {
    k01: number; k10?: number;
    k02: number; k20?: number;
    k12: number; k21?: number;
}

export interface AzeotropeResult { // Consistent result type
    T_K: number;
    x: number[];
    y: number[];
    K_values: number[];
    ZL: number;
    ZV: number;
}

// --- Helper: Ensure this is the ONLY definition in this file ---
function normalizeMoleFractions(x: number[]): number[] {
    const sum = x.reduce((acc, val) => acc + val, 0);
    if (sum === 0) return x.map(() => 1 / x.length);
    const upper_bound_val = 1.0 - MIN_MOLE_FRACTION * (x.length > 1 ? (x.length - 1) : 0);
    return x.map(val => Math.max(MIN_MOLE_FRACTION, Math.min(upper_bound_val, val / sum)));
}

// --- START: Embedded Newton-Raphson Utilities ---
function calculateNumericalJacobian(
    f: (vars: number[]) => number[],
    vars: number[],
    eps: number = NUMERICAL_JACOBIAN_EPS
): number[][] | null {
    const n = vars.length;
    const jacobian: number[][] = Array.from({ length: n }, () => Array(n).fill(0));
    const f0 = f(vars);

    for (let j = 0; j < n; j++) {
        const vars_j_plus_eps = [...vars];
        vars_j_plus_eps[j] += eps;
        const f_j_plus_eps = f(vars_j_plus_eps);

        for (let i = 0; i < n; i++) {
            jacobian[i][j] = (f_j_plus_eps[i] - f0[i]) / eps;
        }
    }
    return jacobian;
}

function solveLinearSystem(A_orig: number[][], b_orig: number[]): number[] | null {
    const n = b_orig.length;
    const A = A_orig.map((row, i) => [...row, b_orig[i]]);
    for (let i = 0; i < n; i++) {
        // Pivoting
        let maxRow = i;
        for (let k = i + 1; k < n; k++) {
            if (Math.abs(A[k][i]) > Math.abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        [A[i], A[maxRow]] = [A[maxRow], A[i]];

        // Elimination
        for (let k = i + 1; k < n; k++) {
            const factor = A[k][i] / A[i][i];
            for (let j = i; j < n + 1; j++) {
                A[k][j] -= A[i][j] * factor;
            }
        }
    }

    // Back substitution
    const x = Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
        x[i] = A[i][n] / A[i][i];
        for (let k = i - 1; k >= 0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;
}

function solveNewtonRaphsonSystem(
    f: (vars: number[]) => number[] | null, // Changed to accept null
    jacobian: (vars: number[]) => number[][], // Analytical Jacobian function
    vars0: number[],
    options?: { maxIterations?: number, tolerance?: number }
): { vars: number[], converged: boolean, iterations?: number, errorNorm?: number, message?: string } {
    const maxIterations = options?.maxIterations || 100;
    const tolerance = options?.tolerance || 1e-6;
    let currentVars = [...vars0];
    let converged = false;
    let iterCount = 0;
    let lastErrorNorm: number | undefined;

    for (iterCount = 0; iterCount < maxIterations; iterCount++) {
        const F_current = f(currentVars);
        if (!F_current) {
            return { vars: currentVars, converged: false, iterations: iterCount, errorNorm: Infinity, message: "System function returned null during iteration." };
        }

        lastErrorNorm = Math.sqrt(F_current.reduce((sum, fi) => sum + fi ** 2, 0));

        if (lastErrorNorm < tolerance) {
            converged = true;
            break;
        }

        const J_current = jacobian(currentVars); // Correctly call the provided jacobian function
        // Assuming jacobian function is robust and doesn't return null/invalid for expected inputs.
        // If J_current could be null, a check would be needed here.

        const neg_F_current = F_current.map(v => -v);
        const delta = solveLinearSystem(J_current, neg_F_current);

        if (!delta) {
            return { vars: currentVars, converged: false, iterations: iterCount, errorNorm: lastErrorNorm, message: "Linear system solution failed." };
        }

        currentVars = currentVars.map((v, i) => v - delta[i]);
    }

    // If loop finished due to max iterations, calculate final error norm for the last state of currentVars
    if (!converged && iterCount === maxIterations) {
        const F_final = f(currentVars);
        if (F_final) {
            lastErrorNorm = Math.sqrt(F_final.reduce((sum, fi) => sum + fi ** 2, 0));
            if (lastErrorNorm < tolerance) { // One last check
                converged = true;
            }
        } else {
            lastErrorNorm = Infinity; // f failed at the very end
        }
    }

    return { 
        vars: currentVars, 
        converged, 
        iterations: iterCount, 
        errorNorm: lastErrorNorm, 
        message: converged ? "Converged" : (iterCount === maxIterations ? "Max iterations reached" : "Did not converge") 
    };
}
// --- END: Embedded Newton-Raphson Utilities ---


// --- Helper functions for SRK ---
function getKijSrk(idx1: number, idx2: number, params: TernarySrkParams): number {
    const key = `${Math.min(idx1, idx2)}${Math.max(idx1, idx2)}`;
    if (idx1 === idx2) return 0;
    if (key === '01') return params.k01;
    if (key === '02') return params.k02;
    if (key === '12') return params.k12;
    if (idx1 === 1 && idx2 === 0 && params.k10 !== undefined) return params.k10;
    if (idx1 === 2 && idx2 === 0 && params.k20 !== undefined) return params.k20;
    if (idx1 === 2 && idx2 === 1 && params.k21 !== undefined) return params.k21;
    return 0;
}

function calculateSrkPureComponentParams(
    T_K: number,
    srkParams: SrkPureComponentParams
): { ac_i: number, b_i: number, m_i: number, alpha_i: number, a_i: number } {
    const Tr_i = T_K / srkParams.Tc_K;
    const m_i = 0.480 + 1.574 * srkParams.omega - 0.176 * srkParams.omega * srkParams.omega;
    const alpha_i = Math.pow(1 + m_i * (1 - Math.sqrt(Tr_i)), 2);
    const ac_i = 0.42748 * R_gas_const_J_molK * R_gas_const_J_molK * srkParams.Tc_K * srkParams.Tc_K / srkParams.Pc_Pa;
    const a_i = ac_i * alpha_i;
    const b_i = 0.08664 * R_gas_const_J_molK * srkParams.Tc_K / srkParams.Pc_Pa;
    return { ac_i, b_i, m_i, alpha_i, a_i };
}

function calculateSrkMixtureParamsTernary(
    phaseCompositions: number[],
    T_K: number,
    componentsData: CompoundData[],
    ternaryKij: TernarySrkParams
): { mixture_a: number, mixture_b: number, component_a_values: number[], component_b_values: number[] } {
    const numComponents = componentsData.length;
    const component_a_values = componentsData.map(c => calculateSrkPureComponentParams(T_K, c.srkParams!).a_i);
    const component_b_values = componentsData.map(c => calculateSrkPureComponentParams(T_K, c.srkParams!).b_i);

    let mixture_a = 0;
    for (let i = 0; i < numComponents; i++) {
        for (let j = 0; j < numComponents; j++) {
            const kij = getKijSrk(i, j, ternaryKij);
            mixture_a += phaseCompositions[i] * phaseCompositions[j] * Math.sqrt(component_a_values[i] * component_a_values[j]) * (1 - kij);
        }
    }
    const mixture_b = phaseCompositions.reduce((sum, x_i, i) => sum + x_i * component_b_values[i], 0);
    return { mixture_a, mixture_b, component_a_values, component_b_values };
}

function calculateSrkCoefficients(mixture_a: number, mixture_b: number, T_K: number, P_Pa: number): { p: number, q: number, r: number, A_mix_eos: number, B_mix_eos: number } {
    const A_mix_eos = mixture_a * P_Pa / Math.pow(R_gas_const_J_molK * T_K, 2);
    const B_mix_eos = mixture_b * P_Pa / (R_gas_const_J_molK * T_K);
    const p_coeff = -1.0;
    const q_coeff = A_mix_eos - B_mix_eos - B_mix_eos * B_mix_eos;
    const r_coeff = -A_mix_eos * B_mix_eos;
    return { p: p_coeff, q: q_coeff, r: r_coeff, A_mix_eos, B_mix_eos };
}

function calculateSrkFugacityCoefficientsTernary(
    phaseCompositions: number[],
    T_K: number, P_Pa: number, Z_phase: number,
    mixture_a: number, mixture_b: number,
    component_a_values: number[], component_b_values: number[],
    ternaryKij: TernarySrkParams
): number[] {
    const numComponents = phaseCompositions.length;
    if (numComponents !== 3) return [NaN, NaN, NaN];

    const A_mix_eos = mixture_a * P_Pa / Math.pow(R_gas_const_J_molK * T_K, 2);
    const B_mix_eos = mixture_b * P_Pa / (R_gas_const_J_molK * T_K);
    const fugacity_coeffs: number[] = [];

    for (let k = 0; k < numComponents; k++) {
        let sum_term_a_ki = 0;
        for (let i = 0; i < numComponents; i++) {
            const kij_ki = getKijSrk(k, i, ternaryKij);
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
        if (Math.abs(B_mix_eos) > 1e-12) {
            term3_factor = A_mix_eos / B_mix_eos;
        } else if (Math.abs(A_mix_eos) < 1e-12) {
            term3_factor = 0;
        } else {
            fugacity_coeffs.push(NaN);
            continue;
        }
        
        const term3_main_num = 2 * sum_term_a_ki;
        let term3_main_frac1 = 0;
        if (Math.abs(mixture_a) > 1e-12) term3_main_frac1 = term3_main_num / mixture_a;

        let term3_main_frac2 = 0;
        if (Math.abs(mixture_b) > 1e-12) term3_main_frac2 = component_b_values[k] / mixture_b;

        const term3_main = term3_main_frac1 - term3_main_frac2;
        
        let term3_log = 0;
        if (Z_phase > MIN_MOLE_FRACTION && (1 + B_mix_eos / Z_phase) > MIN_MOLE_FRACTION) {
            term3_log = Math.log(1 + B_mix_eos / Z_phase);
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


// --- Bubble Point and ODE functions (findBubblePointTemperatureTernarySrk, calculateVaporCompositionEos, residue_curve_ode_dxdxi_srk) ---
export function findBubblePointTemperatureTernarySrk(
    x_liq: number[], P_system_Pa: number, components: CompoundData[],
    ternarySrkKijParams: TernarySrkParams, initialTempGuess_K: number
): { T_K: number, y_vap: number[], K_values: number[], ZL: number, ZV: number } | null {
    // Ensure x_liq is validated before calling normalizeMoleFractions
    if (!x_liq || x_liq.length !== components.length || x_liq.some(isNaN)) return null;
    const x = normalizeMoleFractions(x_liq);
    let T_K_current = initialTempGuess_K; // Use a distinct name for the loop variable

    if (components.some(c => !c.antoine || !c.srkParams) || components.length !== 3) {
        // console.error("SRK Bubble T: Invalid component data or not ternary system.");
        return null;
    }

    let initial_psats_all_valid = true;
    let K_values = components.map(c => {
        const psat = calculatePsat_Pa(c.antoine!, T_K_current); // Use T_K_current
        if (isNaN(psat)) {
            // console.warn(`SRK Bubble T: Psat is NaN for component ${c.name || c.cas_number} at T=${T_K_current.toFixed(1)}K. Antoine: A=${c.antoine?.A}, B=${c.antoine?.B}, C=${c.antoine?.C}`);
            initial_psats_all_valid = false;
        }
        return isNaN(psat) || P_system_Pa === 0 ? 1.0 : psat / P_system_Pa;
    });

    if (!initial_psats_all_valid) {
        // console.warn(`SRK Bubble T: One or more initial Psat values were NaN at T=${T_K_current.toFixed(1)}K. K-values were set to 1.0 for these components. Initial K-values: [${K_values.map(k=>k.toFixed(3)).join(', ')}]`);
        if (K_values.every(k => Math.abs(k - 1.0) < 1e-9)) {
            // console.warn(`SRK Bubble T: All initial K-values defaulted to 1.0 due to Psat failures. This may lead to trivial convergence.`);
        }
    }

    let ZL_final: number | null = null;
    let ZV_final: number | null = null;
    let sum_Kx_for_T_loop = 0; // Initialize here
    // sum_Kx_for_T_loop_final was removed, using sum_Kx_for_T_loop directly.

    for (let T_iter = 0; T_iter < MAX_BUBBLE_T_ITERATIONS; T_iter++) {
        // sum_Kx_for_T_loop is updated inside K_iter loop

        for (let K_iter = 0; K_iter < MAX_K_VALUE_ITERATIONS; K_iter++) {
            const y_vap = normalizeMoleFractions(K_values.map((Ki, i) => Ki * x[i]));

            const liqMix = calculateSrkMixtureParamsTernary(x, T_K_current, components, ternarySrkKijParams);
            const vapMix = calculateSrkMixtureParamsTernary(y_vap, T_K_current, components, ternarySrkKijParams);

            const liqEosCoeffs = calculateSrkCoefficients(liqMix.mixture_a, liqMix.mixture_b, T_K_current, P_system_Pa);
            const vapEosCoeffs = calculateSrkCoefficients(vapMix.mixture_a, vapMix.mixture_b, T_K_current, P_system_Pa);

            const ZL_all_roots = solveCubicEOS(1, liqEosCoeffs.p, liqEosCoeffs.q, liqEosCoeffs.r);
            const ZV_all_roots = solveCubicEOS(1, vapEosCoeffs.p, vapEosCoeffs.q, vapEosCoeffs.r);
            
            ZL_final = null;
            if (ZL_all_roots && ZL_all_roots.length > 0) {
                const valid_ZL = ZL_all_roots.filter(z => z > liqEosCoeffs.B_mix_eos + MIN_MOLE_FRACTION);
                if (valid_ZL.length > 0) ZL_final = valid_ZL[0]; // Smallest positive root for liquid
            }

            ZV_final = null;
            if (ZV_all_roots && ZV_all_roots.length > 0) {
                const valid_ZV = ZV_all_roots.filter(z => z > vapEosCoeffs.B_mix_eos + MIN_MOLE_FRACTION);
                if (valid_ZV.length > 0) ZV_final = valid_ZV[valid_ZV.length - 1]; // Largest positive root for vapor
            }

            if (ZL_final === null || ZV_final === null) {
                break; 
            }

            const phi_L = calculateSrkFugacityCoefficientsTernary(x, T_K_current, P_system_Pa, ZL_final, liqMix.mixture_a, liqMix.mixture_b, liqMix.component_a_values, liqMix.component_b_values, ternarySrkKijParams);
            const phi_V = calculateSrkFugacityCoefficientsTernary(y_vap, T_K_current, P_system_Pa, ZV_final, vapMix.mixture_a, vapMix.mixture_b, vapMix.component_a_values, vapMix.component_b_values, ternarySrkKijParams);

            if (phi_L.some(isNaN) || phi_V.some(isNaN)) {
                ZL_final = null; ZV_final = null;
                break; 
            }

            const new_K_values = phi_L.map((phiL_i, i) => phi_V[i] === 0 ? Infinity : phiL_i / phi_V[i]);
            if (new_K_values.some(k => isNaN(k) || !isFinite(k))) {
                ZL_final = null; ZV_final = null;
                if (new_K_values.some(k => k === Infinity) && T_iter > 0) {
                    T_K_current *= 1.01; 
                }
                break;
            }
            
            const K_diff_sq_sum = K_values.reduce((sum, oldK, i) => sum + (new_K_values[i] - oldK)**2, 0);
            K_values = new_K_values;
            sum_Kx_for_T_loop = K_values.reduce((sum, Ki, i) => sum + Ki * x[i], 0);

            if (Math.sqrt(K_diff_sq_sum) < K_VALUE_TOLERANCE && K_iter > 0) {
                break;
            }
        } // End K-value loop

        if (ZL_final === null || ZV_final === null) {
            T_K_current *= (T_iter % 2 === 0 ? 1.03 : 0.97); // Try to break out of problematic T
            if (T_K_current <= 100 || isNaN(T_K_current)) T_K_current = initialTempGuess_K * (0.8 + Math.random() * 0.4); // Reset if T is too low/NaN
            continue;
        }

        const error_T_loop = sum_Kx_for_T_loop - 1.0;

        if (Math.abs(error_T_loop) < BUBBLE_T_CONVERGENCE_TOLERANCE) {
            const final_y_vap = normalizeMoleFractions(K_values.map((Ki, i) => Ki * x[i]));
            // ... (trivial solution check, if any, remains here) ...
            return { T_K: T_K_current, y_vap: final_y_vap, K_values, ZL: ZL_final!, ZV: ZV_final! }; // Return T_K_current
        }

        const target_T_ratio = 1.0 / sum_Kx_for_T_loop;
        let damped_change = (target_T_ratio - 1.0) * (T_iter < 5 ? 0.1 : 0.25);
        let T_ratio = 1 + damped_change;
        T_ratio = Math.max(0.90, Math.min(1.10, T_ratio)); // Limit change
        const T_K_new = T_K_current * T_ratio;

        if (Math.abs(T_K_new - T_K_current) < BUBBLE_T_TOLERANCE_K && T_iter > 2) {
            if (Math.abs(error_T_loop) < BUBBLE_T_CONVERGENCE_TOLERANCE * 20) { // Looser tolerance if T is stable
                const final_y_vap = normalizeMoleFractions(K_values.map((Ki, i) => Ki * x[i]));
                return { T_K: T_K_current, y_vap: final_y_vap, K_values, ZL: ZL_final!, ZV: ZV_final! }; // Return T_K_current
            }
        }
        T_K_current = Math.max(150, Math.min(800, T_K_new)); // Ensure T_K_current stays within reasonable bounds & Update T_K_current
        if (isNaN(T_K_current) || T_K_current <= 0) {
            T_K_current = initialTempGuess_K * (0.7 + Math.random() * 0.6); // More aggressive reset if T becomes invalid
            if (T_K_current <= 0) T_K_current = 273.15; // Absolute fallback
        }
    }
    // console.warn(`SRK Bubble T: Max T-iterations reached for x=${x.map(f=>f.toFixed(3))}. Last T_K=${T_K_current.toFixed(1)}, Sum Kx=${sum_Kx_for_T_loop?.toFixed(4)}`);
    return null;
}

function calculateVaporCompositionEos(x_liq: number[], K_values: number[]): number[] | null {
    const y_unnormalized = x_liq.map((xi: number, i: number) => K_values[i] * xi);
    const sum_y_unnormalized = y_unnormalized.reduce((s: number, val: number) => s + val, 0);
    if (sum_y_unnormalized === 0 || isNaN(sum_y_unnormalized)) return null;
    const y_vap = y_unnormalized.map((val: number) => val / sum_y_unnormalized);
    return normalizeMoleFractions(y_vap);
}

function residue_curve_ode_dxdxi_srk(
    x_liq_current: number[], P_system_Pa: number, components: CompoundData[],
    ternarySrkParams: TernarySrkParams, initialTempGuess_K: number
): { dxdxi: number[], T_bubble_K: number, y_vap: number[] } | null {
    if (!x_liq_current || x_liq_current.length !== components.length || x_liq_current.some(isNaN)) return null;
    const x_norm = normalizeMoleFractions(x_liq_current);
    
    const bubbleResult = findBubblePointTemperatureTernarySrk(x_norm, P_system_Pa, components, ternarySrkParams, initialTempGuess_K);
    if (!bubbleResult) return null;

    const { T_K: T_bubble_K, K_values } = bubbleResult; // y_vap is also in bubbleResult if needed directly

    const y_vap = calculateVaporCompositionEos(x_norm, K_values);
    if (!y_vap) return null;

    const dxdxi = [x_norm[0] - y_vap[0], x_norm[1] - y_vap[1], x_norm[2] - y_vap[2]];
    return { dxdxi, T_bubble_K, y_vap };
}

// --- Direct Azeotrope Solver for SRK ---
export function findAzeotropeSrk(
    x: number[],
    P_system_Pa: number,
    components: CompoundData[],
    ternarySrkParams: TernarySrkParams,
    initialGuess: number[]
): AzeotropeResult | null {
    const x_norm = normalizeMoleFractions(x);

    const f = (vars: number[]): number[] | null => { // This function already returns number[] | null
        const K_values = vars;
        const y_vap = calculateVaporCompositionEos(x_norm, K_values);
        if (!y_vap) return null;
        return [y_vap[0] - K_values[0], y_vap[1] - K_values[1], y_vap[2] - K_values[2]];
    };

    const jacobian = (vars: number[]): number[][] => { // This function returns number[][]
        const K_values = vars;
        const dFdK: number[][] = [];
        for (let i = 0; i < 3; i++) {
            dFdK[i] = [0, 0, 0];
            const y_vap_i = calculateVaporCompositionEos(x_norm.map((val, idx) => (idx === i ? val + 1e-6 : val)), K_values);
            if (y_vap_i) {
                for (let j = 0; j < 3; j++) {
                    dFdK[i][j] = (y_vap_i[j] - K_values[j]) / 1e-6;
                }
            }
        }
        return dFdK;
    };

    const { vars, converged, iterations, errorNorm, message } = solveNewtonRaphsonSystem(f, jacobian, initialGuess); // Call is now compatible

    if (!converged) {
        console.warn(`Azeotrope solver did not converge: ${message}`);
        return null;
    }

    const K_values = vars;
    const y_vap = calculateVaporCompositionEos(x_norm, K_values);
    if (!y_vap) return null;

    const ZL_final = 1; // Placeholder, calculate if needed
    const ZV_final = 1; // Placeholder, calculate if needed

    return { T_K: 0, x, y: y_vap, K_values, ZL: ZL_final, ZV: ZV_final }; // T_K is not calculated in this context
}


export async function simulateSingleResidueCurveODE_Srk(
    initial_x: number[], P_system_Pa: number, componentsData: CompoundData[],
    ternarySrkParams: TernarySrkParams, d_xi_step: number,
    max_steps_per_direction: number, initialTemp_K: number
): Promise<ResidueCurveODE | null> {

    const run_direction = async (start_x: number[], forward: boolean, current_T_guess: number): Promise<ResidueCurveODE> => {
        const curve: ResidueCurveODE = [];
        let current_x = [...start_x];
        let last_T_K = current_T_guess;

        for (let loop_step = 0; loop_step < max_steps_per_direction; loop_step++) {
            const ode_result = residue_curve_ode_dxdxi_srk(current_x, P_system_Pa, componentsData, ternarySrkParams, last_T_K);
            if (!ode_result) break;
            const { dxdxi, T_bubble_K } = ode_result; // y_vap is also available if needed
            last_T_K = T_bubble_K;

            const currentPointData: ResidueCurvePointODE = {
                x: normalizeMoleFractions(current_x), 
                T_K: T_bubble_K,
                step: (forward ? 1 : -1) * loop_step * d_xi_step,
            };
            
            if (curve.length > 0) {
                const lastPoint = curve[curve.length - 1];
                const x_diff_sq = lastPoint.x.reduce((sum, val, i) => sum + (val - currentPointData.x[i])**2, 0);
                if (Math.sqrt(x_diff_sq) < 1e-7 && Math.abs(lastPoint.T_K - currentPointData.T_K) < 1e-4) {
                    // Stagnation detected
                    break; 
                } else {
                    curve.push(currentPointData);
                }
            } else {
                curve.push(currentPointData);
            }

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
                    });
                }
                break;
            }
            if (loop_step > 0 && curve.length > 0) {
                const lastAddedPoint = curve[curve.length - 1];
                const x_before_euler_step = curve.length > 1 ? curve[curve.length - 2].x : start_x;
                
                const diff_sq_from_prev_step = lastAddedPoint.x.reduce((sum, px, i) => sum + (px - x_before_euler_step[i]) ** 2, 0);

                if (Math.sqrt(diff_sq_from_prev_step) < d_xi_step * 1e-5) { 
                    break; 
                }
            }
        }
        return curve;
    };

    if (!initial_x || initial_x.length !== componentsData.length || initial_x.some(isNaN)) return null;
    const initial_x_norm = normalizeMoleFractions(initial_x);
    
    const initialBubbleResult = findBubblePointTemperatureTernarySrk(initial_x_norm, P_system_Pa, componentsData, ternarySrkParams, initialTemp_K);
    if (!initialBubbleResult) return null;
    const T_start_K = initialBubbleResult.T_K;

    const curve_forward = await run_direction(initial_x_norm, true, T_start_K);
    const curve_backward = await run_direction(initial_x_norm, false, T_start_K);

    const combined_curve: ResidueCurveODE = [ ...curve_backward.reverse().slice(1), ...curve_forward ];
    
    if (combined_curve.length < 2) return curve_forward.length > 0 ? curve_forward : (curve_backward.length > 0 ? curve_backward.reverse() : null);
    return combined_curve;
}
