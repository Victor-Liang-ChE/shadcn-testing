import type { SupabaseClient } from '@supabase/supabase-js';
import type { CompoundData, PrPureComponentParams, BubbleDewResult, AntoineParams } from './vle-types'; // Import shared types
import { calculatePsat_Pa as importedCalculatePsat_Pa } from './vle-calculations-unifac'; // Re-use Psat calculation for initial guesses if needed

// Define R_gas_const_J_molK locally or import from a valid common source if available
export const R_gas_const_J_molK = 8.31446261815324; // J/mol·K

// Re-export calculatePsat_Pa if it's used by other modules importing from this file
export { importedCalculatePsat_Pa as calculatePsat_Pa };

// Solves Z^3 + p*Z^2 + q*Z + r = 0
// Replaces the placeholder solveCubicEOS
function solveCubicRealRoots(p: number, q: number, r: number): number[] {
    // Convert to depressed cubic: y^3 + Ay + B = 0
    // Z = y - p/3
    const A_dep = q - (p * p) / 3.0;
    const B_dep = (2.0 * p * p * p) / 27.0 - (p * q) / 3.0 + r;

    let y_roots: number[] = [];

    const discriminant_term1 = (B_dep / 2.0) * (B_dep / 2.0);
    const discriminant_term2 = (A_dep / 3.0) * (A_dep / 3.0) * (A_dep / 3.0);
    let Delta = discriminant_term1 + discriminant_term2;

    if (Math.abs(Delta) < 1e-12) Delta = 0; // Treat very small Delta as zero

    if (Delta > 0) {
        // One real root
        const sqrtDelta = Math.sqrt(Delta);
        const term1_val = -B_dep / 2.0 + sqrtDelta;
        const term2_val = -B_dep / 2.0 - sqrtDelta;
        
        const y1 = Math.cbrt(term1_val) + Math.cbrt(term2_val);
        y_roots.push(y1);
    } else if (Delta === 0) {
        // Multiple real roots (all real, at least two are equal)
        if (Math.abs(A_dep) < 1e-9 && Math.abs(B_dep) < 1e-9) { // A = B = 0 implies y=0 is a triple root
            y_roots.push(0);
        } else {
            const y1_val = -2.0 * Math.cbrt(B_dep / 2.0);
            const y2_val = Math.cbrt(B_dep / 2.0);
            y_roots.push(y1_val);
            // Ensure distinct root before pushing. If y1 and y2 are very close, it's effectively a triple root (or single distinct root if A=0).
            // If A != 0, y1 is a single root and y2 is a double root.
            if (Math.abs(y1_val - y2_val) > 1e-7) { 
                 y_roots.push(y2_val); 
            }
        }
    } else { // Delta < 0
        // Three distinct real roots (trigonometric solution)
        const term_cos_arg_num = (-B_dep / 2.0);
        const term_cos_arg_den = Math.sqrt(-discriminant_term2); // sqrt(-(A_dep/3)^3)
        
        let cos_arg = term_cos_arg_num / term_cos_arg_den;
        // Clamp cos_arg to [-1, 1] due to potential floating point inaccuracies
        if (cos_arg > 1.0) cos_arg = 1.0;
        if (cos_arg < -1.0) cos_arg = -1.0;

        const angle_third = Math.acos(cos_arg) / 3.0;
        const factor = 2.0 * Math.sqrt(-A_dep / 3.0);

        y_roots.push(factor * Math.cos(angle_third));
        y_roots.push(factor * Math.cos(angle_third - (2.0 * Math.PI) / 3.0));
        y_roots.push(factor * Math.cos(angle_third + (2.0 * Math.PI) / 3.0));
    }

    // Convert y roots back to Z roots: Z = y - p/3
    let Z_roots = y_roots.map(y_root => y_root - p / 3.0);
    
    // Filter out non-physical (Z <= 0 for practical purposes, though strictly Z can be small positive)
    // or duplicate roots (within tolerance)
    Z_roots = Z_roots.filter(z => z > 1e-9); // Must be positive

    // Remove duplicates that might arise from numerical precision
    if (Z_roots.length > 1) {
        Z_roots.sort((a, b) => a - b);
        Z_roots = Z_roots.filter((z, index, self) => index === 0 || Math.abs(z - self[index-1]) > 1e-7);
    }
    
    Z_roots.sort((a, b) => a - b); // Final sort

    if (Z_roots.length === 0) {
        console.warn("Cubic EOS solver found no positive real roots after transformation.", {p, q, r, A_dep, B_dep, Delta});
    }
    return Z_roots;
}

// The old solveCubicEOS placeholder is now replaced by solveCubicRealRoots.
// We can keep the name solveCubicEOS if preferred and adapt its signature,
// or directly use solveCubicRealRoots where needed.
// For minimal changes to calculatePrFugacityCoefficients, let's alias or adapt.
// Let's rename solveCubicRealRoots to solveCubicEOS and adjust its signature slightly.

export function solveCubicEOS(coeff_a: number, coeff_b: number, coeff_c: number, coeff_d: number): number[] | null {
    if (Math.abs(coeff_a) < 1e-9) { // Should not happen for EOS where coeff_a is 1
        console.error("solveCubicEOS: Coefficient 'a' is near zero. Cannot solve.");
        return null;
    }
    // Normalize to Z^3 + pZ^2 + qZ + r = 0
    const p = coeff_b / coeff_a;
    const q = coeff_c / coeff_a;
    const r = coeff_d / coeff_a;
    return solveCubicRealRoots(p, q, r);
}

// --- Interfaces ---
export interface PrInteractionParams {
    k_ij: number; // Binary interaction parameter for PR
    k_ji: number; // Usually k_ji = k_ij, but allow for asymmetry if DB stores it
    casn1?: string;
    casn2?: string;
}

// Placeholder for Peng-Robinson fugacity coefficient calculation
export function calculatePengRobinsonPhi(
    T_K: number,
    P_Pa: number,
    compositions: number[], // liquid (x) or vapor (y) mole fractions
    compDataArray: CompoundData[], // Array of component data (each must have prParams)
    interactionParams: PrInteractionParams | any, // Should be TernaryPrParams for ternary systems
                                                // For binary, it's PrInteractionParams {k_ij, k_ji}
    phaseType: 'liquid' | 'vapor'
): number[] | null { // Returns array of fugacity coefficients [phi_1, phi_2, ...] or null on error
    console.warn("calculatePengRobinsonPhi is a placeholder. Returning ideal phi=1 for all components.");
    if (!compDataArray || compDataArray.length === 0) return null;
    if (compositions.length !== compDataArray.length) return null;

    // This is a highly simplified placeholder.
    // A real implementation would:
    // 1. Calculate mixture parameters A_mix, B_mix using mixing rules and k_ij.
    // 2. Formulate the cubic equation for Z: Z^3 + c2*Z^2 + c1*Z + c0 = 0
    //    where c2, c1, c0 are functions of A_mix, B_mix.
    //    (c2 = B_mix - 1; c1 = A_mix - 2*B_mix - 3*B_mix^2; c0 = B_mix^3 + B_mix^2 - A_mix*B_mix)
    // 3. Solve for Z using solveCubicEOS(1, c2, c1, c0).
    // 4. Select the appropriate Z root (smallest for liquid, largest for vapor).
    // 5. Calculate ln(phi_i) for each component using the PR equation for fugacity coefficients.

    // Placeholder: just to demonstrate call to solveCubicEOS
    const B_mix_placeholder = 0.05; // Example value
    const A_mix_placeholder = 0.1;  // Example value
    const c2 = B_mix_placeholder - 1;
    const c1 = A_mix_placeholder - 2 * B_mix_placeholder - 3 * B_mix_placeholder * B_mix_placeholder;
    const c0 = B_mix_placeholder * B_mix_placeholder * B_mix_placeholder + B_mix_placeholder * B_mix_placeholder - A_mix_placeholder * B_mix_placeholder;

    const Z_roots = solveCubicEOS(1, c2, c1, c0); // Removed phaseType

    if (!Z_roots || Z_roots.length === 0) {
        console.error("PR Phi: Failed to get Z roots from placeholder solver.");
        return null;
    }

    let Z_phase;
    if (phaseType === 'liquid') {
        Z_phase = Z_roots[0]; // Smallest root
    } else {
        Z_phase = Z_roots[Z_roots.length - 1]; // Largest root
    }
    // console.log(`PR Phi Placeholder: Phase Z = ${Z_phase}`);

    return compositions.map(() => 1.0); // Return ideal phi
}

/**
 * Calculates Peng-Robinson fugacity coefficients for a binary mixture.
 */
export function calculatePrFugacityCoefficients(
    components: CompoundData[], // Uses imported CompoundData
    x_or_y: number[], // mole fractions [x1, x2] of the phase being calculated
    T_K: number,
    P_Pa: number,
    interactionParams: PrInteractionParams,
    phase: 'liquid' | 'vapor'
): [number, number] | null { // Returns [phi1, phi2] or null on error

    if (components.length !== 2 || x_or_y.length !== 2) {
        console.error("PR Fugacity calculation currently supports only binary mixtures.");
        return null;
    }
    if (!components[0].prParams || !components[1].prParams) {
        console.error("PR parameters (Tc, Pc, omega) missing for one or both components.");
        return null;
    }

    const comp1 = components[0];
    const comp2 = components[1];
    const x1 = x_or_y[0];
    const x2 = x_or_y[1];

    const { k_ij } = interactionParams;

    // Pure component PR parameters
    const prPure = components.map(c => {
        const { Tc_K, Pc_Pa, omega } = c.prParams!;
        const Tr = T_K / Tc_K;
        let kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
        if (omega > 0.49) { // More accurate for higher acentric factors
             kappa = 0.379642 + 1.485030 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega;
        }
        const alpha = Math.pow(1 + kappa * (1 - Math.sqrt(Tr)), 2);
        const ac_i = 0.45723553 * Math.pow(R_gas_const_J_molK * Tc_K, 2) / Pc_Pa * alpha;
        const b_i = 0.07779607 * R_gas_const_J_molK * Tc_K / Pc_Pa;
        return { ac_i, b_i };
    });

    const ac1 = prPure[0].ac_i;
    const b1 = prPure[0].b_i;
    const ac2 = prPure[1].ac_i;
    const b2 = prPure[1].b_i;

    // Mixture parameters
    const b_mix = x1 * b1 + x2 * b2;
    const a12_mix = Math.sqrt(ac1 * ac2) * (1 - k_ij);
    const a_mix = x1 * x1 * ac1 + x2 * x2 * ac2 + 2 * x1 * x2 * a12_mix;

    const A_mix = a_mix * P_Pa / Math.pow(R_gas_const_J_molK * T_K, 2);
    const B_mix = b_mix * P_Pa / (R_gas_const_J_molK * T_K);

    // Solve cubic EOS for Z: Z^3 - (1-B)Z^2 + (A-3B^2-2B)Z - (AB-B^2-B^3) = 0
    const p_coeff = -(1 - B_mix);
    const q_coeff = A_mix - 3 * B_mix * B_mix - 2 * B_mix;
    const r_coeff = -(A_mix * B_mix - B_mix * B_mix - B_mix * B_mix * B_mix);

    const Z_roots = solveCubicEOS(1, p_coeff, q_coeff, r_coeff); // Call the new solver

    let Z: number;
    if (!Z_roots || Z_roots.length === 0) {
        console.warn(`PR EOS: No positive real roots for Z at T=${T_K.toFixed(1)}, P=${(P_Pa/1000).toFixed(1)}kPa. Phase: ${phase}`, {A_mix, B_mix, x1, p_coeff, q_coeff, r_coeff});
        return null;
    } else if (Z_roots.length === 1) {
        Z = Z_roots[0];
        if (phase === 'liquid' && Z <= B_mix) {
            return null; // No distinct liquid phase if only root is non-physical for liquid
        }
    } else { // Multiple real roots (typically 2 or 3 after filtering positive)
        if (phase === 'liquid') {
            const liquidCandidates = Z_roots.filter(r => r > B_mix);
            if (liquidCandidates.length > 0) {
                Z = liquidCandidates[0]; // Smallest root that is > B_mix
            } else {
                return null; // No physical liquid root
            }
        } else { // Vapor phase
            Z = Z_roots[Z_roots.length - 1]; // Largest root for vapor
        }
    }
    
    if (Z <= B_mix && phase === 'liquid') {
        return null; 
    }
    if (Z <= B_mix && phase === 'vapor' && Z_roots.length > 0) {
        return null;
    }

    // Fugacity coefficients
    const term_common = A_mix / (2 * Math.sqrt(2) * B_mix) * Math.log((Z + (1 + Math.sqrt(2)) * B_mix) / (Z + (1 - Math.sqrt(2)) * B_mix));

    const d_a_mix_d_n1 = 2 * x1 * ac1 + 2 * x2 * a12_mix; // Partial derivative of (n_total * a_mix) wrt n1, then divide by n_total
    const d_a_mix_d_n2 = 2 * x2 * ac2 + 2 * x1 * a12_mix; // Partial derivative of (n_total * a_mix) wrt n2, then divide by n_total

    const ln_phi1 = (b1 / b_mix) * (Z - 1) - Math.log(Z - B_mix) - term_common * ( (d_a_mix_d_n1 / a_mix) - (b1 / b_mix) );
    const ln_phi2 = (b2 / b_mix) * (Z - 1) - Math.log(Z - B_mix) - term_common * ( (d_a_mix_d_n2 / a_mix) - (b2 / b_mix) );
    
    const term1_phi = (2 * (x1*ac1 + x2*a12_mix) / a_mix) - (b1/b_mix);
    const term2_phi = (2 * (x1*a12_mix + x2*ac2) / a_mix) - (b2/b_mix);

    const ln_phi1_corrected = (b1 / b_mix) * (Z - 1) - Math.log(Z - B_mix) - term_common * term1_phi;
    const ln_phi2_corrected = (b2 / b_mix) * (Z - 1) - Math.log(Z - B_mix) - term_common * term2_phi;

    if (isNaN(ln_phi1_corrected) || isNaN(ln_phi2_corrected) || !isFinite(ln_phi1_corrected) || !isFinite(ln_phi2_corrected) || Z <= B_mix) { // Added Z <= B_mix check here for robustness
        console.warn(`PR EOS: NaN/Inf fugacity coefficient or Z <= B_mix. Z=${Z.toFixed(4)}, B_mix=${B_mix.toFixed(4)}, T=${T_K.toFixed(1)}, P=${(P_Pa/1000).toFixed(1)}kPa, x1=${x_or_y[0].toFixed(3)}`, {ln_phi1_corrected, ln_phi2_corrected, A_mix, B_mix, phase});
        return null;
    }

    return [Math.exp(ln_phi1_corrected), Math.exp(ln_phi2_corrected)];
}


export async function fetchPrInteractionParams(
    supabase: SupabaseClient,
    casn1_query: string,
    casn2_query: string
): Promise<PrInteractionParams> {
    if (!casn1_query || !casn2_query) {
        console.warn("PR interaction param fetch: CAS numbers for both compounds are required. Defaulting k_ij = k_ji = 0.");
        return { k_ij: 0, k_ji: 0, casn1: casn1_query, casn2: casn2_query };
    }
    if (casn1_query === casn2_query) {
        // For a pure component, interaction parameter is 0
        return { k_ij: 0, k_ji: 0, casn1: casn1_query, casn2: casn2_query };
    }

    console.log(`Fetching PR interaction parameters from "peng - robinson parameters" for query pair: comp1_cas='${casn1_query}', comp2_cas='${casn2_query}'`);
    
    // Assuming table name is "peng - robinson parameters" and columns are "CASN1", "CASN2" (case-sensitive), and k12 (lowercase).
    // Double quotes are necessary for table/column names with spaces, hyphens, or enforced case sensitivity.
    const query = supabase
        .from('peng-robinson parameters') 
        .select('"CASN1", "CASN2", k12') // Select k12 (lowercase as per DDL), "CASN1", "CASN2" (uppercase as per DDL)
        // Ensure "CASN1" and "CASN2" are quoted as they are case-sensitive in the DB (from DDL).
        .or(`and("CASN1".eq.${casn1_query},"CASN2".eq.${casn2_query}),and("CASN1".eq.${casn2_query},"CASN2".eq.${casn1_query})`)
        .limit(1);

    const { data, error } = await query;

    if (error) {
        console.error(`Supabase PR interaction parameter query error from "peng - robinson parameters": ${error.message}`, error);
        throw new Error(`Supabase PR interaction parameter query error: ${error.message}`);
    }

    if (!data || data.length === 0) {
        console.warn(`PR interaction parameter (k12) not found in "peng - robinson parameters" for pair ${casn1_query}/${casn2_query}. Defaulting to k_ij = k_ji = 0.`);
        return { k_ij: 0, k_ji: 0, casn1: casn1_query, casn2: casn2_query };
    }

    const dbRow = data[0];
    // Assuming k12 is the column name for the interaction parameter.
    // If k12 can be null in the database, provide a default (e.g., 0).
    const k_value = typeof dbRow.k12 === 'number' ? dbRow.k12 : 0; 
    
    console.log(`PR_PARAM_FETCHED: For query ${casn1_query}/${casn2_query} from "peng - robinson parameters", DB row (${dbRow.CASN1}/${dbRow.CASN2}) k12=${dbRow.k12}. Assigned k_ij=${k_value.toFixed(4)}, k_ji=${k_value.toFixed(4)}`);

    return {
        k_ij: k_value,
        k_ji: k_value, // Assuming k_ij = k_ji for Peng-Robinson from a single k12 value
        casn1: casn1_query,
        casn2: casn2_query
    };
}

// --- Bubble Point Calculations using Peng-Robinson ---

// Helper function for numerical derivative (central difference)
function numericalDerivative(
    func: (t: number) => number,
    t: number,
    h: number = 1e-3 // for temperature, can be tuned
): number {
    const f_plus_h = func(t + h);
    const f_minus_h = func(t - h);

    if (isNaN(f_plus_h) || isNaN(f_minus_h)) {
        console.warn(`Numerical derivative: NaN encountered at T=${t.toFixed(2)}`);
        return NaN;
    }
    if (Math.abs(2 * h) < 1e-9) { // Avoid division by zero if h is too small
        console.warn(`Numerical derivative: Step size h is too small or zero.`);
        return NaN;
    }
    return (f_plus_h - f_minus_h) / (2 * h);
}


export function calculateBubbleTemperaturePr(
    components: CompoundData[],
    x1_feed: number,
    P_system_Pa: number,
    prInteractionParams: PrInteractionParams,
    initialTempGuess_K: number,
    maxIter: number = 100,
    tolerance_f: number = 1e-6, // for the objective function |f(T)|
    tolerance_T_step: number = 1e-4 // for |T_new - T_old|
): BubbleDewResult | null {
    console.log(`PR Bubble T (Newton): x1=${x1_feed.toFixed(3)}, P=${(P_system_Pa / 1000).toFixed(1)}kPa, Init_T=${initialTempGuess_K.toFixed(1)}K`);
    let T_K = initialTempGuess_K;
    const x = [x1_feed, 1 - x1_feed]; // Liquid phase composition is fixed for bubble point

    if (!components[0]?.prParams || !components[1]?.prParams) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: NaN, P_Pa: P_system_Pa, error: "PR: Pure component PR parameters missing.", calculationType: 'bubbleT' };
    }

    let K_values_iter: number[] = [1.0, 1.0]; // K-values from current iteration
    let y_iter: number[] = [...x];       // y-values from current iteration

    // Objective function: f(T) = sum(K_i * x_i) - 1
    const objectiveFunction = (T_current_K: number): number => {
        // Calculate K-values at T_current_K, P_system_Pa, and for liquid phase x
        // This requires calculating phi_L and phi_V

        // Calculate phi_L for fixed liquid composition x
        const phi_L = calculatePrFugacityCoefficients(components, x, T_current_K, P_system_Pa, prInteractionParams, 'liquid');
        if (!phi_L) {
            console.warn(`PR Bubble T ObjFunc: Failed to get phi_L at T=${T_current_K.toFixed(2)}K`);
            return NaN; // Indicate failure
        }

        // Estimate y based on current K_values_iter (or initial guess for first pass)
        // For the objective function, y is an intermediate needed to get phi_V.
        // The approach is to use K-values from the *previous* successful T step,
        // or initialize y ~ x if it's the very first T.
        // We'll use K_values_iter which gets updated.
        let y_for_phi_V = [K_values_iter[0] * x[0], K_values_iter[1] * x[1]];
        const sumY_phi_V = y_for_phi_V[0] + y_for_phi_V[1];
        if (sumY_phi_V > 1e-9) {
            y_for_phi_V = [y_for_phi_V[0] / sumY_phi_V, y_for_phi_V[1] / sumY_phi_V];
        } else {
            // If K-values are problematic, use x as a guess for y's structure
            y_for_phi_V = [...x];
        }
        
        // Calculate phi_V for the estimated vapor composition y_for_phi_V
        const phi_V = calculatePrFugacityCoefficients(components, y_for_phi_V, T_current_K, P_system_Pa, prInteractionParams, 'vapor');
        if (!phi_V) {
            console.warn(`PR Bubble T ObjFunc: Failed to get phi_V at T=${T_current_K.toFixed(2)}K for y=[${y_for_phi_V[0].toFixed(3)},${y_for_phi_V[1].toFixed(3)}]`);
            return NaN; // Indicate failure
        }

        // Calculate K-values for this specific T_current_K
        const K_at_T = [phi_L[0] / phi_V[0], phi_L[1] / phi_V[1]];
        if (K_at_T.some(k => isNaN(k) || !isFinite(k))) {
            console.warn(`PR Bubble T ObjFunc: Invalid K-values at T=${T_current_K.toFixed(2)}K`);
            return NaN;
        }

        // Store K-values and the y they imply for the next Newton step's phi_V calc or for final output
        (objectiveFunction as any)._last_K_values = K_at_T; // Attach to function object for access
        const sum_Kx_for_y = K_at_T[0]*x[0] + K_at_T[1]*x[1];
        if (sum_Kx_for_y > 1e-9) {
            (objectiveFunction as any)._last_y_values = [K_at_T[0]*x[0]/sum_Kx_for_y, K_at_T[1]*x[1]/sum_Kx_for_y];
        } else {
             (objectiveFunction as any)._last_y_values = [...x];  // Fallback
        }

        return K_at_T[0] * x[0] + K_at_T[1] * x[1] - 1.0;
    };

    // Initialize K_values_iter and y_iter for the first call to objectiveFunction
    // Using simplified K = Psat/P for initial guess
    if (components[0].antoine && components[1].antoine) {
        const Psat1_init = importedCalculatePsat_Pa(components[0].antoine, T_K); // Removed name, assuming calculatePsat_Pa doesn't need it
        const Psat2_init = importedCalculatePsat_Pa(components[1].antoine, T_K); // Removed name
        if (!isNaN(Psat1_init) && !isNaN(Psat2_init) && Psat1_init > 0 && Psat2_init > 0) {
            K_values_iter = [Psat1_init / P_system_Pa, Psat2_init / P_system_Pa];
            const sumY_init = K_values_iter[0]*x[0] + K_values_iter[1]*x[1];
            if(sumY_init > 1e-9) {
                y_iter = [K_values_iter[0]*x[0]/sumY_init, K_values_iter[1]*x[1]/sumY_init];
            }
        }
    }


    for (let iter = 0; iter < maxIter; iter++) {
        const fT = objectiveFunction(T_K);

        // Update K_values_iter and y_iter from the properties calculated within objectiveFunction
        if ((objectiveFunction as any)._last_K_values) {
            K_values_iter = (objectiveFunction as any)._last_K_values;
        }
        if ((objectiveFunction as any)._last_y_values) {
            y_iter = (objectiveFunction as any)._last_y_values;
        }

        if (isNaN(fT)) {
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, iterations: iter, error: "PR BubbleT: Objective function NaN", calculationType: 'bubbleT' };
        }

        if (Math.abs(fT) < tolerance_f) {
            return {
                comp1_feed: x1_feed,
                comp1_equilibrium: Math.min(Math.max(y_iter[0], 0), 1), // y_comp1 from last successful K calc
                T_K: T_K,
                P_Pa: P_system_Pa,
                iterations: iter + 1,
                calculationType: 'bubbleT'
            };
        }

        const dfdT = numericalDerivative(objectiveFunction, T_K, T_K * 1e-4); // Relative step for h

        if (isNaN(dfdT) || Math.abs(dfdT) < 1e-9) { // Avoid division by zero or small derivative
            T_K += (fT > 0 ? -0.1 : 0.1); // Heuristic kick
            if (T_K <= 0) T_K = initialTempGuess_K * (Math.random() * 0.2 + 0.9); // Reset if T goes bad
            continue; // Skip update and retry
        }

        let deltaT = -fT / dfdT;

        // Damping and bounding deltaT
        const max_deltaT_abs = T_K * 0.20; // Max 20% change of current T
        if (Math.abs(deltaT) > max_deltaT_abs) {
            deltaT = Math.sign(deltaT) * max_deltaT_abs;
        }
        if (iter < 3) { // Stronger damping for initial iterations
            deltaT *= 0.5;
        }

        const T_K_new = T_K + deltaT;

        if (T_K_new <= 0) {
            deltaT /= 2; // Try smaller step
            T_K = T_K + deltaT; // Re-calculate T_K with halved deltaT
            if (T_K <= 0) T_K = (T_K + initialTempGuess_K) / 2; // Fallback if still bad
            continue;
        }
        
        // Check for convergence on T_step (absolute and relative)
        if (Math.abs(deltaT) < tolerance_T_step && Math.abs(deltaT) < tolerance_T_step * Math.abs(T_K_new) ) {
             objectiveFunction(T_K_new); // This updates internal _last_K_values and _last_y_values
             K_values_iter = (objectiveFunction as any)._last_K_values || K_values_iter;
             y_iter        = (objectiveFunction as any)._last_y_values || y_iter;

             return {
                 comp1_feed: x1_feed,
                 comp1_equilibrium: Math.min(Math.max(y_iter[0], 0), 1),
                 T_K: T_K_new,
                 P_Pa: P_system_Pa,
                 iterations: iter + 1,
                 calculationType: 'bubbleT'
             };
        }
        T_K = T_K_new;
    }

    return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: "PR: Bubble T max iterations reached", calculationType: 'bubbleT', iterations: maxIter };
}

export function calculateBubblePressurePr(
    components: CompoundData[], // Uses imported CompoundData
    x1_feed: number,
    T_system_K: number,
    prInteractionParams: PrInteractionParams,
    initialPressureGuess_Pa: number,
    maxIter: number = 50,
    tolerance: number = 1e-5 // Relative tolerance for pressure
): BubbleDewResult | null {
    console.log(`PR Bubble P: x1=${x1_feed.toFixed(3)}, T=${T_system_K.toFixed(1)}K, Init_P=${(initialPressureGuess_Pa/1000).toFixed(1)}kPa`);
    let P = initialPressureGuess_Pa;
    const x = [x1_feed, 1 - x1_feed];
    let y = [...x]; // Initial guess for vapor composition

    for (let iter = 0; iter < maxIter; iter++) {
        const phi_L = calculatePrFugacityCoefficients(components, x, T_system_K, P, prInteractionParams, 'liquid');
        if (!phi_L) return { comp1_feed: x1_feed, comp1_equilibrium: NaN, error: `PR: Failed to get phi_L at P=${(P/1000).toFixed(1)}kPa`, calculationType: 'bubbleP', T_K: T_system_K };

        // Inner loop for y_i convergence
        let K = [1,1];
        for (let j = 0; j < 10; j++) { // Max 10 inner iterations for y
            const phi_V = calculatePrFugacityCoefficients(components, y, T_system_K, P, prInteractionParams, 'vapor');
            if (!phi_V) return { comp1_feed: x1_feed, comp1_equilibrium: NaN, error: `PR: Failed to get phi_V at P=${(P/1000).toFixed(1)}kPa, y1=${y[0].toFixed(3)}`, calculationType: 'bubbleP', T_K: T_system_K };
            
            K = [phi_L[0] / phi_V[0], phi_L[1] / phi_V[1]];
            const sumKx = K[0] * x[0] + K[1] * x[1];
             if (sumKx === 0) return { comp1_feed: x1_feed, comp1_equilibrium: NaN, error: "PR: sumKx is zero in Bubble P", calculationType: 'bubbleP', T_K: T_system_K };

            const y_new = [K[0] * x[0] / sumKx, K[1] * x[1] / sumKx];
             if (Math.abs(y_new[0] - y[0]) < 1e-4 && Math.abs(y_new[1] - y[1]) < 1e-4) {
                y = y_new;
                break;
            }
            y = y_new;
        }
        
        const P_calc = (K[0] * x[0] + K[1] * x[1]) * P; // P_calc = sum(y_i_unnormalized) * P_old / sum(y_i_normalized_implicit=1)

        if (Math.abs(P_calc - P) / P < tolerance || Math.abs(P_calc - P) < 1) { // Relative or absolute tolerance
             return {
                comp1_feed: x1_feed,
                comp1_equilibrium: y[0],
                T_K: T_system_K,
                P_Pa: P_calc,
                iterations: iter + 1,
                calculationType: 'bubbleP'
            };
        }
        P = P_calc; // Direct substitution
        if (P <= 0) {
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, error: "PR: Bubble P calculated non-positive pressure", calculationType: 'bubbleP', T_K: T_system_K, iterations: iter + 1 };
        }
    }
    return { comp1_feed: x1_feed, comp1_equilibrium: NaN, error: "PR: Bubble P max iterations reached", calculationType: 'bubbleP', T_K: T_system_K, iterations: maxIter };
}

export async function calculateDewTemperaturePr(
    components: CompoundData[],
    y1_feed: number,
    P_system_Pa: number,
    prInteractionParams: PrInteractionParams,
    initialTempGuess_K: number,
    maxIter: number = 50,
    tolerance: number = 1e-5
): Promise<BubbleDewResult | null> {
    console.log(`PR Dew T: y1=${y1_feed.toFixed(3)}, P=${(P_system_Pa/1000).toFixed(1)}kPa, Init_T=${initialTempGuess_K.toFixed(1)}K`);
    let T_K = initialTempGuess_K;
    const y = [y1_feed, 1 - y1_feed];

    if (!components[0]?.prParams || !components[1]?.prParams) {
        return { comp1_feed: y1_feed, comp1_equilibrium: NaN, T_K: NaN, P_Pa: P_system_Pa, error: "PR: Pure component PR parameters missing.", calculationType: 'dewT' };
    }
    if (!components[0]?.antoine || !components[1]?.antoine) {
        return { comp1_feed: y1_feed, comp1_equilibrium: NaN, T_K: NaN, P_Pa: P_system_Pa, error: "PR: Antoine parameters missing.", calculationType: 'dewT' };
    }

    for (let iter = 0; iter < maxIter; iter++) {
        const Psat1 = importedCalculatePsat_Pa(components[0].antoine!, T_K);
        const Psat2 = importedCalculatePsat_Pa(components[1].antoine!, T_K);
        if (isNaN(Psat1) || isNaN(Psat2)) {
            return { comp1_feed: y1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: "PR: Invalid Psat values.", calculationType: 'dewT' };
        }

        const K1 = Psat1 / P_system_Pa;
        const K2 = Psat2 / P_system_Pa;
        if (K1 <= 0 || K2 <= 0) {
            return { comp1_feed: y1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: "PR: Invalid K values.", calculationType: 'dewT' };
        }

        const sum_y_div_K = y[0] / K1 + y[1] / K2;
        const error_func = sum_y_div_K - 1.0;

        if (Math.abs(error_func) < tolerance) {
            const x1 = y[0] / K1 / sum_y_div_K;
            return { comp1_feed: y1_feed, comp1_equilibrium: x1, T_K: T_K, P_Pa: P_system_Pa, iterations: iter + 1, calculationType: 'dewT' };
        }

        let T_change = -error_func * (T_K * 0.05);
        const max_T_step = 5.0;
        if (Math.abs(T_change) > max_T_step) T_change = Math.sign(T_change) * max_T_step;
        if (iter < 3) T_change *= 0.5;
        let T_K_new = T_K + T_change;
        if (T_K_new <= 0) T_K_new = T_K * 0.9;
        if (T_K_new < 150) T_K_new = 150;
        if (T_K_new > 800) T_K_new = 800;
        T_K = T_K_new;
    }

    return { comp1_feed: y1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: "PR: Dew T max iterations reached", calculationType: 'dewT', iterations: maxIter };
}

export async function calculateDewPressurePr(
    components: CompoundData[],
    y1_feed: number,
    T_system_K: number,
    prInteractionParams: PrInteractionParams,
    initialPressureGuess_Pa: number,
    maxIter: number = 10,
    tolerance: number = 1e-5
): Promise<BubbleDewResult | null> {
    console.log(`PR Dew P: y1=${y1_feed.toFixed(3)}, T=${T_system_K.toFixed(1)}K`);
    const y = [y1_feed, 1 - y1_feed];

    if (!components[0]?.prParams || !components[1]?.prParams) {
        return { comp1_feed: y1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "PR: Pure component PR parameters missing.", calculationType: 'dewP' };
    }
    if (!components[0]?.antoine || !components[1]?.antoine) {
        return { comp1_feed: y1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "PR: Antoine parameters missing.", calculationType: 'dewP' };
    }

    const Psat1 = importedCalculatePsat_Pa(components[0].antoine!, T_system_K);
    const Psat2 = importedCalculatePsat_Pa(components[1].antoine!, T_system_K);
    if (isNaN(Psat1) || isNaN(Psat2)) {
        return { comp1_feed: y1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "PR: Invalid Psat values.", calculationType: 'dewP' };
    }

    const sum_y_div_Psat = y[0] / Psat1 + y[1] / Psat2;
    if (sum_y_div_Psat <= 0) {
        return { comp1_feed: y1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "PR: Invalid sum_y_div_Psat.", calculationType: 'dewP' };
    }

    const P_calc = 1.0 / sum_y_div_Psat;
    if (isNaN(P_calc) || P_calc <= 0) {
        return { comp1_feed: y1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: NaN, error: "PR: Invalid P_calc.", calculationType: 'dewP' };
    }

    const x1 = y[0] * P_calc / Psat1;

    return {
        comp1_feed: y1_feed,
        comp1_equilibrium: x1,
        T_K: T_system_K,
        P_Pa: P_calc,
        iterations: 1,
        calculationType: 'dewP'
    };
}
