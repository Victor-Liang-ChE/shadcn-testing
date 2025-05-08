import type { SupabaseClient } from '@supabase/supabase-js';
// MOVED CompoundData, PrPureComponentParams, BubbleDewResult, AntoineParams to vle-types.ts
import type { CompoundData, PrPureComponentParams, BubbleDewResult, AntoineParams } from './vle-types'; // Import shared types
import { calculatePsat_Pa } from './vle-calculations-unifac'; // Re-use Psat calculation for initial guesses if needed

// --- Constants ---
const R_J_molK = 8.31446261815324; // Gas constant in J/mol·K or Pa·m³/mol·K

// --- Interfaces ---
export interface PrInteractionParams {
    k12: number;
    casn1?: string;
    casn2?: string;
}

// Helper to solve cubic equation: Z^3 + p*Z^2 + q*Z + r = 0
// Returns real roots, sorted. For PR EOS, we expect 1 or 3 real roots.
function solveCubicRealRoots(p: number, q: number, r: number): number[] {
    // Convert to depressed cubic: y^3 + Ay + B = 0
    // Z = y - p/3
    const A = q - (p * p) / 3.0;
    const B = (2.0 * p * p * p) / 27.0 - (p * q) / 3.0 + r;

    let y_roots: number[] = [];

    const discriminant_term1 = (B / 2.0) * (B / 2.0);
    const discriminant_term2 = (A / 3.0) * (A / 3.0) * (A / 3.0);
    let Delta = discriminant_term1 + discriminant_term2; // Changed const to let

    if (Math.abs(Delta) < 1e-12) Delta = 0; // Treat very small Delta as zero

    if (Delta > 0) {
        // One real root
        const sqrtDelta = Math.sqrt(Delta);
        const term1_val = -B / 2.0 + sqrtDelta;
        const term2_val = -B / 2.0 - sqrtDelta;
        
        const y1 = Math.cbrt(term1_val) + Math.cbrt(term2_val);
        y_roots.push(y1);
    } else if (Delta === 0) {
        // Multiple real roots (all real, at least two are equal)
        if (Math.abs(A) < 1e-9 && Math.abs(B) < 1e-9) { // A = B = 0 implies y=0 is a triple root
            y_roots.push(0);
        } else {
            const y1 = -2.0 * Math.cbrt(B / 2.0);
            const y2 = Math.cbrt(B / 2.0);
            y_roots.push(y1);
            if (Math.abs(y1 - y2) > 1e-7) { // Ensure distinct root before pushing
                 y_roots.push(y2); // y2 is a double root if y1 != y2
            } else { // y1 is a triple root (or y2 is, same value)
                // Only one distinct value, already pushed
            }
        }
    } else { // Delta < 0
        // Three distinct real roots (trigonometric solution)
        const term_cos_arg_num = (-B / 2.0);
        const term_cos_arg_den = Math.sqrt(-discriminant_term2); // sqrt(-(A/3)^3)
        
        let cos_arg = term_cos_arg_num / term_cos_arg_den;
        // Clamp cos_arg to [-1, 1] due to potential floating point inaccuracies
        if (cos_arg > 1.0) cos_arg = 1.0;
        if (cos_arg < -1.0) cos_arg = -1.0;

        const angle_third = Math.acos(cos_arg) / 3.0;
        const factor = 2.0 * Math.sqrt(-A / 3.0);

        y_roots.push(factor * Math.cos(angle_third));
        y_roots.push(factor * Math.cos(angle_third - (2.0 * Math.PI) / 3.0));
        y_roots.push(factor * Math.cos(angle_third + (2.0 * Math.PI) / 3.0));
    }

    // Convert y roots back to Z roots: Z = y - p/3
    let Z_roots = y_roots.map(y_root => y_root - p / 3.0);
    
    // Filter out non-physical (Z <= 0) or duplicate roots (within tolerance)
    Z_roots = Z_roots.filter(z => z > 1e-9); // Z must be positive
    
    // Remove duplicates that might arise from numerical precision
    if (Z_roots.length > 1) {
        Z_roots.sort((a, b) => a - b);
        Z_roots = Z_roots.filter((z, index, self) => index === 0 || Math.abs(z - self[index-1]) > 1e-7);
    }
    
    Z_roots.sort((a, b) => a - b); // Final sort

    if (Z_roots.length === 0) {
        // console.warn("Cubic EOS solver found no positive real roots after transformation.", {p, q, r, A_dep: A, B_dep: B, Delta});
    }
    return Z_roots;
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

    const { k12 } = interactionParams;

    // Pure component PR parameters
    const prPure = components.map(c => {
        const { Tc_K, Pc_Pa, omega } = c.prParams!;
        const Tr = T_K / Tc_K;
        let kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
        if (omega > 0.49) { // More accurate for higher acentric factors
             kappa = 0.379642 + 1.485030 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega;
        }
        const alpha = Math.pow(1 + kappa * (1 - Math.sqrt(Tr)), 2);
        const ac_i = 0.45723553 * Math.pow(R_J_molK * Tc_K, 2) / Pc_Pa * alpha;
        const b_i = 0.07779607 * R_J_molK * Tc_K / Pc_Pa;
        return { ac_i, b_i };
    });

    const ac1 = prPure[0].ac_i;
    const b1 = prPure[0].b_i;
    const ac2 = prPure[1].ac_i;
    const b2 = prPure[1].b_i;

    // Mixture parameters
    const b_mix = x1 * b1 + x2 * b2;
    const a12_mix = Math.sqrt(ac1 * ac2) * (1 - k12);
    const a_mix = x1 * x1 * ac1 + x2 * x2 * ac2 + 2 * x1 * x2 * a12_mix;

    const A_mix = a_mix * P_Pa / Math.pow(R_J_molK * T_K, 2);
    const B_mix = b_mix * P_Pa / (R_J_molK * T_K);

    // Solve cubic EOS for Z: Z^3 - (1-B)Z^2 + (A-3B^2-2B)Z - (AB-B^2-B^3) = 0
    const p_coeff = -(1 - B_mix);
    const q_coeff = A_mix - 3 * B_mix * B_mix - 2 * B_mix;
    const r_coeff = -(A_mix * B_mix - B_mix * B_mix - B_mix * B_mix * B_mix);

    const Z_roots = solveCubicRealRoots(p_coeff, q_coeff, r_coeff);

    let Z: number;
    if (Z_roots.length === 0) {
        console.warn(`PR EOS: No positive real roots for Z at T=${T_K.toFixed(1)}, P=${(P_Pa/1000).toFixed(1)}kPa. Phase: ${phase}`, {A_mix, B_mix, x1});
        return null;
    } else if (Z_roots.length === 1) {
        Z = Z_roots[0];
        if (phase === 'liquid' && Z <= B_mix) {
            // console.warn(`PR EOS: Single root Z=${Z.toFixed(4)} <= B_mix=${B_mix.toFixed(4)} for liquid phase. May indicate single (vapor) phase or supercritical.`);
            return null; // No distinct liquid phase if only root is non-physical for liquid
        }
    } else { // Multiple real roots (typically 2 or 3 after filtering positive)
        if (phase === 'liquid') {
            const liquidCandidates = Z_roots.filter(r => r > B_mix);
            if (liquidCandidates.length > 0) {
                Z = liquidCandidates[0]; // Smallest root that is > B_mix
            } else {
                // console.warn(`PR EOS: No roots Z > B_mix for liquid phase. Roots: ${Z_roots.join(', ')}, B_mix=${B_mix.toFixed(4)}`);
                return null; // No physical liquid root
            }
        } else { // Vapor phase
            Z = Z_roots[Z_roots.length - 1]; // Largest root for vapor
        }
    }
    
    // This check is now more integrated into the selection above, but an explicit final check is good.
    if (Z <= B_mix && phase === 'liquid') {
        // This case should ideally be caught by the logic above.
        // If it's reached, it implies an issue or a very specific state.
        // console.warn(`PR EOS: Final selected liquid Z=${Z.toFixed(4)} is still <= B_mix=${B_mix.toFixed(4)}. Conditions might be single phase vapor or supercritical.`);
        return null; 
    }
    // For vapor phase, Z must also be > B_mix. This is generally true for the largest root if B_mix is physical.
    if (Z <= B_mix && phase === 'vapor' && Z_roots.length > 0) {
        // This would be unusual if B_mix is correctly calculated and positive.
        // console.warn(`PR EOS: Selected vapor Z=${Z.toFixed(4)} is <= B_mix=${B_mix.toFixed(4)}. This is unexpected.`);
        // It might imply that even the largest root is non-physical, which points to issues.
        return null;
    }


    // Fugacity coefficients
    const term_common = A_mix / (2 * Math.sqrt(2) * B_mix) * Math.log((Z + (1 + Math.sqrt(2)) * B_mix) / (Z + (1 - Math.sqrt(2)) * B_mix));

    const d_a_mix_d_n1 = 2 * x1 * ac1 + 2 * x2 * a12_mix; // Partial derivative of (n_total * a_mix) wrt n1, then divide by n_total
    const d_a_mix_d_n2 = 2 * x2 * ac2 + 2 * x1 * a12_mix; // Partial derivative of (n_total * a_mix) wrt n2, then divide by n_total

    const ln_phi1 = (b1 / b_mix) * (Z - 1) - Math.log(Z - B_mix) - term_common * ( (d_a_mix_d_n1 / a_mix) - (b1 / b_mix) );
    const ln_phi2 = (b2 / b_mix) * (Z - 1) - Math.log(Z - B_mix) - term_common * ( (d_a_mix_d_n2 / a_mix) - (b2 / b_mix) );
    
    // Correction for the derivative term in PR fugacity expression:
    // The term ( (2 * sum_j(x_j * A_ij_param)) / A_param - (b_i/b_m) )
    // For component 1: (2 * (x1*ac1 + x2*a12_mix) / a_mix) - (b1/b_mix)
    // For component 2: (2 * (x1*a12_mix + x2*ac2) / a_mix) - (b2/b_mix)

    const term1_phi = (2 * (x1*ac1 + x2*a12_mix) / a_mix) - (b1/b_mix);
    const term2_phi = (2 * (x1*a12_mix + x2*ac2) / a_mix) - (b2/b_mix);

    const ln_phi1_corrected = (b1 / b_mix) * (Z - 1) - Math.log(Z - B_mix) - term_common * term1_phi;
    const ln_phi2_corrected = (b2 / b_mix) * (Z - 1) - Math.log(Z - B_mix) - term_common * term2_phi;


    if (isNaN(ln_phi1_corrected) || isNaN(ln_phi2_corrected) || !isFinite(ln_phi1_corrected) || !isFinite(ln_phi2_corrected)) {
        console.warn(`PR EOS: NaN/Inf fugacity coefficient calculated. Z=${Z.toFixed(4)}, T=${T_K.toFixed(1)}, P=${(P_Pa/1000).toFixed(1)}kPa, x1=${x_or_y[0].toFixed(3)}`, {ln_phi1_corrected, ln_phi2_corrected, A_mix, B_mix, phase});
        return null;
    }

    return [Math.exp(ln_phi1_corrected), Math.exp(ln_phi2_corrected)];
}


export async function fetchPrInteractionParams(
    supabase: SupabaseClient,
    casn1: string,
    casn2: string
): Promise<PrInteractionParams> {
    if (!casn1 || !casn2) {
        console.warn("PR k_ij fetch: CAS numbers for both compounds are required. Defaulting k_ij = 0.");
        return { k12: 0 };
    }
    // If CAS numbers are the same, it's a pure component, k_ij is not applicable in the same way.
    // For pure components, interaction parameters are not used.
    if (casn1 === casn2) {
        // console.log("PR k_ij fetch: Same CAS numbers, k_ij not applicable for pure component interaction.");
        return { k12: 0 }; // Or handle as appropriate for your system
    }

    console.log(`Fetching PR interaction parameters for pair: casn1='${casn1}', casn2='${casn2}'`);

    // Ensure the table name here EXACTLY matches your database table name.
    // If your table name is "peng - robinson parameters" (with spaces and hyphen), use that.
    // It's generally better to use underscores_in_table_names.
    const query = supabase
        .from('peng-robinson parameters') // Corrected table name
        .select('"k12", "CASN1", "CASN2"') // Ensure column names match DB
        .or(`and(CASN1.eq.${casn1},CASN2.eq.${casn2}),and(CASN1.eq.${casn2},CASN2.eq.${casn1})`)
        .limit(1);

    const { data, error } = await query;

    console.log("Supabase PR interaction query response - data:", JSON.stringify(data, null, 2));
    console.log("Supabase PR interaction query response - error:", JSON.stringify(error, null, 2));

    if (error) {
        console.error(`Supabase PR interaction query raw error object:`, error);
        throw new Error(`Supabase PR interaction query error: ${error.message}`);
    }

    if (!data || data.length === 0) {
        // Default to k12 = 0 if not found, with a warning.
        // Or throw error: throw new Error(`PR interaction parameter k12 not found for pair ${casn1}/${casn2}.`);
        console.warn(`PR interaction parameter k12 not found for pair ${casn1}/${casn2}. Defaulting to k12 = 0.`);
        return { k12: 0, casn1, casn2 };
    }
    
    const params = data[0];
    return {
        k12: typeof params.k12 === 'number' ? params.k12 : 0, // Default to 0 if k12 is null/undefined
        casn1: params.CASN1,
        casn2: params.CASN2,
    };
}

// --- Bubble Point Calculations using Peng-Robinson ---

export function calculateBubbleTemperaturePr(
    components: CompoundData[],
    x1_feed: number,
    P_system_Pa: number,
    prInteractionParams: PrInteractionParams,
    initialTempGuess_K: number,
    maxIter: number = 100, // Increased maxIter
    tolerance: number = 1e-5 // Tolerance for sum(Ki*xi) - 1
): BubbleDewResult | null {
    console.log(`PR Bubble T: x1=${x1_feed.toFixed(3)}, P=${(P_system_Pa/1000).toFixed(1)}kPa, Init_T=${initialTempGuess_K.toFixed(1)}K`);
    let T_K = initialTempGuess_K;
    const x = [x1_feed, 1 - x1_feed];
    let y = [...x]; // Initial guess for vapor composition

    if (!components[0]?.prParams || !components[1]?.prParams) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: NaN, P_Pa: P_system_Pa, error: "PR: Pure component PR parameters missing.", calculationType: 'bubbleT' };
    }

    let K_values = [1.0, 1.0]; // K_i = y_i / x_i

    for (let iter = 0; iter < maxIter; iter++) {
        // Calculate fugacity coefficients for liquid phase (phi_L)
        const phi_L = calculatePrFugacityCoefficients(components, x, T_K, P_system_Pa, prInteractionParams, 'liquid');
        if (!phi_L) {
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: `PR: Failed to calculate liquid fugacities at T=${T_K.toFixed(1)}K.`, calculationType: 'bubbleT', iterations: iter };
        }

        // Estimate vapor phase composition y_i = K_i * x_i
        // Initial K can be from Raoult's law or previous iteration
        // y_i_unnormalized = K_i_prev * x_i
        // For the first iteration, y can be x or from a simple K_i = Psat_i/P
        if (iter === 0) { // Initial y guess based on simplified K
            const Psat1_est = components[0].antoine ? calculatePsat_Pa(components[0].antoine, T_K) : P_system_Pa;
            const Psat2_est = components[1].antoine ? calculatePsat_Pa(components[1].antoine, T_K) : P_system_Pa;
            const K1_est = Psat1_est / P_system_Pa;
            const K2_est = Psat2_est / P_system_Pa;
            y = [K1_est * x[0], K2_est * x[1]];
            const sumY_est = y[0] + y[1];
            if (sumY_est > 0) { y[0] /= sumY_est; y[1] /= sumY_est; }
            else { y = [...x]; } // Fallback if Psat estimates are bad
        } else {
             y = [K_values[0] * x[0], K_values[1] * x[1]];
             const sumY_iter = y[0] + y[1];
             if (sumY_iter > 0) { y[0] /= sumY_iter; y[1] /= sumY_iter; }
             else { /* problem with K values leading to zero sumY */ }
        }


        // Calculate fugacity coefficients for vapor phase (phi_V) using current y
        const phi_V = calculatePrFugacityCoefficients(components, y, T_K, P_system_Pa, prInteractionParams, 'vapor');
        if (!phi_V) {
            // If vapor phase calculation fails, it might mean T is too low (subcooled liquid)
            // or too high (superheated vapor if y is far off).
            // Try adjusting T: if sum(Ki*xi) was > 1, T might be too high. If < 1, T might be too low.
            // This is a tricky spot. For now, let's assume an error.
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: `PR: Failed to calculate vapor fugacities at T=${T_K.toFixed(1)}K with y=[${y[0].toFixed(3)},${y[1].toFixed(3)}]`, calculationType: 'bubbleT', iterations: iter };
        }

        // Update K-values: K_i = phi_i^L / phi_i^V
        K_values = [phi_L[0] / phi_V[0], phi_L[1] / phi_V[1]];

        if (isNaN(K_values[0]) || isNaN(K_values[1]) || !isFinite(K_values[0]) || !isFinite(K_values[1])) {
             return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_K, P_Pa: P_system_Pa, error: `PR: K-value calculation failed (NaN/Inf) at T=${T_K.toFixed(1)}K. K1=${K_values[0]}, K2=${K_values[1]}`, calculationType: 'bubbleT', iterations: iter };
        }

        // Check for convergence: sum(K_i * x_i) = 1
        const sum_Ki_xi = K_values[0] * x[0] + K_values[1] * x[1];
        const error_func = sum_Ki_xi - 1.0;

        // console.debug(`PR Bubble T iter ${iter}: T=${T_K.toFixed(2)}, K1=${K_values[0].toFixed(3)}, K2=${K_values[1].toFixed(3)}, sum(Kx)=${sum_Ki_xi.toFixed(5)}, err=${error_func.toExponential(3)}`);


        if (Math.abs(error_func) < tolerance) {
            // Converged: y_i = K_i * x_i (already normalized from K calculation if sum_Ki_xi is 1)
            // Re-calculate y based on final K values and normalize
            y = [K_values[0] * x[0], K_values[1] * x[1]]; // sum of these should be 1.0
            // Normalization for safety, though sum_Ki_xi is already ~1
            const sumY_final = y[0] + y[1];
            if (sumY_final > 0) { y[0] /= sumY_final; y[1] /= sumY_final; }


            return {
                comp1_feed: x1_feed,
                comp1_equilibrium: y[0],
                T_K: T_K,
                P_Pa: P_system_Pa,
                iterations: iter + 1,
                calculationType: 'bubbleT'
            };
        }

        // Update Temperature (T_K)
        // This is a simplified update. A more robust method (e.g., Brent's) would be better.
        // If sum_Ki_xi > 1, K values are too high, implies T is too low (vapor pressure too high relative to fugacity effects). Increase T.
        // If sum_Ki_xi < 1, K values are too low, implies T is too high. Decrease T.
        // The sensitivity d(sum Kx)/dT is complex.
        // Using a proportional step based on the error, bounded.
        let T_change = -error_func * 10; // Heuristic: error_func is (sumKx - 1). If sumKx > 1, error_func > 0, T_change is negative (decrease T). This seems reversed.
                                        // If sumKx > 1, means system is "too volatile" at current T, so T is too LOW. Need to INCREASE T.
                                        // So, T_change should be proportional to error_func.
        T_change = error_func * (T_K * 0.1); // Proportional step, scaled by current T, and a factor (e.g., 0.1 or smaller)
                                             // If sumKx > 1 (error_func > 0), T_change > 0, so T increases. This is correct.

        // Bound the temperature change
        const max_T_step = 10.0; // Max 10K change
        const min_T_step = 0.01; // Min 0.01K change if not zero
        
        if (Math.abs(T_change) > max_T_step) {
            T_change = Math.sign(T_change) * max_T_step;
        }
        if (Math.abs(T_change) < min_T_step && error_func !== 0) {
             //T_change = Math.sign(T_change) * min_T_step; // Don't force min step if T_change is naturally small
        }
        
        // Damping factor for stability, especially in early iterations or large errors
        if (iter < 5 || Math.abs(error_func) > 0.1) {
            T_change *= 0.5; // Damp early steps or large error steps
        }


        let T_K_new = T_K + T_change;

        // Prevent excessive jumps or non-physical temperatures
        if (T_K_new <= 0) {
            T_K_new = T_K * 0.9; // Fallback if T goes too low
            console.warn(`PR Bubble T: New T was <=0. Resetting T from ${T_K.toFixed(1)} to ${T_K_new.toFixed(1)}`);
        }
        if (Math.abs(T_K_new - T_K) < tolerance * 0.1 && iter > 3) { // Stalled
             // console.debug(`PR Bubble T: Stalled at T=${T_K.toFixed(3)}, error=${error_func.toExponential(3)}. Perturbing slightly.`);
             // T_K_new = T_K + Math.sign(error_func) * min_T_step * 5; // Try a slightly larger kick
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
