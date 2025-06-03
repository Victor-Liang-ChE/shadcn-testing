// --- Interfaces for Data Structures ---
// MOVED AntoineParams, UnifacGroupComposition, PrPureComponentParams, CompoundData, BubbleDewResult to vle-types.ts
import type {
    AntoineParams,
    UnifacGroupComposition, // Assuming this is now in vle-types.ts or defined locally if very specific
    CompoundData,
    BubbleDewResult
} from './vle-types'; // Import shared types

export interface UnifacParameters {
    Rk: { [subgroupId: number]: number };
    Qk: { [subgroupId: number]: number };
    mainGroupMap: { [subgroupId: number]: number };
    a_mk: Map<string, number>; // Key: `${mainGroup_m}-${mainGroup_k}`
}

// REMOVED: VleResultPoint interface (BubbleDewResult is more comprehensive)


export function calculatePsat_Pa(
    params: AntoineParams,
    T_kelvin: number,
    componentName?: string
): number {
    if (!params) {
        // console.debug(`DEBUG: calculatePsat_Pa: No Antoine params for ${componentName || 'unknown component'} at T=${T_kelvin.toFixed(2)}K`);
        return NaN;
    }
    // console.debug(`DEBUG: calculatePsat_Pa for ${componentName || 'unknown component'} at T=${T_kelvin.toFixed(2)}K with params:`, JSON.parse(JSON.stringify(params)));

    // Basic range check (consider stricter checks or warnings)
    // if (T_kelvin < params.Tmin_K || T_kelvin > params.Tmax_K) {
    //     console.warn(`Temperature ${T_kelvin.toFixed(1)} K potentially out of Antoine range (${params.Tmin_K}-${params.Tmax_K} K) for ${componentName}`);
    // }

    const conversionFactor = params.Units?.toLowerCase() === 'kpa' ? 1000 : 1;
    let P_sat: number;

    if (params.EquationNo === 1 || params.EquationNo === '1') { // log10 form
        const log10P = params.A - params.B / (T_kelvin + params.C);
        P_sat = Math.pow(10, log10P);
    } else {  // Assume ln form (natural logarithm)
        const logP = params.A - params.B / (T_kelvin + params.C);
        P_sat = Math.exp(logP);
    }
    P_sat *= conversionFactor; // Apply conversion factor after pow/exp

    if (isNaN(P_sat) || P_sat < 0) { // Added check for negative Psat
        console.warn(`calculatePsat_Pa for ${componentName || 'unknown'} resulted in NaN or negative. T=${T_kelvin}K, Antoine: A=${params.A}, B=${params.B}, C=${params.C}`);
        return NaN;
    }
    // console.debug(`DEBUG: calculatePsat_Pa (${componentName}): P_sat_final_Pa=${P_sat.toExponential(3)} (factor: ${conversionFactor})`);
    return P_sat;
}


/** Calculates UNIFAC activity coefficients gamma_i for a mixture */
export function calculateUnifacGamma(
    componentData: CompoundData[], // Uses imported CompoundData
    x: number[], // mole fractions [x1, x2, ..., xn]
    T_kelvin: number,
    params: UnifacParameters
): number[] | null { // Changed return type to number[] | null
    console.log("!!!!!!!!!! EXECUTING calculateUnifacGamma (GENERALIZED VERSION) !!!!!!!!!!");

    try {
        const numComponents = componentData.length;
        if (numComponents === 0) throw new Error("UNIFAC calculation requires at least one component.");
        if (x.length !== numComponents) throw new Error("Mismatch between componentData and mole fraction array lengths.");

        if (!params || !params.Rk || !params.Qk || !params.mainGroupMap || !params.a_mk) {
            console.error("DEBUG: Incomplete UNIFAC parameters provided to calculateUnifacGamma:", { paramsExists: !!params, RkExists: !!params?.Rk, QkExists: !!params?.Qk, mainGroupMapExists: !!params?.mainGroupMap, a_mkExists: !!params?.a_mk });
            throw new Error("Incomplete UNIFAC parameters provided.");
        }

        const Z = 10.0; // Coordination number

        // --- Ensure r_i, q_i are calculated ---
        for (let i = 0; i < componentData.length; i++) {
            const comp = componentData[i];
            if (comp.r_i === undefined || comp.q_i === undefined) {
                comp.r_i = 0;
                comp.q_i = 0;
                if (!comp.unifacGroups) {
                    console.error(`DEBUG: Missing UNIFAC groups for ${comp.name} when r_i/q_i calculation needed.`);
                    throw new Error(`Missing UNIFAC groups for ${comp.name}`);
                }
                for (const subgroupIdStr in comp.unifacGroups) {
                    const subgroupId = parseInt(subgroupIdStr);
                    const count = comp.unifacGroups[subgroupId];
                    if (params.Rk[subgroupId] === undefined || params.Qk[subgroupId] === undefined) {
                        console.error(`DEBUG: Missing Rk/Qk for subgroup ${subgroupId} (needed by ${comp.name}).`);
                        throw new Error(`Missing Rk/Qk for subgroup ${subgroupId} needed by ${comp.name}`);
                    }
                    comp.r_i += count * params.Rk[subgroupId];
                    comp.q_i += count * params.Qk[subgroupId];
                }
            }
        }

        const r = componentData.map(c => {
            if (c.r_i === undefined) throw new Error(`r_i not defined for ${c.name}`);
            return c.r_i;
        });
        const q = componentData.map(c => {
            if (c.q_i === undefined) throw new Error(`q_i not defined for ${c.name}`);
            return c.q_i;
        });

        const Z_by_2 = Z / 2.0;

        // --- 1. Combinatorial Part ---
        const ln_gamma_C: number[] = Array(numComponents).fill(0.0); // Generalized initialization

        // Calculate sum_xr and sum_xq for the mixture
        let sum_xr_mix = 0;
        let sum_xq_mix = 0;
        for (let k = 0; k < numComponents; k++) {
            sum_xr_mix += x[k] * r[k];
            sum_xq_mix += x[k] * q[k];
        }

        if (sum_xr_mix < 1e-12 || sum_xq_mix < 1e-12) {
            console.warn("UNIFAC: sum_xr_mix or sum_xq_mix is zero. Gammas might be unreliable.");
        }

        const L_component_values = r.map((ri_val, idx) => Z_by_2 * (ri_val - q[idx]) - (ri_val - 1.0));
        console.log(`UNIFAC_DETAILS: L_comp = [${L_component_values.map(val => val.toFixed(4)).join(', ')}]`);

        for (let i = 0; i < numComponents; i++) {
            if (Math.abs(x[i] - 1.0) < 1e-9 && numComponents > 1) { // Component i is pure in a mixture context (others are zero)
                ln_gamma_C[i] = 0.0;
            } else if (numComponents === 1) { // Single component system
                ln_gamma_C[i] = 0.0;
            } else { // Mixture case
                const Phi_i = (sum_xr_mix > 1e-12) ? (x[i] * r[i]) / sum_xr_mix : 0;
                const Theta_i = (sum_xq_mix > 1e-12) ? (x[i] * q[i]) / sum_xq_mix : 0;

                if (x[i] < 1e-9 || Phi_i < 1e-9 || Theta_i < 1e-9 || sum_xr_mix < 1e-9 || sum_xq_mix < 1e-9) {
                    console.warn(`UNIFAC Combinatorial: Very small x, Phi, or Theta for component ${i} in mixture. x_i=${x[i]}, Phi_i=${Phi_i}, Theta_i=${Theta_i}`);
                    let sum_xL_component = 0;
                    for (let comp_idx = 0; comp_idx < numComponents; comp_idx++) {
                        sum_xL_component += x[comp_idx] * L_component_values[comp_idx];
                    }
                    const term1 = (Phi_i / x[i] > 1e-9) ? Math.log(Phi_i / x[i]) : 0;
                    const term2 = (Theta_i / Phi_i > 1e-9 && q[i] > 1e-9) ? Z_by_2 * q[i] * Math.log(Theta_i / Phi_i) : 0;
                    const term3 = L_component_values[i];
                    const term4 = (Phi_i / x[i]) * sum_xL_component;
                    ln_gamma_C[i] = (isNaN(term1) ? 0 : term1) + (isNaN(term2) ? 0 : term2) + term3 - (isNaN(term4) ? 0 : term4);
                } else {
                    const Phi_i_div_xi = Phi_i / x[i];
                    const Phi_i_div_Theta_i = (Theta_i > 1e-12) ? Phi_i / Theta_i : 0;

                    const term_A_log = Math.log(Phi_i_div_xi);
                    const term_A_val = 1.0 - Phi_i_div_xi;

                    const term_B_log = (Phi_i_div_Theta_i > 1e-9) ? Math.log(Phi_i_div_Theta_i) : 0;
                    const term_B_val = 1.0 - Phi_i_div_Theta_i;

                    ln_gamma_C[i] = term_A_log + term_A_val - Z_by_2 * q[i] * (term_B_log + term_B_val);
                }
            }
        }

        console.log(`UNIFAC_DETAILS: x = [${x.map(val => val.toFixed(4)).join(', ')}]`);
        console.log(`UNIFAC_DETAILS: r = [${r.map(val => val.toFixed(4)).join(', ')}]`);
        console.log(`UNIFAC_DETAILS: q = [${q.map(val => val.toFixed(4)).join(', ')}]`);
        console.log(`UNIFAC_DETAILS: ln_gamma_C = [${ln_gamma_C.map(val => val.toFixed(4)).join(', ')}]`);

        // --- 2. Residual Part ---
        const all_subgroup_ids_present = new Set<number>();
        componentData.forEach(comp => {
            if (comp.unifacGroups) {
                Object.keys(comp.unifacGroups).forEach(id => all_subgroup_ids_present.add(parseInt(id)));
            }
        });

        if (all_subgroup_ids_present.size === 0) {
            console.warn("No UNIFAC groups found for residual part. Returning combinatorial gammas.");
            const comb_gamma_only: number[] = ln_gamma_C.map(lng_c => Math.exp(lng_c));
            if (comb_gamma_only.some(isNaN)) {
                console.error("NaN in combinatorial-only gamma calculation."); return null;
            }
            return comb_gamma_only;
        }

        const subgroupIds = Array.from(all_subgroup_ids_present).sort((a, b) => a - b);
        const numGroups = subgroupIds.length;
        console.log(`UNIFAC_DETAILS: SubgroupIDs for residual:`, subgroupIds);

        const v_ki: number[][] = Array(numComponents).fill(0).map(() => Array(numGroups).fill(0));
        subgroupIds.forEach((sgId, k_idx) => {
            for (let i = 0; i < numComponents; i++) {
                v_ki[i][k_idx] = componentData[i].unifacGroups?.[sgId] || 0;
            }
        });

        let sum_x_vk_total = 0;
        for (let i = 0; i < numComponents; ++i) {
            for (let k_idx = 0; k_idx < numGroups; ++k_idx) {
                sum_x_vk_total += x[i] * v_ki[i][k_idx];
            }
        }

        const X_m: number[] = Array(numGroups).fill(0);
        if (sum_x_vk_total > 1e-12) {
            for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                let sum_x_vm = 0;
                for (let i = 0; i < numComponents; i++) {
                    sum_x_vm += x[i] * v_ki[i][m_idx];
                }
                X_m[m_idx] = sum_x_vm / sum_x_vk_total;
            }
        }

        let sum_XQ = 0;
        for (let n_idx = 0; n_idx < numGroups; ++n_idx) {
            sum_XQ += X_m[n_idx] * params.Qk[subgroupIds[n_idx]];
        }

        const theta_m: number[] = Array(numGroups).fill(0);
        if (sum_XQ > 1e-12) {
            for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                theta_m[m_idx] = (X_m[m_idx] * params.Qk[subgroupIds[m_idx]]) / sum_XQ;
            }
        }

        const psi_matrix: number[][] = Array(numGroups).fill(0).map(() => Array(numGroups).fill(0));
        for (let n_idx = 0; n_idx < numGroups; n_idx++) {
            for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                const mainG_n = params.mainGroupMap[subgroupIds[n_idx]];
                const mainG_m = params.mainGroupMap[subgroupIds[m_idx]];
                const interactionKey_nm = `${mainG_n}-${mainG_m}`;
                const a_nm = params.a_mk.get(interactionKey_nm) ?? 0.0;
                psi_matrix[n_idx][m_idx] = Math.exp(-a_nm / T_kelvin);
            }
        }

        const ln_Gamma_k: number[] = Array(numGroups).fill(0);
        for (let k_idx = 0; k_idx < numGroups; k_idx++) {
            let sum_theta_psi_mk = 0;
            for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                sum_theta_psi_mk += theta_m[m_idx] * psi_matrix[m_idx][k_idx];
            }
            const term1 = (sum_theta_psi_mk > 1e-12) ? Math.log(sum_theta_psi_mk) : -Infinity;

            let term2_sum = 0;
            for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                let sum_theta_psi_nm = 0;
                for (let n_idx = 0; n_idx < numGroups; n_idx++) {
                    sum_theta_psi_nm += theta_m[n_idx] * psi_matrix[n_idx][m_idx];
                }
                if (sum_theta_psi_nm > 1e-12) {
                    term2_sum += (theta_m[m_idx] * psi_matrix[k_idx][m_idx]) / sum_theta_psi_nm;
                }
            }
            ln_Gamma_k[k_idx] = params.Qk[subgroupIds[k_idx]] * (1.0 - (isFinite(term1) ? term1 : 0) - term2_sum);
        }

        const ln_Gamma_k_pure: number[][] = Array(numComponents).fill(0).map(() => Array(numGroups).fill(0));
        for (let i = 0; i < numComponents; i++) {
            let sum_vk_pure = 0;
            for (let k_idx = 0; k_idx < numGroups; ++k_idx) sum_vk_pure += v_ki[i][k_idx];

            const X_m_pure: number[] = Array(numGroups).fill(0);
            if (sum_vk_pure > 1e-12) {
                for (let m_idx = 0; m_idx < numGroups; ++m_idx) X_m_pure[m_idx] = v_ki[i][m_idx] / sum_vk_pure;
            }

            let sum_XQ_pure = 0;
            for (let n_idx = 0; n_idx < numGroups; ++n_idx) {
                sum_XQ_pure += X_m_pure[n_idx] * params.Qk[subgroupIds[n_idx]];
            }

            const theta_m_pure: number[] = Array(numGroups).fill(0);
            if (sum_XQ_pure > 1e-12) {
                for (let m_idx = 0; m_idx < numGroups; ++m_idx) {
                    theta_m_pure[m_idx] = (X_m_pure[m_idx] * params.Qk[subgroupIds[m_idx]]) / sum_XQ_pure;
                }
            }

            for (let k_idx = 0; k_idx < numGroups; k_idx++) {
                let sum_theta_psi_mk_pure = 0;
                for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                    sum_theta_psi_mk_pure += theta_m_pure[m_idx] * psi_matrix[m_idx][k_idx];
                }
                const term1_pure = (sum_theta_psi_mk_pure > 1e-12) ? Math.log(sum_theta_psi_mk_pure) : -Infinity;

                let term2_sum_pure = 0;
                for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                    let sum_theta_psi_nm_pure = 0;
                    for (let n_idx = 0; n_idx < numGroups; n_idx++) {
                        sum_theta_psi_nm_pure += theta_m_pure[n_idx] * psi_matrix[n_idx][m_idx];
                    }
                    if (sum_theta_psi_nm_pure > 1e-12) {
                        term2_sum_pure += (theta_m_pure[m_idx] * psi_matrix[k_idx][m_idx]) / sum_theta_psi_nm_pure;
                    }
                }
                ln_Gamma_k_pure[i][k_idx] = params.Qk[subgroupIds[k_idx]] * (1.0 - (isFinite(term1_pure) ? term1_pure : 0) - term2_sum_pure);
            }
        }

        const ln_gamma_R: number[] = Array(numComponents).fill(0.0);
        for (let i = 0; i < numComponents; i++) {
            for (let k_idx = 0; k_idx < numGroups; k_idx++) {
                if (v_ki[i][k_idx] > 0) {
                    ln_gamma_R[i] += v_ki[i][k_idx] * (ln_Gamma_k[k_idx] - ln_Gamma_k_pure[i][k_idx]);
                }
            }
        }
        console.log(`UNIFAC_DETAILS: ln_gamma_R = [${ln_gamma_R.map(val => val.toFixed(4)).join(', ')}]`);

        const gamma_final: number[] = Array(numComponents).fill(0.0);
        for (let i = 0; i < numComponents; i++) {
            gamma_final[i] = Math.exp(ln_gamma_C[i] + ln_gamma_R[i]);
        }
        console.log(`DEBUG: Final gamma = [${gamma_final.map(g => g?.toFixed(4)).join(', ')}]`);

        if (gamma_final.some(g => isNaN(g) || !isFinite(g))) {
            console.error("UNIFAC calculation resulted in NaN or Infinity final gamma values.", { x, T_kelvin, ln_gamma_C, ln_gamma_R });
            return null;
        }
        return gamma_final;

    } catch (error) {
        console.error("Error during UNIFAC calculation:", error);
        return null;
    }
}

/** Bubble Temperature Calculation using Newton-Raphson */
export function calculateBubbleTemperature(
    components: CompoundData[], // Uses imported CompoundData
    x1_feed: number,
    P_system_Pa: number,
    unifacParams: UnifacParameters,
    initialTempGuess_K: number,
    maxIter: number = 50,
    tolerance: number = 1e-4
): BubbleDewResult | null { // Uses imported BubbleDewResult
    console.debug(`DEBUG: calculateBubbleTemperature called for:`);
    console.debug(`DEBUG:   Component 1: ${components[0]?.name} (x1=${x1_feed.toFixed(4)})`);
    console.debug(`DEBUG:   Component 2: ${components[1]?.name} (x2=${(1.0 - x1_feed).toFixed(4)})`);
    console.debug(`DEBUG:   Target Pressure: ${P_system_Pa.toFixed(0)} Pa`);
    console.debug(`DEBUG:   Initial Temp Guess: ${initialTempGuess_K.toFixed(2)} K`);

    const x = [x1_feed, 1.0 - x1_feed];

    const objectiveFunction = (T_kelvin: number): number => {
        const psat = [
            calculatePsat_Pa(components[0].antoine!, T_kelvin, components[0].name),
            calculatePsat_Pa(components[1].antoine!, T_kelvin, components[1].name)
        ];
        const gammaResult = calculateUnifacGamma(components, x, T_kelvin, unifacParams);

        if (!gammaResult || isNaN(psat[0]) || isNaN(psat[1])) {
            return NaN; // Indicate failure to evaluate objective
        }
        const gamma = gammaResult;
        const P_calc_Pa = x[0] * gamma[0] * psat[0] + x[1] * gamma[1] * psat[1];
        return P_calc_Pa - P_system_Pa;
    };

    const derivativeObjectiveFunction = (T_kelvin: number): number => {
        const h = 1e-4; // Small step for numerical differentiation
        const f_plus_h = objectiveFunction(T_kelvin + h);
        const f_minus_h = objectiveFunction(T_kelvin - h);
        if (isNaN(f_plus_h) || isNaN(f_minus_h)) {
            return NaN;
        }
        return (f_plus_h - f_minus_h) / (2 * h); // Corrected derivative
    };

    let T_current = initialTempGuess_K;
    for (let i = 0; i < maxIter; i++) {
        const fT = objectiveFunction(T_current);
        if (isNaN(fT)) {
            console.error(`BubbleT: Objective function returned NaN at T=${T_current.toFixed(2)}K, iter=${i}`);
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_current, P_Pa: P_system_Pa, iterations: i, error: "Objective function NaN", calculationType: 'bubbleT' };
        }

        if (Math.abs(fT) < P_system_Pa * tolerance) { // Relative tolerance
            const gammaConverged = calculateUnifacGamma(components, x, T_current, unifacParams);
            const psatConverged = [
                calculatePsat_Pa(components[0].antoine!, T_current, components[0].name),
                calculatePsat_Pa(components[1].antoine!, T_current, components[1].name)
            ];
            if (!gammaConverged || isNaN(psatConverged[0]) || isNaN(psatConverged[1])) {
                return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_current, P_Pa: P_system_Pa, iterations: i, error: "Final property calculation failed", calculationType: 'bubbleT' };
            }
            const P_final_calc = x[0] * gammaConverged[0] * psatConverged[0] + x[1] * gammaConverged[1] * psatConverged[1];
            const y_comp1_final = (x[0] * gammaConverged[0] * psatConverged[0]) / P_final_calc; // Corrected calculation
            return {
                comp1_feed: x1_feed,
                comp1_equilibrium: Math.min(Math.max(y_comp1_final, 0), 1), // y_comp1
                T_K: T_current,
                P_Pa: P_system_Pa, // P_final_calc could also be used here if more accurate
                iterations: i,
                calculationType: 'bubbleT'
            };
        }

        const dfT = derivativeObjectiveFunction(T_current);
        if (isNaN(dfT) || Math.abs(dfT) < 1e-10) {
            console.error(`BubbleT: Derivative is NaN or too small at T=${T_current.toFixed(2)}K, dfT=${dfT?.toExponential(2)}, iter=${i}`);
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_current, P_Pa: P_system_Pa, iterations: i, error: "Derivative issue", calculationType: 'bubbleT' };
        }

        const T_next = T_current - fT / dfT; // Corrected Newton step
        if (Math.abs(T_next - T_current) < tolerance * Math.abs(T_current) + tolerance) { // Combined relative/absolute step tolerance
             T_current = T_next; // One last update before final calc
             const gammaConverged = calculateUnifacGamma(components, x, T_current, unifacParams);
             const psatConverged = [
                 calculatePsat_Pa(components[0].antoine!, T_current, components[0].name),
                 calculatePsat_Pa(components[1].antoine!, T_current, components[1].name)
             ];
             if (!gammaConverged || isNaN(psatConverged[0]) || isNaN(psatConverged[1])) {
                 return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_current, P_Pa: P_system_Pa, iterations: i, error: "Final property calculation failed on T_step convergence", calculationType: 'bubbleT' };
             }
             const P_final_calc = x[0] * gammaConverged[0] * psatConverged[0] + x[1] * gammaConverged[1] * psatConverged[1];
             const y_comp1_final = (x[0] * gammaConverged[0] * psatConverged[0]) / P_final_calc; // Corrected calculation
             return {
                 comp1_feed: x1_feed,
                 comp1_equilibrium: Math.min(Math.max(y_comp1_final, 0), 1),
                 T_K: T_current,
                 P_Pa: P_system_Pa,
                 iterations: i,
                 calculationType: 'bubbleT'
             };
        }
        T_current = T_next;

        if (T_current < 150 || T_current > 700) { // Basic Temperature bounds
             console.warn(`BubbleT: Temperature ${T_current.toFixed(2)}K out of bounds [150-700K] at iter=${i}.`);
            return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_current, P_Pa: P_system_Pa, iterations: i, error: "Temperature out of bounds", calculationType: 'bubbleT' };
        }
    }
    console.warn(`BubbleT: Failed to converge for x1=${x1_feed.toFixed(4)} after ${maxIter} iterations. Last T=${T_current.toFixed(2)}K, f(T)=${objectiveFunction(T_current)?.toExponential(2)}`);
    return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_current, P_Pa: P_system_Pa, iterations: maxIter, error: "Max iterations reached", calculationType: 'bubbleT' };
}

/** Bubble Pressure Calculation */
export function calculateBubblePressure(
    components: CompoundData[], // Uses imported CompoundData
    x1_feed: number,
    T_system_K: number,
    unifacParams: UnifacParameters,
    initialPressureGuess_Pa: number,
    maxIter: number = 20, // Reduced iterations for P, often more stable
    tolerance: number = 1e-5 // Relative tolerance for pressure
): BubbleDewResult | null { // Uses imported BubbleDewResult
    const x = [x1_feed, 1.0 - x1_feed];

    const psat = [ // Psat is constant for a given T_target_K
        calculatePsat_Pa(components[0].antoine!, T_system_K, components[0].name),
        calculatePsat_Pa(components[1].antoine!, T_system_K, components[1].name)
    ];
    const gammaResult = calculateUnifacGamma(components, x, T_system_K, unifacParams);

    if (!gammaResult || isNaN(psat[0]) || isNaN(psat[1]) || psat[0] <=0 || psat[1] <=0) {
        return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, error: "Initial property calculation failed (Psat/Gamma)", calculationType: 'bubbleP' };
    }
    const gamma = gammaResult;

    // Direct calculation for P_bubble
    const P_bubble_Pa = x[0] * gamma[0] * psat[0] + x[1] * gamma[1] * psat[1];

    if (isNaN(P_bubble_Pa) || P_bubble_Pa <= 0) {
         return { comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: T_system_K, P_Pa: P_bubble_Pa, error: "P_bubble calculation resulted in NaN or non-positive value", calculationType: 'bubbleP' };
    }

    const y_comp1_final = (x[0] * gamma[0] * psat[0]) / P_bubble_Pa;
    return {
        comp1_feed: x1_feed,
        comp1_equilibrium: Math.min(Math.max(y_comp1_final, 0), 1), // y_comp1
        T_K: T_system_K,
        P_Pa: P_bubble_Pa,
        iterations: 1, // Direct calculation
        calculationType: 'bubbleP'
    };
}


/** Dew Temperature Calculation */
export function calculateDewTemperature(
    componentData: CompoundData[],
    y_comp1: number, // Mole fraction of component 1 in the vapor phase
    P_target_Pa: number,
    unifacParams: UnifacParameters,
    initialTempGuess_K: number,
    max_iter: number = 100,
    tolerance: number = 1e-6 // Tolerance for T
): BubbleDewResult | null {
    const y = [y_comp1, 1.0 - y_comp1];
    let T_current = initialTempGuess_K;
    let x_iter = [...y]; // Initial guess for liquid composition

    const objectiveFunctionBuilder = (current_y: number[]) => (T_k: number): number => {
        // Inner loop for x and gamma consistency at T_k
        let x_local = [...current_y]; // Start with y as guess for x
        let K_values_local: number[] = [1.0, 1.0];

        for (let j = 0; j < 20; j++) { // Max 20 inner iterations for x
            const gamma_local = calculateUnifacGamma(componentData, x_local, T_k, unifacParams);
            const psat_local = [
                calculatePsat_Pa(componentData[0].antoine!, T_k, componentData[0].name),
                calculatePsat_Pa(componentData[1].antoine!, T_k, componentData[1].name)
            ];

            if (!gamma_local || psat_local.some(p => isNaN(p) || p <= 0)) return NaN; // Cannot evaluate

            K_values_local = [
                (gamma_local[0] * psat_local[0]) / P_target_Pa,
                (gamma_local[1] * psat_local[1]) / P_target_Pa
            ];

            if (K_values_local.some(k => k <= 0 || isNaN(k))) return NaN;

            const x_new_temp = [current_y[0] / K_values_local[0], current_y[1] / K_values_local[1]];
            const sum_x_new_temp = x_new_temp[0] + x_new_temp[1];

            if (sum_x_new_temp <= 1e-9 || isNaN(sum_x_new_temp)) return NaN;

            const x_next_local = [x_new_temp[0] / sum_x_new_temp, x_new_temp[1] / sum_x_new_temp];
            const diff_x = Math.abs(x_next_local[0] - x_local[0]);
            x_local = x_next_local;
            if (diff_x < 1e-5) break; // x converged
        }
        // Store converged x for this T_k to be used by outer loop if needed
        (objectiveFunctionBuilder as any).last_x_iter = x_local; // Attach to function object for access
        return x_local[0] + x_local[1] - 1.0; // sum(x_i) - 1 = 0
    };
    
    const objectiveFuncInstance = objectiveFunctionBuilder(y);

    const derivativeObjectiveFunction = (T_k_deriv: number): number => {
        const h = 1e-3; // Adjusted h for potentially sensitive function
        const f_plus_h = objectiveFuncInstance(T_k_deriv + h);
        const f_minus_h = objectiveFuncInstance(T_k_deriv - h);
        if (isNaN(f_plus_h) || isNaN(f_minus_h)) {
            return NaN;
        }
        return (f_plus_h - f_minus_h) / (2 * h); // Corrected derivative
    };

    for (let i = 0; i < max_iter; i++) {
        const fT = objectiveFuncInstance(T_current);
         x_iter = (objectiveFunctionBuilder as any).last_x_iter || x_iter; // Get the x used for fT calculation

        if (isNaN(fT)) {
            return { comp1_feed: y_comp1, comp1_equilibrium: NaN, T_K: T_current, P_Pa: P_target_Pa, iterations: i, error: "DewT: Objective function NaN (fT)", calculationType: 'dewT' };
        }

        if (Math.abs(fT) < tolerance) { // Absolute tolerance for sum(x_i) - 1
            return {
                comp1_feed: y_comp1,
                comp1_equilibrium: Math.min(Math.max(x_iter[0], 0), 1), // x_comp1
                T_K: T_current,
                P_Pa: P_target_Pa,
                iterations: i,
                calculationType: 'dewT'
            };
        }

        const dfT = derivativeObjectiveFunction(T_current);
        if (isNaN(dfT) || Math.abs(dfT) < 1e-10) {
             return { comp1_feed: y_comp1, comp1_equilibrium: NaN, T_K: T_current, P_Pa: P_target_Pa, iterations: i, error: "DewT: Derivative issue", calculationType: 'dewT' };
        }

        const T_next = T_current - fT / dfT; // Corrected Newton step
        if (Math.abs(T_next - T_current) < tolerance * Math.abs(T_current) + tolerance) {
             T_current = T_next; // Last update
             // Recalculate x_iter for the final T_current
             objectiveFuncInstance(T_current); // This will update .last_x_iter
             x_iter = (objectiveFunctionBuilder as any).last_x_iter || x_iter;

            return {
                comp1_feed: y_comp1,
                comp1_equilibrium: Math.min(Math.max(x_iter[0], 0), 1),
                T_K: T_current,
                P_Pa: P_target_Pa,
                iterations: i +1, // i is 0-indexed, so i+1 iterations
                calculationType: 'dewT'
            };
        }
        T_current = T_next;

        if (T_current < 150 || T_current > 700) {
            return { comp1_feed: y_comp1, comp1_equilibrium: NaN, T_K: T_current, P_Pa: P_target_Pa, iterations: i, error: "DewT: Temperature out of bounds", calculationType: 'dewT' };
        }
    }
    return { comp1_feed: y_comp1, comp1_equilibrium: NaN, T_K: T_current, P_Pa: P_target_Pa, iterations: max_iter, error: "DewT: Max iterations reached", calculationType: 'dewT' };
}


/** Dew Pressure Calculation */
export function calculateDewPressure(
    componentData: CompoundData[],
    y_comp1: number, // Mole fraction of component 1 in the vapor phase
    T_target_K: number,
    unifacParams: UnifacParameters,
    initialPressureGuess_Pa: number,
    max_iter: number = 50,
    tolerance: number = 1e-6
): BubbleDewResult | null {
    const y = [y_comp1, 1.0 - y_comp1];
    let P_current_Pa = initialPressureGuess_Pa;
    let x_current = [...y]; // Initial guess for liquid composition

    const psat_targetT = [ // Psat is constant for a given T_target_K
        calculatePsat_Pa(componentData[0].antoine!, T_target_K, componentData[0].name),
        calculatePsat_Pa(componentData[1].antoine!, T_target_K, componentData[1].name)
    ];
     if (psat_targetT.some(p => isNaN(p) || p <=0)) {
        return { comp1_feed: y_comp1, comp1_equilibrium: NaN, T_K: T_target_K, P_Pa: P_current_Pa, error: "Initial Psat failed for DewP", calculationType: 'dewP' };
    }

    for (let i = 0; i < max_iter; i++) {
        // Inner loop to get consistent x and gamma for current P_current_Pa
        let gamma_current_iter = [1.0, 1.0]; // Start with ideal gamma for inner loop
        for (let j = 0; j < 20; j++) { // Max 20 inner iterations
            const K_iter = [
                (gamma_current_iter[0] * psat_targetT[0]) / P_current_Pa, // Corrected K-value
                (gamma_current_iter[1] * psat_targetT[1]) / P_current_Pa  // Corrected K-value
            ];
             if (K_iter.some(k => k <= 0 || isNaN(k))) {
                // Attempt to recover if P is way off, making K invalid
                if (j === 0 && P_current_Pa > 1E7) P_current_Pa *= 0.5; // P too high
                else if (j === 0 && P_current_Pa < 1E3) P_current_Pa *= 2.0; // P too low
                else return { comp1_feed: y_comp1, comp1_equilibrium: NaN, T_K: T_target_K, P_Pa: P_current_Pa, iterations: i, error: "Invalid K in DewP inner loop", calculationType: 'dewP' };
                continue; // Retry inner loop with adjusted P or fail
            }

            const x_iter_temp = [y[0] / K_iter[0], y[1] / K_iter[1]]; // Corrected x calculation
            const sum_x_iter_temp = x_iter_temp[0] + x_iter_temp[1];

             if (sum_x_iter_temp <= 1e-9 || isNaN(sum_x_iter_temp)) {
                 return { comp1_feed: y_comp1, comp1_equilibrium: NaN, T_K: T_target_K, P_Pa: P_current_Pa, iterations: i, error: "Invalid sum_x_iter in DewP inner loop", calculationType: 'dewP' };
            }

            const x_next_iter = [x_iter_temp[0] / sum_x_iter_temp, x_iter_temp[1] / sum_x_iter_temp]; // Corrected normalization
            const gamma_next = calculateUnifacGamma(componentData, x_next_iter, T_target_K, unifacParams);

            if (!gamma_next) {
                return { comp1_feed: y_comp1, comp1_equilibrium: NaN, T_K: T_target_K, P_Pa: P_current_Pa, iterations: i, error: "UNIFAC failed in DewP inner loop", calculationType: 'dewP' };
            }

            const diff_gamma = Math.abs(gamma_next[0] - gamma_current_iter[0]) + Math.abs(gamma_next[1] - gamma_current_iter[1]);
            gamma_current_iter = gamma_next;
            x_current = x_next_iter; // Update x_current based on this gamma

            if (diff_gamma < 1e-4) break; // Gamma (and x) converged for inner loop
        }
        // gamma_current_iter and x_current are now consistent for P_current_Pa

        const term1 = y[0] / (gamma_current_iter[0] * psat_targetT[0]);
        const term2 = y[1] / (gamma_current_iter[1] * psat_targetT[1]);
        if (isNaN(term1) || isNaN(term2) || (term1 + term2) <= 1e-12) { // Denominator close to zero or NaN
             return { comp1_feed: y_comp1, comp1_equilibrium: NaN, T_K: T_target_K, P_Pa: P_current_Pa, iterations: i, error: "DewP: Division by zero or NaN in P_calc", calculationType: 'dewP' };
        }
        const P_calc_Pa = 1.0 / (term1 + term2); // Corrected P_calc formula

        if (isNaN(P_calc_Pa)) {
            return { comp1_feed: y_comp1, comp1_equilibrium: NaN, T_K: T_target_K, P_Pa: P_current_Pa, iterations: i, error: "DewP: P_calc became NaN", calculationType: 'dewP' };
        }

        if (Math.abs(P_calc_Pa - P_current_Pa) < tolerance * P_current_Pa + tolerance) {
            // x_current is already the liquid composition consistent with P_calc_Pa (or very close P_current_Pa)
            return {
                comp1_feed: y_comp1,
                comp1_equilibrium: Math.min(Math.max(x_current[0], 0), 1), // x_comp1
                T_K: T_target_K,
                P_Pa: P_calc_Pa,
                iterations: i,
                calculationType: 'dewP'
            };
        }
        P_current_Pa = P_calc_Pa; // Successive substitution

        if (P_current_Pa < 0) {
             return { comp1_feed: y_comp1, comp1_equilibrium: NaN, T_K: T_target_K, P_Pa: P_current_Pa, iterations: i, error: "DewP: Pressure became negative", calculationType: 'dewP' };
        }
    }
    return { comp1_feed: y_comp1, comp1_equilibrium: NaN, T_K: T_target_K, P_Pa: P_current_Pa, iterations: max_iter, error: "DewP: Max iterations reached", calculationType: 'dewP' };
}

/** Fetch UNIFAC interaction parameters from Supabase */
import type { SupabaseClient } from '@supabase/supabase-js'; // Add this import

export async function fetchUnifacInteractionParams(
    supabase: SupabaseClient,
    allSubgroupIds: number[]
): Promise<UnifacParameters> {
    if (!supabase) { throw new Error("Supabase client not initialized."); }
    if (allSubgroupIds.length === 0) return { Rk: {}, Qk: {}, mainGroupMap: {}, a_mk: new Map() };

    console.log("Fetching UNIFAC parameters for subgroups (from library):", allSubgroupIds);

    try {
        // 1. Fetch Rk, Qk, and Main Group Mapping
        const { data: groupData, error: groupError } = await supabase
            .from('UNIFAC - Rk and Qk')
            .select('"Subgroup #", "Main Group #", "Rk", "Qk"')
            .in('"Subgroup #"', allSubgroupIds);

        if (groupError) throw new Error(`Supabase UNIFAC subgroup query error: ${groupError.message}`);
        if (!groupData || groupData.length < allSubgroupIds.length) {
            const foundIds = new Set(groupData?.map(g => g["Subgroup #"]) ?? []);
            const missingIds = allSubgroupIds.filter(id => !foundIds.has(id));
            throw new Error(`Missing UNIFAC parameters for subgroup ID(s): ${missingIds.join(', ')}`);
        }

        const Rk: { [id: number]: number } = {};
        const Qk: { [id: number]: number } = {};
        const mainGroupMap: { [id: number]: number } = {};
        const mainGroupIds = new Set<number>();

        groupData.forEach(g => {
            const subgroupId = g["Subgroup #"];
            const mainGroupId = g["Main Group #"];
            const rkVal = g["Rk"];
            const qkVal = g["Qk"];

            if (subgroupId == null || mainGroupId == null || rkVal == null || qkVal == null) {
                console.warn(`Incomplete data for subgroup ${subgroupId} in UNIFAC - Rk and Qk table. Skipping.`);
                // Or throw new Error(`Incomplete data for subgroup ${subgroupId}.`);
                return;
            }
            Rk[subgroupId] = rkVal;
            Qk[subgroupId] = qkVal;
            mainGroupMap[subgroupId] = mainGroupId;
            mainGroupIds.add(mainGroupId);
        });

        for (const reqId of allSubgroupIds) {
            if (mainGroupMap[reqId] === undefined) {
                throw new Error(`Failed to retrieve parameters or main group mapping for subgroup ID: ${reqId}. Check UNIFAC - Rk and Qk table.`);
            }
        }

        console.log("Fetched Rk, Qk, Main Group Map. Required Main Groups:", Array.from(mainGroupIds));

        // 2. Fetch Interaction Parameters (a_mk)
        const mainGroupArray = Array.from(mainGroupIds);
        if (mainGroupArray.length === 0) {
            console.warn("No main groups identified, cannot fetch interaction parameters.");
            return { Rk, Qk, mainGroupMap, a_mk: new Map() };
        }

        const { data: interactionData, error: interactionError } = await supabase
            .from('UNIFAC - a(ij)')
            .select('i, j, "a(ij)"')
            .in('i', mainGroupArray)
            .in('j', mainGroupArray);

        if (interactionError) throw new Error(`Supabase UNIFAC interaction query error: ${interactionError.message}`);
        if (!interactionData) throw new Error("Failed to query UNIFAC interaction parameters (interactionData is null/undefined).");

        const a_mk = new Map<string, number>();
        interactionData.forEach(interaction => {
            const main_group_m = interaction.i;
            const main_group_k = interaction.j;
            const a_mk_value = interaction["a(ij)"];

            if (main_group_m != null && main_group_k != null && a_mk_value != null) {
                a_mk.set(`${main_group_m}-${main_group_k}`, a_mk_value);
            } else {
                console.warn(`Incomplete interaction data found: i=${main_group_m}, j=${main_group_k} in UNIFAC - a(ij) table. Skipping.`);
            }
        });

        // --- Start of new logging ---
        console.log("DEBUG: Fetched a_mk map contents (all entries):");
        for (const [key, value] of a_mk.entries()) {
            console.log(`DEBUG: a_mk[${key}] = ${value}`);
        }

        // Specifically log the critical pairs if you know the main group IDs:
        // Replace with actual main group IDs for Methanol and Water from DWSIM unifac.txt / unifac_sg.txt
        // Example: CH3OH (Main Group 6), H2O (Main Group 7)
        const mg_meoh = 6; // Example Main Group ID for Methanol (e.g., from CH3OH group)
        const mg_h2o = 7;  // Example Main Group ID for Water (H2O group)
        
        console.log(`DEBUG: Specific a_mk for mg${mg_meoh}-mg${mg_h2o} (e.g. MeOH-H2O): ${a_mk.get(`${mg_meoh}-${mg_h2o}`)}`);
        console.log(`DEBUG: Specific a_mk for mg${mg_h2o}-mg${mg_meoh} (e.g. H2O-MeOH): ${a_mk.get(`${mg_h2o}-${mg_meoh}`)}`);
        // Log self-interaction parameters if needed for debugging
        console.log(`DEBUG: Specific a_mk for mg${mg_meoh}-mg${mg_meoh} (e.g. MeOH-MeOH): ${a_mk.get(`${mg_meoh}-${mg_meoh}`)}`);
        console.log(`DEBUG: Specific a_mk for mg${mg_h2o}-mg${mg_h2o} (e.g. H2O-H2O): ${a_mk.get(`${mg_h2o}-${mg_h2o}`)}`);
        // --- End of new logging ---

        for (const mg1 of mainGroupArray) {
            for (const mg2 of mainGroupArray) {
                const key = `${mg1}-${mg2}`;
                if (!a_mk.has(key)) {
                    if (mg1 === mg2) {
                        console.warn(`Missing diagonal interaction parameter for main group ${mg1} (a_mm). Assuming 0.`);
                        a_mk.set(key, 0.0);
                    } else {
                        console.warn(`Missing UNIFAC interaction parameter for main groups: ${key} (a_mk). Assuming 0. Check UNIFAC - a(ij) table.`);
                        a_mk.set(key, 0.0);
                        // Alternatively, could throw an error:
                        // throw new Error(`Missing UNIFAC interaction parameter for main groups: ${key}. Ensure all pairs are in UNIFAC - a(ij) table.`);
                    }
                }
            }
        }

        console.log("Fetched UNIFAC interaction parameters (a_mk) from library.");
        return { Rk, Qk, mainGroupMap, a_mk };
    } catch (err: any) {
        console.error("Error fetching UNIFAC parameters in library:", err.message);
        // Re-throw the error to be caught by the calling function in page.tsx
        throw err; 
    }
}