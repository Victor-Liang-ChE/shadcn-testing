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


/** Calculates UNIFAC activity coefficients gamma_i for a binary mixture */
export function calculateUnifacGamma(
    componentData: CompoundData[], // Uses imported CompoundData
    x: number[], // mole fractions [x1, x2]
    T_kelvin: number,
    params: UnifacParameters
): [number, number] | null { // Return null on error
    console.debug(`DEBUG: calculateUnifacGamma called for:`);
    console.debug(`DEBUG:   Component 1: ${componentData[0]?.name}, x1=${x[0]?.toFixed(4)}`);
    console.debug(`DEBUG:   Component 2: ${componentData[1]?.name}, x2=${x[1]?.toFixed(4)}`);
    console.debug(`DEBUG:   Temperature: ${T_kelvin.toFixed(2)} K`);
    // For Map, JSON.stringify won't work well directly, so we convert a_mk to an object for logging
    const a_mk_object: { [key: string]: number } = {};
    if (params?.a_mk) {
        for (const [key, value] of params.a_mk) {
            a_mk_object[key] = value;
        }
    }
    console.debug("DEBUG:   UNIFAC Rk params:", JSON.parse(JSON.stringify(params?.Rk || {})));
    console.debug("DEBUG:   UNIFAC Qk params:", JSON.parse(JSON.stringify(params?.Qk || {})));
    console.debug("DEBUG:   UNIFAC MainGroupMap:", JSON.parse(JSON.stringify(params?.mainGroupMap || {})));
    console.debug("DEBUG:   UNIFAC Interaction (a_mk) params:", JSON.parse(JSON.stringify(a_mk_object)));


    try {
        const numComponents = componentData.length;
        if (numComponents !== 2) throw new Error("UNIFAC calculation currently supports only binary mixtures.");
        if (!params || !params.Rk || !params.Qk || !params.mainGroupMap || !params.a_mk) {
            console.error("DEBUG: Incomplete UNIFAC parameters provided to calculateUnifacGamma:", { paramsExists: !!params, RkExists: !!params?.Rk, QkExists: !!params?.Qk, mainGroupMapExists: !!params?.mainGroupMap, a_mkExists: !!params?.a_mk });
            throw new Error("Incomplete UNIFAC parameters provided.");
        }

        const Z = 10.0; // Coordination number

        // --- Ensure r_i, q_i are calculated ---
        console.debug("DEBUG: Calculating r_i, q_i for components:");
        for (let i = 0; i < componentData.length; i++) {
            const comp = componentData[i];
            // console.debug(`DEBUG:   Processing component: ${comp.name}`);
            // console.debug(`DEBUG:     Initial r_i: ${comp.r_i}, q_i: ${comp.q_i}`);
            // console.debug(`DEBUG:     UNIFAC Groups for ${comp.name}:`, JSON.parse(JSON.stringify(comp.unifacGroups || {})));

            if (comp.r_i === undefined || comp.q_i === undefined) {
                // console.debug(`DEBUG:     r_i or q_i undefined for ${comp.name}, calculating now.`);
                comp.r_i = 0;
                comp.q_i = 0;
                if (!comp.unifacGroups) {
                    console.error(`DEBUG: Missing UNIFAC groups for ${comp.name} when r_i/q_i calculation needed.`);
                    throw new Error(`Missing UNIFAC groups for ${comp.name}`);
                }
                for (const subgroupIdStr in comp.unifacGroups) {
                    const subgroupId = parseInt(subgroupIdStr);
                    const count = comp.unifacGroups[subgroupId];
                    // console.debug(`DEBUG:     Subgroup ID: ${subgroupId}, Count: ${count}`);
                    if (params.Rk[subgroupId] === undefined || params.Qk[subgroupId] === undefined) {
                        console.error(`DEBUG: Missing Rk/Qk for subgroup ${subgroupId} (needed by ${comp.name}). Rk map: ${JSON.stringify(params.Rk)}, Qk map: ${JSON.stringify(params.Qk)}`);
                        throw new Error(`Missing Rk/Qk for subgroup ${subgroupId} needed by ${comp.name}`);
                    }
                    // console.debug(`DEBUG:       Rk[${subgroupId}] = ${params.Rk[subgroupId]}, Qk[${subgroupId}] = ${params.Qk[subgroupId]}`);
                    comp.r_i += count * params.Rk[subgroupId];
                    comp.q_i += count * params.Qk[subgroupId];
                }
            }
            console.debug(`DEBUG:   Final for ${comp.name}: r_i = ${comp.r_i?.toFixed(4)}, q_i = ${comp.q_i?.toFixed(4)}`);
        }

        // --- 1. Combinatorial Part ---
        const r = componentData.map(c => c.r_i!);
        const q = componentData.map(c => c.q_i!);
        const sum_xr = x[0] * r[0] + x[1] * r[1];
        const sum_xq = x[0] * q[0] + x[1] * q[1];
        // console.debug(`DEBUG: Combinatorial: r = [${r[0]?.toFixed(4)}, ${r[1]?.toFixed(4)}], q = [${q[0]?.toFixed(4)}, ${q[1]?.toFixed(4)}]`);
        // console.debug(`DEBUG: Combinatorial: sum_xr = ${sum_xr?.toFixed(4)}, sum_xq = ${sum_xq?.toFixed(4)}`);


        if (sum_xr === 0 || sum_xq === 0) {
            console.warn("Sum of x*r or x*q is zero. Returning ideal gammas.");
            return [1.0, 1.0];
        }

        const Phi = [x[0] * r[0] / sum_xr, x[1] * r[1] / sum_xr];
        const Theta = [x[0] * q[0] / sum_xq, x[1] * q[1] / sum_xq];
        const L = r.map((ri, i) => (Z / 2.0) * (ri - q[i]) - (ri - 1.0));
        const sum_xL = x[0] * L[0] + x[1] * L[1];
        // console.debug(`DEBUG: Combinatorial: Phi = [${Phi[0]?.toFixed(4)}, ${Phi[1]?.toFixed(4)}], Theta = [${Theta[0]?.toFixed(4)}, ${Theta[1]?.toFixed(4)}]`);
        // console.debug(`DEBUG: Combinatorial: L = [${L[0]?.toFixed(4)}, ${L[1]?.toFixed(4)}], sum_xL = ${sum_xL?.toFixed(4)}`);

        const ln_gamma_C: number[] = [0.0, 0.0];
        for (let i = 0; i < numComponents; i++) {
            const term1 = Math.log(Phi[i] / x[i]);
            const term2 = (Z / 2.0) * q[i] * Math.log(Theta[i] / Phi[i]);
            const term3 = L[i];
            const term4 = (Phi[i] / x[i]) * sum_xL;
            ln_gamma_C[i] = (isNaN(term1) ? 0 : term1) + (isNaN(term2) ? 0 : term2) + term3 - (isNaN(term4) ? 0 : term4);
        }
        console.debug(`DEBUG: ln_gamma_C = [${ln_gamma_C[0]?.toFixed(4)}, ${ln_gamma_C[1]?.toFixed(4)}]`);


        // --- 2. Residual Part ---
        const all_subgroup_ids_present = new Set<number>();
        componentData.forEach(comp => {
            if (comp.unifacGroups) {
                Object.keys(comp.unifacGroups).forEach(id => all_subgroup_ids_present.add(parseInt(id)));
            }
        });
        // console.debug("DEBUG: All subgroup IDs present in mixture:", Array.from(all_subgroup_ids_present));

        if (all_subgroup_ids_present.size === 0) {
            console.warn("No UNIFAC groups found for residual part. Returning combinatorial gammas.");
            return [Math.exp(ln_gamma_C[0]), Math.exp(ln_gamma_C[1])];
        }

        const subgroupIds = Array.from(all_subgroup_ids_present).sort((a, b) => a - b);
        const numGroups = subgroupIds.length;
        // console.debug("DEBUG: Sorted unique subgroup IDs for residual part:", subgroupIds);


        // Check if all required main groups and interactions are present
        const requiredMainGroups = new Set<number>();
        subgroupIds.forEach(sgId => {
            const mainGroupId = params.mainGroupMap[sgId];
            if (mainGroupId === undefined) {
                console.error(`DEBUG: Missing main group mapping for subgroup ${sgId}. MainGroupMap: ${JSON.stringify(params.mainGroupMap)}`);
                throw new Error(`Missing main group mapping for subgroup ${sgId}`);
            }
            requiredMainGroups.add(mainGroupId);
        });
        // console.debug("DEBUG: Required main groups:", Array.from(requiredMainGroups));

        for (const mg1 of requiredMainGroups) {
            for (const mg2 of requiredMainGroups) {
                const interactionKey = `${mg1}-${mg2}`;
                // console.debug(`DEBUG: Checking interaction parameter for main groups ${mg1}-${mg2} (key: ${interactionKey})`);
                if (!params.a_mk.has(interactionKey)) {
                    console.error(`DEBUG: Missing interaction parameter a_mk for main groups ${mg1}-${mg2} (key: ${interactionKey}). Available a_mk keys: ${Array.from(params.a_mk.keys())}`);
                    throw new Error(`Missing interaction parameter a_mk for main groups ${mg1}-${mg2}`);
                }
            }
        }


        const v_ki: number[][] = Array(numComponents).fill(0).map(() => Array(numGroups).fill(0));
        subgroupIds.forEach((sgId, k_idx) => {
            for (let i = 0; i < numComponents; i++) {
                v_ki[i][k_idx] = componentData[i].unifacGroups?.[sgId] || 0;
            }
        });
        // console.debug("DEBUG: Group counts matrix v_ki (comp x group_idx):", JSON.parse(JSON.stringify(v_ki)));


        let sum_x_vk_total = 0;
        for (let i = 0; i < numComponents; ++i) {
            for (let k_idx = 0; k_idx < numGroups; ++k_idx) {
                sum_x_vk_total += x[i] * v_ki[i][k_idx];
            }
        }

        const X_m: number[] = Array(numGroups).fill(0);
        if (sum_x_vk_total > 0) {
            for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                let sum_x_vm = 0;
                for (let i = 0; i < numComponents; i++) {
                    sum_x_vm += x[i] * v_ki[i][m_idx];
                }
                X_m[m_idx] = sum_x_vm / sum_x_vk_total;
            }
        }
        // console.debug("DEBUG: Group mole fractions X_m (mixture):", JSON.parse(JSON.stringify(X_m.map(val => val.toFixed(4)))));

        let sum_XQ = 0;
        for (let n_idx = 0; n_idx < numGroups; ++n_idx) {
            sum_XQ += X_m[n_idx] * params.Qk[subgroupIds[n_idx]];
        }

        const theta_m: number[] = Array(numGroups).fill(0);
        if (sum_XQ > 0) {
            for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                theta_m[m_idx] = (X_m[m_idx] * params.Qk[subgroupIds[m_idx]]) / sum_XQ;
            }
        }
        // console.debug("DEBUG: Group surface area fractions theta_m (mixture):", JSON.parse(JSON.stringify(theta_m.map(val => val.toFixed(4)))));


        const psi_matrix: number[][] = Array(numGroups).fill(0).map(() => Array(numGroups).fill(0));
        console.debug("DEBUG: Calculating Psi_nm matrix:");
        for (let n_idx = 0; n_idx < numGroups; n_idx++) {
            for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                const sg_n = subgroupIds[n_idx];
                const sg_m = subgroupIds[m_idx];
                const mainG_n = params.mainGroupMap[sg_n];
                const mainG_m = params.mainGroupMap[sg_m];
                const interactionKey = `${mainG_n}-${mainG_m}`;
                const a_nm = params.a_mk.get(interactionKey) ?? 0.0;
                psi_matrix[n_idx][m_idx] = Math.exp(-a_nm / T_kelvin);
                // console.debug(`DEBUG:   Psi[sg${sg_n}(mg${mainG_n}) - sg${sg_m}(mg${mainG_m})]: a_mk(${interactionKey}) = ${a_nm.toFixed(2)}, Psi_nm = ${psi_matrix[n_idx][m_idx].toFixed(4)} at T=${T_kelvin.toFixed(2)}K`);
            }
        }
        // console.debug("DEBUG: Psi_matrix (sg_idx x sg_idx):", JSON.parse(JSON.stringify(psi_matrix.map(row => row.map(val => val.toFixed(4))))));


        const ln_Gamma_k: number[] = Array(numGroups).fill(0);
        for (let k_idx = 0; k_idx < numGroups; k_idx++) {
            let sum_theta_psi_mk = 0;
            for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                sum_theta_psi_mk += theta_m[m_idx] * psi_matrix[m_idx][k_idx];
            }
            const term1 = Math.log(sum_theta_psi_mk);

            let term2_sum = 0;
            for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                let sum_theta_psi_nm = 0;
                for (let n_idx = 0; n_idx < numGroups; n_idx++) {
                    sum_theta_psi_nm += theta_m[n_idx] * psi_matrix[n_idx][m_idx];
                }
                if (sum_theta_psi_nm > 0) {
                    term2_sum += (theta_m[m_idx] * psi_matrix[k_idx][m_idx]) / sum_theta_psi_nm;
                }
            }
            ln_Gamma_k[k_idx] = params.Qk[subgroupIds[k_idx]] * (1.0 - (isNaN(term1) ? 0 : term1) - term2_sum);
        }
        // console.debug("DEBUG: ln_Gamma_k (mixture groups):", JSON.parse(JSON.stringify(ln_Gamma_k.map((val, idx) => `sg${subgroupIds[idx]}: ${val.toFixed(4)}`))));


        const ln_Gamma_k_pure: number[][] = Array(numComponents).fill(0).map(() => Array(numGroups).fill(0));
        for (let i = 0; i < numComponents; i++) {
            // console.debug(`DEBUG: Calculating ln_Gamma_k_pure for component ${componentData[i].name}`);
            let sum_vk_pure = 0;
            for (let k_idx = 0; k_idx < numGroups; ++k_idx) sum_vk_pure += v_ki[i][k_idx];

            const X_m_pure: number[] = Array(numGroups).fill(0);
            if (sum_vk_pure > 0) {
                for (let m_idx = 0; m_idx < numGroups; ++m_idx) X_m_pure[m_idx] = v_ki[i][m_idx] / sum_vk_pure;
            }
            // console.debug(`DEBUG:   X_m_pure for ${componentData[i].name}:`, JSON.parse(JSON.stringify(X_m_pure.map(val => val.toFixed(4)))));


            let sum_XQ_pure = 0;
            for (let n_idx = 0; n_idx < numGroups; ++n_idx) sum_XQ_pure += X_m_pure[n_idx] * params.Qk[subgroupIds[n_idx]];

            const theta_m_pure: number[] = Array(numGroups).fill(0);
            if (sum_XQ_pure > 0) {
                for (let m_idx = 0; m_idx < numGroups; ++m_idx) {
                    theta_m_pure[m_idx] = (X_m_pure[m_idx] * params.Qk[subgroupIds[m_idx]]) / sum_XQ_pure;
                }
            }
            // console.debug(`DEBUG:   theta_m_pure for ${componentData[i].name}:`, JSON.parse(JSON.stringify(theta_m_pure.map(val => val.toFixed(4)))));


            for (let k_idx = 0; k_idx < numGroups; k_idx++) {
                let sum_theta_psi_mk_pure = 0;
                for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                    sum_theta_psi_mk_pure += theta_m_pure[m_idx] * psi_matrix[m_idx][k_idx];
                }
                const term1_pure = Math.log(sum_theta_psi_mk_pure);

                let term2_sum_pure = 0;
                for (let m_idx = 0; m_idx < numGroups; m_idx++) {
                    let sum_theta_psi_nm_pure = 0;
                    for (let n_idx = 0; n_idx < numGroups; n_idx++) {
                        sum_theta_psi_nm_pure += theta_m_pure[n_idx] * psi_matrix[n_idx][m_idx];
                    }
                    if (sum_theta_psi_nm_pure > 0) {
                        term2_sum_pure += (theta_m_pure[m_idx] * psi_matrix[k_idx][m_idx]) / sum_theta_psi_nm_pure;
                    }
                }
                ln_Gamma_k_pure[i][k_idx] = params.Qk[subgroupIds[k_idx]] * (1.0 - (isNaN(term1_pure) ? 0 : term1_pure) - term2_sum_pure);
            }
            // console.debug(`DEBUG:   ln_Gamma_k_pure for ${componentData[i].name}:`, JSON.parse(JSON.stringify(ln_Gamma_k_pure[i].map((val, idx) => `sg${subgroupIds[idx]}: ${val.toFixed(4)}` ))));
        }

        const ln_gamma_R: number[] = [0.0, 0.0];
        for (let i = 0; i < numComponents; i++) {
            for (let k_idx = 0; k_idx < numGroups; k_idx++) {
                if (v_ki[i][k_idx] > 0) {
                    ln_gamma_R[i] += v_ki[i][k_idx] * (ln_Gamma_k[k_idx] - ln_Gamma_k_pure[i][k_idx]);
                }
            }
        }
        console.debug(`DEBUG: ln_gamma_R = [${ln_gamma_R[0]?.toFixed(4)}, ${ln_gamma_R[1]?.toFixed(4)}]`);

        const gamma: [number, number] = [
            Math.exp(ln_gamma_C[0] + ln_gamma_R[0]),
            Math.exp(ln_gamma_C[1] + ln_gamma_R[1])
        ];
        console.debug(`DEBUG: Final gamma = [${gamma[0]?.toFixed(4)}, ${gamma[1]?.toFixed(4)}]`);


        if (isNaN(gamma[0]) || isNaN(gamma[1]) || !isFinite(gamma[0]) || !isFinite(gamma[1])) {
             console.error("UNIFAC calculation resulted in NaN or Infinity gamma values.", {x, T_kelvin, ln_gamma_C, ln_gamma_R});
             throw new Error("UNIFAC calculation failed (NaN/Infinity gamma).");
        }
        return gamma;

    } catch (error) {
        console.error("Error during UNIFAC calculation:", error);
        return null; // Indicate failure
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