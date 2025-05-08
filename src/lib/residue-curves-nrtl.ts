import type { SupabaseClient } from '@supabase/supabase-js';
import type { CompoundData, AntoineParams } from './vle-types'; // Assuming vle-types.ts is in the same directory
import { fetchNrtlParameters, NrtlInteractionParams } from './vle-calculations-nrtl'; // Assuming vle-calculations-nrtl.ts is in the same directory
import { calculatePsat_Pa } from './vle-calculations-unifac'; // Assuming vle-calculations-unifac.ts is in the same directory

// Gas constant in cal/mol·K, consistent with vle-calculations-nrtl.ts
const R_cal_molK = 1.98720425864083;

export interface TernaryNrtlParameters {
    params12: NrtlInteractionParams;
    params13: NrtlInteractionParams;
    params23: NrtlInteractionParams;
    // Store CAS numbers for clarity, matching component order 1, 2, 3
    casn1: string;
    casn2: string;
    casn3: string;
}

export interface ResidueCurvePoint {
    x1: number; // Mole fraction of component 1 in liquid
    x2: number; // Mole fraction of component 2 in liquid
    x3: number; // Mole fraction of component 3 in liquid
    y1?: number; // Mole fraction of component 1 in vapor
    y2?: number; // Mole fraction of component 2 in vapor
    y3?: number; // Mole fraction of component 3 in vapor
    T_K: number;  // Temperature (assumed constant for the curve)
    P_Pa?: number; // Total pressure at this point
    xi?: number; // Dimensionless extent of evaporation (optional)
}

export type ResidueCurve = ResidueCurvePoint[];

/**
 * Fetches all necessary binary NRTL parameters for a ternary system.
 */
export async function fetchTernaryNrtlParameters(
    supabase: SupabaseClient,
    casn1: string,
    casn2: string,
    casn3: string
): Promise<TernaryNrtlParameters> {
    if (!casn1 || !casn2 || !casn3) {
        throw new Error("CAS numbers for all three compounds are required.");
    }
    if (new Set([casn1, casn2, casn3]).size !== 3) {
        throw new Error("All three CAS numbers must be unique.");
    }

    const params12 = await fetchNrtlParameters(supabase, casn1, casn2);
    const params13 = await fetchNrtlParameters(supabase, casn1, casn3);
    const params23 = await fetchNrtlParameters(supabase, casn2, casn3);

    return { params12, params13, params23, casn1, casn2, casn3 };
}

/**
 * Calculates NRTL activity coefficients (gamma_i) for a ternary mixture.
 * @param x Array of liquid mole fractions [x1, x2, x3].
 * @param T_K Temperature in Kelvin.
 * @param nrtlParams TernaryNrtlParameters object containing all binary interaction params.
 * @returns A tuple [gamma1, gamma2, gamma3] or null if calculation fails.
 */
export function calculateTernaryNrtlGamma(
    x: number[],
    T_K: number,
    nrtlParams: TernaryNrtlParameters
): [number, number, number] | null {
    if (x.length !== 3) {
        console.error("Mole fraction array must contain 3 components.");
        return null;
    }
    if (Math.abs(x[0] + x[1] + x[2] - 1.0) > 1e-6) {
        console.error("Mole fractions do not sum to 1:", x);
        return null;
    }

    const { params12, params13, params23 } = nrtlParams;

    // Create matrices for tau and G for easier indexing
    // tau_ij = A_ij / (R*T), G_ij = exp(-alpha_ij * tau_ij)
    // A_ii = 0, alpha_ii = (any, typically not used as G_ii = 1)

    const A = [
        [0, params12.A12, params13.A12],
        [params12.A21, 0, params23.A12], // A23 from params23.A12 (comp2 as first in pair 2-3)
        [params13.A21, params23.A21, 0]  // A32 from params23.A21 (comp3 as second in pair 2-3)
    ];

    const alpha = [
        [0, params12.alpha, params13.alpha],
        [params12.alpha, 0, params23.alpha],
        [params13.alpha, params23.alpha, 0]
    ];

    const tau: number[][] = [[0,0,0],[0,0,0],[0,0,0]];
    const G: number[][] = [[1,0,0],[0,1,0],[0,0,1]]; // G_ii = 1

    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            if (i === j) continue;
            tau[i][j] = A[i][j] / (R_cal_molK * T_K);
            G[i][j] = Math.exp(-alpha[i][j] * tau[i][j]);
        }
    }

    const lnGamma: number[] = [0, 0, 0];

    for (let i = 0; i < 3; i++) {
        let term1_num = 0;
        let term1_den = 0;
        for (let j = 0; j < 3; j++) {
            term1_num += x[j] * tau[j][i] * G[j][i];
            term1_den += x[j] * G[j][i];
        }
        const term1 = term1_den === 0 ? 0 : term1_num / term1_den;

        let term2_sum = 0;
        for (let j = 0; j < 3; j++) {
            if (i === j) continue; // Original NRTL formulation sums over j != i for this part, but modern forms sum over all j and rely on G_ii=1, tau_ii=0
            
            let term2_sub_num = 0;
            let term2_sub_den = 0;
            for (let m = 0; m < 3; m++) {
                term2_sub_num += x[m] * tau[m][j] * G[m][j];
                term2_sub_den += x[m] * G[m][j];
            }
            const term2_sub = term2_sub_den === 0 ? 0 : term2_sub_num / term2_sub_den;
            
            let G_ij_x_k_sum_den = 0;
             for (let k = 0; k < 3; k++) {
                G_ij_x_k_sum_den += x[k] * G[k][j];
            }
            if (G_ij_x_k_sum_den === 0) continue;

            term2_sum += (x[j] * G[i][j] / G_ij_x_k_sum_den) * (tau[i][j] - term2_sub);
        }
        lnGamma[i] = term1 + term2_sum;
    }
    
    const gammas: [number, number, number] = [Math.exp(lnGamma[0]), Math.exp(lnGamma[1]), Math.exp(lnGamma[2])];

    if (gammas.some(isNaN) || gammas.some(g => !isFinite(g))) {
        console.warn(`Ternary NRTL gamma calculation resulted in NaN/Infinity. x=[${x.map(xi=>xi.toFixed(3)).join(',')}] T=${T_K.toFixed(1)}K`);
        return null;
    }
    return gammas;
}

/**
 * Defines the system of differential equations for ternary residue curves.
 * dx_i/d(xi_evap) = x_i - y_i
 * We only need to solve for x1 and x2, as x3 = 1 - x1 - x2.
 * @param x_liq_pair Array of [x1, x2] liquid mole fractions.
 * @param _xi_evap Dummy dimensionless evaporation extent (equations are autonomous).
 * @param components Array of three CompoundData objects.
 * @param T_K Constant temperature in Kelvin.
 * @param nrtlParams TernaryNrtlParameters object.
 * @returns [dx1/d(xi_evap), dx2/d(xi_evap)] or null if error.
 */
export function residueCurveDifferentials(
    x_liq_pair: [number, number],
    _xi_evap: number, // Placeholder for ODE solver compatibility
    components: CompoundData[],
    T_K: number,
    nrtlParams: TernaryNrtlParameters
): [number, number] | null {
    const x1 = x_liq_pair[0];
    const x2 = x_liq_pair[1];
    const x3 = 1 - x1 - x2;

    if (x1 < -1e-6 || x2 < -1e-6 || x3 < -1e-6 || x1 > 1 + 1e-6 || x2 > 1 + 1e-6 || x3 > 1 + 1e-6) {
        console.warn("Mole fraction out of bounds in differentials:", [x1, x2, x3]);
        return null; // Stop integration or handle error
    }
    // Clamp to valid range to prevent small numerical errors from propagating
    const x_liq: [number, number, number] = [
        Math.max(0, Math.min(1, x1)),
        Math.max(0, Math.min(1, x2)),
        Math.max(0, Math.min(1, 1 - Math.max(0, Math.min(1, x1)) - Math.max(0, Math.min(1, x2))))
    ];
    
    if (components.length !== 3) return null;
    if (!components.every(c => c.antoine)) return null;

    const Psat: number[] = components.map(c => calculatePsat_Pa(c.antoine!, T_K, c.name));
    if (Psat.some(isNaN)) return null;

    const gammas = calculateTernaryNrtlGamma(x_liq, T_K, nrtlParams);
    if (!gammas) return null;

    const P_total_terms = x_liq.map((x_i, i) => gammas[i] * x_i * Psat[i]);
    const P_total = P_total_terms.reduce((sum, term) => sum + term, 0);

    if (P_total <= 1e-9) { // Avoid division by zero if pressure is effectively zero
        console.warn("Total pressure near zero in residueCurveDifferentials.");
        return [0,0]; // Or handle as an error/stop condition
    }

    const y_vap: number[] = P_total_terms.map(term => term / P_total);

    const dx_d_xi: [number, number] = [
        x_liq[0] - y_vap[0],
        x_liq[1] - y_vap[1]
    ];
    // dx3_d_xi would be x_liq[2] - y_vap[2]

    return dx_d_xi;
}

/**
 * Placeholder for tracing a single residue curve using numerical integration (e.g., Euler or RK4).
 * @param x_initial_pair Initial liquid mole fractions [x1, x2].
 * @param components Array of three CompoundData objects.
 * @param T_K Constant temperature in Kelvin.
 * @param nrtlParams TernaryNrtlParameters object.
 * @param stepSize_xi Step size for the dimensionless evaporation extent.
 * @param numSteps Maximum number of integration steps.
 * @returns A ResidueCurve (array of ResidueCurvePoint).
 */
export function traceResidueCurve(
    x_initial_pair: [number, number],
    components: CompoundData[],
    T_K: number,
    nrtlParams: TernaryNrtlParameters,
    stepSize_xi: number = 0.01,
    numSteps: number = 1000
): ResidueCurve {
    const curve: ResidueCurve = [];
    let current_x_pair: [number, number] = [...x_initial_pair];
    let current_xi = 0;

    for (let step = 0; step < numSteps; step++) {
        const x1 = current_x_pair[0];
        const x2 = current_x_pair[1];
        const x3 = 1 - x1 - x2;

        // Store current point
        // For full point data, recalculate gammas, Psats, P_total, y_vap here
        const Psats = components.map(c => calculatePsat_Pa(c.antoine!, T_K, c.name));
        const gammas = calculateTernaryNrtlGamma([x1,x2,x3], T_K, nrtlParams);
        let P_total_val: number | undefined = undefined;
        let y_vals: [number, number, number] | undefined = undefined;

        if (gammas && !Psats.some(isNaN)) {
            const P_terms = [x1,x2,x3].map((xi, idx) => gammas[idx] * xi * Psats[idx]);
            P_total_val = P_terms.reduce((s, t) => s + t, 0);
            if (P_total_val > 1e-9) {
                 y_vals = P_terms.map(pt => pt / P_total_val!) as [number,number,number];
            }
        }
        
        curve.push({ x1, x2, x3, y1: y_vals?.[0], y2: y_vals?.[1], y3: y_vals?.[2], T_K, P_Pa: P_total_val, xi: current_xi });

        // Calculate differentials
        const differentials = residueCurveDifferentials(current_x_pair, current_xi, components, T_K, nrtlParams);
        if (!differentials) {
            console.warn("Differential calculation failed. Stopping trace.");
            break;
        }

        // Euler step
        const next_x1 = x1 + differentials[0] * stepSize_xi;
        const next_x2 = x2 + differentials[1] * stepSize_xi;
        
        current_xi += stepSize_xi;
        current_x_pair = [next_x1, next_x2];

        // Stopping conditions
        if (next_x1 < 0 || next_x2 < 0 || (next_x1 + next_x2) > 1.0) {
             // Went out of bounds, try to project last valid point or just stop
            console.log("Residue curve went out of bounds. Stopping.");
            break;
        }
        if (Math.abs(differentials[0]) < 1e-5 && Math.abs(differentials[1]) < 1e-5) {
            console.log("Residue curve reached a stationary point (azeotrope/pure component). Stopping.");
            // Add the stationary point
            const final_x3 = 1 - next_x1 - next_x2;
            curve.push({ x1: next_x1, x2: next_x2, x3: final_x3, T_K, xi: current_xi });
            break;
        }
    }
    return curve;
}

// Example usage (conceptual, would be called from a UI or test script)
/*
async function example(supabase: SupabaseClient) {
    const casNumbers = {
        comp1: "67-64-1",  // Acetone
        comp2: "67-56-1",  // Methanol
        comp3: "7732-18-5" // Water (example, ensure Antoine & NRTL params exist)
    };
    const antoineComp1: AntoineParams = { A: 10.0, B: 1000, C: -50 }; // Placeholder
    const antoineComp2: AntoineParams = { A: 10.1, B: 1100, C: -40 }; // Placeholder
    const antoineComp3: AntoineParams = { A: 10.2, B: 1200, C: -30 }; // Placeholder

    const components: CompoundData[] = [
        { name: "Acetone", casn: casNumbers.comp1, antoine: antoineComp1, unifacGroups: [] },
        { name: "Methanol", casn: casNumbers.comp2, antoine: antoineComp2, unifacGroups: [] },
        { name: "Water", casn: casNumbers.comp3, antoine: antoineComp3, unifacGroups: [] }
    ];

    try {
        const ternaryNrtlParams = await fetchTernaryNrtlParameters(supabase, casNumbers.comp1, casNumbers.comp2, casNumbers.comp3);
        const T_System_K = 333.15; // Example temperature

        const initial_x_pair: [number, number] = [0.1, 0.1]; // Starting x1, x2 (x3 = 0.8)
        const curve = traceResidueCurve(initial_x_pair, components, T_System_K, ternaryNrtlParams, 0.01, 500);
        console.log("Generated Residue Curve:", curve);
    } catch (error) {
        console.error("Error in residue curve example:", error);
    }
}
*/

