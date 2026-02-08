// --- Shared Data Structures for VLE Calculations ---
// Updated for HYSYS Supabase tables (7-param extended Antoine, name-based lookups)

/**
 * 7-parameter extended Antoine equation (DIPPR 101 / HYSYS form):
 *   ln(P [kPa]) = A + B/(T+C) + D·ln(T) + E·T^F
 * where T is in Kelvin and P is in kPa.
 * Multiply result by 1000 to get Pa.
 * Note: G is stored but unused in the HYSYS equation form.
 */
export interface AntoineParams {
    A: number;
    B: number;
    C: number;
    D: number;
    E: number;
    F: number;
    G: number;
    Tmin_K: number;
    Tmax_K: number;
}

export interface PrPureComponentParams {
    Tc_K: number;    // Critical Temperature (K)
    Pc_Pa: number;   // Critical Pressure (Pa)
    omega: number;   // Acentric factor
}

export interface SrkPureComponentParams {
    Tc_K: number;    // Critical Temperature (K)
    Pc_Pa: number;   // Critical Pressure (Pa)
    omega: number;   // Acentric factor
}

export interface UniquacPureComponentParams {
    r: number; // UNIQUAC r parameter (volume)
    q: number; // UNIQUAC q parameter (surface area)
}

export interface WilsonPureComponentParams {
    V_L_m3mol: number; // Liquid molar volume in m^3/mol (converted from m³/kmol ÷ 1000)
}

export interface UnifacGroupComposition {
    [subgroupId: number]: number; // e.g. {1: 2, 2: 1} for CH3-CH2-OH (2x CH3, 1x CH2, 1x OH)
}

/**
 * HYSYS NRTL interaction parameters.
 * τ_ij = Aij/T + Bij, G_ij = exp(-Cij_Alpha · τ_ij)
 * Aij in K, Bij dimensionless, Cij_Alpha dimensionless.
 */
export interface HysysNrtlParams {
    Aij: number;  // K
    Aji: number;  // K
    Bij: number;  // dimensionless
    Bji: number;  // dimensionless
    Cij_Alpha: number; // non-randomness parameter (alpha)
    Cji_Alpha: number; // usually same as Cij_Alpha
}

/**
 * HYSYS Wilson interaction parameters.
 * Λ_12 = (V2/V1) · exp(-(A12/T + B12))
 * Aij in K, Bij dimensionless.
 */
export interface HysysWilsonParams {
    Aij: number; // K
    Aji: number; // K
    Bij: number; // dimensionless
    Bji: number; // dimensionless
}

/**
 * HYSYS UNIQUAC interaction parameters.
 * τ_ij = exp(-(Aij/T + Bij))
 * Aij in K, Bij dimensionless.
 */
export interface HysysUniquacParams {
    Aij: number; // K
    Aji: number; // K
    Bij: number; // dimensionless
    Bji: number; // dimensionless
}

/**
 * HYSYS PR SRK shared Kij table.
 * Binary interaction parameter, dimensionless.
 */
export interface HysysPrSrkParams {
    Kij: number; // dimensionless
}

export interface CompoundData {
    name: string;
    cas_number?: string | null;
    molecularWeight?: number | null;
    antoine: AntoineParams | null;
    unifacGroups?: UnifacGroupComposition | null;
    prParams?: PrPureComponentParams | null;
    srkParams?: SrkPureComponentParams | null;
    uniquacParams?: UniquacPureComponentParams | null;
    wilsonParams?: WilsonPureComponentParams | null;
    // For UNIFAC pre-calculation, can be added dynamically
    r_i?: number; // Sum of Rk for the molecule
    q_i?: number; // Sum of Qk for the molecule
}

export interface BubbleDewResult {
    comp1_feed: number;        // Input mole fraction of component 1 (either x or y)
    comp1_equilibrium: number; // Calculated mole fraction of component 1 in the other phase
    T_K?: number;              // Calculated Temperature (for BubbleT, DewT)
    P_Pa?: number;             // Calculated Pressure (for BubbleP, DewP)
    iterations?: number;
    error?: string;            // Error message if calculation failed
    calculationType?: 'bubbleT' | 'bubbleP' | 'dewT' | 'dewP';
}
