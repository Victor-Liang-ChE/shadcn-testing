// --- Shared Data Structures for VLE Calculations ---
// Updated for Aspen Plus (NIST140 / APV140) Supabase tables.

/**
 * 7-parameter extended Antoine equation (Aspen PLXANT / Riedel form):
 *   ln(P* [Pa]) = C1 + C2/(T+C3) + C4·T + C5·ln(T) + C6·T^C7
 * where T is in Kelvin and P is in Pa.
 * Columns map directly to apv140_plxant_wide CSV columns.
 */
export interface AntoineParams {
    C1: number;
    C2: number;
    C3: number;  // shift in denominator (often 0)
    C4: number;  // linear T coefficient
    C5: number;  // ln(T) coefficient
    C6: number;  // power-law coefficient (usually very small)
    C7: number;  // power-law exponent
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
    V_L_m3mol: number; // Liquid molar volume at 25°C in m³/mol (Spencer-Danner/Rackett reference)
    // Optional critical properties for temperature-dependent Rackett VL calculation
    Tc_K?: number;
    Pc_Pa?: number;
    omega?: number;
    RKTZRA?: number; // Aspen Rackett compressibility factor (RKTZRA column); if absent Yamada-Gunn is used
}

export interface UnifacGroupComposition {
    [subgroupId: number]: number; // e.g. {1: 2, 2: 1} for CH3-CH2-OH (2x CH3, 1x CH2, 1x OH)
}

/**
 * Aspen NRTL interaction parameters (stored as R_cal-scaled values for
 * compatibility with the existing calculateNrtlGamma function).
 * Raw Aspen DB cols: aij (dim'less), bij (K), cij (alpha).
 * Stored here as: Aij = bij * R_cal, Bij = aij * R_cal so that
 *   τ_ij = Aij/(R_cal·T) + Bij/R_cal = bij/T + aij  (standard Aspen form).
 */
export interface HysysNrtlParams {
    Aij: number;
    Aji: number;
    Bij: number;
    Bji: number;
    Cij_Alpha: number;
    Cji_Alpha: number;
}

/**
 * Aspen Wilson interaction parameters (R_cal-scaled).
 * Raw Aspen: aij (dim'less), bij (K).
 * Stored as: Aij = -bij * R_cal, Bij = -aij * R_cal so that
 *   Λ_12 = (V2/V1)·exp(-(Aij/(R_cal·T) + Bij/R_cal)) = (V2/V1)·exp(aij + bij/T).
 */
export interface HysysWilsonParams {
    Aij: number;
    Aji: number;
    Bij: number;
    Bji: number;
}

/**
 * Aspen UNIQUAC interaction parameters (R_cal-scaled).
 * Raw Aspen: aij (dim'less), bij (K).
 * Stored as: Aij = -bij * R_cal, Bij = -aij * R_cal so that
 *   τ_ij = exp(-(Aij/(R_cal·T) + Bij/R_cal)) = exp(aij + bij/T).
 */
export interface HysysUniquacParams {
    Aij: number;
    Aji: number;
    Bij: number;
    Bji: number;
}

/**
 * Aspen PR/SRK binary interaction parameter.
 * kij1 from prkbv_wide / rkskbv_wide tables (temperature-independent Kij).
 */
export interface HysysPrSrkParams {
    Kij: number;
}

/**
 * Hayden-O'Connell (1975) pure-component parameters for vapor-phase association.
 * Used to compute second virial coefficients for fugacity corrections.
 *   μ [Debye]  = MUP_stored  / 3.33564e-25
 *   r_d [Å]    = RGYR_m × 1e10
 * eta_self: self-association parameter (0 = non-associating, >0 = H-bonding/polar)
 */
export interface HocPureComponentParams {
    MUP_stored: number;  // raw DB value; convert to Debye: / 3.33564e-25
    RGYR_m: number;      // raw DB value in m; convert to Å: × 1e10
    eta_self: number;    // self-association parameter (HOC_ETA_SELF, null→0 for non-associating)
    // Convenience values pre-computed at fetch time:
    mu_D: number;        // dipole moment in Debye
    rd_A: number;        // radius of gyration in Angstroms
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
    hocProps?: HocPureComponentParams | null;   // HOC vapor-phase association params
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
