// --- Shared Data Structures for VLE Calculations ---

export interface AntoineParams {
    A: number;
    B: number;
    C: number;
    Tmin_K: number;
    Tmax_K: number;
    Units: 'Pa' | 'kPa' | string; // Allow string for flexibility
    EquationNo?: number | string;
}

export interface PrPureComponentParams {
    Tc_K: number;    // Critical Temperature (K)
    Pc_Pa: number;   // Critical Pressure (Pa)
    omega: number;   // Acentric factor
}

export interface SrkPureComponentParams { // ADDED for SRK
    Tc_K: number;    // Critical Temperature (K)
    Pc_Pa: number;   // Critical Pressure (Pa)
    omega: number;   // Acentric factor
}

export interface UniquacPureComponentParams { // ADDED for UNIQUAC
    r: number; // Van der Waals volume parameter
    q: number; // Van der Waals surface area parameter
}

export interface WilsonPureComponentParams { // ADDED for Wilson
    V_L_m3mol: number; // Liquid molar volume in m^3/mol
}

export interface UnifacGroupComposition {
    [subgroupId: number]: number; // e.g. {1: 2, 2: 1} for CH3-CH2-OH (2x CH3, 1x CH2, 1x OH)
}

export interface CompoundData {
    name: string;
    cas_number?: string | null; // CAS number is crucial for NRTL/PR/SRK interaction params
    antoine: AntoineParams | null;
    unifacGroups?: UnifacGroupComposition | null; // Optional: only needed for UNIFAC
    prParams?: PrPureComponentParams | null;       // Optional: only needed for Peng-Robinson
    srkParams?: SrkPureComponentParams | null;      // ADDED: Optional: only needed for SRK
    uniquacParams?: UniquacPureComponentParams | null; // ADDED for UNIQUAC
    wilsonParams?: WilsonPureComponentParams | null; // ADDED for Wilson
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
