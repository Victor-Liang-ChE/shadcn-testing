import type { SupabaseClient } from '@supabase/supabase-js';
import type { AntoineParams, UnifacGroupComposition, UniquacPureComponentParams, CompoundData } from './vle-types';

// ─── Compound name helpers (compounds_master table) ────────────────────────────

/** An entry from the compounds_master table. */
export interface CompoundAlias {
    fullName: string;
    simName: string;
}

/**
 * Format compound names for display using smart title-case rules that mirror
 * the Python smart_title() script used to pre-process the Aspen data.
 *
 * Rules applied (in order):
 *  1. Replace underscores with spaces.
 *  2. Convert to title-case (capitalise first letter of every "word",
 *     where any non-letter character acts as a word boundary – same
 *     behaviour as Python's str.title()).
 *  3. Lower-case standard chemical iso/n/sec/tert/cis/trans/o/m/p prefixes.
 *  4. Lower-case common prepositions and conjunctions (of, and, in, with).
 */
export function formatCompoundName(name: string): string {
    if (!name) return name;

    // 1. Underscores → spaces
    let t = name.replace(/_/g, ' ');

    // 2. Title-case: lower-case everything, then capitalise after any non-letter
    t = t.toLowerCase().replace(/\b\w/g, c => c.toUpperCase());

    // 3. Fix chemical prefixes that should remain lower-case
    const prefixPairs: [RegExp, string][] = [
        [/^N-/,      'n-'],    [/ N-/g,      ' n-'],
        [/^Iso-/,    'iso-'],  [/ Iso-/g,    ' iso-'],
        [/^Sec-/,    'sec-'],  [/ Sec-/g,    ' sec-'],
        [/^Tert-/,   'tert-'], [/ Tert-/g,   ' tert-'],
        [/^Cis-/,    'cis-'],  [/ Cis-/g,    ' cis-'],
        [/^Trans-/,  'trans-'],[/ Trans-/g,  ' trans-'],
        [/^O-/,      'o-'],    [/ O-/g,      ' o-'],
        [/^M-/,      'm-'],    [/ M-/g,      ' m-'],
        [/^P-/,      'p-'],    [/ P-/g,      ' p-'],
    ];
    for (const [pattern, replacement] of prefixPairs) {
        t = t.replace(pattern, replacement);
    }

    // 4. Lower-case common prepositions / conjunctions mid-name
    for (const word of [' Of ', ' And ', ' In ', ' With ']) {
        t = t.split(word).join(word.toLowerCase());
    }

    return t;
}

/**
 * Search compounds_master for compounds whose Name matches the user-typed prefix.
 * Deduplicates on Name since compounds_master stores one row per (Name, Alias) pair.
 * Returns {fullName, simName} pairs where both equal the canonical Aspen Name.
 */
export async function fetchCompoundSuggestions(
    supabase: SupabaseClient,
    query: string,
    limit = 8,
): Promise<CompoundAlias[]> {
    if (!query || query.trim().length < 1) return [];

    // Fetch extra rows to account for Name duplicates (one per alias)
    const { data, error } = await supabase
        .from('compounds_master')
        .select('"Name"')
        .ilike('Name', `${query.trim()}%`)
        .limit(limit * 5);

    if (error || !data) {
        console.error('fetchCompoundSuggestions error:', error);
        return [];
    }

        const seen = new Set<string>();
        const unique: CompoundAlias[] = [];
        for (const row of data as any[]) {
            const name = row.Name as string;
            if (!seen.has(name)) {
                seen.add(name);
                unique.push({ fullName: name, simName: name });
            }
            if (unique.length >= limit) break;
        }
        return unique;
}

/**
 * Given a compound name typed by the user, resolve the canonical Aspen
 * CompoundName used in all parameter tables (apv140_*).
 * Checks compounds_master Name first, then Alias (for nicknames like "MeOH").
 */
export async function resolveSimName(
    supabase: SupabaseClient,
    name: string,
): Promise<string | null> {
    if (!name) return null;

    // 1. Try exact Name match in compounds_master (case-insensitive)
    const { data: byName } = await supabase
        .from('compounds_master')
        .select('"Name"')
        .ilike('Name', name.trim())
        .limit(1);
    if (byName && byName.length > 0) return (byName[0] as any).Name as string;

    // 2. Try Alias match (e.g., user typed "MeOH", "CH3OH", a CAS number, etc.)
    const { data: byAlias } = await supabase
        .from('compounds_master')
        .select('"Name"')
        .ilike('Alias', name.trim())
        .limit(1);
    if (byAlias && byAlias.length > 0) return (byAlias[0] as any).Name as string;

    // 3. Fall back – assume user typed the exact Aspen name already
    return name.trim().toUpperCase();
}

// ─── Core therm-data helpers ────────────────────────────────────────────────────

export interface CriticalProperties {
    Tc_K: number;
    Pc_Pa: number;
    omega: number;
}

export interface FetchedCompoundThermData {
    compoundName: string;
    antoine: AntoineParams | null;
    nbp_K: number | null;
    molecularWeight: number | null;
    V_L_m3mol?: number;
    unifacGroups?: UnifacGroupComposition | null;
    criticalProperties?: CriticalProperties | null;
    uniquacParams?: UniquacPureComponentParams | null;
    criticalVolume_m3kmol?: number | null;
    RKTZRA?: number | null; // Aspen Rackett compressibility factor (from RKTZRA column)
    hocProps?: import('./vle-types').HocPureComponentParams | null; // HOC vapor-phase association params
}

/**
 * Fetches compound thermodynamic data from Aspen NIST/APV Supabase tables.
 *
 * Tables queried:
 *   - "apv140_pure_props_wide" → TC (K), PC (Pa), VC (m³/kmol), OMEGA, GMUQR, GMUQQ
 *   - "apv140_plxant_wide"     → 7-parameter extended Antoine (C1–C7, Tmin/Tmax in K)
 *
 * Aspen PLXANT equation:
 *   ln(P* [Pa]) = C1 + C2/(T+C3) + C4·T + C5·ln(T) + C6·T^C7   (T in K)
 */
export async function fetchAndConvertThermData(
    supabase: SupabaseClient,
    compoundName: string
): Promise<FetchedCompoundThermData> {
    if (!compoundName) {
        throw new Error("Compound name cannot be empty for fetchAndConvertThermData.");
    }

    // 1. Fetch pure component properties from Aspen apv140 table
    const { data: pureRows, error: pureError } = await supabase
        .from('apv140_pure_props_wide')
        .select('*')
        .ilike('CompoundName', compoundName)
        .limit(1);

    if (pureError || !pureRows || pureRows.length === 0) {
        console.error(`Error fetching Aspen pure comp data for "${compoundName}":`, pureError?.message);
        throw new Error(`Compound "${compoundName}" not found in apv140_pure_props_wide.`);
    }
    const pureData = pureRows[0];

    // 2. Fetch Antoine parameters from Aspen PLXANT table
    const { data: antoineRows, error: antoineError } = await supabase
        .from('apv140_plxant_wide')
        .select('*')
        .ilike('CompoundName', compoundName)
        .limit(1);

    const antoineData = antoineRows?.[0] ?? null;
    let antoineParams: AntoineParams | null = null;
    if (!antoineError && antoineData) {
        const C1 = antoineData.C1;
        const C2 = antoineData.C2;
        const C3 = antoineData.C3 ?? 0;
        const C4 = antoineData.C4 ?? 0;
        const C5 = antoineData.C5 ?? 0;
        const C6 = antoineData.C6 ?? 0;
        const C7 = antoineData.C7 ?? 1;
        // Aspen PLXANT stores Tmin/Tmax in Kelvin directly
        const Tmin_K = antoineData.Tmin ?? 0;
        const Tmax_K = antoineData.Tmax ?? 10000;

        if (C1 != null && C2 != null) {
            antoineParams = { C1, C2, C3, C4, C5, C6, C7, Tmin_K, Tmax_K };
        } else {
            console.warn(`Antoine parameters incomplete for "${compoundName}"`);
        }
    } else {
        console.warn(`Antoine data not found for "${compoundName}":`, antoineError?.message);
    }

    const Tc_K   = pureData['TC']     ?? null;
    const Pc_Pa  = pureData['PC']     ?? null;
    const omega  = pureData['OMEGA']  ?? null;
    const RKTZRA = pureData['RKTZRA'] ?? null; // Aspen Rackett compressibility factor
    const VC_m3kmol = pureData['VC']  ?? null;  // m³/kmol

    // GMUQR = UNIQUAC r parameter, GMUQQ = UNIQUAC q parameter
    const uniquac_r = pureData['GMUQR'] ?? null;
    const uniquac_q = pureData['GMUQQ'] ?? null;

    // HOC vapor-phase association parameters
    const MUP_stored  = pureData['MUP']          ?? null;   // raw DB value (C·m)
    const RGYR_m      = pureData['RGYR']          ?? null;   // raw DB value (m)
    const HOC_ETA_RAW = pureData['HOC_ETA_SELF']  ?? null;   // null for non-associating

    // Build critical properties
    let criticalProperties: CriticalProperties | null = null;
    if (Tc_K != null && Pc_Pa != null && omega != null) {
        criticalProperties = { Tc_K, Pc_Pa, omega };
    }

    // Build UNIQUAC params (from GMUQR/GMUQQ)
    let uniquacParams: UniquacPureComponentParams | null = null;
    if (uniquac_r != null && uniquac_q != null) {
        uniquacParams = { r: uniquac_r, q: uniquac_q };
    }

    // Spencer-Danner Rackett equation: VL(T) = (R·Tc/Pc)·ZRA^[1+(1-Tr)^(2/7)]
    // Use stored RKTZRA when available; fall back to Yamada-Gunn: ZRA = 0.29056 − 0.08775·ω
    // Reference temperature: 25°C (298.15 K) — consistent with APV Wilson param regression.
    let V_L_m3mol: number | undefined = undefined;
    if (Tc_K != null && Pc_Pa != null && Tc_K > 0 && Pc_Pa > 0) {
        const ZRA = RKTZRA ?? (omega != null ? 0.29056 - 0.08775 * omega : null);
        if (ZRA != null) {
            const R_gas = 8.314_462_618; // J·mol⁻¹·K⁻¹
            const Tr_ref = 298.15 / Tc_K;
            V_L_m3mol = (R_gas * Tc_K / Pc_Pa) * Math.pow(ZRA, 1 + Math.pow(1 - Tr_ref, 2 / 7));
        }
    } else if (VC_m3kmol != null) {
        // Fallback if critical props unavailable: crude estimate from Vc
        V_L_m3mol = (VC_m3kmol / 1000.0) * 0.3;
    }

    // 4. Fetch UNIFAC groups (still from compound_properties legacy table)
    let unifacGroups: UnifacGroupComposition | null = null;
    try {
        const { data: compoundRows, error: compoundError } = await supabase
            .from('compound_properties')
            .select('properties')
            .ilike('name', compoundName)
            .limit(1);

        const compoundData = compoundRows?.[0] ?? null;
        if (!compoundError && compoundData) {
            const props = compoundData.properties as any;
            const unifacData = props?.UnifacVLE?.group || props?.UNIFAC?.group || props?.unifac?.group;
            if (unifacData) {
                const groups: UnifacGroupComposition = {};
                if (Array.isArray(unifacData)) {
                    for (const group of unifacData) {
                        const subgroupId = parseInt(group.id);
                        const count = parseInt(group.value);
                        if (!isNaN(subgroupId) && !isNaN(count) && count > 0) {
                            groups[subgroupId] = count;
                        }
                    }
                } else if (typeof unifacData === 'object' && unifacData.id && unifacData.value) {
                    const subgroupId = parseInt(unifacData.id);
                    const count = parseInt(unifacData.value);
                    if (!isNaN(subgroupId) && !isNaN(count) && count > 0) {
                        groups[subgroupId] = count;
                    }
                }
                if (Object.keys(groups).length > 0) {
                    unifacGroups = groups;
                }
            }
        }
    } catch (e) {
        console.warn(`UNIFAC groups lookup failed for "${compoundName}" (non-fatal):`, e);
    }

    return {
        compoundName,
        antoine: antoineParams,
        nbp_K: null, // Not directly available; can be estimated from Psat if needed
        molecularWeight: null, // Not available in Aspen pure_props_wide
        V_L_m3mol,
        unifacGroups,
        criticalProperties,
        uniquacParams,
        criticalVolume_m3kmol: VC_m3kmol,
        RKTZRA,
        hocProps: (MUP_stored != null && RGYR_m != null) ? {
            MUP_stored,
            RGYR_m,
            eta_self: HOC_ETA_RAW ?? 0,
            mu_D: MUP_stored / 3.33564e-25,
            rd_A: RGYR_m * 1e10,
        } : null,
    };
}

/**
 * Convenience function to build a full CompoundData object from Aspen NIST/APV tables.
 */
export async function fetchCompoundDataFromHysys(
    supabase: SupabaseClient,
    compoundName: string
): Promise<CompoundData> {
    const therm = await fetchAndConvertThermData(supabase, compoundName);

    return {
        name: compoundName,
        cas_number: null,
        molecularWeight: therm.molecularWeight,
        antoine: therm.antoine,
        unifacGroups: therm.unifacGroups,
        prParams: therm.criticalProperties ? {
            Tc_K: therm.criticalProperties.Tc_K,
            Pc_Pa: therm.criticalProperties.Pc_Pa,
            omega: therm.criticalProperties.omega,
        } : null,
        srkParams: therm.criticalProperties ? {
            Tc_K: therm.criticalProperties.Tc_K,
            Pc_Pa: therm.criticalProperties.Pc_Pa,
            omega: therm.criticalProperties.omega,
        } : null,
        uniquacParams: therm.uniquacParams,
        wilsonParams: therm.V_L_m3mol != null ? {
            V_L_m3mol: therm.V_L_m3mol,
            Tc_K: therm.criticalProperties?.Tc_K,
            Pc_Pa: therm.criticalProperties?.Pc_Pa,
            omega: therm.criticalProperties?.omega,
            RKTZRA: therm.RKTZRA ?? undefined,
        } : null,
        hocProps: therm.hocProps ?? null,
    };
}

