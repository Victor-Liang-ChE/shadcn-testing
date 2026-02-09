import type { SupabaseClient } from '@supabase/supabase-js';
import type { AntoineParams, UnifacGroupComposition, UniquacPureComponentParams, CompoundData } from './vle-types';

// ─── Alias helpers (HYSYS ALIASES table) ───────────────────────────────────────

/** An entry from the HYSYS ALIASES table. */
export interface CompoundAlias {
    fullName: string;
    simName: string;
}

/**
 * Format compound names for display by replacing underscores with spaces.
 * E.g., "Water_Vapor" becomes "Water Vapor".
 */
export function formatCompoundName(name: string): string {
    return name.replace(/_/g, ' ');
}

/**
 * Search the HYSYS ALIASES table for compounds whose FullName matches the
 * user-typed prefix.  Returns `{fullName, simName}` pairs so the UI can
 * display the human-readable name while keeping the internal SimName for
 * look-ups in the other HYSYS tables.
 */
export async function fetchCompoundSuggestions(
    supabase: SupabaseClient,
    query: string,
    limit = 5,
): Promise<CompoundAlias[]> {
    if (!query || query.trim().length < 1) return [];

    const { data, error } = await supabase
        .from('HYSYS ALIASES')
        .select('"FullName","SimName"')
        .ilike('FullName', `${query.trim()}%`)
        .limit(limit);

    if (error || !data) {
        console.error('fetchCompoundSuggestions error:', error?.message);
        return [];
    }

    return data.map((row: any) => ({
        fullName: row.FullName as string,
        simName: row.SimName as string,
    }));
}

/**
 * Given a name (could be FullName _or_ SimName), resolve the SimName that the
 * rest of the HYSYS tables expect.  Tries FullName first, then SimName, so
 * that users can type either.
 */
export async function resolveSimName(
    supabase: SupabaseClient,
    name: string,
): Promise<string | null> {
    if (!name) return null;

    // 1. Try FullName → SimName
    const { data: byFull } = await supabase
        .from('HYSYS ALIASES')
        .select('"SimName"')
        .ilike('FullName', name.trim())
        .limit(1);

    if (byFull && byFull.length > 0) return byFull[0].SimName as string;

    // 2. Try SimName directly (backward compat – user typed a SimName)
    const { data: bySim } = await supabase
        .from('HYSYS ALIASES')
        .select('"SimName"')
        .ilike('SimName', name.trim())
        .limit(1);

    if (bySim && bySim.length > 0) return bySim[0].SimName as string;

    // 3. Fall back – assume the caller typed the exact SimName already
    return name.trim();
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
    V_L_m3mol?: number; // Liquid Molar Volume in m^3/mol for Wilson
    unifacGroups?: UnifacGroupComposition | null; // For UNIFAC
    criticalProperties?: CriticalProperties | null; // For PR and SRK
    uniquacParams?: UniquacPureComponentParams | null; // For UNIQUAC
    criticalVolume_m3kmol?: number | null; // For Kij estimation
}

/**
 * Fetches compound thermodynamic data from HYSYS Supabase tables by component name.
 * 
 * Tables queried:
 *   - "HYSYS PROPS PARAMS" → Tc, Pc, ω, MW, NBP, UNIQUAC r/q, Wilson V, Vc
 *   - "HYSYS ANTOINE" → 7-parameter extended Antoine coefficients (A–G, Tmin, Tmax)
 *   - UNIFAC tables (unchanged) → UNIFAC group composition via compound_properties fallback
 */
export async function fetchAndConvertThermData(
    supabase: SupabaseClient,
    compoundName: string
): Promise<FetchedCompoundThermData> {
    if (!compoundName) {
        throw new Error("Compound name cannot be empty for fetchAndConvertThermData.");
    }

    // 1. Fetch pure component properties from HYSYS table
    const { data: pureRows, error: pureError } = await supabase
        .from('HYSYS PROPS PARAMS')
        .select('*')
        .ilike('Component', compoundName)
        .limit(1);

    if (pureError || !pureRows || pureRows.length === 0) {
        console.error(`Error fetching HYSYS pure comp data for "${compoundName}":`, pureError?.message);
        throw new Error(`Compound "${compoundName}" not found in HYSYS PROPS PARAMS.`);
    }
    const pureData = pureRows[0];

    // 2. Fetch Antoine parameters from HYSYS ANTOINE table
    const { data: antoineRows, error: antoineError } = await supabase
        .from('HYSYS ANTOINE')
        .select('*')
        .ilike('Component', compoundName)
        .limit(1);

    const antoineData = antoineRows?.[0] ?? null;
    let antoineParams: AntoineParams | null = null;
    if (!antoineError && antoineData) {
        const A = antoineData.A;
        const B = antoineData.B;
        const C = antoineData.C;
        const D = antoineData.D ?? 0;
        const E = antoineData.E ?? 0;
        const F = antoineData.F ?? 0;
        const G = antoineData.G ?? 0;
        const Tmin = antoineData.Tmin ?? 0;
        const Tmax = antoineData.Tmax ?? 10000;

        if (A != null && B != null && C != null) {
            // HYSYS stores Tmin/Tmax in °C; convert to K for downstream use
            antoineParams = { A, B, C, D, E, F, G, Tmin_K: Tmin + 273.15, Tmax_K: Tmax + 273.15 };
        } else {
            console.warn(`Antoine parameters incomplete for "${compoundName}"`);
        }
    } else {
        console.warn(`Antoine data not found for "${compoundName}":`, antoineError?.message);
    }

    // 3. Extract pure component properties
    // HYSYS stores temperatures in °C — convert to K
    const Tc_C = pureData['CriticalTemperature'] ?? null;
    const Tc_K = Tc_C != null ? Tc_C + 273.15 : null; // °C → K
    const Pc_kPa = pureData['CriticalPressure'] ?? null;  // kPa
    const omega = pureData['Acentricity'] ?? null;
    const MW = pureData['MolecularWeight'] ?? null;
    const nbp_C = pureData['NormalBoilingPoint'] ?? null;
    const nbp_K = nbp_C != null ? nbp_C + 273.15 : null; // °C → K
    const uniquac_r = pureData['UNIQUAC_r'] ?? null;
    const uniquac_q = pureData['UNIQUAC_q'] ?? null;
    const wilsonVol_m3kmol = pureData['WilsonVolume'] ?? null; // m³/kmol
    const critVol = pureData['CriticalVolume'] ?? null; // m³/kmol

    // Build critical properties (Pc: kPa → Pa)
    let criticalProperties: CriticalProperties | null = null;
    if (Tc_K != null && Pc_kPa != null && omega != null) {
        criticalProperties = {
            Tc_K,
            Pc_Pa: Pc_kPa * 1000, // kPa → Pa
            omega,
        };
    }

    // Build UNIQUAC params
    let uniquacParams: UniquacPureComponentParams | null = null;
    if (uniquac_r != null && uniquac_q != null) {
        uniquacParams = { r: uniquac_r, q: uniquac_q };
    }

    // Wilson molar volume: m³/kmol → m³/mol (÷ 1000)
    let V_L_m3mol: number | undefined = undefined;
    if (wilsonVol_m3kmol != null) {
        V_L_m3mol = wilsonVol_m3kmol / 1000.0;
    }

    // 4. Fetch UNIFAC groups (still from compound_properties table since UNIFAC tables are unchanged)
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
        nbp_K: nbp_K,
        molecularWeight: MW,
        V_L_m3mol,
        unifacGroups,
        criticalProperties,
        uniquacParams,
        criticalVolume_m3kmol: critVol,
    };
}

/**
 * Convenience function to build a full CompoundData object from HYSYS tables.
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
        wilsonParams: therm.V_L_m3mol != null ? { V_L_m3mol: therm.V_L_m3mol } : null,
    };
}
