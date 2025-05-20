import type { SupabaseClient } from '@supabase/supabase-js';
import type { AntoineParams, UnifacGroupComposition, UniquacPureComponentParams } from './vle-types';

export interface CriticalProperties { // Generic structure for Tc, Pc, omega
    Tc_K: number;
    Pc_Pa: number;
    omega: number;
}

export interface FetchedCompoundThermData {
    casNumber: string;
    antoine: AntoineParams | null;
    nbp_K: number | null;
    V_L_m3mol?: number; // Liquid Molar Volume in m^3/mol for Wilson
    unifacGroups?: UnifacGroupComposition | null; // For UNIFAC
    criticalProperties?: CriticalProperties | null; // For PR and SRK
    uniquacParams?: UniquacPureComponentParams | null; // For UNIQUAC
    // any other properties you might fetch
}

export async function fetchAndConvertThermData(
    supabase: SupabaseClient,
    casn: string
): Promise<FetchedCompoundThermData> {
    if (!casn) {
        throw new Error("CAS number cannot be empty for fetchAndConvertThermData.");
    }

    let antoineParams: AntoineParams | null = null;
    let nbp_K_value: number | null = null;
    let V_L_m3mol_value: number | undefined = undefined;
    let unifacGroups_value: UnifacGroupComposition | null = null;
    let criticalProperties_value: CriticalProperties | null = null; // Changed from prSrkParams_value
    let uniquacParams_value: UniquacPureComponentParams | null = null;

    // 1. Fetch compound_id from CAS number
    const { data: compoundIdData, error: compoundIdError } = await supabase
        .from('compounds')
        .select('id, name') // Fetch name too for logging if needed
        .eq('cas_number', casn)
        .single();

    if (compoundIdError || !compoundIdData) {
        console.error(`Error fetching compound ID for CAS ${casn}:`, compoundIdError?.message);
        throw new Error(`Compound with CAS ${casn} not found.`);
    }
    const compound_id = compoundIdData.id;

    // 2. Fetch properties from compound_properties table
    const { data: propertiesData, error: propertiesError } = await supabase
        .from('compound_properties')
        .select('properties')
        .eq('compound_id', compound_id)
        .eq('source', 'chemsep1') // Example: Filter by a specific source if necessary
        .maybeSingle();

    if (propertiesError) {
        console.warn(`Error fetching compound_properties for compound_id ${compound_id} (CAS ${casn}): ${propertiesError.message}`);
    }

    if (propertiesData && propertiesData.properties) {
        const props = propertiesData.properties as any;

        // Extract Antoine Parameters
        if (props.Antoine && typeof props.Antoine.A === 'number' && typeof props.Antoine.B === 'number' && typeof props.Antoine.C === 'number') {
            const rawUnits = props.Antoine.units;
            let resolvedUnits: string;

            if (rawUnits === 'Pa' || rawUnits === 'mmHg' || rawUnits === 'bar') {
                resolvedUnits = rawUnits;
            } else {
                if (rawUnits && typeof rawUnits === 'string') {
                    console.warn(`Unrecognized Antoine units '${rawUnits}' for CAS ${casn}. Defaulting to 'Pa'.`);
                } else {
                    console.warn(`Antoine units not specified or invalid for CAS ${casn}. Defaulting to 'Pa'.`);
                }
                resolvedUnits = 'Pa';
            }

            antoineParams = {
                A: props.Antoine.A,
                B: props.Antoine.B,
                C: props.Antoine.C,
                Tmin_K: props.Antoine.Tmin?.value,
                Tmax_K: props.Antoine.Tmax?.value,
                Units: resolvedUnits,
            };
        } else {
            console.warn(`Antoine parameters missing or incomplete for CAS ${casn}`);
        }

        // Extract Normal Boiling Point
        if (props["Normal boiling point"] && typeof props["Normal boiling point"].value === 'number') {
            nbp_K_value = props["Normal boiling point"].value;
        } else {
            console.warn(`Normal Boiling Point (nbp_K) missing for CAS ${casn}`);
        }

        // Extract Liquid Molar Volume (V_L_m3mol)
        if (props["Wilson volume"] && typeof props["Wilson volume"].value === 'number' && props["Wilson volume"].units === "m3/kmol") {
            V_L_m3mol_value = props["Wilson volume"].value / 1000.0;
        } else if (props["Liquid molar volume at normal boiling point"] && typeof props["Liquid molar volume at normal boiling point"].value === 'number' && props["Liquid molar volume at normal boiling point"].units === "m3/kmol") {
            V_L_m3mol_value = props["Liquid molar volume at normal boiling point"].value / 1000.0;
        } else if (props["Liquid density"] && typeof props["Liquid density"].A === 'number' && props["Liquid density"].units === "kmol/m3") {
            console.warn(`Direct V_L not found for CAS ${casn}. "Liquid density" might require T-dependent calculation not yet implemented.`);
        }

        if (V_L_m3mol_value === undefined) {
            console.warn(`Liquid Molar Volume (V_L_m3mol) could not be extracted for CAS ${casn}`);
        }

        // Extract UNIFAC Groups
        if (props.elements_composition?.UNIFAC && typeof props.elements_composition.UNIFAC === 'object') {
            const groups: UnifacGroupComposition = {};
            for (const key in props.elements_composition.UNIFAC) {
                const subgroupId = parseInt(key);
                const count = parseInt(props.elements_composition.UNIFAC[key]);
                if (!isNaN(subgroupId) && !isNaN(count) && count > 0) {
                    groups[subgroupId] = count;
                }
            }
            if (Object.keys(groups).length > 0) {
                unifacGroups_value = groups;
            } else {
                console.warn(`UNIFAC groups found but empty or invalid for CAS ${casn}`);
            }
        } else {
            console.warn(`UNIFAC groups (props.elements_composition.UNIFAC) missing for CAS ${casn}`);
        }

        // Extract Critical Properties (Tc, Pc, omega)
        const tcProp = props["Critical temperature"];
        const pcProp = props["Critical pressure"];
        const omegaProp = props["Acentric factor"];

        if (tcProp && typeof tcProp.value === 'number' &&
            pcProp && typeof pcProp.value === 'number' && typeof pcProp.units === 'string' &&
            omegaProp && typeof omegaProp.value === 'number') {
            
            const Tc_K = tcProp.value;
            let Pc_Pa = pcProp.value;
            const pcUnitsLower = pcProp.units.toLowerCase();

            if (pcUnitsLower === 'kpa') Pc_Pa *= 1000;
            else if (pcUnitsLower === 'bar') Pc_Pa *= 100000;
            else if (pcUnitsLower !== 'pa') {
                console.warn(`Unrecognized critical pressure units '${pcProp.units}' for CAS ${casn}. Assuming Pa if not kPa or bar.`);
            }
            
            const omega = omegaProp.value;
            criticalProperties_value = { Tc_K, Pc_Pa, omega };
        } else {
            console.warn(`Critical properties (Tc, Pc, omega) missing or incomplete for CAS ${casn}`);
        }

        // Extract UNIQUAC Parameters (r, q)
        const rPropUQ = props["UNIQUAC r"] || props["Van der Waals volume"];
        const qPropUQ = props["UNIQUAC q"] || props["Van der Waals area"];

        if (rPropUQ && typeof rPropUQ.value === 'number' &&
            qPropUQ && typeof qPropUQ.value === 'number') {
            uniquacParams_value = { r: rPropUQ.value, q: qPropUQ.value };
        } else {
            console.warn(`UNIQUAC parameters (r, q) missing or incomplete for CAS ${casn}`);
        }

    } else {
        console.warn(`No properties data found for compound_id ${compound_id} (CAS ${casn}).`);
    }

    return {
        casNumber: casn,
        antoine: antoineParams,
        nbp_K: nbp_K_value,
        V_L_m3mol: V_L_m3mol_value,
        unifacGroups: unifacGroups_value,
        criticalProperties: criticalProperties_value, // Use generic criticalProperties
        uniquacParams: uniquacParams_value,
    };
}
