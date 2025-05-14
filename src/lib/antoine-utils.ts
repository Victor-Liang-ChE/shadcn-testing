import type { SupabaseClient } from '@supabase/supabase-js';
import type { AntoineParams } from './vle-types';

export interface FetchedCompoundThermData {
    casNumber: string;
    antoine: AntoineParams | null;
    nbp_K: number | null;
    V_L_m3mol?: number; // Liquid Molar Volume in m^3/mol
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
    // Assuming 'DDB_Primary' or a similar source contains the most relevant data
    // Adjust 'source' if your primary data source for these properties is different
    const { data: propertiesData, error: propertiesError } = await supabase
        .from('compound_properties')
        .select('properties')
        .eq('compound_id', compound_id)
        .eq('source', 'chemsep1') // Example: Filter by a specific source if necessary
        .maybeSingle(); // Use maybeSingle as properties might not exist for all compounds/sources

    if (propertiesError) {
        console.warn(`Error fetching compound_properties for compound_id ${compound_id} (CAS ${casn}): ${propertiesError.message}`);
        // Continue, as some data might still be usable or defaults applied
    }

    if (propertiesData && propertiesData.properties) {
        const props = propertiesData.properties as any; // Cast to any for easier access to dynamic keys

        // Extract Antoine Parameters
        if (props.Antoine && typeof props.Antoine.A === 'number' && typeof props.Antoine.B === 'number' && typeof props.Antoine.C === 'number') {
            const rawUnits = props.Antoine.units;
            let resolvedUnits: string;

            if (rawUnits === 'Pa' || rawUnits === 'mmHg' || rawUnits === 'bar') {
                resolvedUnits = rawUnits;
            } else {
                if (rawUnits && typeof rawUnits === 'string') { // If units were provided but not recognized
                    console.warn(`Unrecognized Antoine units '${rawUnits}' for CAS ${casn}. Defaulting to 'Pa'.`);
                } else { // If units were undefined, null, or not a string
                    console.warn(`Antoine units not specified or invalid for CAS ${casn}. Defaulting to 'Pa'.`);
                }
                resolvedUnits = 'Pa'; // Default to 'Pa'
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
        // Priority: "Wilson volume", then "Liquid molar volume at normal boiling point"
        // Both are expected in m3/kmol, convert to m3/mol
        if (props["Wilson volume"] && typeof props["Wilson volume"].value === 'number' && props["Wilson volume"].units === "m3/kmol") {
            V_L_m3mol_value = props["Wilson volume"].value / 1000.0; // kmol to mol
        } else if (props["Liquid molar volume at normal boiling point"] && typeof props["Liquid molar volume at normal boiling point"].value === 'number' && props["Liquid molar volume at normal boiling point"].units === "m3/kmol") {
            V_L_m3mol_value = props["Liquid molar volume at normal boiling point"].value / 1000.0; // kmol to mol
        } else if (props["Liquid density"] && typeof props["Liquid density"].A === 'number' && props["Liquid density"].units === "kmol/m3") {
            // Fallback: Calculate from liquid density if available (e.g. using the first coefficient 'A' if it represents density directly)
            // This is a simplification; density might be temperature-dependent.
            // For this example, let's assume 'A' is a direct density value or a primary term in a correlation.
            // If props["Liquid density"].eqno === 105 (as in example), it's a complex equation.
            // For simplicity, we'll only use it if it's a single value or simple to interpret.
            // The provided example for "Liquid density" has eqno: 105, which is complex.
            // A simpler structure might be: "Liquid density": { "value": 10.0, "units": "kmol/m3" }
            // For now, we'll prefer the direct volume values.
            // If you need to calculate from density equations, that logic would go here.
            console.warn(`Direct V_L not found for CAS ${casn}. "Liquid density" might require T-dependent calculation not yet implemented.`);
        }
        
        if (V_L_m3mol_value === undefined) {
             console.warn(`Liquid Molar Volume (V_L_m3mol) could not be extracted for CAS ${casn}`);
        }

    } else {
        console.warn(`No properties data found for compound_id ${compound_id} (CAS ${casn}).`);
    }

    return {
        casNumber: casn,
        antoine: antoineParams,
        nbp_K: nbp_K_value,
        V_L_m3mol: V_L_m3mol_value,
    };
}
