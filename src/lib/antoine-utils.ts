import type { SupabaseClient } from '@supabase/supabase-js';
import type { AntoineParams } from './vle-types';

interface DbAntoineValue {
    value: number;
    units: string;
}

interface DbNbpValue { // Added for NBP
    value: number;
    units: string; // Expecting K
}

interface DbAntoineDetails {
    A: number;
    B: number;
    C: number;
    Tmin: DbAntoineValue;
    Tmax: DbAntoineValue;
    units: 'Pa' | string; // Expecting Pa from DB for ln(P_Pa) form
}

interface DbCompoundProperties {
    properties: {
        Antoine?: DbAntoineDetails;
        "Normal boiling point"?: DbNbpValue; // Added for NBP
        // ... other properties
    };
}

/**
 * New return type for fetched compound thermo data.
 */
export interface FetchedCompoundThermData {
    antoine: AntoineParams;
    nbp_K: number | null; // Normal Boiling Point in Kelvin
}

/**
 * Fetches Antoine parameters and Normal Boiling Point for a given CAS number from Supabase.
 * Antoine params are converted to log10(P_mmHg) = A - B / (T_°C + C).
 */
export async function fetchAndConvertThermData( // Renamed function
    supabase: SupabaseClient,
    casn: string
): Promise<FetchedCompoundThermData> {
    if (!casn) {
        throw new Error('CAS number is required to fetch thermo data.');
    }

    // Step 1: Get compound_id from CAS number
    const { data: compoundData, error: compoundError } = await supabase
        .from('compounds')
        .select('id')
        .eq('cas_number', casn)
        .single();

    if (compoundError || !compoundData) {
        console.error(`Error fetching compound ID for CAS ${casn}:`, compoundError);
        throw new Error(`Compound with CAS ${casn} not found. ${compoundError?.message || ''}`);
    }
    const compoundId = compoundData.id;

    // Step 2: Get properties using compound_id
    const { data: propertiesData, error: propertiesError } = await supabase
        .from('compound_properties')
        .select('properties')
        .eq('compound_id', compoundId)
        .limit(1) // Take the first one if multiple sources exist for a compound
        .single(); // Expect a single record or null

    if (propertiesError || !propertiesData) {
        console.error(`Error fetching properties for compound ID ${compoundId} (CAS ${casn}):`, propertiesError);
        throw new Error(`Properties not found for compound CAS ${casn}. ${propertiesError?.message || ''}`);
    }

    const dbProps = propertiesData as DbCompoundProperties;
    const dbAntoine = dbProps.properties?.Antoine;
    const dbNbp = dbProps.properties?.["Normal boiling point"]; // Access NBP

    if (!dbAntoine || typeof dbAntoine.A !== 'number' || typeof dbAntoine.B !== 'number' || typeof dbAntoine.C !== 'number') {
        throw new Error(`Antoine parameters not found or incomplete in DB for CAS ${casn}.`);
    }

    let nbp_K: number | null = null;
    if (dbNbp && typeof dbNbp.value === 'number') {
        if (dbNbp.units === 'K') {
            nbp_K = dbNbp.value;
        } else {
            console.warn(`NBP for CAS ${casn} is in units ${dbNbp.units}, expected K. NBP not used.`);
            // Potentially add conversion here if other units are common, e.g., °C
        }
    } else {
        console.warn(`Normal Boiling Point data not found or invalid in DB for CAS ${casn}.`);
    }

    // Antoine Conversion (same as before)
    // DB: ln(P_Pa) = A_db - B_db / (C_db + T_K)
    // Target: log10(P_mmHg) = A_ui - B_ui / (C_ui + T_°C)
    const A_converted = (dbAntoine.A / Math.log(10)) - Math.log10(133.3223684);
    const B_converted = dbAntoine.B / Math.log(10);
    // Assuming C_db is for (T_K + C_db) form, convert to (T_C + C_converted)
    // T_K + C_db = (T_C + 273.15) + C_db
    // We want T_C + C_converted, so C_converted = C_db + 273.15
    const C_converted = dbAntoine.C + 273.15;

    const antoineParams: AntoineParams = {
        A: A_converted,
        B: B_converted,
        C: C_converted,
        Tmin_K: dbAntoine.Tmin?.value || 0, // Use DB Tmin if available
        Tmax_K: dbAntoine.Tmax?.value || 1000, // Use DB Tmax if available
        Units: 'mmHg', // The converted parameters are for this unit system
    };

    return {
        antoine: antoineParams,
        nbp_K: nbp_K,
    };
}
