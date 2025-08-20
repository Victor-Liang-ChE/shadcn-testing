import type { SupabaseClient } from '@supabase/supabase-js';
import type { AntoineParams, UnifacGroupComposition, UniquacPureComponentParams } from './vle-types';
import { parseCoefficient } from './property-equations';

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

    // 1. Fetch compound data directly from compound_properties table using CAS number
    const { data: compoundData, error: compoundError } = await supabase
        .from('compound_properties')
        .select('name, properties')
        .eq('properties->CAS->>value', casn)
        .single();

    if (compoundError || !compoundData) {
        console.error(`Error fetching compound data for CAS ${casn}:`, compoundError?.message);
        throw new Error(`Compound with CAS ${casn} not found.`);
    }

    const props = compoundData.properties as any;

    // Extract Antoine Parameters
    const antoineData = props.AntoineVaporPressure || props.Antoine;
    if (antoineData) {
        const A = parseCoefficient(antoineData.A);
        const B = parseCoefficient(antoineData.B);
        const C = parseCoefficient(antoineData.C);
        const Tmin_K = parseCoefficient(antoineData.Tmin?.value ?? antoineData.Tmin);
        const Tmax_K = parseCoefficient(antoineData.Tmax?.value ?? antoineData.Tmax);
        const eqno = parseCoefficient(antoineData.eqno?.value ?? antoineData.eqno);

        const rawUnits = typeof antoineData.units === 'string' ? antoineData.units : undefined;
        // Accept common pressure units; default to 'Pa' if unspecified/unknown
        let resolvedUnits: string = 'Pa';
        if (rawUnits) {
            const u = rawUnits.toLowerCase();
            if (u === 'pa') resolvedUnits = 'Pa';
            else if (u === 'kpa') resolvedUnits = 'kPa';
            else if (u === 'bar') resolvedUnits = 'bar';
            else if (u === 'mmhg' || u === 'torr') resolvedUnits = 'mmHg';
            else if (u === 'atm') resolvedUnits = 'atm';
            else {
                console.warn(`Unrecognized Antoine units '${rawUnits}' for CAS ${casn}. Defaulting to 'Pa'.`);
            }
        }

        if (A !== null && B !== null && C !== null) {
            antoineParams = {
                A,
                B,
                C,
                Tmin_K: Tmin_K ?? undefined as any,
                Tmax_K: Tmax_K ?? undefined as any,
                Units: resolvedUnits,
                EquationNo: eqno ?? undefined,
            };
        } else {
            console.warn(`Antoine parameters missing or incomplete for CAS ${casn}`);
        }
    }

    // Extract Normal Boiling Point
    if (props.NormalBoilingPointTemperature && typeof props.NormalBoilingPointTemperature.value === 'number') {
        nbp_K_value = props.NormalBoilingPointTemperature.value;
    } else {
        console.warn(`Normal Boiling Point (nbp_K) missing for CAS ${casn}`);
    }

    // Extract Liquid Molar Volume (V_L_m3mol)
    if (props.WilsonVolume && props.WilsonVolume.units === "m3/kmol") {
        const wilsonVol = parseCoefficient(props.WilsonVolume.value);
        if (wilsonVol !== null) {
            V_L_m3mol_value = wilsonVol / 1000.0; // Convert from m3/kmol to m3/mol
        }
    } else if (props.LiquidVolumeAtNormalBoilingPoint && props.LiquidVolumeAtNormalBoilingPoint.units === "m3/kmol") {
        const liqVol = parseCoefficient(props.LiquidVolumeAtNormalBoilingPoint.value);
        if (liqVol !== null) {
            V_L_m3mol_value = liqVol / 1000.0; // Convert from m3/kmol to m3/mol
        }
    } else {
        console.warn(`Direct V_L not found for CAS ${casn}.`);
    }

    if (V_L_m3mol_value === undefined) {
        console.warn(`Liquid Molar Volume (V_L_m3mol) could not be extracted for CAS ${casn}`);
    }

    // Extract UNIFAC Groups
    const unifacData = props.UnifacVLE?.group || props.UNIFAC?.group || props.unifac?.group;
    if (unifacData) {
        const groups: UnifacGroupComposition = {};
        
        if (Array.isArray(unifacData)) {
            // Handle array format: [{"id": "1", "value": "2"}, ...]
            for (const group of unifacData) {
                const subgroupId = parseInt(group.id);
                const count = parseInt(group.value);
                if (!isNaN(subgroupId) && !isNaN(count) && count > 0) {
                    groups[subgroupId] = count;
                }
            }
        } else if (typeof unifacData === 'object' && unifacData.id && unifacData.value) {
            // Handle single object format: {"id": "17", "value": "1"}
            const subgroupId = parseInt(unifacData.id);
            const count = parseInt(unifacData.value);
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
        console.warn(`UNIFAC groups missing for CAS ${casn}`);
    }

    // Extract Critical Properties (Tc, Pc, omega)
    const tcProp = props.CriticalTemperature;
    const pcProp = props.CriticalPressure;
    const omegaProp = props.AcentricityFactor;

    if (tcProp && pcProp && omegaProp) {
        const Tc_K = parseCoefficient(tcProp.value);
        const Pc_raw = parseCoefficient(pcProp.value);
        const omega = parseCoefficient(omegaProp.value);
        
        if (Tc_K !== null && Pc_raw !== null && omega !== null && typeof pcProp.units === 'string') {
            let Pc_Pa = Pc_raw;
            const pcUnitsLower = pcProp.units.toLowerCase();

            if (pcUnitsLower === 'kpa') Pc_Pa *= 1000;
            else if (pcUnitsLower === 'bar') Pc_Pa *= 100000;
            else if (pcUnitsLower !== 'pa') {
                console.warn(`Unrecognized critical pressure units '${pcProp.units}' for CAS ${casn}. Assuming Pa if not kPa or bar.`);
            }

            criticalProperties_value = { Tc_K, Pc_Pa, omega };
        } else {
            console.warn(`Critical properties (Tc, Pc, omega) could not be parsed for CAS ${casn}`);
        }
    } else {
        console.warn(`Critical properties (Tc, Pc, omega) missing or incomplete for CAS ${casn}`);
    }

    // Extract UNIQUAC Parameters (r, q)
    const rPropUQ = props.UniquacR || props.VanDerWaalsVolume;
    const qPropUQ = props.UniquacQ || props.VanDerWaalsArea;

    if (rPropUQ && qPropUQ) {
        const r = parseCoefficient(rPropUQ.value);
        const q = parseCoefficient(qPropUQ.value);
        if (r !== null && q !== null) {
            uniquacParams_value = { r, q };
        } else {
            console.warn(`UNIQUAC parameters (r, q) could not be parsed for CAS ${casn}`);
        }
    } else {
        console.warn(`UNIQUAC parameters (r, q) missing or incomplete for CAS ${casn}`);
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
