'use client';

import React, { useState, useEffect } from 'react';
import { createClient, SupabaseClient } from '@supabase/supabase-js';
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Terminal } from "lucide-react";
import { useTheme } from "next-themes";

import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsCoreOption } from 'echarts/core';
import { LineChart, LineSeriesOption } from 'echarts/charts';
import {
    TitleComponent,
    TooltipComponent,
    GridComponent,
    LegendComponent,
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

// Import UNIFAC calculation logic and types
import {
    calculatePsat_Pa as libCalculatePsat_Pa,
    calculateUnifacGamma as libCalculateUnifacGamma,
    calculateBubbleTemperature as calculateBubbleTemperatureUnifac, // Renamed for clarity
    calculateBubblePressure as calculateBubblePressureUnifac,   // Renamed for clarity
    fetchUnifacInteractionParams,
    // type AntoineParams, // Moved
    // type UnifacGroupComposition, // Moved or local to unifac.ts
    // type CompoundData, // Moved
    type UnifacParameters, // Stays in unifac.ts
    // type BubbleDewResult, // Moved
    // type PrPureComponentParams // Moved
} from '@/lib/vle-calculations-unifac';

// Import NRTL calculation logic and types
import {
    fetchNrtlParameters,
    calculateNrtlGamma,
    calculateBubbleTemperatureNrtl,
    calculateBubblePressureNrtl,
    // calculateDewTemperatureNrtl, // Keep if you implement
    // calculateDewPressureNrtl,     // Keep if you implement
    type NrtlInteractionParams
} from '@/lib/vle-calculations-nrtl';

// Import PR calculation logic and types
import {
    fetchPrInteractionParams,
    calculateBubbleTemperaturePr,
    calculateBubblePressurePr,
    type PrInteractionParams
} from '@/lib/vle-calculations-pr';

// Import SRK calculation logic and types // ADDED
import {
    fetchSrkInteractionParams,
    calculateBubbleTemperatureSrk,
    calculateBubblePressureSrk,
    type SrkInteractionParams
} from '@/lib/vle-calculations-srk';

// Import UNIQUAC calculation logic and types // ADDED
import {
    fetchUniquacInteractionParams,
    calculateBubbleTemperatureUniquac,
    calculateBubblePressureUniquac,
    type UniquacInteractionParams
} from '@/lib/vle-calculations-uniquac';

// Import Wilson calculation logic and types // ADDED
import {
    fetchWilsonInteractionParams,
    calculateBubbleTemperatureWilson,
    calculateBubblePressureWilson,
    type WilsonInteractionParams
} from '@/lib/vle-calculations-wilson';

// Import Shared VLE Types
import type {
    AntoineParams,
    PrPureComponentParams,
    SrkPureComponentParams,
    UniquacPureComponentParams,
    WilsonPureComponentParams, // ADDED
    CompoundData,
    BubbleDewResult
} from '@/lib/vle-types';


// Register ECharts components
echarts.use([
    TitleComponent,
    TooltipComponent,
    GridComponent,
    LegendComponent,
    LineChart,
    CanvasRenderer
]);

// --- Supabase Client Setup ---
const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL || 'YOUR_SUPABASE_URL';
const supabaseAnonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY || 'YOUR_SUPABASE_ANON_KEY';

let supabase: SupabaseClient;
try {
    supabase = createClient(supabaseUrl, supabaseAnonKey);
    console.log("Supabase client initialized.");
} catch (error) {
    console.error("Error initializing Supabase client:", error);
}


// --- React Component ---
type DiagramType = "Txy" | "Pxy" | "xy_constP" | "xy_constT";
type AzeotropeScanType = 'vs_P_find_T' | 'vs_T_find_P';
type FluidPackageType = 'unifac' | 'nrtl' | 'pr' | 'srk' | 'uniquac' | 'wilson'; // ADDED: 'wilson'

export default function VleCalculatorPage() {
    const [comp1Name, setComp1Name] = useState<string>("Methanol"); // Changed default for NRTL testing
    const [comp2Name, setComp2Name] = useState<string>("Water");   // Changed default for NRTL testing
    const [pressureKPa, setPressureKPa] = useState<number>(101.325);
    const [temperatureK, setTemperatureK] = useState<number>(343.15); // e.g. 70C for Ethanol/Water
    const [isLoading, setIsLoading] = useState<boolean>(false);
    const [error, setError] = useState<string | null>(null);
    const [results, setResults] = useState<BubbleDewResult[]>([]);
    const [logMessages, setLogMessages] = useState<string[]>([]);
    const { resolvedTheme } = useTheme();
    const [vleChartOptions, setVleChartOptions] = useState<EChartsCoreOption>({});
    const [diagramType, setDiagramType] = useState<DiagramType>("Txy");
    const [fluidPackage, setFluidPackage] = useState<FluidPackageType>('unifac'); // Default to unifac

    // --- Azeotrope Finder States ---
    const [azeotropeFinderActive, setAzeotropeFinderActive] = useState<boolean>(false);
    const [azeotropeScanType, setAzeotropeScanType] = useState<AzeotropeScanType>('vs_P_find_T');
    const [scanSteps, setScanSteps] = useState<number>(100); // MODIFIED DEFAULT VALUE
    const [azeotropeScanData, setAzeotropeScanData] = useState<{ scanVal: number, x_az: number, dependentVal: number }[]>([]);


    // --- Logging Hook (remains the same) ---
    useEffect(() => {
        const originalLog = console.log;
        const originalWarn = console.warn;
        const originalError = console.error;

        const addToLog = (level: string, ...args: any[]) => {
            const message = args.map(arg => typeof arg === 'object' ? JSON.stringify(arg) : arg).join(' ');
            setLogMessages(prev => [...prev, `[${level.toUpperCase()}] ${new Date().toLocaleTimeString()}: ${message}`].slice(-100)); // Keep last 100 messages
        };

        console.log = (...args) => {
            originalLog.apply(console, args);
            addToLog('info', ...args);
        };
        console.warn = (...args) => {
            originalWarn.apply(console, args);
            addToLog('warn', ...args);
        };
        console.error = (...args) => {
            originalError.apply(console, args);
            addToLog('error', ...args);
        };

        // Cleanup function to restore original console methods
        return () => {
            console.log = originalLog;
            console.warn = originalWarn;
            console.error = originalError;
        };
    }, []);


    // --- Data Fetching Functions (fetchCompoundData) ---
    async function fetchCompoundData(compoundName: string): Promise<CompoundData | null> { // Uses imported CompoundData
        if (!supabase) { throw new Error("Supabase client not initialized."); }
        console.log(`Fetching data for ${compoundName}...`);
        setError(null);

        try {
            // 1. Find compound ID and CAS Number
            const { data: compoundDbData, error: compoundError } = await supabase
                .from('compounds')
                .select('id, name, cas_number') // Fetch cas_number
                .ilike('name', compoundName)
                .limit(1);

            if (compoundError) throw new Error(`Supabase compound query error: ${compoundError.message}`);
            if (!compoundDbData || compoundDbData.length === 0) {
                throw new Error(`Compound '${compoundName}' not found.`);
            }
            const compoundId = compoundDbData[0].id;
            const foundName = compoundDbData[0].name;
            const casNumber = compoundDbData[0].cas_number; // Get CAS number
            console.log(`Found ${foundName} (ID: ${compoundId}, CAS: ${casNumber})`);

            if (fluidPackage === 'nrtl' && !casNumber) {
                // For NRTL, CAS number is essential for fetching parameters.
                // UNIFAC might proceed if groups are found even without CAS, but good to have.
                console.warn(`CAS number for ${foundName} is missing, which might be required for the selected fluid package.`);
                // Depending on strictness, you might throw an error here for NRTL.
            }

            // 2. Get properties (try specific sources first, e.g., 'chemsep1', 'DWSIM')
            // Adjust sources based on your database content
            const sourcesToTry = ['chemsep1', 'chemsep2', 'DWSIM', 'biod_db'];
            let properties: any = null;
            let foundSource: string | null = null;

            for (const source of sourcesToTry) {
                 const { data: propsData, error: propsError } = await supabase
                    .from('compound_properties')
                    .select('properties')
                    .eq('compound_id', compoundId)
                    .eq('source', source)
                    .single(); // Expect only one entry per compound/source pair

                 if (propsError && propsError.code !== 'PGRST116') { // Ignore 'not found' error
                     console.warn(`Supabase properties query error for source ${source}: ${propsError.message}`);
                     // Continue trying other sources
                 } else if (propsData) {
                     properties = propsData.properties;
                     foundSource = source;
                     console.log(`Found properties from source: ${source}`);
                     break; // Stop searching once properties are found
                 }
            }


            if (!properties) {
                // Fallback: Get any available property set if specific sources failed
                 const { data: anyPropsData, error: anyPropsError } = await supabase
                    .from('compound_properties')
                    .select('properties, source')
                    .eq('compound_id', compoundId)
                    .limit(1); // Get the first available one

                 if (anyPropsError) {
                     console.warn(`Supabase fallback properties query error: ${anyPropsError.message}`);
                 } else if (anyPropsData && anyPropsData.length > 0) {
                     properties = anyPropsData[0].properties;
                     foundSource = anyPropsData[0].source;
                     console.log(`Found properties from fallback source: ${foundSource}`);
                 } else {
                    throw new Error(`No properties found for compound ID ${compoundId} ('${foundName}').`);
                 }
            }

             if (typeof properties !== 'object' || properties === null) {
                 throw new Error(`Invalid properties format found for ${foundName} from source ${foundSource}. Expected an object.`);
             }


            // 3. Extract Antoine Parameters (adapt based on common patterns)
            let antoine: AntoineParams | null = null; // Uses imported AntoineParams
            // Try ChemSep style first
            const antoineChemsep = properties.Antoine || properties.AntoineVaporPressure;
            if (antoineChemsep && typeof antoineChemsep === 'object' && antoineChemsep.A && antoineChemsep.B && antoineChemsep.C) {
                antoine = {
                    A: parseFloat(antoineChemsep.A?.value ?? antoineChemsep.A),
                    B: parseFloat(antoineChemsep.B?.value ?? antoineChemsep.B),
                    C: parseFloat(antoineChemsep.C?.value ?? antoineChemsep.C),
                    Tmin_K: parseFloat(antoineChemsep.Tmin?.value ?? antoineChemsep.Tmin ?? 0),
                    Tmax_K: parseFloat(antoineChemsep.Tmax?.value ?? antoineChemsep.Tmax ?? 10000),
                    Units: antoineChemsep.units || 'Pa', // Default to Pa if missing
                    EquationNo: antoineChemsep.eqno
                };
                console.log(`Extracted ChemSep-style Antoine for ${foundName}`);
            }
            // Add extraction logic for other styles (biod_db, DIPPR) if needed here
            // else if (properties.Vapor_Pressure_kPa_a) { ... }
            // else if (properties.DIPPR_Vapor_Pressure_Constant_A) { ... }

            if (!antoine) {
                 console.warn(`Could not extract Antoine parameters for ${foundName} from source ${foundSource}.`);
                 // Allow proceeding without Antoine if only UNIFAC needed? Or throw error?
                 // For VLE, Antoine is essential.
                 throw new Error(`Failed to extract Antoine parameters for ${foundName}.`);
            }
             // Validate extracted Antoine params
             if (isNaN(antoine.A) || isNaN(antoine.B) || isNaN(antoine.C)) {
                 throw new Error(`Invalid Antoine parameters (NaN) extracted for ${foundName}. A=${antoine.A}, B=${antoine.B}, C=${antoine.C}`);
             }


            // 4. Extract UNIFAC Group Composition (only if UNIFAC is potentially used or data is available)
            let unifacGroups: import('@/lib/vle-types').UnifacGroupComposition | null = null; // Explicit import if UnifacGroupComposition is in vle-types
            if (properties.elements_composition?.UNIFAC && typeof properties.elements_composition.UNIFAC === 'object') {
                unifacGroups = {};
                for (const key in properties.elements_composition.UNIFAC) {
                    if (Object.prototype.hasOwnProperty.call(properties.elements_composition.UNIFAC, key)) {
                        const subgroupId = parseInt(key);
                        const count = parseInt(properties.elements_composition.UNIFAC[key]);
                        if (!isNaN(subgroupId) && !isNaN(count) && count > 0) {
                            unifacGroups[subgroupId] = count;
                        }
                    }
                }
                if (Object.keys(unifacGroups).length > 0) {
                    console.log(`Extracted UNIFAC groups for ${foundName}:`, unifacGroups);
                } else {
                    console.warn(`UNIFAC group data found for ${foundName} but was empty or invalid.`);
                    unifacGroups = null;
                }
            } else {
                console.warn(`No UNIFAC group data found in properties.elements_composition.UNIFAC for ${foundName}.`);
            }

             if (fluidPackage === 'unifac' && !unifacGroups) {
                 // Throw error only if UNIFAC is selected and groups are missing
                 throw new Error(`Failed to extract valid UNIFAC groups for ${foundName}, which are required for UNIFAC model.`);
             }

            // 5. Extract Peng-Robinson Parameters (Tc, Pc, omega)
            let prParams: PrPureComponentParams | null = null;
            // Use correct lowercase, space-separated keys as per the database JSON structure
            const tcPropObj = properties["Critical temperature"];
            const pcPropObj = properties["Critical pressure"];
            const omegaPropObj = properties["Acentric factor"];

            if (tcPropObj && pcPropObj && omegaPropObj) {
                const Tc_K_val = (tcPropObj && typeof tcPropObj === 'object' && tcPropObj.value !== undefined) ? parseFloat(tcPropObj.value) : NaN;
                
                const pcValue = (pcPropObj && typeof pcPropObj === 'object' && pcPropObj.value !== undefined) ? parseFloat(pcPropObj.value) : NaN;
                const pcUnits = (pcPropObj && typeof pcPropObj === 'object' && pcPropObj.units !== undefined) ? String(pcPropObj.units).toLowerCase() : '';
                const Pc_Pa_val = pcValue * (pcUnits === 'kpa' ? 1000 : pcUnits === 'bar' ? 100000 : 1);
                
                const omega_val = (omegaPropObj && typeof omegaPropObj === 'object' && omegaPropObj.value !== undefined) ? parseFloat(omegaPropObj.value) : NaN;

                if (!isNaN(Tc_K_val) && !isNaN(Pc_Pa_val) && !isNaN(omega_val)) {
                    prParams = { Tc_K: Tc_K_val, Pc_Pa: Pc_Pa_val, omega: omega_val };
                    console.log(`Extracted Peng-Robinson parameters for ${foundName}: Tc=${Tc_K_val.toFixed(1)}K, Pc=${(Pc_Pa_val/1000).toFixed(1)}kPa, omega=${omega_val.toFixed(3)}`);
                } else {
                    console.warn(`One or more PR parameters (Tc, Pc, omega) are invalid or have unexpected structure for ${foundName}. Tc=${Tc_K_val}, Pc_Pa=${Pc_Pa_val}, omega=${omega_val}`);
                    console.warn("PR Property Objects:", {tcPropObj, pcPropObj, omegaPropObj});
                }
            } else {
                 console.warn(`Peng-Robinson parameters ("Critical temperature", "Critical pressure", "Acentric factor") not found or incomplete in properties for ${foundName} from source ${foundSource}.`);
            }

            if (fluidPackage === 'pr' && !prParams) {
                throw new Error(`Failed to extract Peng-Robinson parameters for ${foundName}, which are required for PR model.`);
            }

            // 6. Extract SRK Parameters (Tc, Pc, omega) - same as PR for now // ADDED
            let srkParams: SrkPureComponentParams | null = null;
            // Assuming SRK uses the same critical properties and acentric factor keys as PR
            if (tcPropObj && pcPropObj && omegaPropObj) { // Re-use already fetched property objects
                const Tc_K_val = (tcPropObj && typeof tcPropObj === 'object' && tcPropObj.value !== undefined) ? parseFloat(tcPropObj.value) : NaN;
                const pcValue = (pcPropObj && typeof pcPropObj === 'object' && pcPropObj.value !== undefined) ? parseFloat(pcPropObj.value) : NaN;
                const pcUnits = (pcPropObj && typeof pcPropObj === 'object' && pcPropObj.units !== undefined) ? String(pcPropObj.units).toLowerCase() : '';
                const Pc_Pa_val = pcValue * (pcUnits === 'kpa' ? 1000 : pcUnits === 'bar' ? 100000 : 1);
                const omega_val = (omegaPropObj && typeof omegaPropObj === 'object' && omegaPropObj.value !== undefined) ? parseFloat(omegaPropObj.value) : NaN;

                if (!isNaN(Tc_K_val) && !isNaN(Pc_Pa_val) && !isNaN(omega_val)) {
                    srkParams = { Tc_K: Tc_K_val, Pc_Pa: Pc_Pa_val, omega: omega_val };
                    console.log(`Extracted SRK parameters for ${foundName}: Tc=${Tc_K_val.toFixed(1)}K, Pc=${(Pc_Pa_val/1000).toFixed(1)}kPa, omega=${omega_val.toFixed(3)}`);
                } else {
                    // Warning already issued for PR, avoid duplicate if keys are identical
                    if (fluidPackage === 'srk') { // Only warn specifically if SRK is selected and PR didn't already warn
                        console.warn(`One or more SRK parameters (Tc, Pc, omega) are invalid or have unexpected structure for ${foundName}.`);
                    }
                }
            } else {
                 if (fluidPackage === 'srk') { // Only warn specifically if SRK is selected
                    console.warn(`SRK parameters ("Critical temperature", "Critical pressure", "Acentric factor") not found or incomplete in properties for ${foundName} from source ${foundSource}.`);
                 }
            }

            if (fluidPackage === 'srk' && !srkParams) {
                throw new Error(`Failed to extract SRK parameters for ${foundName}, which are required for SRK model.`);
            }

            // 7. Extract UNIQUAC Parameters (r, q) // ADDED
            let uniquacParams: UniquacPureComponentParams | null = null;
            // Attempt to find UNIQUAC r and q. Keys might vary.
            // Example keys: "UNIQUAC r", "UNIQUAC q", "Van der Waals volume", "Van der Waals area"
            const rPropObj = properties["UNIQUAC r"] || properties["Van der Waals volume"];
            const qPropObj = properties["UNIQUAC q"] || properties["Van der Waals area"];

            if (rPropObj && qPropObj) {
                const r_val = (rPropObj && typeof rPropObj === 'object' && rPropObj.value !== undefined) ? parseFloat(rPropObj.value) : parseFloat(rPropObj); // Allow direct value
                const q_val = (qPropObj && typeof qPropObj === 'object' && qPropObj.value !== undefined) ? parseFloat(qPropObj.value) : parseFloat(qPropObj); // Allow direct value

                if (!isNaN(r_val) && !isNaN(q_val)) {
                    uniquacParams = { r: r_val, q: q_val };
                    console.log(`Extracted UNIQUAC parameters for ${foundName}: r=${r_val.toFixed(4)}, q=${q_val.toFixed(4)}`);
                } else {
                    console.warn(`One or more UNIQUAC parameters (r, q) are invalid for ${foundName}. r=${r_val}, q=${q_val}`);
                }
            } else {
                if (fluidPackage === 'uniquac') { // Only warn if UNIQUAC is selected
                    console.warn(`UNIQUAC parameters (r, q) not found in properties for ${foundName} from source ${foundSource}.`);
                }
            }
            if (fluidPackage === 'uniquac' && !uniquacParams) {
                throw new Error(`Failed to extract UNIQUAC parameters (r,q) for ${foundName}, which are required for UNIQUAC model.`);
            }

            // 8. Extract Wilson Parameters (V_L) // ADDED
            let wilsonParams: WilsonPureComponentParams | null = null;
            // Attempt to find Liquid Molar Volume. Key might be "Liquid molar volume", "Molar volume", etc.
            // Ensure units are converted to m^3/mol if necessary.
            const vLPropObj = properties["Liquid molar volume"] || properties["Molar volume"] || properties["Wilson volume"]; // ADDED "Wilson volume"
            if (vLPropObj) {
                const vL_val_any_unit = (vLPropObj && typeof vLPropObj === 'object' && vLPropObj.value !== undefined) ? parseFloat(vLPropObj.value) : parseFloat(vLPropObj as any);
                const vL_units = (vLPropObj && typeof vLPropObj === 'object' && vLPropObj.units !== undefined) ? String(vLPropObj.units).toLowerCase() : 'cm3/mol'; // Default assumption

                let vL_m3mol: number | undefined;
                if (!isNaN(vL_val_any_unit)) {
                    if (vL_units === 'cm3/mol' || vL_units === 'cm^3/mol') {
                        vL_m3mol = vL_val_any_unit * 1e-6; // cm^3 to m^3
                    } else if (vL_units === 'm3/mol' || vL_units === 'm^3/mol') {
                        vL_m3mol = vL_val_any_unit;
                    } else if (vL_units === 'm3/kmol' || vL_units === 'm^3/kmol') { // ADDED m3/kmol case
                        vL_m3mol = vL_val_any_unit / 1000; // kmol to mol
                    } else if (vL_units === 'l/mol' || vL_units === 'dm3/mol' || vL_units === 'dm^3/mol') {
                        vL_m3mol = vL_val_any_unit * 1e-3; // L or dm^3 to m^3
                    } else {
                        console.warn(`Wilson: Unknown units for molar volume ('${vL_units}') for ${foundName}. Assuming cm3/mol if value looks reasonable, otherwise parameter will be invalid.`);
                        // Heuristic: if value is small (e.g. < 0.1) it might be m3/mol already. If large (e.g. > 1) it's likely cm3/mol or L/mol.
                        // This is risky, better to have consistent units in DB or more explicit unit handling.
                        // For now, if unit is unknown, we might default to cm3/mol conversion or fail.
                        // Let's be conservative and fail if units are not recognized among common ones.
                         if (vL_val_any_unit > 0.01 && vL_val_any_unit < 1000) { // typical range for cm3/mol
                            vL_m3mol = vL_val_any_unit * 1e-6; // Tentatively assume cm3/mol
                            console.warn(`Wilson: Assumed cm3/mol for V_L of ${foundName} due to unrecognized unit '${vL_units}'. Value: ${vL_val_any_unit}`);
                         } else if (vL_val_any_unit > 1 && vL_val_any_unit < 200000 && vL_units === 'm3/kmol') { // typical range for m3/kmol (0.001 to 0.2 m3/mol)
                            vL_m3mol = vL_val_any_unit / 1000; // Tentatively assume m3/kmol
                            console.warn(`Wilson: Assumed m3/kmol for V_L of ${foundName} due to unrecognized unit '${vL_units}'. Value: ${vL_val_any_unit}`);
                         } else {
                            console.error(`Wilson: Molar volume for ${foundName} has unrecognized unit '${vL_units}' and value ${vL_val_any_unit} is outside typical cm3/mol or m3/kmol range. Cannot proceed with Wilson.`);
                         }
                    }
                }

                if (vL_m3mol !== undefined && !isNaN(vL_m3mol) && vL_m3mol > 0) {
                    wilsonParams = { V_L_m3mol: vL_m3mol };
                    console.log(`Extracted Wilson parameter for ${foundName}: V_L=${vL_m3mol.toExponential(4)} m^3/mol (Original: ${vL_val_any_unit} ${vL_units})`);
                } else {
                    console.warn(`Wilson parameter (V_L) is invalid or has unrecognized units for ${foundName}. V_L_val=${vL_val_any_unit}, units=${vL_units}`);
                }
            } else {
                if (fluidPackage === 'wilson') {
                    console.warn(`Wilson parameter (Liquid molar volume) not found in properties for ${foundName} from source ${foundSource}.`);
                }
            }
            if (fluidPackage === 'wilson' && !wilsonParams) {
                throw new Error(`Failed to extract Wilson parameter (V_L) for ${foundName}, which is required for Wilson model.`);
            }


            return { name: foundName, antoine, unifacGroups, cas_number: casNumber, prParams, srkParams, uniquacParams, wilsonParams }; // Include wilsonParams
        } catch (err: any) {
            console.error(`Error fetching data for ${compoundName}:`, err.message);
            setError(`Failed to fetch data for ${compoundName}: ${err.message}`);
            return null;
        }
    }


    // --- Bisection Solver Helper ---
    function bisectionSolve(
        func: (x: number) => number | null,
        a: number,
        b: number,
        tolerance: number = 1e-5,
        maxIter: number = 50
    ): number | null {
        let fa = func(a);
        let fb = func(b);

        if (fa === null || fb === null) {
            console.warn("Bisection: func evaluation failed at initial endpoints.", { a, fa, b, fb });
            return null;
        }
        if (fa * fb >= 0) {
            if (Math.abs(fa) < tolerance) return a;
            if (Math.abs(fb) < tolerance) return b;
            // console.warn("Bisection: Root not bracketed or multiple roots.", { a, fa, b, fb });
            return null;
        }

        let c = a;
        for (let i = 0; i < maxIter; i++) {
            c = (a + b) / 2;
            const fc = func(c);

            if (fc === null) {
                console.warn(`Bisection: func evaluation failed at midpoint c=${c.toFixed(4)}, iter=${i}. Trying to shrink interval.`);
                // Attempt to recover by checking which half might still be valid
                // This is a simple recovery; more sophisticated logic might be needed
                if (fa !== null && func((a+c)/2) !== null) { // if left sub-interval seems okay
                    b = c; // shrink from right
                } else if (fb !== null && func((c+b)/2) !== null) { // if right sub-interval seems okay
                    a = c; // shrink from left
                } else {
                    console.error("Bisection: Cannot recover from null function evaluation at midpoint.");
                    return null; // Give up if recovery fails
                }
                continue; // Retry with shrunk interval
            }
            if (Math.abs(fc) < tolerance || (b - a) / 2 < tolerance) {
                return c;
            }

            if (fa * fc < 0) {
                b = c;
                // fb = fc; // Not needed for standard bisection
            } else {
                a = c;
                fa = fc; // Update fa because the new 'a' is 'c'
            }
        }
        // console.warn("Bisection: Max iterations reached. Last c=", c, "fc=", func(c));
        const final_fc = func(c);
        return (final_fc !== null && Math.abs(final_fc) < tolerance * 10) ? c : null;
    }


    // --- MODIFIED Calculation Handler ---
    const handleVleDiagramCalculation = async () => {
        setIsLoading(true); setError(null); setResults([]); setLogMessages([]);
        console.log(`\n--- Starting VLE Calculation (${diagramType}, ${fluidPackage.toUpperCase()}) ---`);
        console.log(`Compounds: ${comp1Name}/${comp2Name}`);
        console.log(`Target Pressure for T-calcs: ${pressureKPa} kPa`);
        console.log(`Target Temperature for P-calcs: ${temperatureK} K`);


        if (!comp1Name || !comp2Name) {
            setError("Please enter valid compound names."); setIsLoading(false); return;
        }
        // Validate both pressure and temperature as they are always used now
        if (!pressureKPa || pressureKPa <= 0) {
            setError("Please enter a positive pressure value."); setIsLoading(false); return;
        }
        if (!temperatureK || temperatureK <= 0) {
            setError("Please enter a positive temperature value."); setIsLoading(false); return;
        }
        if (!supabase) { setError("Supabase client is not available."); setIsLoading(false); return; }

        try {
            const [data1, data2] = await Promise.all([fetchCompoundData(comp1Name), fetchCompoundData(comp2Name)]);
            if (!data1 || !data2 || !data1.antoine || !data2.antoine) { // unifacGroups check moved
                setIsLoading(false); return;
            }
            
            // Reset UNIFAC specific pre-calculated values if any
            data1.r_i = undefined; data1.q_i = undefined; 
            data2.r_i = undefined; data2.q_i = undefined; 

            const components: CompoundData[] = [data1, data2];
            let activityParameters: UnifacParameters | NrtlInteractionParams | PrInteractionParams | SrkInteractionParams | UniquacInteractionParams | WilsonInteractionParams; // Updated type

            if (fluidPackage === 'unifac') {
                if (!data1.unifacGroups || !data2.unifacGroups) {
                    throw new Error("UNIFAC groups missing for one or both compounds.");
                }
                const allSubgroupIds = new Set<number>();
                components.forEach(comp => { if (comp.unifacGroups) Object.keys(comp.unifacGroups).forEach(id => allSubgroupIds.add(parseInt(id))); });
                activityParameters = await fetchUnifacInteractionParams(supabase, Array.from(allSubgroupIds));
            } else if (fluidPackage === 'nrtl') { // nrtl
                if (!data1.cas_number || !data2.cas_number) {
                    throw new Error("CAS numbers missing for one or both compounds, required for NRTL.");
                }
                // ADDED: Debug log before calling fetchNrtlParameters
                console.log(`[DEBUG page.tsx] Attempting to call fetchNrtlParameters with: cas1='${data1.cas_number}', cas2='${data2.cas_number}'`);
                activityParameters = await fetchNrtlParameters(supabase, data1.cas_number, data2.cas_number);
            } else if (fluidPackage === 'pr') { 
                if (!data1.prParams || !data2.prParams) {
                     throw new Error("Peng-Robinson parameters (Tc, Pc, omega) missing for one or both compounds.");
                }
                if (!data1.cas_number || !data2.cas_number) {
                    console.warn("CAS numbers missing for one or both compounds; PR k_ij might default to 0 if not found by other means.");
                     activityParameters = { k12: 0 }; 
                } else {
                    activityParameters = await fetchPrInteractionParams(supabase, data1.cas_number, data2.cas_number);
                }
            } else if (fluidPackage === 'srk') { // srk  // CORRECTED
                if (!data1.srkParams || !data2.srkParams) {
                    throw new Error("SRK parameters (Tc, Pc, omega) missing for one or both compounds.");
                }
                if (!data1.cas_number || !data2.cas_number) {
                    console.warn("CAS numbers missing for one or both compounds; SRK k_ij might default to 0.");
                    activityParameters = { k12: 0 };
                } else {
                    activityParameters = await fetchSrkInteractionParams(supabase, data1.cas_number, data2.cas_number);
                }
            } else if (fluidPackage === 'uniquac') { 
                if (!data1.uniquacParams || !data2.uniquacParams) {
                    throw new Error("UNIQUAC parameters (r, q) missing for one or both compounds.");
                }
                if (!data1.cas_number || !data2.cas_number) {
                    console.warn("CAS numbers missing for one or both compounds; UNIQUAC Aij might default to 0.");
                    activityParameters = { A12: 0, A21: 0 };
                } else {
                    activityParameters = await fetchUniquacInteractionParams(supabase, data1.cas_number, data2.cas_number);
                }
            } else if (fluidPackage === 'wilson') { // ADDED WILSON CASE
                if (!data1.wilsonParams || !data2.wilsonParams) {
                    throw new Error("Wilson parameters (V_L) missing for one or both compounds.");
                }
                if (!data1.cas_number || !data2.cas_number) {
                    console.warn("CAS numbers missing for one or both compounds; Wilson a_ij might default to 0.");
                    activityParameters = { a12_J_mol: 0, a21_J_mol: 0 };
                } else {
                    activityParameters = await fetchWilsonInteractionParams(supabase, data1.cas_number, data2.cas_number);
                }
            } else {
                // This case should ideally not be reached if fluidPackage is correctly typed and handled
                throw new Error(`Unsupported fluid package selected: ${fluidPackage}`);
            }

            const targetPressurePa = pressureKPa * 1000; 
            const targetTemperatureK = temperatureK;   

            const x_feed_values = Array.from({ length: 21 }, (_, i) => parseFloat((i * 0.05).toFixed(3))); 
            const calculatedResults: BubbleDewResult[] = [];

            let T_bp1_est_at_targetP: number | undefined;
            let T_bp2_est_at_targetP: number | undefined;
            let P_sat1_est_at_targetT: number | undefined;
            let P_sat2_est_at_targetT: number | undefined;

            if (data1.antoine) T_bp1_est_at_targetP = antoineBoilingPointSolver(data1.antoine, targetPressurePa) ?? undefined;
            if (data2.antoine) T_bp2_est_at_targetP = antoineBoilingPointSolver(data2.antoine, targetPressurePa) ?? undefined;
            if (data1.antoine) P_sat1_est_at_targetT = libCalculatePsat_Pa(data1.antoine, targetTemperatureK);
            if (data2.antoine) P_sat2_est_at_targetT = libCalculatePsat_Pa(data2.antoine, targetTemperatureK);

            for (const x1_feed of x_feed_values) {
                console.log(`\nCalculating for feed mole fraction ${comp1Name} (x₁) = ${x1_feed.toFixed(3)} using ${fluidPackage.toUpperCase()}...`);
                let bubbleT_resultPoint: BubbleDewResult | null = null;
                let bubbleP_resultPoint: BubbleDewResult | null = null;

                // --- Calculate Bubble Temperature at targetPressurePa ---
                const initialTempGuessForBubbleT = (T_bp1_est_at_targetP && T_bp2_est_at_targetP 
                    ? (x1_feed * T_bp1_est_at_targetP + (1 - x1_feed) * T_bp2_est_at_targetP) 
                    : targetTemperatureK) || 350;
                
                if (fluidPackage === 'unifac') {
                    bubbleT_resultPoint = calculateBubbleTemperatureUnifac(components, x1_feed, targetPressurePa, activityParameters as UnifacParameters, initialTempGuessForBubbleT);
                } else if (fluidPackage === 'nrtl') { // nrtl
                    bubbleT_resultPoint = calculateBubbleTemperatureNrtl(components, x1_feed, targetPressurePa, activityParameters as NrtlInteractionParams, initialTempGuessForBubbleT);
                } else if (fluidPackage === 'pr') { 
                    bubbleT_resultPoint = calculateBubbleTemperaturePr(components, x1_feed, targetPressurePa, activityParameters as PrInteractionParams, initialTempGuessForBubbleT);
                } else if (fluidPackage === 'srk') { 
                    bubbleT_resultPoint = calculateBubbleTemperatureSrk(components, x1_feed, targetPressurePa, activityParameters as SrkInteractionParams, initialTempGuessForBubbleT);
                } else if (fluidPackage === 'uniquac') { 
                    bubbleT_resultPoint = calculateBubbleTemperatureUniquac(components, x1_feed, targetPressurePa, activityParameters as UniquacInteractionParams, initialTempGuessForBubbleT);
                } else if (fluidPackage === 'wilson') { // ADDED WILSON CASE
                    bubbleT_resultPoint = calculateBubbleTemperatureWilson(components, x1_feed, targetPressurePa, activityParameters as WilsonInteractionParams, initialTempGuessForBubbleT);
                }

                if (bubbleT_resultPoint && !bubbleT_resultPoint.error) {
                    calculatedResults.push(bubbleT_resultPoint);
                } else {
                    console.error(`   Bubble Temperature calculation (${fluidPackage.toUpperCase()}) failed for x₁ = ${x1_feed.toFixed(3)}. Error: ${bubbleT_resultPoint?.error || 'Unknown error'}`);
                    // Ensure T_K is present, even if NaN, for bubbleT type errors
                    if (bubbleT_resultPoint) {
                        calculatedResults.push({...bubbleT_resultPoint, T_K: bubbleT_resultPoint.T_K ?? NaN });
                    } else {
                        calculatedResults.push({comp1_feed: x1_feed, comp1_equilibrium: NaN, T_K: NaN, error: "BubbleT calc object null", calculationType: 'bubbleT', P_Pa: targetPressurePa});
                    }
                }

                // --- Calculate Bubble Pressure at targetTemperatureK ---
                const initialPressureGuessForBubbleP = (P_sat1_est_at_targetT && P_sat2_est_at_targetT 
                    ? (x1_feed * P_sat1_est_at_targetT + (1 - x1_feed) * P_sat2_est_at_targetT) 
                    : targetPressurePa) || 101325;

                if (fluidPackage === 'unifac') {
                    bubbleP_resultPoint = calculateBubblePressureUnifac(components, x1_feed, targetTemperatureK, activityParameters as UnifacParameters, initialPressureGuessForBubbleP);
                } else if (fluidPackage === 'nrtl') { // nrtl
                    bubbleP_resultPoint = calculateBubblePressureNrtl(components, x1_feed, targetTemperatureK, activityParameters as NrtlInteractionParams, initialPressureGuessForBubbleP);
                } else if (fluidPackage === 'pr') { 
                     bubbleP_resultPoint = calculateBubblePressurePr(components, x1_feed, targetTemperatureK, activityParameters as PrInteractionParams, initialPressureGuessForBubbleP);
                } else if (fluidPackage === 'srk') { 
                    bubbleP_resultPoint = calculateBubblePressureSrk(components, x1_feed, targetTemperatureK, activityParameters as SrkInteractionParams, initialPressureGuessForBubbleP);
                } else if (fluidPackage === 'uniquac') { 
                    bubbleP_resultPoint = calculateBubblePressureUniquac(components, x1_feed, targetTemperatureK, activityParameters as UniquacInteractionParams, initialPressureGuessForBubbleP);
                } else if (fluidPackage === 'wilson') { // ADDED WILSON CASE
                    bubbleP_resultPoint = calculateBubblePressureWilson(components, x1_feed, targetTemperatureK, activityParameters as WilsonInteractionParams, initialPressureGuessForBubbleP);
                }
                
                if (bubbleP_resultPoint && !bubbleP_resultPoint.error) {
                    calculatedResults.push(bubbleP_resultPoint);
                } else {
                    console.error(`   Bubble Pressure calculation (${fluidPackage.toUpperCase()}) failed for x₁ = ${x1_feed.toFixed(3)}. Error: ${bubbleP_resultPoint?.error || 'Unknown error'}`);
                    // Ensure P_Pa is present, even if NaN, for bubbleP type errors
                    if (bubbleP_resultPoint) {
                        calculatedResults.push({...bubbleP_resultPoint, P_Pa: bubbleP_resultPoint.P_Pa ?? NaN});
                    } else {
                        calculatedResults.push({comp1_feed: x1_feed, comp1_equilibrium: NaN, P_Pa: NaN, error: "BubbleP calc object null", calculationType: 'bubbleP', T_K: targetTemperatureK});
                    }
                }
            }
            
            // Add pure component points
            // For Bubble T type results (at targetPressurePa)
            if (T_bp2_est_at_targetP) calculatedResults.push({comp1_feed: 0.0, comp1_equilibrium: 0.0, T_K: T_bp2_est_at_targetP, P_Pa: targetPressurePa, iterations:0, calculationType: 'bubbleT'});
            if (T_bp1_est_at_targetP) calculatedResults.push({comp1_feed: 1.0, comp1_equilibrium: 1.0, T_K: T_bp1_est_at_targetP, P_Pa: targetPressurePa, iterations:0, calculationType: 'bubbleT'});
            
            // For Bubble P type results (at targetTemperatureK)
            if (P_sat2_est_at_targetT) calculatedResults.push({comp1_feed: 0.0, comp1_equilibrium: 0.0, P_Pa: P_sat2_est_at_targetT, T_K: targetTemperatureK, iterations:0, calculationType: 'bubbleP'});
            if (P_sat1_est_at_targetT) calculatedResults.push({comp1_feed: 1.0, comp1_equilibrium: 1.0, P_Pa: P_sat1_est_at_targetT, T_K: targetTemperatureK, iterations:0, calculationType: 'bubbleP'});


            calculatedResults.sort((a, b) => { // Sort primarily by feed, then by calculation type (optional, but good for consistency)
                if (a.comp1_feed !== b.comp1_feed) return a.comp1_feed - b.comp1_feed;
                if (a.calculationType && b.calculationType) return a.calculationType.localeCompare(b.calculationType);
                return 0;
            });
            setResults(calculatedResults);
            console.log(`\n--- Calculation Complete (Both T and P variations, ${fluidPackage.toUpperCase()}) ---`);
        } catch (err: any) {
            let errorMessage = "An unknown error occurred during VLE diagram calculation.";
            if (err instanceof Error) {
                errorMessage = err.message;
            } else if (typeof err === 'string') {
                errorMessage = err;
            } else {
                try {
                    errorMessage = JSON.stringify(err);
                } catch (stringifyError) {
                    errorMessage = "An un-stringifiable error object was thrown.";
                }
            }
            console.error("VLE Diagram Calculation failed. Raw error object:", err); // Log the raw error
            console.error("Processed error message for VLE Diagram:", errorMessage);
            setError(`VLE Diagram Calculation failed: ${errorMessage}`);
        } finally {
            setIsLoading(false);
        }
    };

    const handleAzeotropeScan = async () => {
        setIsLoading(true); setError(null); setAzeotropeScanData([]); setLogMessages([]);
        console.log(`\n--- Starting Azeotrope Scan (${azeotropeScanType}, ${fluidPackage.toUpperCase()}) ---`);
        console.log(`Compounds: ${comp1Name}/${comp2Name}`);
        // console.log(`Scan Range: ${scanStartValue} to ${scanEndValue}, Steps: ${scanSteps}`); // Commented out as range is now internal

        if (!comp1Name || !comp2Name) {
            setError("Please enter valid compound names."); setIsLoading(false); return;
        }
        if (!supabase) { setError("Supabase client is not available."); setIsLoading(false); return; }

        try {
            const [data1, data2] = await Promise.all([fetchCompoundData(comp1Name), fetchCompoundData(comp2Name)]);
            if (!data1 || !data2 || !data1.antoine || !data2.antoine) {
                setIsLoading(false); return;
            }
            data1.r_i = undefined; data2.r_i = undefined; // Reset for UNIFAC
            const components: CompoundData[] = [data1, data2];
            let activityParameters: UnifacParameters | NrtlInteractionParams | PrInteractionParams | SrkInteractionParams | UniquacInteractionParams | WilsonInteractionParams; // Updated type

            if (fluidPackage === 'unifac') {
                if (!data1.unifacGroups || !data2.unifacGroups) {
                    throw new Error("UNIFAC groups missing for one or both compounds.");
                }
                const allSubgroupIds = new Set<number>();
                components.forEach(comp => { if (comp.unifacGroups) Object.keys(comp.unifacGroups).forEach(id => allSubgroupIds.add(parseInt(id))); });
                activityParameters = await fetchUnifacInteractionParams(supabase, Array.from(allSubgroupIds));
            } else if (fluidPackage === 'nrtl') { // nrtl
                if (!data1.cas_number || !data2.cas_number) {
                    throw new Error("CAS numbers missing for one or both compounds, required for NRTL.");
                }
                // ADDED: Debug log before calling fetchNrtlParameters
                console.log(`[DEBUG page.tsx] Attempting to call fetchNrtlParameters (Azeotrope Scan) with: cas1='${data1.cas_number}', cas2='${data2.cas_number}'`);
                activityParameters = await fetchNrtlParameters(supabase, data1.cas_number, data2.cas_number);
            } else if (fluidPackage === 'pr') { 
                 if (!data1.prParams || !data2.prParams) {
                     throw new Error("Peng-Robinson parameters (Tc, Pc, omega) missing for one or both compounds for Azeotrope scan.");
                }
                 if (!data1.cas_number || !data2.cas_number) {
                    console.warn("CAS numbers missing for Azeotrope scan with PR; k_ij might default to 0.");
                    activityParameters = { k12: 0 };
                } else {
                    activityParameters = await fetchPrInteractionParams(supabase, data1.cas_number, data2.cas_number);
                }
            } else if (fluidPackage === 'srk') { // srk // CORRECTED
                if (!data1.srkParams || !data2.srkParams) {
                    throw new Error("SRK parameters (Tc, Pc, omega) missing for one or both compounds for Azeotrope scan.");
                }
                if (!data1.cas_number || !data2.cas_number) {
                   console.warn("CAS numbers missing for Azeotrope scan with SRK; k_ij might default to 0.");
                   activityParameters = { k12: 0 };
               } else {
                   activityParameters = await fetchSrkInteractionParams(supabase, data1.cas_number, data2.cas_number);
               }
           } else if (fluidPackage === 'uniquac') { // uniquac // CORRECTED
                if (!data1.uniquacParams || !data2.uniquacParams) {
                    throw new Error("UNIQUAC parameters (r, q) missing for one or both compounds for Azeotrope scan.");
                }
                if (!data1.cas_number || !data2.cas_number) {
                    console.warn("CAS numbers missing for Azeotrope scan with UNIQUAC; Aij might default to 0.");
                    activityParameters = { A12: 0, A21: 0 };
                } else {
                    activityParameters = await fetchUniquacInteractionParams(supabase, data1.cas_number, data2.cas_number);
                }
           } else if (fluidPackage === 'wilson') { // ADDED WILSON CASE
                if (!data1.wilsonParams || !data2.wilsonParams) {
                    throw new Error("Wilson parameters (V_L) missing for one or both compounds for Azeotrope scan.");
                }
                if (!data1.cas_number || !data2.cas_number) {
                    console.warn("CAS numbers missing for Azeotrope scan with Wilson; a_ij might default to 0.");
                    activityParameters = { a12_J_mol: 0, a21_J_mol: 0 };
                } else {
                    activityParameters = await fetchWilsonInteractionParams(supabase, data1.cas_number, data2.cas_number);
                }
           } else {
                throw new Error(`Unsupported fluid package selected for azeotrope scan: ${fluidPackage}`);
           }

            const scanResultsArray: { scanVal: number, x_az: number, dependentVal: number }[] = [];
            
            // --- Define internal wide scan ranges ---
            let currentScanStartValue: number;
            let currentScanEndValue: number;

            if (azeotropeScanType === 'vs_P_find_T') {
                currentScanStartValue = 1; // Default P_start (kPa) - very low
                currentScanEndValue = 2000;  // Default P_end (kPa) - quite high
                console.log(`Automatic Scan Range (Pressure): ${currentScanStartValue} kPa to ${currentScanEndValue} kPa, Steps: ${scanSteps}`);
            } else { // vs_T_find_P
                currentScanStartValue = 230; // Default T_start (K) - below typical ambient
                currentScanEndValue = 550;   // Default T_end (K) - towards higher boiling points / critical region
                console.log(`Automatic Scan Range (Temperature): ${currentScanStartValue} K to ${currentScanEndValue} K, Steps: ${scanSteps}`);
            }

            const stepValue = (currentScanEndValue - currentScanStartValue) / (scanSteps > 1 ? scanSteps - 1 : 1);

            for (let i = 0; i < scanSteps; i++) {
                const currentScanVal = currentScanStartValue + i * stepValue;
                let x_az_found: number | null = null;
                let dependent_val_found: number | null = null;

                if (azeotropeScanType === 'vs_P_find_T') {
                    const P_system_Pa = currentScanVal * 1000;
                    console.log(`Scanning for azeotrope at P = ${currentScanVal.toFixed(2)} kPa using ${fluidPackage.toUpperCase()}...`);

                    const objective_h = (x1_guess: number): number | null => {
                        if (x1_guess <= 1e-4 || x1_guess >= 1.0 - 1e-4) return 1.0;
                        const T_bp1_est_at_P = antoineBoilingPointSolver(components[0].antoine, P_system_Pa) || 350;
                        const T_bp2_est_at_P = antoineBoilingPointSolver(components[1].antoine, P_system_Pa) || 350;
                        const initialTempGuess = x1_guess * T_bp1_est_at_P + (1 - x1_guess) * T_bp2_est_at_P;
                        
                        let bubbleT_result: BubbleDewResult | null = null;
                        if (fluidPackage === 'unifac') {
                            bubbleT_result = calculateBubbleTemperatureUnifac(components, x1_guess, P_system_Pa, activityParameters as UnifacParameters, initialTempGuess, 30, 1e-4);
                        } else if (fluidPackage === 'nrtl') {
                            bubbleT_result = calculateBubbleTemperatureNrtl(components, x1_guess, P_system_Pa, activityParameters as NrtlInteractionParams, initialTempGuess, 30, 1e-4);
                        } else if (fluidPackage === 'pr') { 
                            bubbleT_result = calculateBubbleTemperaturePr(components, x1_guess, P_system_Pa, activityParameters as PrInteractionParams, initialTempGuess, 30, 1e-4);
                        } else if (fluidPackage === 'srk') { // CORRECTED
                            bubbleT_result = calculateBubbleTemperatureSrk(components, x1_guess, P_system_Pa, activityParameters as SrkInteractionParams, initialTempGuess, 30, 1e-4);
                        } else if (fluidPackage === 'uniquac') { // CORRECTED
                            bubbleT_result = calculateBubbleTemperatureUniquac(components, x1_guess, P_system_Pa, activityParameters as UniquacInteractionParams, initialTempGuess, 30, 1e-4);
                        } else if (fluidPackage === 'wilson') { // ADDED WILSON CASE
                            bubbleT_result = calculateBubbleTemperatureWilson(components, x1_guess, P_system_Pa, activityParameters as WilsonInteractionParams, initialTempGuess, 30, 1e-4);
                        }

                        if (bubbleT_result && bubbleT_result.error === undefined && typeof bubbleT_result.comp1_equilibrium === 'number' && !isNaN(bubbleT_result.comp1_equilibrium)) {
                            return x1_guess - bubbleT_result.comp1_equilibrium;
                        }
                        return null;
                    };
                    x_az_found = bisectionSolve(objective_h, 0.001, 0.999, 1e-4, 30);

                    if (x_az_found !== null) {
                        const T_bp1_est_final = antoineBoilingPointSolver(components[0].antoine, P_system_Pa) || 350;
                        const T_bp2_est_final = antoineBoilingPointSolver(components[1].antoine, P_system_Pa) || 350;
                        const initialTempGuessFinal = x_az_found * T_bp1_est_final + (1 - x_az_found) * T_bp2_est_final;
                        let final_bubble_res: BubbleDewResult | null = null;
                        if (fluidPackage === 'unifac') {
                            final_bubble_res = calculateBubbleTemperatureUnifac(components, x_az_found, P_system_Pa, activityParameters as UnifacParameters, initialTempGuessFinal, 30, 1e-4);
                        } else if (fluidPackage === 'nrtl') {
                            final_bubble_res = calculateBubbleTemperatureNrtl(components, x_az_found, P_system_Pa, activityParameters as NrtlInteractionParams, initialTempGuessFinal, 30, 1e-4);
                        } else if (fluidPackage === 'pr') { 
                            final_bubble_res = calculateBubbleTemperaturePr(components, x_az_found, P_system_Pa, activityParameters as PrInteractionParams, initialTempGuessFinal, 30, 1e-4);
                        } else if (fluidPackage === 'srk') { // CORRECTED
                            final_bubble_res = calculateBubbleTemperatureSrk(components, x_az_found, P_system_Pa, activityParameters as SrkInteractionParams, initialTempGuessFinal, 30, 1e-4);
                        } else if (fluidPackage === 'uniquac') { // CORRECTED
                            final_bubble_res = calculateBubbleTemperatureUniquac(components, x_az_found, P_system_Pa, activityParameters as UniquacInteractionParams, initialTempGuessFinal, 30, 1e-4);
                        } else if (fluidPackage === 'wilson') { // ADDED WILSON CASE
                            final_bubble_res = calculateBubbleTemperatureWilson(components, x_az_found, P_system_Pa, activityParameters as WilsonInteractionParams, initialTempGuessFinal, 30, 1e-4);
                        }
                        if (final_bubble_res && final_bubble_res.error === undefined && final_bubble_res.T_K) {
                            dependent_val_found = final_bubble_res.T_K;
                            console.log(`  Found: P=${currentScanVal.toFixed(2)} kPa -> x_az=${x_az_found.toFixed(4)}, T_az=${dependent_val_found.toFixed(2)} K`);
                        } else { x_az_found = null; }
                    }
                } else { // vs_T_find_P
                    const T_system_K = currentScanVal;
                    console.log(`Scanning for azeotrope at T = ${T_system_K.toFixed(2)} K using ${fluidPackage.toUpperCase()}...`);
                    const objective_k = (x1_guess: number): number | null => {
                        if (x1_guess <= 1e-4 || x1_guess >= 1.0 - 1e-4) return 1.0;
                        const P_sat1_est = libCalculatePsat_Pa(components[0].antoine!, T_system_K);
                        const P_sat2_est = libCalculatePsat_Pa(components[1].antoine!, T_system_K);
                        if (isNaN(P_sat1_est) || isNaN(P_sat2_est)) return null;
                        const initialPressureGuess = x1_guess * P_sat1_est + (1 - x1_guess) * P_sat2_est;

                        let bubbleP_result: BubbleDewResult | null = null;
                        if (fluidPackage === 'unifac') {
                            bubbleP_result = calculateBubblePressureUnifac(components, x1_guess, T_system_K, activityParameters as UnifacParameters, initialPressureGuess || 101325, 10, 1e-5);
                        } else if (fluidPackage === 'nrtl') {
                            bubbleP_result = calculateBubblePressureNrtl(components, x1_guess, T_system_K, activityParameters as NrtlInteractionParams, initialPressureGuess || 101325, 10, 1e-5);
                        } else if (fluidPackage === 'pr') { 
                            bubbleP_result = calculateBubblePressurePr(components, x1_guess, T_system_K, activityParameters as PrInteractionParams, initialPressureGuess || 101325, 10, 1e-5);
                        } else if (fluidPackage === 'srk') { // CORRECTED
                            bubbleP_result = calculateBubblePressureSrk(components, x1_guess, T_system_K, activityParameters as SrkInteractionParams, initialPressureGuess || 101325, 10, 1e-5);
                        } else if (fluidPackage === 'uniquac') { 
                            bubbleP_result = calculateBubblePressureUniquac(components, x1_guess, T_system_K, activityParameters as UniquacInteractionParams, initialPressureGuess || 101325, 10, 1e-5);
                        } else if (fluidPackage === 'wilson') { // ADDED WILSON CASE
                            bubbleP_result = calculateBubblePressureWilson(components, x1_guess, T_system_K, activityParameters as WilsonInteractionParams, initialPressureGuess || 101325, 10, 1e-5);
                        }
                        
                        if (bubbleP_result && bubbleP_result.error === undefined && typeof bubbleP_result.comp1_equilibrium === 'number' && !isNaN(bubbleP_result.comp1_equilibrium)) {
                            return x1_guess - bubbleP_result.comp1_equilibrium;
                        }
                        return null;
                    };
                    x_az_found = bisectionSolve(objective_k, 0.001, 0.999, 1e-4, 30);

                    if (x_az_found !== null) {
                        const P_sat1_est_final = libCalculatePsat_Pa(components[0].antoine!, T_system_K);
                        const P_sat2_est_final = libCalculatePsat_Pa(components[1].antoine!, T_system_K);
                        if (isNaN(P_sat1_est_final) || isNaN(P_sat2_est_final)) { x_az_found = null; }
                        else {
                            const initialPressureGuessFinal = x_az_found * P_sat1_est_final + (1 - x_az_found) * P_sat2_est_final;
                            let final_bubble_res: BubbleDewResult | null = null;
                            if (fluidPackage === 'unifac') {
                                final_bubble_res = calculateBubblePressureUnifac(components, x_az_found, T_system_K, activityParameters as UnifacParameters, initialPressureGuessFinal || 101325, 10, 1e-5);
                            } else if (fluidPackage === 'nrtl') {
                                final_bubble_res = calculateBubblePressureNrtl(components, x_az_found, T_system_K, activityParameters as NrtlInteractionParams, initialPressureGuessFinal || 101325, 10, 1e-5);
                            } else if (fluidPackage === 'pr') { 
                                final_bubble_res = calculateBubblePressurePr(components, x_az_found, T_system_K, activityParameters as PrInteractionParams, initialPressureGuessFinal || 101325, 10, 1e-5);
                            } else if (fluidPackage === 'srk') { // CORRECTED
                                final_bubble_res = calculateBubblePressureSrk(components, x_az_found, T_system_K, activityParameters as SrkInteractionParams, initialPressureGuessFinal || 101325, 10, 1e-5);
                            } else if (fluidPackage === 'uniquac') { 
                                final_bubble_res = calculateBubblePressureUniquac(components, x_az_found, T_system_K, activityParameters as UniquacInteractionParams, initialPressureGuessFinal || 101325, 10, 1e-5);
                            } else if (fluidPackage === 'wilson') { // ADDED WILSON CASE
                                final_bubble_res = calculateBubblePressureWilson(components, x_az_found, T_system_K, activityParameters as WilsonInteractionParams, initialPressureGuessFinal || 101325, 10, 1e-5);
                            }
                            if (final_bubble_res && final_bubble_res.error === undefined && final_bubble_res.P_Pa) {
                                dependent_val_found = final_bubble_res.P_Pa / 1000;
                                console.log(`  Found: T=${T_system_K.toFixed(2)} K -> x_az=${x_az_found.toFixed(4)}, P_az=${dependent_val_found.toFixed(2)} kPa`);
                            } else { x_az_found = null; }
                        }
                    }
                }

                if (x_az_found !== null && dependent_val_found !== null) {
                    scanResultsArray.push({ scanVal: currentScanVal, x_az: x_az_found, dependentVal: dependent_val_found });
                } else {
                    console.log(`  No azeotrope found for scan value: ${currentScanVal.toFixed(2)}`);
                }
            }
            setAzeotropeScanData(scanResultsArray);
            console.log(`\n--- Azeotrope Scan Complete (${fluidPackage.toUpperCase()}) ---`);
        } catch (err: any) {
            let errorMessage = "An unknown error occurred during Azeotrope scan.";
            if (err instanceof Error) {
                errorMessage = err.message;
            } else if (typeof err === 'string') {
                errorMessage = err;
            } else {
                try {
                    errorMessage = JSON.stringify(err);
                } catch (stringifyError) {
                    errorMessage = "An un-stringiable error object was thrown during azeotrope scan.";
                }
            }
            console.error("Azeotrope Scan failed. Raw error object:", err); // Log the raw error
            console.error("Processed error message for Azeotrope Scan:", errorMessage);
            setError(`Azeotrope Scan failed: ${errorMessage}`);
        } finally {
            setIsLoading(false);
        }
    };
    
    const handleCalculate = async () => {
        if (azeotropeFinderActive) {
            await handleAzeotropeScan();
        } else {
            await handleVleDiagramCalculation();
        }
    };


    // --- MODIFIED useEffect hook for ECharts Plotting ---
    useEffect(() => {
        const isDark = resolvedTheme === 'dark';
        const textColor = isDark ? '#E5E7EB' : '#1F2937';
        const lineColor = isDark ? '#4B5563' : '#9CA3AF';
        const liquidColor = isDark ? '#60A5FA' : '#3B82F6'; // For x_az line
        const dependentParamColor = isDark ? '#FBBF24' : '#F59E0B'; // For T_az/P_az line
        const chartFontFamily = '"Merriweather Sans", sans-serif'; // Define font family

        if (azeotropeFinderActive) {
            if (azeotropeScanData.length > 0) {
                const scanParamName = azeotropeScanType === 'vs_P_find_T' ? "Pressure (kPa)" : "Temperature (K)";
                const dependentParamName = azeotropeScanType === 'vs_P_find_T' ? "Azeotropic Temperature (K)" : "Azeotropic Pressure (kPa)";
                const titleText = `Azeotrope Locus for ${comp1Name}/${comp2Name}`;

                const x_az_data = azeotropeScanData.map(d => [d.scanVal, d.x_az]);
                const dependent_val_data = azeotropeScanData.map(d => [d.scanVal, d.dependentVal]);

                setVleChartOptions({
                    title: { text: titleText, left: 'center', textStyle: { color: textColor, fontSize: 16, fontFamily: chartFontFamily } },
                    tooltip: { 
                        trigger: 'axis', 
                        axisPointer: { type: 'cross' },
                        textStyle: { fontFamily: chartFontFamily },
                        formatter: (params: any) => {
                            if (!params || params.length === 0) return '';
                            let tooltipString = `${params[0].axisValueLabel} ${scanParamName.split(' ')[0] === 'Pressure' ? 'kPa' : 'K'}<br/>`;
                            
                            const compParam = params.find((p: any) => p.seriesName === 'Azeotropic Composition (x₁)');
                            const depParam = params.find((p: any) => p.seriesName === dependentParamName);

                            if (compParam) {
                                tooltipString += `${compParam.marker} ${compParam.seriesName}: ${parseFloat(compParam.value[1]).toFixed(3)}<br/>`;
                            }
                            if (depParam) {
                                tooltipString += `${depParam.marker} ${depParam.seriesName}: ${Math.round(parseFloat(depParam.value[1]))} ${dependentParamName.includes('(K)') ? 'K' : 'kPa'}`;
                            }
                            return tooltipString;
                        }
                    },
                    legend: { data: ['Azeotropic Composition (x₁)', dependentParamName], bottom: 5, textStyle: { color: textColor, fontFamily: chartFontFamily } },
                    grid: { left: '10%', right: '12%', bottom: '15%', top: '15%', containLabel: true }, // Adjusted for dual y-axis
                    xAxis: {
                        type: 'value', name: scanParamName, nameLocation: 'middle', nameGap: 30,
                        axisLabel: { color: textColor, fontSize: 12, fontFamily: chartFontFamily }, 
                        nameTextStyle: { color: textColor, fontSize: 14, fontWeight: 'bold', fontFamily: chartFontFamily },
                        axisLine: { lineStyle: { color: lineColor } }, 
                        splitLine: { show: false },
                        scale: true,
                    },
                    yAxis: [
                        {
                            type: 'value', name: 'Azeotropic Composition (x₁)', nameLocation: 'middle', nameGap: 55, min: 0, max: 1,
                            axisLabel: { color: textColor, fontSize: 12, formatter: (v:number) => v.toFixed(3), fontFamily: chartFontFamily },
                            nameTextStyle: { color: textColor, fontSize: 14, fontWeight: 'bold', fontFamily: chartFontFamily },
                            axisLine: { lineStyle: { color: liquidColor } }, 
                            splitLine: { show: false },
                        },
                        {
                            type: 'value', name: dependentParamName, nameLocation: 'middle', nameGap: 70, position: 'right',
                            axisLabel: { color: textColor, fontSize: 12, formatter: (v:number) => Math.round(v).toString(), fontFamily: chartFontFamily }, 
                            nameTextStyle: { color: textColor, fontSize: 14, fontWeight: 'bold', fontFamily: chartFontFamily },
                            axisLine: { lineStyle: { color: dependentParamColor } }, 
                            splitLine: { show: false },
                            scale: true,
                        }
                    ],
                    series: [
                        { name: 'Azeotropic Composition (x₁)', type: 'line', yAxisIndex: 0, symbol: 'circle', symbolSize: 6, lineStyle: { color: liquidColor, width: 2 }, data: x_az_data },
                        { name: dependentParamName, type: 'line', yAxisIndex: 1, symbol: 'triangle', symbolSize: 6, lineStyle: { color: dependentParamColor, width: 2 }, data: dependent_val_data }
                    ]
                });
            } else { // Azeotrope finder active, but no data
                setVleChartOptions({
                    title: {
                        text: `Azeotrope Scan for ${comp1Name}/${comp2Name}`,
                        subtext: 'No azeotropes found within the scanned range.',
                        left: 'center',
                        top: 'center', // Center vertically as well
                        textStyle: { color: textColor, fontSize: 18, fontWeight: 'bold', fontFamily: chartFontFamily },
                        subtextStyle: { color: textColor, fontSize: 14, fontFamily: chartFontFamily }
                    },
                    xAxis: { show: false }, // Hide axes
                    yAxis: { show: false },
                    series: [] // No data series
                });
            }
        } else if (!azeotropeFinderActive && results.length > 0 && comp1Name && comp2Name) {
            // --- Existing VLE Diagram Plotting Logic ---
            let titleText = "";
            let xAxisName = `Mole Fraction ${comp1Name} (x₁, y₁)`;
            let yAxisName = "";
            let series: LineSeriesOption[] = []; 

            // Filter results based on diagramType
            let relevantResults: BubbleDewResult[];
            if (diagramType === "Txy" || diagramType === "xy_constP") {
                relevantResults = results.filter(p => p.calculationType === 'bubbleT' && p.T_K !== undefined);
                yAxisName = "Temperature (K)";
                if (diagramType === "Txy") titleText = `Txy Diagram: ${comp1Name}/${comp2Name} at ${pressureKPa.toFixed(1)} kPa`;
                else titleText = `xy Diagram (T-basis): ${comp1Name}/${comp2Name} at ${pressureKPa.toFixed(1)} kPa`;
            } else { // Pxy or xy_constT
                relevantResults = results.filter(p => p.calculationType === 'bubbleP' && p.P_Pa !== undefined);
                yAxisName = "Pressure (Pa)";
                if (diagramType === "Pxy") titleText = `Pxy Diagram: ${comp1Name}/${comp2Name} at ${temperatureK.toFixed(1)} K`;
                else titleText = `xy Diagram (P-basis): ${comp1Name}/${comp2Name} at ${temperatureK.toFixed(1)} K`;
            }
            
            // Sort relevant Results by feed composition for plotting
            relevantResults.sort((a,b) => a.comp1_feed - b.comp1_feed);


            const liquidLineData: [number, number][] = [];
            const vaporLineDataRaw: [number, number][] = []; 

            if (diagramType === "Txy" || diagramType === "xy_constP") {
                relevantResults.forEach(p => {
                    if (p.T_K !== undefined && p.comp1_equilibrium !== undefined && !isNaN(p.comp1_equilibrium)) {
                        liquidLineData.push([p.comp1_feed, p.T_K]);
                        vaporLineDataRaw.push([p.comp1_equilibrium, p.T_K]);
                    }
                });
            } else { // Pxy or xy_constT
                 relevantResults.forEach(p => {
                    if (p.P_Pa !== undefined && p.comp1_equilibrium !== undefined && !isNaN(p.comp1_equilibrium)) {
                        liquidLineData.push([p.comp1_feed, p.P_Pa]);
                        vaporLineDataRaw.push([p.comp1_equilibrium, p.P_Pa]);
                    }
                });
            }
            
            // liquidLineData is already sorted by comp1_feed due to relevantResults sort
            const vaporLineData = [...vaporLineDataRaw].sort((a,b) => a[0] - b[0]);


            const isDark = resolvedTheme === 'dark';
            const textColor = isDark ? '#E5E7EB' : '#1F2937';
            const lineColor = isDark ? '#4B5563' : '#9CA3AF';
            const liquidColor = isDark ? '#60A5FA' : '#3B82F6';
            const vaporColor = isDark ? '#F87171' : '#EF4444';
            const chartFontFamily = '"Merriweather Sans", sans-serif'; // Define font family for VLE diagrams too

            if (diagramType === "xy_constP" || diagramType === "xy_constT") {
                // For xy diagrams, y-axis is y1, x-axis is x1
                xAxisName = `Liquid Mole Fraction ${comp1Name} (x₁)`;
                yAxisName = `Vapor Mole Fraction ${comp1Name} (y₁)`;
                const xyData = relevantResults.map(p => [p.comp1_feed, p.comp1_equilibrium] as [number,number]).filter(p => !isNaN(p[1]));
                // xyData is already sorted by comp1_feed due to relevantResults sort
                const diagonalData = [[0,0], [1,1]];

                const equilibriumSeries: LineSeriesOption = {
                    name: 'Equilibrium Line (y₁ vs x₁)',
                    type: 'line',
                    smooth: true,
                    symbol: 'none',
                    lineStyle: { color: liquidColor, width: 2 },
                    data: xyData
                };

                const azeotropeTolerance = 1e-3; // Tolerance for x1 approx y1
                const azeotropesFound: any[] = [];

                // Check for points directly on/near the diagonal
                relevantResults.forEach(p => {
                    if (p.comp1_feed > 1e-6 && p.comp1_feed < (1.0 - 1e-6) && 
                        typeof p.comp1_equilibrium === 'number' && !isNaN(p.comp1_equilibrium)) {
                        if (Math.abs(p.comp1_feed - p.comp1_equilibrium) < azeotropeTolerance) {
                            const x_az_candidate = p.comp1_feed; 
                            // Avoid duplicate if already found by interpolation or very close points
                            if (!azeotropesFound.some(az => Math.abs(az.coord[0] - x_az_candidate) < azeotropeTolerance)) {
                                azeotropesFound.push({
                                    name: 'Azeotrope',
                                    coord: [x_az_candidate, x_az_candidate], // Plot on the y=x diagonal
                                    value: x_az_candidate.toFixed(3)
                                });
                            }
                        }
                    }
                });

                // Check for crossings between points
                for (let i = 0; i < relevantResults.length - 1; i++) {
                    const p1 = relevantResults[i];
                    const p2 = relevantResults[i+1];

                    if (typeof p1.comp1_equilibrium !== 'number' || isNaN(p1.comp1_equilibrium) ||
                        typeof p2.comp1_equilibrium !== 'number' || isNaN(p2.comp1_equilibrium)) {
                        continue; // Skip if equilibrium data is invalid for the segment
                    }

                    const x1_i = p1.comp1_feed;
                    const y1_i = p1.comp1_equilibrium;
                    const x1_j = p2.comp1_feed;
                    const y1_j = p2.comp1_equilibrium;

                    // Ensure points are distinct to avoid division by zero in interpolation
                    if (Math.abs(x1_i - x1_j) < 1e-9) continue;

                    const func_i = x1_i - y1_i;
                    const func_j = x1_j - y1_j;

                    // Check for crossing: func_i and func_j have opposite signs
                    if (func_i * func_j < 0) {
                        // Linear interpolation for x_azeotrope where x_azeotrope = y_azeotrope(x_azeotrope)
                        // Solve for x_az in: (x_az - x1_i) / (x1_j - x1_i) = (0 - func_i) / (func_j - func_i)
                        const x_azeotrope = x1_i - func_i * (x1_j - x1_i) / (func_j - func_i);

                        // Ensure interpolated point is within the segment bounds and not at pure component ends
                        // (handles cases where x1_i > x1_j or vice-versa if data wasn't perfectly sorted, though it should be)
                        const min_x_segment = Math.min(x1_i, x1_j);
                        const max_x_segment = Math.max(x1_i, x1_j);

                        if (x_azeotrope > min_x_segment - 1e-9 && x_azeotrope < max_x_segment + 1e-9 && // Check within segment (with small tolerance)
                            x_azeotrope > 1e-6 && x_azeotrope < (1.0 - 1e-6)) { // Exclude pure component ends
                            
                            // Avoid duplicate if already found or very close
                            if (!azeotropesFound.some(az => Math.abs(az.coord[0] - x_azeotrope) < azeotropeTolerance)) {
                                azeotropesFound.push({
                                    name: 'Azeotrope',
                                    coord: [x_azeotrope, x_azeotrope], // Plot on the y=x diagonal
                                    value: x_azeotrope.toFixed(3)
                                });
                            }
                        }
                    }
                }
                
                // Sort azeotropes by x-coordinate if multiple are found
                azeotropesFound.sort((a,b) => a.coord[0] - b.coord[0]);

                if (azeotropesFound.length > 0) {
                    equilibriumSeries.markPoint = {
                        symbol: 'circle', 
                        symbolSize: 8,    
                        itemStyle: {
                            color: 'red', 
                            borderColor: isDark ? '#333' : '#fff',
                            borderWidth: 1,
                        },
                        label: {
                            show: true,
                            formatter: '{b}\nx₁={c}', 
                            position: 'top',       
                            color: textColor,
                            fontSize: 10,
                            fontFamily: chartFontFamily, // Apply font to markPoint label
                            backgroundColor: isDark ? 'rgba(50,50,50,0.8)' : 'rgba(255,255,255,0.8)',
                            padding: [3, 5],
                            borderRadius: 3,
                        },
                        data: azeotropesFound
                    };
                }

                series = [
                    equilibriumSeries,
                    { name: 'Diagonal (y₁=x₁)', type: 'line', symbol: 'none', lineStyle: { color: lineColor, width: 1, type: 'dashed' }, data: diagonalData }
                ];
            } else { // Txy or Pxy
                // Split Vapor Line Data Logic for Txy/Pxy (min/max point splitting)
                let splitPointIndex = -1;
                if (vaporLineData.length > 1) {
                    if (diagramType === "Txy") { // Find min T for Txy
                        let minVal = Infinity;
                        vaporLineData.forEach((p, index) => { if (p[1] < minVal) { minVal = p[1]; splitPointIndex = index; } });
                    } else if (diagramType === "Pxy") { // Find max P for Pxy
                        let maxVal = -Infinity;
                        vaporLineData.forEach((p, index) => { if (p[1] > maxVal) { maxVal = p[1]; splitPointIndex = index; } });
                    }
                }


                let vaporPart1: [number, number][] = [];
                let vaporPart2: [number, number][] = [];

                if (splitPointIndex !== -1 && splitPointIndex < vaporLineData.length -1 && vaporLineData.length > 1 ) {
                    vaporPart1 = vaporLineData.slice(0, splitPointIndex + 1);
                    vaporPart2 = vaporLineData.slice(splitPointIndex);
                    // Already sorted by x-coordinate (y1)
                    if (vaporPart1.length > 0 && vaporPart2.length > 0 &&
                        Math.abs(vaporPart1[vaporPart1.length-1][0] - vaporPart2[0][0]) < 1e-9 && // Compare x
                        Math.abs(vaporPart1[vaporPart1.length-1][1] - vaporPart2[0][1]) < 1e-9 ) { // Compare y
                       vaporPart2[0][0] += 1e-9; // Nudge x if identical to avoid rendering issues
                    }
                } else {
                    vaporPart1 = [...vaporLineData]; // Plot as single sorted line if no clear split
                    vaporPart2 = [];
                }


                series = [
                    { name: diagramType === "Txy" ? 'Liquid Line (T-x₁)' : 'Liquid Line (P-x₁)', type: 'line', smooth: true, symbol: 'none', lineStyle: { color: liquidColor, width: 2 }, data: liquidLineData },
                    { name: diagramType === "Txy" ? 'Vapor Line (T-y₁)' : 'Vapor Line (P-y₁)', type: 'line', smooth: false, symbol: 'none', lineStyle: { color: vaporColor, width: 2 }, data: vaporPart1, legendHoverLink: false, },
                     ...(vaporPart2.length > 0 ? [{ name: diagramType === "Txy" ? 'Vapor Line (T-y₁)' : 'Vapor Line (P-y₁)', type: 'line' as const, smooth: false, symbol: 'none' as const, lineStyle: { color: vaporColor, width: 2 }, data: vaporPart2, legendHoverLink: false }] : []) 
                ];
            }


            setVleChartOptions({
                title: { text: titleText, left: 'center', textStyle: { color: textColor, fontSize: 16, fontFamily: chartFontFamily } },
                tooltip: { 
                    trigger: 'axis', 
                    axisPointer: { type: 'cross' }, 
                    textStyle: { fontFamily: chartFontFamily },
                    formatter: (params: any) => { // Custom formatter for VLE diagrams
                        if (!params || params.length === 0) return '';
                        let tooltipString = `${params[0].axisValueLabel}<br/>`; // x-axis value (mole fraction or T/P)
                        params.forEach((param: any) => {
                            let seriesName = param.seriesName;
                            let value = parseFloat(param.value[1]); // y-axis value
                            
                            // Customize series name for legend/tooltip if needed
                            if (seriesName.includes('Liquid Line')) seriesName = `Liquid (${yAxisName.split(' ')[0]})`;
                            if (seriesName.includes('Vapor Line')) seriesName = `Vapor (${yAxisName.split(' ')[0]})`;
                            if (seriesName.includes('Equilibrium Line')) seriesName = `y₁`;
                            if (seriesName.includes('Diagonal')) seriesName = `y₁ = x₁`;

                            let formattedValue = '';
                            if (yAxisName.includes("Temperature")) {
                                formattedValue = value.toFixed(1) + " K";
                            } else if (yAxisName.includes("Pressure")) {
                                formattedValue = (value / 1000).toPrecision(4) + " kPa";
                            } else { // Mole fraction for y-axis (xy diagram)
                                formattedValue = value.toFixed(3);
                            }
                            tooltipString += `${param.marker} ${seriesName}: ${formattedValue}<br/>`;
                        });
                        return tooltipString;
                    }
                },
                legend: { data: series.map((s: LineSeriesOption) => s.name as string).filter((name?: string) => name !== undefined && name !== null), bottom: 5, textStyle: { color: textColor, fontFamily: chartFontFamily } }, 
                grid: { left: '10%', right: '5%', bottom: '15%', top: '15%', containLabel: true }, 
                xAxis: {
                     type: 'value', name: xAxisName, nameLocation: 'middle', nameGap: 30, min: 0, max: 1,
                     axisLabel: { color: textColor, fontSize: 12, fontFamily: chartFontFamily }, 
                     nameTextStyle: { color: textColor, fontSize: 14, fontWeight: 'bold', fontFamily: chartFontFamily },
                     axisLine: { lineStyle: { color: lineColor } }, 
                     splitLine: { show: false }, // REMOVE GRID LINES
                },
                yAxis: {
                    type: 'value', name: yAxisName, nameLocation: 'middle', nameGap: 55, 
                    scale: true,
                    axisLabel: { 
                        color: textColor, 
                        fontSize: 12, 
                        fontFamily: chartFontFamily,
                        formatter: (value: number) => {
                            if (diagramType === "Pxy") { // Y-axis is Pressure (Pa), show as kPa
                                return (value / 1000).toPrecision(4);
                            } else if (diagramType === "Txy") { // Y-axis is Temperature (K)
                                return value.toFixed(1);
                            } else { // xy_constP or xy_constT, Y-axis is Vapor Mole Fraction (y1)
                                return value.toFixed(3); // Mole fractions typically 0-1
                            }
                        } 
                    },
                    nameTextStyle: { color: textColor, fontSize: 14, fontWeight: 'bold', fontFamily: chartFontFamily },
                    axisLine: { lineStyle: { color: lineColor } }, 
                    splitLine: { show: false }, // REMOVE GRID LINES
                 },
                series: series.filter((s: LineSeriesOption) => s && s.data && (s.data as any[]).length > 0) 
            });
        } else {
            setVleChartOptions({}); // Clear chart if no data or mode switch
        }
    }, [results, comp1Name, comp2Name, pressureKPa, temperatureK, diagramType, resolvedTheme, azeotropeFinderActive, azeotropeScanData, azeotropeScanType]);


    // --- Local Helper Functions (dPsat_dT_Pa, newtonRaphsonSolve, antoineBoilingPointSolver) ---
    /** Derivative of Saturation Pressure (Pa/K) wrt Temperature (K) - Local Helper */
    function dPsat_dT_Pa(params: AntoineParams, T_kelvin: number): number { // Uses imported AntoineParams
        if (!params || (T_kelvin + params.C === 0)) return NaN;
        // Uses imported libCalculatePsat_Pa for consistency
        const Psat_Pa_at_T = libCalculatePsat_Pa(params, T_kelvin); // Imported
        if (isNaN(Psat_Pa_at_T)) return NaN;
        let factor = 1.0; if (params.EquationNo === 1 || params.EquationNo === '1') factor = Math.log(10);
        return Psat_Pa_at_T * factor * params.B / Math.pow(T_kelvin + params.C, 2);
    }
    /** Newton-Raphson Root-finder - Local Helper for pure BP estimation */
    function newtonRaphsonSolve(f: (x: number) => number, df: (x: number) => number, initialGuess: number, tolerance: number = 1e-7, maxIter: number = 100): number | null {
        let x = initialGuess;
        for (let iter = 0; iter < maxIter; iter++) {
            const fx = f(x); if (isNaN(fx)) { console.error(`NR: f(x) NaN at x=${x}`); return null; }
            if (Math.abs(fx) < tolerance) return x;
            const dfx = df(x); if (isNaN(dfx) || Math.abs(dfx) < 1e-10) { console.error(`NR: df(x) NaN or too small (${dfx}) at x=${x}`); return null; }
            const x_next = x - fx / dfx; if (isNaN(x_next) || !isFinite(x_next)) { console.error(`NR: x_next NaN/Inf at x=${x}`); return null; }
            if (iter > 0 && Math.abs(x_next - x) < (tolerance * Math.abs(x_next) + tolerance)) return x_next;
            x = x_next;
        }
        // console.warn("NR failed for initial ${initialGuess}, last x=${x}, f(x)=${f(x)}");
        return null;
    }
    // Solver for pure component boiling point (used for initial temp guesses)
    const antoineBoilingPointSolver = (antoineParams: AntoineParams | null, P_target: number): number | null => { // Uses imported AntoineParams
       
        if(!antoineParams) return null;
        const func = (T: number) => libCalculatePsat_Pa(antoineParams, T) - P_target; // Imported
        const dfunc = (T: number) => dPsat_dT_Pa(antoineParams, T); // Local dPsat_dT_Pa which internally uses libCalculatePsat_Pa
        let T_initial_guess = (antoineParams.Tmin_K > 0 && antoineParams.Tmax_K < 10000 && antoineParams.Tmax_K > antoineParams.Tmin_K)
            ? (antoineParams.Tmin_K + antoineParams.Tmax_K) / 2
            : 373.15;
        T_initial_guess = Math.max(150, Math.min(T_initial_guess, 700)); // Bound initial guess
        return newtonRaphsonSolve(func, dfunc, T_initial_guess, 1e-5, 50);
    };


    // --- JSX Rendering ---
    return (
        <div className="container mx-auto p-4">
            <Card className="mb-4">
                <CardHeader><CardTitle>VLE Phase Diagram Calculator</CardTitle></CardHeader>
                <CardContent className="space-y-6">
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                        <div><label htmlFor="comp1" className="block text-sm font-medium mb-1">Compound 1</label><Input id="comp1" value={comp1Name} onChange={(e) => setComp1Name(e.target.value)} placeholder="e.g., Acetone" /></div>
                        <div><label htmlFor="comp2" className="block text-sm font-medium mb-1">Compound 2</label><Input id="comp2" value={comp2Name} onChange={(e) => setComp2Name(e.target.value)} placeholder="e.g., Methanol" /></div>
                    </div>

                    {/* Fluid Package Selector */}
                    <div>
                        <label htmlFor="fluidPackage" className="block text-sm font-medium mb-1">Fluid Package</label>
                        <Select value={fluidPackage} onValueChange={(value) => setFluidPackage(value as FluidPackageType)}>
                            <SelectTrigger><SelectValue placeholder="Select fluid package" /></SelectTrigger>
                            <SelectContent>
                                <SelectItem value="unifac">UNIFAC</SelectItem>
                                <SelectItem value="nrtl">NRTL</SelectItem>
                                <SelectItem value="pr">Peng-Robinson</SelectItem>
                                <SelectItem value="srk">SRK</SelectItem>
                                <SelectItem value="uniquac">UNIQUAC</SelectItem>
                                <SelectItem value="wilson">Wilson</SelectItem> {/* ADDED Wilson Option */}
                            </SelectContent>
                        </Select>
                    </div>

                    <Button onClick={() => setAzeotropeFinderActive(!azeotropeFinderActive)} variant="outline" className="w-full">
                        {azeotropeFinderActive ? "Switch to VLE Diagrams" : "Switch to Azeotrope Finder"}
                    </Button>

                    {azeotropeFinderActive ? (
                        <div className="space-y-4 p-4 border rounded-md">
                            <h3 className="text-lg font-semibold">Azeotrope Finder Settings</h3>
                            <div>
                                <label htmlFor="azeotropeScanType" className="block text-sm font-medium mb-1">Scan Type</label>
                                <Select value={azeotropeScanType} onValueChange={(value) => setAzeotropeScanType(value as AzeotropeScanType)}>
                                    <SelectTrigger><SelectValue /></SelectTrigger>
                                    <SelectContent>
                                        <SelectItem value="vs_P_find_T">Scan Pressure (find x_az, T_az)</SelectItem>
                                        <SelectItem value="vs_T_find_P">Scan Temperature (find x_az, P_az)</SelectItem>
                                    </SelectContent>
                                </Select>
                            </div>
                            {/* REMOVED Scan Start and Scan End inputs */}
                            <div>
                                <label htmlFor="scanSteps" className="block text-sm font-medium mb-1">Number of Scan Points</label>
                                <Input id="scanSteps" type="number" value={scanSteps} onChange={(e) => setScanSteps(parseInt(e.target.value))} min="5" max="200" />
                                <p className="text-xs text-muted-foreground mt-1">
                                    Defines the resolution of the automatic wide-range scan. Recommended: 30-100.
                                </p>
                            </div>
                        </div>
                    ) : (
                        <>
                            <div>
                                <label htmlFor="diagramType" className="block text-sm font-medium mb-1">Diagram Type</label>
                                <Select value={diagramType} onValueChange={(value) => setDiagramType(value as DiagramType)}>
                                    <SelectTrigger><SelectValue placeholder="Select diagram type" /></SelectTrigger>
                                    <SelectContent>
                                        <SelectItem value="Txy">T-xy Diagram (Constant P)</SelectItem>
                                        <SelectItem value="Pxy">P-xy Diagram (Constant T)</SelectItem>
                                        <SelectItem value="xy_constP">x-y Diagram (Constant P)</SelectItem>
                                        <SelectItem value="xy_constT">x-y Diagram (Constant T)</SelectItem>
                                    </SelectContent>
                                </Select>
                            </div>
                            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                                <div>
                                    <label htmlFor="temperature" className="block text-sm font-medium mb-1">System Temperature (K)</label>
                                    <Input 
                                        id="temperature" type="number" value={temperatureK} 
                                        onChange={(e) => setTemperatureK(parseFloat(e.target.value))} placeholder="e.g., 323.15"
                                        disabled={diagramType === "Txy" || diagramType === "xy_constP"}
                                       
                                        className={(diagramType === "Txy" || diagramType === "xy_constP") ? "opacity-50 cursor-not-allowed" : ""}
                                    />
                                </div>
                                <div>
                                    <label htmlFor="pressure" className="block text-sm font-medium mb-1">System Pressure (kPa)</label>
                                    <Input 
                                        id="pressure" type="number" value={pressureKPa}
                                        onChange={(e) => setPressureKPa(parseFloat(e.target.value))} placeholder="e.g., 101.325"
                                        disabled={diagramType === "Pxy" || diagramType === "xy_constT"}
                                        className={(diagramType === "Pxy" || diagramType === "xy_constT") ? "opacity-50 cursor-not-allowed" : ""}
                                    />
                                </div>
                            </div>
                        </>
                    )}

                    <Button onClick={handleCalculate} disabled={isLoading} className="w-full md:w-auto">
                        {isLoading ? 'Calculating...' : (azeotropeFinderActive ? 'Find Azeotropes' : `Calculate ${diagramType} Diagram`)}
                    </Button>
                </CardContent>
            </Card>
            {error && (<Alert variant="destructive" className="mb-4"><Terminal className="h-4 w-4" /><AlertTitle>Error</AlertTitle><AlertDescription>{error}</AlertDescription></Alert>)}
            {/* Loading Message */}
            {isLoading && <p className="text-center my-4">Loading and calculating, please wait...</p>}
            {/* Chart Area or "No Results" Message - Rendered when not loading */}
            {!isLoading && (
                <>
                    {/* Chart Card: Render if in azeotropeFinder mode OR in VLE diagram mode with results */}
                    {(azeotropeFinderActive || (!azeotropeFinderActive && results.length > 0)) && (
                        <Card className="mb-4">
                            <CardContent className="pt-6">
                                {Object.keys(vleChartOptions).length > 0 ? (
                                    <div className="relative h-[500px] md:h-[600px] rounded-md border bg-card">
                                        <ReactECharts
                                            echarts={echarts}
                                            option={vleChartOptions}
                                            style={{ height: '100%', width: '100%' }}
                                            notMerge={true}
                                            lazyUpdate={true}
                                            theme={resolvedTheme === 'dark' ? 'dark' : undefined}
                                        />
                                    </div>
                                ) : (
                                    // Fallback if vleChartOptions is empty for some reason.
                                    // useEffect should ensure vleChartOptions is set with a "no data" message if needed.
                                    (<p className="text-muted-foreground text-center h-[500px] md:h-[600px] flex items-center justify-center">
                                        {azeotropeFinderActive && azeotropeScanData.length === 0 && !error ? 
                                            "No azeotropes found within the scanned range." : 
                                            "Generating chart..."
                                        }
                                    </p>)
                                )}
                                <p className="text-xs text-muted-foreground mt-2 text-center">
                                    Note: Accuracy depends on database parameters and model limitations.
                                </p>
                            </CardContent>
                        </Card>
                    )}

                    {/* "Calculate VLE diagram" Message: Render if in VLE mode, no results, not loading, and no error */}
                    {!azeotropeFinderActive && results.length === 0 && !error && (
                        <Card className="mb-4">
                            <CardContent className="pt-6">
                                <p className="text-muted-foreground text-center h-[500px] md:h-[600px] flex items-center justify-center">
                                    Calculate a VLE diagram to see results.
                                </p>
                            </CardContent>
                        </Card>
                    )}
                </>
            )}
            <Card>
                <CardHeader><CardTitle>Calculation Log</CardTitle></CardHeader>
                <CardContent><div className="h-64 overflow-y-auto bg-muted p-2 rounded text-xs font-mono">{logMessages.length === 0 ? <p>No messages yet...</p> : logMessages.map((msg, index) => <div key={index}>{msg}</div>)}</div></CardContent>
            </Card>
        </div>
    );
}

// REMOVE local interfaces and calculation functions that are now imported from vle-calculations.ts library
// Keep only essential local helpers if any (like the pure component BP solver if you prefer it local)