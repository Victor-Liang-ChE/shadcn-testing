'use client';

import React, { useState, useEffect } from 'react';
import { createClient, SupabaseClient } from '@supabase/supabase-js';
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"; // Added for diagram type
import { Terminal } from "lucide-react";
import { useTheme } from "next-themes";

import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsCoreOption } from 'echarts/core'; // CHANGED: EChartsOption to EChartsCoreOption, REMOVED: SeriesOption
import { LineChart, LineSeriesOption } from 'echarts/charts'; // LineSeriesOption is already imported
import {
    TitleComponent,
    TooltipComponent,
    GridComponent,
    LegendComponent,
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

// Import calculation logic and types from the library
import {
    calculatePsat_Pa as libCalculatePsat_Pa, // Renamed to avoid conflict
    calculateUnifacGamma as libCalculateUnifacGamma, // Renamed
    // NEW: Import all four calculation functions
    calculateBubbleTemperature,
    calculateBubblePressure,
    calculateDewTemperature,
    calculateDewPressure,
    // Keep necessary type imports:
    type AntoineParams,
    type UnifacGroupComposition,
    type CompoundData,
    type UnifacParameters,
    type BubbleDewResult // Use the new result type
} from '@/lib/vle-calculations'; // Adjust path as needed

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
type AzeotropeScanType = 'vs_P_find_T' | 'vs_T_find_P'; // Vary P, find T_az OR Vary T, find P_az

export default function VleCalculatorPage() {
    const [comp1Name, setComp1Name] = useState<string>("Acetone");
    const [comp2Name, setComp2Name] = useState<string>("Methanol");
    const [pressureKPa, setPressureKPa] = useState<number>(101.325); // Used for Txy and xy_constP
    const [temperatureK, setTemperatureK] = useState<number>(323.15); // Used for Pxy and xy_constT
    const [isLoading, setIsLoading] = useState<boolean>(false);
    const [error, setError] = useState<string | null>(null);
    const [results, setResults] = useState<BubbleDewResult[]>([]); // For VLE diagrams
    const [logMessages, setLogMessages] = useState<string[]>([]);
    const { resolvedTheme } = useTheme();
    const [vleChartOptions, setVleChartOptions] = useState<EChartsCoreOption>({});
    const [diagramType, setDiagramType] = useState<DiagramType>("Txy");

    // --- MODIFIED: Azeotrope Finder States ---
    const [azeotropeFinderActive, setAzeotropeFinderActive] = useState<boolean>(false);
    const [azeotropeScanType, setAzeotropeScanType] = useState<AzeotropeScanType>('vs_P_find_T');
    // REMOVED: scanStartValue, scanEndValue states
    const [scanSteps, setScanSteps] = useState<number>(30); // Increased default
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


    // --- Data Fetching Functions (fetchCompoundData, fetchUnifacInteractionParams) ---
    // Ensure they use the types imported from vle-calculations
    async function fetchCompoundData(compoundName: string): Promise<CompoundData | null> {
        if (!supabase) { throw new Error("Supabase client not initialized."); }
        console.log(`Fetching data for ${compoundName}...`);
        setError(null); // Clear previous errors specific to this fetch

        try {
            // 1. Find compound ID
            const { data: compoundData, error: compoundError } = await supabase
                .from('compounds')
                .select('id, name')
                .ilike('name', compoundName)
                .limit(1); // Take the first match

            if (compoundError) throw new Error(`Supabase compound query error: ${compoundError.message}`);
            if (!compoundData || compoundData.length === 0) {
                throw new Error(`Compound '${compoundName}' not found.`);
            }
            const compoundId = compoundData[0].id;
            const foundName = compoundData[0].name;
            console.log(`Found ${foundName} (ID: ${compoundId})`);

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
            let antoine: AntoineParams | null = null; // Use imported type
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


            // 4. Extract UNIFAC Group Composition
            let unifacGroups: UnifacGroupComposition | null = null; // Use imported type
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
                    unifacGroups = null; // Set back to null if empty/invalid
                }
            } else {
                console.warn(`No UNIFAC group data found in properties.elements_composition.UNIFAC for ${foundName}.`);
            }

             if (!unifacGroups) {
                 throw new Error(`Failed to extract valid UNIFAC groups for ${foundName}.`);
             }


            return { name: foundName, antoine, unifacGroups };
        } catch (err: any) {
            console.error(`Error fetching data for ${compoundName}:`, err.message);
            setError(`Failed to fetch data for ${compoundName}: ${err.message}`);
            return null;
        }
    }
    async function fetchUnifacInteractionParams(subgroupIds: number[]): Promise<UnifacParameters | null> {
        if (!supabase) { throw new Error("Supabase client not initialized."); }
        if (subgroupIds.length === 0) return { Rk: {}, Qk: {}, mainGroupMap: {}, a_mk: new Map() }; // Return empty valid structure
        console.log("Fetching UNIFAC parameters for subgroups:", subgroupIds);
        setError(null);

        try {
            // 1. Fetch Rk, Qk, and Main Group Mapping for required subgroups
            // ---- MODIFICATION: REMOVE explicit quotes from table name string ----
            const { data: groupData, error: groupError } = await supabase
                .from('UNIFAC - Rk and Qk') // Pass exact name without extra quotes
                // Select statement still needs quotes for columns with special chars/caps
                .select('"Subgroup #", "Main Group #", "Rk", "Qk"')
                .in('"Subgroup #"', subgroupIds); // Filter column still needs quotes

            if (groupError) throw new Error(`Supabase UNIFAC subgroup query error: ${groupError.message}`);
            // ---- MODIFICATION: Improved check for missing data ----
            // Ensure you handle the case where groupData might be empty or shorter than subgroupIds length
             if (!groupData || groupData.length < subgroupIds.length) {
                const foundIds = new Set(groupData?.map(g => g["Subgroup #"]) ?? []);
                const missingIds = subgroupIds.filter(id => !foundIds.has(id));
                 throw new Error(`Missing UNIFAC parameters for subgroup ID(s): ${missingIds.join(', ')}`);
             }


            const Rk: { [id: number]: number } = {};
            const Qk: { [id: number]: number } = {};
            const mainGroupMap: { [id: number]: number } = {};
            const mainGroupIds = new Set<number>();

            groupData.forEach(g => {
                 // ---- MODIFICATION: Access data using correct quoted column names ----
                const subgroupId = g["Subgroup #"];
                const mainGroupId = g["Main Group #"];
                const rkVal = g["Rk"]; // Using quoted names for consistency/safety
                const qkVal = g["Qk"];

                if (subgroupId == null || mainGroupId == null || rkVal == null || qkVal == null) {
                     console.warn(`Incomplete data for subgroup ${subgroupId}. Skipping.`);
                     return; // Skip incomplete entries
                }
                Rk[subgroupId] = rkVal;
                Qk[subgroupId] = qkVal;
                mainGroupMap[subgroupId] = mainGroupId;
                mainGroupIds.add(mainGroupId);
            });

             // Check if all requested subgroups were found and mapped (already partially done above)
             for (const reqId of subgroupIds) {
                 if (mainGroupMap[reqId] === undefined) {
                     // This condition should ideally be caught by the check after the query now
                     throw new Error(`Failed to retrieve parameters or main group mapping for subgroup ID: ${reqId}`);
                 }
             }

            console.log("Fetched Rk, Qk, Main Group Map. Required Main Groups:", Array.from(mainGroupIds));

            // 2. Fetch Interaction Parameters (a_mk) for all pairs of involved main groups
            const mainGroupArray = Array.from(mainGroupIds);
            if (mainGroupArray.length === 0) {
                 console.warn("No main groups identified, cannot fetch interaction parameters.");
                 return { Rk, Qk, mainGroupMap, a_mk: new Map() };
            }

            // ---- MODIFICATION: REMOVE explicit quotes from table name string ----
            const { data: interactionData, error: interactionError } = await supabase
                .from('UNIFAC - a(ij)') // Pass exact name without extra quotes
                 // Select statement still needs quotes for columns with special chars/caps
                .select('i, j, "a(ij)"')
                .in('i', mainGroupArray)
                .in('j', mainGroupArray);


            if (interactionError) throw new Error(`Supabase UNIFAC interaction query error: ${interactionError.message}`);
            if (!interactionData) throw new Error("Failed to query UNIFAC interaction parameters."); // Should return empty array, not null

            const a_mk = new Map<string, number>();
            interactionData.forEach(interaction => {
                 // ---- MODIFICATION: Access data using correct column names ----
                const main_group_m = interaction.i;
                const main_group_k = interaction.j;
                const a_mk_value = interaction["a(ij)"]; // Access column with quotes

                if (main_group_m != null && main_group_k != null && a_mk_value != null) {
                    // Still use the conceptual main group IDs for the map key
                    a_mk.set(`${main_group_m}-${main_group_k}`, a_mk_value);
                } else {
                    console.warn(`Incomplete interaction data found: i=${main_group_m}, j=${main_group_k}. Skipping.`);
                }
            });

            // Verify all required pairs are present (logic remains the same, uses mainGroupArray)
            // ---- MODIFICATION: Warn and default missing interactions to 0 ----
            for (const mg1 of mainGroupArray) {
                for (const mg2 of mainGroupArray) {
                    const key = `${mg1}-${mg2}`;
                    if (!a_mk.has(key)) {
                         if (mg1 === mg2) {
                             console.warn(`Missing diagonal interaction parameter for main group ${mg1}. Assuming 0.`);
                             a_mk.set(key, 0.0);
                         } else {
                             // Consider if missing off-diagonal parameters should default to 0 or throw error
                              console.warn(`Missing UNIFAC interaction parameter for main groups: ${key}. Assuming 0 (check database).`);
                              a_mk.set(key, 0.0); // Defaulting to 0, adjust if needed
                            // throw new Error(`Missing UNIFAC interaction parameter for main groups: ${key}`);
                         }
                    }
                }
            }

            console.log("Fetched UNIFAC interaction parameters (a_mk).");
            return { Rk, Qk, mainGroupMap, a_mk };
        } catch (err: any) {
            console.error("Error fetching UNIFAC parameters:", err.message);
            setError(`Failed to fetch UNIFAC parameters: ${err.message}`);
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
        // This function now contains the original logic of handleCalculate
        setIsLoading(true); setError(null); setResults([]); setLogMessages([]);
        console.log(`\n--- Starting VLE Calculation (${diagramType}) ---`);
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
            if (!data1 || !data2 || !data1.antoine || !data2.antoine || !data1.unifacGroups || !data2.unifacGroups) {
                setError(error || "Missing essential data for one or both compounds."); setIsLoading(false); return;
            }
            
            data1.r_i = undefined; data1.q_i = undefined; 
            data2.r_i = undefined; data2.q_i = undefined; 

            const components: CompoundData[] = [data1, data2];
            const allSubgroupIds = new Set<number>();
            components.forEach(comp => { if (comp.unifacGroups) Object.keys(comp.unifacGroups).forEach(id => allSubgroupIds.add(parseInt(id))); });
            const unifacParams = await fetchUnifacInteractionParams(Array.from(allSubgroupIds));
            if (!unifacParams) { setIsLoading(false); return; }

            const targetPressurePa = pressureKPa * 1000; 
            const targetTemperatureK = temperatureK;   

            const x_feed_values = Array.from({ length: 21 }, (_, i) => parseFloat((i * 0.05).toFixed(3))); 
            const calculatedResults: BubbleDewResult[] = [];

            // Estimate pure component saturation properties
            let T_bp1_est_at_targetP: number | undefined;
            let T_bp2_est_at_targetP: number | undefined;
            let P_sat1_est_at_targetT: number | undefined;
            let P_sat2_est_at_targetT: number | undefined;

            if (data1.antoine) T_bp1_est_at_targetP = antoineBoilingPointSolver(data1.antoine, targetPressurePa) ?? undefined;
            if (data2.antoine) T_bp2_est_at_targetP = antoineBoilingPointSolver(data2.antoine, targetPressurePa) ?? undefined;
            if (data1.antoine) P_sat1_est_at_targetT = libCalculatePsat_Pa(data1.antoine, targetTemperatureK);
            if (data2.antoine) P_sat2_est_at_targetT = libCalculatePsat_Pa(data2.antoine, targetTemperatureK);

            for (const x1_feed of x_feed_values) {
                console.log(`\nCalculating for feed mole fraction ${comp1Name} (x₁) = ${x1_feed.toFixed(3)}...`);

                // --- Calculate Bubble Temperature at targetPressurePa ---
                const initialTempGuessForBubbleT = (T_bp1_est_at_targetP && T_bp2_est_at_targetP 
                    ? (x1_feed * T_bp1_est_at_targetP + (1 - x1_feed) * T_bp2_est_at_targetP) 
                    : targetTemperatureK) || 350; // Fallback to targetTemperatureK or 350K
                
                const bubbleT_resultPoint = calculateBubbleTemperature(components, x1_feed, targetPressurePa, unifacParams, initialTempGuessForBubbleT);
                if (bubbleT_resultPoint && !bubbleT_resultPoint.error) {
                    calculatedResults.push(bubbleT_resultPoint);
                } else {
                    console.error(`   Bubble Temperature calculation failed for x₁ = ${x1_feed.toFixed(3)}. Error: ${bubbleT_resultPoint?.error || 'Unknown error'}`);
                    // Push error result to still have a placeholder if needed, or handle differently
                    if (bubbleT_resultPoint) calculatedResults.push(bubbleT_resultPoint); else calculatedResults.push({comp1_feed: x1_feed, comp1_equilibrium: NaN, error: "BubbleT calc object null", calculationType: 'bubbleT', P_Pa: targetPressurePa});
                }

                // --- Calculate Bubble Pressure at targetTemperatureK ---
                // initialPressureGuess for calculateBubblePressure is not strictly used by its current direct calculation logic, but good to provide
                const initialPressureGuessForBubbleP = (P_sat1_est_at_targetT && P_sat2_est_at_targetT 
                    ? (x1_feed * P_sat1_est_at_targetT + (1 - x1_feed) * P_sat2_est_at_targetT) 
                    : targetPressurePa) || 101325; // Fallback to targetPressurePa or 101325 Pa

                const bubbleP_resultPoint = calculateBubblePressure(components, x1_feed, targetTemperatureK, unifacParams, initialPressureGuessForBubbleP);
                if (bubbleP_resultPoint && !bubbleP_resultPoint.error) {
                    calculatedResults.push(bubbleP_resultPoint);
                } else {
                    console.error(`   Bubble Pressure calculation failed for x₁ = ${x1_feed.toFixed(3)}. Error: ${bubbleP_resultPoint?.error || 'Unknown error'}`);
                    if (bubbleP_resultPoint) calculatedResults.push(bubbleP_resultPoint); else calculatedResults.push({comp1_feed: x1_feed, comp1_equilibrium: NaN, error: "BubbleP calc object null", calculationType: 'bubbleP', T_K: targetTemperatureK});
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
            console.log(`\n--- Calculation Complete (Both T and P variations) ---`);
        } catch (err: any) {
            console.error("Calculation failed:", err);
            setError(`Calculation failed: ${err.message}`);
        } finally {
            setIsLoading(false);
        }
    };

    const handleAzeotropeScan = async () => {
        setIsLoading(true); setError(null); setAzeotropeScanData([]); setLogMessages([]);
        console.log(`\n--- Starting Azeotrope Scan (${azeotropeScanType}) ---`);
        console.log(`Compounds: ${comp1Name}/${comp2Name}`);
        // console.log(`Scan Range: ${scanStartValue} to ${scanEndValue}, Steps: ${scanSteps}`); // Commented out as range is now internal

        if (!comp1Name || !comp2Name) {
            setError("Please enter valid compound names."); setIsLoading(false); return;
        }
        if (!supabase) { setError("Supabase client is not available."); setIsLoading(false); return; }

        try {
            const [data1, data2] = await Promise.all([fetchCompoundData(comp1Name), fetchCompoundData(comp2Name)]);
            if (!data1 || !data2 || !data1.antoine || !data2.antoine || !data1.unifacGroups || !data2.unifacGroups) {
                setError(error || "Missing essential data for one or both compounds."); setIsLoading(false); return;
            }
            data1.r_i = undefined; data2.r_i = undefined; // Reset for UNIFAC
            const components: CompoundData[] = [data1, data2];
            const allSubgroupIds = new Set<number>();
            components.forEach(comp => { if (comp.unifacGroups) Object.keys(comp.unifacGroups).forEach(id => allSubgroupIds.add(parseInt(id))); });
            const unifacParams = await fetchUnifacInteractionParams(Array.from(allSubgroupIds));
            if (!unifacParams) { setIsLoading(false); return; }

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

                if (azeotropeScanType === 'vs_P_find_T') { // Vary P (kPa), find x_az and T_az (K)
                    const P_system_Pa = currentScanVal * 1000;
                    console.log(`Scanning for azeotrope at P = ${currentScanVal.toFixed(2)} kPa...`);

                    const objective_h = (x1_guess: number): number | null => {
                        if (x1_guess <= 1e-4 || x1_guess >= 1.0 - 1e-4) return 1.0; // Push away from edges

                        // Estimate boiling points at current P_system_Pa for a better initial T guess
                        const T_bp1_est_at_P = antoineBoilingPointSolver(components[0].antoine, P_system_Pa) || 350; // Default if solver fails
                        const T_bp2_est_at_P = antoineBoilingPointSolver(components[1].antoine, P_system_Pa) || 350;
                        const initialTempGuess = x1_guess * T_bp1_est_at_P + (1 - x1_guess) * T_bp2_est_at_P;
                        
                        const bubbleT_result = calculateBubbleTemperature(components, x1_guess, P_system_Pa, unifacParams, initialTempGuess, 30, 1e-4);
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
                        const final_bubble_res = calculateBubbleTemperature(components, x_az_found, P_system_Pa, unifacParams, initialTempGuessFinal, 30, 1e-4);
                        if (final_bubble_res && final_bubble_res.error === undefined && final_bubble_res.T_K) {
                            dependent_val_found = final_bubble_res.T_K;
                            console.log(`  Found: P=${currentScanVal.toFixed(2)} kPa -> x_az=${x_az_found.toFixed(4)}, T_az=${dependent_val_found.toFixed(2)} K`);
                        } else { x_az_found = null; }
                    }
                } else { // vs_T_find_P: Vary T (K), find x_az and P_az (kPa)
                    const T_system_K = currentScanVal;
                    console.log(`Scanning for azeotrope at T = ${T_system_K.toFixed(2)} K...`);
                    const objective_k = (x1_guess: number): number | null => {
                        if (x1_guess <= 1e-4 || x1_guess >= 1.0 - 1e-4) return 1.0;

                        const P_sat1_est = libCalculatePsat_Pa(components[0].antoine!, T_system_K);
                        const P_sat2_est = libCalculatePsat_Pa(components[1].antoine!, T_system_K);
                        if (isNaN(P_sat1_est) || isNaN(P_sat2_est)) return null;
                        const initialPressureGuess = x1_guess * P_sat1_est + (1 - x1_guess) * P_sat2_est;

                        const bubbleP_result = calculateBubblePressure(components, x1_guess, T_system_K, unifacParams, initialPressureGuess || 101325, 10, 1e-5);
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
                            const final_bubble_res = calculateBubblePressure(components, x_az_found, T_system_K, unifacParams, initialPressureGuessFinal || 101325, 10, 1e-5);
                            if (final_bubble_res && final_bubble_res.error === undefined && final_bubble_res.P_Pa) {
                                dependent_val_found = final_bubble_res.P_Pa / 1000; // Convert to kPa
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
            console.log(`\n--- Azeotrope Scan Complete ---`);
        } catch (err: any) {
            console.error("Azeotrope Scan failed:", err);
            setError(`Azeotrope Scan failed: ${err.message}`);
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

        if (azeotropeFinderActive && azeotropeScanData.length > 0) {
            const scanParamName = azeotropeScanType === 'vs_P_find_T' ? "Pressure (kPa)" : "Temperature (K)";
            const dependentParamName = azeotropeScanType === 'vs_P_find_T' ? "Azeotropic Temperature (K)" : "Azeotropic Pressure (kPa)";
            const titleText = `Azeotrope Locus for ${comp1Name}/${comp2Name}`;

            const x_az_data = azeotropeScanData.map(d => [d.scanVal, d.x_az]);
            const dependent_val_data = azeotropeScanData.map(d => [d.scanVal, d.dependentVal]);

            setVleChartOptions({
                title: { text: titleText, left: 'center', textStyle: { color: textColor, fontSize: 16 } },
                tooltip: { trigger: 'axis', axisPointer: { type: 'cross' } },
                legend: { data: ['Azeotropic Composition (x₁)', dependentParamName], bottom: 5, textStyle: { color: textColor } },
                grid: { left: '10%', right: '12%', bottom: '15%', top: '15%', containLabel: true }, // Adjusted for dual y-axis
                xAxis: {
                    type: 'value', name: scanParamName, nameLocation: 'middle', nameGap: 30,
                    axisLabel: { color: textColor, fontSize: 12 }, nameTextStyle: { color: textColor, fontSize: 14, fontWeight: 'bold' },
                    axisLine: { lineStyle: { color: lineColor } }, splitLine: { lineStyle: { color: isDark ? '#374151' : '#E5E7EB' } },
                    scale: true,
                },
                yAxis: [
                    {
                        type: 'value', name: 'Azeotropic Composition (x₁)', nameLocation: 'middle', nameGap: 55, min: 0, max: 1,
                        axisLabel: { color: textColor, fontSize: 12, formatter: (v:number) => v.toFixed(3) },
                        nameTextStyle: { color: textColor, fontSize: 14, fontWeight: 'bold' },
                        axisLine: { lineStyle: { color: liquidColor } }, splitLine: { show: false },
                    },
                    {
                        type: 'value', name: dependentParamName, nameLocation: 'middle', nameGap: 70, position: 'right',
                        axisLabel: { color: textColor, fontSize: 12, formatter: (v:number) => v.toFixed(azeotropeScanType === 'vs_T_find_P' ? 1 : 2) }, // kPa or K
                        nameTextStyle: { color: textColor, fontSize: 14, fontWeight: 'bold' },
                        axisLine: { lineStyle: { color: dependentParamColor } }, splitLine: { show: false },
                        scale: true,
                    }
                ],
                series: [
                    { name: 'Azeotropic Composition (x₁)', type: 'line', yAxisIndex: 0, symbol: 'circle', symbolSize: 6, lineStyle: { color: liquidColor, width: 2 }, data: x_az_data },
                    { name: dependentParamName, type: 'line', yAxisIndex: 1, symbol: 'triangle', symbolSize: 6, lineStyle: { color: dependentParamColor, width: 2 }, data: dependent_val_data }
                ]
            });

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
                title: { text: titleText, left: 'center', textStyle: { color: textColor, fontSize: 16 } },
                tooltip: { trigger: 'axis', axisPointer: { type: 'cross' }, /* TODO: Add formatter if needed */ },
                legend: { data: series.map((s: LineSeriesOption) => s.name as string).filter((name?: string) => name !== undefined && name !== null), bottom: 5, textStyle: { color: textColor } }, 
                grid: { left: '10%', right: '5%', bottom: '15%', top: '15%', containLabel: true }, 
                xAxis: {
                     type: 'value', name: xAxisName, nameLocation: 'middle', nameGap: 30, min: 0, max: 1,
                     axisLabel: { color: textColor, fontSize: 12 }, nameTextStyle: { color: textColor, fontSize: 14, fontWeight: 'bold' },
                     axisLine: { lineStyle: { color: lineColor } }, splitLine: { lineStyle: { color: isDark ? '#374151' : '#E5E7EB' } },
                },
                yAxis: {
                    type: 'value', name: yAxisName, nameLocation: 'middle', nameGap: 55, 
                    scale: true,
                    axisLabel: { 
                        color: textColor, 
                        fontSize: 12, 
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
                    nameTextStyle: { color: textColor, fontSize: 14, fontWeight: 'bold' },
                    axisLine: { lineStyle: { color: lineColor } }, splitLine: { lineStyle: { color: isDark ? '#374151' : '#E5E7EB' } },
                 },
                series: series.filter((s: LineSeriesOption) => s && s.data && (s.data as any[]).length > 0) 
            });
        } else {
            setVleChartOptions({}); // Clear chart if no data or mode switch
        }
    }, [results, comp1Name, comp2Name, pressureKPa, temperatureK, diagramType, resolvedTheme, azeotropeFinderActive, azeotropeScanData, azeotropeScanType]);


    // --- Local Helper Functions (dPsat_dT_Pa, newtonRaphsonSolve, antoineBoilingPointSolver) ---
    // ... (Keep your existing local helper functions if they are still used for pure component BP estimation)
    /** Derivative of Saturation Pressure (Pa/K) wrt Temperature (K) - Local Helper */
    function dPsat_dT_Pa(params: AntoineParams, T_kelvin: number): number { // AntoineParams is imported
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
        console.warn(`NR failed for initial ${initialGuess}, last x=${x}, f(x)=${f(x)}`);
        return null;
    }
    // Solver for pure component boiling point (used for initial temp guesses)
    const antoineBoilingPointSolver = (antoineParams: AntoineParams | null, P_target: number): number | null => {
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
            {isLoading && <p className="text-center my-4">Loading and calculating, please wait...</p>}

            {/* --- Chart Card --- */}
            {(results.length > 0 || azeotropeScanData.length > 0) && !isLoading && (
                 <Card className="mb-4">
                     <CardContent className="pt-6">
                         {Object.keys(vleChartOptions).length > 0 ? (
                              <div className="relative h-[500px] md:h-[600px] rounded-md border bg-card">
                                 {/* --- REPLACE PLOTLY WITH ECHARTS --- */}
                                 <ReactECharts
                                     echarts={echarts}
                                     option={vleChartOptions}
                                     style={{ height: '100%', width: '100%' }}
                                     notMerge={true} // Recommended for dynamic data
                                     lazyUpdate={true} // Recommended for performance
                                     theme={resolvedTheme === 'dark' ? 'dark' : undefined} // Optional: ECharts theme
                                 />
                              </div>
                          ) : (
                             <p className="text-muted-foreground text-center h-[500px] md:h-[600px] flex items-center justify-center">Generating chart...</p>
                          )}
                          <p className="text-xs text-muted-foreground mt-2 text-center">
                              Note: Accuracy depends on database parameters and model limitations.
                          </p>
                      </CardContent>
                  </Card>
             )}
             {/* --- END Chart Card --- */}

            <Card>
                <CardHeader><CardTitle>Calculation Log</CardTitle></CardHeader>
                <CardContent><div className="h-64 overflow-y-auto bg-muted p-2 rounded text-xs font-mono">{logMessages.length === 0 ? <p>No messages yet...</p> : logMessages.map((msg, index) => <div key={index}>{msg}</div>)}</div></CardContent>
            </Card>
        </div>
    );
}

// REMOVE local interfaces and calculation functions that are now imported from vle-calculations.ts library
// Keep only essential local helpers if any (like the pure component BP solver if you prefer it local)