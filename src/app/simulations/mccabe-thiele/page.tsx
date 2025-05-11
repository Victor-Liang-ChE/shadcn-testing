'use client';

import React, { useState, useEffect, useCallback, useRef } from 'react';
import { createClient, SupabaseClient } from '@supabase/supabase-js';

// Import ECharts components
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
// Import specific types from the main echarts package
import type { EChartsOption, SeriesOption } from 'echarts';
import { LineChart, ScatterChart } from 'echarts/charts';
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  MarkLineComponent,
  MarkPointComponent,
  ToolboxComponent
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

// Register necessary ECharts components
echarts.use([
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  MarkLineComponent,
  MarkPointComponent,
  ToolboxComponent,
  LineChart,
  ScatterChart,
  CanvasRenderer
]);

import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Slider } from "@/components/ui/slider";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Tooltip as ShadTooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { ArrowLeftRight } from 'lucide-react'; // Import swap icon

// Import VLE calculation logic and types
import {
    calculatePsat_Pa as libCalculatePsat_Pa, // Shared Psat
    // UNIFAC
    calculateUnifacGamma,
    calculateBubbleTemperature as calculateBubbleTemperatureUnifac,
    calculateBubblePressure as calculateBubblePressureUnifac,
    fetchUnifacInteractionParams,
    type UnifacParameters,
} from '@/lib/vle-calculations-unifac';

import {
    // NRTL
    fetchNrtlParameters,
    calculateNrtlGamma,
    calculateBubbleTemperatureNrtl,
    calculateBubblePressureNrtl,
    type NrtlInteractionParams
} from '@/lib/vle-calculations-nrtl';

import {
    // PR
    fetchPrInteractionParams,
    calculateBubbleTemperaturePr,
    calculateBubblePressurePr,
    type PrInteractionParams
} from '@/lib/vle-calculations-pr';

import {
    // SRK
    fetchSrkInteractionParams,
    calculateBubbleTemperatureSrk,
    calculateBubblePressureSrk,
    type SrkInteractionParams
} from '@/lib/vle-calculations-srk';

import {
    // UNIQUAC
    fetchUniquacInteractionParams,
    calculateBubbleTemperatureUniquac,
    calculateBubblePressureUniquac,
    type UniquacInteractionParams as LibUniquacInteractionParams // Alias to avoid conflict if local type exists
} from '@/lib/vle-calculations-uniquac';

import {
    // Wilson
    fetchWilsonInteractionParams,
    calculateBubbleTemperatureWilson,
    calculateBubblePressureWilson,
    type WilsonInteractionParams as LibWilsonInteractionParams // Alias
} from '@/lib/vle-calculations-wilson';

// Import Shared VLE Types
import type {
    AntoineParams,
    PrPureComponentParams,
    SrkPureComponentParams,
    UniquacPureComponentParams,
    WilsonPureComponentParams,
    CompoundData, // This will be the primary type for component properties
    BubbleDewResult,
    UnifacGroupComposition
} from '@/lib/vle-types';


type EChartsPoint = [number, number] | [number | null, number | null];
type FluidPackageTypeMcCabe = 'unifac' | 'nrtl' | 'pr' | 'srk' | 'uniquac' | 'wilson';

// --- Supabase Client Setup ---
const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL || 'YOUR_SUPABASE_URL';
const supabaseAnonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY || 'YOUR_SUPABASE_ANON_KEY';

let supabase: SupabaseClient;
try {
    supabase = createClient(supabaseUrl, supabaseAnonKey);
} catch (error) {
    console.error("Error initializing Supabase client for McCabe-Thiele:", error);
}

export default function McCabeThielePage() {
  // Input States - Updated defaults for initial load
  const [comp1Name, setComp1Name] = useState('methanol');
  const [comp2Name, setComp2Name] = useState('water');
  const [temperatureC, setTemperatureC] = useState<number | null>(60);
  const [pressureBar, setPressureBar] = useState<number | null>(1); // Default, but will be overridden by useTemperature=true
  const [useTemperature, setUseTemperature] = useState(true);
  const [fluidPackage, setFluidPackage] = useState<FluidPackageTypeMcCabe>('uniquac');

  // Parameter States
  const [xd, setXd] = useState(0.9);
  const [xb, setXb] = useState(0.1);
  const [xf, setXf] = useState(0.5);
  const [q, setQ] = useState(1.0);
  const [r, setR] = useState(1.5);

  const buffer = 0.01;

  // Data & Control States
  const [equilibriumData, setEquilibriumData] = useState<{ x: number[], y: number[] } | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Result States
  const [stages, setStages] = useState<number | null>(null);
  const [feedStage, setFeedStage] = useState<number | null>(null);

  // State for ECharts options - Use the imported EChartsOption type
  const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
  const echartsRef = useRef<ReactECharts | null>(null); // Ref for ECharts instance
  // State for displayed parameters in the title
  const [displayedComp1, setDisplayedComp1] = useState('');
  const [displayedComp2, setDisplayedComp2] = useState('');
  const [displayedTemp, setDisplayedTemp] = useState<number | null>(null);
  const [displayedPressure, setDisplayedPressure] = useState<number | null>(null);
  const [displayedUseTemp, setDisplayedUseTemp] = useState(true);
  const [displayedFluidPackage, setDisplayedFluidPackage] = useState<FluidPackageTypeMcCabe>('unifac');


  // --- useEffect for initial graph load ---
  useEffect(() => {
    // Automatically calculate and display the graph on initial component mount
    // with the default states (Methanol/Water, 60C, UNIQUAC).
    console.log("McCabe: Initial load effect triggered.");
    calculateEquilibriumCurve();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []); // Empty dependency array ensures this runs only once on mount

  // --- Data Fetching (Adapted from test/page.tsx) ---
  async function fetchCompoundDataLocal(compoundName: string): Promise<CompoundData | null> {
    if (!supabase) { throw new Error("Supabase client not initialized."); }
    console.log(`McCabe: Fetching data for ${compoundName}...`);

    try {
        const { data: compoundDbData, error: compoundError } = await supabase
            .from('compounds')
            .select('id, name, cas_number')
            .ilike('name', compoundName)
            .limit(1);

        if (compoundError) throw new Error(`Supabase compound query error: ${compoundError.message}`);
        if (!compoundDbData || compoundDbData.length === 0) throw new Error(`Compound '${compoundName}' not found.`);
        
        const compoundId = compoundDbData[0].id;
        const foundName = compoundDbData[0].name;
        const casNumber = compoundDbData[0].cas_number;

        const sourcesToTry = ['chemsep1', 'chemsep2', 'DWSIM', 'biod_db'];
        let properties: any = null;
        let foundSource: string | null = null;

        for (const source of sourcesToTry) {
            const { data: propsData, error: propsError } = await supabase
                .from('compound_properties')
                .select('properties')
                .eq('compound_id', compoundId)
                .eq('source', source)
                .single();
            if (propsError && propsError.code !== 'PGRST116') console.warn(`Supabase props query error (${source}): ${propsError.message}`);
            else if (propsData) { properties = propsData.properties; foundSource = source; break; }
        }
        if (!properties) {
            const { data: anyPropsData, error: anyPropsError } = await supabase
                .from('compound_properties').select('properties, source').eq('compound_id', compoundId).limit(1);
            if (anyPropsError) console.warn(`Supabase fallback props query error: ${anyPropsError.message}`);
            else if (anyPropsData && anyPropsData.length > 0) { properties = anyPropsData[0].properties; foundSource = anyPropsData[0].source; }
            else throw new Error(`No properties found for ${foundName}.`);
        }
        if (typeof properties !== 'object' || properties === null) throw new Error(`Invalid props format for ${foundName}.`);

        let antoine: AntoineParams | null = null;
        const antoineChemsep = properties.Antoine || properties.AntoineVaporPressure;
        if (antoineChemsep?.A && antoineChemsep.B && antoineChemsep.C) {
            antoine = {
                A: parseFloat(antoineChemsep.A?.value ?? antoineChemsep.A),
                B: parseFloat(antoineChemsep.B?.value ?? antoineChemsep.B),
                C: parseFloat(antoineChemsep.C?.value ?? antoineChemsep.C),
                Tmin_K: parseFloat(antoineChemsep.Tmin?.value ?? antoineChemsep.Tmin ?? 0),
                Tmax_K: parseFloat(antoineChemsep.Tmax?.value ?? antoineChemsep.Tmax ?? 10000),
                Units: antoineChemsep.units || 'Pa',
                EquationNo: antoineChemsep.eqno
            };
        }
        if (!antoine || isNaN(antoine.A) || isNaN(antoine.B) || isNaN(antoine.C)) throw new Error(`Failed to extract valid Antoine params for ${foundName}.`);

        let unifacGroups: UnifacGroupComposition | null = null;
        if (properties.elements_composition?.UNIFAC) {
            unifacGroups = {};
            for (const key in properties.elements_composition.UNIFAC) {
                const subgroupId = parseInt(key); const count = parseInt(properties.elements_composition.UNIFAC[key]);
                if (!isNaN(subgroupId) && !isNaN(count) && count > 0) unifacGroups[subgroupId] = count;
            }
            if (Object.keys(unifacGroups).length === 0) unifacGroups = null;
        }
        if (fluidPackage === 'unifac' && !unifacGroups) throw new Error(`UNIFAC groups required for ${foundName} with UNIFAC model.`);

        let prParams: PrPureComponentParams | null = null;
        let srkParams: SrkPureComponentParams | null = null;
        const tcPropObj = properties["Critical temperature"];
        const pcPropObj = properties["Critical pressure"];
        const omegaPropObj = properties["Acentric factor"];
        if (tcPropObj && pcPropObj && omegaPropObj) {
            const Tc_K_val = parseFloat(tcPropObj.value);
            const pcValue = parseFloat(pcPropObj.value);
            const pcUnits = String(pcPropObj.units).toLowerCase();
            const Pc_Pa_val = pcValue * (pcUnits === 'kpa' ? 1000 : pcUnits === 'bar' ? 100000 : 1);
            const omega_val = parseFloat(omegaPropObj.value);
            if (!isNaN(Tc_K_val) && !isNaN(Pc_Pa_val) && !isNaN(omega_val)) {
                prParams = { Tc_K: Tc_K_val, Pc_Pa: Pc_Pa_val, omega: omega_val };
                srkParams = { Tc_K: Tc_K_val, Pc_Pa: Pc_Pa_val, omega: omega_val };
            }
        }
        if ((fluidPackage === 'pr' && !prParams) || (fluidPackage === 'srk' && !srkParams)) throw new Error(`EOS params (Tc,Pc,omega) required for ${foundName} with ${fluidPackage.toUpperCase()} model.`);

        let uniquacParams: UniquacPureComponentParams | null = null;
        const rPropObj = properties["UNIQUAC r"] || properties["Van der Waals volume"];
        const qPropObj = properties["UNIQUAC q"] || properties["Van der Waals area"];
        if (rPropObj && qPropObj) {
            const r_val = parseFloat(rPropObj.value ?? rPropObj);
            const q_val = parseFloat(qPropObj.value ?? qPropObj);
            if (!isNaN(r_val) && !isNaN(q_val)) uniquacParams = { r: r_val, q: q_val };
        }
        if (fluidPackage === 'uniquac' && !uniquacParams) throw new Error(`UNIQUAC params (r,q) required for ${foundName}.`);
        
        let wilsonParams: WilsonPureComponentParams | null = null;
        const vLPropObj = properties["Liquid molar volume"] || properties["Molar volume"] || properties["Wilson volume"];
        if (vLPropObj) {
            const vL_val_any_unit = parseFloat(vLPropObj.value ?? vLPropObj);
            const vL_units = String(vLPropObj.units).toLowerCase() || 'cm3/mol';
            let vL_m3mol: number | undefined;
            if (!isNaN(vL_val_any_unit)) {
                if (vL_units === 'cm3/mol' || vL_units === 'cm^3/mol') vL_m3mol = vL_val_any_unit * 1e-6;
                else if (vL_units === 'm3/mol' || vL_units === 'm^3/mol') vL_m3mol = vL_val_any_unit;
                else if (vL_units === 'm3/kmol' || vL_units === 'm^3/kmol') vL_m3mol = vL_val_any_unit / 1000;
                else if (vL_units === 'l/mol' || vL_units === 'dm3/mol' || vL_units === 'dm^3/mol') vL_m3mol = vL_val_any_unit * 1e-3;
                else console.warn(`Wilson: Unknown units for V_L ('${vL_units}') for ${foundName}.`);
            }
            if (vL_m3mol !== undefined && !isNaN(vL_m3mol) && vL_m3mol > 0) wilsonParams = { V_L_m3mol: vL_m3mol };
        }
        if (fluidPackage === 'wilson' && !wilsonParams) throw new Error(`Wilson param (V_L) required for ${foundName}.`);

        return { name: foundName, antoine, unifacGroups, cas_number: casNumber, prParams, srkParams, uniquacParams, wilsonParams };
    } catch (err: any) {
        console.error(`Error fetching data for ${compoundName} in McCabe:`, err.message);
        setError(`Data fetch failed for ${compoundName}: ${err.message}`);
        return null;
    }
  }

  const calculateEquilibriumCurve = useCallback(async () => {
    console.log(`McCabe: calculateEquilibriumCurve called. Comp1: ${comp1Name}, Comp2: ${comp2Name}, Package: ${fluidPackage}`);
    setLoading(true);
    setError(null);
    setEquilibriumData(null);

    if (!comp1Name || !comp2Name) {
        setError("Please enter valid compound names.");
        setLoading(false);
        return;
    }
    if (!supabase) {
        setError("Supabase client not available.");
        setLoading(false);
        return;
    }

    // Validate temperature or pressure based on mode
    if (useTemperature && temperatureC === null) {
        setError("Temperature must be provided to calculate.");
        setLoading(false);
        return;
    }
    if (!useTemperature && pressureBar === null) {
        setError("Pressure must be provided to calculate.");
        setLoading(false);
        return;
    }

    try {
        const [data1, data2] = await Promise.all([fetchCompoundDataLocal(comp1Name), fetchCompoundDataLocal(comp2Name)]);
        if (!data1 || !data2) {
            // Error state already set by fetchCompoundDataLocal
            setLoading(false);
            return;
        }
        
        const components: CompoundData[] = [data1, data2];
        let activityParameters: any; // To hold UnifacParameters, NrtlInteractionParams, etc.

        // Fetch interaction parameters based on fluidPackage
        if (fluidPackage === 'unifac') {
            if (!data1.unifacGroups || !data2.unifacGroups) throw new Error("UNIFAC groups missing.");
            const allSubgroupIds = new Set<number>();
            components.forEach(comp => { if (comp.unifacGroups) Object.keys(comp.unifacGroups).forEach(id => allSubgroupIds.add(parseInt(id))); });
            activityParameters = await fetchUnifacInteractionParams(supabase, Array.from(allSubgroupIds));
        } else if (fluidPackage === 'nrtl') {
            if (!data1.cas_number || !data2.cas_number) throw new Error("CAS numbers required for NRTL.");
            activityParameters = await fetchNrtlParameters(supabase, data1.cas_number, data2.cas_number);
        } else if (fluidPackage === 'pr') {
            if (!data1.prParams || !data2.prParams) throw new Error("PR pure component params missing.");
            activityParameters = await fetchPrInteractionParams(supabase, data1.cas_number!, data2.cas_number!);
        } else if (fluidPackage === 'srk') {
            if (!data1.srkParams || !data2.srkParams) throw new Error("SRK pure component params missing.");
            activityParameters = await fetchSrkInteractionParams(supabase, data1.cas_number!, data2.cas_number!);
        } else if (fluidPackage === 'uniquac') {
            if (!data1.uniquacParams || !data2.uniquacParams) throw new Error("UNIQUAC pure component params missing.");
            activityParameters = await fetchUniquacInteractionParams(supabase, data1.cas_number!, data2.cas_number!);
        } else if (fluidPackage === 'wilson') {
            if (!data1.wilsonParams || !data2.wilsonParams) throw new Error("Wilson pure component params missing.");
            activityParameters = await fetchWilsonInteractionParams(supabase, data1.cas_number!, data2.cas_number!);
        } else {
            throw new Error(`Unsupported fluid package: ${fluidPackage}`);
        }

        const x_feed_values = Array.from({ length: 51 }, (_, i) => parseFloat((i * 0.02).toFixed(3))); // 51 points for x1 from 0 to 1 (step 0.02)
        const calculated_x_values: number[] = [];
        const calculated_y_values: number[] = [];

        // No top-level let fixedTempK, fixedPressurePa needed for the loop itself

        for (const x1_val of x_feed_values) {
            let resultPoint: BubbleDewResult | null = null;
            if (useTemperature) { // Constant Temperature, calculate P_bubble and y1
                // temperatureC is guaranteed not null here due to the check above
                const currentFixedTempK = temperatureC! + 273.15;
                const initialPressureGuess = (libCalculatePsat_Pa(data1.antoine!, currentFixedTempK) + libCalculatePsat_Pa(data2.antoine!, currentFixedTempK)) / 2 || 101325;
                
                if (fluidPackage === 'unifac') resultPoint = calculateBubblePressureUnifac(components, x1_val, currentFixedTempK, activityParameters as UnifacParameters, initialPressureGuess);
                else if (fluidPackage === 'nrtl') resultPoint = calculateBubblePressureNrtl(components, x1_val, currentFixedTempK, activityParameters as NrtlInteractionParams, initialPressureGuess);
                else if (fluidPackage === 'pr') resultPoint = calculateBubblePressurePr(components, x1_val, currentFixedTempK, activityParameters as PrInteractionParams, initialPressureGuess);
                else if (fluidPackage === 'srk') resultPoint = calculateBubblePressureSrk(components, x1_val, currentFixedTempK, activityParameters as SrkInteractionParams, initialPressureGuess);
                else if (fluidPackage === 'uniquac') resultPoint = calculateBubblePressureUniquac(components, x1_val, currentFixedTempK, activityParameters as LibUniquacInteractionParams, initialPressureGuess);
                else if (fluidPackage === 'wilson') resultPoint = calculateBubblePressureWilson(components, x1_val, currentFixedTempK, activityParameters as LibWilsonInteractionParams, initialPressureGuess);
            } else { // Constant Pressure, calculate T_bubble and y1
                // pressureBar is guaranteed not null here
                const currentFixedPressurePa = pressureBar! * 1e5;

                const Tbp1 = antoineBoilingPointSolverLocal(data1.antoine, currentFixedPressurePa);
                const Tbp2 = antoineBoilingPointSolverLocal(data2.antoine, currentFixedPressurePa);
                
                let initialTempGuess: number; // For the Newton solver
                if (Tbp1 !== null && Tbp2 !== null) {
                    initialTempGuess = x1_val * Tbp1 + (1 - x1_val) * Tbp2;
                } else if (Tbp1 !== null) {
                    initialTempGuess = Tbp1;
                } else if (Tbp2 !== null) {
                    initialTempGuess = Tbp2;
                } else {
                    initialTempGuess = 373.15; // Default if Antoine fails for both
                }
                initialTempGuess = Math.max(150, Math.min(initialTempGuess, 700)); // Bound it

                if (fluidPackage === 'unifac') resultPoint = calculateBubbleTemperatureUnifac(components, x1_val, currentFixedPressurePa, activityParameters as UnifacParameters, initialTempGuess);
                else if (fluidPackage === 'nrtl') resultPoint = calculateBubbleTemperatureNrtl(components, x1_val, currentFixedPressurePa, activityParameters as NrtlInteractionParams, initialTempGuess);
                else if (fluidPackage === 'pr') resultPoint = calculateBubbleTemperaturePr(components, x1_val, currentFixedPressurePa, activityParameters as PrInteractionParams, initialTempGuess);
                else if (fluidPackage === 'srk') resultPoint = calculateBubbleTemperatureSrk(components, x1_val, currentFixedPressurePa, activityParameters as SrkInteractionParams, initialTempGuess);
                else if (fluidPackage === 'uniquac') resultPoint = calculateBubbleTemperatureUniquac(components, x1_val, currentFixedPressurePa, activityParameters as LibUniquacInteractionParams, initialTempGuess);
                else if (fluidPackage === 'wilson') resultPoint = calculateBubbleTemperatureWilson(components, x1_val, currentFixedPressurePa, activityParameters as LibWilsonInteractionParams, initialTempGuess);
            }

            if (resultPoint && resultPoint.error === undefined && typeof resultPoint.comp1_equilibrium === 'number' && !isNaN(resultPoint.comp1_equilibrium)) {
                calculated_x_values.push(x1_val);
                calculated_y_values.push(resultPoint.comp1_equilibrium);
            } else {
                console.warn(`McCabe: Calculation failed for x1=${x1_val} with ${fluidPackage}. Error: ${resultPoint?.error}`);
                // Optionally push NaN or skip point
            }
        }
        
        // Add pure component points if not already perfectly covered
        if (!calculated_x_values.includes(0.0)) { calculated_x_values.unshift(0.0); calculated_y_values.unshift(0.0); }
        if (!calculated_x_values.includes(1.0)) { calculated_x_values.push(1.0); calculated_y_values.push(1.0); }
        // Sort just in case
        const sortedPairs = calculated_x_values.map((x_val, i) => ({x: x_val, y: calculated_y_values[i]}))
            .sort((a,b) => a.x - b.x);

        setEquilibriumData({ x: sortedPairs.map(p => p.x), y: sortedPairs.map(p => p.y) });

        setDisplayedComp1(comp1Name);
        setDisplayedComp2(comp2Name);
        setDisplayedUseTemp(useTemperature);
        setDisplayedFluidPackage(fluidPackage);
        if (useTemperature) { setDisplayedTemp(temperatureC); setDisplayedPressure(null); }
        else { setDisplayedTemp(null); setDisplayedPressure(pressureBar); }

    } catch (err: any) {
        console.error("McCabe: Error calculating equilibrium curve:", err);
        setError(`Calculation failed: ${err.message}`);
        setEquilibriumData(null);
    } finally {
        setLoading(false);
    }
  }, [comp1Name, comp2Name, useTemperature, temperatureC, pressureBar, fluidPackage]);

  // REMOVED: useEffect that automatically called calculateEquilibriumCurve
  // useEffect(() => {
  //   calculateEquilibriumCurve();
  // }, [calculateEquilibriumCurve]); 

  // Initial load with default values can be triggered manually by the user via the button.
  // If an initial automatic load with defaults is strictly required *once* on page load,
  // a separate useEffect with an empty dependency array could call calculateEquilibriumCurve(),
  // but the current request is to stop fetching until the button is pressed.


  const generateEChartsOptions = useCallback((xValues: number[], yValues: number[]) => {
    if (!xValues || !yValues || xValues.length === 0) {
        setEchartsOptions({});
        return;
    }
    const series: SeriesOption[] = [];
    series.push({
      name: 'Equilibrium Line',
      type: 'line',
      data: xValues.map((x, i) => [x, yValues[i]]),
      color: 'yellow', symbol: 'none', lineStyle: { width: 2.5 }, z: 5, animation: false,
    });
    series.push({
      name: 'y = x Line', type: 'line', data: [[0, 0], [1, 1]],
      color: 'white', symbol: 'none', lineStyle: { width: 1.5, type: 'dotted' }, animation: false,
    });

    const rectifyingSlope = r / (r + 1);
    const rectifyingIntercept = xd / (r + 1);
    let feedSlope: number;
    if (Math.abs(q - 1) < 1e-6) { feedSlope = Infinity; } else { feedSlope = q / (q - 1); }
    const feedIntercept = (feedSlope === Infinity) ? NaN : xf - feedSlope * xf;
    let xIntersect: number, yIntersect: number;
    if (feedSlope === Infinity) { xIntersect = xf; yIntersect = rectifyingSlope * xf + rectifyingIntercept; }
    else { xIntersect = (feedIntercept - rectifyingIntercept) / (rectifyingSlope - feedSlope); yIntersect = rectifyingSlope * xIntersect + rectifyingIntercept; }
    xIntersect = Math.max(0, Math.min(1, xIntersect));
    yIntersect = Math.max(0, Math.min(1, yIntersect));
    const strippingSlope = (xIntersect === xb) ? Infinity : (yIntersect - xb) / (xIntersect - xb);
    const strippingIntercept = (strippingSlope === Infinity) ? NaN : xb - strippingSlope * xb;

    series.push({
        name: 'Rectifying Section', type: 'line', data: [[xd, xd], [xIntersect, yIntersect]],
        color: 'orange', symbol: 'none', lineStyle: { width: 2.5 }, animation: false,
    });
    const feedLineData: EChartsPoint[] = [[xf, xf]];
    if (feedSlope === Infinity) { feedLineData.push([xf, yIntersect]); } else { feedLineData.push([xIntersect, yIntersect]); }
    series.push({
        name: 'Feed Section', type: 'line', data: feedLineData,
        color: 'red', symbol: 'none', lineStyle: { width: 2.5 }, animation: false,
    });
    series.push({
        name: 'Stripping Section', type: 'line', data: [[xIntersect, yIntersect], [xb, xb]],
        color: 'green', symbol: 'none', lineStyle: { width: 2.5 }, animation: false,
    });
    series.push({
        name: 'Key Points', type: 'scatter',
        data: [
            { value: [xd, xd], itemStyle: { color: 'orange' }, name: 'Distillate' },
            { value: [xb, xb], itemStyle: { color: 'green' }, name: 'Bottoms' },
            { value: [xf, xf], itemStyle: { color: 'red' }, name: 'Feed' }
        ],
        symbolSize: 8, label: { show: false },
        tooltip: { formatter: (params: any) => `${params.name}: (${params.value[0].toFixed(3)}, ${params.value[1].toFixed(3)})` },
        animation: false,
    });

    let stageCount = 0; let feedStageCount = 0; let currentX = xd; let currentY = xd;
    const stageLineData: EChartsPoint[] = [];
    let previousSectionIsRectifying = true;

    while (currentX > xb + 0.005 && stageCount < 25) {
        let intersectX = NaN;
        for (let i = 0; i < xValues.length - 1; i++) {
            if ((yValues[i] <= currentY && yValues[i+1] >= currentY) || (yValues[i] >= currentY && yValues[i+1] <= currentY)) {
                 if (Math.abs(yValues[i+1] - yValues[i]) < 1e-9) { intersectX = (yValues[i] === currentY) ? xValues[i] : xValues[i+1]; }
                 else { const fraction = (currentY - yValues[i]) / (yValues[i+1] - yValues[i]); intersectX = xValues[i] + fraction * (xValues[i+1] - xValues[i]); }
                 break;
            }
        }
        if (isNaN(intersectX)) {
            if (currentY <= yValues[0]) intersectX = xValues[0];
            else if (currentY >= yValues[yValues.length - 1]) intersectX = xValues[xValues.length - 1];
            else break;
        }
        intersectX = Math.max(0, Math.min(1, intersectX));
        stageLineData.push([currentX, currentY]); stageLineData.push([intersectX, currentY]); stageLineData.push([null, null]);

        let nextY: number;
        const currentSectionIsRectifying = intersectX > xIntersect;
        if (currentSectionIsRectifying) { nextY = rectifyingSlope * intersectX + rectifyingIntercept; }
        else { nextY = (strippingSlope === Infinity) ? xb : strippingSlope * intersectX + strippingIntercept; }
        nextY = Math.max(0, Math.min(1, nextY));
        if (feedStageCount === 0 && currentSectionIsRectifying !== previousSectionIsRectifying) { feedStageCount = stageCount + 1; }
        previousSectionIsRectifying = currentSectionIsRectifying;

        const isLastStage = (intersectX <= xb + 0.005) || (stageCount >= 24);
        if (isLastStage) {
            stageLineData.push([intersectX, currentY]); stageLineData.push([intersectX, intersectX]); stageLineData.push([null, null]);
            currentX = intersectX;
        } else {
            stageLineData.push([intersectX, currentY]); stageLineData.push([intersectX, nextY]); stageLineData.push([null, null]);
            currentX = intersectX; currentY = nextY;
        }
        stageCount++;
    }
    if (feedStageCount === 0 && stageCount > 0) { feedStageCount = stageCount; }
    series.push({
        name: 'Stages', type: 'line', data: stageLineData, color: 'white',
        symbol: 'none', lineStyle: { width: 2 }, connectNulls: false, legendHoverLink: false, animation: false,
    });
    setStages(stageCount); setFeedStage(feedStageCount);

    const capitalizeFirst = (str: string) => str ? str.charAt(0).toUpperCase() + str.slice(1) : '';
    const dispComp1Cap = capitalizeFirst(displayedComp1);
    const dispComp2Cap = capitalizeFirst(displayedComp2);
    const titleCondition = displayedUseTemp ? 
        (displayedTemp !== null ? `${displayedTemp.toFixed(0)} °C` : 'N/A Temp') : 
        (displayedPressure !== null ? `${displayedPressure.toFixed(2)} bar` : 'N/A Pressure');
    const titleText = displayedComp1 && displayedComp2 ?
        `McCabe-Thiele: ${dispComp1Cap}/${dispComp2Cap} at ${titleCondition}` : // Removed fluid package
        'McCabe-Thiele Diagram';
    
    // Axis label will refer to comp1
    const axisCompLabel = dispComp1Cap;

    setEchartsOptions({
        backgroundColor: '#08306b',
        title: { text: titleText, left: 'center', textStyle: { color: 'white', fontSize: 18, fontFamily: 'Merriweather Sans' } },
        grid: { left: '5%', right: '5%', bottom: '5%', top: '5%', containLabel: true },
        xAxis: {
            type: 'value', min: 0, max: 1, interval: 0.1,
            name: `Liquid Mole Fraction ${axisCompLabel} (x)`, // Updated label
            nameLocation: 'middle', nameGap: 30, nameTextStyle: { color: 'white', fontSize: 15, fontFamily: 'Merriweather Sans' },
            axisLine: { lineStyle: { color: 'white' } }, axisTick: { lineStyle: { color: 'white' }, length: 5, inside: false },
            axisLabel: { color: 'white', fontSize: 16, fontFamily: 'Merriweather Sans', formatter: '{value}' }, splitLine: { show: false },
        },
        yAxis: {
            type: 'value', min: 0, max: 1, interval: 0.1,
            name: `Vapor Mole Fraction ${axisCompLabel} (y)`, // Updated label
            nameLocation: 'middle', nameGap: 40, nameTextStyle: { color: 'white', fontSize: 15, fontFamily: 'Merriweather Sans' },
            axisLine: { lineStyle: { color: 'white' } }, axisTick: { lineStyle: { color: 'white' }, length: 5, inside: false },
            axisLabel: { color: 'white', fontSize: 16, fontFamily: 'Merriweather Sans', formatter: '{value}' }, splitLine: { show: false },
        },
        legend: {
            orient: 'vertical', right: '2%', top: 'center',
            data: series.map(s => s.name).filter(name => name !== 'Key Points' && name !== 'Stages') as string[],

            textStyle: { color: 'white', fontSize: 12, fontFamily: 'Merriweather Sans' },
            itemWidth: 12, itemHeight: 12, icon: 'rect',
        },
        tooltip: { show: true, trigger: 'axis', axisPointer: { type: 'cross' } }, // Enabled tooltip
        animationDuration: 300, animationEasing: 'cubicInOut',
        toolbox: {
            show: true, orient: 'vertical', right: 0, top: 'bottom',
            feature: { saveAsImage: { show: true, title: 'Save as Image', name: `mccabe-thiele-${displayedComp1}-${displayedComp2}`, backgroundColor: '#08306b', pixelRatio: 2 } },
            iconStyle: { borderColor: '#fff' }
        },
        series: series,
    });
  }, [xd, xb, xf, q, r, displayedComp1, displayedComp2, displayedTemp, displayedPressure, displayedUseTemp, displayedFluidPackage, buffer]);

  useEffect(() => {
    if (equilibriumData?.x && equilibriumData?.y) {
      generateEChartsOptions(equilibriumData.x, equilibriumData.y);
    } else {
        setEchartsOptions({}); setStages(null); setFeedStage(null);
    }
  }, [equilibriumData, generateEChartsOptions]);

  const handleUpdateGraphClick = () => {
    calculateEquilibriumCurve();
  };

  // Local helper for Antoine boiling point (simplified from test/page.tsx)
  const antoineBoilingPointSolverLocal = (antoineParams: AntoineParams | null, P_target: number): number | null => {
    if (!antoineParams) return null;
    // Simplified: T = B / (A - logP_base_e_or_10) - C
    // This is an approximation and might not be accurate for all Antoine forms.
    // For a more robust solution, a root finder (Newton-Raphson) is needed.
    // For now, using a rough estimate or assuming it's provided if critical.
    // This is primarily for initial guesses in bubbleT calculations.
    try {
        let logP: number;
        const P_converted_to_antoine_units = P_target / (antoineParams.Units?.toLowerCase() === 'kpa' ? 1000 : 1);

        if (antoineParams.EquationNo === 1 || antoineParams.EquationNo === '1') { // log10 form
            logP = Math.log10(P_converted_to_antoine_units);
        } else { // Assume ln form
            logP = Math.log(P_converted_to_antoine_units);
        }
        if (antoineParams.A - logP === 0) return null; // Avoid division by zero
        const T_K = antoineParams.B / (antoineParams.A - logP) - antoineParams.C;
        return (T_K > 0 && T_K < 1000) ? T_K : null; // Basic validity check
    } catch {
        return null;
    }
  };
  
  const getFeedQualityState = () => {
    if (q === 1) return "Saturated Liquid";
    if (q === 0) return "Saturated Vapor";
    if (q > 1) return "Subcooled Liquid";
    if (q > 0 && q < 1) return "Partially Vaporized";
    return "Super Heated Vapor";
  };

  const handleKeyDown = (event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key === 'Enter') {
      event.preventDefault();
      handleUpdateGraphClick(); // Trigger graph update on Enter key
    }
  };
  
  const handleSwapComponents = () => {
    const tempName = comp1Name;
    setComp1Name(comp2Name);
    setComp2Name(tempName);
  };

  const updateCompositions = (type: 'xd' | 'xf' | 'xb', value: number) => {
    // Get the current state values to check against
    const currentXd = xd;
    const currentXf = xf;
    const currentXb = xb;

    let targetValue = value; // Start with the value from the slider event

    // Check constraints and adjust the targetValue if needed (snap to boundary)
    if (type === 'xd') {
        // Prevent xd from going below or equal to xf + buffer
        if (targetValue <= currentXf + buffer) {
            targetValue = currentXf + buffer; // Snap to the boundary
        }
        // Ensure targetValue doesn't exceed absolute max
        targetValue = Math.min(targetValue, 0.99);
        setXd(targetValue);
    } else if (type === 'xf') {
        // Prevent xf from going above or equal to xd - buffer
        if (targetValue >= currentXd - buffer) {
            targetValue = currentXd - buffer; // Snap to the boundary
        }
        // Prevent xf from going below or equal to xb + buffer
        if (targetValue <= currentXb + buffer) {
            targetValue = currentXb + buffer; // Snap to the boundary
        }
         // Ensure targetValue is within absolute bounds
        targetValue = Math.max(0.01, Math.min(targetValue, 0.99));
        setXf(targetValue);
    } else if (type === 'xb') {
        // Prevent xb from going above or equal to xf - buffer
        if (targetValue >= currentXf - buffer) {
            targetValue = currentXf - buffer; // Snap to the boundary
        }
        // Ensure targetValue doesn't go below absolute min
        targetValue = Math.max(targetValue, 0.01);
        setXb(targetValue);
    }

    // Note: We directly set the state for the changed slider.
    // The constraints prevent invalid values based on the *current* state of the other sliders.
    // This avoids complex cascading updates and relies on the user not being able to move sliders past the calculated boundaries.
  };

  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8 px-8 md:px-32">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          <div className="lg:col-span-1 space-y-6">
            <Card>
              <CardContent className="space-y-6 p-4">
                <div className="space-y-4">
                  {/* Temperature/Pressure Toggle and Input */}
                  <div className="flex items-center gap-2">
                    <Tabs value={useTemperature ? "temperature" : "pressure"} onValueChange={(value) => setUseTemperature(value === "temperature")}>
                        <TabsList className="grid grid-cols-2">
                            <TabsTrigger value="temperature">Temperature</TabsTrigger> {/* Full word */}
                            <TabsTrigger value="pressure">Pressure</TabsTrigger>    {/* Full word */}
                        </TabsList>
                    </Tabs>
                    {useTemperature ? (
                        <Input
                           id="temperature" type="text"
                           value={temperatureC === null ? '' : String(temperatureC)}
                           onChange={(e) => {
                              const valStr = e.target.value;
                              if (valStr === '') { setTemperatureC(null); }
                              else {
                                  const num = parseFloat(valStr);
                                  if (!isNaN(num)) { setTemperatureC(num); }
                                  else if (valStr === '-' || (valStr.endsWith('.') && !valStr.substring(0, valStr.length -1).includes('.'))) { /* Allow partial input */ }
                                  else { setTemperatureC(null); }
                              }
                           }}
                           onKeyDown={handleKeyDown} // Attach Enter key handler
                           placeholder="Enter value" // Updated placeholder
                           required={useTemperature} disabled={!useTemperature} className="flex-1"
                        />
                    ) : (
                        <Input
                           id="pressure" type="text"
                           value={pressureBar === null ? '' : String(pressureBar)}
                           onChange={(e) => {
                              const valStr = e.target.value;
                              if (valStr === '') { setPressureBar(null); }
                              else {
                                  const num = parseFloat(valStr);
                                  if (!isNaN(num)) { setPressureBar(num); }
                                  else if (valStr === '-' || (valStr.endsWith('.') && !valStr.substring(0, valStr.length -1).includes('.'))) { /* Allow partial input */ }
                                  else { setPressureBar(null); }
                              }
                           }}
                           onKeyDown={handleKeyDown} // Attach Enter key handler
                           placeholder="Enter value" // Updated placeholder
                           required={!useTemperature} disabled={useTemperature} className="flex-1"
                        />
                    )}
                    {/* Unit Indicator */}
                    <span className="ml-1 text-sm text-muted-foreground w-8 text-left">
                        {useTemperature ? "°C" : "bar"}
                    </span>
                   </div>

                  {/* Component Inputs - Side by side with placeholders and swap button */}
                  <div className="flex items-center gap-2">
                    <Input
                      id="comp1Name"
                      value={comp1Name}
                      onChange={(e) => setComp1Name(e.target.value)}
                      onKeyDown={handleKeyDown} // Attach Enter key handler
                      placeholder="Methanol" // Specific placeholder
                      required
                      className="flex-1"
                    />
                    <Button variant="ghost" size="icon" onClick={handleSwapComponents} title="Swap Components">
                        <ArrowLeftRight className="h-4 w-4" />
                    </Button>
                    <Input
                      id="comp2Name"
                      value={comp2Name}
                      onChange={(e) => setComp2Name(e.target.value)}
                      onKeyDown={handleKeyDown} // Attach Enter key handler
                      placeholder="Water" // Specific placeholder
                      required
                      className="flex-1"
                    />
                  </div>
                  {/* Fluid Package Selector - Label side-by-side, reordered */}
                  <div className="flex items-center gap-2">
                      <Label htmlFor="fluidPackageMcCabe" className="text-sm font-medium whitespace-nowrap">Fluid Package:</Label>
                      <Select value={fluidPackage} onValueChange={(value) => setFluidPackage(value as FluidPackageTypeMcCabe)}>
                          <SelectTrigger id="fluidPackageMcCabe" className="flex-1"><SelectValue placeholder="Select fluid package" /></SelectTrigger>
                          <SelectContent>
                              <SelectItem value="uniquac">UNIQUAC</SelectItem>
                              <SelectItem value="pr">Peng-Robinson</SelectItem>
                              <SelectItem value="wilson">Wilson</SelectItem>
                              <SelectItem value="nrtl">NRTL</SelectItem>
                              <SelectItem value="srk">SRK</SelectItem>
                              <SelectItem value="unifac">UNIFAC</SelectItem>
                          </SelectContent>
                      </Select>
                  </div>
                </div>
                <Button onClick={handleUpdateGraphClick} disabled={loading} className="w-full">
                  {loading ? 'Calculating...' : 'Update Graph & Calculate XY'}
                </Button>
                {error && <p className="text-sm text-red-500 mt-2">{error}</p>}
              </CardContent>
            </Card>
            {/* Slider Card */}
            <Card>
               <CardContent className="space-y-8 pt-6 pb-6">
                 {/* xd Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="xd">
                     <span dangerouslySetInnerHTML={{ __html: `Distillate Comp (x<sub>D</sub>): ${xd.toFixed(2)}` }} />
                   </Label>
                   <Slider id="xd" min={0.01} max={0.99} step={0.01} value={[xd]} onValueChange={(value) => updateCompositions('xd', value[0])} style={{ '--primary': 'hsl(38 92% 50%)' } as React.CSSProperties}/>
                 </div>
                 {/* xf Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="xf">
                     <span dangerouslySetInnerHTML={{ __html: `Feed Comp (x<sub>F</sub>): ${xf.toFixed(2)}` }} />
                   </Label>
                   <Slider id="xf" min={0.01} max={0.99} step={0.01} value={[xf]} onValueChange={(value) => updateCompositions('xf', value[0])} style={{ '--primary': 'hsl(0 84% 60%)' } as React.CSSProperties}/>
                 </div>
                 {/* xb Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="xb">
                     <span dangerouslySetInnerHTML={{ __html: `Bottoms Comp (x<sub>B</sub>): ${xb.toFixed(2)}` }} />
                   </Label>
                   <Slider id="xb" min={0.01} max={0.99} step={0.01} value={[xb]} onValueChange={(value) => updateCompositions('xb', value[0])} style={{ '--primary': 'hsl(142 71% 45%)' } as React.CSSProperties}/>
                 </div>
                 {/* q Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="q" className="flex items-center">Feed Quality (q): {q.toFixed(2)}
                       <ShadTooltip><TooltipTrigger asChild><Button variant="ghost" size="icon" className="h-5 w-5 rounded-full"><span className="text-xs">ⓘ</span></Button></TooltipTrigger><TooltipContent><p>{getFeedQualityState()}</p></TooltipContent></ShadTooltip>
                   </Label>
                   <Slider id="q" min={-1} max={2} step={0.05} value={[q]} onValueChange={(value) => setQ(value[0])} style={{ '--primary': 'hsl(0 0% 98%)' } as React.CSSProperties}/>
                 </div>
                 {/* r Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="r">Reflux Ratio (R): {r.toFixed(2)}</Label>
                   <Slider id="r" min={0.1} max={10} step={0.05} value={[r]} onValueChange={(value) => setR(value[0])} style={{ '--primary': 'hsl(262 84% 58%)' } as React.CSSProperties}/>
                 </div>
               </CardContent>
            </Card>
          </div>
          {/* Column 2: Plot and Results Cards */}
          <div className="lg:col-span-2 space-y-6">
            <Card>
              <CardContent className="py-2">
                <div className="relative h-[500px] md:h-[600px] rounded-md" style={{ backgroundColor: '#08306b' }}>
                   {loading && ( <div className="absolute inset-0 flex items-center justify-center text-white"><div className="text-center"><div className="mb-2">Loading & Calculating XY Data...</div><div className="text-sm text-gray-300">Using {fluidPackage.toUpperCase()} model.</div></div></div> )}
                   {!loading && !equilibriumData && !error && ( <div className="absolute inset-0 flex items-center justify-center text-white">Please provide inputs and update graph.</div> )}
                   {error && !loading && ( <div className="absolute inset-0 flex items-center justify-center text-red-400">Error: {error}</div> )}
                  {!loading && equilibriumData && Object.keys(echartsOptions).length > 0 && (
                    <ReactECharts ref={echartsRef} echarts={echarts} option={echartsOptions} style={{ height: '100%', width: '100%', borderRadius: '0.375rem', overflow: 'hidden' }} notMerge={false} lazyUpdate={false} />
                  )}
                </div>
                {stages !== null && feedStage !== null && (
                  <div className="grid grid-cols-2 gap-4 text-center mt-4 pt-2">
                    <div className="p-4 bg-muted rounded-md"><p className="text-sm font-medium">Number of Stages</p><p className="text-2xl font-bold">{stages}</p></div>
                    <div className="p-4 bg-muted rounded-md"><p className="text-sm font-medium">Feed Stage</p><p className="text-2xl font-bold">{feedStage}</p></div>
                  </div>
                )}
              </CardContent>
            </Card>
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}