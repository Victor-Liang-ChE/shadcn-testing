'use client';

import React, { useState, useEffect, useCallback, useRef } from 'react';
import { createClient, SupabaseClient } from '@supabase/supabase-js';

// ECharts imports (ensure all necessary components are registered as in other pages)
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsOption, SeriesOption as EChartsSeriesOption } from 'echarts'; // Use EChartsOption
import { LineChart } from 'echarts/charts';
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  ToolboxComponent,
  DataZoomComponent,
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

echarts.use([
  TitleComponent, TooltipComponent, GridComponent, LegendComponent, ToolboxComponent, DataZoomComponent, LineChart, CanvasRenderer
]);

// Shadcn UI imports
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { TooltipProvider } from "@/components/ui/tooltip"; // Assuming you might use Shadcn tooltips
import { ArrowLeftRight, Terminal } from 'lucide-react';
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"; // Added Tabs

// VLE Calculation Libraries (ensure all are imported as in test/page.tsx)
import {
    calculatePsat_Pa as libCalculatePsat_Pa,
    calculateUnifacGamma, calculateBubbleTemperature as calculateBubbleTemperatureUnifac, calculateBubblePressure as calculateBubblePressureUnifac, fetchUnifacInteractionParams, type UnifacParameters,
} from '@/lib/vle-calculations-unifac';
import {
    fetchNrtlParameters, calculateNrtlGamma, calculateBubbleTemperatureNrtl, calculateBubblePressureNrtl, type NrtlInteractionParams
} from '@/lib/vle-calculations-nrtl';
import {
    fetchPrInteractionParams, calculateBubbleTemperaturePr, calculateBubblePressurePr, type PrInteractionParams
} from '@/lib/vle-calculations-pr';
import {
    fetchSrkInteractionParams, calculateBubbleTemperatureSrk, calculateBubblePressureSrk, type SrkInteractionParams
} from '@/lib/vle-calculations-srk';
import {
    fetchUniquacInteractionParams, calculateBubbleTemperatureUniquac, calculateBubblePressureUniquac, type UniquacInteractionParams as LibUniquacInteractionParams
} from '@/lib/vle-calculations-uniquac';
import {
    fetchWilsonInteractionParams, calculateBubbleTemperatureWilson, calculateBubblePressureWilson, type WilsonInteractionParams as LibWilsonInteractionParams
} from '@/lib/vle-calculations-wilson';

// Shared VLE Types
import type {
    AntoineParams, PrPureComponentParams, SrkPureComponentParams, UniquacPureComponentParams, WilsonPureComponentParams,
    CompoundData, BubbleDewResult, UnifacGroupComposition
} from '@/lib/vle-types';

type FluidPackageTypeAzeotrope = 'unifac' | 'nrtl' | 'pr' | 'srk' | 'uniquac' | 'wilson';
type AzeotropeScanType = 'vs_P_find_T' | 'vs_T_find_P'; // Scan P, find T_az OR Scan T, find P_az

// Supabase Client Setup
const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL || '';
const supabaseAnonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY || '';
let supabase: SupabaseClient;
if (supabaseUrl && supabaseAnonKey) {
    try {
        supabase = createClient(supabaseUrl, supabaseAnonKey);
    } catch (error) {
        console.error("Error initializing Supabase client for Azeotrope Finder:", error);
    }
} else {
    console.error("Supabase URL or Anon Key is missing for Azeotrope Finder.");
}

// Debounce function (copy from mccabe-thiele or test page)
function debounce<F extends (...args: any[]) => any>(func: F, waitFor: number) {
  let timeout: ReturnType<typeof setTimeout> | null = null;
  const debounced = (...args: Parameters<F>) => {
    if (timeout !== null) { clearTimeout(timeout); timeout = null; }
    timeout = setTimeout(() => func(...args), waitFor);
  };
  return debounced as (...args: Parameters<F>) => void;
}

// formatNumberToPrecision function (copy from mccabe-thiele)
const formatNumberToPrecision = (num: any, precision: number = 3): string => {
    if (typeof num === 'number') {
      if (num === 0) return '0';
      const fixed = num.toPrecision(precision);
      if (fixed.includes('.')) { return parseFloat(fixed).toString(); }
      return fixed;
    }
    return String(num);
};


export default function AzeotropeFinderPage() {
  // Input States
  const [comp1Name, setComp1Name] = useState('Acetone'); // Default to acetone
  const [comp2Name, setComp2Name] = useState('Water');   // Default to water
  const [fluidPackage, setFluidPackage] = useState<FluidPackageTypeAzeotrope>('uniquac'); // Default to uniquac
  const [azeotropeScanType, setAzeotropeScanType] = useState<AzeotropeScanType>('vs_P_find_T'); // Default to Scan Pressure
  const scanSteps = 50; // Always use 50 scan points, no UI to change this

  // Data & Control States
  const [azeotropeScanData, setAzeotropeScanData] = useState<{ scanVal: number, x_az: number, dependentVal: number }[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  // const [logMessages, setLogMessages] = useState<string[]>([]); // Removed logMessages

  // ECharts State
  const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
  const echartsRef = useRef<ReactECharts | null>(null);

  // Suggestion states (copy from mccabe-thiele)
  const [comp1Suggestions, setComp1Suggestions] = useState<string[]>([]);
  const [comp2Suggestions, setComp2Suggestions] = useState<string[]>([]);
  const [showComp1Suggestions, setShowComp1Suggestions] = useState(false);
  const [showComp2Suggestions, setShowComp2Suggestions] = useState(false);
  const [activeSuggestionInput, setActiveSuggestionInput] = useState<'comp1' | 'comp2' | null>(null);
  const input1Ref = useRef<HTMLInputElement>(null);
  const input2Ref = useRef<HTMLInputElement>(null);
  const suggestions1Ref = useRef<HTMLDivElement>(null);
  const suggestions2Ref = useRef<HTMLDivElement>(null);
  
  // Displayed parameters for title
  const [displayedComp1, setDisplayedComp1] = useState('');
  const [displayedComp2, setDisplayedComp2] = useState('');
  const [displayedFluidPackage, setDisplayedFluidPackage] = useState<FluidPackageTypeAzeotrope | ''>('');
  const [displayedScanType, setDisplayedScanType] = useState<AzeotropeScanType | ''>('');


  // --- Logging Hook (copy from test/page.tsx or mccabe-thiele) ---
  // Removed logging useEffect and related state

  // --- useEffect for initial graph load ---
  useEffect(() => {
    if (supabase && comp1Name && comp2Name) { // Ensure supabase is initialized and default components are set
        // console.log("AzeoFinder: Initial load effect triggered for acetone/water."); // Logging removed
        handleAzeotropeScan();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [supabase]); // Run once when supabase client is ready

  // --- Data Fetching (fetchCompoundDataLocal - adapt from mccabe-thiele/test page) ---
  // This function needs to be robust to fetch all params for all models.
  async function fetchCompoundDataLocal(compoundName: string): Promise<CompoundData | null> {
    if (!supabase) { throw new Error("Supabase client not initialized."); }
    console.log(`AzeotropeFinder: Fetching data for ${compoundName}...`);
    try {
        const { data: compoundDbData, error: compoundError } = await supabase.from('compounds').select('id, name, cas_number').ilike('name', compoundName).limit(1).single();
        if (compoundError || !compoundDbData) throw new Error(compoundError?.message || `Compound '${compoundName}' not found.`);
        
        const compoundId = compoundDbData.id;
        const foundName = compoundDbData.name;
        const casNumber = compoundDbData.cas_number;

        const sourcesToTry = ['chemsep1', 'chemsep2', 'DWSIM', 'biod_db'];
        let properties: any = null; let foundSource: string | null = null;
        for (const source of sourcesToTry) {
            const { data: propsData, error: propsError } = await supabase.from('compound_properties').select('properties').eq('compound_id', compoundId).eq('source', source).single();
            if (!propsError && propsData) { properties = propsData.properties; foundSource = source; break; }
        }
        if (!properties) {
            const { data: anyPropsData, error: anyPropsError } = await supabase.from('compound_properties').select('properties, source').eq('compound_id', compoundId).limit(1).single();
            if (anyPropsError || !anyPropsData) throw new Error(`No properties found for ${foundName}.`);
            properties = anyPropsData.properties; foundSource = anyPropsData.source;
        }
        if (typeof properties !== 'object' || properties === null) throw new Error(`Invalid props format for ${foundName} from source ${foundSource}.`);

        let antoine: AntoineParams | null = null;
        const antoineChemsep = properties.Antoine || properties.AntoineVaporPressure;
        if (antoineChemsep?.A && antoineChemsep.B && antoineChemsep.C) {
            antoine = {
                A: parseFloat(antoineChemsep.A?.value ?? antoineChemsep.A), B: parseFloat(antoineChemsep.B?.value ?? antoineChemsep.B), C: parseFloat(antoineChemsep.C?.value ?? antoineChemsep.C),
                Tmin_K: parseFloat(antoineChemsep.Tmin?.value ?? antoineChemsep.Tmin ?? 0), Tmax_K: parseFloat(antoineChemsep.Tmax?.value ?? antoineChemsep.Tmax ?? 10000),
                Units: antoineChemsep.units || 'Pa', EquationNo: antoineChemsep.eqno
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
        // Validation for UNIFAC groups will be done in handleAzeotropeScan

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
        // Validation for PR/SRK params will be done in handleAzeotropeScan

        let uniquacParams: UniquacPureComponentParams | null = null;
        const rPropObjUQ = properties["UNIQUAC r"] || properties["Van der Waals volume"]; // Renamed to avoid conflict
        const qPropObjUQ = properties["UNIQUAC q"] || properties["Van der Waals area"]; // Renamed to avoid conflict
        if (rPropObjUQ && qPropObjUQ) {
            const r_val = parseFloat(rPropObjUQ.value ?? rPropObjUQ);
            const q_val = parseFloat(qPropObjUQ.value ?? qPropObjUQ);
            if (!isNaN(r_val) && !isNaN(q_val)) uniquacParams = { r: r_val, q: q_val };
        }
        // Validation for UNIQUAC params will be done in handleAzeotropeScan
        
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
                else console.warn(`Wilson (AzeotropeFinder): Unknown units for V_L ('${vL_units}') for ${foundName}.`);
            }
            if (vL_m3mol !== undefined && !isNaN(vL_m3mol) && vL_m3mol > 0) wilsonParams = { V_L_m3mol: vL_m3mol };
        }
        // Validation for Wilson params will be done in handleAzeotropeScan

        return { name: foundName, antoine, unifacGroups, cas_number: casNumber, prParams, srkParams, uniquacParams, wilsonParams };
    } catch (err: any) {
        console.error(`Error fetching data for ${compoundName} in AzeotropeFinder:`, err.message);
        setError(`Data fetch failed for ${compoundName}: ${err.message}`);
        return null;
    }
  }

  // --- Suggestion Handlers (copy from mccabe-thiele) ---
  const fetchSuggestions = useCallback(async (inputValue: string, inputTarget: 'comp1' | 'comp2') => {
    if (!inputValue || inputValue.length < 2 || !supabase) {
      if (inputTarget === 'comp1') { setComp1Suggestions([]); setShowComp1Suggestions(false); }
      else { setComp2Suggestions([]); setShowComp2Suggestions(false); }
      return;
    }
    try {
      const { data, error } = await supabase
        .from('compounds')
        .select('name')
        .ilike('name', `${inputValue}%`) // Changed from %${inputValue}% to prioritize prefix matches
        .limit(5);
      if (error) { console.error("Azeo: Supabase suggestion fetch error:", error); if (inputTarget === 'comp1') setComp1Suggestions([]); else setComp2Suggestions([]); return; }
      const suggestions = data ? data.map(item => item.name) : [];
      if (inputTarget === 'comp1') { setComp1Suggestions(suggestions); setShowComp1Suggestions(suggestions.length > 0); }
      else { setComp2Suggestions(suggestions); setShowComp2Suggestions(suggestions.length > 0); }
    } catch (err) { console.error("Azeo: Error fetching suggestions:", err); if (inputTarget === 'comp1') setComp1Suggestions([]); else setComp2Suggestions([]); }
  }, [supabase]);

  const debouncedFetchSuggestions = useCallback(debounce(fetchSuggestions, 300), [fetchSuggestions]);

  const handleComp1NameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = e.target.value;
    setComp1Name(newValue);
    setActiveSuggestionInput('comp1');
    if (newValue.trim() === "") { setShowComp1Suggestions(false); setComp1Suggestions([]); }
    else { debouncedFetchSuggestions(newValue, 'comp1'); }
  };

  const handleComp2NameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = e.target.value;
    setComp2Name(newValue);
    setActiveSuggestionInput('comp2');
    if (newValue.trim() === "") { setShowComp2Suggestions(false); setComp2Suggestions([]); }
    else { debouncedFetchSuggestions(newValue, 'comp2'); }
  };

  const handleSuggestionClick = (suggestion: string, inputTarget: 'comp1' | 'comp2') => {
    if (inputTarget === 'comp1') {
      setComp1Name(suggestion); setShowComp1Suggestions(false); setComp1Suggestions([]); input1Ref.current?.focus();
    } else {
      setComp2Name(suggestion); setShowComp2Suggestions(false); setComp2Suggestions([]); input2Ref.current?.focus();
    }
  };

  useEffect(() => {
    function handleClickOutside(event: MouseEvent) {
      const target = event.target as Node;
      if (activeSuggestionInput === 'comp1' && suggestions1Ref.current && !suggestions1Ref.current.contains(target) && input1Ref.current && !input1Ref.current.contains(target)) {
        setShowComp1Suggestions(false);
      }
      if (activeSuggestionInput === 'comp2' && suggestions2Ref.current && !suggestions2Ref.current.contains(target) && input2Ref.current && !input2Ref.current.contains(target)) {
        setShowComp2Suggestions(false);
      }
    }
    document.addEventListener("mousedown", handleClickOutside);
    return () => { document.removeEventListener("mousedown", handleClickOutside); };
  }, [activeSuggestionInput]);

  const handleSwapComponents = () => { const temp = comp1Name; setComp1Name(comp2Name); setComp2Name(temp); };


  // --- CSV Export Function ---
  const handleExportToCSV = () => {
    if (azeotropeScanData.length === 0) {
      alert("No data to export.");
      return;
    }

    let csvContent = "";
    let scanValHeader = "";
    let xAzHeader = `Azeotropic ${displayedComp1 || 'Comp1'} Composition (x1)`;
    let dependentValHeader = "";

    if (displayedScanType === 'vs_P_find_T') {
      scanValHeader = "Pressure (bar)";
      dependentValHeader = "Azeotropic Temperature (°C)";
    } else { // vs_T_find_P
      scanValHeader = "Temperature (K)";
      dependentValHeader = "Azeotropic Pressure (bar)";
    }

    const headers = [scanValHeader, xAzHeader, dependentValHeader];
    csvContent += headers.join(",") + "\r\n";

    azeotropeScanData.forEach(row => {
      const scanValFormatted = formatNumberToPrecision(row.scanVal, 4);
      const xAzFormatted = formatNumberToPrecision(row.x_az, 4); // Increased precision for CSV
      const dependentValFormatted = (displayedScanType === 'vs_P_find_T')
        ? formatNumberToPrecision(row.dependentVal, 2) // Temp to 2 decimal places
        : formatNumberToPrecision(row.dependentVal, 4); // Pressure to 4 sig figs

      const csvRow = [scanValFormatted, xAzFormatted, dependentValFormatted].join(",");
      csvContent += csvRow + "\r\n";
    });

    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement("a");
    if (link.download !== undefined) {
      const url = URL.createObjectURL(blob);
      const filename = `azeotrope_locus_${displayedComp1 || 'Comp1'}_${displayedComp2 || 'Comp2'}_${displayedFluidPackage || 'model'}.csv`;
      link.setAttribute("href", url);
      link.setAttribute("download", filename);
      link.style.visibility = 'hidden';
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    }
  };

  // --- Bisection Solver (copy from test/page.tsx) ---
  function bisectionSolve(
    func: (x: number) => number | null, a: number, b: number, tolerance: number = 1e-5, maxIter: number = 50
  ): number | null {
    // ... (Full implementation from test/page.tsx)
    let fa = func(a); let fb = func(b);
    if (fa === null || fb === null || fa * fb >=0) {
        if (fa !== null && Math.abs(fa) < tolerance) return a;
        if (fb !== null && Math.abs(fb) < tolerance) return b;
        return null;
    }
    let c = a;
    for (let i = 0; i < maxIter; i++) {
        c = (a+b)/2; const fc = func(c);
        if (fc === null || Math.abs(fc) < tolerance || (b-a)/2 < tolerance) return (fc === null && Math.abs(func(c-tolerance) || 1) > tolerance && Math.abs(func(c+tolerance) || 1) > tolerance) ? null : c; // check if fc is null due to bad point
        if (fa * fc < 0) { b = c; } else { a = c; fa = fc; }
    }
    return null;
  }

  // --- Antoine Boiling Point Solver (for initial guesses - copy from mccabe-thiele) ---
  const antoineBoilingPointSolverLocal = (antoineParams: AntoineParams | null, P_target_Pa: number): number | null => {
    if (!antoineParams) return null;
    try {
        let logP: number;
        // Convert P_target_Pa to the pressure units expected by Antoine equation (Pa or kPa)
        const P_converted_to_antoine_units = antoineParams.Units?.toLowerCase() === 'kpa' ? P_target_Pa / 1000 : P_target_Pa;

        if (antoineParams.EquationNo === 1 || antoineParams.EquationNo === '1') { // log10 form
            logP = Math.log10(P_converted_to_antoine_units);
        } else { // Assume ln form (default for ChemSep if not specified as eqno 1)
            logP = Math.log(P_converted_to_antoine_units);
        }
        if (antoineParams.A - logP === 0) return null; 
        const T_K = antoineParams.B / (antoineParams.A - logP) - antoineParams.C;
        
        // Basic validity check against Tmin/Tmax if available and reasonable
        const Tmin = antoineParams.Tmin_K ?? 0;
        const Tmax = antoineParams.Tmax_K ?? 10000;
        if (T_K >= Tmin && T_K <= Tmax) return T_K;
        
        // Fallback if outside range but still a number (might be due to poor Tmin/Tmax in DB)
        if (T_K > 0 && T_K < 10000 && !isNaN(T_K)) {
            console.warn(`Azeo: Antoine BP for ${P_target_Pa} Pa is ${T_K}K, outside DB range [${Tmin}-${Tmax}]. Using it anyway.`);
            return T_K;
        }
        return null;
    } catch (e) {
        console.error("Azeo: Error in antoineBoilingPointSolverLocal", e);
        return null;
    }
  };


  // --- Main Azeotrope Scan Logic ---
  const handleAzeotropeScan = useCallback(async () => {
    setLoading(true); setError(null); setAzeotropeScanData([]); 
    // setLogMessages(prev => ["--- Starting Azeotrope Scan ---", ...prev]); // Logging removed

    if (!comp1Name || !comp2Name) { setError("Please enter compound names."); setLoading(false); return; }
    if (!supabase) { setError("Supabase not initialized."); setLoading(false); return; }

    try {
        const [data1, data2] = await Promise.all([fetchCompoundDataLocal(comp1Name), fetchCompoundDataLocal(comp2Name)]);
        if (!data1 || !data2) { setLoading(false); return; } // Error set in fetchCompoundDataLocal

        // Validate required parameters based on fluidPackage
        if (fluidPackage === 'unifac' && (!data1.unifacGroups || !data2.unifacGroups)) throw new Error("UNIFAC groups missing.");
        if (fluidPackage === 'pr' && (!data1.prParams || !data2.prParams)) throw new Error("PR parameters missing.");
        if (fluidPackage === 'srk' && (!data1.srkParams || !data2.srkParams)) throw new Error("SRK parameters missing.");
        if (fluidPackage === 'uniquac' && (!data1.uniquacParams || !data2.uniquacParams)) throw new Error("UNIQUAC parameters missing.");
        if (fluidPackage === 'wilson' && (!data1.wilsonParams || !data2.wilsonParams)) throw new Error("Wilson parameters missing.");
        // NRTL params are fetched as interaction params, CAS numbers checked in fetchNrtlParameters

        const components: CompoundData[] = [data1, data2];
        let activityParameters: any; // UnifacParameters | NrtlInteractionParams | etc.

        // Fetch interaction parameters
        if (fluidPackage === 'unifac') {
            const allSubgroupIds = new Set<number>();
            components.forEach(comp => { if (comp.unifacGroups) Object.keys(comp.unifacGroups).forEach(id => allSubgroupIds.add(parseInt(id))); });
            activityParameters = await fetchUnifacInteractionParams(supabase, Array.from(allSubgroupIds));
        } else if (fluidPackage === 'nrtl') {
            activityParameters = await fetchNrtlParameters(supabase, data1.cas_number!, data2.cas_number!);
        } else if (fluidPackage === 'pr') {
            activityParameters = await fetchPrInteractionParams(supabase, data1.cas_number!, data2.cas_number!);
        } else if (fluidPackage === 'srk') {
            activityParameters = await fetchSrkInteractionParams(supabase, data1.cas_number!, data2.cas_number!);
        } else if (fluidPackage === 'uniquac') {
            activityParameters = await fetchUniquacInteractionParams(supabase, data1.cas_number!, data2.cas_number!);
        } else if (fluidPackage === 'wilson') {
            activityParameters = await fetchWilsonInteractionParams(supabase, data1.cas_number!, data2.cas_number!);
        } else {
            throw new Error(`Unsupported fluid package: ${fluidPackage}`);
        }

        const scanResultsArray: { scanVal: number, x_az: number, dependentVal: number }[] = [];
        
        let currentScanStartValue: number, currentScanEndValue: number;
        if (azeotropeScanType === 'vs_P_find_T') {
            currentScanStartValue = 0.01; // bar (changed from 1 kPa)
            currentScanEndValue = 20; // bar (changed from 2000 kPa)
        } else { // vs_T_find_P
            currentScanStartValue = 230; // K
            currentScanEndValue = 550; // K
        }
        const stepValue = (currentScanEndValue - currentScanStartValue) / (scanSteps > 1 ? scanSteps - 1 : 1);

        for (let i = 0; i < scanSteps; i++) {
            const currentScanVal = currentScanStartValue + i * stepValue;
            let x_az_found: number | null = null;
            let dependent_val_found: number | null = null;

            const objectiveFunction = (x1_guess: number): number | null => {
                if (x1_guess <= 1e-5 || x1_guess >= 1.0 - 1e-5) return 1.0; // Avoid edges for stability

                let bubbleResult: BubbleDewResult | null = null;
                if (azeotropeScanType === 'vs_P_find_T') {
                    const P_system_Pa = currentScanVal * 100000; // Convert bar to Pa
                    const T_bp1_est = antoineBoilingPointSolverLocal(data1.antoine, P_system_Pa) || 350;
                    const T_bp2_est = antoineBoilingPointSolverLocal(data2.antoine, P_system_Pa) || 350;
                    const initialTempGuess = x1_guess * T_bp1_est + (1 - x1_guess) * T_bp2_est;

                    if (fluidPackage === 'unifac') bubbleResult = calculateBubbleTemperatureUnifac(components, x1_guess, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'nrtl') bubbleResult = calculateBubbleTemperatureNrtl(components, x1_guess, P_system_Pa, activityParameters, initialTempGuess);
                    // ... other fluid packages for bubbleT
                    else if (fluidPackage === 'pr') bubbleResult = calculateBubbleTemperaturePr(components, x1_guess, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'srk') bubbleResult = calculateBubbleTemperatureSrk(components, x1_guess, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'uniquac') bubbleResult = calculateBubbleTemperatureUniquac(components, x1_guess, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'wilson') bubbleResult = calculateBubbleTemperatureWilson(components, x1_guess, P_system_Pa, activityParameters, initialTempGuess);


                } else { // vs_T_find_P
                    const T_system_K = currentScanVal;
                    const P_sat1_est = libCalculatePsat_Pa(data1.antoine!, T_system_K);
                    const P_sat2_est = libCalculatePsat_Pa(data2.antoine!, T_system_K);
                    const initialPressureGuess = (P_sat1_est && P_sat2_est) ? (x1_guess * P_sat1_est + (1 - x1_guess) * P_sat2_est) : 101325;
                    
                    if (fluidPackage === 'unifac') bubbleResult = calculateBubblePressureUnifac(components, x1_guess, T_system_K, activityParameters, initialPressureGuess);
                    else if (fluidPackage === 'nrtl') bubbleResult = calculateBubblePressureNrtl(components, x1_guess, T_system_K, activityParameters, initialPressureGuess);
                    // ... other fluid packages for bubbleP
                    else if (fluidPackage === 'pr') bubbleResult = calculateBubblePressurePr(components, x1_guess, T_system_K, activityParameters, initialPressureGuess);
                    else if (fluidPackage === 'srk') bubbleResult = calculateBubblePressureSrk(components, x1_guess, T_system_K, activityParameters, initialPressureGuess);
                    else if (fluidPackage === 'uniquac') bubbleResult = calculateBubblePressureUniquac(components, x1_guess, T_system_K, activityParameters, initialPressureGuess);
                    else if (fluidPackage === 'wilson') bubbleResult = calculateBubblePressureWilson(components, x1_guess, T_system_K, activityParameters, initialPressureGuess);
                }

                if (bubbleResult && bubbleResult.error === undefined && typeof bubbleResult.comp1_equilibrium === 'number' && !isNaN(bubbleResult.comp1_equilibrium)) {
                    return x1_guess - bubbleResult.comp1_equilibrium;
                }
                return null; // Indicate failure to solver
            };

            x_az_found = bisectionSolve(objectiveFunction, 0.0001, 0.9999, 1e-4, 30);

            if (x_az_found !== null) {
                // Recalculate the dependent variable at the found x_az
                let finalBubbleResult: BubbleDewResult | null = null;
                 if (azeotropeScanType === 'vs_P_find_T') {
                    const P_system_Pa = currentScanVal * 100000; // Convert bar to Pa
                    const T_bp1_est = antoineBoilingPointSolverLocal(data1.antoine, P_system_Pa) || 350;
                    const T_bp2_est = antoineBoilingPointSolverLocal(data2.antoine, P_system_Pa) || 350;
                    const initialTempGuess = x_az_found * T_bp1_est + (1 - x_az_found) * T_bp2_est;

                    if (fluidPackage === 'unifac') finalBubbleResult = calculateBubbleTemperatureUnifac(components, x_az_found, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'nrtl') finalBubbleResult = calculateBubbleTemperatureNrtl(components, x_az_found, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'pr') finalBubbleResult = calculateBubbleTemperaturePr(components, x_az_found, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'srk') finalBubbleResult = calculateBubbleTemperatureSrk(components, x_az_found, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'uniquac') finalBubbleResult = calculateBubbleTemperatureUniquac(components, x_az_found, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'wilson') finalBubbleResult = calculateBubbleTemperatureWilson(components, x_az_found, P_system_Pa, activityParameters, initialTempGuess);
                    
                    if (finalBubbleResult && finalBubbleResult.T_K) {
                        dependent_val_found = finalBubbleResult.T_K - 273.15; // Convert to Celsius
                    }
                } else { // vs_T_find_P
                    const T_system_K = currentScanVal;
                    const P_sat1_est = libCalculatePsat_Pa(data1.antoine!, T_system_K);
                    const P_sat2_est = libCalculatePsat_Pa(data2.antoine!, T_system_K);
                    const initialPressureGuess = (P_sat1_est && P_sat2_est) ? (x_az_found * P_sat1_est + (1 - x_az_found) * P_sat2_est) : 101325;

                    if (fluidPackage === 'unifac') finalBubbleResult = calculateBubblePressureUnifac(components, x_az_found, T_system_K, activityParameters, initialPressureGuess);
                    else if (fluidPackage === 'nrtl') finalBubbleResult = calculateBubblePressureNrtl(components, x_az_found, T_system_K, activityParameters, initialPressureGuess);
                    else if (fluidPackage === 'pr') finalBubbleResult = calculateBubblePressurePr(components, x_az_found, T_system_K, activityParameters, initialPressureGuess);
                    else if (fluidPackage === 'srk') finalBubbleResult = calculateBubblePressureSrk(components, x_az_found, T_system_K, activityParameters, initialPressureGuess);
                    else if (fluidPackage === 'uniquac') finalBubbleResult = calculateBubblePressureUniquac(components, x_az_found, T_system_K, activityParameters, initialPressureGuess);
                    else if (fluidPackage === 'wilson') finalBubbleResult = calculateBubblePressureWilson(components, x_az_found, T_system_K, activityParameters, initialPressureGuess);

                    if (finalBubbleResult && finalBubbleResult.P_Pa) dependent_val_found = finalBubbleResult.P_Pa / 100000; // Convert to bar
                }
            }
            
            if (x_az_found !== null && dependent_val_found !== null) {
                scanResultsArray.push({ scanVal: currentScanVal, x_az: x_az_found, dependentVal: dependent_val_found });
                const dependentUnit = azeotropeScanType === 'vs_P_find_T' ? "°C" : "bar";
                const dependentPrecision = azeotropeScanType === 'vs_P_find_T' ? 1 : 3; 
                const dependentValFormatted = azeotropeScanType === 'vs_P_find_T' 
                    ? dependent_val_found.toFixed(dependentPrecision) 
                    : dependent_val_found.toPrecision(dependentPrecision);

                // console.log(`Azeotrope found: ScanVal=${currentScanVal.toFixed(2)}, x_az=${x_az_found.toFixed(4)}, DependentVal=${dependentValFormatted} ${dependentUnit}`); // Logging removed
            }
        }
        setAzeotropeScanData(scanResultsArray);
        setDisplayedComp1(comp1Name);
        setDisplayedComp2(comp2Name);
        setDisplayedFluidPackage(fluidPackage);
        setDisplayedScanType(azeotropeScanType);

    } catch (err: any) {
        console.error("Azeotrope Scan failed:", err);
        setError(`Scan failed: ${err.message}`);
    } finally {
        setLoading(false);
        // setLogMessages(prev => ["--- Azeotrope Scan Complete ---", ...prev]); // Logging removed
    }
  }, [comp1Name, comp2Name, fluidPackage, azeotropeScanType, supabase]);

  // --- ECharts Plotting ---
  const generateAzeotropeEChartsOptions = useCallback(() => {
    if (azeotropeScanData.length === 0) {
      setEchartsOptions({
        title: { show: false }, // Hide ECharts title when no data; JSX placeholders will be used
        grid: { show: false }, 
        xAxis: { show: false }, 
        yAxis: { show: false }
      });
      return;
    }

    const comp1DisplayCap = displayedComp1 ? displayedComp1.charAt(0).toUpperCase() + displayedComp1.slice(1) : "Comp1";
    const scanParamName = displayedScanType === 'vs_P_find_T' ? "Azeotropic Pressure (bar)" : "Azeotropic Temperature (K)"; // Changed kPa to bar
    const dependentParamName = displayedScanType === 'vs_P_find_T' ? `Azeotropic Temperature (°C)` : `Azeotropic Pressure (bar)`;
    const compositionLegendName = `Azeotropic ${comp1DisplayCap} Comp. (x₁)`;
    const compositionAxisName = `Azeotropic ${comp1DisplayCap} Composition (x₁)`;
    const titleText = displayedComp1 && displayedComp2 ? `Azeotrope Locus: ${displayedComp1}/${displayedComp2}` : 'Azeotrope Locus';


    const x_az_data = azeotropeScanData.map(d => [d.scanVal, d.x_az]);
    const dependent_val_data = azeotropeScanData.map(d => [d.scanVal, d.dependentVal]);

    const options: EChartsOption = {
      title: { 
        text: titleText, 
        left: 'center', 
        textStyle: { color: '#E5E7EB', fontSize: 18, fontFamily: 'Merriweather Sans' }, // Match McCabe-Thiele fontSize
      },
      backgroundColor: 'transparent',
      tooltip: {
        trigger: 'axis',
        axisPointer: { type: 'cross' },
        backgroundColor: '#1e293b',
        borderColor: '#3b82f6',
        textStyle: { color: '#e5e7eb', fontFamily: 'Merriweather Sans' },
        formatter: (params: any) => {
            if (!Array.isArray(params) || params.length === 0) return '';
            let tooltipHtml = `${scanParamName}: <b>${formatNumberToPrecision(params[0].axisValue, 4)}</b><br/>`;
            params.forEach((param: any) => {
                const seriesName = param.seriesName;
                const value = param.value[1];
                let precision = 3;
                if (seriesName === compositionLegendName) precision = 3;
                else if (displayedScanType === 'vs_T_find_P') precision = 4; 
                else if (displayedScanType === 'vs_P_find_T') precision = 3; 
                
                const formattedValue = (displayedScanType === 'vs_P_find_T' && seriesName === dependentParamName)
                    ? value.toFixed(1) 
                    : formatNumberToPrecision(value, precision);

                tooltipHtml += `${param.marker} ${seriesName}: <b>${formattedValue}</b><br/>`;
            });
            return tooltipHtml;
        }
      },
      // legend: { // Legend removed as per request
      //   data: [compositionLegendName, dependentParamName],
      //   bottom: 10,
      //   textStyle: { color: '#9ca3af', fontFamily: 'Merriweather Sans' }, 
      //   inactiveColor: '#4b5563'
      // },
      grid: { left: '8%', right: '10%', bottom: '8%', top: '8%', containLabel: true }, // Adjusted grid.bottom from 15% to 8%
      xAxis: {
        show: true, 
        type: 'value', name: scanParamName, nameLocation: 'middle', nameGap: 30,
        axisLabel: { color: '#e7e7eb', fontFamily: 'Merriweather Sans', fontSize: 16 }, // Changed color, Added fontSize
        nameTextStyle: { color: '#e7e7eb', fontSize: 15, fontFamily: 'Merriweather Sans' }, // Added fontSize
        axisLine: { lineStyle: { color: '#4b5563', width: 2.5 } }, // Added width
        splitLine: { show: false }, // Removed grid lines
        scale: true, 
      },
      yAxis: [
        {
          show: true, 
          type: 'value', name: compositionAxisName, nameLocation: 'middle', nameGap: 50, min: 0, max: 1,
          position: 'left', 
          axisLabel: { color: '#e7e7eb', formatter: (v: number) => formatNumberToPrecision(v,3), fontFamily: 'Merriweather Sans', fontSize: 16 }, // Changed color, Added fontSize
          nameTextStyle: { color: '#e7e7eb', fontSize: 15, fontFamily: 'Merriweather Sans' }, // Added fontSize
          axisLine: { lineStyle: { color: '#22c55e', width: 2.5 } }, // Changed color to green for composition axis, added width
          splitLine: { show: false }, // Removed grid lines
        },
        {
          show: true, 
          type: 'value', name: dependentParamName, nameLocation: 'middle', nameGap: 60, position: 'right',
          axisLabel: { 
            color: '#e7e7eb', // Changed color
            formatter: (v: number) => (displayedScanType === 'vs_P_find_T') ? v.toFixed(1) : formatNumberToPrecision(v, 3), 
            fontFamily: 'Merriweather Sans',
            fontSize: 16 // Added fontSize
          },
          nameTextStyle: { color: '#e7e7eb', fontSize: 15, fontFamily: 'Merriweather Sans' }, // Added fontSize
          axisLine: { lineStyle: { color: '#f59e0b', width: 2.5 } }, // Added width
          splitLine: { show: false }, // Removed grid lines
          scale: true,
        }
      ],
      series: [
        { name: compositionLegendName, type: 'line', yAxisIndex: 0, data: x_az_data, smooth: true, lineStyle: { color: '#22c55e', width: 2.5 }, itemStyle: { color: '#22c55e' }, symbolSize: 6, symbol: 'circle' }, // Changed color to green
        { name: dependentParamName, type: 'line', yAxisIndex: 1, data: dependent_val_data, smooth: true, lineStyle: { color: '#f59e0b', width: 2.5 }, itemStyle: { color: '#f59e0b' }, symbolSize: 6, symbol: 'triangle' }
      ],
      // dataZoom removed
      toolbox: {
            show: true, 
            orient: 'horizontal', // Changed from vertical
            right: 10, 
            bottom: 10, // Changed from top: 'center'
            feature: { saveAsImage: { show: true, title: 'Save as Image', backgroundColor: '#0f172a' } },
            iconStyle: { borderColor: '#9ca3af' }
      },
    };
    setEchartsOptions(options);
  }, [azeotropeScanData, displayedComp1, displayedComp2, displayedFluidPackage, displayedScanType]);

  useEffect(() => {
    generateAzeotropeEChartsOptions();
  }, [azeotropeScanData, generateAzeotropeEChartsOptions]);


  // --- JSX Structure ---
  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Left Column: Inputs */}
          <div className="lg:col-span-1 space-y-6">
            <Card>
              {/* <CardHeader>{/* <CardTitle>Azeotrope Finder Inputs</CardTitle> * /}</CardHeader> */} {/* Removed CardHeader */}
              <CardContent className="p-4 space-y-4"> {/* Added p-4 for consistent padding */}
                <Tabs value={azeotropeScanType} onValueChange={(value) => setAzeotropeScanType(value as AzeotropeScanType)} className="mb-4">
                    <TabsList className="grid w-full grid-cols-2">
                        <TabsTrigger value="vs_P_find_T">Scan Pressure</TabsTrigger>
                        <TabsTrigger value="vs_T_find_P">Scan Temperature</TabsTrigger>
                    </TabsList>
                </Tabs>
                {/* Component Inputs with Suggestions (similar to McCabe-Thiele) */}
                <div className="flex items-center gap-2">
                    <div className="relative flex-1">
                        {/* <Label htmlFor="comp1NameAzeo">Component 1</Label> */}
                        <Input ref={input1Ref} id="comp1NameAzeo" value={comp1Name} onChange={handleComp1NameChange} 
                               onFocus={() => { setActiveSuggestionInput('comp1'); if (comp1Name.trim()) debouncedFetchSuggestions(comp1Name, 'comp1');}}
                               placeholder="e.g., Ethanol" autoComplete="off" />
                        {showComp1Suggestions && comp1Suggestions.length > 0 && (
                            <div ref={suggestions1Ref} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1 max-h-48 overflow-y-auto">
                                {comp1Suggestions.map((s, i) => <div key={i} onClick={() => handleSuggestionClick(s, 'comp1')} className="px-3 py-2 hover:bg-accent cursor-pointer text-sm">{s}</div>)}
                            </div>
                        )}
                    </div>
                    <Button variant="ghost" size="icon" onClick={handleSwapComponents} title="Swap Components" className="self-end mb-1"><ArrowLeftRight className="h-4 w-4" /></Button>
                    <div className="relative flex-1">
                        {/* <Label htmlFor="comp2NameAzeo">Component 2</Label> */}
                        <Input ref={input2Ref} id="comp2NameAzeo" value={comp2Name} onChange={handleComp2NameChange} 
                               onFocus={() => { setActiveSuggestionInput('comp2'); if (comp2Name.trim()) debouncedFetchSuggestions(comp2Name, 'comp2');}}
                               placeholder="e.g., Water" autoComplete="off" />
                        {showComp2Suggestions && comp2Suggestions.length > 0 && (
                            <div ref={suggestions2Ref} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1 max-h-48 overflow-y-auto">
                                {comp2Suggestions.map((s, i) => <div key={i} onClick={() => handleSuggestionClick(s, 'comp2')} className="px-3 py-2 hover:bg-accent cursor-pointer text-sm">{s}</div>)}
                            </div>
                        )}
                    </div>
                </div>

                <div className="flex items-center gap-2">
                  <Label htmlFor="fluidPackageAzeo" className="whitespace-nowrap">Fluid Package:</Label>
                  <Select value={fluidPackage} onValueChange={(value) => setFluidPackage(value as FluidPackageTypeAzeotrope)}>
                    <SelectTrigger id="fluidPackageAzeo" className="flex-1"><SelectValue /></SelectTrigger>
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
                {/* Scan Type Select removed */}
                {/* 
                <div>
                  <Label htmlFor="azeotropeScanTypeAzeo">Scan Type:</Label>
                  <Select value={azeotropeScanType} onValueChange={(value) => setAzeotropeScanType(value as AzeotropeScanType)}>
                    <SelectTrigger id="azeotropeScanTypeAzeo"><SelectValue /></SelectTrigger>
                    <SelectContent>
                      <SelectItem value="vs_P_find_T">Scan Pressure (find x_az, T_az)</SelectItem>
                      <SelectItem value="vs_T_find_P">Scan Temperature (find x_az, P_az)</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
                */}
                {/* Scan Steps input removed - defaults to 50 */}
                {/* 
                <div className="flex items-center gap-2">
                  <Label htmlFor="scanStepsAzeo" className="whitespace-nowrap">Scan Points:</Label>
                  <Input id="scanStepsAzeo" type="number" value={scanSteps} onChange={(e) => setScanSteps(Math.max(5, parseInt(e.target.value)))} min="5" max="200" className="flex-1"/>
                </div>
                */}
                <Button onClick={handleAzeotropeScan} disabled={loading} className="w-full">
                  {loading ? 'Scanning...' : 'Find Azeotropes'}
                </Button>
                <Button 
                  onClick={handleExportToCSV} 
                  disabled={loading || azeotropeScanData.length === 0} 
                  className="w-full mt-2"
                  variant="outline"
                >
                  Export to CSV
                </Button>
                {error && <Alert variant="destructive" className="mt-2"><Terminal className="h-4 w-4" /><AlertTitle>Error</AlertTitle><AlertDescription>{error}</AlertDescription></Alert>}
              </CardContent>
            </Card>
            {/* Log Card Removed */}
            {/* 
            <Card>
                <CardHeader><CardTitle>Calculation Log</CardTitle></CardHeader>
                <CardContent><div className="h-48 overflow-y-auto bg-muted p-2 rounded text-xs font-mono">{logMessages.length === 0 ? <p>No messages yet...</p> : logMessages.map((msg, index) => <div key={index}>{msg}</div>)}</div></CardContent>
            </Card>
            */}
          </div>

          {/* Right Column: Plot */}
          <div className="lg:col-span-2 space-y-6">
            <Card>
              <CardContent className="py-2"> {/* Minimal padding for graph card content */}
                <div className="relative h-[600px] md:h-[700px] rounded-md" style={{ backgroundColor: '#08306b' }}> {/* Dark background for plot area */}
                  {loading && (<div className="absolute inset-0 flex items-center justify-center text-white"><div className="text-center"><div className="mb-2">Scanning for Azeotropes...</div><div className="text-sm text-gray-300">Using {fluidPackage.toUpperCase()} model.</div></div></div>)}
                  {!loading && !displayedComp1 && !error && (<div className="absolute inset-0 flex items-center justify-center text-white">Please provide inputs and find azeotropes.</div>)}
                  {error && !loading && (<div className="absolute inset-0 flex items-center justify-center text-red-400">Error: {error}</div>)}
                  {!loading && Object.keys(echartsOptions).length > 0 && (
                    <ReactECharts ref={echartsRef} echarts={echarts} option={echartsOptions} style={{ height: '100%', width: '100%' }} notMerge={true} lazyUpdate={true} />
                  )}
                </div>
                 <p className="text-xs text-muted-foreground mt-2 text-center">
                    Note: Accuracy depends on database parameters and model limitations. Results are from a wide scan range.
                </p>
              </CardContent>
            </Card>
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}
