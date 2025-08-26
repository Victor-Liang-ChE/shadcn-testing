'use client';

import React, { useState, useEffect, useCallback, useRef } from 'react';
import { createClient, SupabaseClient } from '@supabase/supabase-js';
import { useTheme } from "next-themes";

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
  calculatePsat_Pa as libCalculatePsat_Pa,
  // UNIFAC
  calculateUnifacGamma,
  calculateBubbleTemperatureUnifac,
  calculateBubblePressureUnifac,
  fetchUnifacInteractionParams,
  type UnifacParameters,
  // NRTL
  fetchNrtlParameters,
  calculateNrtlGamma,
  calculateBubbleTemperatureNrtl,
  calculateBubblePressureNrtl,
  type NrtlInteractionParams,
  // Peng–Robinson
  fetchPrInteractionParams,
  calculateBubbleTemperaturePr,
  calculateBubblePressurePr,
  type PrInteractionParams,
  // SRK
  fetchSrkInteractionParams,
  calculateBubbleTemperatureSrk,
  calculateBubblePressureSrk,
  type SrkInteractionParams,
  // UNIQUAC
  fetchUniquacInteractionParams,
  calculateBubbleTemperatureUniquac,
  calculateBubblePressureUniquac,
  type UniquacInteractionParams as LibUniquacInteractionParams,
  // Wilson
  fetchWilsonInteractionParams,
  calculateBubbleTemperatureWilson,
  calculateBubblePressureWilson,
  type WilsonInteractionParams as LibWilsonInteractionParams
} from '@/lib/vle-calculations';

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

// Debounce function
function debounce<F extends (...args: any[]) => any>(func: F, waitFor: number) {
  let timeout: ReturnType<typeof setTimeout> | null = null;

  const debounced = (...args: Parameters<F>) => {
    if (timeout !== null) {
      clearTimeout(timeout);
      timeout = null;
    }
    timeout = setTimeout(() => func(...args), waitFor);
  };

  return debounced as (...args: Parameters<F>) => void;
}

// Function to fetch CAS number by name from the new compound_properties table
const fetchCasNumberByName = async (supabaseClient: SupabaseClient, name: string): Promise<string | null> => {
    if (!name || !name.trim()) {
        console.error("fetchCasNumberByName (McCabe-Thiele): Name is empty.");
        return null;
    }
    const trimmedName = name.trim();
    const { data, error } = await supabaseClient
        .from('compound_properties')
        .select('properties')
        .ilike('name', trimmedName)
        .limit(1)
        .single();

    if (error) {
        console.error(`fetchCasNumberByName (McCabe-Thiele): Error fetching CAS for "${trimmedName}":`, error);
        return null;
    }

    if (!data || !data.properties?.CAS?.value) {
        console.warn(`fetchCasNumberByName (McCabe-Thiele): No CAS number found for "${trimmedName}". Data received:`, data);
        return null;
    }

    return data.properties.CAS.value;
};

export default function McCabeThielePage() {
  const { resolvedTheme } = useTheme(); // Get the resolved theme ('light' or 'dark')
  
  // Input States - Updated defaults for initial load
  const [comp1Name, setComp1Name] = useState('methanol');
  const [comp2Name, setComp2Name] = useState('water');
  const [temperatureC, setTemperatureC] = useState<number | null>(60);
  const [pressureBar, setPressureBar] = useState<number | null>(1); // Default, but will be overridden by useTemperature=true
  
  // New string states for input fields
  const [temperatureInput, setTemperatureInput] = useState<string>(temperatureC !== null ? String(temperatureC) : '');
  const [pressureInput, setPressureInput] = useState<string>(pressureBar !== null ? String(pressureBar) : '');
  const [tempMax, setTempMax] = useState<string>('100');
  const [pressureMax, setPressureMax] = useState<string>('10');
  const [qMinSlider, setQMinSlider] = useState<string>('-1');
  const [qMaxSlider, setQMaxSlider] = useState<string>('2');
  const [rMinSlider, setRMinSlider] = useState<string>('0.1');
  const [rMaxSlider, setRMaxSlider] = useState<string>('10');

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

  // --- NEW: Caching states for component data and parameters ---
  const [comp1Data, setComp1Data] = useState<CompoundData | null>(null);
  const [comp2Data, setComp2Data] = useState<CompoundData | null>(null);
  const [interactionParams, setInteractionParams] = useState<any>(null);

  // Result States
  const [stages, setStages] = useState<number | null>(null);
  const [feedStage, setFeedStage] = useState<number | null>(null);
  const [rMin, setRMin] = useState<number | null>(null);
  const [murphreeEfficiency, setMurphreeEfficiency] = useState(1.0);

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

  // Suggestion states
  const [comp1Suggestions, setComp1Suggestions] = useState<string[]>([]);
  const [comp2Suggestions, setComp2Suggestions] = useState<string[]>([]);
  const [showComp1Suggestions, setShowComp1Suggestions] = useState(false);
  const [showComp2Suggestions, setShowComp2Suggestions] = useState(false);
  const [activeSuggestionInput, setActiveSuggestionInput] = useState<'comp1' | 'comp2' | null>(null);

  const input1Ref = useRef<HTMLInputElement>(null);
  const input2Ref = useRef<HTMLInputElement>(null);
  const activeComponentRef = useRef<HTMLDivElement>(null);

  // Flag to trigger automatic graph regeneration after UI actions like swap, toggle, or fluid-package change
  const [autoGeneratePending, setAutoGeneratePending] = useState(false);

  // --- useEffect for initial graph load ---
  useEffect(() => {
    // Automatically calculate and display the graph on initial component mount
    // with the default states (Methanol/Water, 60C, UNIQUAC).
    if (supabase) { // Ensure supabase is initialized before first calculation
        calculateEquilibriumCurve();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [supabase]); // Add supabase as a dependency if it can be initially undefined

  // --- Data Fetching (Updated for new compound_properties table structure) ---
  async function fetchCompoundDataLocal(compoundName: string): Promise<CompoundData | null> {
    if (!supabase) { throw new Error("Supabase client not initialized."); }
    try {
        // Fetch compound data from the new compound_properties table
        const { data: compoundDbData, error: compoundError } = await supabase
            .from('compound_properties')
            .select('name, properties')
            .ilike('name', compoundName)
            .limit(1)
            .single();

        if (compoundError) throw new Error(`Supabase compound query error: ${compoundError.message}`);
        if (!compoundDbData || !compoundDbData.properties) throw new Error(`Compound '${compoundName}' not found.`);
        
        const foundName = compoundDbData.name;
        const properties = compoundDbData.properties;
        
        // Extract CAS number from properties
        const casNumber = properties.CAS?.value || null;

        if (typeof properties !== 'object' || properties === null) throw new Error(`Invalid props format for ${foundName}.`);

        let antoine: AntoineParams | null = null;
        const antoineChemsep = properties.Antoine || properties.AntoineVaporPressure || properties.VaporPressure;
        if (antoineChemsep?.A && antoineChemsep.B && antoineChemsep.C) {
            antoine = {
                A: parseFloat(antoineChemsep.A?.value ?? antoineChemsep.A),
                B: parseFloat(antoineChemsep.B?.value ?? antoineChemsep.B),
                C: parseFloat(antoineChemsep.C?.value ?? antoineChemsep.C),
                Tmin_K: parseFloat(antoineChemsep.Tmin?.value ?? antoineChemsep.Tmin ?? 0),
                Tmax_K: parseFloat(antoineChemsep.Tmax?.value ?? antoineChemsep.Tmax ?? 10000),
                Units: antoineChemsep.units || 'Pa',
                EquationNo: antoineChemsep.eqno?.value ?? antoineChemsep.eqno
            };
        }
        if (!antoine || isNaN(antoine.A) || isNaN(antoine.B) || isNaN(antoine.C)) throw new Error(`Failed to extract valid Antoine params for ${foundName}.`);

    let unifacGroups: UnifacGroupComposition | null = null;
    // UNIFAC groups may be stored either as an array of {id,value} or as an object map {"1": "2", ...}
    const unifacRaw = properties.UnifacVLE || properties.UnifacVLE;
    if (unifacRaw) {
      unifacGroups = {};
      try {
        const groupData = unifacRaw.group ?? unifacRaw;
        if (Array.isArray(groupData)) {
          for (const groupItem of groupData) {
            const subgroupId = parseInt(String(groupItem?.id ?? groupItem?.group));
            const count = parseInt(String(groupItem?.value ?? groupItem?.count ?? 0));
            if (!isNaN(subgroupId) && !isNaN(count) && count > 0) {
              unifacGroups[subgroupId] = count;
            }
          }
        } else if (groupData && typeof groupData === 'object') {
          for (const [key, val] of Object.entries(groupData)) {
            const subgroupId = parseInt(key);
            let count = NaN;
            if (val && typeof val === 'object') {
              count = parseInt(String((val as any).value ?? (val as any).count ?? NaN));
            } else {
              count = parseInt(String(val as any));
            }
            if (!isNaN(subgroupId) && !isNaN(count) && count > 0) {
              unifacGroups[subgroupId] = count;
            }
          }
        } else {
          // Unexpected format: leave unifacGroups as null
          console.warn(`UNIFAC group data in unexpected format for ${foundName}.`, groupData);
          unifacGroups = null;
        }
      } catch (err) {
        console.warn(`Error parsing UNIFAC groups for ${foundName}:`, err);
        unifacGroups = null;
      }
      if (unifacGroups && Object.keys(unifacGroups).length === 0) unifacGroups = null;
    }
        if (fluidPackage === 'unifac' && !unifacGroups) throw new Error(`UNIFAC groups required for ${foundName} with UNIFAC model.`);

        let prParams: PrPureComponentParams | null = null;
        let srkParams: SrkPureComponentParams | null = null;
        const tcPropObj = properties["CriticalTemperature"] || properties["Critical temperature"];
        const pcPropObj = properties["CriticalPressure"] || properties["Critical pressure"];
        const omegaPropObj = properties["AcentricityFactor"] || properties["Acentric factor"];
        if (tcPropObj && pcPropObj && omegaPropObj) {
            const Tc_K_val = parseFloat(tcPropObj.value);
            const pcValue = parseFloat(pcPropObj.value);
            const pcUnits = String(pcPropObj.units).toLowerCase();
            let Pc_Pa_val: number;
            if (pcUnits === 'pa') Pc_Pa_val = pcValue;
            else if (pcUnits === 'kpa') Pc_Pa_val = pcValue * 1e3;
            else if (pcUnits === 'mpa') Pc_Pa_val = pcValue * 1e6;
            else if (pcUnits === 'bar') Pc_Pa_val = pcValue * 1e5;
            else {
                Pc_Pa_val = pcValue; // Assume Pa if unknown; log once
                console.warn(`Unknown Pc units ('${pcUnits}') for ${foundName}. Assuming Pa.`);
            }
            const omega_val = parseFloat(omegaPropObj.value);
            if (!isNaN(Tc_K_val) && !isNaN(Pc_Pa_val) && !isNaN(omega_val)) {
                prParams = { Tc_K: Tc_K_val, Pc_Pa: Pc_Pa_val, omega: omega_val };
                srkParams = { Tc_K: Tc_K_val, Pc_Pa: Pc_Pa_val, omega: omega_val };
            }
        }
        if ((fluidPackage === 'pr' && !prParams) || (fluidPackage === 'srk' && !srkParams)) throw new Error(`EOS params (Tc,Pc,omega) required for ${foundName} with ${fluidPackage.toUpperCase()} model.`);

        let uniquacParams: UniquacPureComponentParams | null = null;
        const rPropObj = properties["UniquacR"] || properties["UNIQUAC r"];
        const qPropObj = properties["UniquacQ"] || properties["UNIQUAC q"];
        if (rPropObj && qPropObj) {
            const r_val = parseFloat(rPropObj.value ?? rPropObj);
            const q_val = parseFloat(qPropObj.value ?? qPropObj);
            if (!isNaN(r_val) && !isNaN(q_val)) uniquacParams = { r: r_val, q: q_val };
        }
        if (fluidPackage === 'uniquac' && !uniquacParams) throw new Error(`UNIQUAC params (r,q) required for ${foundName}.`);
        
        let wilsonParams: WilsonPureComponentParams | null = null;
        const vLPropObj = properties["WilsonVolume"] || properties["Wilson volume"] || properties["LiquidVolumeAtNormalBoilingPoint"];
        if (vLPropObj) {
            const vL_val_any_unit = parseFloat(vLPropObj.value ?? vLPropObj);
            const vL_units = String(vLPropObj.units || 'm3/kmol').toLowerCase();
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

  // --- Fetch Suggestions ---
  const fetchSuggestions = useCallback(async (inputValue: string, inputTarget: 'comp1' | 'comp2') => {
    if (!inputValue || inputValue.length < 2 || !supabase) {
      if (inputTarget === 'comp1') {
        setComp1Suggestions([]);
        setShowComp1Suggestions(false);
      } else {
        setComp2Suggestions([]);
        setShowComp2Suggestions(false);
      }
      return;
    }

    try {
      const { data, error } = await supabase
        .from('compound_properties')
        .select('name')
        .ilike('name', `${inputValue}%`) // Changed from %${inputValue}% to prioritize prefix matches
        .limit(5);

      if (error) {
        console.error("Supabase suggestion fetch error:", error);
        if (inputTarget === 'comp1') setComp1Suggestions([]); else setComp2Suggestions([]);
        return;
      }

      const suggestions = data ? data.map(item => item.name) : [];
      if (inputTarget === 'comp1') {
        setComp1Suggestions(suggestions);
        setShowComp1Suggestions(suggestions.length > 0);
      } else {
        setComp2Suggestions(suggestions);
        setShowComp2Suggestions(suggestions.length > 0);
      }
    } catch (err) {
      console.error("Error fetching suggestions:", err);
      if (inputTarget === 'comp1') setComp1Suggestions([]); else setComp2Suggestions([]);
    }
  }, [supabase]);

  const handleComp1NameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = e.target.value;
    setComp1Name(newValue);
    setActiveSuggestionInput('comp1');
    if (newValue.trim() === "") {
      setShowComp1Suggestions(false);
      setComp1Suggestions([]);
    } else {
      fetchSuggestions(newValue, 'comp1');
    }
  };

  const handleComp2NameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = e.target.value;
    setComp2Name(newValue);
    setActiveSuggestionInput('comp2');
    if (newValue.trim() === "") {
      setShowComp2Suggestions(false);
      setComp2Suggestions([]);
    } else {
      fetchSuggestions(newValue, 'comp2');
    }
  };

  const handleSuggestionClick = (suggestion: string, inputTarget: 'comp1' | 'comp2') => {
    if (inputTarget === 'comp1') {
      setComp1Name(suggestion);
      setShowComp1Suggestions(false);
      setComp1Suggestions([]);
    } else {
      setComp2Name(suggestion);
      setShowComp2Suggestions(false);
      setComp2Suggestions([]);
    }
    // setActiveSuggestionInput(null); // Don't reset active input here, let blur handle it or next focus
  };
  
  // Effect to handle clicks outside of the suggestion boxes
  useEffect(() => {
    function handleClickOutside(event: MouseEvent) {
      const target = event.target as Node;
      
      // If there's an active component with suggestions showing, check if click is outside
      if (activeComponentRef.current && !activeComponentRef.current.contains(target)) {
        setShowComp1Suggestions(false);
        setShowComp2Suggestions(false);
      }
    }
    
    document.addEventListener("mousedown", handleClickOutside);
    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, [activeSuggestionInput]);

  const calculateEquilibriumCurve = useCallback(async (showLoading: boolean = true, pointsCount: number = 101) => {
    if (showLoading) {
      setLoading(true);
      // setEquilibriumData(null); // This line is removed to keep the old graph visible.
    }
    setError(null);

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
        let data1 = comp1Data;
        let data2 = comp2Data;
        let activityParameters = interactionParams;

        // Check if we need to re-fetch data
        const needsFetching = !data1 || !data2 || !activityParameters || data1.name !== comp1Name || data2.name !== comp2Name || displayedFluidPackage !== fluidPackage;

        if (needsFetching) {
            const [fetchedData1, fetchedData2] = await Promise.all([fetchCompoundDataLocal(comp1Name), fetchCompoundDataLocal(comp2Name)]);
            if (!fetchedData1 || !fetchedData2) {
                setLoading(false);
                return;
            }
            
            const fetchedComponents: CompoundData[] = [fetchedData1, fetchedData2];
            let fetchedActivityParameters: any;

            if (fluidPackage === 'unifac') {
                if (!fetchedData1.unifacGroups || !fetchedData2.unifacGroups) throw new Error("UNIFAC groups missing.");
                const allSubgroupIds = new Set<number>();
                fetchedComponents.forEach(comp => { if (comp.unifacGroups) Object.keys(comp.unifacGroups).forEach(id => allSubgroupIds.add(parseInt(id))); });
                fetchedActivityParameters = await fetchUnifacInteractionParams(supabase, Array.from(allSubgroupIds));
            } else if (fluidPackage === 'nrtl') {
                if (!fetchedData1.cas_number || !fetchedData2.cas_number) throw new Error("CAS numbers required for NRTL.");
                fetchedActivityParameters = await fetchNrtlParameters(supabase, fetchedData1.cas_number, fetchedData2.cas_number);
            } else if (fluidPackage === 'pr') {
                if (!fetchedData1.prParams || !fetchedData2.prParams) throw new Error("PR pure component params missing.");
                fetchedActivityParameters = await fetchPrInteractionParams(supabase, fetchedData1.cas_number!, fetchedData2.cas_number!);
            } else if (fluidPackage === 'srk') {
                if (!fetchedData1.srkParams || !fetchedData2.srkParams) throw new Error("SRK pure component params missing.");
                fetchedActivityParameters = await fetchSrkInteractionParams(supabase, fetchedData1.cas_number!, fetchedData2.cas_number!);
            } else if (fluidPackage === 'uniquac') {
                if (!fetchedData1.uniquacParams || !fetchedData2.uniquacParams) throw new Error("UNIQUAC pure component params missing.");
                fetchedActivityParameters = await fetchUniquacInteractionParams(supabase, fetchedData1.cas_number!, fetchedData2.cas_number!);
            } else if (fluidPackage === 'wilson') {
                if (!fetchedData1.wilsonParams || !fetchedData2.wilsonParams) throw new Error("Wilson pure component params missing.");
                fetchedActivityParameters = await fetchWilsonInteractionParams(supabase, fetchedData1.cas_number!, fetchedData2.cas_number!);
            } else {
                throw new Error(`Unsupported fluid package: ${fluidPackage}`);
            }
            
            // Update cache and use the new data for this run
            setComp1Data(fetchedData1);
            setComp2Data(fetchedData2);
            setInteractionParams(fetchedActivityParameters);
            data1 = fetchedData1;
            data2 = fetchedData2;
            activityParameters = fetchedActivityParameters;
        }

        if (!data1 || !data2 || !activityParameters) {
             setError("Component data or parameters not available for calculation.");
             if (showLoading) setLoading(false);
             return;
        }

        const components: CompoundData[] = [data1, data2];

        const stepSize = 1 / (pointsCount - 1);
        const x_feed_values = Array.from({ length: pointsCount }, (_, i) => parseFloat((i * stepSize).toFixed(4)));
        const calculated_x_values: number[] = [0.0];
        const calculated_y_values: number[] = [0.0];

        // Loop over interior points only to avoid issues at pure component conditions
        for (const x1_val of x_feed_values.slice(1, -1)) {
            let resultPoint: BubbleDewResult | null = null;
            if (useTemperature) {
                const currentFixedTempK = temperatureC! + 273.15;
                const initialPressureGuess = (libCalculatePsat_Pa(data1.antoine!, currentFixedTempK) + libCalculatePsat_Pa(data2.antoine!, currentFixedTempK)) / 2 || 101325;
                
                if (fluidPackage === 'unifac') resultPoint = calculateBubblePressureUnifac(components, x1_val, currentFixedTempK, activityParameters as UnifacParameters);
                else if (fluidPackage === 'nrtl') resultPoint = calculateBubblePressureNrtl(components, x1_val, currentFixedTempK, activityParameters as NrtlInteractionParams);
                else if (fluidPackage === 'pr') resultPoint = calculateBubblePressurePr(components, x1_val, currentFixedTempK, activityParameters as PrInteractionParams, initialPressureGuess);
                else if (fluidPackage === 'srk') resultPoint = calculateBubblePressureSrk(components, x1_val, currentFixedTempK, activityParameters as SrkInteractionParams, initialPressureGuess);
                else if (fluidPackage === 'uniquac') resultPoint = calculateBubblePressureUniquac(components, x1_val, currentFixedTempK, activityParameters as LibUniquacInteractionParams);
                else if (fluidPackage === 'wilson') resultPoint = calculateBubblePressureWilson(components, x1_val, currentFixedTempK, activityParameters as LibWilsonInteractionParams);
            } else {
                const currentFixedPressurePa = pressureBar! * 1e5;

                const Tbp1 = antoineBoilingPointSolverLocal(data1.antoine, currentFixedPressurePa);
                const Tbp2 = antoineBoilingPointSolverLocal(data2.antoine, currentFixedPressurePa);
                
                let initialTempGuess: number;
                if (Tbp1 !== null && Tbp2 !== null) {
                    initialTempGuess = x1_val * Tbp1 + (1 - x1_val) * Tbp2;
                } else if (Tbp1 !== null) {
                    initialTempGuess = Tbp1;
                } else if (Tbp2 !== null) {
                    initialTempGuess = Tbp2;
                } else {
                    initialTempGuess = 373.15;
                }
                initialTempGuess = Math.max(150, Math.min(initialTempGuess, 700));

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
            }
        }
        
        calculated_x_values.push(1.0);
        calculated_y_values.push(1.0);

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
        // if (showLoading) setEquilibriumData(null); // This line is removed to keep the graph visible on error.
    } finally {
        if (showLoading) setLoading(false);
    }
  }, [comp1Name, comp2Name, useTemperature, temperatureC, pressureBar, fluidPackage, comp1Data, comp2Data, interactionParams, displayedFluidPackage]);

  // REMOVED: useEffect that automatically called calculateEquilibriumCurve
  // useEffect(() => {
  //   calculateEquilibriumCurve();
  // }, [calculateEquilibriumCurve]); 

  // Initial load with default values can be triggered manually by the user via the button.
  // If an initial automatic load with defaults is strictly required *once* on page load,
  // a separate useEffect with an empty dependency array could call calculateEquilibriumCurve(),
  // but the current request is to stop fetching until the button is pressed.

  const interpolateY = (xTarget: number, xData: number[], yData: number[]): number | null => {
    if (!xData || !yData || xData.length === 0 || xData.length !== yData.length) return null;
    
    // Handle cases where xTarget is outside the data range (clamp to ends)
    if (xTarget <= xData[0]) return yData[0];
    if (xTarget >= xData[xData.length - 1]) return yData[yData.length - 1];
  
    for (let i = 0; i < xData.length - 1; i++) {
      if (xTarget >= xData[i] && xTarget <= xData[i + 1]) {
        const x1 = xData[i], y1 = yData[i];
        const x2 = xData[i + 1], y2 = yData[i + 1];
        if (x2 === x1) return y1; // Avoid division by zero if points are identical
        return y1 + (y2 - y1) * (xTarget - x1) / (x2 - x1);
      }
    }
    return null; // Should ideally be covered by boundary checks
  };

  const interpolateX = (yTarget: number, yData: number[], xData: number[]): number | null => {
    if (!xData || !yData || yData.length === 0 || yData.length !== xData.length) return null;

    // Attempt to find the correct segment for interpolation.
    // This assumes yData is reasonably monotonic for typical VLE.
    let bestSegment = -1;

    // First pass: direct match or exact boundary
    if (yTarget <= yData[0]) { // Using the sorted nature of xData, yData should correspond
        if (yData[0] < yData[yData.length-1]) return xData[0]; // y is increasing
        else return xData[xData.length-1]; // y is decreasing (less common for typical x-y)
    }
    if (yTarget >= yData[yData.length - 1]) {
        if (yData[0] < yData[yData.length-1]) return xData[xData.length-1];
        else return xData[0];
    }


    for (let i = 0; i < yData.length - 1; i++) {
        const y1 = yData[i];
        const y2 = yData[i+1];
        if ((y1 <= yTarget && yTarget <= y2) || (y2 <= yTarget && yTarget <= y1)) {
            bestSegment = i;
            break;
        }
    }

    if (bestSegment !== -1) {
        const y1 = yData[bestSegment], x1 = xData[bestSegment];
        const y2 = yData[bestSegment+1], x2 = xData[bestSegment+1];
        if (y2 === y1) return (yTarget === y1) ? x1 : null; // Horizontal segment
        return x1 + (x2 - x1) * (yTarget - y1) / (y2 - y1);
    }
    
    // Fallback if no segment found (e.g. yTarget outside range not caught by initial checks)
    // This part might need more robust handling depending on VLE data characteristics
    console.warn("interpolateX: yTarget might be out of yData range or yData not suitable for simple interpolation.");
    return null;
  };

  const formatNumberToPrecision = (num: any, precision: number = 3): string => {
    if (typeof num === 'number') {
      if (num === 0) return '0';
      const fixed = num.toPrecision(precision);
      // Remove trailing zeros after decimal point, but keep integer part if it ends in zero
      if (fixed.includes('.')) {
        return parseFloat(fixed).toString(); 
      }
      return fixed;
    }
    return String(num); // Return as string if not a number
  };


    const generateEChartsOptions = useCallback((xValues: number[], yValues: number[]) => {
    if (!xValues || !yValues || xValues.length === 0) {
        setEchartsOptions({});
        setRMin(null);
        return;
    }
    const series: SeriesOption[] = [];
    
    // --- Define Operating Lines ---
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

    // --- Define Data for Calculations ---
    let yDataForStepping = yValues;
    if (murphreeEfficiency < 1.0) {
        yDataForStepping = xValues.map((x, i) => {
            const y_eq = yValues[i];
            let y_op: number;
            if (x >= xIntersect) {
                y_op = rectifyingSlope * x + rectifyingIntercept;
            } else {
                y_op = (strippingSlope === Infinity) ? xb : strippingSlope * x + strippingIntercept;
            }
            return y_op + murphreeEfficiency * (y_eq - y_op);
        });
    }

    // --- PASS 1: Calculate Stage Data and Find Final X Position ---
    let stageCount = 0;
    let feedStageCount = 0;
    let currentX = xd;
    let currentY = xd;
    let finalStageX = xb;
    const stageLineData: EChartsPoint[] = [];
    let previousSectionIsRectifying = true;

    for(let i = 0; i < 50; i++) {
        if (currentX <= xb + buffer) break;

        let intersectX = NaN;
        for (let j = 0; j < xValues.length - 1; j++) {
            if (yDataForStepping[j] === null || yDataForStepping[j+1] === null) continue;
            if ((yDataForStepping[j] <= currentY && yDataForStepping[j+1] >= currentY) || (yDataForStepping[j] >= currentY && yDataForStepping[j+1] <= currentY)) {
                 if (Math.abs(yDataForStepping[j+1] - yDataForStepping[j]) < 1e-9) {
                     intersectX = (yDataForStepping[j] === currentY) ? xValues[j] : xValues[j+1];
                 } else {
                     const fraction = (currentY - yDataForStepping[j]) / (yDataForStepping[j+1] - yDataForStepping[j]);
                     intersectX = xValues[j] + fraction * (xValues[j+1] - xValues[j]);
                 }
                 break;
            }
        }
        if (isNaN(intersectX)) break;

        finalStageX = intersectX;
        stageCount++;

        stageLineData.push([currentX, currentY]);
        stageLineData.push([intersectX, currentY]);
        stageLineData.push([null, null]);

        const isFinalStage = intersectX <= xb + buffer;
        if (isFinalStage) {
            stageLineData.push([intersectX, currentY]);
            stageLineData.push([intersectX, intersectX]);
            stageLineData.push([null, null]);
            break;
        } else {
            const currentSectionIsRectifying = intersectX >= xIntersect;
            let nextY;
            if (currentSectionIsRectifying) {
                nextY = rectifyingSlope * intersectX + rectifyingIntercept;
            } else {
                nextY = (strippingSlope === Infinity) ? xb : strippingSlope * intersectX + strippingIntercept;
            }
            nextY = Math.max(0, Math.min(1, nextY));

            stageLineData.push([intersectX, currentY]);
            stageLineData.push([intersectX, nextY]);
            stageLineData.push([null, null]);
            
            if (feedStageCount === 0 && currentSectionIsRectifying !== previousSectionIsRectifying) {
                feedStageCount = stageCount;
            }
            previousSectionIsRectifying = currentSectionIsRectifying;
            currentX = intersectX;
            currentY = nextY;
        }
    }
    
    setStages(stageCount);
    if (feedStageCount === 0 && stageCount > 0) setFeedStage(stageCount);
    else setFeedStage(feedStageCount);
    
    // --- PASS 2: Assemble All Chart Series With Finalized Data ---
    series.push({ name: 'Equilibrium Line', type: 'line', data: xValues.map((x, i) => [x, yValues[i]]), color: 'cyan', symbol: 'none', lineStyle: { width: 3.5 }, z: 5, animation: false });
    series.push({ name: 'y = x Line', type: 'line', data: [[0, 0], [1, 1]], color: 'blue', symbol: 'none', lineStyle: { width: 3.5, type: 'dotted' }, animation: false });

    if (murphreeEfficiency < 1.0) {
        const plotData = xValues
            .map((x, i) => (x >= finalStageX && x <= xd) ? [x, yDataForStepping[i]] as EChartsPoint : null)
            .filter(p => p !== null);

        const yAtFinalStageX = interpolateY(finalStageX, xValues, yDataForStepping);
        if (yAtFinalStageX !== null) {
            plotData.unshift([finalStageX, yAtFinalStageX]);
        }
        
        series.push({ name: 'Effective Equilibrium', type: 'line', data: plotData, connectNulls: false, color: 'cyan', symbol: 'none', lineStyle: { width: 3.5, type: 'dashed' }, z: 4, animation: false });
    }

    const feedLineData: EChartsPoint[] = [[xf, xf]];
    if (feedSlope === Infinity) { feedLineData.push([xf, yIntersect]); } else { feedLineData.push([xIntersect, yIntersect]); }
    series.push({ name: 'Rectifying Section', type: 'line', data: [[xd, xd], [xIntersect, yIntersect]], color: 'orange', symbol: 'none', lineStyle: { width: 3.5 }, animation: false });
    series.push({ name: 'Feed Section', type: 'line', data: feedLineData, color: 'red', symbol: 'none', lineStyle: { width: 3.5 }, animation: false });
    series.push({ name: 'Stripping Section', type: 'line', data: [[xIntersect, yIntersect], [xb, xb]], color: 'green', symbol: 'none', lineStyle: { width: 3.5 }, animation: false });
    series.push({ name: 'Key Points', type: 'scatter', data: [{ value: [xd, xd], itemStyle: { color: 'orange' }, name: 'Distillate' },{ value: [xb, xb], itemStyle: { color: 'green' }, name: 'Bottoms' },{ value: [xf, xf], itemStyle: { color: 'red' }, name: 'Feed' }], symbolSize: 8, label: { show: false }, tooltip: { formatter: (params: any) => `${params.name}: (${params.value[0].toFixed(3)}, ${params.value[1].toFixed(3)})` }, animation: false });
    series.push({ name: 'Stages', type: 'line', data: stageLineData, color: 'black', symbol: 'none', lineStyle: { width: 3 }, connectNulls: false, legendHoverLink: false, animation: false });

    // --- FINAL Rmin Calculation ---
    let rMinValue: number | null = null;
    let maxOperatingSlope = -Infinity;

    // Step 1: Calculate q-line/equilibrium intersection
    let intersectionPinch: { x: number; y: number } | null = null;
    if (Math.abs(q - 1) < 1e-9) {
        const y_val = interpolateY(xf, xValues, yValues);
        if (y_val !== null) {
            intersectionPinch = { x: xf, y: y_val };
        }
    } else {
        for (let i = 0; i < xValues.length - 1; i++) {
            const x1 = xValues[i], y1 = yValues[i];
            const x2 = xValues[i + 1], y2 = yValues[i + 1];
            if (x1 === x2) continue;
            const m_eq = (y2 - y1) / (x2 - x1);
            const c_eq = y1 - m_eq * x1;
            const denominator = (q - 1) * m_eq - q;
            if (Math.abs(denominator) < 1e-9) continue;
            const x_int = (-(q - 1) * c_eq - xf) / denominator;
            if (x_int >= Math.min(x1, x2) - 1e-6 && x_int <= Math.max(x1, x2) + 1e-6) {
                const y_int = m_eq * x_int + c_eq;
                intersectionPinch = { x: x_int, y: y_int };
                break;
            }
        }
    }

    if (intersectionPinch && intersectionPinch.x <= xd) {
        if (Math.abs(xd - intersectionPinch.x) > 1e-9) {
            maxOperatingSlope = (xd - intersectionPinch.y) / (xd - intersectionPinch.x);
        }
    }

    // Step 2: Check for a more limiting tangency
    if (q >= 1 || q <= 0) {
        // ✅ FIXED: Start loop at i = 1 to skip the invalid (0, 0) point.
        for (let i = 1; i < xValues.length; i++) {
            const x_eq = xValues[i];
            const y_eq = yValues[i];

            if (x_eq > xd + 1e-6) continue;

            let isCandidate = false;
            if (q >= 1) {
                if (x_eq <= xf + 1e-6) isCandidate = true;
            } else { // q <= 0
                if (x_eq >= xf - 1e-6) isCandidate = true;
            }
            
            if (isCandidate) {
                if (Math.abs(xd - x_eq) < 1e-9) continue;
                const slope = (xd - y_eq) / (xd - x_eq);
                if (slope > maxOperatingSlope) {
                    maxOperatingSlope = slope;
                }
            }
        }
    }

    // Step 3: Final Calculation
    if (maxOperatingSlope > -Infinity && maxOperatingSlope < 1) {
        rMinValue = maxOperatingSlope / (1 - maxOperatingSlope);
    }
    setRMin(rMinValue !== null && rMinValue >= 0 ? rMinValue : null);

    const capitalizeFirst = (str: string) => str ? str.charAt(0).toUpperCase() + str.slice(1) : '';
    const dispComp1Cap = capitalizeFirst(displayedComp1);
    const dispComp2Cap = capitalizeFirst(displayedComp2);
    const titleCondition = displayedUseTemp ? 
        (displayedTemp !== null ? `${displayedTemp.toFixed(1)} °C` : 'N/A Temp') : 
        (displayedPressure !== null ? `${displayedPressure.toFixed(2)} bar` : 'N/A Pressure');
    const titleText = displayedComp1 && displayedComp2 ?
        `McCabe-Thiele: ${dispComp1Cap}/${dispComp2Cap} at ${titleCondition}` :
        'McCabe-Thiele Diagram';
    
    const axisCompLabel = dispComp1Cap;
    const isDark = resolvedTheme === 'dark';
    const textColor = isDark ? 'white' : '#000000';

    setEchartsOptions({
        backgroundColor: 'transparent',
        title: { text: titleText, left: 'center', textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } },
        grid: { left: '5%', right: '5%', bottom: '5%', top: '5%', containLabel: true },
        xAxis: {
            type: 'value', min: 0, max: 1, interval: 0.1,
            name: `Liquid Mole Fraction ${axisCompLabel} (x)`,
            nameLocation: 'middle', nameGap: 30, nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
            axisLine: { lineStyle: { color: textColor } }, axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
            axisLabel: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans', formatter: '{value}' }, splitLine: { show: false },
        },
        yAxis: {
            type: 'value', min: 0, max: 1, interval: 0.1,
            name: `Vapor Mole Fraction ${axisCompLabel} (y)`,
            nameLocation: 'middle', nameGap: 40, nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
            axisLine: { lineStyle: { color: textColor } }, axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
            axisLabel: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans', formatter: '{value}' }, splitLine: { show: false },
        },
        legend: {
            orient: 'vertical', right: '2%', top: 'center',
            data: (() => {
                const names = series.map(s => s.name).filter((n): n is string => typeof n === 'string');
                const filtered = names.filter(name => name !== 'Key Points' && name !== 'Stages');
                let ordered: string[] = [];
                if (filtered.includes('Equilibrium Line')) ordered.push('Equilibrium Line');
                if (filtered.includes('Effective Equilibrium') && murphreeEfficiency < 1.0) ordered.push('Effective Equilibrium');
                filtered.forEach(name => {
                  if (name !== 'Equilibrium Line' && name !== 'Effective Equilibrium') ordered.push(name);
                });
                return ordered as string[];
            })(),
            textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
            itemWidth: 25, itemHeight: 2,
        },        tooltip: {
            show: true,
            trigger: 'axis',
            triggerOn: 'mousemove',
            backgroundColor: isDark ? '#08306b' : '#ffffff',
            borderColor: isDark ? '#55aaff' : '#333333',
            borderWidth: 1,
            textStyle: {
                color: textColor,
                fontSize: 12,
                fontFamily: 'Merriweather Sans'
            },
            axisPointer: {
                type: 'cross',
                snap: true,
                label: {
                    show: true,
                    backgroundColor: isDark ? '#08306b' : '#ffffff',
                    color: isDark ? 'white' : 'black',
                    borderColor: isDark ? '#55aaff' : '#333333',
                    borderWidth: 1,
                    shadowBlur: 0,
                    shadowColor: 'transparent',
                    fontFamily: 'Merriweather Sans',
                    precision: 2,
                    formatter: function (params: any) {
                        if (params.axisDimension === 'x') {
                            return `x: ${formatNumberToPrecision(params.value, 3)}`;
                        }
                        if (params.axisDimension === 'y') {
                            return `y: ${formatNumberToPrecision(params.value, 3)}`;
                        }
                        return formatNumberToPrecision(params.value, 3);
                    }
                },
                crossStyle: {
                    color: isDark ? '#999' : '#666'
                }
            },
            formatter: function (params: any) {
                let tooltipHtml = '';
                if (Array.isArray(params) && params.length > 0) {
                    const xAxisValue = params[0].axisValue;
                    if (xAxisValue !== undefined) {
                         tooltipHtml += `x: ${formatNumberToPrecision(xAxisValue, 3)}<br/>`;
                    }

                    if (xAxisValue !== undefined) {
                        const operatingLineInfo = [];
                        
                        if (xAxisValue >= xIntersect && xAxisValue <= xd) {
                            operatingLineInfo.push('<span style="color: orange;">Rectifying Section</span>');
                        }
                        
                        if (xAxisValue >= xb && xAxisValue <= xIntersect) {
                            operatingLineInfo.push('<span style="color: green;">Stripping Section</span>');
                        }
                        
                        if (operatingLineInfo.length > 0) {
                            tooltipHtml += operatingLineInfo.join(' & ') + '<br/>';
                        }
                        
                        if (stages !== null && stages > 0) {
                            const isOnVerticalStageLine = params.some((param: any) => {
                                return param.seriesName === 'Stages' && param.value && 
                                       Array.isArray(param.value) && param.value.length >= 2 &&
                                       Math.abs(param.value[0] - xAxisValue) < 0.001;
                            });
                            
                            if (!isOnVerticalStageLine) {
                                const stageWidth = (xd - xb) / stages;
                                const stageNumber = Math.ceil((xd - xAxisValue) / stageWidth);
                                if (stageNumber >= 1 && stageNumber <= stages) {
                                    tooltipHtml += `<span style="color: purple;">Stage ${stageNumber}</span><br/>`;
                                }
                            }
                        }
                    }

                    params.forEach((param: any, index: number) => {
                        const seriesName = param.seriesName;
                        if (seriesName === 'Key Points') {
                            if (param.data && param.data.name) {
                                tooltipHtml += `<span style="color: ${param.color};">${param.data.name}: (${formatNumberToPrecision(param.value[0], 3)}, ${formatNumberToPrecision(param.value[1], 3)})</span><br/>`;
                            }
                        } else if (seriesName === 'Stages') {
                            const isVerticalLine = param.value && Array.isArray(param.value) && param.value.length >= 2 &&
                                                  Math.abs(param.value[0] - xAxisValue) < 0.001;
                            if (!isVerticalLine) {
                                tooltipHtml += `<span style="color: ${param.color};">${seriesName}: ${formatNumberToPrecision(param.value[1], 3)}</span><br/>`;
                            }
                        } else if (seriesName && param.value && Array.isArray(param.value) && param.value.length >= 2) {
                            tooltipHtml += `<span style="color: ${param.color};">${seriesName}: ${formatNumberToPrecision(param.value[1], 3)}</span><br/>`;
                        }
                    });
                }
                return tooltipHtml;
            }
        }, 
        animationDuration: 300, animationEasing: 'cubicInOut',
        toolbox: {
            show: false
        },
        series: series,
    });
  }, [xd, xb, xf, q, r, murphreeEfficiency, displayedComp1, displayedComp2, displayedTemp, displayedPressure, displayedUseTemp, displayedFluidPackage, buffer, resolvedTheme]);

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

  const handleKeyDown = (event: React.KeyboardEvent<HTMLInputElement>, inputTarget?: 'comp1' | 'comp2') => {
    if (event.key === 'Enter') {
      event.preventDefault();
      if (inputTarget === 'comp1' && showComp1Suggestions && comp1Suggestions.length > 0) {
        // Potentially select first suggestion or handle differently
        // For now, just close suggestions and update graph
        setShowComp1Suggestions(false);
      } else if (inputTarget === 'comp2' && showComp2Suggestions && comp2Suggestions.length > 0) {
        setShowComp2Suggestions(false);
      }
      // Trigger graph update on Enter key, which now internally syncs numeric states
      handleUpdateGraphClick();
    } else if (event.key === 'Escape') {
        if (inputTarget === 'comp1') setShowComp1Suggestions(false);
        if (inputTarget === 'comp2') setShowComp2Suggestions(false);
    }
  };
  
  const handleSwapComponents = () => {
    const tempName = comp1Name;
    setComp1Name(comp2Name);
    setComp2Name(tempName);
    setAutoGeneratePending(true); // trigger auto regenerate after swap
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

  // Helper to compute a 'nice' slider step given the maximum value
  const computeStep = (maxVal: number): number => {
    const desiredSteps = 100;
    const raw = maxVal / desiredSteps;
    const exponent = Math.floor(Math.log10(raw));
    const pow10 = Math.pow(10, exponent);
    const mant = raw / pow10;
    let niceMant: number;
    if (mant <= 1) niceMant = 1;
    else if (mant <= 2) niceMant = 2;
    else if (mant <= 5) niceMant = 5;
    else niceMant = 10;
    return niceMant * pow10;
  };

  // Silent recalculation when fixed condition slider changes
  const hasMounted = useRef(false);
  useEffect(() => {
      if (!hasMounted.current) {
          hasMounted.current = true;
          return;
      }
      // Use a lighter resolution (21 points) for real-time slider updates
      calculateEquilibriumCurve(false, 21);
  }, [temperatureC, pressureBar]);

  // When autoGeneratePending is set, wait for the relevant state to update (next render) then regenerate the graph
  useEffect(() => {
    if (autoGeneratePending) {
      handleUpdateGraphClick();
      setAutoGeneratePending(false);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [comp1Name, comp2Name, useTemperature, fluidPackage, autoGeneratePending]);

  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8 px-4 md:px-16">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          <div className="lg:col-span-1 space-y-6">
            <Card>
              <CardContent className="space-y-6 p-4">
                <div className="space-y-4">
                  {/* Temperature/Pressure Toggle and Input */}
                  <div>
                    <Tabs value={useTemperature ? "temperature" : "pressure"} onValueChange={(value) => { 
                        setUseTemperature(value === "temperature");
                        setAutoGeneratePending(true); // trigger auto regenerate after mode change
                    }}>
                      <TabsList className="flex w-full">
                        <TabsTrigger value="temperature" className="flex-1 text-center">Temperature</TabsTrigger>
                        <TabsTrigger value="pressure" className="flex-1 text-center">Pressure</TabsTrigger>
                      </TabsList>
                    </Tabs>
                  </div>

                  {/* Component Inputs - Side by side with placeholders and swap button */}
                  <div className="flex items-center gap-2">
                    <div className="relative flex-1">
                      <Input
                        ref={input1Ref}
                        id="comp1Name"
                        value={comp1Name}
                        onChange={handleComp1NameChange}
                        onKeyDown={(e) => handleKeyDown(e, 'comp1')}
                        onFocus={() => {
                            setActiveSuggestionInput('comp1');
                            if (comp1Name.trim() !== "" && comp1Suggestions.length > 0) setShowComp1Suggestions(true);
                            else if (comp1Name.trim() !== "") fetchSuggestions(comp1Name, 'comp1');
                        }}
                        placeholder="Methanol"
                        required
                        className="w-full"
                        autoComplete="off"
                      />
                      {showComp1Suggestions && comp1Suggestions.length > 0 && (
                        <div ref={activeSuggestionInput === 'comp1' ? activeComponentRef : null} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1">
                          {comp1Suggestions.map((suggestion, index) => (
                            <div
                              key={index}
                              onClick={() => handleSuggestionClick(suggestion, 'comp1')}
                              className="px-3 py-2 hover:bg-accent cursor-pointer text-sm"
                            >
                              {suggestion}
                            </div>
                          ))}
                        </div>
                      )}
                    </div>
                    <Button variant="ghost" size="icon" onClick={handleSwapComponents} title="Swap Components">
                        <ArrowLeftRight className="h-4 w-4" />
                    </Button>
                    <div className="relative flex-1">
                      <Input
                        ref={input2Ref}
                        id="comp2Name"
                        value={comp2Name}
                        onChange={handleComp2NameChange}
                        onKeyDown={(e) => handleKeyDown(e, 'comp2')}
                        onFocus={() => {
                            setActiveSuggestionInput('comp2');
                            if (comp2Name.trim() !== "" && comp2Suggestions.length > 0) setShowComp2Suggestions(true);
                            else if (comp2Name.trim() !== "") fetchSuggestions(comp2Name, 'comp2');
                        }}
                        placeholder="Water"
                        required
                        className="w-full"
                        autoComplete="off"
                      />
                      {showComp2Suggestions && comp2Suggestions.length > 0 && (
                        <div ref={activeSuggestionInput === 'comp2' ? activeComponentRef : null} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1">
                          {comp2Suggestions.map((suggestion, index) => (
                            <div
                              key={index}
                              onClick={() => handleSuggestionClick(suggestion, 'comp2')}
                              className="px-3 py-2 hover:bg-accent cursor-pointer text-sm"
                            >
                              {suggestion}
                            </div>
                          ))}
                        </div>
                      )}
                    </div>
                  </div>
                  {/* Fluid Package Selector - Label side-by-side, reordered */}
                  <div className="flex items-center gap-2">
                      <Label htmlFor="fluidPackageMcCabe" className="text-sm font-medium whitespace-nowrap">Fluid Package:</Label>
                      <Select value={fluidPackage} onValueChange={(value) => {
                          setFluidPackage(value as FluidPackageTypeMcCabe);
                          setAutoGeneratePending(true); // trigger auto regenerate after fluid package change
                      }}>
                          <SelectTrigger id="fluidPackageMcCabe" className="flex-1"><SelectValue placeholder="Select fluid package" /></SelectTrigger>
                          <SelectContent>
                              <SelectItem value="wilson">Wilson</SelectItem>
                              <SelectItem value="uniquac">UNIQUAC</SelectItem>
                              <SelectItem value="nrtl">NRTL</SelectItem>
                              <SelectItem value="unifac">UNIFAC</SelectItem>
                              <SelectItem value="pr">Peng-Robinson</SelectItem>
                              <SelectItem value="srk">SRK</SelectItem>
                          </SelectContent>
                      </Select>
                  </div>
                </div>
                <Button onClick={handleUpdateGraphClick} disabled={loading} className="w-full">
                  {loading ? 'Calculating...' : 'Update Graph'}
                </Button>
                {error && <p className="text-sm text-red-500 mt-2">{error}</p>}
              </CardContent>
            </Card>
            {/* Slider Card */}
            <Card>
               <CardContent className="space-y-8 pt-2 pb-2">
                {/* Fixed condition slider */}
                <div className="space-y-2">
                  <div className="flex items-center justify-between">
                    <Label htmlFor="fixedCondition">
                      {useTemperature ? `Temperature (°C): ${temperatureC}` : `Pressure (bar): ${pressureBar}`}
                    </Label>
                    <div className="flex items-center gap-1">
                      <Label className="text-xs text-muted-foreground">Max:</Label>
                      <Input
                        type="text"
                        value={useTemperature ? tempMax : pressureMax}
                        onChange={(e: React.ChangeEvent<HTMLInputElement>) => {
                          const raw = e.target.value;
                          const num = parseFloat(raw);
                          const CLAMP_MAX = 2000;
                          if (useTemperature) {
                            if (raw === "" || isNaN(num)) {
                              setTempMax(raw);
                            } else {
                              setTempMax(String(Math.min(num, CLAMP_MAX)));
                            }
                          } else {
                            if (raw === "" || isNaN(num)) {
                              setPressureMax(raw);
                            } else {
                              setPressureMax(String(Math.min(num, CLAMP_MAX)));
                            }
                          }
                        }}
                        className="w-16 h-8 text-xs"
                      />
                    </div>
                  </div>
                  <Slider
                    id="fixedCondition"
                    min={0}
                    max={useTemperature ? parseFloat(tempMax) : parseFloat(pressureMax)}
                    step={computeStep(useTemperature ? parseFloat(tempMax) : parseFloat(pressureMax))}
                    value={[useTemperature ? (temperatureC || 0) : (pressureBar || 0)]}
                    onValueChange={([v]: number[]) => {
                      if (useTemperature) setTemperatureC(v);
                      else setPressureBar(v);
                    }}
                    className="w-full"
                  />
                </div>
                 {/* xd Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="xd">
                     <span dangerouslySetInnerHTML={{ __html: `Distillate Composition (x<sub>D</sub>): ${xd.toFixed(2)}` }} />
                   </Label>
                   <Slider id="xd" min={0.01} max={0.99} step={0.01} value={[xd]} onValueChange={(value) => updateCompositions('xd', value[0])} style={{ '--primary': 'hsl(38 92% 50%)' } as React.CSSProperties}/>
                 </div>
                 {/* xf Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="xf">
                     <span dangerouslySetInnerHTML={{ __html: `Feed Composition (x<sub>F</sub>): ${xf.toFixed(2)}` }} />
                   </Label>
                   <Slider id="xf" min={0.01} max={0.99} step={0.01} value={[xf]} onValueChange={(value) => updateCompositions('xf', value[0])} style={{ '--primary': 'hsl(0 84% 60%)' } as React.CSSProperties}/>
                 </div>
                 {/* xb Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="xb">
                     <span dangerouslySetInnerHTML={{ __html: `Bottoms Composition (x<sub>B</sub>): ${xb.toFixed(2)}` }} />
                   </Label>
                   <Slider id="xb" min={0.01} max={0.99} step={0.01} value={[xb]} onValueChange={(value) => updateCompositions('xb', value[0])} style={{ '--primary': 'hsl(142 71% 45%)' } as React.CSSProperties}/>
                 </div>
                 {/* q Slider */}
                 <div className="space-y-3">
                   <div className="flex items-center justify-between">
                     <Label htmlFor="q" className="flex items-center">Feed Quality (q): {q.toFixed(2)}
                         <ShadTooltip><TooltipTrigger asChild><Button variant="ghost" size="icon" className="h-5 w-5 rounded-full"><span className="text-xs">ⓘ</span></Button></TooltipTrigger><TooltipContent><p>{getFeedQualityState()}</p></TooltipContent></ShadTooltip>
                     </Label>
                     <div className="flex items-center gap-1">
                       <Label className="text-xs text-muted-foreground">Min:</Label>
                       <Input
                         type="text"
                         value={qMinSlider}
                         onChange={(e: React.ChangeEvent<HTMLInputElement>) => setQMinSlider(e.target.value)}
                         className="w-12 h-8 text-xs"
                       />
                       <Label className="text-xs text-muted-foreground">Max:</Label>
                       <Input
                         type="text"
                         value={qMaxSlider}
                         onChange={(e: React.ChangeEvent<HTMLInputElement>) => setQMaxSlider(e.target.value)}
                         className="w-12 h-8 text-xs"
                       />
                     </div>
                   </div>
                   <Slider id="q" min={parseFloat(qMinSlider)} max={parseFloat(qMaxSlider)} step={0.05} value={[q]} onValueChange={(value) => setQ(value[0])} style={{ '--primary': 'hsl(60 100% 50%)' } as React.CSSProperties}/>
                 </div>
                 {/* r Slider */}
                 <div className="space-y-3">
                   <div className="flex justify-between items-center">
                     <Label htmlFor="r" className="flex justify-between items-center">
                       <span>Reflux Ratio (R): {r.toFixed(2)}</span>
                       {rMin !== null && (
                          <span className="text-xs text-muted-foreground ml-2">
                              R<sub>min</sub>: {rMin.toFixed(2)}
                          </span>
                       )}
                     </Label>
                     <div className="flex items-center gap-1">
                       <Label className="text-xs text-muted-foreground">Min:</Label>
                       <Input
                         type="text"
                         value={rMinSlider}
                         onChange={(e: React.ChangeEvent<HTMLInputElement>) => setRMinSlider(e.target.value)}
                         className="w-12 h-8 text-xs"
                       />
                       <Label className="text-xs text-muted-foreground">Max:</Label>
                       <Input
                         type="text"
                         value={rMaxSlider}
                         onChange={(e: React.ChangeEvent<HTMLInputElement>) => setRMaxSlider(e.target.value)}
                         className="w-12 h-8 text-xs"
                       />
                     </div>
                   </div>
                   <Slider id="r" min={parseFloat(rMinSlider)} max={parseFloat(rMaxSlider)} step={0.05} value={[r]} onValueChange={(value) => setR(value[0])} style={{ '--primary': 'hsl(262 84% 58%)' } as React.CSSProperties}/>
                 </div>
                 {/* Murphree Efficiency Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="murphree" className="flex justify-between items-center">
                     <span dangerouslySetInnerHTML={{ __html: `Murphree Efficiency (E<sub>M</sub>): ${(murphreeEfficiency * 100).toFixed(0)}%` }} />
                   </Label>
                   <Slider
                     id="murphree"
                     min={0.01}
                     max={1}
                     step={0.01}
                     value={[murphreeEfficiency]}
                     onValueChange={(value) => setMurphreeEfficiency(value[0])}
                     style={{ '--primary': 'hsl(180 100% 50%)' } as React.CSSProperties}
                   />
                 </div>
               </CardContent>
            </Card>
            {stages !== null && feedStage !== null && (
              <Card>
                <CardContent>
                  <div className="grid grid-cols-2 gap-4 text-center">
                    <div className="p-4 bg-muted rounded-md">
                      <p className="text-sm font-medium">Number of Stages</p>
                      <p className="text-2xl font-bold">{stages}</p>
                    </div>
                    <div className="p-4 bg-muted rounded-md">
                      <p className="text-sm font-medium">Feed Stage</p>
                      <p className="text-2xl font-bold">{feedStage}</p>
                    </div>
                  </div>
                </CardContent>
              </Card>
            )}
          </div>
          {/* Column 2: Plot and Results Cards */}
          <div className="lg:col-span-2 space-y-6">
            <Card>
              <CardContent className="py-2">
                <div className="relative aspect-square rounded-md">
                   {/* The loading overlay has been removed. The button shows "Calculating..." instead. */}

                   {!equilibriumData && !error && ( <div className="absolute inset-0 flex items-center justify-center text-muted-foreground">Please provide inputs and update graph.</div> )}
                   
                   {/* The error message will now overlay the graph if a calculation fails. */}
                   {error && !loading && ( <div className="absolute inset-0 flex items-center justify-center bg-background/80 text-red-400 p-4 text-center rounded-md">{error}</div> )}

                  {/* The graph is shown as long as equilibriumData exists, even during reloads. */}
                  {equilibriumData && Object.keys(echartsOptions).length > 0 && (
                    <ReactECharts 
                  key={`echarts-${murphreeEfficiency === 1.0 ? 'no-eff' : 'with-eff'}`}
                  ref={echartsRef} 
                  echarts={echarts} 
                  option={echartsOptions} 
                  style={{ height: '100%', width: '100%', borderRadius: '0.375rem', overflow: 'hidden' }} 
                  notMerge={true} 
                  lazyUpdate={false} 
                />
                  )}
                </div>
              </CardContent>
            </Card>
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}