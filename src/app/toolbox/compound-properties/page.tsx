'use client';

import React, { useState, useEffect, useCallback, useRef } from 'react';
import { createClient, SupabaseClient } from '@supabase/supabase-js';

import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsOption, SeriesOption as EChartsSeriesOption } from 'echarts';
import { LineChart } from 'echarts/charts';
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  ToolboxComponent,
  DataZoomComponent, // Keep DataZoomComponent if you plan to re-add zoom later
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

echarts.use([
  TitleComponent, TooltipComponent, GridComponent, LegendComponent, ToolboxComponent, DataZoomComponent, LineChart, CanvasRenderer
]);

import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Terminal, Trash2, PlusCircle } from 'lucide-react'; 
import { TooltipProvider } from "@/components/ui/tooltip";
import { Label } from "@/components/ui/label"; 

import {
    Select,
    SelectContent,
    SelectItem,
    SelectTrigger,
    SelectValue,
} from "@/components/ui/select";

import {
    calculateEq101,
    calculateEq105,
    calculatePolynomial,
    calculateEq106,
    calculateEq102_conductivity_viscosity,
    calculateEq104_virial,
    calculateEq16Complex, // Added import
    parseCoefficient
} from '@/lib/property-equations';

// Supabase Client Setup
const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL || '';
const supabaseAnonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY || '';
let supabase: SupabaseClient;
if (supabaseUrl && supabaseAnonKey) {
    try {
        supabase = createClient(supabaseUrl, supabaseAnonKey);
    } catch (error) {
        console.error("Error initializing Supabase client for Compound Properties:", error);
    }
} else {
    console.error("Supabase URL or Anon Key is missing for Compound Properties.");
}

// Updated formatNumberToPrecision function (assuming 2-argument version from previous state)
const formatNumberToPrecision = (num: any, precision: number = 4): string => {
    if (typeof num !== 'number' || isNaN(num)) return String(num);
    if (num === 0) return '0';

    const s = num.toPrecision(precision);

    if (s.includes('e')) {
        const numFromScientific = Number(s);
        const plainString = numFromScientific.toString();
        if (!plainString.includes('e') && (plainString !== "0" || numFromScientific === 0)) {
            if (plainString.includes('.')) {
                return parseFloat(plainString).toString();
            }
            return plainString;
        }
        return s;
    } else {
        if (s.includes('.')) {
            return parseFloat(s).toString();
        }
        return s;
    }
};


interface PropertyDefinition {
  displayName: string;
  jsonKey: string;
  equationType: 'eq101' | 'eq105' | 'polynomial' | 'eq106' | 'eq102_cv' | 'eq104_virial' | 'eq16_complex'; // Added eq16_complex
  yAxisIndex: number;
  targetUnitName: string;
  color: string; // This will be used as a base, might need array of colors for multiple compounds
  coeffs: string[];
  requiresTc?: boolean;
  requiresMolarMass?: boolean;
  conversionFactor?: number;
  equationTemplate?: string;
}

const baseColors = ['#5470C6', '#91CC75', '#FAC858', '#EE6666', '#73C0DE', '#3BA272', '#FC8452', '#9A60B4', '#EA7CCC'];
const propertiesToPlotConfig: PropertyDefinition[] = [
  { displayName: "Vapor Pressure", jsonKey: "Vapour pressure", equationType: "eq101", yAxisIndex: 0, targetUnitName: "Pa", color: baseColors[0], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "P = exp(A + B/T + C ln(T) + D T<sup>E</sup>)" },
  { displayName: "Liquid Density", jsonKey: "Liquid density", equationType: "eq105", yAxisIndex: 0, targetUnitName: "kg/m³", color: baseColors[1], coeffs: ['A', 'B', 'C', 'D'], requiresMolarMass: true, requiresTc: true, equationTemplate: "ρ = (A / B<sup>(1+(1-T/T<sub>c</sub>)<sup>D</sup>)</sup>) MW" },
  { displayName: "Liquid Heat Capacity", jsonKey: "Liquid heat capacity", equationType: "eq16_complex", yAxisIndex: 0, targetUnitName: "J/mol·K", color: baseColors[2], coeffs: ['A', 'B', 'C', 'D', 'E'], conversionFactor: 1/1000, equationTemplate: "Cp = A + exp(B/T + C + D T + E T<sup>2</sup>)" },
  { displayName: "Liquid Viscosity", jsonKey: "Liquid viscosity", equationType: "eq101", yAxisIndex: 0, targetUnitName: "Pa·s", color: baseColors[3], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "μ = exp(A + B/T + C ln(T) + D T<sup>E</sup>)" },
  { displayName: "Heat of Vaporization", jsonKey: "Heat of vaporization", equationType: "eq106", yAxisIndex: 0, targetUnitName: "J/mol", color: baseColors[4], coeffs: ['A', 'B', 'C', 'D', 'E'], requiresTc: true, equationTemplate: "ΔHᵥ = A(1-T/T<sub>c</sub>)<sup>(B+C(T/T<sub>c</sub>)+D(T/T<sub>c</sub>)<sup>2</sup>+E(T/T<sub>c</sub>)<sup>3</sup>)</sup>" },
  { 
    displayName: "Ideal Gas Heat Capacity", 
    jsonKey: "Ideal gas heat capacity", 
    equationType: "eq16_complex",
    yAxisIndex: 0, 
    targetUnitName: "J/mol·K", 
    color: baseColors[5], 
    coeffs: ['A', 'B', 'C', 'D', 'E'], 
    conversionFactor: 1/1000, 
    equationTemplate: "Cp⁰ = A + exp(B/T + C + D·T + E·T<sup>2</sup>)" 
  },
  // { displayName: "Ideal Gas Heat Capacity (RPP)", jsonKey: "Ideal gas heat capacity (RPP)", equationType: "polynomial", yAxisIndex: 0, targetUnitName: "J/mol·K", color: baseColors[6], coeffs: ['A', 'B', 'C', 'D', 'E'], conversionFactor: 1/1000, equationTemplate: "Cp⁰(RPP) = A + B T + C T² + D T³ + E T⁴" }, 
  { displayName: "Liquid Thermal Conductivity", jsonKey: "Liquid thermal conductivity", equationType: "eq16_complex", yAxisIndex: 0, targetUnitName: "W/m·K", color: baseColors[7], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "kL = A + exp(B/T + C + D T + E T<sup>2</sup>)" },
  { displayName: "Second Virial Coefficient", jsonKey: "Second virial coefficient", equationType: "eq104_virial", yAxisIndex: 0, targetUnitName: "m³/kmol", color: baseColors[8], coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "Bᵥ = A + B/T + C/T<sup>2</sup> + D/T<sup>8</sup> + E/T<sup>9</sup>" },
  { displayName: "Solid Density", jsonKey: "Solid density", equationType: "polynomial", yAxisIndex: 0, targetUnitName: "kg/m³", color: '#808080', coeffs: ['A', 'B'], requiresMolarMass: true, equationTemplate: "ρS = (A + B T) MW" },
  { displayName: "Solid Heat Capacity", jsonKey: "Solid heat capacity", equationType: "polynomial", yAxisIndex: 0, targetUnitName: "J/mol·K", color: '#FFD700', coeffs: ['A', 'B', 'C', 'D', 'E'], conversionFactor: 1/1000, equationTemplate: "CpS = A + B T + C T<sup>2</sup> + D T<sup>3</sup> + E T<sup>4</sup>" },
  { displayName: "Surface Tension", jsonKey: "Surface tension", equationType: "eq16_complex", yAxisIndex: 0, targetUnitName: "N/m", color: '#00CED1', coeffs: ['A', 'B', 'C', 'D', 'E'], equationTemplate: "σ = A + exp(B/T + C + D T + E T<sup>2</sup>)" },
  { displayName: "Vapour Thermal Conductivity", jsonKey: "Vapour thermal conductivity", equationType: "eq102_cv", yAxisIndex: 0, targetUnitName: "W/m·K", color: '#DA70D6', coeffs: ['A', 'B', 'C', 'D'], equationTemplate: "kV = (A T<sup>B</sup>) / (1 + C/T + D/T<sup>2</sup>)" },
  { displayName: "Vapour Viscosity", jsonKey: "Vapour viscosity", equationType: "eq102_cv", yAxisIndex: 0, targetUnitName: "Pa·s", color: '#6A5ACD', coeffs: ['A', 'B', 'C', 'D'], equationTemplate: "μV = (A T<sup>B</sup>) / (1 + C/T + D/T<sup>2</sup>)" },
];

const compoundColors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']; // Colors for up to 4 compounds

interface FetchedCompoundData {
    properties: any;
    molarWeight: number | null;
    criticalTemp: number | null;
    name: string;
}

interface CompoundInputState {
    id: string;
    name: string;
    suggestions: string[];
    showSuggestions: boolean;
    data: FetchedCompoundData | null;
    error: string | null;
    inputRef: React.RefObject<HTMLInputElement | null>;         // was RefObject<HTMLInputElement>
    suggestionsRef: React.RefObject<HTMLDivElement | null>;  // <-- allow null here
}


// Debounce function
function debounce<F extends (...args: any[]) => any>(func: F, waitFor: number) {
  let timeout: ReturnType<typeof setTimeout> | null = null;
  const debounced = (...args: Parameters<F>) => {
    if (timeout !== null) { clearTimeout(timeout); timeout = null; }
    timeout = setTimeout(() => func(...args), waitFor);
  };
  return debounced as (...args: Parameters<F>) => void;
}

export default function CompoundPropertiesPage() {
  const nextCompoundId = useRef(0); 

  const createNewCompoundState = (name: string = ''): CompoundInputState => {
    const id = `compound-${nextCompoundId.current++}`;
    return {
      id, name, suggestions: [], showSuggestions: false, data: null, error: null,
      inputRef: React.createRef<HTMLInputElement | null>(),
      suggestionsRef: React.createRef<HTMLDivElement | null>()
    };
  };

  const [compounds, setCompounds] = useState<CompoundInputState[]>(() => [
    createNewCompoundState('Ethanol') // Use the new function for initial state
  ]);
  
  const [loading, setLoading] = useState(false);
  const [overallError, setOverallError] = useState<string | null>(null); // For global errors like "select a property"
  const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
  const echartsRef = useRef<ReactECharts | null>(null);

  const [selectedPropertyKey, setSelectedPropertyKey] = useState<string>(propertiesToPlotConfig[0].jsonKey);
  const [availablePropertiesForSelection, setAvailablePropertiesForSelection] = useState<PropertyDefinition[]>(propertiesToPlotConfig);
 


  const fetchCompoundPropertiesLocal = async (compoundName: string): Promise<FetchedCompoundData | null> => {
    if (!supabase) throw new Error("Supabase client not initialized.");
    // setLoading(true); // Loading state handled by caller
    // setError(null); // Error state handled by caller (per compound)
    console.log(`CompoundProperties: Fetching data for ${compoundName}...`);

    try {
      const { data: compoundDbData, error: compoundError } = await supabase
        .from('compounds').select('id, name, cas_number').ilike('name', compoundName).limit(1).single();
      if (compoundError || !compoundDbData) throw new Error(compoundError?.message || `Compound '${compoundName}' not found.`);
      
      const compoundId = compoundDbData.id;
      const foundName = compoundDbData.name;

      const sourcesToTry = ['chemsep1', 'DWSIM', 'chemsep2', 'biod_db'];
      let properties: any = null;
      for (const source of sourcesToTry) {
        const { data: propsData, error: propsError } = await supabase
          .from('compound_properties').select('properties').eq('compound_id', compoundId).eq('source', source).single();
        if (!propsError && propsData) { properties = propsData.properties; break; }
      }
      if (!properties) {
        const { data: anyPropsData, error: anyPropsError } = await supabase
          .from('compound_properties').select('properties').eq('compound_id', compoundId).limit(1).single();
        if (anyPropsError || !anyPropsData) throw new Error(`No properties found for ${foundName}.`);
        properties = anyPropsData.properties;
      }
      if (typeof properties !== 'object' || properties === null) throw new Error(`Invalid properties format for ${foundName}.`);

      const molarWeight = parseCoefficient(properties["Molecular weight"]?.value ?? properties["Molecular weight"]);
      const criticalTemp = parseCoefficient(properties["Critical temperature"]?.value ?? properties["Critical temperature"]);

      return { properties, molarWeight, criticalTemp, name: foundName };
    } catch (err: any) {
      console.error(`Error fetching data for ${compoundName}:`, err.message);
      // setError(`Data fetch failed for ${compoundName}: ${err.message}`); // Caller will set this per compound
      throw err; // Re-throw to be caught by caller
    }
  };

  const processAndPlotProperties = useCallback((
    allCompoundsData: CompoundInputState[],
    propertyKey: string | null
  ) => {
    setOverallError(null);
    if (!propertyKey) {
      setEchartsOptions({});
      setOverallError("Please select a property to plot.");
      return;
    }

    const propDef = propertiesToPlotConfig.find(p => p.jsonKey === propertyKey);
    if (!propDef) {
      console.warn(`Property definition for key ${propertyKey} not found.`);
      setEchartsOptions({});
      setOverallError(`Property definition for ${propertyKey} not found.`);
      return;
    }

    const seriesData: EChartsSeriesOption[] = [];
    let commonTmin: number | null = null;
    let commonTmax: number | null = null;
    let atLeastOneCompoundHasData = false;
    let titleCompoundNames = allCompoundsData.filter(c => c.data).map(c => c.data!.name).join(' vs ');
    if (!titleCompoundNames && allCompoundsData.length > 0 && allCompoundsData[0].name) {
        titleCompoundNames = allCompoundsData[0].name; 
    }


    allCompoundsData.forEach((compoundState, compoundIndex) => {
        if (!compoundState.data || !compoundState.data.properties) {
            // console.warn(`No data for compound ${compoundState.name} or property ${propDef.displayName}`);
            return;
        }
        atLeastOneCompoundHasData = true;
        const currentCompoundData = compoundState.data;
        const propData = currentCompoundData.properties[propDef.jsonKey];

        if (!propData) {
            console.warn(`Property ${propDef.displayName} (${propDef.jsonKey}) not found for ${currentCompoundData.name}.`);
            return;
        }

        const Tmin = parseCoefficient(propData.Tmin?.value ?? propData.Tmin);
        const Tmax = parseCoefficient(propData.Tmax?.value ?? propData.Tmax);

        if (Tmin === null || Tmax === null || Tmin >= Tmax) {
            console.warn(`Invalid Tmin/Tmax for ${propDef.displayName} for ${currentCompoundData.name}. Tmin: ${Tmin}, Tmax: ${Tmax}`);
            return;
        }
        
        // Update common temperature range
        if (commonTmin === null || Tmin > commonTmin) commonTmin = Tmin;
        if (commonTmax === null || Tmax < commonTmax) commonTmax = Tmax;


        const rawCoeffs = propDef.coeffs.map(cKey => ({ key: cKey, value: parseCoefficient(propData[cKey]) }));
        
        // --- Debugging: Log Ideal Gas Heat Capacity Coefficients ---
        if (propDef.jsonKey === "Ideal gas heat capacity") {
            console.log(`DEBUG: Ideal Gas Heat Capacity coeffs for ${currentCompoundData.name}:`, 
                JSON.stringify(rawCoeffs.reduce((acc, coeff) => {
                    acc[coeff.key] = coeff.value;
                    return acc;
                }, {} as Record<string, number | null>), null, 2)
            );
        }
        // --- End Debugging ---

        let Tc_K: number | null = null;
        if (propDef.requiresTc) {
            Tc_K = currentCompoundData.criticalTemp ?? null;
            if (Tc_K === null) {
                console.warn(`Critical Temperature (Tc) required for ${propDef.displayName} for ${currentCompoundData.name} but not found.`);
                return;
            }
        }
        
        let mw_kg_kmol: number | null = null;
        if (propDef.requiresMolarMass) {
            mw_kg_kmol = currentCompoundData.molarWeight ?? null;
            if (mw_kg_kmol === null) {
                console.warn(`Molar Mass required for ${propDef.displayName} for ${currentCompoundData.name} but not found.`);
                return;
            }
        }

        let currentEquation = propDef.equationTemplate || propDef.displayName;
        rawCoeffs.forEach(c => {
            if (c.value !== null) {
                currentEquation = currentEquation.replace(new RegExp(`\\b${c.key}\\b`, 'g'), formatNumberToPrecision(c.value, 2));
            } else {
                currentEquation = currentEquation.replace(new RegExp(`\\b${c.key}\\b`, 'g'), "0");
            }
        });
        if (Tc_K !== null) {
            // Display Tc in Celsius in the equation string
            const tcCelsiusDisplay = formatNumberToPrecision(Tc_K - 273.15, 2);
            currentEquation = currentEquation.replace(/\bTc\b/g, `${tcCelsiusDisplay}°C`); // Or T<sub>c</sub> later
        }
        if (mw_kg_kmol !== null) currentEquation = currentEquation.replace(/\bMW\b/g, formatNumberToPrecision(mw_kg_kmol, 4));
        
        // Apply HTML for subscripts/superscripts and clean up "+ -"
        currentEquation = currentEquation.replace(/T²|T\^2/g, 'T<sup>2</sup>')
                                     .replace(/T³|T\^3/g, 'T<sup>3</sup>')
                                     .replace(/T⁴|T\^4/g, 'T<sup>4</sup>')
                                     .replace(/\+ -/g, '- '); // Replace "+ -" with "- "
        currentEquation = currentEquation.replace(/\bTc\b/g, 'T<sub>c</sub>');


        const legendEquation = currentEquation;

        const points: [number, number][] = [];
        // Use the determined commonTmin and commonTmax for plotting range if available, otherwise individual Tmin/Tmax
        const plotTmin = commonTmin ?? Tmin;
        const plotTmax = commonTmax ?? Tmax;

        if (plotTmin >= plotTmax) { // If common range is invalid (e.g. no overlap)
            console.warn(`No valid common temperature range for plotting ${propDef.displayName} for ${currentCompoundData.name}. Individual range: ${Tmin}-${Tmax}`);
            // Fallback to individual range or skip plotting this compound for this property
            // For now, let's try with individual range if common range is bad.
            // This part might need refinement based on desired behavior for non-overlapping ranges.
            // A better approach might be to only plot if there's *some* overlap.
            // For now, let's assume we want to plot on the *intersection* of valid ranges.
            // The commonTmin/commonTmax logic above tries to find this intersection.
            // If after processing all compounds, commonTmin >= commonTmax, then no valid intersection.
            // This check is now done before the loop for points.
        }


        const numSteps = 50;
        // This loop for points will be done *after* all compounds are processed to determine final commonTmin/commonTmax
        // For now, this structure is per-compound, which is fine, but the T range needs to be consistent.
        // Let's defer point generation until after commonTmin/commonTmax are finalized.
        // This part is just to gather data for each series.
        
        // Store necessary info to generate points later
        seriesData.push({
            name: `${currentCompoundData.name}`, // Legend shows only compound name
            type: 'line',
            data: [], 
            yAxisIndex: 0,
            smooth: true,
            lineStyle: { color: compoundColors[compoundIndex % compoundColors.length], width: 2.5 },
            itemStyle: { color: compoundColors[compoundIndex % compoundColors.length] },
            symbolSize: 6,
            // Store extra info for point generation
            _internal_prop_data: propData,
            _internal_compound_data: currentCompoundData,
            _internal_legend_equation: legendEquation,
        } as any); // Use 'as any' for custom internal properties temporarily
    });

    if (!atLeastOneCompoundHasData) {
        setEchartsOptions({});
        setOverallError("No data available for the selected property and compounds.");
        return;
    }
    
    if (commonTmin === null || commonTmax === null || commonTmin >= commonTmax) {
        setEchartsOptions({});
        setOverallError(`No overlapping temperature range found for property '${propDef.displayName}' across selected compounds.`);
        return;
    }

    // Convert common Kelvin limits to Celsius
    const commonTminCelsius = commonTmin - 273.15;
    const commonTmaxCelsius = commonTmax - 273.15;
    const rangeCelsius = commonTmaxCelsius - commonTminCelsius;

    let celsiusStep;
    if (rangeCelsius <= 10) {
        celsiusStep = 2; // Finer step for small ranges
    } else if (rangeCelsius <= 25) {
        celsiusStep = 5;
    } else if (rangeCelsius <= 50) {
        celsiusStep = 10;
    } else if (rangeCelsius <= 200) {
        celsiusStep = 20;
    } else if (rangeCelsius <= 500) {
        celsiusStep = 50;
    } else {
        celsiusStep = 100;
    }

    let finalPlotTminCelsius = Math.floor(commonTminCelsius / celsiusStep) * celsiusStep;
    let finalPlotTmaxCelsius = Math.ceil(commonTmaxCelsius / celsiusStep) * celsiusStep;

    if (finalPlotTminCelsius === finalPlotTmaxCelsius) {
        finalPlotTmaxCelsius = finalPlotTminCelsius + celsiusStep;
    }
    
    // Ensure at least one step if the range is very small, even after ceiling.
    // This can happen if commonTmin and commonTmax are very close and fall into the same step.
    if (finalPlotTmaxCelsius <= finalPlotTminCelsius) { // Should be strictly greater
        finalPlotTmaxCelsius = finalPlotTminCelsius + celsiusStep;
    }


    // Convert nice Celsius plot limits back to Kelvin for ECharts axis and data generation
    const finalPlotTminKelvin = finalPlotTminCelsius + 273.15;
    const finalPlotTmaxKelvin = finalPlotTmaxCelsius + 273.15;


    // Now generate points for each series using the commonTmin and commonTmax
    seriesData.forEach(series => {
        const s = series as any; // Cast to access internal properties
        const propData = s._internal_prop_data;
        const currentCompoundData = s._internal_compound_data as FetchedCompoundData;
        const points: [number, number][] = [];
        const numSteps = 100;
        const tempStep = (finalPlotTmaxKelvin! - finalPlotTminKelvin!) / numSteps;

        const rawCoeffs = propDef.coeffs.map(cKey => ({ key: cKey, value: parseCoefficient(propData[cKey]) }));
        const passedCoeffs = rawCoeffs.map(c => c.value);

        // --- START DEBUGGING BLOCK ---
        if (propDef.jsonKey === "Ideal gas heat capacity" || propDef.jsonKey === "Second virial coefficient") {
            console.log(`DEBUG: Processing "${propDef.displayName}" for "${currentCompoundData.name}"`);
            console.log(`DEBUG:   Equation Type: ${propDef.equationType}`);
            console.log(`DEBUG:   Coefficients (A,B,C,D,E order for func):`, passedCoeffs);
            console.log(`DEBUG:   Target Unit: ${propDef.targetUnitName}, Conversion Factor: ${propDef.conversionFactor}`);
            console.log(`DEBUG:   Raw Coeffs from DB for this property:`, propData); // Log the whole propData for this key
        }
        // --- END DEBUGGING BLOCK ---

        let Tc_K_series: number | null = null;
        if (propDef.requiresTc) Tc_K_series = currentCompoundData.criticalTemp ?? null;
        
        let mw_kg_kmol_series: number | null = null;
        if (propDef.requiresMolarMass) mw_kg_kmol_series = currentCompoundData.molarWeight ?? null;

        let displayUnit = propDef.targetUnitName;
        let conversionFactorForDisplay = 1;

        // This conversion to 'bar' is specific to vapor pressure in your yAxisConfig.
        // For calculation, we use the raw value and apply targetUnitName conversion later if needed.
        if (propDef.targetUnitName === "Pa" && propDef.jsonKey === "Vapour pressure") { // Make sure this condition is specific if needed
            displayUnit = "bar";
            conversionFactorForDisplay = 1 / 100000;
        }
        s._internal_display_unit = displayUnit;

        for (let i = 0; i <= numSteps; i++) {
            const T = finalPlotTminKelvin! + i * tempStep;
            const compoundSpecificTmin = parseCoefficient(propData.Tmin?.value ?? propData.Tmin);
            const compoundSpecificTmax = parseCoefficient(propData.Tmax?.value ?? propData.Tmax);
            
            if (T < (compoundSpecificTmin ?? -Infinity) - 1e-6 || T > (compoundSpecificTmax ?? Infinity) + 1e-6 ) {
                continue;
            }

            let rawValue: number | null = null; // Value before conversionFactor and displayFactor
            
            // --- START DEBUGGING BLOCK for a specific temperature point ---
            // You can pick a T value that you have reference data for, e.g., 373.15 K for 100°C
            const isDebugTempPoint = (Math.abs(T - 373.15) < tempStep / 2); // Check if T is close to 100C
            if ((propDef.jsonKey === "Ideal gas heat capacity" || propDef.jsonKey === "Second virial coefficient") && isDebugTempPoint) {
                console.log(`DEBUG:   Calculating at T=${T.toFixed(2)} K`);
            }
            // --- END DEBUGGING BLOCK ---

            switch (propDef.equationType) {
              case 'eq101':
                if (passedCoeffs.length >= 5) rawValue = calculateEq101(T, passedCoeffs[0]!, passedCoeffs[1]!, passedCoeffs[2]!, passedCoeffs[3]!, passedCoeffs[4] ?? undefined);
                break;
              case 'eq105':
                if (passedCoeffs.length >= 4 && Tc_K_series !== null && mw_kg_kmol_series !== null) {
                     rawValue = calculateEq105(T, passedCoeffs[0]!, passedCoeffs[1]!, Tc_K_series, passedCoeffs[3]!, mw_kg_kmol_series);
                } else if (passedCoeffs.length >=3 && Tc_K_series === null && parseCoefficient(propData.C) !== null && mw_kg_kmol_series !== null) {
                     const tcFromCoeffC = parseCoefficient(propData.C)!; // Assuming 'C' holds Tc if dedicated Tc not available
                     rawValue = calculateEq105(T, passedCoeffs[0]!, passedCoeffs[1]!, tcFromCoeffC, passedCoeffs[2]!, mw_kg_kmol_series);
                }
                break;
              case 'polynomial': // Used by Ideal Gas Heat Capacity if config is not changed
                rawValue = calculatePolynomial(T, passedCoeffs[0]!, passedCoeffs[1] ?? undefined, passedCoeffs[2] ?? undefined, passedCoeffs[3] ?? undefined, passedCoeffs[4] ?? undefined); 
                break;
              case 'eq16_complex': // Should be used by Ideal Gas Heat Capacity if config IS changed
                if (passedCoeffs.length >= 5) rawValue = calculateEq16Complex(T, passedCoeffs[0]!, passedCoeffs[1]!, passedCoeffs[2]!, passedCoeffs[3]!, passedCoeffs[4]!);
                break;
              case 'eq106':
                if (passedCoeffs.length >= 5 && Tc_K_series !== null) rawValue = calculateEq106(T, passedCoeffs[0]!, passedCoeffs[1]!, passedCoeffs[2]!, passedCoeffs[3]!, passedCoeffs[4]!, Tc_K_series);
                break;
              case 'eq102_cv':
                if (passedCoeffs.length >= 4) rawValue = calculateEq102_conductivity_viscosity(T, passedCoeffs[0]!, passedCoeffs[1]!, passedCoeffs[2]!, passedCoeffs[3]!);
                break;
              case 'eq104_virial': // Used by Second Virial Coefficient
                if (passedCoeffs.length >= 5) rawValue = calculateEq104_virial(T, passedCoeffs[0]!, passedCoeffs[1]!, passedCoeffs[2]!, passedCoeffs[3]!, passedCoeffs[4]!);
                break;
            }

            let finalValue = rawValue;
            if (finalValue !== null && isFinite(finalValue)) {
              if (propDef.conversionFactor) {
                // --- START DEBUGGING BLOCK ---
                if ((propDef.jsonKey === "Ideal gas heat capacity" || propDef.jsonKey === "Second virial coefficient") && isDebugTempPoint) {
                    console.log(`DEBUG:     Raw value (from eq): ${rawValue}`);
                    console.log(`DEBUG:     Applying conversionFactor: ${propDef.conversionFactor}`);
                }
                // --- END DEBUGGING BLOCK ---
                finalValue *= propDef.conversionFactor;
              }
              // --- START DEBUGGING BLOCK ---
              if ((propDef.jsonKey === "Ideal gas heat capacity" || propDef.jsonKey === "Second virial coefficient") && isDebugTempPoint) {
                  console.log(`DEBUG:     Value after propDef.conversionFactor: ${finalValue}`);
                  console.log(`DEBUG:     Applying displayConversionFactor: ${conversionFactorForDisplay} (for display unit ${displayUnit})`);
              }
              // --- END DEBUGGING BLOCK ---
              finalValue *= conversionFactorForDisplay; // This is mainly for Pa to bar for vapor pressure display
              // --- START DEBUGGING BLOCK ---
              if ((propDef.jsonKey === "Ideal gas heat capacity" || propDef.jsonKey === "Second virial coefficient") && isDebugTempPoint) {
                  console.log(`DEBUG:     Final value for plot point: ${finalValue}`);
              }
              // --- END DEBUGGING BLOCK ---
              points.push([T, finalValue]);
            }
        }
        s.data = points;
    });

    // --- Debugging: Log all calculated points ---
    const allCalculatedPointsForDebugging = seriesData.map(s => ({
        name: s.name,
        data: (s.data as [number, number][]).map(p => ({ 
            temperatureKelvin: p[0], 
            temperatureCelsius: parseFloat((p[0] - 273.15).toPrecision(5)), // Match tooltip precision idea
            value: p[1] 
        }))
    })).filter(s => s.data.length > 0);
    console.log("Calculated Series Data for Plotting:", JSON.stringify(allCalculatedPointsForDebugging, null, 2));
    // --- End Debugging ---


    // --- determine axis unit dynamically ---
    const yAxisDisplayName = propDef.displayName;
    let yAxisUnit = propDef.targetUnitName;
    if (seriesData.length > 0 && (seriesData[0] as any)._internal_display_unit) {
      yAxisUnit = (seriesData[0] as any)._internal_display_unit;
    }

    const yAxisConfig: any = {
      type: 'value' as const,
      name: `${yAxisDisplayName} (${yAxisUnit})`,
      nameLocation: 'middle' as const,
      nameGap: 65,
      position: 'left' as const,
      axisLine: { show: true, lineStyle: { color: propDef.color, width: 2 } },
      axisLabel: { formatter: (val: number) => formatNumberToPrecision(val, 3), color: '#e0e6f1', fontSize: 14, fontFamily: 'Merriweather Sans' },
      nameTextStyle: { color: '#e0e6f1', fontSize: 16, fontFamily: 'Merriweather Sans' },
      splitLine: { show: false },
      scale: true,
    };
    


    setEchartsOptions({
      backgroundColor: '#08306b', // Changed background color
      title: {
        text: `${propDef.displayName}: ${titleCompoundNames || 'Selected Compounds'}`,
        left: 'center',
        textStyle: { color: '#E5E7EB', fontSize: 18, fontFamily: 'Merriweather Sans' },
      },
      tooltip: {
        trigger: 'axis',
        axisPointer: { type: 'cross' },
        backgroundColor: '#1e293b',
        borderColor: '#3b82f6',
        textStyle: { color: '#e5e7eb', fontFamily: 'Merriweather Sans' },
        formatter: (params: any) => {
            if (!Array.isArray(params) || params.length === 0) return '';
            // Temperature conversion for tooltip display
            const tempInCelsius = params[0].axisValue - 273.15;
            let tooltipHtml = `Temperature: <b>${formatNumberToPrecision(tempInCelsius, 3)} °C</b><br/>`;
            params.forEach((param: any) => {
                const seriesFullname = param.seriesName; 
                tooltipHtml += `${param.marker} ${seriesFullname}: <b>${formatNumberToPrecision(param.value[1], 4)}</b><br/>`;
            });
            return tooltipHtml;
        }
      },
      legend: {
        data: seriesData.map(s => s.name as string), // Legend data is just compound names
        bottom: 5, // Reduced gap to x-axis
        textStyle: { color: '#9ca3af', fontFamily: 'Merriweather Sans', fontSize: 11 },
        inactiveColor: '#4b5563',
        type: 'scroll',
        formatter: (name) => {
            // Find the series to get its internal legend equation
            const series = seriesData.find(s => s.name === name) as any;
            if (series && series._internal_legend_equation) {
                 // Truncate name if too long for legend display, full name in tooltip
                const displayName = name.length > 40 ? name.substring(0, 37) + "..." : name;
                return `${displayName}`; // Equation shown in tooltip
            }
            return name;
        },
        tooltip: {
            show: true,
            formatter: (params: any) => { // params here is { name: string }
                const series = seriesData.find(s => s.name === params.name) as any;
                if (series && series._internal_legend_equation) {
                    return `<strong>${params.name}</strong><br/>Equation: ${series._internal_legend_equation}`;
                }
                return params.name;
            }
        }
      },
      grid: { left: '8%', right: '5%', bottom: '12%', top: '10%', containLabel: true }, // Adjusted grid.bottom
      xAxis: {
        type: 'value',
        name: 'Temperature (°C)', 
        nameLocation: 'middle',
        nameGap: 30,
        axisLabel: { 
            color: '#e0e6f1', 
            fontFamily: 'Merriweather Sans', 
            formatter: (val: number) => formatNumberToPrecision(val - 273.15, 1) // Convert K to °C for display, adjust precision
        },
        nameTextStyle: { color: '#e0e6f1', fontSize: 15, fontFamily: 'Merriweather Sans' },
        axisLine: { lineStyle: { color: '#4b5563', width: 2 } },
        splitLine: { show: false },
        scale: false, // Set to false when min, max, and interval are explicitly defined
        min: finalPlotTminKelvin, // Set x-axis range based on adjusted Kelvin Tmin
        max: finalPlotTmaxKelvin, // Set x-axis range based on adjusted Kelvin Tmax
        interval: celsiusStep, // Set the interval (in Kelvin units, which corresponds to Celsius step)
      },
      yAxis: yAxisConfig,
      series: seriesData.filter(s => (s.data as any[]).length > 0), // Only plot series with data
      toolbox: {
        show: true, orient: 'vertical', right: 10, top: 'center',
        feature: { saveAsImage: { show: true, title: 'Save as Image', backgroundColor: '#0f172a' } },
        iconStyle: { borderColor: '#9ca3af' }
      },
    });

  }, [/* …other deps… */]);  // ← removed useLogScale here

  const handleFetchAndPlot = useCallback(async () => {
    setLoading(true);
    setOverallError(null);
    let allDataFetchedSuccessfully = true;

    const updatedCompounds = await Promise.all(compounds.map(async (compound) => {
        if (!compound.name.trim()) {
            return { ...compound, data: null, error: "Compound name is empty." };
        }
        try {
            const data = await fetchCompoundPropertiesLocal(compound.name);
            return { ...compound, data, error: null };
        } catch (err: any) {
            allDataFetchedSuccessfully = false;
            return { ...compound, data: null, error: err.message || "Failed to fetch data." };
        }
    }));
    setCompounds(updatedCompounds);
    
    // Filter available properties based on ALL fetched compounds
    if (allDataFetchedSuccessfully && updatedCompounds.every(c => c.data)) {
        const commonProps = propertiesToPlotConfig.filter(propDef => {
            return updatedCompounds.every(compoundState => {
                if (!compoundState.data) return false;
                const propData = compoundState.data.properties[propDef.jsonKey];
                if (!propData) return false;
                const Tmin = parseCoefficient(propData.Tmin?.value ?? propData.Tmin);
                const Tmax = parseCoefficient(propData.Tmax?.value ?? propData.Tmax);
                if (Tmin === null || Tmax === null || Tmin >= Tmax) return false;
                // Add more checks if needed (coeffs, Tc, MW)
                if (propDef.requiresTc && compoundState.data.criticalTemp === null) return false;
                if (propDef.requiresMolarMass && compoundState.data.molarWeight === null) return false;
                // Check if essential coefficients are present
                if (propDef.coeffs.length > 0 && parseCoefficient(propData[propDef.coeffs[0]]) === null && propDef.equationType !== 'polynomial') return false; // Basic check for first coeff
                if (propDef.equationType === 'polynomial' && propDef.coeffs.length > 0 && parseCoefficient(propData[propDef.coeffs[0]]) === null) {
                     let firstCoeffPresent = false;
                     for(const cKey of propDef.coeffs){ if(parseCoefficient(propData[cKey]) !== null){ firstCoeffPresent = true; break;}}
                     if(!firstCoeffPresent) return false;
                }
                return true;
            });
        });
        setAvailablePropertiesForSelection(commonProps);
        if (commonProps.length > 0 && !commonProps.find(p => p.jsonKey === selectedPropertyKey)) {
            setSelectedPropertyKey(commonProps[0].jsonKey); // Auto-select first common property if current is not available
            processAndPlotProperties(updatedCompounds, commonProps[0].jsonKey);
        } else if (commonProps.length === 0) {
            setAvailablePropertiesForSelection([]);
            setSelectedPropertyKey('');
            processAndPlotProperties(updatedCompounds, null); // No common properties
            setOverallError("No common plottable properties found for the selected compounds.");
        } else {
             processAndPlotProperties(updatedCompounds, selectedPropertyKey);
        }
    } else {
        // If not all data fetched, plot what's available but indicate errors
        processAndPlotProperties(updatedCompounds, selectedPropertyKey);
        // availablePropertiesForSelection might need to be conservative or cleared
        setAvailablePropertiesForSelection(propertiesToPlotConfig); // Or filter based on what *did* load
    }
    setLoading(false);
  }, [compounds, selectedPropertyKey, processAndPlotProperties]);


  const fetchSuggestionsForCompound = useCallback(debounce(async (compoundId: string, inputValue: string) => {
    if (!inputValue || inputValue.length < 2 || !supabase) {
      setCompounds(prev => prev.map(c => c.id === compoundId ? { ...c, suggestions: [], showSuggestions: false } : c));
      return;
    }
    try {
      const { data, error: fetchError } = await supabase.from('compounds').select('name').ilike('name', `%${inputValue}%`).limit(5);
      if (fetchError) { console.error("Supabase suggestion fetch error:", fetchError); return; }
      setCompounds(prev => prev.map(c => c.id === compoundId ? { ...c, suggestions: data ? data.map(item => item.name) : [], showSuggestions: data && data.length > 0 } : c));
    } catch (err) { console.error("Error fetching suggestions:", err); }
  }, 300), [supabase]); // supabase dependency for useCallback

  const handleCompoundNameChange = (id: string, value: string) => {
    setCompounds(prev => prev.map(c => c.id === id ? { ...c, name: value, data: null, error: null } : c)); // Clear data and error on name change
    if (value.trim() === "") {
      setCompounds(prev => prev.map(c => c.id === id ? { ...c, suggestions: [], showSuggestions: false } : c));
    } else {
      fetchSuggestionsForCompound(id, value);
    }
  };

  const handleSuggestionClick = (compoundId: string, suggestion: string) => {
    setCompounds(prev => prev.map(c => {
      if (c.id === compoundId) {
        c.inputRef.current?.focus();
        return { ...c, name: suggestion, suggestions: [], showSuggestions: false, data: null, error: null };
      }
      return c;
    }));
  };
  
  const addCompoundInput = () => {
    if (compounds.length < 4) {
        setCompounds(prev => [
            ...prev,
            createNewCompoundState() // Use the new function to add compounds
        ]);
    }
  };

  const removeCompoundInput = (idToRemove: string) => {
    setCompounds(prev => prev.filter(c => c.id !== idToRemove));
    // After removing, if a plot was based on this compound, it might need re-evaluation
    // handleFetchAndPlot(); // Optionally re-plot, or let user do it. For now, manual re-plot.
    // If the removed compound was the only one, or crucial for selectedPropertyKey, clear plot or re-evaluate available properties.
    // This logic can be complex, for now, let's assume user will re-plot.
    // A simple re-plot if data existed:
    const remainingCompounds = compounds.filter(c => c.id !== idToRemove);
    if (remainingCompounds.some(c => c.data)) {
        processAndPlotProperties(remainingCompounds, selectedPropertyKey);
    } else {
        setEchartsOptions({}); // Clear chart if no data left
        setAvailablePropertiesForSelection(propertiesToPlotConfig); // Reset available props
    }
  };


  useEffect(() => {
    function handleClickOutside(event: MouseEvent) {
      compounds.forEach(c => {
        const target = event.target as Node;
        if (c.suggestionsRef.current && !c.suggestionsRef.current.contains(target) && c.inputRef.current && !c.inputRef.current.contains(target)) {
          setCompounds(prev => prev.map(pc => pc.id === c.id ? { ...pc, showSuggestions: false } : pc));
        }
      });
    }
    document.addEventListener("mousedown", handleClickOutside);
    return () => { document.removeEventListener("mousedown", handleClickOutside); };
  }, [compounds]);

  // Initial fetch for default compound
  useEffect(() => {
    if (supabase && compounds.length > 0 && compounds[0].name && !compounds[0].data) {
        handleFetchAndPlot();
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [supabase]); // Trigger once when supabase is ready for the initial compound

  // Re-plot if selectedPropertyKey changes and data is available
  useEffect(() => {
    if (compounds.some(c => c.data) && selectedPropertyKey) { // Check if any compound has data
        processAndPlotProperties(compounds, selectedPropertyKey);
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selectedPropertyKey]); // Removed useLogScale dependency

  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          <div className="lg:col-span-1 space-y-6">
            <Card>
              <CardContent className="space-y-4 pt-6">
                {compounds.map((compound, index) => (
                    <div key={compound.id} className="space-y-1">
                        <div className="flex items-center gap-2"> {/* Flex container for label and input group */}
                            <Label htmlFor={compound.id} className="text-sm whitespace-nowrap w-28"> {/* Label on the left */}
                                Compound {index + 1}:
                            </Label>
                            <div className="relative flex items-center gap-2 flex-grow"> {/* Input and trash icon group */}
                                <Input
                                    ref={compound.inputRef}
                                    id={compound.id}
                                    value={compound.name}
                                    onChange={(e) => handleCompoundNameChange(compound.id, e.target.value)}
                                    onFocus={() => { setCompounds(prev => prev.map(c => c.id === compound.id ? { ...c, showSuggestions: c.suggestions.length > 0 } : c));}}
                                    placeholder={`e.g., ${index === 0 ? 'Ethanol' : 'Water'}`}
                                    autoComplete="off"
                                    className="flex-grow"
                                />
                                {index > 0 && (
                                    <Button variant="ghost" size="icon" onClick={() => removeCompoundInput(compound.id)} title="Remove Compound">
                                        <Trash2 className="h-4 w-4 text-red-500" />
                                    </Button>
                                )}
                                {compound.showSuggestions && compound.suggestions.length > 0 && (
                                    <div ref={compound.suggestionsRef} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1 top-full max-h-48 overflow-y-auto">
                                    {compound.suggestions.map((s, i) => <div key={i} onClick={() => handleSuggestionClick(compound.id, s)} className="px-3 py-2 hover:bg-accent cursor-pointer text-sm">{s}</div>)}
                                    </div>
                                )}
                            </div>
                        </div>
                        {compound.error && <p className="text-xs text-red-500 mt-1 pl-32">{compound.error}</p>} {/* Indent error to align with input */}
                    </div>
                ))}

                {compounds.length < 4 && (
                    <Button onClick={addCompoundInput} variant="outline" className="w-full mt-2">
                        <PlusCircle className="mr-2 h-4 w-4" /> Add Compound to Compare
                    </Button>
                )}
                
                <div className="flex items-center gap-2 pt-4">
                  <Label htmlFor="propertySelect" className="block text-sm font-medium text-gray-300 whitespace-nowrap">
                    Property:
                  </Label>
                  <Select
                    value={selectedPropertyKey}
                    onValueChange={(value) => {
                        setSelectedPropertyKey(value);
                    }}
                    disabled={availablePropertiesForSelection.length === 0}
                  >
                    <SelectTrigger id="propertySelect" className="flex-1">
                      <SelectValue placeholder="Select a property" />
                    </SelectTrigger>
                    <SelectContent>
                      {availablePropertiesForSelection.length > 0 ? availablePropertiesForSelection.map(prop => (
                        <SelectItem key={prop.jsonKey} value={prop.jsonKey}>
                          {prop.displayName}
                        </SelectItem>
                      )) : <SelectItem value="no-props" disabled>No common properties available</SelectItem>}
                    </SelectContent>
                  </Select>
                </div>

                <Button onClick={handleFetchAndPlot} disabled={loading} className="w-full">
                  {loading ? 'Fetching...' : 'Fetch & Plot Properties'}
                </Button>
                {overallError && !loading && <Alert variant="destructive" className="mt-2"><Terminal className="h-4 w-4" /><AlertTitle>Error</AlertTitle><AlertDescription>{overallError}</AlertDescription></Alert>}
              </CardContent>
            </Card>
          </div>

          <div className="lg:col-span-2 space-y-6">
            <Card>
              <CardContent className="py-2">
                <div className="relative h-[600px] md:h-[700px] rounded-md bg-[#0f172a]">
                  {loading && (<div className="absolute inset-0 flex items-center justify-center text-white"><div className="text-center"><div className="mb-2">Loading Properties...</div></div></div>)}
                  {!loading && compounds.every(c => !c.data && !c.error) && !overallError && (<div className="absolute inset-0 flex items-center justify-center text-white">Please enter compound(s) and fetch properties.</div>)}
                  {overallError && !loading && (<div className="absolute inset-0 flex items-center justify-center text-red-400 px-4 text-center">Error: {overallError}</div>)}
                  {!loading && Object.keys(echartsOptions).length > 0 && compounds.some(c=>c.data) && !overallError && (
                    <ReactECharts ref={echartsRef} echarts={echarts} option={echartsOptions} style={{ height: '100%', width: '100%' }} notMerge={true} lazyUpdate={true} />
                  )}
                </div>
                <p className="text-xs text-muted-foreground mt-2 text-center">
                    Note: Property correlations are from database (ChemSep, DWSIM, etc.) and have specific temperature ranges and accuracies.
                    <br/>Hover over legend items for the full equation.
                </p>
              </CardContent>
            </Card>
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}

