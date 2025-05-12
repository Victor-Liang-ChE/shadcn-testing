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
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"; // Added for plot mode toggle

import {
    calculateEq101,
    calculateEq105,
    calculatePolynomial,
    calculateEq106,
    calculateEq102_conductivity_viscosity,
    calculateEq104_virial,
    calculateEq16Complex, // Added import
    calculateEq105_molar, // Added import
    parseCoefficient,
    propertiesToPlotConfig, // Import from property-equations
    type PropertyDefinition, // Import type from property-equations
    calculateEq121, // Added import for eq121
    calculateEq13 // Added import for eq13
} from '@/lib/property-equations';
import {
    type ConstantPropertyDefinition,
    constantPropertiesConfig
} from '@/lib/constant-properties'; // Added import for constants

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


const compoundColors = ['#FAC858', '#ff7f0e', '#2ca02c', '#d62728']; // Colors for up to 4 compounds (Yellow, Orange, Green, Red)

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

// Helper function to render property name with optional symbol and subscripts/superscripts
const renderPropertyName = (displayName: string, symbol?: string): React.ReactNode => {
  let symbolSpan: React.ReactNode | null = null;
  if (symbol) {
    // Aggressively trim leading/trailing whitespace from the initial symbol string
    const baseSymbol = symbol.replace(/^\s+|\s+$/g, '');

    // Regex to find parts like _X or ^Y or normal characters
    const parts = baseSymbol.match(/([^_^+]+|_[^_^+]+|\^[^_^+]+)/g);

    if (parts) {
      symbolSpan = (
        <span style={{ padding: '0px', margin: '0px' }}> {/* Wrap all symbol parts in a single span */}
          {parts.map((part, index) => {
            // Each part from regex might still have spaces if baseSymbol had internal spaces
            // that became leading/trailing for a part. e.g. baseSymbol = "A _B" -> parts ["A ", "_B"]
            // So, aggressively trim each part as well.
            const cleanedPart = part.replace(/^\s+|\s+$/g, '');

            if (cleanedPart.startsWith('_') && cleanedPart.length > 1) {
              // Aggressively trim the content *after* substring
              return <sub key={index} style={{ padding: '0px', margin: '0px' }}>{cleanedPart.substring(1).replace(/^\s+|\s+$/g, '')}</sub>;
            } else if (cleanedPart.startsWith('^') && cleanedPart.length > 1) {
              // Aggressively trim the content *after* substring here as well
              return <sup key={index} style={{ padding: '0px', margin: '0px' }}>{cleanedPart.substring(1).replace(/^\s+|\s+$/g, '')}</sup>;
            }
            // cleanedPart is already aggressively trimmed
            return cleanedPart; 
          })}
        </span>
      );
    } else if (baseSymbol) { // Fallback if regex doesn't match but baseSymbol has content
      // baseSymbol is already aggressively trimmed from the start
      symbolSpan = <span style={{ padding: '0px', margin: '0px' }}>{baseSymbol}</span>;
    }
  }

  return (
    <>
      {displayName}
      {symbolSpan && (
        <>,{symbolSpan}</> 
      )}
    </>
  );
};

// Helper function to convert symbol string to HTML string
const renderSymbolToHtml = (symbol?: string): string => {
  if (!symbol) return '';
  const baseSymbol = symbol.replace(/^\s+|\s+$/g, '');
  const parts = baseSymbol.match(/([^_^+]+|_[^_^+]+|\^[^_^+]+)/g);

  if (parts) {
    return parts.map(part => {
      const cleanedPart = part.replace(/^\s+|\s+$/g, '');
      if (cleanedPart.startsWith('_') && cleanedPart.length > 1) {
        return `<sub>${cleanedPart.substring(1).replace(/^\s+|\s+$/g, '')}</sub>`;
      } else if (cleanedPart.startsWith('^') && cleanedPart.length > 1) {
        return `<sup>${cleanedPart.substring(1).replace(/^\s+|\s+$/g, '')}</sup>`;
      }
      return cleanedPart;
    }).join('');
  }
  return baseSymbol;
};


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
    createNewCompoundState('Water') // Use the new function for initial state
  ]);
  
  const [loading, setLoading] = useState(false);
  // Corrected useState syntax for overallError
  const [overallError, setOverallError] = useState<string | null>(null); 
  const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
  const echartsRef = useRef<ReactECharts | null>(null);
  const yAxisUnitRef = useRef<string | null>(null); // Ref to store yAxisUnit for CSV export

  const [selectedPropertyKey, setSelectedPropertyKey] = useState<string>(propertiesToPlotConfig[0].jsonKey);
  const [selectedUnit, setSelectedUnit] = useState<string>(() => {
    const initialPropDef = propertiesToPlotConfig.find(p => p.jsonKey === propertiesToPlotConfig[0].jsonKey);
    return initialPropDef?.availableUnits?.[0]?.unit || initialPropDef?.targetUnitName || '';
  });
  const [availablePropertiesForSelection, setAvailablePropertiesForSelection] = useState<PropertyDefinition[]>(propertiesToPlotConfig);
 
  // New state for plot mode and constants
  const [plotMode, setPlotMode] = useState<'tempDependent' | 'constants'>('tempDependent');
  const [selectedConstantKey, setSelectedConstantKey] = useState<string>(() => constantPropertiesConfig[0]?.jsonKey || '');
  const [selectedConstantUnit, setSelectedConstantUnit] = useState<string>(() => {
    const initialConstDef = constantPropertiesConfig.find(c => c.jsonKey === (constantPropertiesConfig[0]?.jsonKey || ''));
    return initialConstDef?.availableUnits?.[0]?.unit || initialConstDef?.targetUnitName || '';
  });
  const [availableConstantsForSelection, setAvailableConstantsForSelection] = useState<ConstantPropertyDefinition[]>(constantPropertiesConfig);


  const fetchCompoundPropertiesLocal = async (compoundName: string): Promise<FetchedCompoundData | null> => {
    if (!supabase) throw new Error("Supabase client not initialized.");
    // setLoading(true); // Loading state handled by caller
    // setError(null); // Error state handled by caller (per compound)
    console.log(`CompoundProperties: Fetching data for ${compoundName}...`);

    try {
      const { data: compoundDbData, error: compoundError } = await supabase
        .from('compounds').select('id, name, cas_number, molecular_weight').ilike('name', compoundName).limit(1).single(); // Added molecular_weight
      
      console.log(`DEBUG_FETCH_MW: For compound "${compoundName}", compoundDbData:`, JSON.stringify(compoundDbData)); // Log entire compoundDbData

      if (compoundError || !compoundDbData) {
        console.error(`DEBUG_FETCH_MW: Error or no data for compound "${compoundName}" from 'compounds' table. Error:`, compoundError);
        throw new Error(compoundError?.message || `Compound '${compoundName}' not found in 'compounds' table.`);
      }
      
      const compoundId = compoundDbData.id;
      const foundName = compoundDbData.name;
      const directMolarWeight = compoundDbData.molecular_weight; 
      console.log(`DEBUG_FETCH_MW: For compound "${foundName}" (ID: ${compoundId}), directMolarWeight from DB: ${directMolarWeight} (type: ${typeof directMolarWeight})`);


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

      // Prioritize directMolarWeight, then fallback to properties JSON
      let molarWeight: number | null = null;
      if (directMolarWeight !== null && typeof directMolarWeight === 'number') {
        molarWeight = directMolarWeight;
        console.log(`DEBUG_FETCH_MW: Using directMolarWeight for "${foundName}": ${molarWeight}`);
      } else {
        const mwFromProps = properties["Molecular weight"];
        console.log(`DEBUG_FETCH_MW: directMolarWeight for "${foundName}" is null or not a number (value: ${directMolarWeight}). Fallback to properties["Molecular weight"]:`, mwFromProps);
        molarWeight = parseCoefficient(mwFromProps?.value ?? mwFromProps);
        console.log(`DEBUG_FETCH_MW: Parsed molarWeight from properties for "${foundName}": ${molarWeight}`);
      }
      
      const criticalTemp = parseCoefficient(properties["Critical temperature"]?.value ?? properties["Critical temperature"]);
      console.log(`DEBUG_FETCH_MW: Final molarWeight for "${foundName}" being returned: ${molarWeight}, CriticalTemp: ${criticalTemp}`);

      return { properties, molarWeight, criticalTemp, name: foundName };
    } catch (err: any) {
      console.error(`Error fetching data for ${compoundName}:`, err.message);
      // setError(`Data fetch failed for ${compoundName}: ${err.message}`); // Caller will set this per compound
      throw err; // Re-throw to be caught by caller
    }
  };

  const processAndPlotProperties = useCallback((
    allCompoundsData: CompoundInputState[],
    propertyKey: string | null, // This will be for temp-dependent or constant key based on plotMode
    currentSelectedUnit: string // Pass selectedUnit as an argument
  ) => {
    setOverallError(null);
    if (!propertyKey) {
      setEchartsOptions({});
      setOverallError(plotMode === 'tempDependent' ? "Please select a temperature-dependent property to plot." : "Please select a constant property to plot.");
      return;
    }

    if (plotMode === 'tempDependent') {
        const propDef = propertiesToPlotConfig.find(p => p.jsonKey === propertyKey);
        if (!propDef) {
          console.warn(`Property definition for key ${propertyKey} not found.`);
          setEchartsOptions({});
          setOverallError(`Property definition for ${propertyKey} not found.`);
          return;
        }

        const unitDefToUse = propDef.availableUnits?.find(u => u.unit === currentSelectedUnit) || 
                             propDef.availableUnits?.[0] || 
                             { unit: propDef.targetUnitName, conversionFactorFromBase: 1 };
        const displayUnit = unitDefToUse.unit;

        const seriesData: EChartsSeriesOption[] = [];
        let commonTmin: number | null = null;
        let commonTmax: number | null = null;
        let atLeastOneCompoundHasData = false;
        let titleCompoundNames = allCompoundsData.filter(c => c.data).map(c => c.data!.name).join(' vs ');
        if (!titleCompoundNames && allCompoundsData.length > 0 && allCompoundsData[0].name) {
            titleCompoundNames = allCompoundsData[0].name; 
        }

        allCompoundsData.forEach((compoundState, compoundIndex) => {
            if (!compoundState.data || !compoundState.data.properties) return;
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
            if (commonTmin === null || Tmin > commonTmin) commonTmin = Tmin;
            if (commonTmax === null || Tmax < commonTmax) commonTmax = Tmax;

            const rawCoeffs = propDef.coeffs.map(cKey => ({ key: cKey, value: parseCoefficient(propData[cKey]) }));
            let Tc_K: number | null = null;
            if (propDef.requiresTc) {
                Tc_K = currentCompoundData.criticalTemp ?? null;
                if (Tc_K === null) {
                    console.warn(`Critical Temperature (Tc) required for ${propDef.displayName} for ${currentCompoundData.name} but not found.`);
                    return;
                }
            }
            const mw_kg_kmol_series: number | null = currentCompoundData.molarWeight ?? null;
            let equationBody = propDef.equationTemplate || "";
            rawCoeffs.forEach(c => {
                if (c.value !== null) {
                    equationBody = equationBody.replace(new RegExp(`\\b${c.key}\\b`, 'g'), formatNumberToPrecision(c.value, 2));
                } else {
                    equationBody = equationBody.replace(new RegExp(`\\b${c.key}\\b`, 'g'), "0"); 
                }
            });
            if (Tc_K !== null) {
                const tcCelsiusDisplay = formatNumberToPrecision(Tc_K - 273.15, 2);
                equationBody = equationBody.replace(/\bTc\b/g, `${tcCelsiusDisplay}°C`); 
            }
            if (mw_kg_kmol_series !== null) {
                equationBody = equationBody.replace(/\bMW\b/g, formatNumberToPrecision(mw_kg_kmol_series, 4));
            }
            equationBody = equationBody.replace(/T\^(\d+)/g, 'T<sup>$1</sup>');
            equationBody = equationBody.replace(/T\^([A-Za-z]+)/g, 'T<sup>$1</sup>');
            equationBody = equationBody.replace(/\b0\s*\*?\s*T<sup>[^<]+<\/sup>/g, "0");
            equationBody = equationBody.replace(/\b0\s*\*?\s*ln\(T\)/g, "0");
            equationBody = equationBody.replace(/\b0\/T\b/g, "0");
            equationBody = equationBody.replace(/\b0\s*\*?\s*T\b(?!\w)/g, "0");
            equationBody = equationBody.replace(/\+\s*0\b(?![.\d])/g, ""); 
            equationBody = equationBody.replace(/-\s*0\b(?![.\d])/g, ""); 
            equationBody = equationBody.replace(/([(,=])\s*0\s*\+\s*/g, "$1"); 
            equationBody = equationBody.replace(/\+\s*-\s*/g, '- '); 
            equationBody = equationBody.replace(/-\s*-\s*/g, '+ '); 
            equationBody = equationBody.replace(/\+\s*\+\s*/g, '+ '); 
            equationBody = equationBody.replace(/([(,=])\s*\+\s*/g, '$1'); 
            equationBody = equationBody.trim();
            equationBody = equationBody.replace(/\(\s+/g, '(');
            equationBody = equationBody.replace(/\s+\)/g, ')');
            if (equationBody === "" || equationBody === "()") equationBody = "0";
            equationBody = equationBody.replace(/exp\(\s*\)/g, "exp(0)");
            equationBody = equationBody.replace(/ln\(\s*\)/g, "ln(0)");
            const legendEquation = renderSymbolToHtml(propDef.symbol) + " = " + equationBody;

            seriesData.push({
                name: `${currentCompoundData.name}`,
                type: 'line',
                data: [], 
                yAxisIndex: 0,
                smooth: true,
                lineStyle: { color: compoundColors[compoundIndex % compoundColors.length], width: 2.5 },
                itemStyle: { color: compoundColors[compoundIndex % compoundColors.length] },
                symbol: 'none', 
                _internal_prop_data: propData,
                _internal_compound_data: currentCompoundData,
                _internal_legend_equation: legendEquation,
            } as any); 
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

        const commonTminCelsius = commonTmin - 273.15;
        const commonTmaxCelsius = commonTmax - 273.15;
        const rangeCelsius = commonTmaxCelsius - commonTminCelsius;
        let celsiusStep;
        if (rangeCelsius <= 20) celsiusStep = 2; 
        else if (rangeCelsius <= 50) celsiusStep = 5;
        else if (rangeCelsius <= 100) celsiusStep = 10;
        else if (rangeCelsius <= 250) celsiusStep = 25;
        else if (rangeCelsius <= 500) celsiusStep = 50;
        else celsiusStep = 100;
        if (rangeCelsius > 0 && celsiusStep <= 0) celsiusStep = 1; 
        let finalPlotTminCelsius = Math.floor(commonTminCelsius / celsiusStep) * celsiusStep;
        let finalPlotTmaxCelsius = Math.ceil(commonTmaxCelsius / celsiusStep) * celsiusStep;
        if (finalPlotTminCelsius === finalPlotTmaxCelsius) finalPlotTmaxCelsius = finalPlotTminCelsius + celsiusStep;
        if (finalPlotTmaxCelsius <= finalPlotTminCelsius) finalPlotTmaxCelsius = finalPlotTminCelsius + celsiusStep;
        const finalPlotTminKelvin = finalPlotTminCelsius + 273.15;
        const finalPlotTmaxKelvin = finalPlotTmaxCelsius + 273.15;

        seriesData.forEach(series => {
            const s = series as any;
            const propData = s._internal_prop_data;
            const currentCompoundData = s._internal_compound_data as FetchedCompoundData;
            const points: [number, number][] = [];
            const dataPointCelsiusStep = 0.5;
            const rawCoeffs = propDef.coeffs.map(cKey => ({ key: cKey, value: parseCoefficient(propData[cKey]) }));
            const passedCoeffs = rawCoeffs.map(c => c.value);
            let Tc_K_series: number | null = null;
            if (propDef.requiresTc) Tc_K_series = currentCompoundData.criticalTemp ?? null;
            const mw_kg_kmol_series: number | null = currentCompoundData.molarWeight ?? null;
            s._internal_display_unit = displayUnit;

            for (let currentC = finalPlotTminCelsius; currentC <= finalPlotTmaxCelsius + 1e-9; currentC += dataPointCelsiusStep) {
                const T = currentC + 273.15;
                const compoundSpecificTmin = parseCoefficient(propData.Tmin?.value ?? propData.Tmin);
                const compoundSpecificTmax = parseCoefficient(propData.Tmax?.value ?? propData.Tmax);
                if (T < (compoundSpecificTmin ?? -Infinity) - 1e-6 || T > (compoundSpecificTmax ?? Infinity) + 1e-6 ) continue;
                let rawValue: number | null = null;
                switch (propDef.equationType) {
                  case 'eq101':
                    if (passedCoeffs.length >= 5) rawValue = calculateEq101(T, passedCoeffs[0]!, passedCoeffs[1]!, passedCoeffs[2]!, passedCoeffs[3]!, passedCoeffs[4] ?? undefined);
                    break;
                  case 'eq105':
                    if (passedCoeffs.length >= 4 && Tc_K_series !== null && mw_kg_kmol_series !== null) {
                         rawValue = calculateEq105(T, passedCoeffs[0]!, passedCoeffs[1]!, Tc_K_series, passedCoeffs[3]!, mw_kg_kmol_series);
                    } else if (passedCoeffs.length >=3 && Tc_K_series === null && parseCoefficient(propData.C) !== null && mw_kg_kmol_series !== null) {
                         const tcFromCoeffC = parseCoefficient(propData.C)!; 
                         rawValue = calculateEq105(T, passedCoeffs[0]!, passedCoeffs[1]!, tcFromCoeffC, passedCoeffs[2]!, mw_kg_kmol_series);
                    }
                    break;
                  case 'eq105_molar':
                    if (passedCoeffs.length >= 4 && Tc_K_series !== null) { 
                        const coeffA = passedCoeffs[0]!
                        const coeffB = passedCoeffs[1]!
                        const coeffD_for_eq105_molar = passedCoeffs[3];
                        if (coeffD_for_eq105_molar !== null && coeffD_for_eq105_molar !== undefined) {
                             rawValue = calculateEq105_molar(T, coeffA, coeffB, Tc_K_series, coeffD_for_eq105_molar);
                        } else {
                            console.warn(`eq105_molar: Coefficient D is missing or null/undefined for ${currentCompoundData.name}. D_val: ${coeffD_for_eq105_molar}`);
                        }
                    }
                    break;
                  case 'polynomial':
                    rawValue = calculatePolynomial(T, passedCoeffs[0]!, passedCoeffs[1] ?? undefined, passedCoeffs[2] ?? undefined, passedCoeffs[3] ?? undefined, passedCoeffs[4] ?? undefined); 
                    break;
                  case 'eq16_complex': 
                    if (passedCoeffs.length >= 5) rawValue = calculateEq16Complex(T, passedCoeffs[0]!, passedCoeffs[1]!, passedCoeffs[2]!, passedCoeffs[3]!, passedCoeffs[4]!);
                    break;
                  case 'eq106':
                    if (passedCoeffs.length >= 5 && Tc_K_series !== null) rawValue = calculateEq106(T, passedCoeffs[0]!, passedCoeffs[1]!, passedCoeffs[2]!, passedCoeffs[3]!, passedCoeffs[4]!, Tc_K_series);
                    break;
                  case 'eq102_cv':
                    if (passedCoeffs.length >= 4) rawValue = calculateEq102_conductivity_viscosity(T, passedCoeffs[0]!, passedCoeffs[1]!, passedCoeffs[2]!, passedCoeffs[3]!);
                    break;
                  case 'eq104_virial':
                    if (passedCoeffs.length >= 5) rawValue = calculateEq104_virial(T, passedCoeffs[0]!, passedCoeffs[1]!, passedCoeffs[2]!, passedCoeffs[3]!, passedCoeffs[4]!);
                    break;
                  case 'eq121':
                    if (passedCoeffs.length >= 4) rawValue = calculateEq121(T, passedCoeffs[0]!, passedCoeffs[1]!, passedCoeffs[2]!, passedCoeffs[3]!);
                    break;
                  case 'eq13':
                    if (passedCoeffs.length >= 3) rawValue = calculateEq13(T, passedCoeffs[0]!, passedCoeffs[1]!, passedCoeffs[2]!);
                    break;
                }

                let valueInBaseUnit = rawValue;
                if (valueInBaseUnit !== null && isFinite(valueInBaseUnit)) {
                  if (propDef.conversionFactor && propDef.conversionFactor !== 1) {
                    valueInBaseUnit *= propDef.conversionFactor;
                  }
                  let finalDisplayValue = valueInBaseUnit;
                  const cFactorBase = unitDefToUse.conversionFactorFromBase;
                  if (typeof cFactorBase === 'number') {
                    if (cFactorBase !== 1) finalDisplayValue *= cFactorBase;
                  } else if (typeof cFactorBase === 'object' && cFactorBase !== null) {
                    if (mw_kg_kmol_series === null && (cFactorBase.operation === 'divide_by_mw' || cFactorBase.operation === 'multiply_by_mw')) {
                        console.warn(`Molar mass required for unit conversion to ${displayUnit} for ${currentCompoundData.name}, but not available. Skipping point.`);
                        continue; 
                    }
                    if (mw_kg_kmol_series !== null) {
                        if (cFactorBase.operation === 'divide_by_mw') {
                            if (mw_kg_kmol_series === 0) { console.warn(`Molar mass is zero for ${currentCompoundData.name}, cannot divide. Skipping point for unit ${displayUnit}.`); continue; }
                            finalDisplayValue /= mw_kg_kmol_series;
                        } else if (cFactorBase.operation === 'multiply_by_mw') {
                            finalDisplayValue *= mw_kg_kmol_series;
                        }
                    }
                    if (typeof cFactorBase.factor === 'number' && cFactorBase.factor !== 1) {
                        finalDisplayValue *= cFactorBase.factor;
                    }
                  }
                  points.push([T, finalDisplayValue]);
                }
            }
            s.data = points;
        });

        const finalSeriesToPlot = seriesData.filter(s => (s.data as any[]).length > 0);
        if (seriesData.length > 0 && finalSeriesToPlot.length === 0 && atLeastOneCompoundHasData) {
            console.warn("All series ended up with no data points after processing. This might be due to Tmin/Tmax issues or all points being skipped during calculation/conversion.");
        }

        // --- determine axis unit dynamically ---
        const yAxisDisplayName = propDef.displayName;
        let yAxisUnit = displayUnit; // Use the selected display unit
        yAxisUnitRef.current = yAxisUnit; // Store for CSV export

        let yAxisTickFormatter = (val: number) => formatNumberToPrecision(val, 7); // Increased precision from 5 to 7

        const yAxisConfig: any = {
          type: 'value' as const,
          name: `${yAxisDisplayName} (${yAxisUnit})`,
          nameLocation: 'middle' as const,
          nameGap: propDef.jsonKey === "Liquid viscosity (RPS)" ? 60 : 80, // Adjusted nameGap for RPS
          position: 'left' as const,
          axisLine: { show: true, lineStyle: { color: propDef.color, width: 2 } },
          axisLabel: { 
            formatter: yAxisTickFormatter, // Use the conditional formatter
            color: '#e0e6f1', 
            fontSize: 16, 
            fontFamily: 'Merriweather Sans' 
          },
          nameTextStyle: { color: '#e0e6f1', fontSize: 15, fontFamily: 'Merriweather Sans' }, // Font size 15
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
            axisPointer: { 
              type: 'cross',
              label: {
                backgroundColor: '#ffffff', // Background of the label box
                formatter: function (params: { value: number | string | Date; axisDimension: string; }) { 
                  let valueAsNumber: number;
                  const paramValue = params.value; 
                  const axisDim = params.axisDimension;

                  if (typeof paramValue === 'number') {
                    valueAsNumber = paramValue;
                  } else if (typeof paramValue === 'string') {
                    valueAsNumber = parseFloat(paramValue);
                  } else if (paramValue instanceof Date) { // paramValue can be Date
                    valueAsNumber = paramValue.getTime(); 
                  } else {
                    const valStr = String(paramValue);
                    return valStr.length > 10 ? valStr.substring(0,7) + '...' : valStr;
                  }

                  if (isNaN(valueAsNumber)) {
                      const valStr = String(paramValue); 
                      return valStr.length > 10 ? valStr.substring(0,7) + '...' : valStr;
                  }

                  if (axisDim === 'x') {
                    return formatNumberToPrecision(valueAsNumber - 273.15, 3) + ' °C';
                  } else if (axisDim === 'y') {
                    // Removed conditional formatting for y-axis pointer based on property
                    // if (propDef.jsonKey === "Relative static permittivity") {
                    //   return (valueAsNumber === 0 ? '0' : valueAsNumber.toExponential(2)) + (yAxisUnit ? ' ' + yAxisUnit : '');
                    // }
                    return formatNumberToPrecision(valueAsNumber, 3) + (yAxisUnit ? ' ' + yAxisUnit : '');
                  }
                  return formatNumberToPrecision(valueAsNumber, 3); // Fallback
                },
                // Moved text styling properties directly under label
                fontFamily: 'Merriweather Sans', // Font for the text inside the label
                color: '#333' // Color of the text inside the label
              }
            },
            backgroundColor: '#1e293b', // Background of the main tooltip box
            borderColor: '#3b82f6',
            textStyle: { color: '#e5e7eb', fontFamily: 'Merriweather Sans' },
            formatter: (params: any) => {
                if (!Array.isArray(params) || params.length === 0) return '';
                // Temperature conversion for tooltip display
                const tempInCelsius = params[0].axisValue - 273.15;
                let tooltipHtml = `Temperature: <b>${formatNumberToPrecision(tempInCelsius, 3)} °C</b><br/>`;
                params.forEach((param: any) => {
                    const seriesFullname = param.seriesName; 
                    // Add yAxisUnit to the tooltip for the dependent variable
                    tooltipHtml += `${param.marker} ${seriesFullname}: <b>${formatNumberToPrecision(param.value[1], 4)} ${yAxisUnit || ''}</b><br/>`;
                });
                return tooltipHtml;
            }
          },
          legend: {
            data: finalSeriesToPlot.map(s => s.name as string), // Use filtered series for legend data
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
          grid: { left: '8%', right: '5%', bottom: '12%', top: '10%', containLabel: true },
          xAxis: {
            type: 'value',
            name: 'Temperature (°C)', 
            nameLocation: 'middle',
            nameGap: 30,
            axisLabel: { 
                color: '#e0e6f1', 
                fontFamily: 'Merriweather Sans', 
                fontSize: 16, // Font size 16
                formatter: (val: number) => (val - 273.15).toFixed(0) // Changed to toFixed(0) for integer display
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
          series: finalSeriesToPlot, // Only plot series with data
          toolbox: {
            show: true, orient: 'vertical', right: 10, top: 'center',
            feature: { 
              saveAsImage: { show: true, title: 'Save as Image', backgroundColor: '#0f172a' },
              // TODO: Add data view later if needed
            },
            iconStyle: { borderColor: '#9ca3af' }
          },
        });

    } else if (plotMode === 'constants') {
        const constDef = constantPropertiesConfig.find(c => c.jsonKey === propertyKey);
        if (!constDef) {
          console.warn(`Constant definition for key ${propertyKey} not found.`);
          setEchartsOptions({});
          setOverallError(`Constant definition for ${propertyKey} not found.`);
          return;
        }

        const unitDefToUse = constDef.availableUnits?.find(u => u.unit === currentSelectedUnit) ||
                             constDef.availableUnits?.[0] ||
                             { unit: constDef.targetUnitName, conversionFactorFromBase: 1 };
        const displayUnit = unitDefToUse.unit;
        const conversionFactorForDisplay = unitDefToUse.conversionFactorFromBase; // This will be 1 for K, °C, °F from temp props
        yAxisUnitRef.current = displayUnit; // Store for CSV export

        const barChartData: { name: string, value: number, itemStyle: { color: string } }[] = [];
        const compoundNames: string[] = [];
        let atLeastOneCompoundHasData = false;

        allCompoundsData.forEach((compoundState, compoundIndex) => {
            if (!compoundState.data || !compoundState.data.properties) return;
            
            const currentCompoundData = compoundState.data;
            let rawValue = parseCoefficient(currentCompoundData.properties[constDef.jsonKey]);

            // Special case for Molecular Weight: prefer direct value from `compounds` table if available
            if (constDef.jsonKey === "Molecular weight" && currentCompoundData.molarWeight !== null) {
                rawValue = currentCompoundData.molarWeight;
            }

            if (rawValue !== null && isFinite(rawValue)) {
                atLeastOneCompoundHasData = true;
                let displayValue = rawValue; // Start with raw value (which is in targetUnitName, e.g., K for temperatures)
                
                // Apply a general numeric conversionFactorFromBase if it's provided and not 1
                // For K, °C, °F temperature units, conversionFactorForDisplay is set to 1 in constant-properties.ts,
                // so this block won't apply to them, specific logic below will.
                if (typeof conversionFactorForDisplay === 'number' && conversionFactorForDisplay !== 1) {
                    displayValue *= conversionFactorForDisplay;
                }

                // Special handling for temperature units if targetUnitName is K
                if (constDef.targetUnitName === "K") {
                    if (displayUnit === "°C") {
                        displayValue = rawValue - 273.15; // rawValue is in K
                    } else if (displayUnit === "°F") {
                        displayValue = (rawValue - 273.15) * 9/5 + 32; // rawValue is in K
                    }
                    // If displayUnit is "K", displayValue (initially rawValue) is already correct.
                }


                barChartData.push({
                    name: currentCompoundData.name,
                    value: displayValue,
                    itemStyle: { color: compoundColors[compoundIndex % compoundColors.length] }
                });
                compoundNames.push(currentCompoundData.name);
            } else {
                 console.warn(`Constant ${constDef.displayName} not found or invalid for ${currentCompoundData.name}. Value: ${rawValue}`);
            }
        });
        
        if (!atLeastOneCompoundHasData) {
            setEchartsOptions({});
            setOverallError(`No data available for constant '${constDef.displayName}' for the selected compounds.`);
            return;
        }
        
        const yAxisDisplayName = constDef.displayName;
        const yAxisConfig: any = {
            type: 'value' as const,
            name: `${yAxisDisplayName} (${displayUnit})`,
            nameLocation: 'middle' as const,
            nameGap: 50, // Adjust as needed
            axisLine: { show: true, lineStyle: { color: constDef.color || '#EE6666', width: 2 } },
            axisLabel: { formatter: (val: number) => formatNumberToPrecision(val, 4), color: '#e0e6f1', fontSize: 16, fontFamily: 'Merriweather Sans' },
            nameTextStyle: { color: '#e0e6f1', fontSize: 15, fontFamily: 'Merriweather Sans' },
            splitLine: { show: false },
            scale: true,
        };

        setEchartsOptions({
            backgroundColor: '#08306b',
            title: {
                text: `${constDef.displayName} for ${allCompoundsData.filter(c => c.data).map(c => c.data!.name).join(', ')}`,
                left: 'center',
                textStyle: { color: '#E5E7EB', fontSize: 18, fontFamily: 'Merriweather Sans' },
            },
            tooltip: {
                trigger: 'item', // Changed to item for bar chart
                formatter: (params: any) => {
                    return `${params.name}<br/>${constDef.displayName}: <b>${formatNumberToPrecision(params.value, 4)} ${displayUnit}</b>`;
                },
                backgroundColor: '#1e293b',
                borderColor: '#3b82f6',
                textStyle: { color: '#e5e7eb', fontFamily: 'Merriweather Sans' },
            },
            legend: { show: false }, // No legend needed for single series bar chart by default
            grid: { left: '3%', right: '4%', bottom: '3%', containLabel: true },
            xAxis: {
                type: 'category',
                data: compoundNames,
                axisLabel: { color: '#e0e6f1', fontFamily: 'Merriweather Sans', fontSize: 14, interval: 0, rotate: compoundNames.length > 3 ? 30 : 0 },
                axisLine: { lineStyle: { color: '#4b5563', width: 2 } },
            },
            yAxis: yAxisConfig,
            series: [{
                name: constDef.displayName,
                type: 'bar',
                data: barChartData.map(d => d.value), // Echarts expects values directly for simple bar
                itemStyle: {
                    // If you want each bar to have the color defined in barChartData
                     color: (params: any) => barChartData[params.dataIndex].itemStyle.color
                },
                label: {
                    show: true,
                    position: 'top',
                    formatter: (params: any) => formatNumberToPrecision(params.value, 3),
                    color: '#e0e6f1'
                }
            }],
            toolbox: { show: true, orient: 'vertical', right: 10, top: 'center', feature: { saveAsImage: { show: true, title: 'Save as Image', backgroundColor: '#0f172a' }, }, iconStyle: { borderColor: '#9ca3af' } },
        });
    }
  }, [plotMode, /* …other deps from temp dependent mode… */]);

  const handleFetchAndPlot = useCallback(async () => {
    setLoading(true);
    setOverallError(null);
    console.log("handleFetchAndPlot: Starting fetch and plot process.");
    const fetchedCompoundStates = await Promise.all(compounds.map(async (compound) => {
        if (!compound.name.trim()) {
            return { ...compound, data: null, error: compound.error || "Compound name is empty." };
        }
        try {
            const data = await fetchCompoundPropertiesLocal(compound.name);
            console.log(`handleFetchAndPlot: Fetched data for ${compound.name}:`, data ? 'Success' : 'Failed/Null');
            return { ...compound, data, error: null };
        } catch (err: any) {
            console.log(`handleFetchAndPlot: Error fetching data for ${compound.name}: ${err.message}`);
            return { ...compound, data: null, error: err.message || "Failed to fetch data." };
        }
    }));
    setCompounds(fetchedCompoundStates);
    console.log("handleFetchAndPlot: All compound data fetched (or attempted):", fetchedCompoundStates.map(c => ({ name: c.name, hasData: !!c.data, error: c.error })));
    
    const activeCompoundsWithData = fetchedCompoundStates.filter(c => c.data && c.name.trim());
    
    if (plotMode === 'tempDependent') {
        let commonProps: PropertyDefinition[] = [];
        if (activeCompoundsWithData.length > 0) {
            commonProps = propertiesToPlotConfig.filter(propDef => {
                return activeCompoundsWithData.every(compoundState => {
                    const propData = compoundState.data!.properties[propDef.jsonKey];
                    if (!propData) return false;
                    const Tmin = parseCoefficient(propData.Tmin?.value ?? propData.Tmin);
                    const Tmax = parseCoefficient(propData.Tmax?.value ?? propData.Tmax);
                    if (Tmin === null || Tmax === null || Tmin >= Tmax) return false;
                    if (propDef.requiresTc && compoundState.data!.criticalTemp === null) return false;
                    if (propDef.requiresMolarMass && compoundState.data!.molarWeight === null) return false;
                    const hasMwDependentUnit = propDef.availableUnits?.some(unit => typeof unit.conversionFactorFromBase === 'object' && (unit.conversionFactorFromBase.operation === 'divide_by_mw' || unit.conversionFactorFromBase.operation === 'multiply_by_mw'));
                    if (hasMwDependentUnit && compoundState.data!.molarWeight === null) return false;
                    if (propDef.coeffs.length > 0) {
                        if (propDef.equationType !== 'polynomial') {
                            const firstCoeffKey = propDef.coeffs[0];
                            const firstCoeffVal = parseCoefficient(propData[firstCoeffKey]);
                            if (firstCoeffVal === null) return false;
                            if (propDef.equationType === 'eq105_molar') {
                                const coeffBVal = parseCoefficient(propData['B']);
                                if (coeffBVal === null || coeffBVal <= 0) return false;
                            }
                        } else {
                            if (!propDef.coeffs.some(cKey => parseCoefficient(propData[cKey]) !== null)) return false;
                        }
                    }
                    return true;
                });
            });
        }
        setAvailablePropertiesForSelection(commonProps);
        // let currentPropertyKeyForPlot = selectedPropertyKey; // No longer needed here
        // let currentUnitForPlot = selectedUnit; // No longer needed here
        if (commonProps.length > 0) {
            const currentSelectionStillValid = commonProps.find(p => p.jsonKey === selectedPropertyKey);
            if (!currentSelectionStillValid) {
                const newPropDef = commonProps[0];
                // currentPropertyKeyForPlot = newPropDef.jsonKey; // No longer needed here
                // currentUnitForPlot = newPropDef.availableUnits?.[0]?.unit || newPropDef.targetUnitName || ''; // No longer needed here
                setSelectedPropertyKey(newPropDef.jsonKey);
                setSelectedUnit(newPropDef.availableUnits?.[0]?.unit || newPropDef.targetUnitName || '');
            }
        } else {
            setSelectedPropertyKey(''); setSelectedUnit(''); // currentPropertyKeyForPlot = ''; currentUnitForPlot = ''; // No longer needed here
            if (fetchedCompoundStates.filter(c => c.name.trim()).length > 0) setOverallError("No common plottable temperature-dependent properties found.");
        }
        // REMOVED: processAndPlotProperties(fetchedCompoundStates, currentPropertyKeyForPlot || null, currentUnitForPlot);

    } else { // plotMode === 'constants'
        let commonConsts: ConstantPropertyDefinition[] = [];
        if (activeCompoundsWithData.length > 0) {
            commonConsts = constantPropertiesConfig.filter(constDef => {
                return activeCompoundsWithData.every(compoundState => {
                    // Ensure compoundState.data and compoundState.data.properties exist
                    if (!compoundState.data || !compoundState.data.properties) {
                        return false; 
                    }
                    // Special case for Molecular Weight: prefer direct value from `compounds` table
                    if (constDef.jsonKey === "Molecular weight" && compoundState.data.molarWeight !== null) {
                        // Assuming molarWeight, if present, is a valid finite number.
                        return true;
                    }
                    const constData = compoundState.data.properties[constDef.jsonKey];
                    const parsedValue = parseCoefficient(constData);
                    // Ensure the parsed value is a non-null finite number to be considered common/plottable
                    return parsedValue !== null && isFinite(parsedValue);
                });
            });
        }
        setAvailableConstantsForSelection(commonConsts);
        // let currentConstantKeyForPlot = selectedConstantKey; // No longer needed here
        // let currentConstantUnitForPlot = selectedConstantUnit; // No longer needed here

        if (commonConsts.length > 0) {
            const currentSelectionStillValid = commonConsts.find(c => c.jsonKey === selectedConstantKey);
            if (!currentSelectionStillValid) {
                const newConstDef = commonConsts[0];
                // currentConstantKeyForPlot = newConstDef.jsonKey; // No longer needed here
                // currentConstantUnitForPlot = newConstDef.availableUnits?.[0]?.unit || newConstDef.targetUnitName || ''; // No longer needed here
                setSelectedConstantKey(newConstDef.jsonKey);
                setSelectedConstantUnit(newConstDef.availableUnits?.[0]?.unit || newConstDef.targetUnitName || '');
            }
        } else {
            setSelectedConstantKey(''); setSelectedConstantUnit(''); // currentConstantKeyForPlot = ''; currentConstantUnitForPlot = ''; // No longer needed here
             if (fetchedCompoundStates.filter(c => c.name.trim()).length > 0) setOverallError("No common constant properties found for the selected compounds.");
        }
        // REMOVED: processAndPlotProperties(fetchedCompoundStates, currentConstantKeyForPlot || null, currentConstantUnitForPlot);
    }
    
    setLoading(false);
  }, [compounds, selectedPropertyKey, selectedUnit, fetchCompoundPropertiesLocal, plotMode, selectedConstantKey, selectedConstantUnit, /* removed processAndPlotProperties from here if it was, but it's used by the effect */]);


  const handleExportCSV = useCallback(() => {
    if (!echartsOptions || !echartsOptions.series || (echartsOptions.series as EChartsSeriesOption[]).length === 0) {
      console.warn("No data to export.");
      return;
    }
    const seriesForCSV = echartsOptions.series as EChartsSeriesOption[];
    if (!seriesForCSV.some(s => s.data && (s.data as any[]).length > 0)) {
        console.warn("No data points in series to export.");
        return;
    }

    let filename = "exported_data.csv";
    let csvContent = '\uFEFF'; // UTF-8 BOM

    if (plotMode === 'tempDependent') {
        const currentPropDef = propertiesToPlotConfig.find(p => p.jsonKey === selectedPropertyKey);
        const propDisplayName = currentPropDef ? currentPropDef.displayName : 'Property';
        filename = `${propDisplayName.replace(/\s+/g, '_')}_vs_Temperature.csv`;
        const unitSuffix = yAxisUnitRef.current ? ` (${yAxisUnitRef.current})` : '';
        const headers = ['Temperature (°C)'];
        seriesForCSV.forEach(s => headers.push(`${s.name} ${propDisplayName}${unitSuffix}`));
        csvContent += headers.join(',') + '\r\n';
        const allKelvinTemps = new Set<number>();
        seriesForCSV.forEach(s => { if (s.data) { (s.data as [number, number][]).forEach(p => allKelvinTemps.add(p[0])); } });
        const sortedKelvinTemps = Array.from(allKelvinTemps).sort((a, b) => a - b);
        sortedKelvinTemps.forEach(tk => {
            const tc = tk - 273.15;
            const row = [formatNumberToPrecision(tc, 3)];
            seriesForCSV.forEach(s => {
                let valueFound = '';
                if (s.data) {
                    const point = (s.data as [number, number][]).find(p => Math.abs(p[0] - tk) < 1e-9);
                    if (point) valueFound = formatNumberToPrecision(point[1], 4);
                }
                row.push(valueFound);
            });
            csvContent += row.join(',') + '\r\n';
        });
    } else { // plotMode === 'constants'
        const currentConstDef = constantPropertiesConfig.find(c => c.jsonKey === selectedConstantKey);
        const constDisplayName = currentConstDef ? currentConstDef.displayName : 'Constant';
        filename = `${constDisplayName.replace(/\s+/g, '_')}_Data.csv`;
        const unitSuffix = yAxisUnitRef.current ? ` (${yAxisUnitRef.current})` : '';
        
        // Assuming seriesForCSV[0].data contains objects like { name: string, value: number } for bar chart
        // Or, if it's simpler, xAxis.data for names and series[0].data for values
        const xAxisData = echartsOptions.xAxis as any; // Assuming single xAxis
        const compoundNames = xAxisData?.data as string[] || [];
        const constantValues = (seriesForCSV[0]?.data as any[])?.map(item => typeof item === 'object' ? item.value : item) || [];

        const headers = ['Compound', `${constDisplayName}${unitSuffix}`];
        csvContent += headers.join(',') + '\r\n';

        compoundNames.forEach((name, index) => {
            const value = constantValues[index] !== undefined ? formatNumberToPrecision(constantValues[index], 4) : '';
            csvContent += `"${name.replace(/"/g, '""')}",${value}\r\n`;
        });
    }

    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement("a");
    if (link.download !== undefined) {
        const url = URL.createObjectURL(blob);
        link.setAttribute("href", url);
        link.setAttribute("download", filename);
        link.style.visibility = 'hidden';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }
  }, [echartsOptions, selectedPropertyKey, yAxisUnitRef]);


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
        processAndPlotProperties(remainingCompounds, selectedPropertyKey, selectedUnit);
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
    if (supabase && compounds.length > 0 && compounds[0].name && !compounds[0].data && !loading) { // Added !loading to prevent re-trigger if already in progress
        handleFetchAndPlot();
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [supabase, compounds[0]?.name, compounds[0]?.data]); // Refined dependencies for initial fetch

  // Re-plot if selectedPropertyKey/selectedConstantKey changes and data is available
  useEffect(() => {
    if (loading) return; // Do not re-plot if a fetch operation is in progress

    const currentKey = plotMode === 'tempDependent' ? selectedPropertyKey : selectedConstantKey;
    const currentUnit = plotMode === 'tempDependent' ? selectedUnit : selectedConstantUnit;
    
    let keyIsValidForCurrentAvailableList = false;
    if (currentKey) {
      if (plotMode === 'tempDependent') {
        keyIsValidForCurrentAvailableList = !!availablePropertiesForSelection.find(p => p.jsonKey === currentKey);
      } else { // plotMode === 'constants'
        keyIsValidForCurrentAvailableList = !!availableConstantsForSelection.find(c => c.jsonKey === currentKey);
      }
    }

    // Only plot if data exists for at least one compound, a key is selected, and that key is valid for the current available list.
    // If no key is selected (currentKey is empty), but data exists, processAndPlotProperties will handle showing an appropriate message.
    if (compounds.some(c => c.data)) {
        if (currentKey && keyIsValidForCurrentAvailableList) {
            processAndPlotProperties(compounds, currentKey, currentUnit);
        } else if (!currentKey) { 
            // If no key is selected (e.g. after filtering, no common props found, or initial state before selection)
            // Call processAndPlotProperties with null key to clear chart or show "select property"
            processAndPlotProperties(compounds, null, currentUnit);
        }
    } else if (!compounds.some(c => c.data) && Object.keys(echartsOptions).length > 0 && !overallError) {
        // If there's no compound data (e.g., all compounds removed or failed to load),
        // and there was a chart, clear it.
        // This case might be covered if selectedKey becomes null/empty.
        // processAndPlotProperties(compounds, null, currentUnit); // This will show "Please select property"
        // Or more directly:
        // setEchartsOptions({});
        // if (!overallError) setOverallError("No data to plot. Please add compounds and fetch properties.");
        // Let processAndPlotProperties(..., null, ...) handle this.
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [compounds, selectedPropertyKey, selectedUnit, selectedConstantKey, selectedConstantUnit, plotMode, loading, availablePropertiesForSelection, availableConstantsForSelection, processAndPlotProperties, overallError]); // Added overallError

  const canExportCSV = echartsOptions && 
                       echartsOptions.series && 
                       (echartsOptions.series as EChartsSeriesOption[]).length > 0 &&
                       (echartsOptions.series as EChartsSeriesOption[]).some(s => s.data && (s.data as any[]).length > 0);

  const handleInputKeyDown = (event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key === 'Enter') {
      event.preventDefault(); // Prevent default form submission or other Enter key behaviors
      handleFetchAndPlot();
    }
  };

  // Generate a key for the Select component to force re-render when options change
  const propertySelectKey = availablePropertiesForSelection.map(p => p.jsonKey).join(',') || 'no-props';
  const constantSelectKey = availableConstantsForSelection.map(c => c.jsonKey).join(',') || 'no-const-props';


  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                   <div className="lg:col-span-1 space-y-6">
            <Card>
              <CardContent className="space-y-4 pt-6">
                <div className="flex items-center space-x-2 mb-4">
                  {/* Replace Switch with Tabs */}
                  <Tabs value={plotMode} onValueChange={(value) => setPlotMode(value as 'tempDependent' | 'constants')} className="w-full">
                    <TabsList className="grid w-full grid-cols-2">
                      <TabsTrigger value="tempDependent">Temp. Dependent</TabsTrigger>
                      <TabsTrigger value="constants">Constants</TabsTrigger>
                    </TabsList>
                  </Tabs>
                </div>

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
                                    onKeyDown={handleInputKeyDown} // Added event handler
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
                
                <div className="flex items-start gap-2 pt-4">
                  {plotMode === 'tempDependent' ? (
                    <>
                      <div className="flex-grow space-y-2">
                        <Select
                          key={propertySelectKey}
                          value={selectedPropertyKey}
                          onValueChange={(value) => {
                              setSelectedPropertyKey(value);
                              const newPropDef = availablePropertiesForSelection.find(p => p.jsonKey === value);
                              const newUnit = newPropDef?.availableUnits?.[0]?.unit || newPropDef?.targetUnitName || '';
                              setSelectedUnit(newUnit);
                          }}
                          disabled={availablePropertiesForSelection.length === 0}
                        >
                          <SelectTrigger id="propertySelect" className="w-full">
                            <SelectValue placeholder="Select a property" />
                          </SelectTrigger>
                          <SelectContent>
                            {availablePropertiesForSelection.length > 0 ? availablePropertiesForSelection.map(prop => (
                              <SelectItem key={prop.jsonKey} value={prop.jsonKey}>
                                {renderPropertyName(prop.displayName, prop.symbol)}
                              </SelectItem>
                            )) : <SelectItem value="no-props" disabled>No common properties</SelectItem>}
                          </SelectContent>
                        </Select>
                      </div>
                      
                      {(() => {
                        const currentPropDef = propertiesToPlotConfig.find(p => p.jsonKey === selectedPropertyKey);
                        const unitsForCurrentProp = currentPropDef?.availableUnits;
                        if (unitsForCurrentProp && unitsForCurrentProp.length > 1) {
                          return (
                            <div className="flex-grow space-y-2 max-w-[150px]">
                              <Select
                                value={selectedUnit}
                                onValueChange={(value) => setSelectedUnit(value)}
                              >
                                <SelectTrigger id="unitSelect" className="w-full">
                                  <SelectValue placeholder="Unit" />
                                </SelectTrigger>
                                <SelectContent>
                                  {unitsForCurrentProp.map(u => (
                                    <SelectItem key={u.unit} value={u.unit}>
                                      {u.displayName || u.unit}
                                    </SelectItem>
                                  ))}
                                </SelectContent>
                              </Select>
                            </div>
                          );
                        }
                        return null;
                      })()}
                    </>
                  ) : ( // plotMode === 'constants'
                    <>
                      <div className="flex-grow space-y-2">
                        <Select
                          key={constantSelectKey}
                          value={selectedConstantKey}
                          onValueChange={(value) => {
                              setSelectedConstantKey(value);
                              const newConstDef = availableConstantsForSelection.find(c => c.jsonKey === value);
                              const newUnit = newConstDef?.availableUnits?.[0]?.unit || newConstDef?.targetUnitName || '';
                              setSelectedConstantUnit(newUnit);
                          }}
                          disabled={availableConstantsForSelection.length === 0}
                        >
                          <SelectTrigger id="constantSelect" className="w-full">
                            <SelectValue placeholder="Select a constant" />
                          </SelectTrigger>
                          <SelectContent>
                            {availableConstantsForSelection.length > 0 ? availableConstantsForSelection.map(prop => (
                              <SelectItem key={prop.jsonKey} value={prop.jsonKey}>
                                {renderPropertyName(prop.displayName, prop.symbol)}
                              </SelectItem>
                            )) : <SelectItem value="no-const-props" disabled>No common constants</SelectItem>}
                          </SelectContent>
                        </Select>
                      </div>
                      
                      {(() => {
                        const currentConstDef = constantPropertiesConfig.find(c => c.jsonKey === selectedConstantKey);
                        const unitsForCurrentConst = currentConstDef?.availableUnits;
                        // Only show unit selector if there's more than one unit OR if the single unit is not "-" (dimensionless)
                        if (unitsForCurrentConst && (unitsForCurrentConst.length > 1 || (unitsForCurrentConst.length === 1 && unitsForCurrentConst[0].unit !== "-"))) {
                          return (
                            <div className="flex-grow space-y-2 max-w-[150px]">
                              <Select
                                value={selectedConstantUnit}
                                onValueChange={(value) => setSelectedConstantUnit(value)}
                              >
                                <SelectTrigger id="constantUnitSelect" className="w-full">
                                  <SelectValue placeholder="Unit" />
                                </SelectTrigger>
                                <SelectContent>
                                  {unitsForCurrentConst.map(u => (
                                    <SelectItem key={u.unit} value={u.unit}>
                                      {u.displayName || u.unit}
                                    </SelectItem>
                                  ))}
                                </SelectContent>
                              </Select>
                            </div>
                          );
                        }
                        return null;
                      })()}
                    </>
                  )}
                </div>

                <Button onClick={handleFetchAndPlot} disabled={loading} className="w-full">
                  {loading ? 'Fetching...' : 'Fetch & Plot Properties'}
                </Button>
                {canExportCSV && !loading && (
                    <Button onClick={handleExportCSV} variant="outline" className="w-full mt-2">
                        Export Data as CSV
                    </Button>
                )}
                {overallError && !loading && <Alert variant="destructive" className="mt-2"><Terminal className="h-4 w-4" /><AlertTitle>Error</AlertTitle><AlertDescription>{overallError}</AlertDescription></Alert>}
              </CardContent>
            </Card>
          </div>

          <div className="lg:col-span-2 space-y-6">
            <Card>
              <CardContent className="py-2">
                <div className="relative h-[600px] md:h-[700px] rounded-md bg-[#08306b]">
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

