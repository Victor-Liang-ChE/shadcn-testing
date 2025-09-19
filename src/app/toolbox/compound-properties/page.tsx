'use client';

import React, { useState, useEffect, useCallback, useRef } from 'react';
import { createClient, SupabaseClient } from '@supabase/supabase-js';
import { useTheme } from "next-themes";

import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import { LineChart, LineSeriesOption } from 'echarts/charts'; // Added LineSeriesOption
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  ToolboxComponent,
  DataZoomComponent, // Keep DataZoomComponent if you plan to re-add zoom later
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers'; // Added CanvasRenderer
import type { EChartsOption } from 'echarts'; 

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
    calculatePropertyByEquation, // The only direct calculator you need now
    calculateBoilingPointEq101, // Keep for the special boiling point solver
    parseCoefficient,
    propertiesToPlotConfig,
    type PropertyDefinition,
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
// Format numbers with default rounding precision
const formatNumberToPrecision = (num: any, precision: number = 4): string => {
  if (num === null || num === undefined || isNaN(num)) return '';
  const factor = Math.pow(10, precision);
  return String(Math.round(num * factor) / factor);
};
// Format with significant figures and scientific notation for extreme values
const formatWithSigFigs = (num: number, sigFigs: number = 3): string => {
  const absVal = Math.abs(num);
  if ((absVal !== 0 && absVal < 1e-6) || absVal >= 1e6) {
    return num.toExponential(sigFigs - 1);
  }
  return num.toPrecision(sigFigs);
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

// Debounce function removed for instantaneous suggestions

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

// Helper function (can be placed inside CompoundPropertiesPage or at an accessible scope)
function calculateDynamicAxisParams(
    dataMinY: number,
    dataMaxY: number,
    baseTotalGapForName: number // Represents the original desired distance from axis line to axis name text
) {
    let dynamicAxisLabelMargin = 8; // Default ECharts margin for axisLabel if no data
    let maxAbsVal = -1;

    if (isFinite(dataMinY) && isFinite(dataMaxY)) {
        maxAbsVal = Math.max(Math.abs(dataMinY), Math.abs(dataMaxY));
        if (dataMinY <= 0 && dataMaxY >= 0) { // If range crosses zero
             maxAbsVal = Math.max(maxAbsVal, 0); // Consider 0 if it's the max abs value
        }
    } else if (isFinite(dataMinY)) {
        maxAbsVal = Math.abs(dataMinY);
    } else if (isFinite(dataMaxY)) {
        maxAbsVal = Math.abs(dataMaxY);
    }
    // else maxAbsVal remains -1, dynamicAxisLabelMargin uses default.

    if (maxAbsVal !== -1) { // If we have a valid max absolute value
        if (maxAbsVal >= 1000000) {      // For values like 1,000,000s or millions
            dynamicAxisLabelMargin = 5;
        } else if (maxAbsVal >= 100000) { // For values in the 100,000s range
            dynamicAxisLabelMargin = 5;
        } else if (maxAbsVal >= 10000) { // For values in the 10,000s range
            dynamicAxisLabelMargin = 10;
        } else if (maxAbsVal >= 1000) {  // For values in the 1,000s range
            dynamicAxisLabelMargin = 15;
        } else if (maxAbsVal >= 100) {   // For values in the 100s range
            dynamicAxisLabelMargin = 20;
        } else if (maxAbsVal >= 10) {    // For values in the 10s range
            dynamicAxisLabelMargin = 25;
        } else {                         // For values less than 10 (including 0 and fractionals)
            if (maxAbsVal < 1 && maxAbsVal > 0) { // Fractional values that might be long when formatted (e.g., 0.0000123)
                const order = Math.floor(Math.log10(maxAbsVal));
                if (order <= -4) {       // e.g., 0.000x or smaller
                    dynamicAxisLabelMargin = 0;
                } else if (order <= -2) {  // e.g., 0.0x
                    dynamicAxisLabelMargin = 15;
                } else {                   // e.g., 0.x
                    dynamicAxisLabelMargin = 18;
                }
            } else {                       // For values from 0 up to 9.99...
                dynamicAxisLabelMargin = 22; // Smallest margin
            }
        }
    }

    let calculatedNameGap = baseTotalGapForName - dynamicAxisLabelMargin;
    // Ensure nameGap doesn't become too small or negative
    calculatedNameGap = Math.max(10, calculatedNameGap); // Minimum 10px gap

    return { dynamicAxisLabelMargin, calculatedNameGap };
}

// Define the specific sampling options ECharts uses
type SeriesSamplingOption = 'sum' | 'min' | 'max' | 'none' | 'average' | 'minmax' | 'lttb';

// Define a custom series option type at module scope
interface AppLineSeriesOption extends LineSeriesOption {
  _internal_legend_equation?: string;
  _internal_display_unit?: string;
  large?: boolean; // Explicitly add large
  sampling?: SeriesSamplingOption; // Explicitly add sampling with the correct type
}

export default function CompoundPropertiesPage() {
  const { resolvedTheme } = useTheme(); // Get the resolved theme ('light' or 'dark')
  
  // Add cache reference for efficient data handling
  const propertiesCache = useRef(new Map<string, FetchedCompoundData | null>());
  
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

  const [selectedPropertyKey, setSelectedPropertyKey] = useState<string>(propertiesToPlotConfig[0].displayName); // Changed to displayName
  const [selectedUnit, setSelectedUnit] = useState<string>(() => {
    const initialPropDef = propertiesToPlotConfig[0]; // Simpler: directly use the first definition
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

  // State for manual temperature input and points
  const [manualTempInput, setManualTempInput] = useState<string>('');
  type ManualTempPoint = {
    compoundName: string;
    seriesName: string; // To match ECharts series name
    tempCelsius: number; // For Boiling Point: stores INPUT PRESSURE. For others: INPUT TEMP C.
    tempKelvin: number; // For Boiling Point: not directly used for input. For others: INPUT TEMP K.
    value: number | null; // For Boiling Point: stores CALCULATED TEMP C. For others: CALCULATED PROP VALUE.
    unit: string; // For Boiling Point: stores INPUT PRESSURE UNIT. For others: PROP VALUE UNIT.
    isValid: boolean;
    message?: string;
  };
  const [manualTempPointsData, setManualTempPointsData] = useState<ManualTempPoint[]>([]);

  // State to control when auto-updates should happen
  const [shouldAutoUpdate, setShouldAutoUpdate] = useState(true);
  
  // State to track if user explicitly clicked fetch button (to show empty compound errors)
  const [userClickedFetch, setUserClickedFetch] = useState(false);

  // useEffect for initial data fetching (e.g., for 'Water')
  useEffect(() => {
    compounds.forEach(compound => {
      if (compound.name && !compound.data && !compound.error) {
        // Check if supabase client is available before fetching
        if (!supabase) {
          console.error("Supabase client not initialized. Cannot fetch initial data for", compound.name);
          setCompounds(prev =>
            prev.map(c =>
              c.id === compound.id
                ? { ...c, error: `Supabase client not ready for ${compound.name}` }
                : c
            )
          );
          return;
        }

        const fetchInitialDataForCompound = async (cmpdToFetch: CompoundInputState) => {
          try {
            const fetchedData = await fetchCompoundPropertiesLocal(cmpdToFetch.name);
            setCompounds(prev =>
              prev.map(c =>
                c.id === cmpdToFetch.id
                  ? { ...c, data: fetchedData, error: null }
                  : c
              )
            );
          } catch (err: any) {
            console.error(`Initial fetch failed for ${cmpdToFetch.name}:`, err.message);
            setCompounds(prev =>
              prev.map(c =>
                c.id === cmpdToFetch.id
                  ? { ...c, data: null, error: `Failed to load ${cmpdToFetch.name}: ${err.message}` }
                  : c
              )
            );
          }
        };
        fetchInitialDataForCompound(compound);
      }
    });
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []); // Empty dependency array ensures this runs only once on mount


  const fetchCompoundPropertiesLocal = async (compoundName: string): Promise<FetchedCompoundData | null> => {
    if (!supabase) throw new Error("Supabase client not initialized.");
    
    const trimmedCompoundName = compoundName.trim(); // Trim the compound name
    if (!trimmedCompoundName) { // Check if trimmed name is empty
        console.warn("fetchCompoundPropertiesLocal: Compound name is empty after trimming.");
        throw new Error("Compound name cannot be empty.");
    }

    const cacheKey = trimmedCompoundName.toLowerCase();
    
    // 1. Check cache before fetching
    if (propertiesCache.current.has(cacheKey)) {
        const cachedData = propertiesCache.current.get(cacheKey);
        // Handle the case where a previous fetch for this name failed
        if (cachedData === null) {
            throw new Error(`Compound '${trimmedCompoundName}' not found in the database.`);
        }
        return cachedData || null;
    }

    try {
      const { data: compoundDbData, error: compoundError } = await supabase
        .from('compound_properties').select('name, properties').ilike('name', trimmedCompoundName).limit(1).single();

      if (compoundError || !compoundDbData) {
        throw new Error(compoundError?.message || `Compound '${trimmedCompoundName}' not found in 'compound_properties' table.`);
      }
      
      const foundName = compoundDbData.name;
      const properties = compoundDbData.properties;
      
      if (typeof properties !== 'object' || properties === null) throw new Error(`Invalid properties format for ${foundName}.`);

  // Extract molecular weight and critical temperature from properties
  let molarWeight: number | null = null;
  const mwFromProps = properties.MolecularWeight;
  molarWeight = parseCoefficient(mwFromProps);

  const criticalTemp = parseCoefficient(properties.CriticalTemperature);

      const finalData = { properties, molarWeight, criticalTemp, name: foundName };
      
      propertiesCache.current.set(cacheKey, finalData);
      return finalData;
    } catch (err: any) {
      console.error(`Error fetching data for "${trimmedCompoundName}":`, err.message);
      
      propertiesCache.current.set(cacheKey, null);
      throw err;
    }
  };

  const processAndPlotProperties = useCallback((
    allCompoundsData: CompoundInputState[],
    propertyKey: string | null, // This will be displayName for temp-dependent or jsonKey for constant
    currentSelectedUnit: string, // Pass selectedUnit as an argument
    manualPointsToAnnotate: ManualTempPoint[] // Adapted for annotations
  ) => {
    setOverallError(null);
    if (!propertyKey) {
      setEchartsOptions({});
      
      // Count compounds with data
      const compoundsWithData = allCompoundsData.filter(c => c.data && c.name.trim());
      
      if (compoundsWithData.length === 0) {
        // No compounds with data
        setOverallError(plotMode === 'tempDependent' ? "Please select a temperature-dependent property to plot." : "Please select a constant property to plot.");
      } else if (compoundsWithData.length === 1) {
        // Only one compound
        setOverallError(plotMode === 'tempDependent' ? "No temperature-dependent properties available to display for this compound." : "No constant properties available to display for this compound.");
      } else {
        // Multiple compounds but no common properties
        setOverallError(plotMode === 'tempDependent' ? "No common temperature-dependent properties found between the selected compounds." : "No common constant properties found between the selected compounds.");
      }
      return;
    }

    if (plotMode === 'tempDependent') {
        const propDef = propertiesToPlotConfig.find(p => p.displayName === propertyKey);
        if (!propDef) {
          setEchartsOptions({});
          setOverallError(`Property definition for ${propertyKey} not found.`);
          return;
        }

        const unitDefToUse = propDef.availableUnits?.find(u => u.unit === currentSelectedUnit) || 
                             propDef.availableUnits?.[0] || 
                             { unit: propDef.targetUnitName, conversionFactorFromBase: 1 };
        const displayUnit = unitDefToUse.unit; // For Boiling Point, this is the selected PRESSURE unit for X-axis
        const isBoilingPointPlot = propDef.displayName === "Boiling Point";

        const seriesData: AppLineSeriesOption[] = []; // Use the custom type here
        // let commonTmin: number | null = null; // Removed, replaced by commonTminK_forPlot
        // let commonTmax: number | null = null; // Removed, replaced by commonTmaxK_forPlot
        let atLeastOneCompoundHasData = false;
        let titleCompoundNames = allCompoundsData.filter(c => c.data).map(c => c.data!.name).join(' vs ');
        if (!titleCompoundNames && allCompoundsData.length > 0 && allCompoundsData[0].name) {
            titleCompoundNames = allCompoundsData[0].name; 
        }

        // Variables for common temperature range calculation (for non-Boiling Point plots)
        let commonTminK_forPlot: number | null = null;
        let commonTmaxK_forPlot: number | null = null;
        let finalPlotTminKelvin: number | undefined = undefined;
        let finalPlotTmaxKelvin: number | undefined = undefined;
        let xAxisKelvinInterval: number | undefined = undefined; // Interval in Kelvin for X-axis ticks

        // Define pressure range and step for Boiling Point X-axis (data generation and axis ticks)
        let finalPlotPressureMin: number | undefined = undefined;
        let finalPlotPressureMax: number | undefined = undefined;
        let pressureAxisInterval: number | undefined = undefined; // For X-axis ticks (in displayUnit)

        if (isBoilingPointPlot) {
            // These define the conceptual range of pressures you want to see data for.
            const dataDisplayPressureMin = 0.01; // e.g., 0.01 bar (smallest pressure for the curve)
            const dataDisplayPressureMax = 20;   // e.g., 20 bar (largest pressure for the curve)
            const pressureRangeForStepping = dataDisplayPressureMax - dataDisplayPressureMin;

            // Determine a "nice" step for the pressure axis ticks
            if (pressureRangeForStepping <= 0) {
                pressureAxisInterval = 1;
            } else if (pressureRangeForStepping <= 0.1) pressureAxisInterval = 0.01;
            else if (pressureRangeForStepping <= 0.5) pressureAxisInterval = 0.05;
            else if (pressureRangeForStepping <= 1.0) pressureAxisInterval = 0.1;
            else if (pressureRangeForStepping <= 2.0) pressureAxisInterval = 0.2;
            else if (pressureRangeForStepping <= 5.0) pressureAxisInterval = 0.5;
            else if (pressureRangeForStepping <= 10.0) pressureAxisInterval = 1;
            else if (pressureRangeForStepping <= 20.0) pressureAxisInterval = 2;
            else if (pressureRangeForStepping <= 50.0) pressureAxisInterval = 5;
            else pressureAxisInterval = 10;

            // Calculate axis min/max to align with the interval, similar to temperature axis
            finalPlotPressureMin = Math.floor(dataDisplayPressureMin / pressureAxisInterval) * pressureAxisInterval;
            // If dataDisplayPressureMin is small (e.g. 0.01) and interval is larger (e.g. 0.1),
            // finalPlotPressureMin might become 0. We need data generation to start > 0.
            // The axis itself can start at 0 if appropriate.
            if (finalPlotPressureMin < 0) finalPlotPressureMin = 0; // Ensure axis doesn't start negative

            finalPlotPressureMax = Math.ceil(dataDisplayPressureMax / pressureAxisInterval) * pressureAxisInterval;

            if (finalPlotPressureMin >= finalPlotPressureMax) { // Ensure max > min
                finalPlotPressureMax = finalPlotPressureMin + pressureAxisInterval;
            }
            // `finalPlotPressureMin` and `finalPlotPressureMax` now define the x-axis range for the plot.
        } else { // Not Boiling Point Plot - existing logic for temperature axis
            // First pass to determine common Tmin/Tmax for plotting range
      allCompoundsData.forEach(compoundState => {
                if (compoundState.data && compoundState.data.properties) {
        const propDataForRange = compoundState.data.properties[propDef.jsonKey];
                    if (propDataForRange) {
                        const tMinK_compound = parseCoefficient(propDataForRange.Tmin?.value ?? propDataForRange.Tmin);
                        const tMaxK_compound = parseCoefficient(propDataForRange.Tmax?.value ?? propDataForRange.Tmax);

                        if (tMinK_compound !== null && tMaxK_compound !== null && tMinK_compound < tMaxK_compound) {
                            commonTminK_forPlot = (commonTminK_forPlot === null) ? tMinK_compound : Math.max(commonTminK_forPlot, tMinK_compound);
                            commonTmaxK_forPlot = (commonTmaxK_forPlot === null) ? tMaxK_compound : Math.min(commonTmaxK_forPlot, tMaxK_compound);
                        }
                    }
                }
            });

            if (commonTminK_forPlot === null || commonTmaxK_forPlot === null || commonTminK_forPlot >= commonTmaxK_forPlot) {
                // finalPlotTmin/maxKelvin remain undefined, xAxisKelvinInterval remains undefined. Axis will auto-scale.
            } else {
                const commonTminCelsius_plot = commonTminK_forPlot - 273.15;
                const commonTmaxCelsius_plot = commonTmaxK_forPlot - 273.15;
                const rangeCelsius_plot = commonTmaxCelsius_plot - commonTminCelsius_plot;

                let celsiusStepForAxisTicks: number;
                if (rangeCelsius_plot <= 0) { 
                    celsiusStepForAxisTicks = 10; 
                } else if (rangeCelsius_plot <= 20) celsiusStepForAxisTicks = 2;
                else if (rangeCelsius_plot <= 50) celsiusStepForAxisTicks = 5;
                else if (rangeCelsius_plot <= 100) celsiusStepForAxisTicks = 10;
                else if (rangeCelsius_plot <= 250) celsiusStepForAxisTicks = 25;
                else if (rangeCelsius_plot <= 500) celsiusStepForAxisTicks = 50;
                else celsiusStepForAxisTicks = 100;
                
                if (rangeCelsius_plot > 0 && celsiusStepForAxisTicks <= 0) celsiusStepForAxisTicks = 1;


                let calculatedPlotTminCelsius = Math.floor(commonTminCelsius_plot / celsiusStepForAxisTicks) * celsiusStepForAxisTicks;
                let calculatedPlotTmaxCelsius = Math.ceil(commonTmaxCelsius_plot / celsiusStepForAxisTicks) * celsiusStepForAxisTicks;

                if (calculatedPlotTminCelsius >= calculatedPlotTmaxCelsius) { // Ensure max > min
                    calculatedPlotTmaxCelsius = calculatedPlotTminCelsius + celsiusStepForAxisTicks;
                }
                
                finalPlotTminKelvin = calculatedPlotTminCelsius + 273.15;
                finalPlotTmaxKelvin = calculatedPlotTmaxCelsius + 273.15;
                xAxisKelvinInterval = celsiusStepForAxisTicks; // Interval is in Kelvin, magnitude matches Celsius step
            }
        }

        // Define pressure range for Boiling Point X-axis (data generation)
        // const boilingPointPressureMinDisplay = 0.01; // Replaced by curveDataStartPressure
        // const boilingPointPressureMaxDisplay = 20;   // Replaced by curveDataEndPressure
        // const numPressurePoints = 100; // Replaced by numDataPointsForCurve

        allCompoundsData.forEach((compoundState, compoundIndex) => {
            if (!compoundState.data || !compoundState.data.properties) return;
            atLeastOneCompoundHasData = true;
            const currentCompoundData = compoundState.data;
            const propData = currentCompoundData.properties[propDef.jsonKey];
            if (!propData) {
                return;
            }
            const Tmin_compound_K = parseCoefficient(propData.Tmin?.value ?? propData.Tmin);
            const Tmax_compound_K = parseCoefficient(propData.Tmax?.value ?? propData.Tmax);
            if (Tmin_compound_K === null || Tmax_compound_K === null || Tmin_compound_K >= Tmax_compound_K) {
                return;
            }

            const rawCoeffs = propDef.coeffs.map(cKey => ({ key: cKey, value: parseCoefficient(propData[cKey]) }));
            const passedCoeffs = rawCoeffs.map((c: { key: string, value: number | null }) => c.value);
            let Tc_K_series: number | null = null;
            if (propDef.requiresTc && !isBoilingPointPlot) Tc_K_series = currentCompoundData.criticalTemp ?? null;
            const mw_kg_kmol_series: number | null = currentCompoundData.molarWeight ?? null;
            
            const points: [number, number][] = [];
            const equationString = propDef.equationTemplate 
                ? propDef.equationTemplate.replace(/([A-E])/g, (match, p1) => {
                    const coeffIndex = propDef.coeffs.indexOf(p1);
                    if (coeffIndex !== -1 && passedCoeffs[coeffIndex] !== null) {
                        return formatNumberToPrecision(passedCoeffs[coeffIndex]!, 3);
                    }
                    return match;
                  })
                : 'N/A';

            if (isBoilingPointPlot) {
                if (passedCoeffs.length >= 5 && passedCoeffs.every((c: number | null) => c !== null)) {
                    const [coeffA, coeffB, coeffC, coeffD, coeffE] = passedCoeffs as [number, number, number, number, number | undefined];
                    
                    const xUnitDef = propDef.availableUnits?.find(u => u.unit === displayUnit);
                    const conversionToPaFactor = xUnitDef ? (typeof xUnitDef.conversionFactorFromBase === 'number' ? 1 / xUnitDef.conversionFactorFromBase : 1) : 1;

                    // Determine the actual pressure range for generating curve points
                    // Start from the larger of the calculated plot minimum or a very small positive number (e.g., 0.001 in display units)
                    // to avoid issues with log(0) or P=0 in calculations, especially if finalPlotPressureMin ended up as 0.
                    const curveDataStartPressure = Math.max(finalPlotPressureMin!, 0.001); // Ensure positive start for calc
                    const curveDataEndPressure = finalPlotPressureMax!;

                    const numDataPointsForCurve = 200; // For a smooth curve
                    let dataPointPressureStep = (curveDataEndPressure - curveDataStartPressure) / numDataPointsForCurve;
                    
                    if (dataPointPressureStep <= 0 && curveDataEndPressure > curveDataStartPressure) {
                        // This case should ideally not happen if curveDataEndPressure > curveDataStartPressure
                        dataPointPressureStep = (curveDataEndPressure - curveDataStartPressure) / 2; // Minimal 2 points
                    }
                    
                    if (dataPointPressureStep > 0) {
                        for (let pDisplay = curveDataStartPressure; pDisplay <= curveDataEndPressure + 1e-9; pDisplay += dataPointPressureStep) {
                            const currentPressureInPa = pDisplay * conversionToPaFactor;

                            if (currentPressureInPa <= 0) continue; // Safeguard, should be handled by curveDataStartPressure
                            
                            let initialTempGuess = 373.15; 
                            if (coeffA && coeffB && currentPressureInPa > 0) { // Antoine-like pre-guess
                                const lnP = Math.log(currentPressureInPa);
                                if (coeffA - lnP !== 0) {
                                    const T_antoine_approx = coeffB / (coeffA - lnP);
                                    if (T_antoine_approx > 0 && T_antoine_approx < 1000) { // Plausible temp range
                                        initialTempGuess = T_antoine_approx;
                                    }
                                }
                            }
                            const tempKelvin = calculateBoilingPointEq101(
                                currentPressureInPa,
                                coeffA, coeffB, coeffC, coeffD, coeffE,
                                initialTempGuess, 50, 1e-6
                            );
                            if (tempKelvin !== null && isFinite(tempKelvin) && tempKelvin > 0) {
                                const tempCelsius = tempKelvin - 273.15;
                                points.push([pDisplay, tempCelsius]);
                            }
                        }
                    } else if (curveDataStartPressure === curveDataEndPressure) { // Single point if range is zero
                        const currentPressureInPa = curveDataStartPressure * conversionToPaFactor;
                         if (currentPressureInPa > 0) {
                            let initialTempGuessSingle = 373.15;
                            if (coeffA && coeffB) {
                                const lnPSingle = Math.log(currentPressureInPa);
                                if (coeffA - lnPSingle !== 0) {
                                    const T_antoine_approx_single = coeffB / (coeffA - lnPSingle);
                                    if (T_antoine_approx_single > 0 && T_antoine_approx_single < 1000) {
                                        initialTempGuessSingle = T_antoine_approx_single;
                                    }
                                }
                            }
                            const tempKelvin = calculateBoilingPointEq101(
                                currentPressureInPa,
                                coeffA, coeffB, coeffC, coeffD, coeffE,
                                initialTempGuessSingle, 50, 1e-6
                            );
                            if (tempKelvin !== null && isFinite(tempKelvin) && tempKelvin > 0) {
                                points.push([curveDataStartPressure, tempKelvin - 273.15]);
                            }
                         }
                    }
                }
            } else {
                // Determine loop range for data points for this specific compound
                const loopStartCelsius = (finalPlotTminKelvin !== undefined) ? (finalPlotTminKelvin - 273.15) : (Tmin_compound_K - 273.15);
                const loopEndCelsius = (finalPlotTmaxKelvin !== undefined) ? (finalPlotTmaxKelvin - 273.15) : (Tmax_compound_K - 273.15);
                
                const numDataPointsForLoop = 200; // Desired number of points over the range
                let dataPointCelsiusStep = 1; // Fixed 1°C increment for data points

                for (let currentC = loopStartCelsius; currentC <= loopEndCelsius + 1e-9; currentC += dataPointCelsiusStep) {
                    const T = currentC + 273.15; // T is in Kelvin for calculations
                    
                    // Ensure T is within the specific compound's valid range before calculation
                    if (T < Tmin_compound_K - 1e-6 || T > Tmax_compound_K + 1e-6) {
                        continue;
                    }
                        
                    let rawValue: number | null = null;

                const eqnoParsed = parseCoefficient(propData?.eqno);
                const eqnoStr = eqnoParsed !== null ? String(eqnoParsed) : (propData?.eqno ? String(propData.eqno) : '');
                if (!eqnoStr) {
                    console.warn(`Equation number (eqno) not found for property ${propDef.displayName} for ${currentCompoundData.name}.`);
                    continue; // Skips this iteration of the forEach loop
                }

                const coeffsForCalc = {
                    A: parseCoefficient(propData.A),
                    B: parseCoefficient(propData.B),
                    C: parseCoefficient(propData.C),
                    D: parseCoefficient(propData.D),
                    E: parseCoefficient(propData.E),
                };

                // Tc is required for some equations, so we pass it if available.
                const Tc_K_for_calc = currentCompoundData.criticalTemp ?? undefined;

        rawValue = calculatePropertyByEquation(
          eqnoStr,
                    T, // Temperature in Kelvin
                    coeffsForCalc,
                    Tc_K_for_calc
                );

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
                      // For non-boiling point plots, X is TempK, Y is finalDisplayValue
                      points.push([T, finalDisplayValue]);
                    }
                }
            }
            
            seriesData.push({
                name: currentCompoundData.name,
                type: 'line', // Ensures this is treated as a LineSeriesOption
                data: points,
                lineStyle: { width: 2.5 },
                color: compoundColors[compoundIndex % compoundColors.length],
                symbol: 'none',
                emphasis: { lineStyle: { width: 3.5 } },
                large: true, 
                sampling: 'lttb', 
                _internal_legend_equation: equationString,
                _internal_display_unit: displayUnit,
            });
        });

        const finalSeriesToPlot = seriesData.filter(s => (s.data as any[]).length > 0);

        // Add manual annotations if any
        if (manualPointsToAnnotate && manualPointsToAnnotate.length > 0) {
            finalSeriesToPlot.forEach(series => {
                const seriesName = series.name as string;
                const pointForThisSeries = manualPointsToAnnotate.find(p => p.seriesName === seriesName && p.isValid && p.value !== null);
                if (pointForThisSeries) {
                    if (!series.markPoint) {
                        series.markPoint = { data: [] };
                    }
                    // For Boiling Point: pointForThisSeries.tempCelsius is INPUT PRESSURE, pointForThisSeries.value is CALC TEMP C
                    // For others: pointForThisSeries.tempCelsius is INPUT TEMP C, pointForThisSeries.value is CALC PROP VALUE
                    const xCoord = isBoilingPointPlot ? pointForThisSeries.tempCelsius : (pointForThisSeries.tempCelsius + 273.15);
                    const yCoord = pointForThisSeries.value;

                    (series.markPoint as any).data.push({
                        name: isBoilingPointPlot 
                            ? `Point at ${formatNumberToPrecision(xCoord, 2)} ${pointForThisSeries.unit}`
                            : `Point at ${formatNumberToPrecision(pointForThisSeries.tempCelsius, 2)}°C`,
                        coord: [xCoord, yCoord], 
                        itemStyle: { color: 'red' },
                        symbol: 'circle',
                        symbolSize: 8,
                        label: {
                            show: false, // Changed to false to hide label on graph
                        }
                    });
                }
            });
        }


        // --- Calculate Y-axis data range for dynamic margin ---
        let dataMinY = Infinity;
        let dataMaxY = -Infinity;
        if (finalSeriesToPlot && finalSeriesToPlot.length > 0) {
            finalSeriesToPlot.forEach(series => {
                if (series.data) {
                    (series.data as [number, number][]).forEach(point => {
                        // point[1] is Y-value (TempC for Boiling Point, Prop Value for others)
                        if (point && typeof point[1] === 'number' && isFinite(point[1])) {
                            if (point[1] < dataMinY) dataMinY = point[1];
                            if (point[1] > dataMaxY) dataMaxY = point[1];
                        }
                    });
                }
            });
        }
        // --- End Y-axis data range calculation ---

        const yAxisDisplayName = isBoilingPointPlot ? "Temperature" : propDef.displayName;
        let yAxisUnit = isBoilingPointPlot ? "°C" : displayUnit;
        // yAxisUnitRef.current stores X-axis unit for Boiling Point (pressure unit), or Y-axis unit for others.
        yAxisUnitRef.current = isBoilingPointPlot ? displayUnit : yAxisUnit; 

        // Theme-dependent colors
        const isDark = resolvedTheme === 'dark';
        const textColor = isDark ? 'white' : '#000000';

  // Format Y-axis ticks: scientific notation for extremes, default for <100, 3 sig figs otherwise
  const yAxisTickFormatter = (val: number) => {
    const absVal = Math.abs(val);
    // Scientific notation for very small or large values
    if ((absVal !== 0 && absVal < 1e-6) || absVal >= 1e6) {
      return formatWithSigFigs(val, 3);
    }
    // Default 2-decimal precision for values under 100
    if (absVal < 100) {
      return formatNumberToPrecision(val, 2);
    }
    // Use 3 significant figures for mid-range values
    return formatWithSigFigs(val, 3);
  };

        // --- Determine dynamic axis label margin and nameGap ---
        // For Boiling Point, Y-axis is Temperature (°C).
        // For other plots, Y-axis is the property value.
        const yAxisRangeForGapCalc = { min: dataMinY, max: dataMaxY }; // dataMinY/MaxY are already correct for Y-axis

        const originalBaseTotalGap = (propDef.jsonKey === "Liquid viscosity (RPS)" ? 85 : 75);
        const { dynamicAxisLabelMargin, calculatedNameGap } = calculateDynamicAxisParams(yAxisRangeForGapCalc.min, yAxisRangeForGapCalc.max, originalBaseTotalGap);
        // --- End dynamic parameter calculation ---

        const yAxisConfig: any = {
          type: 'value' as const,
          name: yAxisUnit === '-' ? yAxisDisplayName : `${yAxisDisplayName} (${yAxisUnit})`,
          nameLocation: 'middle' as const,
          nameGap: calculatedNameGap, 
          position: 'left' as const,
          axisLine: { show: true, lineStyle: { color: 'white', width: 2 } },
          axisTick: { show: true, lineStyle: { color: 'white' }, length: 12 },
          axisLabel: { 
            formatter: yAxisTickFormatter, 
            color: 'white', 
            fontSize: 16, 
            fontFamily: 'Merriweather Sans',
            margin: dynamicAxisLabelMargin, 
          },
          nameTextStyle: { color: 'white', fontSize: 15, fontFamily: 'Merriweather Sans' }, 
          splitLine: { show: false },
          scale: true, // Y-axis always auto-scales
          min: undefined, 
          max: undefined, 
        };
        
        const xAxisConfig: any = {
            type: 'value',
            name: isBoilingPointPlot ? `Pressure (${displayUnit})` : 'Temperature (°C)',
            nameLocation: 'middle',
            nameGap: 30,
            axisLabel: {
                color: textColor,
                fontFamily: 'Merriweather Sans',
                fontSize: 16,
                formatter: (val: number) => {
                  // For temp-dependent, convert Kelvin to Celsius for display, else use pressure directly
                  const displayVal = isBoilingPointPlot ? val : val - 273.15;
                  return formatWithSigFigs(displayVal, 3);
                }
            },
            nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
            axisLine: { lineStyle: { color: textColor, width: 2 } },
            axisTick: { show: true, lineStyle: { color: textColor }, length: 12 },
            splitLine: { show: false },
            // Apply calculated range and interval for Boiling Point X-axis
            scale: isBoilingPointPlot ? false : (finalPlotTminKelvin !== undefined ? false : true),
            min: isBoilingPointPlot ? finalPlotPressureMin : finalPlotTminKelvin,
            max: isBoilingPointPlot ? finalPlotPressureMax : finalPlotTmaxKelvin,
            interval: isBoilingPointPlot ? pressureAxisInterval : xAxisKelvinInterval,
        };

        setEchartsOptions({
          // Disable all chart animations when swapping properties
          animation: false,
          animationDuration: 0,
          animationDurationUpdate: 0,
          // Base chart styling
          backgroundColor: 'transparent', 
          title: {
            text: isBoilingPointPlot 
                ? `Boiling Point (Temperature vs. Pressure): ${titleCompoundNames || 'Selected Compounds'}` 
                : `${propDef.displayName}: ${titleCompoundNames || 'Selected Compounds'}`,
            left: 'center',
            textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' },
          },
          tooltip: {
            trigger: 'axis',
            axisPointer: { 
              type: 'cross',
              label: {
                backgroundColor: isDark ? '#1e293b' : '#ffffff',
                formatter: function (params: { value: number | string | Date; axisDimension: string; }) { 
                  let valueAsNumber: number;
                  const paramValue = params.value;
                  const axisDim = params.axisDimension;
                  if (typeof paramValue === 'number') {
                    valueAsNumber = paramValue;
                  } else if (typeof paramValue === 'string') {
                    valueAsNumber = parseFloat(paramValue);
                  } else if (paramValue instanceof Date) {
                    valueAsNumber = paramValue.getTime();
                  } else {
                    const valStr = String(paramValue);
                    return valStr.length > 10 ? valStr.substring(0,7) + '...' : valStr;
                  }
                  const displayVal = isBoilingPointPlot
                    ? (axisDim === 'x' ? valueAsNumber : valueAsNumber)
                    : (axisDim === 'x' ? valueAsNumber - 273.15 : valueAsNumber);
                  const suffix = isBoilingPointPlot
                    ? (axisDim === 'x' ? ` ${displayUnit}` : ' °C')
                    : (axisDim === 'y' ? (yAxisUnit ? ` ${yAxisUnit}` : '') : ' °C');
                  return formatWithSigFigs(displayVal, 3) + suffix;
                },
                fontFamily: 'Merriweather Sans',
                color: textColor
              }
            },
            backgroundColor: resolvedTheme === 'dark' ? '#1e293b' : '#ffffff', // Background of the main tooltip box
            borderColor: resolvedTheme === 'dark' ? '#3b82f6' : '#333333',
            textStyle: { color: textColor, fontFamily: 'Merriweather Sans' },
      formatter: (params: any) => {
        if (!Array.isArray(params) || params.length === 0) return '';
        const sortedParams = [...params].sort((a: any, b: any) => b.value[1] - a.value[1]);
        let tooltipHtml = '';
        if (isBoilingPointPlot) {
          const pressureVal = params[0].axisValue;
          tooltipHtml = `Pressure: <b>${formatWithSigFigs(pressureVal,3)} ${displayUnit}</b><br/>`;
          sortedParams.forEach((param: any) => {
            tooltipHtml += `<span style="color: ${param.color};"><b>${param.seriesName}: ${formatWithSigFigs(param.value[1],3)} °C</b></span><br/>`;
          });
        } else {
          const tempInCelsius = params[0].axisValue - 273.15;
          tooltipHtml = `Temperature: <b>${formatWithSigFigs(tempInCelsius,3)} °C</b><br/>`;
          sortedParams.forEach((param: any) => {
            tooltipHtml += `<span style="color: ${param.color};"><b>${param.seriesName}: ${formatWithSigFigs(param.value[1],3)}${yAxisUnit !== '-' ? ' ' + yAxisUnit : ''}</b></span><br/>`;
          });
        }
        return tooltipHtml;
      }
          },
          legend: {
            data: finalSeriesToPlot.map(s => s.name as string), // Use filtered series for legend data
            bottom: 40, // Closer to x-axis
            textStyle: { color: textColor, fontFamily: 'Merriweather Sans', fontSize: 16 },
            inactiveColor: '#4b5563',
            type: 'scroll',
            itemWidth: 25,
            itemHeight: 2,
            formatter: (name) => {
                const series = seriesData.find(s => s.name === name) as any;
                if (series && series._internal_legend_equation) {
                    const displayName = name.length > 40 ? name.substring(0, 37) + "..." : name;
                    return `${displayName}`; 
                }
                return name;
            },
            tooltip: {
                show: true,
                formatter: (params: any) => { 
                    const series = seriesData.find(s => s.name === params.name) as any;
                    if (series && series._internal_legend_equation) {
                        return `<strong>${params.name}</strong><br/>Equation: ${series._internal_legend_equation}`;
                    }
                    return params.name;
                }
            }
          },
          grid: { left: '5%', right: '5%', bottom: '12%', top: '10%', containLabel: true }, 
          xAxis: xAxisConfig,
          yAxis: yAxisConfig,
          series: finalSeriesToPlot as any, // Use type assertion here
          toolbox: {
            show: false,
          },
        });

    } else if (plotMode === 'constants') {
        // Theme-dependent colors
        const isDark = resolvedTheme === 'dark';
        const textColor = isDark ? 'white' : '#000000';
        
        const constDef = constantPropertiesConfig.find(c => c.jsonKey === propertyKey);
        if (!constDef) {
          setEchartsOptions({});
          setOverallError(`Constant definition for ${propertyKey} not found.`);
          return;
        }

        const unitDefToUse = constDef.availableUnits?.find(u => u.unit === currentSelectedUnit) ||
                             constDef.availableUnits?.[0] ||
                             { unit: constDef.targetUnitName, conversionFactorFromBase: 1 };
        const displayUnit = unitDefToUse.unit;
        const conversionFactorForDisplay = unitDefToUse.conversionFactorFromBase; 
        yAxisUnitRef.current = displayUnit; 

        const barChartData: { name: string, value: number, itemStyle: { color: string } }[] = [];
        const compoundNames: string[] = [];
        let atLeastOneCompoundHasData = false;

        allCompoundsData.forEach((compoundState, compoundIndex) => {
            if (!compoundState.data || !compoundState.data.properties) return;
            
            const currentCompoundData = compoundState.data;
            let rawValue = parseCoefficient(currentCompoundData.properties[constDef.jsonKey]);

            if (constDef.jsonKey === "Molecular weight" && currentCompoundData.molarWeight !== null) {
                rawValue = currentCompoundData.molarWeight;
            }

            if (rawValue !== null && isFinite(rawValue)) {
                atLeastOneCompoundHasData = true;
                let displayValue = rawValue; 
                
                if (typeof conversionFactorForDisplay === 'number' && conversionFactorForDisplay !== 1) {
                    displayValue *= conversionFactorForDisplay;
                }

                if (constDef.targetUnitName === "K") {
                    if (displayUnit === "°C") {
                        displayValue = rawValue - 273.15; 
                    } else if (displayUnit === "°F") {
                        displayValue = (rawValue - 273.15) * 9/5 + 32; 
                    }
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
        
        let dataMinY = Infinity;
        let dataMaxY = -Infinity;
        if (barChartData && barChartData.length > 0) {
            barChartData.forEach(bar => {
                if (typeof bar.value === 'number' && isFinite(bar.value)) {
                    if (bar.value < dataMinY) dataMinY = bar.value;
                    if (bar.value > dataMaxY) dataMaxY = bar.value;
                }
            });
        }
        
        const yAxisDisplayName = constDef.displayName;

        const originalBaseTotalGapConst = (constDef.jsonKey === "Liquid viscosity (RPS)" ? 85 : 75);
        const { dynamicAxisLabelMargin, calculatedNameGap } = calculateDynamicAxisParams(dataMinY, dataMaxY, originalBaseTotalGapConst);

        const yAxisConfig: any = {
            type: 'value' as const,
            name: displayUnit === '-' ? yAxisDisplayName : `${yAxisDisplayName} (${displayUnit})`,
            nameLocation: 'middle' as const,
            nameGap: calculatedNameGap, 
            axisLine: { show: true, lineStyle: { color: 'white', width: 2 } },
            axisTick: { show: true, lineStyle: { color: 'white' }, length: 12 },
            axisLabel: { 
              formatter: (val: number) => formatNumberToPrecision(val, 4), 
              color: 'white', 
              fontSize: 16, 
              fontFamily: 'Merriweather Sans',
              margin: dynamicAxisLabelMargin, 
            },
            nameTextStyle: { color: 'white', fontSize: 15, fontFamily: 'Merriweather Sans' },
            splitLine: { show: false },
            scale: true, // Reverted: Ensure auto-scaling
            min: undefined, // Reverted: Let ECharts determine min
            max: undefined, // Reverted: Let ECharts determine max
        };

        setEchartsOptions({
            backgroundColor: 'transparent',
            title: {
                text: `${constDef.displayName} for ${allCompoundsData.filter(c => c.data && c.name.trim()).map(c => c.data!.name).join(' vs ')}`, // Updated title format
                left: 'center',
                textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' },
            },
            tooltip: {
                show: false,
            },
            legend: { show: false }, 
            grid: { left: '5%', right: '5%', bottom: '12%', top: '10%', containLabel: true }, 
            xAxis: {
                type: 'category',
                data: compoundNames,
                axisLabel: { color: textColor, fontFamily: 'Merriweather Sans', fontSize: 14, interval: 0, rotate: compoundNames.length > 4 ? 30 : 0 }, 
                axisLine: { lineStyle: { color: textColor, width: 2 } },
                axisTick: { show: true, lineStyle: { color: textColor }, length: 12 },
            },
            yAxis: yAxisConfig,
            series: [{
                name: constDef.displayName,
                type: 'bar',
                data: barChartData.map(d => d.value), 
                itemStyle: {
                     color: (params: any) => barChartData[params.dataIndex].itemStyle.color
                },
                label: {
                    show: true,
                    position: 'top',
                    formatter: (params: any) => `${formatNumberToPrecision(params.value, 3)}${displayUnit !== '-' ? ' ' + displayUnit : ''}`,
                    color: textColor,
                    fontSize: Math.max(12, 30 - compoundNames.length * 4),
                    fontFamily: 'Merriweather Sans' // Added font family for bar labels
                }
            }] as any, // Use type assertion here for consistency
            toolbox: { 
              show: false,
            },
        });
    }
  }, [plotMode, resolvedTheme]);

  const handleFetchAndPlot = useCallback(async (isUserInitiated = false) => {
    setLoading(true);
    setOverallError(null);
    setManualTempInput('');
    setManualTempPointsData([]);
    setShouldAutoUpdate(true);
    const fetchedCompoundStates = await Promise.all(compounds.map(async (compound) => {
        if (!compound.name.trim()) {
            return { ...compound, data: null, error: (isUserInitiated ? "Compound name is empty." : null) };
        }
        try {
            const data = await fetchCompoundPropertiesLocal(compound.name);
            return { ...compound, data, error: null };
        } catch (err: any) {
            return { ...compound, data: null, error: err.message || "Failed to fetch data." };
        }
    }));
    setCompounds(fetchedCompoundStates);
    
    const activeCompoundsWithData = fetchedCompoundStates.filter(c => c.data && c.name.trim());
    
    if (plotMode === 'tempDependent') {
        let commonProps: PropertyDefinition[] = [];
        if (activeCompoundsWithData.length > 0) {
      commonProps = propertiesToPlotConfig.filter(propDef => {
                const allCompoundsMeetBasicCriteria = activeCompoundsWithData.every(compoundState => {
        const propData = compoundState.data!.properties[propDef.jsonKey];
                    if (!propData) return false;
                    
                    const TminCompound = parseCoefficient(propData.Tmin?.value ?? propData.Tmin);
                    const TmaxCompound = parseCoefficient(propData.Tmax?.value ?? propData.Tmax);
                    if (TminCompound === null || TmaxCompound === null || TminCompound >= TmaxCompound) return false;

                    if (propDef.requiresTc && compoundState.data!.criticalTemp === null) return false;
                    if (propDef.requiresMolarMass && compoundState.data!.molarWeight === null) return false;
                    const hasMwDependentUnit = propDef.availableUnits?.some(unit => typeof unit.conversionFactorFromBase === 'object' && (unit.conversionFactorFromBase.operation === 'divide_by_mw' || unit.conversionFactorFromBase.operation === 'multiply_by_mw'));
                    if (hasMwDependentUnit && compoundState.data!.molarWeight === null) return false;
                    
          if (propDef.coeffs.length > 0) {
                        // Check if at least one coefficient is available for calculation
                        const hasValidCoeff = propDef.coeffs.some(cKey => parseCoefficient(propData?.[cKey]) !== null);
                        if (!hasValidCoeff) return false;
                        
            // Check if eqno is available and parseable
            if (parseCoefficient(propData?.eqno) === null) return false;
                    }
                    return true;
                });

                if (!allCompoundsMeetBasicCriteria) {
                    return false; 
                }

                let overallCommonTminForProp: number | null = null;
                let overallCommonTmaxForProp: number | any = null;

        for (const compoundState of activeCompoundsWithData) {
          const propData = compoundState.data!.properties[propDef.jsonKey];
                    const tminRawValue = propData?.Tmin?.value ?? propData?.Tmin;
                    const tmaxRawValue = propData?.Tmax?.value ?? propData?.Tmax;

                    const Tmin = parseCoefficient(tminRawValue)!;
                    const Tmax = parseCoefficient(tmaxRawValue)!;

                    if (overallCommonTminForProp === null) { 
                        overallCommonTminForProp = Tmin;
                        overallCommonTmaxForProp = Tmax;
                    } else {
                        overallCommonTminForProp = Math.max(overallCommonTminForProp, Tmin);
                        overallCommonTmaxForProp = Math.min(overallCommonTmaxForProp, Tmax);
                    }
                }
                
                if (overallCommonTminForProp !== null && overallCommonTmaxForProp !== null && overallCommonTminForProp < overallCommonTmaxForProp) {
                    return true; 
                }

                return false; 
            });
        }
        setAvailablePropertiesForSelection(commonProps);
        if (commonProps.length > 0) {
            const currentSelectionStillValid = commonProps.find(p => p.displayName === selectedPropertyKey); // Changed to find by displayName
            if (!currentSelectionStillValid) {
                const newPropDef = commonProps[0];
                setSelectedPropertyKey(newPropDef.displayName); // Changed to set displayName
                setSelectedUnit(newPropDef.availableUnits?.[0]?.unit || newPropDef?.targetUnitName || '');
            }
        } else {
            setSelectedPropertyKey(''); setSelectedUnit(''); 
            if (fetchedCompoundStates.filter(c => c.name.trim()).length > 0) setOverallError("No common plottable temperature-dependent properties found.");
        }

    } else { 
        let commonConsts: ConstantPropertyDefinition[] = [];
        if (activeCompoundsWithData.length > 0) {
            commonConsts = constantPropertiesConfig.filter(constDef => {
                return activeCompoundsWithData.every(compoundState => {
                    if (!compoundState.data || !compoundState.data.properties) {
                        return false; 
                    }
                    if (constDef.jsonKey === "Molecular weight" && compoundState.data.molarWeight !== null) {
                        return true;
                    }
                    const constData = compoundState.data.properties[constDef.jsonKey];
                    const parsedValue = parseCoefficient(constData);
                    return parsedValue !== null && isFinite(parsedValue);
                });
            });
        }
        setAvailableConstantsForSelection(commonConsts);

        if (commonConsts.length > 0) {
            const currentSelectionStillValid = commonConsts.find(c => c.jsonKey === selectedConstantKey);
            if (!currentSelectionStillValid) {
                const newConstDef = commonConsts[0];
                setSelectedConstantKey(newConstDef.jsonKey);
                setSelectedConstantUnit(newConstDef.availableUnits?.[0]?.unit || newConstDef?.targetUnitName || '');
            }
        } else {
            setSelectedConstantKey(''); setSelectedConstantUnit(''); 
             if (fetchedCompoundStates.filter(c => c.name.trim()).length > 0) setOverallError("No common constant properties found for the selected compounds.");
        }
    }
    
    setLoading(false);
  }, [compounds, selectedPropertyKey, selectedUnit, fetchCompoundPropertiesLocal, plotMode, selectedConstantKey, selectedConstantUnit]);

  // Separate handler for button click
  const handleFetchButtonClick = () => {
    setUserClickedFetch(true);
    handleFetchAndPlot(true);
  };

  const handleExportCSV = useCallback(() => {
    if (!echartsOptions || !echartsOptions.series || (echartsOptions.series as AppLineSeriesOption[]).length === 0) { // Use AppLineSeriesOption here
      console.warn("No data to export.");
      return;
    }
    const seriesForCSV = echartsOptions.series as AppLineSeriesOption[]; // Use AppLineSeriesOption here
    if (!seriesForCSV.some(s => s.data && (s.data as any[]).length > 0)) {
        console.warn("No data points in series to export.");
        return;
    }

    let filename = "exported_data.csv";
    let csvContent = '\uFEFF'; // UTF-8 BOM

    if (plotMode === 'tempDependent') {
        const currentPropDef = propertiesToPlotConfig.find(p => p.displayName === selectedPropertyKey); // Changed to find by displayName
        const propDisplayName = currentPropDef ? currentPropDef.displayName : 'Property';
        const isBoilingPointPlotCSV = currentPropDef?.displayName === "Boiling Point";
        
        // yAxisUnitRef.current stores pressure unit for Boiling Point, otherwise the Y-axis unit

        const xAxisUnitForCSV = yAxisUnitRef.current || ''; 
        const xAxisLabelForCSV = isBoilingPointPlotCSV ? `Pressure (${xAxisUnitForCSV})` : 'Temperature (°C)';
        
        filename = isBoilingPointPlotCSV 
            ? `Boiling_Point_Temperature_vs_Pressure.csv`
            : `${propDisplayName.replace(/\s+/g, '_')}_vs_Temperature.csv`;

        const headers = [xAxisLabelForCSV];
        
        seriesForCSV.forEach(s => {
            const seriesYAxisLabel = isBoilingPointPlotCSV ? `Temperature (°C)` : `${propDisplayName}${xAxisUnitForCSV && !isBoilingPointPlotCSV ? ` (${xAxisUnitForCSV})` : ''}`;
            headers.push(`${s.name} ${seriesYAxisLabel}`);
        });
        csvContent += headers.join(',') + '\r\n';

        const allXValues = new Set<number>();
        seriesForCSV.forEach(s => { 
            if (s.data) { 
                (s.data as [number, number][]).forEach(p => allXValues.add(p[0])); 
            } 
        });
        const sortedXValues = Array.from(allXValues).sort((a, b) => a - b);

        sortedXValues.forEach(xVal => {
            let currentXFormatted = '';
            if (isBoilingPointPlotCSV) { // xVal is Pressure
                currentXFormatted = formatNumberToPrecision(xVal, 3);
            } else { // xVal is Temp Kelvin
                currentXFormatted = formatNumberToPrecision(xVal - 273.15, 3);
            }
            const row = [currentXFormatted];

            seriesForCSV.forEach(s => {
                let valueFound = '';
                if (s.data) {
                    const point = (s.data as [number, number][]).find(p => Math.abs(p[0] - xVal) < 1e-9);
                    if (point) {
                         // point[1] is TempC for Boiling Point, or Property Value otherwise
                        valueFound = formatNumberToPrecision(point[1], isBoilingPointPlotCSV ? 3 : 4);
                    }
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
  }, [echartsOptions, plotMode, selectedPropertyKey, selectedConstantKey]); // adjust deps accordingly

  const fetchSuggestionsForCompound = useCallback(async (compoundId: string, inputValue: string) => {
    const trimmedInputValue = inputValue.trim(); 
    if (!trimmedInputValue || trimmedInputValue.length < 2 || !supabase) {
      setCompounds(prev => prev.map(c => c.id === compoundId ? { ...c, suggestions: [], showSuggestions: false } : c));
      return;
    }
    try {
      const { data, error: fetchError } = await supabase
        .from('compound_properties')
        .select('name')
        .ilike('name', `${trimmedInputValue}%`) 
        .limit(5);

  if (fetchError) { console.error("Supabase suggestion fetch error:", fetchError); return; }
  setCompounds(prev => prev.map(c => c.id === compoundId ? { ...c, suggestions: data ? data.map((item: any) => item.name) : [], showSuggestions: data && data.length > 0 } : c));
    } catch (err) { console.error("Error fetching suggestions:", err); }
  }, [supabase]); 

  const handleCompoundNameChange = (id: string, value: string) => {
    setCompounds(prev => prev.map(c => c.id === id ? { ...c, name: value, error: null } : c)); 
    // Disable auto-updates when user is typing compound names
    setShouldAutoUpdate(false);
    if (value.trim() === "") { 
      setCompounds(prev => prev.map(c => c.id === id ? { ...c, suggestions: [], showSuggestions: false } : c));
    } else {
      fetchSuggestionsForCompound(id, value); 
    }
  };

  const handleSuggestionClick = (compoundId: string, suggestion: string) => {
    setCompounds(prev => prev.map(c => {
      if (c.id === compoundId) {
        return { ...c, name: suggestion, suggestions: [], showSuggestions: false, error: null };
      }
      return c;
    }));
  };
  
  const addCompoundInput = () => {
    if (compounds.length < 4) {
        setCompounds(prev => [
            ...prev,
            createNewCompoundState() 
        ]);
    }
  };

  const removeCompoundInput = (idToRemove: string) => {
    // Remove the compound from state and re-enable auto-updates
    setCompounds(prev => prev.filter(c => c.id !== idToRemove));
    setShouldAutoUpdate(true);
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

  useEffect(() => {
    if (supabase && compounds.length > 0 && compounds[0].name && !compounds[0].data && !loading && shouldAutoUpdate) {
      handleFetchAndPlot(false);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [supabase, compounds[0]?.name, compounds[0]?.data, shouldAutoUpdate]);

  // useEffect to handle compound removal - only triggers meaningful changes

  useEffect(() => {
    if (plotMode === 'constants' || plotMode === 'tempDependent') {
      if (compounds.some(c => c.name.trim() !== '')) {
        // Clear manual temp data when plot mode changes or initial fetch for mode
        setManualTempInput('');
        setManualTempPointsData([]);
        // Enable auto-updates when plot mode changes
        setShouldAutoUpdate(true);
        handleFetchAndPlot(false);
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [plotMode]); 


  useEffect(() => {
    if (loading || !shouldAutoUpdate) return; 

    const currentKey = plotMode === 'tempDependent' ? selectedPropertyKey : selectedConstantKey;
    const currentUnit = plotMode === 'tempDependent' ? selectedUnit : selectedConstantUnit;
    
    let keyIsValidForCurrentAvailableList = false;
    if (currentKey) {
      if (plotMode === 'tempDependent') {
        keyIsValidForCurrentAvailableList = !!availablePropertiesForSelection.find(p => p.displayName === currentKey);
      } else { 
        keyIsValidForCurrentAvailableList = !!availableConstantsForSelection.find(c => c.jsonKey === currentKey);
      }
    }

    if (compounds.some(c => c.data)) {
        if (currentKey && keyIsValidForCurrentAvailableList) {
           
            processAndPlotProperties(compounds, currentKey, currentUnit, manualTempPointsData);
        } else if (!currentKey) { 
            processAndPlotProperties(compounds, null, currentUnit, manualTempPointsData);
        }
       } else if (!compounds.some(c => c.data) && Object.keys(echartsOptions).length > 0 && !overallError) {
    }
  }, [compounds.filter(c => c.data).map(c => c.data?.name).join(','), selectedPropertyKey, selectedUnit, selectedConstantKey, selectedConstantUnit, plotMode, loading, availablePropertiesForSelection, availableConstantsForSelection, processAndPlotProperties, overallError, manualTempPointsData, resolvedTheme, shouldAutoUpdate]); 

  const canExportCSV = echartsOptions && 
                       echartsOptions.series && 
                       (echartsOptions.series as AppLineSeriesOption[]).length > 0 && // Use AppLineSeriesOption here
                       (echartsOptions.series as AppLineSeriesOption[]).some(s => s.data && (s.data as any[]).length > 0); // Use AppLineSeriesOption here

  const handleInputKeyDown = (event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key === 'Enter') {
      event.preventDefault();
      // Re-enable auto-updates when Enter is pressed
      setShouldAutoUpdate(true); 
      handleFetchAndPlot(true);
    }
  };

  const propertySelectKey = availablePropertiesForSelection.map(p => p.displayName).join(',') || 'no-props'; // Changed to use displayName
  const constantSelectKey = availableConstantsForSelection.map(c => c.jsonKey).join(',') || 'no-const-props';

  const handleCalculateManualPoint = () => {
    if (plotMode !== 'tempDependent' || !selectedPropertyKey) {
      setManualTempPointsData([]);
      return;
    }
    
    const propDef = propertiesToPlotConfig.find(p => p.displayName === selectedPropertyKey);
    if (!propDef) {
      setManualTempPointsData([]);
      return;
    }
    const isBoilingPointSelected = propDef.displayName === "Boiling Point";
    const inputValue = parseFloat(manualTempInput);

    if (isNaN(inputValue)) {
      setManualTempPointsData(compounds.filter(c => c.data).map(c => ({
        compoundName: c.data!.name, seriesName: c.data!.name, tempCelsius: NaN,
        tempKelvin: NaN, value: null,
        unit: isBoilingPointSelected ? selectedUnit : selectedUnit,
        isValid: false, message: "Invalid input."
      })));
      return;
    }

    const unitDefToUse = propDef.availableUnits?.find(u => u.unit === selectedUnit) ||
                         propDef.availableUnits?.[0] ||
                         { unit: propDef.targetUnitName, conversionFactorFromBase: 1 };
    const currentDisplayUnit = unitDefToUse.unit;

    const newPointsData = compounds.filter(c => c.data).map(compoundState => {
  const currentCompoundData = compoundState.data!;
  const propData = currentCompoundData.properties[propDef.jsonKey];
      
      let point: ManualTempPoint = {
        compoundName: currentCompoundData.name, seriesName: currentCompoundData.name,
        tempCelsius: inputValue, tempKelvin: isBoilingPointSelected ? NaN : inputValue + 273.15,
        value: null, unit: currentDisplayUnit, isValid: false,
      };

      if (!propData) {
        point.message = `Property data not found for ${currentCompoundData.name}.`;
        return point;
      }

      const rawCoeffs = propDef.coeffs.map(cKey => ({ key: cKey, value: parseCoefficient(propData[cKey]) }));
      const passedCoeffs = rawCoeffs.map(c => c.value);
      const mw_kg_kmol_series: number | null = currentCompoundData.molarWeight ?? null;
      
      if (isBoilingPointSelected) {
        if (passedCoeffs.length >= 5 && passedCoeffs.every(c => c !== null)) {
          const [coeffA, coeffB, coeffC, coeffD, coeffE] = passedCoeffs as [number, number, number, number, number | undefined];
          const pressUnitDef = propDef.availableUnits?.find(u => u.unit === currentDisplayUnit);
          let pressureInPa = inputValue;
          if (pressUnitDef && typeof pressUnitDef.conversionFactorFromBase === 'number' && pressUnitDef.conversionFactorFromBase !== 0) {
            pressureInPa = inputValue / pressUnitDef.conversionFactorFromBase;
          } else if (pressUnitDef && currentDisplayUnit === "Pa") {
            pressureInPa = inputValue;
          } else {
            point.message = `Pressure unit conversion error for ${currentDisplayUnit}.`;
            return point;
          }

          if (pressureInPa <= 0) {
            point.message = "Pressure must be positive."; return point;
          }
          
          let initialTempGuessBp = 373.15;
          if (coeffA && coeffB && pressureInPa > 0) {
              const lnPVal = Math.log(pressureInPa);
              if (coeffA - lnPVal !== 0) {
                  const T_antoine_approx_bp = coeffB / (coeffA - lnPVal);
                  if (T_antoine_approx_bp > 0 && T_antoine_approx_bp < 1000) {
                      initialTempGuessBp = T_antoine_approx_bp;
                  }
              }
          }

          const calculatedTempKelvin = calculateBoilingPointEq101(pressureInPa, coeffA, coeffB, coeffC, coeffD, coeffE, initialTempGuessBp);
          if (calculatedTempKelvin !== null && isFinite(calculatedTempKelvin) && calculatedTempKelvin > 0) {
            point.value = calculatedTempKelvin - 273.15;
            point.isValid = true;
          } else {
            point.message = "Boiling temperature calculation failed or out of range.";
          }
        } else {
          point.message = "Insufficient coefficients for boiling point calculation.";
        }
      } else {
        const tempKelvinForCalc = inputValue + 273.15;
  const Tmin = parseCoefficient(propData.Tmin?.value ?? propData.Tmin);
  const Tmax = parseCoefficient(propData.Tmax?.value ?? propData.Tmax);

        if (Tmin === null || Tmax === null || tempKelvinForCalc < Tmin - 1e-6 || tempKelvinForCalc > Tmax + 1e-6) {
          // Use Supabase-derived temperature limits
          const TminC = Tmin !== null ? Tmin - 273.15 : NaN;
          const TmaxC = Tmax !== null ? Tmax - 273.15 : NaN;
          point.message = `Temp out of range (${formatWithSigFigs(TminC, 3)} to ${formatWithSigFigs(TmaxC, 3)} °C).`;
          return point;
        }
        
  const eqnoParsed = parseCoefficient(propData?.eqno);
  const eqnoStr = eqnoParsed !== null ? String(eqnoParsed) : (propData?.eqno ? String(propData.eqno) : '');
  if (!eqnoStr) {
          point.message = `Equation number (eqno) not found.`;
          return point;
        }

        const coeffsForCalc = { A: parseCoefficient(propData.A), B: parseCoefficient(propData.B), C: parseCoefficient(propData.C), D: parseCoefficient(propData.D), E: parseCoefficient(propData.E), };
        const Tc_K_for_calc = currentCompoundData.criticalTemp ?? undefined;

        if (propDef.requiresTc && Tc_K_for_calc === undefined) {
          point.message = `Critical Temp required but not found.`; return point;
        }

  let rawValue = calculatePropertyByEquation(eqnoStr, tempKelvinForCalc, coeffsForCalc, Tc_K_for_calc);

        if (rawValue !== null && isFinite(rawValue)) {
          let valueInBaseUnit = rawValue;
          if (propDef.conversionFactor && propDef.conversionFactor !== 1) {
            valueInBaseUnit *= propDef.conversionFactor;
          }
          let finalDisplayValue = valueInBaseUnit;
          const cFactorBase = unitDefToUse.conversionFactorFromBase;
          if (typeof cFactorBase === 'number') {
            if (cFactorBase !== 1) finalDisplayValue *= cFactorBase;
          } else if (typeof cFactorBase === 'object' && cFactorBase !== null) {
            if (mw_kg_kmol_series === null && (cFactorBase.operation === 'divide_by_mw' || cFactorBase.operation === 'multiply_by_mw')) {
              point.message = `Molar mass required for unit conversion to ${currentDisplayUnit}, but not available.`; return point;
            }
            if (mw_kg_kmol_series !== null) {
              if (cFactorBase.operation === 'divide_by_mw') {
                if (mw_kg_kmol_series === 0) { point.message = `Molar mass is zero, cannot divide for unit ${currentDisplayUnit}.`; return point; }
                finalDisplayValue /= mw_kg_kmol_series;
              } else if (cFactorBase.operation === 'multiply_by_mw') {
                finalDisplayValue *= mw_kg_kmol_series;
              }
            }
            if (typeof cFactorBase.factor === 'number' && cFactorBase.factor !== 1) {
              finalDisplayValue *= cFactorBase.factor;
            }
          }
          point.value = finalDisplayValue;
          point.isValid = true;
        } else {
          point.message = "Calculation failed or resulted in non-finite number.";
        }
      }
      return point;
    });
    setManualTempPointsData(newPointsData);
  };

  // Effect to clear manual points if selected property changes
  useEffect(() => {
    setManualTempInput('');
    setManualTempPointsData([]);
  }, [selectedPropertyKey, selectedUnit]);

  const handleManualTempInputKeyDown = (event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key === 'Enter') {
      event.preventDefault();
      if (manualTempInput.trim() && !isNaN(parseFloat(manualTempInput))) {
        handleCalculateManualPoint(); // Call renamed function
      }
    }
  };


  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                   <div className="lg:col-span-1 space-y-6">
            <Card>
              <CardContent className="space-y-4 pt-6">
                <div className="flex items-center space-x-2 mb-4">
                  <Tabs value={plotMode} onValueChange={(value) => {
                    setPlotMode(value as 'tempDependent' | 'constants');
                    // Enable auto-updates when switching plot modes
                    setShouldAutoUpdate(true);
                  }} className="w-full">
                    <TabsList className="grid w-full grid-cols-2">
                      <TabsTrigger value="tempDependent">Temp. Dependent</TabsTrigger>
                      <TabsTrigger value="constants">Constants</TabsTrigger>
                    </TabsList>
                  </Tabs>
                </div>

                {compounds.map((compound, index) => (
                    <div key={compound.id} className="space-y-1">
                        <div className="flex items-center gap-2"> 
                            <Label htmlFor={compound.id} className="text-sm whitespace-nowrap w-28"> 
                                Compound {index + 1}:
                            </Label>
                            <div className="relative flex items-center gap-2 flex-grow"> 
                                <Input
                                    ref={compound.inputRef}
                                    id={compound.id}
                                    value={compound.name}
                                    onChange={(e) => handleCompoundNameChange(compound.id, e.target.value)}
                                    onFocus={() => { setCompounds(prev => prev.map(c => c.id === compound.id ? { ...c, showSuggestions: c.suggestions.length > 0 } : c));}}
                                    onKeyDown={handleInputKeyDown} 
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
                                    <div ref={compound.suggestionsRef} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1 top-full">
                                    {compound.suggestions.map((s, i) => <div key={i} onClick={() => handleSuggestionClick(compound.id, s)} className="px-3 py-2 hover:bg-accent cursor-pointer text-sm">{s}</div>)}
                                    </div>
                                )}
                            </div>
                        </div>
                        {compound.error && <p className="text-xs text-red-500 mt-1 pl-32">{compound.error}</p>} 
                   
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
                              setSelectedPropertyKey(value); // value is now displayName
                              const newPropDef = availablePropertiesForSelection.find(p => p.displayName === value); // Find by displayName
                              const newUnit = newPropDef?.availableUnits?.[0]?.unit || newPropDef?.targetUnitName || '';
                              setSelectedUnit(newUnit);
                              // Enable auto-updates when changing properties
                              setShouldAutoUpdate(true);
                          }}
                          disabled={availablePropertiesForSelection.length === 0}
                        >
                          <SelectTrigger id="propertySelect" className="w-full">
                            <SelectValue placeholder="Select a property" />
                          </SelectTrigger>
                          <SelectContent>
                            {availablePropertiesForSelection.length > 0 ? availablePropertiesForSelection.map(prop => (
                              <SelectItem key={prop.displayName} value={prop.displayName}> {/* Changed key and value to prop.displayName */}
                                {renderPropertyName(prop.displayName, prop.symbol)}
                              </SelectItem>
                            )) : <SelectItem value="no-props" disabled>No common properties</SelectItem>}
                          </SelectContent>
                        </Select>
                      </div>
                      
                      {(() => {
                        const currentPropDef = propertiesToPlotConfig.find(p => p.displayName === selectedPropertyKey); // Changed to find by displayName
                        const unitsForCurrentProp = currentPropDef?.availableUnits;
                        if (unitsForCurrentProp && unitsForCurrentProp.length > 1) {
                          return (
                            <div className="flex-grow space-y-2 max-w-[224px]"> 
                                                              <Select
                                  value={selectedUnit}
                                  onValueChange={(value) => {
                                    setSelectedUnit(value);
                                    // Enable auto-updates when changing units
                                    setShouldAutoUpdate(true);
                                  }}
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
                  ) : ( 
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
                              // Enable auto-updates when changing constants
                              setShouldAutoUpdate(true);
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
                        if (unitsForCurrentConst && (unitsForCurrentConst.length > 1 || (unitsForCurrentConst.length === 1 && unitsForCurrentConst[0].unit !== "-"))) {
                          return (
                            <div className="flex-grow space-y-2 max-w-[224px]"> 
                                                              <Select
                                  value={selectedConstantUnit}
                                  onValueChange={(value) => {
                                    setSelectedConstantUnit(value);
                                    // Enable auto-updates when changing constant units
                                    setShouldAutoUpdate(true);
                                  }}
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

                {/* Second Virial Coefficient warning removed per user request */}

                <Button onClick={handleFetchButtonClick} disabled={loading} className="w-full">
                  {loading ? 'Fetching...' : 'Fetch & Plot Properties'}
                </Button>
                {canExportCSV && !loading && (
                    <Button onClick={handleExportCSV} variant="outline" className="w-full mt-2">
                        Export Data as CSV
                    </Button>
                )}

              </CardContent>
            </Card>

            {plotMode === 'tempDependent' && compounds.some(c => c.data) && selectedPropertyKey && (
              <Card>
                <CardHeader>
                  <CardTitle className="text-base">
                    {propertiesToPlotConfig.find(p => p.displayName === selectedPropertyKey)?.isValueAtPressure
                      ? "Get Boiling Temp at Specific Pressure"
                      : "Get Value at Specific Temperature"}
                  </CardTitle>
                </CardHeader>
                <CardContent className="space-y-3">
                  <div className="flex items-center gap-2">
                    <Input
                      type="number"
                      value={manualTempInput}
                      onChange={(e) => {
                        const raw = e.target.value;
                        const num = parseFloat(raw);
                        const CLAMP_MAX = 2000;
                        if (raw === '' || isNaN(num)) {
                          setManualTempInput(raw);
                        } else {
                          const clamped = Math.min(Math.abs(num), CLAMP_MAX) * Math.sign(num);
                          setManualTempInput(String(clamped));
                        }
                        if (raw === '') {
                          setManualTempPointsData([]);
                        }
                      }}
                      onKeyDown={handleManualTempInputKeyDown}
                      placeholder={
                        propertiesToPlotConfig.find(p => p.displayName === selectedPropertyKey)?.isValueAtPressure
                          ? "Enter Pressure"
                          : "Enter Temp"
                      }
                      className="flex-grow"
                    />
                    <span className="text-sm text-muted-foreground">
                      {propertiesToPlotConfig.find(p => p.displayName === selectedPropertyKey)?.isValueAtPressure
                        ? selectedUnit // This is the selected pressure unit
                        : "°C"}
                    </span>
                    <Button 
                      onClick={handleCalculateManualPoint} // Call renamed function
                      disabled={!manualTempInput.trim() || isNaN(parseFloat(manualTempInput))}
                    >
                      Calculate
                    </Button>
                  </div>
                  {manualTempPointsData.length > 0 && (
                    <div className="space-y-1 text-sm mt-2">
                      {manualTempPointsData.map((point, idx) => {
                        // Find the compound index to get the matching color
                        const compoundIndex = compounds.findIndex(c => c.data?.name === point.compoundName);
                        const compoundColor = compoundIndex >= 0 ? compoundColors[compoundIndex % compoundColors.length] : '#666';
                        
                        return (
                          <div key={idx} className={`pl-2 ${point.isValid ? '' : 'text-red-500'}`} style={{ color: point.isValid ? compoundColor : undefined }}>
                            <strong>{point.compoundName}:</strong>{' '}
                            {point.isValid && point.value !== null
                              ? (propertiesToPlotConfig.find(p => p.displayName === selectedPropertyKey)?.isValueAtPressure
                                  ? `${formatNumberToPrecision(point.value, 2)} °C (at ${formatNumberToPrecision(point.tempCelsius, 3)} ${point.unit})` // value is Temp C, tempCelsius is Pressure, unit is Pressure unit
                                  : `${formatNumberToPrecision(point.value, 3)} ${point.unit} (at ${formatNumberToPrecision(point.tempCelsius, 2)} °C)`) // value is Prop Val, unit is Prop unit, tempCelsius is Temp C
                              : (point.message || "N/A")}
                          </div>
                        );
                      })}
                    </div>
                  )}
                </CardContent>
              </Card>
            )}
          </div>

          <div className="lg:col-span-2 space-y-6">
            <Card>
              <CardContent className="py-2">
                <div className="relative aspect-square rounded-md">
                  {loading && (<div className="absolute inset-0 flex items-center justify-center text-muted-foreground"><div className="text-center"><div className="mb-2">Loading Properties...</div></div></div>)}
                  {!loading && compounds.every(c => !c.data && !c.error) && !overallError && (<div className="absolute inset-0 flex items-center justify-center text-muted-foreground">Please enter compound(s) and fetch properties.</div>)}
                  {overallError && !loading && (<div className="absolute inset-0 flex items-center justify-center text-red-400 px-4 text-center">Error: {overallError}</div>)}
                  {!loading && Object.keys(echartsOptions).length > 0 && compounds.some(c=>c.data) && !overallError && (
                    <ReactECharts ref={echartsRef} echarts={echarts} option={echartsOptions} style={{ height: '100%', width: '100%' }} notMerge={true} lazyUpdate={true} />
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