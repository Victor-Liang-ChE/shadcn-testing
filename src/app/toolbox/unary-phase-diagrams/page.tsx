'use client';

import React, { useState, useEffect, useCallback, useRef } from 'react';
import { createClient, SupabaseClient } from '@supabase/supabase-js';
import { useTheme } from "next-themes";

import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import { LineChart, ScatterChart } from 'echarts/charts';
import type { LineSeriesOption, ScatterSeriesOption } from 'echarts/charts';
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  ToolboxComponent,
  MarkPointComponent,
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';
import type { EChartsOption } from 'echarts';

echarts.use([
  TitleComponent, TooltipComponent, GridComponent, LegendComponent, ToolboxComponent, MarkPointComponent, LineChart, ScatterChart, CanvasRenderer
]);

import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Label } from "@/components/ui/label";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { TooltipProvider } from "@/components/ui/tooltip";


import {
    calculatePropertyByEquation,
    parseCoefficient,
} from '@/lib/property-equations';

const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL || '';
const supabaseAnonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY || '';
let supabase: SupabaseClient;
if (supabaseUrl && supabaseAnonKey) {
    try {
        supabase = createClient(supabaseUrl, supabaseAnonKey);
    } catch (error) {
        console.error("Error initializing Supabase client for Phase Diagram:", error);
    }
} else {
    console.error("Supabase URL or Anon Key is missing for Phase Diagram.");
}

const formatNumber = (num: any, precision: number = 4): string => {
    if (typeof num !== 'number' || isNaN(num)) return String(num);
    return num.toPrecision(precision);
};

// NEW: Temperature conversion functions
const convertTempFromK = (valueK: number, unit: 'C' | 'K' | 'F'): number => {
    if (unit === 'C') return valueK - 273.15;
    if (unit === 'F') return (valueK - 273.15) * 9/5 + 32;
    return valueK; // Return Kelvin
};

const getTempUnitSymbol = (unit: 'C' | 'K' | 'F'): string => {
    if (unit === 'C') return '°C';
    if (unit === 'F') return '°F';
    return 'K';
};


interface FetchedPhaseData {
    name: string;
    criticalTemp: number | null;
    criticalPressure: number | null;
    triplePointTemp: number | null;
    triplePointPressure: number | null;
    normalBoilingPoint: number | null;
    normalMeltingPoint: number | null;
    molarWeight: number | null;
    heatOfFusion: number | null;
    heatOfVaporization: any | null;
    vapourPressure: any | null;
    liquidDensity: any | null;
    solidDensity: any | null;
}

export default function PhaseDiagramPage() {
  const { resolvedTheme } = useTheme();
  const dataCache = useRef(new Map<string, FetchedPhaseData | null>());
  const hasInitialLoaded = useRef(false);
  const [compoundName, setCompoundName] = useState<string>('Water');
  const [compoundSuggestions, setCompoundSuggestions] = useState<string[]>([]);
  const [showSuggestions, setShowSuggestions] = useState(false);
  const [compoundData, setCompoundData] = useState<FetchedPhaseData | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
  const inputRef = useRef<HTMLInputElement>(null);
  const suggestionsRef = useRef<HTMLDivElement>(null);
  
  // NEW: State for temperature unit
  const [tempUnit, setTempUnit] = useState<'C' | 'K' | 'F'>('C');
  const [pressureUnit, setPressureUnit] = useState<'Pa' | 'kPa' | 'bar' | 'atm'>('bar');
  const [logScaleY, setLogScaleY] = useState<boolean>(false);

  const getUnitFactor = useCallback((unit: typeof pressureUnit) => {
    switch (unit) {
        case 'Pa': return 1;
        case 'kPa': return 1 / 1000;
        case 'bar': return 1 / 100000;
        case 'atm': return 1 / 101325;
        default: return 1 / 100000;
    }
  }, []);

  const fetchCompoundData = useCallback(async (name: string): Promise<FetchedPhaseData | null> => {
    if (!supabase) throw new Error("Supabase client not initialized.");
    const trimmedName = name.trim();
    if (!trimmedName) throw new Error("Compound name cannot be empty.");

    const cacheKey = trimmedName.toLowerCase();
    if (dataCache.current.has(cacheKey)) {
        const cached = dataCache.current.get(cacheKey);
        if (cached === null) throw new Error(`Compound '${trimmedName}' not found.`);
        return cached ?? null;
    }

    setLoading(true);
    setError(null);

    try {
        const { data: compoundDbData, error: compoundError } = await supabase
            .from('compounds').select('id, name, molecular_weight').ilike('name', trimmedName).limit(1).single();

        if (compoundError || !compoundDbData) {
            throw new Error(compoundError?.message || `Compound '${trimmedName}' not found.`);
        }
        const compoundId = compoundDbData.id;
        const foundName = compoundDbData.name;

        const { data: propsData, error: propsError } = await supabase
            .from('compound_properties').select('properties').eq('compound_id', compoundId).limit(1).single();

        if (propsError || !propsData) {
            throw new Error(`No properties found for ${foundName}.`);
        }
        const properties = propsData.properties;

        const heatOfFusionKmol = parseCoefficient(properties["Heat of fusion at melting point"]);

        const fetchedData: FetchedPhaseData = {
            name: foundName,
            criticalTemp: parseCoefficient(properties["Critical temperature"]),
            criticalPressure: parseCoefficient(properties["Critical pressure"]),
            triplePointTemp: parseCoefficient(properties["Triple point temperature"]),
            triplePointPressure: parseCoefficient(properties["Triple point pressure"]),
            normalBoilingPoint: parseCoefficient(properties["Normal boiling point"]),
            normalMeltingPoint: parseCoefficient(properties["Melting point"]),
            molarWeight: compoundDbData.molecular_weight,
            heatOfFusion: heatOfFusionKmol ? heatOfFusionKmol / 1000 : null,
            heatOfVaporization: properties["Heat of vaporization"],
            vapourPressure: properties["Vapour pressure"],
            liquidDensity: properties["Liquid density"],
            solidDensity: properties["Solid density"],
        };

        dataCache.current.set(cacheKey, fetchedData);
        return fetchedData;
    } catch (err: any) {
        dataCache.current.set(cacheKey, null);
        throw err;
    } finally {
        setLoading(false);
    }
  }, []);

  const handleFetchAndPlot = useCallback(async () => {
      try {
          const data = await fetchCompoundData(compoundName);
          setCompoundData(data);
      } catch (err: any) {
          setError(err.message);
          setCompoundData(null);
          setEchartsOptions({});
      }
  }, [compoundName, fetchCompoundData]);

  useEffect(() => {
    if (!hasInitialLoaded.current) {
      hasInitialLoaded.current = true;
      handleFetchAndPlot();
    }
  }, [handleFetchAndPlot]);

  const plotPhaseDiagram = useCallback((data: FetchedPhaseData | null) => {
    if (!data) {
        setEchartsOptions({});
        return;
    }
    
    const {
        criticalTemp, criticalPressure,
        triplePointTemp, triplePointPressure,
        vapourPressure, molarWeight,
        heatOfFusion, heatOfVaporization
    } = data;

    const factor = getUnitFactor(pressureUnit);
    const series: (LineSeriesOption | ScatterSeriesOption)[] = [];

    // NEW: bounds to build custom y-axis for log scale
    let minChartY = Infinity;
    let maxChartY = -Infinity;
    const updateBounds = (val: number) => {
        if (val < minChartY) minChartY = val;
        if (val > maxChartY) maxChartY = val;
    };

    // helper to transform pressure to chart value
    const transformPressure = (p: number) => {
        const value = p * factor;
        return logScaleY ? Math.log10(value) : value;
    };
    
    // 1. Vaporization Curve (Liquid-Gas) - Starts at Triple Point
    if (vapourPressure && triplePointTemp && criticalTemp && criticalPressure) {
        const vpCoeffs = { A: parseCoefficient(vapourPressure.A), B: parseCoefficient(vapourPressure.B), C: parseCoefficient(vapourPressure.C), D: parseCoefficient(vapourPressure.D), E: parseCoefficient(vapourPressure.E) };
        const eqno = vapourPressure['eqno']?.toString();

        if (eqno && Object.values(vpCoeffs).every(c => c !== null)) {
            const curveData = [];
            // Start from the triple point
            const p0 = transformPressure(triplePointPressure!);
            updateBounds(p0);
            curveData.push([convertTempFromK(triplePointTemp, tempUnit), p0]);
            
            // Generate points every 1 degree in the selected unit
            const tempRangeInSelectedUnit = convertTempFromK(criticalTemp, tempUnit) - convertTempFromK(triplePointTemp, tempUnit);
            const numPoints = Math.floor(Math.abs(tempRangeInSelectedUnit));
            
            for (let i = 1; i < numPoints; i++) {
                const tempInSelectedUnit = convertTempFromK(triplePointTemp, tempUnit) + i * Math.sign(tempRangeInSelectedUnit);
                // Convert back to Kelvin for calculation
                let TK = tempInSelectedUnit;
                if (tempUnit === 'C') TK = tempInSelectedUnit + 273.15;
                else if (tempUnit === 'F') TK = (tempInSelectedUnit - 32) * 5/9 + 273.15;
                
                if (TK >= triplePointTemp && TK <= criticalTemp) {
                    const pressurePa = calculatePropertyByEquation(eqno, TK, vpCoeffs, criticalTemp);
                    if (pressurePa !== null && pressurePa > 0) {
                        const pVal = transformPressure(pressurePa);
                        updateBounds(pVal);
                        curveData.push([tempInSelectedUnit, pVal]);
                    }
                }
            }
            // End at the critical point
            const pCrit = transformPressure(criticalPressure!);
            updateBounds(pCrit);
            curveData.push([convertTempFromK(criticalTemp, tempUnit), pCrit]);
            
            series.push({
                name: 'Vaporization', type: 'line', color: '#dc2626', data: curveData,
                symbol: 'none', lineStyle: { width: 2.5 }, showSymbol: false, z: 1
            });
        }
    }

    // 2. Fusion Curve (Solid-Liquid) - NEW IMPLEMENTATION
    if (triplePointTemp && triplePointPressure && criticalPressure) {
        const fusionCurveData = [];
        const startTemp = convertTempFromK(triplePointTemp, tempUnit);
        const startPress = triplePointPressure * factor;

        // Define a high pressure point to draw the line up to. Use 1.5x critical pressure as a reasonable max.
        const highPressure = criticalPressure * 1.5 * factor;
        let endTemp = startTemp;

        if (data.name.toLowerCase() === 'water') {
            // For water, slope is negative. Decrease temp slightly as pressure increases.
            endTemp = startTemp - 5; // Negative slope
        } else {
            // For most substances, slope is positive. Increase temp slightly.
            endTemp = startTemp + 20; // Steep positive slope
        }

        // Generate points every 1 degree for better tooltip coverage
        const tempRange = Math.abs(endTemp - startTemp);
        const numPoints = Math.max(Math.floor(tempRange), 2);
        
        for (let i = 0; i <= numPoints; i++) {
            const tempFraction = i / numPoints;
            const currentTemp = startTemp + (endTemp - startTemp) * tempFraction;
            const currentPress = startPress + (highPressure - startPress) * tempFraction;
            const pFusion = transformPressure(currentPress / factor); // currentPress already has factor applied, revert to Pa for transform
            updateBounds(pFusion);
            fusionCurveData.push([currentTemp, pFusion]);
        }
        
        series.push({
            name: 'Fusion', type: 'line', color: '#91CC75', data: fusionCurveData,
            symbol: 'none', lineStyle: { width: 2.5 }, showSymbol: false, z: 1
        });
    }

    // 3. Sublimation Curve (Solid-Gas)
    if (triplePointTemp && triplePointPressure && heatOfFusion && heatOfVaporization) {
        const hvCoeffs = {
            A: parseCoefficient(heatOfVaporization.A), B: parseCoefficient(heatOfVaporization.B),
            C: parseCoefficient(heatOfVaporization.C), D: parseCoefficient(heatOfVaporization.D),
            E: parseCoefficient(heatOfVaporization.E),
        };
        const hvEqno = heatOfVaporization['eqno']?.toString();
        
        let deltaH_vap_at_tp = null;
        if (hvEqno && Object.values(hvCoeffs).every(c => c !== null)) {
            const deltaH_vap_kmol = calculatePropertyByEquation(hvEqno, triplePointTemp, hvCoeffs, criticalTemp);
            // FIX: Convert the result from J/kmol to J/mol
            if (deltaH_vap_kmol !== null) {
                deltaH_vap_at_tp = deltaH_vap_kmol / 1000;
            }
        }

        // Make sure both heat of fusion and vaporization are valid before proceeding
        if (deltaH_vap_at_tp !== null && heatOfFusion !== null) {
            // Now both values are correctly in J/mol
            const deltaH_sub = heatOfFusion + deltaH_vap_at_tp;
            const R = 8.314; // Gas constant in J/(mol·K)
            
            const sublimationCurveData = [];
            const T1 = triplePointTemp;
            const P1 = triplePointPressure;
            
            const startTempK = Math.max(150, T1 - 120);
            const numPoints = 100;

            for (let i = 0; i <= numPoints; i++) {
                const T2 = startTempK + (T1 - startTempK) * i / numPoints;
                if (T2 > 0) {
                    const exponent = (-deltaH_sub / R) * (1 / T2 - 1 / T1);
                    const pressurePa = P1 * Math.exp(exponent);
                    
                    if (pressurePa > 0) {
                        const pSub = transformPressure(pressurePa);
                        updateBounds(pSub);
                        sublimationCurveData.push([convertTempFromK(T2, tempUnit), pSub]);
                    }
                }
            }

            series.push({
                name: 'Sublimation', type: 'line', color: '#1d4ed8',
                data: sublimationCurveData,
                symbol: 'none', lineStyle: { width: 2.5 }, showSymbol: false, z: 1
            });
        }
    }

    // Add key points AFTER curves so they appear on top
    if (criticalTemp && criticalPressure) {
        series.push({
            name: 'Critical Point', type: 'scatter', color: '#d62728', 
            data: (() => { const y = transformPressure(criticalPressure!); updateBounds(y); return [{ name: 'Critical Point', value: [convertTempFromK(criticalTemp, tempUnit), y], itemStyle: { color: '#d62728', opacity: 1 } }]; })(),
            symbolSize: 10, z: 10, itemStyle: { opacity: 1 }
        });
    }
    if (triplePointTemp && triplePointPressure) {
        series.push({
            name: 'Triple Point', type: 'scatter', color: '#ff7f0e', 
            data: (() => { const y = transformPressure(triplePointPressure!); updateBounds(y); return [{ name: 'Triple Point', value: [convertTempFromK(triplePointTemp, tempUnit), y], itemStyle: { color: '#ff7f0e', opacity: 1 } }]; })(),
            symbolSize: 10, z: 10, itemStyle: { opacity: 1 }
        });
    }
    
    
    // ECharts options (same as before)
    const isDark = resolvedTheme === 'dark';
    const textColor = isDark ? 'white' : '#000000';
    setEchartsOptions({
        backgroundColor: 'transparent',
        title: { text: `Phase Diagram for ${data.name}`, left: 'center', textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } },
        grid: { left: '8%', right: '5%', bottom: '10%', top: '5%', containLabel: true },
        tooltip: {
            trigger: 'axis',
            axisPointer: { 
                type: 'cross',
                label: {
                    backgroundColor: isDark ? '#1e293b' : '#ffffff',
                    borderColor: isDark ? '#3b82f6' : '#333333',
                    color: textColor,
                    fontFamily: 'Merriweather Sans',
                    formatter: (params: any) => {
                        if (params.axisDimension === 'x') {
                            return `${formatNumber(params.value, 3)} ${getTempUnitSymbol(tempUnit)}`;
                        } else if (params.axisDimension === 'y') {
                            const realP = logScaleY ? Math.pow(10, params.value) : params.value;
                            return `${formatNumber(realP, 3)} ${pressureUnit}`;
                        }
                        return params.value;
                    }
                }
            },
            backgroundColor: isDark ? '#1e293b' : '#ffffff',
            borderColor: isDark ? '#3b82f6' : '#333333',
            textStyle: { color: textColor, fontFamily: 'Merriweather Sans' },
            formatter: (params: any) => {
                if (!params || params.length === 0) return '';
                const temp = params[0].axisValue;
                let tooltipContent = `Temperature: <b>${formatNumber(temp, 3)} ${getTempUnitSymbol(tempUnit)}</b><br/>`;
                
                params.forEach((param: any) => {
                    if (param.seriesName === 'Vaporization' || param.seriesName === 'Fusion' || param.seriesName === 'Sublimation') {
                        const chartP = param.value[1];
                        const realP = logScaleY ? Math.pow(10, chartP) : chartP;
                        tooltipContent += `${param.seriesName}: <b>${formatNumber(realP, 3)} ${pressureUnit}</b><br/>`;
                    } else if (param.seriesName.includes('Point')) {
                        const chartP2 = param.value[1];
                        const realP2 = logScaleY ? Math.pow(10, chartP2) : chartP2;
                        tooltipContent += `${param.seriesName}: <b>${formatNumber(realP2, 3)} ${pressureUnit}</b><br/>`;
                    }
                });
                
                return tooltipContent.slice(0, -5); // Remove last <br/>
            }
        },
        legend: { bottom: 5, textStyle: { color: textColor, fontFamily: 'Merriweather Sans', fontSize: 14 }, inactiveColor: '#4b5563' },
        xAxis: { type: 'value', name: `Temperature (${getTempUnitSymbol(tempUnit)})`, nameLocation: 'middle', nameGap: 30, axisLabel: { color: textColor, fontFamily: 'Merriweather Sans', fontSize: 16 }, nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' }, axisLine: { show: true, lineStyle: { color: textColor, width: 2 }, onZero: false }, axisTick: { show: true, lineStyle: { color: textColor }, length: 12 }, splitLine: { show: false }, scale: true },
        yAxis: logScaleY ? (() => {
            const minVal = isFinite(minChartY) ? minChartY : 1e-3;
            const maxVal = isFinite(maxChartY) ? maxChartY : 1e3;
            const minLog = Math.floor(minVal);
            const maxLog = Math.ceil(maxVal);
            const decadeCount = maxLog - minLog;
            return {
                type: 'value',
                name: `Pressure (${pressureUnit})`,
                nameLocation: 'middle',
                nameGap: 75,
                min: minLog,
                max: maxLog,
                interval: 1,
                logBase: undefined,
                axisLabel: {
                    color: textColor,
                    fontFamily: 'Merriweather Sans',
                    fontSize: 14,
                    hideOverlap: false,
                    showMaxLabel: true,
                    showMinLabel: true,
                    interval: 0,
                    formatter: (val: number) => `1e${val.toFixed(0)}`
                },
                nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
                axisLine: { show: true, lineStyle: { color: textColor, width: 2 }, onZero: false },
                axisTick: { show: true, lineStyle: { color: textColor }, length: 10 },
                minorTick: { show: true, splitNumber: 9, length: 5 },
                splitLine: { show: false },
                minorSplitLine: { show: false },
                scale: true,
            } as any;
        })() : {
            type: 'value',
            name: `Pressure (${pressureUnit})`,
            nameLocation: 'middle',
            nameGap: 75,
            axisLabel: {
                color: textColor,
                fontFamily: 'Merriweather Sans',
                fontSize: 14,
                formatter: (val: number) => {
                    if (Math.abs(val) >= 1e6 || (Math.abs(val) < 1e-3 && val !== 0)) {
                        return val.toExponential(1);
                    }
                    return val.toString();
                }
            },
            nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
            axisLine: { show: true, lineStyle: { color: textColor, width: 2 }, onZero: false },
            axisTick: { show: true, lineStyle: { color: textColor }, length: 10 },
            splitLine: { show: false },
            scale: true,
        } as any,
        series: series as any,
    });
  }, [pressureUnit, tempUnit, logScaleY, resolvedTheme, getUnitFactor]);

  useEffect(() => {
      plotPhaseDiagram(compoundData);
  }, [compoundData, plotPhaseDiagram]);

  const handleFetchSuggestions = useCallback(async (value: string) => {
    if (!value || value.length < 2 || !supabase) {
        setCompoundSuggestions([]);
        setShowSuggestions(false);
        return;
    }
    try {
        const { data, error } = await supabase.from('compounds').select('name').ilike('name', `${value}%`).limit(5);
        if (error) throw error;
        setCompoundSuggestions(data.map(d => d.name));
        setShowSuggestions(data.length > 0);
    } catch (err) { console.error("Suggestion fetch error:", err); }
  }, []);

  const handleNameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
      const value = e.target.value;
      setCompoundName(value);
      handleFetchSuggestions(value);
  };
  const handleSuggestionClick = (name: string) => {
      setCompoundName(name);
      setShowSuggestions(false);
      inputRef.current?.focus();
  };
  const handleInputKeyDown = (event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key === 'Enter') {
      event.preventDefault();
      handleFetchAndPlot();
    }
  };
  useEffect(() => {
    function handleClickOutside(event: MouseEvent) {
      if (suggestionsRef.current && !suggestionsRef.current.contains(event.target as Node) && inputRef.current && !inputRef.current.contains(event.target as Node)) {
        setShowSuggestions(false);
      }
    }
    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, []);

  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          <div className="lg:col-span-1 space-y-6">
            <Card>
              <CardContent className="space-y-4">
                <div className="space-y-1">
                    <Label htmlFor="compound-input">Compound Name</Label>
                    <div className="relative">
                        <Input id="compound-input" ref={inputRef} value={compoundName} onChange={handleNameChange} onFocus={() => handleFetchSuggestions(compoundName)} onKeyDown={handleInputKeyDown} placeholder="e.g., Water" autoComplete="off" />
                        {showSuggestions && compoundSuggestions.length > 0 && (
                            <div ref={suggestionsRef} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1 top-full">
                                {compoundSuggestions.map((s, i) => (<div key={i} onClick={() => handleSuggestionClick(s)} className="px-3 py-2 hover:bg-accent cursor-pointer text-sm">{s}</div>))}
                            </div>
                        )}
                    </div>
                </div>
                <div className="grid grid-cols-3 gap-4 items-end">
                    <div className="space-y-1">
                        <Label htmlFor="temp-unit">T. Unit</Label>
                        <Select value={tempUnit} onValueChange={(v) => setTempUnit(v as any)}>
                            <SelectTrigger id="temp-unit" className="justify-center"><SelectValue placeholder="Select unit" className="text-center" /></SelectTrigger>
                            <SelectContent>
                                <SelectItem value="C">°C</SelectItem>
                                <SelectItem value="K">K</SelectItem>
                                <SelectItem value="F">°F</SelectItem>
                            </SelectContent>
                        </Select>
                    </div>
                    <div className="space-y-1">
                        <Label htmlFor="pressure-unit">P. Unit</Label>
                        <Select value={pressureUnit} onValueChange={(v) => setPressureUnit(v as any)}>
                            <SelectTrigger id="pressure-unit" className="justify-center"><SelectValue placeholder="Select unit" className="text-center" /></SelectTrigger>
                            <SelectContent>
                                <SelectItem value="bar">bar</SelectItem>
                                <SelectItem value="Pa">Pa</SelectItem>
                                <SelectItem value="kPa">kPa</SelectItem>
                                <SelectItem value="atm">atm</SelectItem>
                            </SelectContent>
                        </Select>
                    </div>
                    <div className="space-y-1">
                        <Label htmlFor="y-axis-scale">Scale</Label>
                        <Select value={logScaleY ? 'log' : 'linear'} onValueChange={(v) => setLogScaleY(v === 'log')}>
                            <SelectTrigger id="y-axis-scale" className="justify-center"><SelectValue placeholder="Select scale" className="text-center" /></SelectTrigger>
                            <SelectContent>
                                <SelectItem value="linear">Linear</SelectItem>
                                <SelectItem value="log">Log</SelectItem>
                            </SelectContent>
                        </Select>
                    </div>
                </div>
                <Button onClick={handleFetchAndPlot} disabled={loading} className="w-full">{loading ? 'Loading...' : 'Generate Diagram'}</Button>
              </CardContent>
            </Card>
            {compoundData && !loading && (
              <Card>
                <CardContent className="space-y-4">
                  {/* Critical Properties */}
                  <div className="space-y-2">
                    <h4 className="font-semibold text-sm text-muted-foreground uppercase tracking-wide">Critical Point</h4>
                    <div className="grid grid-cols-2 gap-3">
                      {compoundData.criticalTemp && (
                        <div className="bg-muted rounded-lg p-3 text-center">
                          <div className="text-xs text-muted-foreground mb-1">Temperature</div>
                          <div className="font-semibold">{formatNumber(convertTempFromK(compoundData.criticalTemp, tempUnit), 3)} {getTempUnitSymbol(tempUnit)}</div>
                        </div>
                      )}
                      {compoundData.criticalPressure && (
                        <div className="bg-muted rounded-lg p-3 text-center">
                          <div className="text-xs text-muted-foreground mb-1">Pressure</div>
                          <div className="font-semibold">{formatNumber(compoundData.criticalPressure * getUnitFactor(pressureUnit), 3)} {pressureUnit}</div>
                        </div>
                      )}
                    </div>
                  </div>

                  {/* Triple Point Properties */}
                  {(compoundData.triplePointTemp || compoundData.triplePointPressure) && (
                    <div className="space-y-2">
                      <h4 className="font-semibold text-sm text-muted-foreground uppercase tracking-wide">Triple Point</h4>
                      <div className="grid grid-cols-2 gap-3">
                        {compoundData.triplePointTemp && (
                          <div className="bg-muted rounded-lg p-3 text-center">
                            <div className="text-xs text-muted-foreground mb-1">Temperature</div>
                            <div className="font-semibold">{formatNumber(convertTempFromK(compoundData.triplePointTemp, tempUnit), 3)} {getTempUnitSymbol(tempUnit)}</div>
                          </div>
                        )}
                        {compoundData.triplePointPressure && (
                          <div className="bg-muted rounded-lg p-3 text-center">
                            <div className="text-xs text-muted-foreground mb-1">Pressure</div>
                            <div className="font-semibold">{formatNumber(compoundData.triplePointPressure * getUnitFactor(pressureUnit), 3)} {pressureUnit}</div>
                          </div>
                        )}
                      </div>
                    </div>
                  )}

                  {/* Phase Transition Points */}
                  {(compoundData.normalBoilingPoint || compoundData.normalMeltingPoint) && (
                    <div className="space-y-2">
                      <h4 className="font-semibold text-sm text-muted-foreground uppercase tracking-wide">Phase Transitions</h4>
                      <div className="grid grid-cols-1 gap-2">
                        {compoundData.normalMeltingPoint && (
                          <div className="bg-muted rounded-lg p-3 flex justify-between items-center">
                            <div className="text-sm font-medium">Normal Melting Point</div>
                            <div className="font-semibold">{formatNumber(convertTempFromK(compoundData.normalMeltingPoint, tempUnit), 3)} {getTempUnitSymbol(tempUnit)}</div>
                          </div>
                        )}
                        {compoundData.normalBoilingPoint && (
                          <div className="bg-muted rounded-lg p-3 flex justify-between items-center">
                            <div className="text-sm font-medium">
                              {(() => {
                                const isNormalSublimationPoint = compoundData.triplePointTemp && compoundData.normalBoilingPoint < compoundData.triplePointTemp;
                                return isNormalSublimationPoint ? 'Normal Sublimation Point' : 'Normal Boiling Point';
                              })()}
                            </div>
                            <div className="font-semibold">{formatNumber(convertTempFromK(compoundData.normalBoilingPoint, tempUnit), 3)} {getTempUnitSymbol(tempUnit)}</div>
                          </div>
                        )}
                      </div>
                    </div>
                  )}
                </CardContent>
              </Card>
            )}
          </div>
          <div className="lg:col-span-2">
            <Card>
              <CardContent className="py-2">
                <div className="relative aspect-square rounded-md">
                    {loading && (<div className="absolute inset-0 flex items-center justify-center text-muted-foreground">Loading Diagram...</div>)}
                    {error && !loading && (<div className="absolute inset-0 flex items-center justify-center text-red-400 px-4 text-center">Error: {error}</div>)}
                    {!loading && !error && Object.keys(echartsOptions).length > 0 && (<ReactECharts echarts={echarts} option={echartsOptions} style={{ height: '100%', width: '100%' }} notMerge={true} lazyUpdate={true} />)}
                    {!loading && !error && Object.keys(echartsOptions).length === 0 && (<div className="absolute inset-0 flex items-center justify-center text-muted-foreground">Enter a compound to generate its phase diagram.</div>)}
                </div>
              </CardContent>
            </Card>
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}