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
  MarkLineComponent,
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';
import type { EChartsOption } from 'echarts';

echarts.use([
  TitleComponent, TooltipComponent, GridComponent, LegendComponent, ToolboxComponent, MarkPointComponent, MarkLineComponent, LineChart, ScatterChart, CanvasRenderer
]);

import { Card, CardContent } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Label } from "@/components/ui/label";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Checkbox } from "@/components/ui/checkbox";
import { TooltipProvider } from "@/components/ui/tooltip";

import {
    calculatePropertyByEquation,
    parseCoefficient,
} from '@/lib/property-equations';

// --- Annotation Placement Helpers ---
interface NormalizationBox { width: number; height: number; xMin: number; yMin: number; }

const distanceSqNormalized = (p1: number[], p2: number[], norm: NormalizationBox): number => {
  if (norm.width === 0 || norm.height === 0) { return (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2; }
  const dx = (p1[0] - p2[0]) / norm.width;
  const dy = (p1[1] - p2[1]) / norm.height;
  return dx**2 + dy**2;
};

const distToSegmentSqNormalized = (p: number[], v: number[], w: number[], norm: NormalizationBox): number => {
  const p_norm_x = (p[0] - norm.xMin) / norm.width;
  const p_norm_y = (p[1] - norm.yMin) / norm.height;
  const v_norm_x = (v[0] - norm.xMin) / norm.width;
  const v_norm_y = (v[1] - norm.yMin) / norm.height;
  const w_norm_x = (w[0] - norm.xMin) / norm.width;
  const w_norm_y = (w[1] - norm.yMin) / norm.height;

  const dx = w_norm_x - v_norm_x;
  const dy = w_norm_y - v_norm_y;
  const l2 = dx*dx + dy*dy;

  if (l2 === 0) return distanceSqNormalized(p, v, norm);
  
  let t = ((p_norm_x - v_norm_x) * dx + (p_norm_y - v_norm_y) * dy) / l2;
  t = Math.max(0, Math.min(1, t));
  
  const closestPoint = [v[0] + t * (w[0] - v[0]), v[1] + t * (w[1] - v[1])];
  return distanceSqNormalized(p, closestPoint, norm);
};

const minDistToCurveSqNormalized = (point: number[], curveData: number[][], norm: NormalizationBox): number => {
  if (curveData.length < 2) return Infinity;
  let minDistance = Infinity;
  for (let i = 0; i < curveData.length - 1; i++) {
    const dist = distToSegmentSqNormalized(point, curveData[i], curveData[i + 1], norm);
    if (dist < minDistance) minDistance = dist;
  }
  return minDistance;
};

const getBoundaryValue = (coord: number, curveData: number[][], solveForY: boolean): number | null => {
  if (curveData.length < 2) return null;
  const xIndex = solveForY ? 0 : 1;
  const yIndex = solveForY ? 1 : 0;
  for (let i = 0; i < curveData.length - 1; i++) {
    const p1 = curveData[i];
    const p2 = curveData[i + 1];
    if ((p1[xIndex] <= coord && p2[xIndex] >= coord) || (p1[xIndex] >= coord && p2[xIndex] <= coord)) {
      if (Math.abs(p1[xIndex] - p2[xIndex]) < 1e-9) return p1[yIndex];
      const fraction = (coord - p1[xIndex]) / (p2[xIndex] - p1[xIndex]);
      return p1[yIndex] + fraction * (p2[yIndex] - p1[yIndex]);
    }
  }
  return null;
};

const getNiceAxisBounds = (dataMin: number, dataMax: number, tickCount: number = 5) => {
    const range = dataMax - dataMin;
    if (range === 0) {
      return { min: dataMin - 1, max: dataMax + 1 };
    }

    const rawStep = range / tickCount;
    const exponent = Math.floor(Math.log10(rawStep));
    const powerOf10 = 10 ** exponent;
    const fraction = rawStep / powerOf10;

    let niceFraction;
    if (fraction <= 1) niceFraction = 1;
    else if (fraction <= 2) niceFraction = 2;
    else if (fraction <= 5) niceFraction = 5;
    else niceFraction = 10;
    
    const niceStep = niceFraction * powerOf10;

    const min = Math.floor(dataMin / niceStep) * niceStep;
    const max = Math.ceil(dataMax / niceStep) * niceStep;

    return { min, max };
};

const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL || '';
const supabaseAnonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY || '';
let supabase: SupabaseClient;
if (supabaseUrl && supabaseAnonKey) {
    try { supabase = createClient(supabaseUrl, supabaseAnonKey); } catch (error) { console.error("Error initializing Supabase client:", error); }
} else { console.error("Supabase URL or Anon Key is missing."); }

const formatNumber = (num: any, precision: number = 4): string => {
    if (typeof num !== 'number' || isNaN(num)) return String(num);
    return num.toPrecision(precision);
};

const convertTempFromK = (valueK: number, unit: 'C' | 'K' | 'F'): number => {
    if (unit === 'C') return valueK - 273.15;
    if (unit === 'F') return (valueK - 273.15) * 9/5 + 32;
    return valueK;
};

const getTempUnitSymbol = (unit: 'C' | 'K' | 'F'): string => {
    if (unit === 'C') return '°C';
    if (unit === 'F') return '°F';
    return 'K';
};

interface FetchedPhaseData { name: string; criticalTemp: number | null; criticalPressure: number | null; triplePointTemp: number | null; triplePointPressure: number | null; normalBoilingPoint: number | null; normalMeltingPoint: number | null; molarWeight: number | null; heatOfFusion: number | null; heatOfVaporization: any | null; vapourPressure: any | null; liquidDensity: any | null; solidDensity: any | null; }

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
  const [tempUnit, setTempUnit] = useState<'C' | 'K' | 'F'>('C');
  const [pressureUnit, setPressureUnit] = useState<'Pa' | 'kPa' | 'bar' | 'atm'>('bar');
  const [logScaleY, setLogScaleY] = useState<boolean>(false);

  const getUnitFactor = useCallback((unit: typeof pressureUnit) => {
    switch (unit) { case 'Pa': return 1; case 'kPa': return 1 / 1000; case 'bar': return 1 / 100000; case 'atm': return 1 / 101325; default: return 1 / 100000; }
  }, []);

  const fetchCompoundData = useCallback(async (name: string): Promise<FetchedPhaseData | null> => {
    if (!supabase) throw new Error("Supabase client not initialized.");
    const trimmedName = name.trim();
    if (!trimmedName) throw new Error("Compound name cannot be empty.");
    const cacheKey = trimmedName.toLowerCase();
    if (dataCache.current.has(cacheKey)) { const cached = dataCache.current.get(cacheKey); if (cached === null) throw new Error(`Compound '${trimmedName}' not found.`); return cached ?? null; }
    setLoading(true); setError(null);
    try {
        const { data: compoundDbData, error: compoundError } = await supabase.from('compounds').select('id, name, molecular_weight').ilike('name', trimmedName).limit(1).single();
        if (compoundError || !compoundDbData) { throw new Error(compoundError?.message || `Compound '${trimmedName}' not found.`); }
        const { data: propsData, error: propsError } = await supabase.from('compound_properties').select('properties').eq('compound_id', compoundDbData.id).limit(1).single();
        if (propsError || !propsData) { throw new Error(`No properties found for ${compoundDbData.name}.`); }
        const heatOfFusionKmol = parseCoefficient(propsData.properties["Heat of fusion at melting point"]);
        const fetchedData: FetchedPhaseData = { name: compoundDbData.name, criticalTemp: parseCoefficient(propsData.properties["Critical temperature"]), criticalPressure: parseCoefficient(propsData.properties["Critical pressure"]), triplePointTemp: parseCoefficient(propsData.properties["Triple point temperature"]), triplePointPressure: parseCoefficient(propsData.properties["Triple point pressure"]), normalBoilingPoint: parseCoefficient(propsData.properties["Normal boiling point"]), normalMeltingPoint: parseCoefficient(propsData.properties["Melting point"]), molarWeight: compoundDbData.molecular_weight, heatOfFusion: heatOfFusionKmol ? heatOfFusionKmol / 1000 : null, heatOfVaporization: propsData.properties["Heat of vaporization"], vapourPressure: propsData.properties["Vapour pressure"], liquidDensity: propsData.properties["Liquid density"], solidDensity: propsData.properties["Solid density"], };
        dataCache.current.set(cacheKey, fetchedData); return fetchedData;
    } catch (err: any) { dataCache.current.set(cacheKey, null); throw err; } finally { setLoading(false); }
  }, []);

  const handleFetchAndPlot = useCallback(async () => {
      try { const data = await fetchCompoundData(compoundName); setCompoundData(data); } catch (err: any) { setError(err.message); setCompoundData(null); setEchartsOptions({}); }
  }, [compoundName, fetchCompoundData]);

  useEffect(() => { if (!hasInitialLoaded.current) { hasInitialLoaded.current = true; handleFetchAndPlot(); } }, [handleFetchAndPlot]);

  const plotPhaseDiagram = useCallback((data: FetchedPhaseData | null) => {
    if (!data || !data.criticalTemp || !data.criticalPressure || !data.triplePointTemp || !data.triplePointPressure) { setEchartsOptions({}); return; }
    
    const { criticalTemp, criticalPressure, triplePointTemp, triplePointPressure, vapourPressure, heatOfFusion, heatOfVaporization } = data;
    const isDark = resolvedTheme === 'dark';
    const textColor = isDark ? 'white' : '#000000';
    const factor = getUnitFactor(pressureUnit);
    const series: (LineSeriesOption | ScatterSeriesOption | any)[] = [];
    let vaporizationData: [number, number][] = [], fusionData: [number, number][] = [], sublimationData: [number, number][] = [];

    const transformPressure = (p: number | null): number => {
        if (p === null) return -Infinity;
        const value = p * factor;
        return logScaleY ? (value > 0 ? Math.log10(value) : -Infinity) : value;
    };

    // 0. Define Axis Ranges Manually using the helper
    const T_c_unit = convertTempFromK(criticalTemp, tempUnit);
    const T_t_unit = convertTempFromK(triplePointTemp, tempUnit);
    const { min: axisXMin, max: axisXMax } = getNiceAxisBounds(T_t_unit - 200, T_c_unit + 150, 8);
    let axisYMin: number, axisYMax: number;

    if (logScaleY) {
        axisYMin = Math.floor(transformPressure(triplePointPressure / 10000));
        axisYMax = Math.ceil(transformPressure(criticalPressure * 2));
    } else {
        const { min: niceYMin, max: niceYMax } = getNiceAxisBounds(0, transformPressure(criticalPressure * 1.5));
        axisYMin = niceYMin;
        axisYMax = niceYMax;
    }
    
    // 1. Generate Curve Data
    if (vapourPressure) {
        const vpCoeffs = { A: parseCoefficient(vapourPressure.A), B: parseCoefficient(vapourPressure.B), C: parseCoefficient(vapourPressure.C), D: parseCoefficient(vapourPressure.D), E: parseCoefficient(vapourPressure.E) };
        if (vapourPressure['eqno'] && Object.values(vpCoeffs).every(c => c !== null)) {
            const numPoints = Math.max(Math.floor(Math.abs(criticalTemp - triplePointTemp)), 2) || 100;
            for (let i = 0; i <= numPoints; i++) {
                const tempK = triplePointTemp + (criticalTemp - triplePointTemp) * (i / numPoints);
                const pPa = calculatePropertyByEquation(vapourPressure['eqno'], tempK, vpCoeffs, criticalTemp);
                if (pPa !== null && pPa > 0) vaporizationData.push([convertTempFromK(tempK, tempUnit), transformPressure(pPa)]);
            }
        }
    }
    // Corrected Fusion Line Logic
    if (triplePointTemp && triplePointPressure && criticalPressure && criticalTemp) {
        const startTempK = triplePointTemp, startPressPa = triplePointPressure;
        let highPressurePa;
        if (logScaleY) {
            highPressurePa = Math.pow(10, axisYMax) / factor;
        } else {
            highPressurePa = axisYMax / factor;
        }
        let endTempK = data.name.toLowerCase() === 'water' ? startTempK - (criticalTemp - startTempK) * 0.1 : startTempK + (criticalTemp - startTempK) * 0.2;
        for (let i = 0; i <= 50; i++) {
            const fraction = i / 50;
            fusionData.push([convertTempFromK(startTempK + (endTempK - startTempK) * fraction, tempUnit), transformPressure(startPressPa + (highPressurePa - startPressPa) * fraction)]);
        }
    }
    if (heatOfFusion && heatOfVaporization) {
        const hvCoeffs = { A: parseCoefficient(heatOfVaporization.A), B: parseCoefficient(heatOfVaporization.B), C: parseCoefficient(heatOfVaporization.C), D: parseCoefficient(heatOfVaporization.D), E: parseCoefficient(heatOfVaporization.E) };
        if (heatOfVaporization['eqno'] && Object.values(hvCoeffs).every(c => c !== null)) {
            const hVapKm = calculatePropertyByEquation(heatOfVaporization['eqno'], triplePointTemp, hvCoeffs, criticalTemp);
            if (hVapKm !== null && heatOfFusion !== null) {
                const deltaH_sub = heatOfFusion + (hVapKm / 1000);
                const startTempK = convertTempFromK(axisXMin, tempUnit) + 273.15;
                for (let i = 0; i <= 100; i++) {
                    const T2 = startTempK + (triplePointTemp - startTempK) * i / 100;
                    if (T2 > 0) {
                        const pPa = triplePointPressure * Math.exp((-deltaH_sub / 8.314) * (1 / T2 - 1 / triplePointTemp));
                        if (pPa > 0) sublimationData.push([convertTempFromK(T2, tempUnit), transformPressure(pPa)]);
                    }
                }
            }
        }
    }
    if (vaporizationData.length > 0) series.push({ name: 'Vaporization', type: 'line', data: vaporizationData, symbol: 'none', lineStyle: { width: 2.5 }, showSymbol: false, z: 2 });
    if (fusionData.length > 0) series.push({ name: 'Fusion', type: 'line', data: fusionData, symbol: 'none', lineStyle: { width: 2.5 }, showSymbol: false, z: 2 });
    if (sublimationData.length > 0) series.push({ name: 'Sublimation', type: 'line', data: sublimationData, symbol: 'none', lineStyle: { width: 2.5 }, showSymbol: false, z: 2 });
    
    // 2. Add Key Points and Supercritical Lines
    const T_c = convertTempFromK(criticalTemp, tempUnit), P_c = transformPressure(criticalPressure);
    const T_t = convertTempFromK(triplePointTemp, tempUnit), P_t = transformPressure(triplePointPressure);
    series.push({ name: 'Critical Point', type: 'scatter', data: [{ value: [T_c, P_c], itemStyle: { color: '#d62728' } }], symbolSize: 10, z: 10 });
    series.push({ name: 'Supercritical Boundary', type: 'line',
        markLine: { silent: true, symbol: 'none', label: { show: false }, lineStyle: { type: 'dashed', color: isDark ? 'rgba(255,255,255,0.4)' : 'rgba(0,0,0,0.4)', width: 1.5 }, data: [[{ coord: [T_c, P_c] }, { yAxis: axisYMax, xAxis: T_c }], [{ coord: [T_c, P_c] }, { xAxis: axisXMax, yAxis: P_c }]], z: 1 },
    });
    series.push({ name: 'Triple Point', type: 'scatter', data: [{ value: [T_t, P_t], itemStyle: { color: '#ff7f0e' } }], symbolSize: 10, z: 10 });

    // 3. Annotation Placement
    const debugSeries: any[] = [];
    const annotationPoints = [];
    const chartBox = { xMin: axisXMin, xMax: axisXMax, yMin: axisYMin, yMax: axisYMax };
    const normBox: NormalizationBox = { ...chartBox, width: chartBox.xMax - chartBox.xMin, height: chartBox.yMax - chartBox.yMin };
            
    const findBestSpot = (regionName: string, searchBox: any, boundaries: any[], isPointInRegion: (x: number, y: number) => boolean, breakTiesByCentroid: boolean) => {
        let bestPoints: number[][] = []; let maxMinDistSq = -1;
        let dbgGrid: number[][] = [], dbgValid: number[][] = [];
        const gridDensity = 100; 
        for (let i = 0; i <= gridDensity; i++) {
            const x = searchBox.xMin + (searchBox.xMax - searchBox.xMin) * i / gridDensity;
            for (let j = 0; j <= gridDensity; j++) {
                const y = searchBox.yMin + (searchBox.yMax - searchBox.yMin) * j / gridDensity;
                if (!isPointInRegion(x, y)) continue;
                const normDistToLeft = ((x - searchBox.xMin) / normBox.width)**2; const normDistToRight = ((searchBox.xMax - x) / normBox.width)**2; const normDistToBottom = ((y - searchBox.yMin) / normBox.height)**2; const normDistToTop = ((searchBox.yMax - y) / normBox.height)**2;
                let minBoundaryDistSq = Math.min(normDistToLeft, normDistToRight, normDistToBottom, normDistToTop);
                for (const curve of boundaries) { minBoundaryDistSq = Math.min(minBoundaryDistSq, minDistToCurveSqNormalized([x,y], curve, normBox)); }
                if (minBoundaryDistSq > maxMinDistSq) { maxMinDistSq = minBoundaryDistSq; bestPoints = [[x, y]]; } else if (Math.abs(minBoundaryDistSq - maxMinDistSq) < 1e-9) { bestPoints.push([x, y]); }
            }
        }
        let finalBestPoint: number[] | null = null;
        if (bestPoints.length > 0) {
            if (breakTiesByCentroid && bestPoints.length > 1) { const sumX = bestPoints.reduce((acc, p) => acc + p[0], 0); const sumY = bestPoints.reduce((acc, p) => acc + p[1], 0); finalBestPoint = [sumX / bestPoints.length, sumY / bestPoints.length]; } 
            else { finalBestPoint = bestPoints[0]; }
        }
        return finalBestPoint;
    };
    
    const verticalCritLine = [[T_c, P_c], [T_c, axisYMax]];
    const horizontalCritLine = [[T_c, P_c], [axisXMax, P_c]];
    const spotSolid = findBestSpot("Solid", chartBox, [fusionData, sublimationData], (x,y) => { if(x > T_t && y < P_t) return false; const fusionT = getBoundaryValue(y, fusionData, false); const subT = getBoundaryValue(y, sublimationData, false); return (fusionT !== null && x < fusionT) || (subT !== null && x < subT); }, false);
    if(spotSolid) annotationPoints.push({ name: 'Solid', value: spotSolid });
    const spotLiquid = findBestSpot("Liquid", chartBox, [fusionData, vaporizationData], (x,y) => { if (x >= T_c || y >= P_c || y < P_t) return false; const fusionT = getBoundaryValue(y, fusionData, false); const vapP = getBoundaryValue(x, vaporizationData, true); return fusionT !== null && vapP !== null && x > fusionT && y > vapP; }, false);
    if(spotLiquid) annotationPoints.push({ name: 'Liquid', value: spotLiquid });
    const spotGas = findBestSpot("Gas", chartBox, [sublimationData, vaporizationData, horizontalCritLine], (x,y) => { if (x < T_t) { const subP = getBoundaryValue(x, sublimationData, true); return subP !== null && y < subP; } else if (x >= T_t && x < T_c) { const vapP = getBoundaryValue(x, vaporizationData, true); return vapP !== null && y < vapP; } else { return y < P_c; } }, false);
    if(spotGas) annotationPoints.push({ name: 'Gas', value: spotGas });
    const spotSuper = findBestSpot("Supercritical", chartBox, [verticalCritLine, horizontalCritLine], (x,y) => x > T_c && y > P_c, true);
    if(spotSuper) annotationPoints.push({ name: 'Supercritical', value: spotSuper });
    series.push({ name: 'Phases', type: 'scatter', data: annotationPoints, symbol: 'circle', symbolSize: 0, label: { show: true, position: 'inside', formatter: '{b}', color: textColor, fontSize: 16, fontWeight: 'bold', fontFamily: 'Merriweather Sans', opacity: 0.6 }, z: 5, silent: true });
    
    // 4. ECharts Options
    setEchartsOptions({
        backgroundColor: 'transparent',
        title: { text: `Phase Diagram for ${data.name}`, left: 'center', textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } },
        grid: { left: '8%', right: '5%', bottom: '10%', top: '5%', containLabel: true },
        tooltip: { trigger: 'axis', axisPointer: { type: 'cross', label: { backgroundColor: isDark ? '#1e293b' : '#ffffff', borderColor: isDark ? '#3b82f6' : '#333333', color: textColor, fontFamily: 'Merriweather Sans', formatter: (params: any) => { if (params.axisDimension === 'x') { return `${formatNumber(params.value, 3)} ${getTempUnitSymbol(tempUnit)}`; } else if (params.axisDimension === 'y') { const realP = logScaleY ? Math.pow(10, params.value) : params.value; return `${formatNumber(realP, 3)} ${pressureUnit}`; } return params.value; } } }, backgroundColor: isDark ? '#1e293b' : '#ffffff', borderColor: isDark ? '#3b82f6' : '#333333', textStyle: { color: textColor, fontFamily: 'Merriweather Sans' }, formatter: (params: any) => { if (!params || params.length === 0) return ''; const temp = params[0].axisValue; let tooltipContent = `Temperature: <b>${formatNumber(temp, 3)} ${getTempUnitSymbol(tempUnit)}</b><br/>`; params.forEach((param: any) => { if (param.seriesName && !['Phases', 'Supercritical Boundary'].includes(param.seriesName) && !param.seriesName.includes('Grid') && !param.seriesName.includes('Valid') && !param.seriesName.includes('Best') && !param.seriesName.includes('Circle')) { const chartP = param.value[1]; const realP = logScaleY ? Math.pow(10, chartP) : chartP; tooltipContent += `${param.seriesName}: <b>${formatNumber(realP, 3)} ${pressureUnit}</b><br/>`; } }); return tooltipContent.slice(0, -5); } },
        legend: { data: ['Vaporization', 'Fusion', 'Sublimation', 'Critical Point', 'Triple Point'], bottom: 5, textStyle: { color: textColor, fontFamily: 'Merriweather Sans', fontSize: 14 }, inactiveColor: '#4b5563' },
        xAxis: { type: 'value', name: `Temperature (${getTempUnitSymbol(tempUnit)})`, nameLocation: 'middle', nameGap: 30, axisLabel: { color: textColor, fontFamily: 'Merriweather Sans', fontSize: 16 }, nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' }, axisLine: { show: true, lineStyle: { color: textColor, width: 2 }, onZero: false }, axisTick: { show: true, lineStyle: { color: textColor }, length: 12 }, splitLine: { show: false }, scale: false, min: axisXMin, max: axisXMax },
        yAxis: logScaleY ? { type: 'value', name: `Pressure (${pressureUnit})`, nameLocation: 'middle', nameGap: 75, min: axisYMin, max: axisYMax, axisLabel: { color: textColor, fontFamily: 'Merriweather Sans', fontSize: 14, hideOverlap: false, showMaxLabel: true, showMinLabel: true, interval: 1, formatter: (val: number) => `1e${val.toFixed(0)}` }, nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' }, axisLine: { show: true, lineStyle: { color: textColor, width: 2 }, onZero: false }, axisTick: { show: true, lineStyle: { color: textColor }, length: 10 }, minorTick: { show: true, splitNumber: 9, length: 5 }, splitLine: { show: false }, minorSplitLine: { show: false }, scale: false, } as any : { type: 'value', name: `Pressure (${pressureUnit})`, nameLocation: 'middle', nameGap: 75, axisLabel: { color: textColor, fontFamily: 'Merriweather Sans', fontSize: 14, formatter: (val: number) => { if (Math.abs(val) >= 1e6 || (Math.abs(val) < 1e-3 && val !== 0)) { return val.toExponential(1); } return val.toString(); } }, nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' }, axisLine: { show: true, lineStyle: { color: textColor, width: 2 }, onZero: false }, axisTick: { show: true, lineStyle: { color: textColor }, length: 10 }, splitLine: { show: false }, scale: false, min: axisYMin, max: axisYMax } as any,
        series: series,
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
              <CardContent className="space-y-4 pt-2">
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
                            <SelectTrigger id="temp-unit"><SelectValue /></SelectTrigger>
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
                            <SelectTrigger id="pressure-unit"><SelectValue /></SelectTrigger>
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
                            <SelectTrigger id="y-axis-scale"><SelectValue /></SelectTrigger>
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
                <CardContent className="space-y-4 pt-2">
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

                  {/* Phase Transitions */}
                  {(() => {
                    // A substance sublimes if its triple point pressure is above 1 atm (101325 Pa).
                    const sublimesAtNormalPressure = compoundData.triplePointPressure && compoundData.triplePointPressure > 101325;
                    // Only show melting point if the substance does NOT sublime at normal pressure.
                    const showMeltingPoint = !sublimesAtNormalPressure && compoundData.normalMeltingPoint;
                    // The boiling/sublimation point is shown if the data exists.
                    const showBoilingOrSublimation = compoundData.normalBoilingPoint;
                    // Don't render the section at all if there's nothing to display.
                    if (!showMeltingPoint && !showBoilingOrSublimation) {
                      return null;
                    }
                    // Count how many are shown
                    const shownCount = [showMeltingPoint, showBoilingOrSublimation].filter(Boolean).length;
                    const gridClass = shownCount === 1 ? 'grid grid-cols-1' : 'grid grid-cols-2 gap-3';
                    return (
                      <div className="space-y-2">
                        <h4 className="font-semibold text-sm text-muted-foreground uppercase tracking-wide">Phase Transitions</h4>
                        <div className={gridClass}>
                          {/* Conditionally render Normal Melting Point */}
                          {showMeltingPoint && (
                            <div className="bg-muted rounded-lg p-3 text-center">
                              <div className="text-xs text-muted-foreground mb-1">Normal Melting Point</div>
                              <div className="font-semibold">{formatNumber(convertTempFromK(compoundData.normalMeltingPoint!, tempUnit), 3)} {getTempUnitSymbol(tempUnit)}</div>
                            </div>
                          )}
                          {/* Conditionally render Normal Boiling or Sublimation Point */}
                          {showBoilingOrSublimation && (
                            <div className="bg-muted rounded-lg p-3 text-center">
                              <div className="text-xs text-muted-foreground mb-1">
                                {sublimesAtNormalPressure ? 'Normal Sublimation Point' : 'Normal Boiling Point'}
                              </div>
                              <div className="font-semibold">{formatNumber(convertTempFromK(compoundData.normalBoilingPoint!, tempUnit), 3)} {getTempUnitSymbol(tempUnit)}</div>
                            </div>
                          )}
                        </div>
                      </div>
                    );
                  })()}
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