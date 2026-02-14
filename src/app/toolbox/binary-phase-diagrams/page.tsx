'use client';

import React, { useState, useEffect, useCallback, useRef, useMemo } from 'react';
import { createClient } from '@supabase/supabase-js';
import { useTheme } from "next-themes";
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsOption, SeriesOption } from 'echarts';
import { LineChart, ScatterChart } from 'echarts/charts';
import { TitleComponent, TooltipComponent, GridComponent, LegendComponent, MarkPointComponent, ToolboxComponent } from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Tabs, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { ArrowLeftRight } from 'lucide-react';
import { TooltipProvider } from "@/components/ui/tooltip";
import { Slider } from "@/components/ui/slider";

import {
    calculatePsat_Pa as libCalculatePsat_Pa,
    fetchUnifacInteractionParams,
    calculateBubbleTemperatureUnifac,
    calculateBubblePressureUnifac,
    calculateDewTemperatureUnifac,
    calculateDewPressureUnifac,
    fetchNrtlParameters,
    calculateBubbleTemperatureNrtl,
    calculateBubblePressureNrtl,
    fetchPrInteractionParams,
    calculateBubbleTemperaturePr,
    calculateBubblePressurePr,
    calculateDewTemperaturePr,
    calculateDewPressurePr,
    fetchSrkInteractionParams,
    calculateBubbleTemperatureSrk,
    calculateBubblePressureSrk,
    calculateDewTemperatureSrk,
    calculateDewPressureSrk,
    fetchUniquacInteractionParams,
    calculateBubbleTemperatureUniquac,
    calculateBubblePressureUniquac,
    fetchWilsonInteractionParams,
    calculateBubbleTemperatureWilson,
    calculateBubblePressureWilson,
    antoineBoilingPointSolverLocal,
} from '@/lib/vle-calculations';

// Import Shared VLE Types
import type {
    CompoundData,
    BubbleDewResult,
    AntoineParams,
    UnifacGroupComposition,
    PrPureComponentParams,
    SrkPureComponentParams,
    UniquacPureComponentParams,
    WilsonPureComponentParams,
} from '@/lib/vle-types';
import { fetchCompoundDataFromHysys, fetchCompoundSuggestions, resolveSimName, formatCompoundName } from '@/lib/antoine-utils';
import type { CompoundAlias } from '@/lib/antoine-utils';

export type FluidPackageType = 'unifac' | 'pr' | 'srk' | 'uniquac' | 'wilson' | 'nrtl';
export type DiagramType = 'txy' | 'pxy' | 'xy';

// Temporary types to satisfy linter
type UnifacParameters = any;
type NrtlInteractionParams = any;
type PrInteractionParams = any;
type SrkInteractionParams = any;
type UniquacInteractionParams = any;
type WilsonInteractionParams = any;

echarts.use([TitleComponent, TooltipComponent, GridComponent, LegendComponent, MarkPointComponent, ToolboxComponent, LineChart, ScatterChart, CanvasRenderer]);

const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL!;
const supabaseAnonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY!;
const supabase = createClient(supabaseUrl, supabaseAnonKey);

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

// Utility to compute a 'nice' slider step given the maximum value
function computeStep(maxVal: number): number {
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
}

export default function VleDiagramPage() {
    const { resolvedTheme } = useTheme();
    const [comp1Name, setComp1Name] = useState('Methanol');
    const [comp2Name, setComp2Name] = useState('H2O');
    const [diagramType, setDiagramType] = useState<DiagramType>('txy');
    const [temperatureC, setTemperatureC] = useState<number | null>(null);
    const [pressureBar, setPressureBar] = useState<number | null>(1);
    const [fluidPackage, setFluidPackage] = useState<FluidPackageType>('uniquac');
    const [temperatureInput, setTemperatureInput] = useState<string>(String(temperatureC));
    const [pressureInput, setPressureInput] = useState<string>(String(pressureBar));
    const [tempMax, setTempMax] = useState<string>('120');
    const [pressureMax, setPressureMax] = useState<string>('10');
    const [useTemperatureForXY, setUseTemperatureForXY] = useState(false);
    const [chartData, setChartData] = useState<any | null>(null); // Using 'any' for now during refactor
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
    const [displayedParams, setDisplayedParams] = useState({ comp1: '', comp2: '', temp: null as number | null, pressure: null as number | null, package: '', type: '' });
    
    // Caching states
    const [comp1Data, setComp1Data] = useState<CompoundData | null>(null);
    const [comp2Data, setComp2Data] = useState<CompoundData | null>(null);
    const [interactionParams, setInteractionParams] = useState<any>(null);

    const [comp1Suggestions, setComp1Suggestions] = useState<CompoundAlias[]>([]);
    const [comp2Suggestions, setComp2Suggestions] = useState<CompoundAlias[]>([]);
    const [showComp1Suggestions, setShowComp1Suggestions] = useState(false);
    const [showComp2Suggestions, setShowComp2Suggestions] = useState(false);
    const [activeSuggestionInput, setActiveSuggestionInput] = useState<'comp1' | 'comp2' | null>(null);
    const [autoGenerateOnCompChange, setAutoGenerateOnCompChange] = useState(true); // Auto-generate on first load

    // Display names (human-readable FullName from HYSYS ALIASES)
    const [comp1DisplayName, setComp1DisplayName] = useState('Methanol');
    const [comp2DisplayName, setComp2DisplayName] = useState('Water');

    const input1Ref = useRef<HTMLInputElement>(null);
    const input2Ref = useRef<HTMLInputElement>(null);
    const activeComponentRef = useRef<HTMLDivElement>(null);
    const echartsRef = useRef<ReactECharts | null>(null);
    const plotDataRef = useRef<{ bubble: [number, number][], dew: [number, number][] }>({ bubble: [], dew: [] });
    const leverRuleDataRef = useRef<{ liquid: number | null, vapor: number | null }>({ liquid: null, vapor: null });
    const yMinRef = useRef<number | undefined>(undefined);

    const generateEchartsOptions = useCallback((data: any, params = displayedParams, themeOverride?: string) => {
        if (!data || data.x.length === 0) return;
        // Reset plot data on each generation
        plotDataRef.current = { bubble: [], dew: [] };
        const series: SeriesOption[] = [];
        let yAxisName = "Mole Fraction", titleConditionText = "";
        const capitalize = (s: string) => s ? s.charAt(0).toUpperCase() + s.slice(1) : '';
        const comp1Label = capitalize(params.comp1);
        const comp2Label = capitalize(params.comp2);

        const interp = (xTarget:number, xArr:number[], yArr:number[]): number|null => {
            if (!xArr || !yArr || xArr.length < 2) return null;
            if (xTarget < Math.min(...xArr) || xTarget > Math.max(...xArr)) return null;
            for (let i=0;i<xArr.length-1;i++) {
                const x1=xArr[i], x2=xArr[i+1];
                if ((x1<=xTarget&&xTarget<=x2)||(x2<=xTarget&&xTarget<=x1)) {
                    const y1=yArr[i], y2=yArr[i+1];
                    if (Math.abs(x2-x1)<1e-12) return y1;
                    return y1 + (y2-y1)*(xTarget-x1)/(x2-x1);
                }
            }
            return null;
        };

        let tC: number[] | undefined;
        let pressBar: number[] | undefined;
        let dewTArray: number[] | undefined;
        let dewPArray: number[] | undefined;
        
        const sortedData = data.x.map((x: number, i: number) => ({
            x: x,
            y: data.y[i],
            t: data.t ? data.t[i] : undefined,
            p: data.p ? data.p[i] : undefined,
        })).sort((a: { x: number; }, b: { x: number; }) => a.x - b.x);
        
        const sortedX = sortedData.map((d: { x: any; }) => d.x);
        const sortedY = sortedData.map((d: { y: any; }) => d.y);

        if (diagramType === 'txy' && data.t) {
            yAxisName = "Temperature (°C)";
            titleConditionText = `at ${formatNumberToPrecision(params.pressure ?? 0)} bar`;
            tC = sortedData.map((d: { t: number; }) => d.t - 273.15);
            
            const dewPointData = sortedData.map((d: {y: number}, i: number) => ({ y: d.y, t: tC![i] })).sort((a: {y: number}, b: {y: number}) => a.y - b.y);
            dewTArray = dewPointData.map((p: {t: number}) => p.t); // for tooltip interpolation
            
            const bubblePoints = sortedX.map((x: any, i: number) => [x, tC![i]]);
            const dewPoints = dewPointData.map((p: {y: number, t: number}) => [p.y, p.t]);
            
            // Store the plot data for tie line calculations
            plotDataRef.current = { bubble: bubblePoints, dew: dewPoints };
            
            series.push({ name: 'Bubble Point', type: 'line', data: bubblePoints, symbol: 'none', color: 'green', lineStyle: { width: 2.5 } });
            series.push({ name: 'Dew Point',    type: 'line', data: dewPoints, symbol: 'none', color: '#3b82f6', lineStyle: { width: 2.5 } });

        } else if (diagramType === 'pxy' && data.p) {
            yAxisName = "Pressure (bar)";
            titleConditionText = `at ${formatNumberToPrecision(params.temp ?? 0)} °C`;
            pressBar = sortedData.map((d: { p: number; }) => d.p/1e5);

            const dewPointData = sortedData.map((d: {y: number}, i: number) => ({ y: d.y, p: pressBar![i] })).sort((a: {y: number}, b: {y: number}) => a.y - b.y);
            dewPArray = dewPointData.map((p: { p: number; }) => p.p); // for tooltip interpolation

            const bubblePoints = sortedX.map((x: any,i: number)=>[x, pressBar![i]]);
            const dewPoints = dewPointData.map((p: {y: number, p: number}) => [p.y, p.p]);
            
            // Store the plot data for tie line calculations
            plotDataRef.current = { bubble: bubblePoints, dew: dewPoints };

            series.push({ name: 'Bubble Point', type: 'line', data: bubblePoints, symbol: 'none', color: 'green', lineStyle: { width: 2.5 } });
            series.push({ name: 'Dew Point',    type: 'line', data: dewPoints, symbol: 'none', color: '#3b82f6', lineStyle: { width: 2.5 } });
        } else {
            yAxisName = `Vapor Mole Fraction ${comp1Label} (y)`;
            titleConditionText = useTemperatureForXY ? `at ${formatNumberToPrecision(params.temp ?? 0)} °C` : `at ${formatNumberToPrecision(params.pressure ?? 0)} bar`;
            series.push({ name: 'Equilibrium', type: 'line', data: sortedX.map((x: number, i: number) => [x, sortedY[i]]), symbol: 'none', color: 'red', lineStyle: { width: 2.5 } });
            series.push({ name: 'y=x', type: 'line', data: [[0, 0], [1, 1]], symbol: 'none', lineStyle: { width: 1.5, type: 'dotted' }, color: 'cyan' });
        }
        
        let yMin: number | undefined;
        let yMax: number | undefined;
        if (diagramType === 'txy' && tC) {
            const allT = tC.concat(dewTArray?.filter(t => !isNaN(t) && isFinite(t)) ?? []);
            if (allT.length > 0) {
                const minVal = Math.min(...allT);
                const maxVal = Math.max(...allT);
                const padding = (maxVal - minVal) * 0.1;
                yMin = Math.round((minVal - padding) * 100) / 100; // Round to 2 decimal places
                yMax = Math.round((maxVal + padding) * 100) / 100; // Round to 2 decimal places
            }
        } else if (diagramType === 'pxy' && pressBar) {
            const allP = pressBar.concat(dewPArray?.filter(p => !isNaN(p) && isFinite(p)) ?? []);
            if (allP.length > 0) {
                const minVal = Math.min(...allP);
                const maxVal = Math.max(...allP);
                const padding = (maxVal - minVal) * 0.1;
                yMin = Math.round((minVal - padding) * 1000) / 1000; // Round to 3 decimal places for pressure
                yMax = Math.round((maxVal + padding) * 1000) / 1000; // Round to 3 decimal places for pressure
            }
        } else if (diagramType === 'xy') {
            yMin = 0;
            yMax = 1;
        }
        // Store the calculated y-axis minimum so we can position labels/vertical lines
        yMinRef.current = yMin;
        
        const themeToUse = themeOverride ?? resolvedTheme;
        const isDark = themeToUse === 'dark';
        const textColor = isDark ? 'white' : '#000000';
        const tooltipBg = isDark ? '#08306b' : '#ffffff';
        const tooltipBorder = isDark ? '#55aaff' : '#333333';
        
        // Add tie-line series for txy and pxy diagrams
        if (diagramType === 'txy' || diagramType === 'pxy') {
            series.push({
                id: 'tie-line',
                type: 'line',
                data: [], // This series doesn't plot main points
                markLine: {
                    symbol: ['none', 'none'],
                    animation: false,
                    label: { show: false },
                    lineStyle: { type: 'dashed' },
                    emphasis: {
                        disabled: true
                    },
                    data: [] // Tie line data will be dynamically inserted here
                },
                markPoint: {
                    animation: false,
                    symbol: 'circle',
                    // This makes the dot invisible but keeps its label
                    symbolSize: 0, 
                    itemStyle: {
                        // This style is for the invisible dot's border
                        borderColor: textColor
                    },
                    label: {
                        show: true,
                        position: 'bottom',
                        backgroundColor: tooltipBg,
                        borderColor: tooltipBorder,
                        color: textColor,
                        padding: [5, 8],
                        borderRadius: 4,
                        fontFamily: 'Merriweather Sans',
                        offset: [0, 5]
                    },
                    data: [],
                }
            });
        }
        
        // ADD THIS NEW SERIES FOR THE CURSOR DOT
        series.push({
            id: 'cursor-dot',
            type: 'scatter',
            symbolSize: 10,
            itemStyle: {
                color: '#ef4444', // A nice red color
                shadowBlur: 0
            },
            data: [],
            z: 100,
            silent: true // Prevents dot from triggering tooltips
        });
        
        const displayDiagramLabel = diagramType === 'txy' ? 'Txy' : (diagramType === 'pxy' ? 'Pxy' : 'xy');
        const titleText = `${comp1Label}-${comp2Label} ${displayDiagramLabel} Diagram ${titleConditionText}`;
        const xAxisName = diagramType === 'xy' ? `Liquid Mole Fraction ${comp1Label} (x)` : `Mole Fraction ${comp1Label} (x/y)`;

        setEchartsOptions({
            backgroundColor: 'transparent',
            title: { text: titleText, left: 'center', textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } },
            grid: { top: '12%', bottom: '15%', left: '15%', right: '5%', containLabel: false },
            xAxis: {
                type: 'value', min: 0, max: 1, interval: 0.1,
                name: xAxisName,
                nameLocation: 'middle', nameGap: 40, nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
                axisLine: { lineStyle: { color: textColor } }, axisTick: { lineStyle: { color: textColor } },
                axisLabel: { color: textColor, fontSize: 14, fontFamily: 'Merriweather Sans' }, splitLine: { show: false },
                // Keep vertical axis pointer active but invisible (labels are from tooltip.axisPointer.label)
                axisPointer: {
                    show: true,
                    type: 'line',
                    lineStyle: {
                        color: 'transparent'
                    },
                    label: {
                        show: true,
                        backgroundColor: tooltipBg,
                        color: textColor,
                        borderColor: tooltipBorder,
                        fontFamily: 'Merriweather Sans',
                        formatter: (params: any) => `x/y: ${formatNumberToPrecision(Number(params.value), 3)}`
                    }
                }
            },
            yAxis: {
                type: 'value', name: yAxisName,
                min: yMin,
                max: yMax,
                interval: diagramType === 'xy' ? 0.1 : undefined,
                nameLocation: 'middle', nameGap: 60, nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
                axisLine: { lineStyle: { color: textColor } }, axisTick: { lineStyle: { color: textColor } },
                axisLabel: { 
                    color: textColor, 
                    fontSize: 14, 
                    fontFamily: 'Merriweather Sans',
                    formatter: (value: number) => {
                        // In XY mode show all labels (including 0 and 1). For other diagrams, hide the min and max values.
                        if (diagramType !== 'xy' && ((yMin !== undefined && Math.abs(value - yMin) < 0.001) || 
                            (yMax !== undefined && Math.abs(value - yMax) < 0.001))) {
                            return '';
                        }
                        return value.toString();
                    }
                }, 
                splitLine: { show: false },
                // Keep a y-axis pointer active but invisible (prevents horizontal crosshair)
                axisPointer: {
                    show: true,
                    type: 'line',
                    lineStyle: { color: 'transparent' },
                    label: {
                        show: true,
                        backgroundColor: tooltipBg,
                        color: textColor,
                        borderColor: tooltipBorder,
                        fontFamily: 'Merriweather Sans',
                        formatter: (params: any) => {
                            // Determine appropriate unit depending on diagram type or XY mode selection
                            const unit = diagramType === 'txy' ? '°C' : (diagramType === 'pxy' ? 'bar' : (useTemperatureForXY ? '°C' : 'bar'));
                            if (diagramType === 'xy') {
                                return `y: ${formatNumberToPrecision(Number(params.value), 3)}`;
                            } else {
                                const prefix = diagramType === 'txy' ? 'T' : 'P';
                                return `${prefix}: ${formatNumberToPrecision(Number(params.value), 3)} ${unit}`;
                            }
                        }
                    }
                }
            },
            legend: {
                orient: 'horizontal', bottom: 10, left: 'center',
                data: series.map(s => s.name).filter((n): n is string => !!n),
                textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
                itemHeight: 2,
            },
            tooltip: {
                trigger: 'axis', backgroundColor: tooltipBg, borderColor: tooltipBorder, borderWidth: 1,
                textStyle: { color: textColor, fontFamily: 'Merriweather Sans' },
                axisPointer: {
                    label: {
                        backgroundColor: tooltipBg,
                        color: textColor,
                        borderColor: tooltipBorder,
                        fontFamily: 'Merriweather Sans',
                        formatter: (params: any) => {
                            if (params.axisDimension === 'x') return `x/y: ${formatNumberToPrecision(params.value, 3)}`;
                            const unit = diagramType === 'txy' ? '°C' : (diagramType === 'pxy' ? 'bar' : (useTemperatureForXY ? '°C' : 'bar'));
                            const prefix = (diagramType === 'txy' || (diagramType === 'xy' && useTemperatureForXY)) ? 'T' : 'P';
                            return `${prefix}: ${formatNumberToPrecision(params.value, 3)} ${unit}`;
                        }
                    }
                },
                formatter: (params: any) => {
                    if (!params || !Array.isArray(params) || params.length === 0) return '';
                    const xVal = params[0].axisValue ?? params[0].value[0];
                    let html = `Overall Comp: ${formatNumberToPrecision(xVal, 3)}<br/>`;

                    if (diagramType === 'txy') {
                        const unit='°C';
                        const bubbleT = tC ? interp(xVal, sortedX, tC) : null;
                        const dewYPoints = sortedData.map((d: {y:number}) => d.y).sort((a: number,b: number)=>a-b);
                        const dewTPoints = sortedData.map((d: {t:number}) => d.t - 273.15).sort((a:number,b:number)=> a - (interp(b,dewYPoints, sortedData.map((d: {t:number})=>d.t-273.15))??0) );
                        const dewT = dewTArray ? interp(xVal, dewYPoints, dewTArray) : null;
                        
                        if(bubbleT!==null) html += `<span style="color: green;">Bubble Temp: ${formatNumberToPrecision(bubbleT,3)} ${unit}</span><br/>`;
                        if(dewT!==null)    html += `<span style="color: #3b82f6;">Dew Temp: ${formatNumberToPrecision(dewT,3)} ${unit}</span><br/>`;
                    } else if (diagramType === 'pxy') {
                        const unit='bar';
                        const bubbleP = pressBar ? interp(xVal, sortedX, pressBar) : null;
                        const dewYPoints = sortedData.map((d: {y: number}) => d.y).sort((a: number,b: number)=>a-b);
                        const dewP = pressBar && dewPArray ? interp(xVal, dewYPoints, dewPArray) : null;
                        if(bubbleP!==null) html += `<span style="color: green;">Bubble P: ${formatNumberToPrecision(bubbleP,3)} ${unit}</span><br/>`;
                        if(dewP!==null)    html += `<span style="color: #3b82f6;">Dew P: ${formatNumberToPrecision(dewP,3)} ${unit}</span><br/>`;
                    } else { // xy
                        html = `x: ${formatNumberToPrecision(xVal, 3)}<br/>`;
                        const yVal = interp(xVal, sortedX, sortedY);
                        if (yVal !== null) {
                            html += `<span style="color: red;">y: ${formatNumberToPrecision(yVal,3)}</span><br/>`;
                        }
                    }

                    // Add Lever Rule Info
                    const leverData = leverRuleDataRef.current;
                    if (leverData && leverData.liquid !== null) {
                        html += `<hr style="border-color: #888; margin-top: 5px; margin-bottom: 5px;" />`;
                        // Match bubble curve color (green) for liquid phase
                        html += `<span style="color: green;">Liquid Phase: <strong style=\"color: green;\">${formatNumberToPrecision(leverData.liquid * 100, 3)}%</strong></span><br/>`;
                        // Match dew curve color (#3b82f6) for vapor phase
                        html += `<span style="color: #3b82f6;">Vapor Phase: <strong style=\"color: #3b82f6;\">${formatNumberToPrecision(leverData.vapor! * 100, 3)}%</strong></span>`;
                    }

                    return html;
                }
            },
            series: series,
            animation: false,
            toolbox: {
                show: false
            },
        });
    }, [diagramType, useTemperatureForXY, displayedParams, resolvedTheme]);

    // ...existing code... (event handlers replaced by memoized versions below)

    async function fetchCompoundDataLocal(compoundName: string): Promise<CompoundData | null> {
        if (!supabase) { throw new Error("Supabase client not initialized."); }
        console.log(`PhaseDiagrams: Fetching data for ${compoundName}...`);
        try {
            const result = await fetchCompoundDataFromHysys(supabase, compoundName);
            return result;
        } catch (err: any) {
            console.error(`Error fetching data for ${compoundName}:`, err.message);
            setError(`Data fetch failed for ${compoundName}: ${err.message}`);
            return null;
        }
    }


    
    const generateDiagram = useCallback(async (showLoading: boolean = true) => {
        if (showLoading) {
            setLoading(true);
            setChartData(null);
        }
        setError(null);
        if (!comp1Name || !comp2Name) {
            setError("Please enter both compound names.");
            if (showLoading) setLoading(false);
            return;
        }

        try {
            let d1 = comp1Data, d2 = comp2Data, params = interactionParams;
            const needsFetching = !d1 || !d2 || !params || d1.name !== comp1Name || d2.name !== comp2Name || displayedParams.package !== fluidPackage;
            
            if (needsFetching) {
                console.log("Cache miss or inputs changed. Fetching new data...");
                // Resolve display names to SimNames
                const [sim1, sim2] = await Promise.all([resolveSimName(supabase, comp1Name), resolveSimName(supabase, comp2Name)]);
                const resolvedComp1 = sim1 || comp1Name;
                const resolvedComp2 = sim2 || comp2Name;
                const [fetched1, fetched2] = await Promise.all([fetchCompoundDataLocal(resolvedComp1), fetchCompoundDataLocal(resolvedComp2)]);
                if (!fetched1 || !fetched2) { if (showLoading) setLoading(false); return; }
                
                let fetchedParams: any;
                if (fluidPackage === 'unifac') {
                    if (!fetched1.unifacGroups || !fetched2.unifacGroups) throw new Error("UNIFAC groups missing.");
                    const allSubgroupIds = new Set<number>();
                    [fetched1, fetched2].forEach(comp => { if (comp.unifacGroups) Object.keys(comp.unifacGroups).forEach(id => allSubgroupIds.add(parseInt(id))); });
                    fetchedParams = await fetchUnifacInteractionParams(supabase, Array.from(allSubgroupIds));
                } else if (fluidPackage === 'nrtl') {
                    fetchedParams = await fetchNrtlParameters(supabase, fetched1.name, fetched2.name);
                    if ((fetchedParams as any)?._usedUnifacFallback) {
                        setError(`Warning: NRTL parameters not found in database for ${fetched1.name}–${fetched2.name}. Using UNIFAC-based estimate.`);
                    }
                } else if (fluidPackage === 'pr') {
                    fetchedParams = await fetchPrInteractionParams(supabase, fetched1.name, fetched2.name);
                } else if (fluidPackage === 'srk') {
                    fetchedParams = await fetchSrkInteractionParams(supabase, fetched1.name, fetched2.name);
                } else if (fluidPackage === 'uniquac') {
                    fetchedParams = await fetchUniquacInteractionParams(supabase, fetched1.name, fetched2.name);
                    if ((fetchedParams as any)?._usedUnifacFallback) {
                        setError(`Warning: UNIQUAC parameters not found in database for ${fetched1.name}–${fetched2.name}. Using UNIFAC-based estimate.`);
                    }
                } else if (fluidPackage === 'wilson') {
                    fetchedParams = await fetchWilsonInteractionParams(supabase, fetched1.name, fetched2.name);
                }
                
                setComp1Data(fetched1); setComp2Data(fetched2); setInteractionParams(fetchedParams);
                d1 = fetched1; d2 = fetched2; params = fetchedParams;
            }
            if (!d1 || !d2 || !params) { throw new Error("Component data or parameters not available."); }
            const components: [CompoundData, CompoundData] = [d1, d2];
            
            const pointsCount = 101; // 101 points => 0.01 spacing for x
            // Use fine resolution matching McCabe-Thiele (0.01 spacing)
            const x_feed = Array.from({ length: pointsCount }, (_, i) => i / (pointsCount - 1));
            
            const results: { x: number[], y: number[], t?: number[], p?: number[] } = { x: [], y: [] };
            if (diagramType === 'txy') results.t = []; else if (diagramType === 'pxy') results.p = [];

            for (const x1 of x_feed) {
                let resultPoint: BubbleDewResult | null = null;
                if (x1 === 0 || x1 === 1) {
                    const pureComp = x1 === 0 ? components[1] : components[0];
                    if (diagramType === 'txy') {
                        const Tbub = antoineBoilingPointSolverLocal(pureComp.antoine, pressureBar! * 1e5);
                        if (Tbub) { resultPoint = { comp1_feed: x1, comp1_equilibrium: x1, T_K: Tbub }; }
                    } else if (diagramType === 'pxy') {
                        const Pbub = libCalculatePsat_Pa(pureComp.antoine!, temperatureC! + 273.15);
                        if (Pbub) { resultPoint = { comp1_feed: x1, comp1_equilibrium: x1, P_Pa: Pbub }; }
                    } else { // xy
                        resultPoint = { comp1_feed: x1, comp1_equilibrium: x1 };
                    }
                } else {
                     if (diagramType === 'txy') {
                        const Tbp1 = antoineBoilingPointSolverLocal(d1.antoine, pressureBar! * 1e5);
                        const Tbp2 = antoineBoilingPointSolverLocal(d2.antoine, pressureBar! * 1e5);
                        const initialTempGuess = Tbp1 && Tbp2 ? (x1 * Tbp1 + (1 - x1) * Tbp2) : 373.15;
                        if(fluidPackage==='unifac') resultPoint=calculateBubbleTemperatureUnifac(components,x1,pressureBar!*1e5,params,initialTempGuess);
                        else if(fluidPackage==='nrtl') resultPoint=calculateBubbleTemperatureNrtl(components,x1,pressureBar!*1e5,params,initialTempGuess);
                        else if(fluidPackage==='pr') resultPoint=calculateBubbleTemperaturePr(components,x1,pressureBar!*1e5,params,initialTempGuess);
                        else if(fluidPackage==='srk') resultPoint=calculateBubbleTemperatureSrk(components,x1,pressureBar!*1e5,params,initialTempGuess);
                        else if(fluidPackage==='uniquac') resultPoint=calculateBubbleTemperatureUniquac(components,x1,pressureBar!*1e5,params,initialTempGuess);
                        else if(fluidPackage==='wilson') resultPoint=calculateBubbleTemperatureWilson(components,x1,pressureBar!*1e5,params,initialTempGuess);
                    } else { // pxy or xy
                        const fixedTempK = useTemperatureForXY || diagramType === 'pxy' ? temperatureC! + 273.15 : null;
                        const fixedPressPa = !useTemperatureForXY && diagramType === 'xy' ? pressureBar! * 1e5 : null;
                        if (fixedTempK) {
                            const initialPressureGuess = (libCalculatePsat_Pa(d1.antoine!, fixedTempK) + libCalculatePsat_Pa(d2.antoine!, fixedTempK)) / 2 || 101325;
                            if(fluidPackage==='unifac') resultPoint=calculateBubblePressureUnifac(components,x1,fixedTempK,params);
                            else if(fluidPackage==='nrtl') resultPoint=calculateBubblePressureNrtl(components,x1,fixedTempK,params);
                            else if(fluidPackage==='pr') resultPoint=calculateBubblePressurePr(components,x1,fixedTempK,params,initialPressureGuess);
                            else if(fluidPackage==='srk') resultPoint=calculateBubblePressureSrk(components,x1,fixedTempK,params,initialPressureGuess);
                            else if(fluidPackage==='uniquac') resultPoint=calculateBubblePressureUniquac(components,x1,fixedTempK,params);
                            else if(fluidPackage==='wilson') resultPoint=calculateBubblePressureWilson(components,x1,fixedTempK,params);
                        } else if (fixedPressPa) {
                            const Tbp1 = antoineBoilingPointSolverLocal(d1.antoine, fixedPressPa);
                            const Tbp2 = antoineBoilingPointSolverLocal(d2.antoine, fixedPressPa);
                            const initialTempGuess = Tbp1 && Tbp2 ? (x1 * Tbp1 + (1 - x1) * Tbp2) : 373.15;
                            if(fluidPackage==='unifac') resultPoint=calculateBubbleTemperatureUnifac(components,x1,fixedPressPa,params,initialTempGuess);
                            else if(fluidPackage==='nrtl') resultPoint=calculateBubbleTemperatureNrtl(components,x1,fixedPressPa,params,initialTempGuess);
                            else if(fluidPackage==='pr') resultPoint=calculateBubbleTemperaturePr(components,x1,fixedPressPa,params,initialTempGuess);
                            else if(fluidPackage==='srk') resultPoint=calculateBubbleTemperatureSrk(components,x1,fixedPressPa,params,initialTempGuess);
                            else if(fluidPackage==='uniquac') resultPoint=calculateBubbleTemperatureUniquac(components,x1,fixedPressPa,params,initialTempGuess);
                            else if(fluidPackage==='wilson') resultPoint=calculateBubbleTemperatureWilson(components,x1,fixedPressPa,params,initialTempGuess);
                        }
                    }
                }
                if (resultPoint?.error === undefined && typeof resultPoint?.comp1_equilibrium === 'number') {
                    results.x.push(x1);
                    results.y.push(resultPoint.comp1_equilibrium);
                    if (diagramType === 'txy' && results.t && resultPoint.T_K) results.t.push(resultPoint.T_K);
                    else if (diagramType === 'pxy' && results.p && resultPoint.P_Pa) results.p.push(resultPoint.P_Pa);
                }
            }
            
            setChartData(results);
            const newParams = {
                comp1: comp1DisplayName || d1.name, comp2: comp2DisplayName || d2.name,
                temp: diagramType === 'pxy' || (diagramType === 'xy' && useTemperatureForXY) ? temperatureC : null,
                pressure: diagramType === 'txy' || (diagramType === 'xy' && !useTemperatureForXY) ? pressureBar : null,
                package: fluidPackage, type: diagramType
            };
            setDisplayedParams(newParams);
            generateEchartsOptions(results, newParams, resolvedTheme);

        } catch (err: any) {
            setError(`Calculation failed: ${err.message}`);
        } finally {
            if (showLoading) setLoading(false);
        }
    }, [comp1Name, comp2Name, diagramType, fluidPackage, pressureBar, temperatureC, useTemperatureForXY, resolvedTheme, comp1Data, comp2Data, interactionParams, displayedParams.package, generateEchartsOptions]);

    useEffect(() => {
        if (supabase && autoGenerateOnCompChange) {
            generateDiagram();
            setAutoGenerateOnCompChange(false);
        }
    }, [supabase, autoGenerateOnCompChange, generateDiagram]);

    // Regenerate chart options on theme change
    useEffect(() => {
        if (chartData) {
            generateEchartsOptions(chartData, displayedParams, resolvedTheme);
        }
    }, [resolvedTheme, chartData, displayedParams, generateEchartsOptions]);

    const debouncedSliderUpdate = useCallback(debounce(() => generateDiagram(false), 200), [generateDiagram]);

    useEffect(() => {
        if (diagramType === 'txy' || (diagramType === 'xy' && !useTemperatureForXY)) debouncedSliderUpdate();
    }, [pressureBar]);

    useEffect(() => {
        if (diagramType === 'pxy' || (diagramType === 'xy' && useTemperatureForXY)) debouncedSliderUpdate();
    }, [temperatureC]);
    
    // Auto-generate when diagram type or fluid package changes
    useEffect(() => {
        generateDiagram();
    }, [diagramType, fluidPackage, useTemperatureForXY]);


    const formatNumberToPrecision = (num: any, precision: number = 4): string => {
        if (typeof num === 'number' && !isNaN(num)) {
          if (num === 0) return '0';
          const fixed = num.toPrecision(precision);
          if (fixed.includes('.')) {
            return parseFloat(fixed).toString(); 
          }
          return fixed;
        }
        return "N/A";
    };
    
    const handleKeyDown = (event: React.KeyboardEvent<HTMLInputElement>) => {
        if (event.key === 'Enter') {
            event.preventDefault();
            generateDiagram();
        }
    };

    const renderConditionalInput = () => {
        switch (diagramType) {
            case 'txy':
                return (
                    <div className="space-y-2">
                        <div className="flex items-center justify-between">
                            <Label htmlFor="pressure-slider">Pressure (bar): {formatNumberToPrecision(pressureBar ?? 0)}</Label>
                            <div className="flex items-center gap-1">
                                <Label className="text-xs text-muted-foreground">Max:</Label>
                                <Input
                                    type="text"
                                    value={pressureMax}
                                    onChange={(e) => {
                                        const raw = e.target.value;
                                        const num = parseFloat(raw);
                                        const CLAMP_MAX = 2000;
                                        if (raw === "" || isNaN(num)) {
                                            setPressureMax(raw);
                                        } else {
                                            setPressureMax(String(Math.min(num, CLAMP_MAX)));
                                        }
                                    }}
                                    className="w-16 h-8 text-xs"
                                />
                            </div>
                        </div>
                        <Slider id="pressure-slider" min={0.1} max={parseFloat(pressureMax) || 10} step={computeStep(parseFloat(pressureMax) || 10)} value={[pressureBar || 1]} onValueChange={([v]) => setPressureBar(v)} className="w-full"/>
                    </div>
                );
            case 'pxy':
                return (
                    <div className="space-y-2">
                        <div className="flex items-center justify-between">
                            <Label htmlFor="temperature-slider">Temperature (°C): {formatNumberToPrecision(temperatureC ?? 0)}</Label>
                            <div className="flex items-center gap-1">
                                <Label className="text-xs text-muted-foreground">Max:</Label>
                                <Input
                                    type="text"
                                    value={tempMax}
                                    onChange={(e) => {
                                        const raw = e.target.value;
                                        const num = parseFloat(raw);
                                        const CLAMP_MAX = 2000;
                                        if (raw === "" || isNaN(num)) {
                                            setTempMax(raw);
                                        } else {
                                            setTempMax(String(Math.min(num, CLAMP_MAX)));
                                        }
                                    }}
                                    className="w-16 h-8 text-xs"
                                />
                            </div>
                        </div>
                        <Slider id="temperature-slider" min={-50} max={parseFloat(tempMax) || 200} step={computeStep(parseFloat(tempMax) || 200)} value={[temperatureC || 60]} onValueChange={([v]) => setTemperatureC(v)} className="w-full"/>
                    </div>
                );
            case 'xy':
                return (
                    <div className="space-y-2">
                        <div className="flex items-center justify-between">
                            <div className="flex items-center gap-2">
                                <Tabs value={useTemperatureForXY ? "temperature" : "pressure"} onValueChange={v => {
                                    setUseTemperatureForXY(v === "temperature");
                                    generateDiagram();
                                }}>
                                    <TabsList><TabsTrigger value="temperature">T</TabsTrigger><TabsTrigger value="pressure">P</TabsTrigger></TabsList>
                                </Tabs>
                                <span className="text-sm">
                                    {useTemperatureForXY
                                        ? `${formatNumberToPrecision(temperatureC ?? 0)} °C`
                                        : `${formatNumberToPrecision(pressureBar ?? 0)} bar`}
                                </span>
                            </div>
                            <div className="flex items-center gap-1">
                                <Label className="text-xs text-muted-foreground">Max:</Label>
                                <Input
                                    type="text"
                                    value={useTemperatureForXY ? tempMax : pressureMax}
                                    onChange={(e) => {
                                        const raw = e.target.value;
                                        const num = parseFloat(raw);
                                        const CLAMP_MAX = 2000;
                                        if (useTemperatureForXY) {
                                            if (raw === "" || isNaN(num)) setTempMax(raw);
                                            else setTempMax(String(Math.min(num, CLAMP_MAX)));
                                        } else {
                                            if (raw === "" || isNaN(num)) setPressureMax(raw);
                                            else setPressureMax(String(Math.min(num, CLAMP_MAX)));
                                        }
                                    }}
                                    className="w-16 h-8 text-xs"
                                />
                            </div>
                        </div>
                        <Slider
                            id="xy-slider"
                            min={useTemperatureForXY ? -50 : 0.1}
                            max={useTemperatureForXY ? parseFloat(tempMax) || 200 : parseFloat(pressureMax) || 10}
                            step={computeStep(useTemperatureForXY ? parseFloat(tempMax) || 200 : parseFloat(pressureMax) || 10)}
                            value={[useTemperatureForXY ? temperatureC ?? 60 : pressureBar ?? 1]}
                            onValueChange={([v]) => useTemperatureForXY ? setTemperatureC(v) : setPressureBar(v)}
                            className="w-full"
                        />
                    </div>
                );
            default: return null;
        }
    };

    const fetchSuggestions = useCallback(async (inputValue: string, inputTarget: 'comp1' | 'comp2') => {
        if (!inputValue || inputValue.length < 2) {
            if (inputTarget === 'comp1') setShowComp1Suggestions(false); else setShowComp2Suggestions(false);
            return;
        }
        try {
            const suggestions = await fetchCompoundSuggestions(supabase, inputValue);
            if (inputTarget === 'comp1') {
                setComp1Suggestions(suggestions);
                setShowComp1Suggestions(suggestions.length > 0);
            } else {
                setComp2Suggestions(suggestions);
                setShowComp2Suggestions(suggestions.length > 0);
            }
        } catch (err) { console.error("Error fetching suggestions:", err); }
    }, []);

    const handleComp1NameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const newValue = e.target.value;
        setComp1Name(newValue);
        setComp1DisplayName(newValue);
        setActiveSuggestionInput('comp1');
        fetchSuggestions(newValue, 'comp1');
    };

    const handleComp2NameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const newValue = e.target.value;
        setComp2Name(newValue);
        setComp2DisplayName(newValue);
        setActiveSuggestionInput('comp2');
        fetchSuggestions(newValue, 'comp2');
    };

    const handleSuggestionClick = (alias: CompoundAlias, inputTarget: 'comp1' | 'comp2') => {
        setAutoGenerateOnCompChange(true);
        if (inputTarget === 'comp1') {
            setComp1Name(alias.simName);
            setComp1DisplayName(formatCompoundName(alias.fullName));
            setShowComp1Suggestions(false);
        } else {
            setComp2Name(alias.simName);
            setComp2DisplayName(formatCompoundName(alias.fullName));
            setShowComp2Suggestions(false);
        }
    };

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
        return () => document.removeEventListener("mousedown", handleClickOutside);
    }, [activeSuggestionInput]);

    // Memoized event handlers for the chart
    const handleAxisPointerUpdate = useCallback((params: any) => {
        const echartsInstance = echartsRef.current?.getEchartsInstance()
        if (!echartsInstance || !params.axesInfo || params.axesInfo.length < 2) {
            return
        }

        const cursorX = params.axesInfo[0].value
        const cursorY = params.axesInfo[1].value

        // Logic for Txy and Pxy diagrams
        if (diagramType === 'txy' || diagramType === 'pxy') {
            const { bubble, dew } = plotDataRef.current
            if (bubble.length < 2 || dew.length < 2) {
                // If plot data is not ready, just update the dot
                echartsInstance.setOption({
                    series: [{ id: 'cursor-dot', data: [[cursorX, cursorY]] }]
                })
                return
            }

            const interpInverse = (yTarget: number, data: [number, number][]): number | null => {
                for (let i = 0; i < data.length - 1; i++) {
                    const [x1, y1] = data[i]
                    const [x2, y2] = data[i + 1]
                    if ((y1 <= yTarget && yTarget <= y2) || (y2 <= yTarget && yTarget <= y1)) {
                        if (Math.abs(y2 - y1) < 1e-9) return x1
                        const t = (yTarget - y1) / (y2 - y1)
                        return x1 + t * (x2 - x1)
                    }
                }
                return null
            }

            const bubbleX = interpInverse(cursorY, bubble)
            const dewX = interpInverse(cursorY, dew)
            const yAxisMin = yMinRef.current ?? 0

            const inRegion = bubbleX !== null && dewX !== null && cursorX >= Math.min(bubbleX, dewX) && cursorX <= Math.max(bubbleX, dewX)

            // *** UNCOMMENT THE LINE BELOW TO DEBUG IN YOUR BROWSER'S CONSOLE ***
            // console.log(`Cursor: (${cursorX.toFixed(3)}, ${cursorY.toFixed(3)}), BubbleX: ${bubbleX?.toFixed(3)}, DewX: ${dewX?.toFixed(3)}, InRegion: ${inRegion}`)

            if (inRegion) {
                // Added Math.abs for robustness in lever rule calculation
                leverRuleDataRef.current = {
                    liquid: Math.abs((dewX - cursorX) / (dewX - bubbleX)),
                    vapor: Math.abs((cursorX - bubbleX) / (dewX - bubbleX))
                }
                
                const tieLineData = [
                    // Horizontal segment from bubble curve to cursor
                    [{ xAxis: bubbleX, yAxis: cursorY, lineStyle: { color: 'green', type: [4, 4] } }, { xAxis: cursorX, yAxis: cursorY }],
                    // Horizontal segment from cursor to dew curve
                    [{ xAxis: cursorX, yAxis: cursorY, lineStyle: { color: '#3b82f6', type: [4, 4] } }, { xAxis: dewX, yAxis: cursorY }],
                    // Vertical segment from bubble curve to x-axis
                    [{ xAxis: bubbleX, yAxis: cursorY, lineStyle: { color: 'green', type: [4, 4] } }, { xAxis: bubbleX, yAxis: yAxisMin }],
                    // Vertical segment from dew curve to x-axis
                    [{ xAxis: dewX, yAxis: cursorY, lineStyle: { color: '#3b82f6', type: [4, 4] } }, { xAxis: dewX, yAxis: yAxisMin }]
                ]
                
                const intersectionPoints = [
                    {
                        name: 'BubbleIntersection',
                        xAxis: bubbleX,
                        yAxis: cursorY,
                        symbol: 'circle',
                        symbolSize: 8,
                        label: { show: false },
                        emphasis: {
                            label: { show: false }
                        },
                        silent: true,
                        itemStyle: { color: 'green', borderColor: 'transparent', borderWidth: 0 },
                        z: 100
                    },
                    {
                        name: 'DewIntersection',
                        xAxis: dewX,
                        yAxis: cursorY,
                        symbol: 'circle',
                        symbolSize: 8,
                        label: { show: false },
                        emphasis: {
                            label: { show: false }
                        },
                        silent: true,
                        itemStyle: { color: '#3b82f6', borderColor: 'transparent', borderWidth: 0 },
                        z: 100
                    }
                ]

                // Add bottom labels at the vertical line intersections with x-axis
                const isDark = resolvedTheme === 'dark'
                const labelBg = isDark ? '#08306b' : '#ffffff'
                const labelBorder = isDark ? '#55aaff' : '#333333'
                const labelColor = isDark ? 'white' : '#000000'
                
                const bottomLabels = [
                    { 
                        name: 'BubbleLabel', 
                        value: `x: ${formatNumberToPrecision(bubbleX, 3)}`, 
                        xAxis: bubbleX, 
                        yAxis: yAxisMin,
                        label: {
                            show: true,
                            position: 'bottom',
                            backgroundColor: labelBg,
                            borderColor: labelBorder,
                            color: labelColor,
                            padding: [5, 8],
                            borderRadius: 4,
                            fontFamily: 'Merriweather Sans',
                            offset: [0, 5]
                        },
                        itemStyle: { color: 'transparent' },
                        symbolSize: 0
                    },
                    { 
                        name: 'DewLabel', 
                        value: `y: ${formatNumberToPrecision(dewX, 3)}`, 
                        xAxis: dewX, 
                        yAxis: yAxisMin,
                        label: {
                            show: true,
                            position: 'bottom',
                            backgroundColor: labelBg,
                            borderColor: labelBorder,
                            color: labelColor,
                            padding: [5, 8],
                            borderRadius: 4,
                            fontFamily: 'Merriweather Sans',
                            offset: [0, 5]
                        },
                        itemStyle: { color: 'transparent' },
                        symbolSize: 0
                    }
                ]

                const markPointData = [ ...intersectionPoints, ...bottomLabels ]

                echartsInstance.setOption({
                    // Hide the x-axis axisPointer label while inside the two-phase region
                    xAxis: { axisPointer: { label: { show: false } } },
                    yAxis: { axisPointer: { label: { show: true } } },
                    series: [
                        { id: 'tie-line', markLine: { data: tieLineData }, markPoint: { data: markPointData } },
                        { id: 'cursor-dot', data: [[cursorX, cursorY]] }
                    ]
                })
            } else {
                leverRuleDataRef.current = { liquid: null, vapor: null }
                echartsInstance.setOption({
                    xAxis: { axisPointer: { label: { show: true } } },
                    yAxis: { axisPointer: { label: { show: true } } },
                    series: [
                        { id: 'tie-line', markLine: { data: [] }, markPoint: { data: [] } },
                        { id: 'cursor-dot', data: [[cursorX, cursorY]] }
                    ]
                })
            }
        // Logic for other diagrams (like xy)
        } else {
            echartsInstance.setOption({
                series: [
                    { id: 'cursor-dot', data: [[cursorX, cursorY]] }
                ]
            })
        }
    }, [diagramType, resolvedTheme])

    const handleChartMouseOut = useCallback(() => {
        const echartsInstance = echartsRef.current?.getEchartsInstance();
        if (!echartsInstance) return;
        leverRuleDataRef.current = { liquid: null, vapor: null };
        echartsInstance.setOption({
            series: [
                { id: 'tie-line', markLine: { data: [] } },
                { id: 'cursor-dot', data: [] }
            ]
        });
    }, []);

    const onEvents = useMemo(() => ({
        'updateAxisPointer': handleAxisPointerUpdate,
        'mouseout': handleChartMouseOut,
    }), [handleAxisPointerUpdate, handleChartMouseOut]);

    return (
        <div className="container mx-auto p-4 md:p-8">
            <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                <div className="lg:col-span-1 space-y-6">
                    <Card>
                        <CardHeader><CardTitle>VLE Diagram Generator</CardTitle></CardHeader>
                        <CardContent className="space-y-6">
                            <Tabs value={diagramType} onValueChange={(v) => setDiagramType(v as DiagramType)} className="w-full">
                                <TabsList className="grid w-full grid-cols-3"><TabsTrigger value="txy">Txy</TabsTrigger><TabsTrigger value="pxy">Pxy</TabsTrigger><TabsTrigger value="xy">xy</TabsTrigger></TabsList>
                            </Tabs>
                            <div className="space-y-4">
                                <div className="flex items-center gap-2">
                                    <div className="relative flex-1">
                                        <Input
                                            ref={input1Ref} id="comp1Name" value={comp1DisplayName} onChange={handleComp1NameChange} onKeyDown={handleKeyDown}
                                            onFocus={() => { setActiveSuggestionInput('comp1'); fetchSuggestions(comp1DisplayName, 'comp1'); }}
                                            placeholder="Methanol" required className="w-full" autoComplete="off"
                                        />
                                        {showComp1Suggestions && comp1Suggestions.length > 0 && (
                                            <div ref={activeSuggestionInput === 'comp1' ? activeComponentRef : null} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1">
                                                {comp1Suggestions.map((alias, index) => (
                                                    <div key={index} onClick={() => handleSuggestionClick(alias, 'comp1')} className="px-3 py-2 hover:bg-accent cursor-pointer text-sm">
                                                        {formatCompoundName(alias.fullName)}
                                                    </div>
                                                ))}
                                            </div>
                                        )}
                                    </div>
                                    <Button variant="ghost" size="icon" onClick={() => {
                                        const temp = comp1Name;
                                        const tempDisplay = comp1DisplayName;
                                        setComp1Name(comp2Name);
                                        setComp1DisplayName(comp2DisplayName);
                                        setComp2Name(temp);
                                        setComp2DisplayName(tempDisplay);
                                        setAutoGenerateOnCompChange(true);
                                        // REMOVED: generateDiagram();
                                        // By removing the direct call, we let the useEffect hook handle
                                        // the diagram generation after the state has been correctly updated.
                                    }} title="Swap Components"><ArrowLeftRight className="h-4 w-4" /></Button>
                                    <div className="relative flex-1">
                                        <Input
                                            ref={input2Ref} id="comp2Name" value={comp2DisplayName} onChange={handleComp2NameChange} onKeyDown={handleKeyDown}
                                            onFocus={() => { setActiveSuggestionInput('comp2'); fetchSuggestions(comp2DisplayName, 'comp2'); }}
                                            placeholder="Water" required className="w-full" autoComplete="off"
                                        />
                                        {showComp2Suggestions && comp2Suggestions.length > 0 && (
                                            <div ref={activeSuggestionInput === 'comp2' ? activeComponentRef : null} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1">
                                                {comp2Suggestions.map((alias, index) => (
                                                    <div key={index} onClick={() => handleSuggestionClick(alias, 'comp2')} className="px-3 py-2 hover:bg-accent cursor-pointer text-sm">
                                                        {formatCompoundName(alias.fullName)}
                                                    </div>
                                                ))}
                                            </div>
                                        )}
                                    </div>
                                </div>
                                {renderConditionalInput()}
                                <div className="flex items-center gap-2 w-full">
                                    <Label htmlFor="fluidPackage" className="text-sm font-medium whitespace-nowrap">Fluid Package:</Label>
                                    <div className="flex-1">
                                        <Select value={fluidPackage} onValueChange={(v) => {
                                            setFluidPackage(v as FluidPackageType);
                                            // The generateDiagram() call is removed from here.
                                            // The useEffect hook will handle the regeneration
                                            // when the fluidPackage state changes.
                                        }}>
                                            <SelectTrigger id="fluidPackage" className="w-full"><SelectValue /></SelectTrigger>
                                            <SelectContent>
                                                <SelectItem value="wilson">Wilson</SelectItem>
                                                <SelectItem value="uniquac">UNIQUAC</SelectItem>
                                                <SelectItem value="nrtl">NRTL</SelectItem>
                                                <SelectItem value="pr">Peng-Robinson</SelectItem>
                                                <SelectItem value="srk">SRK</SelectItem>
                                            </SelectContent>
                                        </Select>
                                    </div>
                                </div>
                            </div>
                            <Button onClick={() => generateDiagram()} disabled={loading} className="w-full">{loading ? 'Calculating...' : 'Generate Diagram'}</Button>
                            {error && <p className="text-sm text-red-500 mt-2">{error}</p>}
                        </CardContent>
                    </Card>
                    {chartData && diagramType === 'txy' && chartData.t && (
                      <Card>
                        <CardHeader><CardTitle>Pure Component Boiling Points</CardTitle></CardHeader>
                        <CardContent>
                                                                                <div className="grid grid-cols-2 gap-4 text-center">
                                                                                    <div className="p-4 bg-muted rounded-md responsive-text-container text-center">
                                                            <p className="text-sm font-medium">{comp1Name.charAt(0).toUpperCase() + comp1Name.slice(1)}</p>
                                                            <p className="font-bold whitespace-nowrap shrinkable-text">{formatNumberToPrecision(chartData.t[chartData.x.indexOf(1)] - 273.15, 3)} °C</p>
                                                                                    </div>
                                                                                    <div className="p-4 bg-muted rounded-md responsive-text-container text-center">
                                                            <p className="text-sm font-medium">{comp2Name.charAt(0).toUpperCase() + comp2Name.slice(1)}</p>
                                                            <p className="font-bold whitespace-nowrap shrinkable-text">{formatNumberToPrecision(chartData.t[chartData.x.indexOf(0)] - 273.15, 3)} °C</p>
                                                                                    </div>
                                                    </div>
                        </CardContent>
                      </Card>
                    )}
                    {chartData && diagramType === 'pxy' && chartData.p && (
                      <Card>
                        <CardHeader><CardTitle>Pure Component Vapor Pressures</CardTitle></CardHeader>
                        <CardContent>
                                                                                <div className="grid grid-cols-2 gap-4 text-center">
                                                                                    <div className="p-4 bg-muted rounded-md responsive-text-container text-center">
                                                            <p className="text-sm font-medium">{comp1Name.charAt(0).toUpperCase() + comp1Name.slice(1)}</p>
                                                            <p className="font-bold whitespace-nowrap shrinkable-text">{formatNumberToPrecision(chartData.p[chartData.x.indexOf(1)] / 1e5, 3)} bar</p>
                                                                                    </div>
                                                                                    <div className="p-4 bg-muted rounded-md responsive-text-container text-center">
                                                            <p className="text-sm font-medium">{comp2Name.charAt(0).toUpperCase() + comp2Name.slice(1)}</p>
                                                            <p className="font-bold whitespace-nowrap shrinkable-text">{formatNumberToPrecision(chartData.p[chartData.x.indexOf(0)] / 1e5, 3)} bar</p>
                                                                                    </div>
                                                    </div>
                        </CardContent>
                      </Card>
                    )}
                </div>
                <div className="lg:col-span-2">
                    <Card className="h-full"><CardContent className="py-2 h-full">
                                <div className="relative aspect-square rounded-md h-full z-10">
                           {loading && ( <div className="absolute inset-0 flex items-center justify-center text-muted-foreground"><div className="text-center"><div className="mb-2">Loading & Calculating...</div><div className="text-sm text-muted-foreground/70">Using { fluidPackage.toUpperCase()} model.</div></div></div> )}
                           {!loading && !chartData && !error && ( <div className="absolute inset-0 flex items-center justify-center text-muted-foreground">Provide inputs and generate a diagram.</div> )}
                           {error && !loading && ( <div className="absolute inset-0 flex items-center justify-center text-red-400">Error: {error}</div> )}
                           {!loading && chartData && Object.keys(echartsOptions).length > 0 && (
                            <ReactECharts 
                                ref={echartsRef}
                                echarts={echarts} 
                                style={{ height: '100%', width: '100%' }} 
                                option={echartsOptions} 
                                notMerge={false} 
                                lazyUpdate={true}
                                onEvents={onEvents}
                            />
                           )}
                        </div>
                    </CardContent></Card>
                </div>
            </div>
        </div>
    );
}