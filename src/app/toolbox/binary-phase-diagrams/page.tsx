'use client';

import React, { useState, useEffect, useCallback, useRef } from 'react';
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
    fetchNrtlParameters,
    calculateBubbleTemperatureNrtl,
    calculateBubblePressureNrtl,
    fetchPrInteractionParams,
    calculateBubbleTemperaturePr,
    calculateBubblePressurePr,
    fetchSrkInteractionParams,
    calculateBubbleTemperatureSrk,
    calculateBubblePressureSrk,
    fetchUniquacInteractionParams,
    calculateBubbleTemperatureUniquac,
    calculateBubblePressureUniquac,
    fetchWilsonInteractionParams,
    calculateBubbleTemperatureWilson,
    calculateBubblePressureWilson,
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
    const [comp1Name, setComp1Name] = useState('methanol');
    const [comp2Name, setComp2Name] = useState('water');
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

    const [comp1Suggestions, setComp1Suggestions] = useState<string[]>([]);
    const [comp2Suggestions, setComp2Suggestions] = useState<string[]>([]);
    const [showComp1Suggestions, setShowComp1Suggestions] = useState(false);
    const [showComp2Suggestions, setShowComp2Suggestions] = useState(false);
    const [activeSuggestionInput, setActiveSuggestionInput] = useState<'comp1' | 'comp2' | null>(null);
    const [autoGenerateOnCompChange, setAutoGenerateOnCompChange] = useState(true); // Auto-generate on first load

    const input1Ref = useRef<HTMLInputElement>(null);
    const input2Ref = useRef<HTMLInputElement>(null);
    const suggestions1Ref = useRef<HTMLDivElement>(null);
    const suggestions2Ref = useRef<HTMLDivElement>(null);

    const generateEchartsOptions = useCallback((data: any, params = displayedParams, themeOverride?: string) => {
        if (!data || data.x.length === 0) return;
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
            
            series.push({ name: 'Bubble Point', type: 'line', data: sortedX.map((x: any, i: number) => [x, tC![i]]), symbol: 'none', color: 'green', lineStyle: { width: 2.5 } });
            series.push({ name: 'Dew Point',    type: 'line', data: dewPointData.map((p: {y: number, t: number}) => [p.y, p.t]), symbol: 'none', color: '#3b82f6', lineStyle: { width: 2.5 } });

        } else if (diagramType === 'pxy' && data.p) {
            yAxisName = "Pressure (bar)";
            titleConditionText = `at ${formatNumberToPrecision(params.temp ?? 0)} °C`;
            pressBar = sortedData.map((d: { p: number; }) => d.p/1e5);

            const dewPointData = sortedData.map((d: {y: number}, i: number) => ({ y: d.y, p: pressBar![i] })).sort((a: {y: number}, b: {y: number}) => a.y - b.y);
            dewPArray = dewPointData.map((p: { p: number; }) => p.p); // for tooltip interpolation

            series.push({ name: 'Bubble Point', type: 'line', data: sortedX.map((x: any,i: number)=>[x, pressBar![i]]), symbol: 'none', color: 'green', lineStyle: { width: 2.5 } });
            series.push({ name: 'Dew Point',    type: 'line', data: dewPointData.map((p: {y: number, p: number}) => [p.y, p.p]), symbol: 'none', color: '#3b82f6', lineStyle: { width: 2.5 } });
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
        
        const themeToUse = themeOverride ?? resolvedTheme;
        const isDark = themeToUse === 'dark';
        const textColor = isDark ? 'white' : '#000000';
        const tooltipBg = isDark ? '#08306b' : '#ffffff';
        const tooltipBorder = isDark ? '#55aaff' : '#333333';
        
        const titleText = `${comp1Label}-${comp2Label} ${diagramType.toUpperCase()} Diagram ${titleConditionText}`;
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
            },
            yAxis: {
                type: 'value', name: yAxisName,
                min: yMin,
                max: yMax,
                nameLocation: 'middle', nameGap: 60, nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
                axisLine: { lineStyle: { color: textColor } }, axisTick: { lineStyle: { color: textColor } },
                axisLabel: { 
                    color: textColor, 
                    fontSize: 14, 
                    fontFamily: 'Merriweather Sans',
                    formatter: (value: number) => {
                        // Hide the min and max values, show others
                        if ((yMin !== undefined && Math.abs(value - yMin) < 0.001) || 
                            (yMax !== undefined && Math.abs(value - yMax) < 0.001)) {
                            return '';
                        }
                        return value.toString();
                    }
                }, 
                splitLine: { show: false },
            },
            legend: {
                orient: 'horizontal', bottom: 10, left: 'center',
                data: series.map(s => s.name).filter((n): n is string => !!n),
                textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
            },
            tooltip: {
                trigger: 'axis', backgroundColor: tooltipBg, borderColor: tooltipBorder, borderWidth: 1,
                textStyle: { color: textColor, fontFamily: 'Merriweather Sans' },
                axisPointer: {
                    type: 'cross',
                    label: {
                        backgroundColor: tooltipBg, color: textColor, borderColor: tooltipBorder,
                        formatter: (params: any) => `${params.axisDimension === 'x' ? 'x/y' : (diagramType==='txy' ? 'T' : 'P')}: ${formatNumberToPrecision(params.value, 3)}`
                    },
                },
                formatter: (params: any) => {
                    if (!params || !Array.isArray(params) || params.length === 0) return '';
                    const xVal = params[0].axisValue ?? params[0].value[0];
                    let html = `x/y: ${formatNumberToPrecision(xVal, 3)}<br/>`;

                    if (diagramType === 'txy') {
                        const unit='°C';
                        const bubbleT = tC ? interp(xVal, sortedX, tC) : null;
                        const dewYPoints = sortedData.map((d: {y:number}) => d.y).sort((a: number,b: number)=>a-b);
                        const dewTPoints = sortedData.map((d: {t:number}) => d.t - 273.15).sort((a:number,b:number)=> a - (interp(b,dewYPoints, sortedData.map((d: {t:number})=>d.t-273.15))??0) );
                        const dewT = dewTArray ? interp(xVal, dewYPoints, dewTArray) : null;
                        
                        if(bubbleT!==null) html += `<span style=\"color: green;\">Bubble T: ${formatNumberToPrecision(bubbleT,3)} ${unit}</span><br/>`;
                        if(dewT!==null)    html += `<span style=\"color: #3b82f6;\">Dew T: ${formatNumberToPrecision(dewT,3)} ${unit}</span><br/>`;
                    } else if (diagramType === 'pxy') {
                        const unit='bar';
                        const bubbleP = pressBar ? interp(xVal, sortedX, pressBar) : null;
                        const dewYPoints = sortedData.map((d: {y: number}) => d.y).sort((a: number,b: number)=>a-b);
                        const dewP = pressBar && dewPArray ? interp(xVal, dewYPoints, dewPArray) : null;
                        if(bubbleP!==null) html += `<span style=\"color: green;\">Bubble P: ${formatNumberToPrecision(bubbleP,3)} ${unit}</span><br/>`;
                        if(dewP!==null)    html += `<span style=\"color: #3b82f6;\">Dew P: ${formatNumberToPrecision(dewP,3)} ${unit}</span><br/>`;
                    } else { // xy
                        html = `x: ${formatNumberToPrecision(xVal, 3)}<br/>`;
                        const yVal = interp(xVal, sortedX, sortedY);
                        if (yVal !== null) {
                            html += `<span style=\"color: red;\">y: ${formatNumberToPrecision(yVal,3)}</span><br/>`;
                        }
                    }
                    return html;
                }
            },
            series: series,
            animation: false,
            toolbox: {
                show: true, orient: 'vertical', right: 0, top: 'bottom',
                feature: { saveAsImage: { show: true, title: 'Save as Image', name: `phase-diagram-${comp1Label}-${comp2Label}`, backgroundColor: isDark ? '#08306b' : '#ffffff', pixelRatio: 2 } },
                iconStyle: { borderColor: textColor }
            },
        });
    }, [diagramType, useTemperatureForXY, displayedParams, resolvedTheme]);

    async function fetchCompoundDataLocal(compoundName: string): Promise<CompoundData | null> {
        if (!supabase) { throw new Error("Supabase client not initialized."); }
        console.log(`PhaseDiagrams: Fetching data for ${compoundName}...`);
    
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
    
            for (const source of sourcesToTry) {
                const { data: propsData, error: propsError } = await supabase
                    .from('compound_properties')
                    .select('properties')
                    .eq('compound_id', compoundId)
                    .eq('source', source)
                    .single();
                if (!propsError && propsData) { properties = propsData.properties; break; }
            }
            if (!properties) {
                const { data: anyPropsData } = await supabase.from('compound_properties').select('properties').eq('compound_id', compoundId).limit(1);
                if (anyPropsData && anyPropsData.length > 0) properties = anyPropsData[0].properties;
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
    
            let uniquacParams: UniquacPureComponentParams | null = null;
            const rPropObj = properties["UNIQUAC r"] || properties["Van der Waals volume"];
            const qPropObj = properties["UNIQUAC q"] || properties["Van der Waals area"];
            if (rPropObj && qPropObj) {
                const r_val = parseFloat(rPropObj.value ?? rPropObj);
                const q_val = parseFloat(qPropObj.value ?? qPropObj);
                if (!isNaN(r_val) && !isNaN(q_val)) uniquacParams = { r: r_val, q: q_val };
            }
            
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
                }
                if (vL_m3mol !== undefined && !isNaN(vL_m3mol) && vL_m3mol > 0) wilsonParams = { V_L_m3mol: vL_m3mol };
            }
    
            return { name: foundName, antoine, unifacGroups, cas_number: casNumber, prParams, srkParams, uniquacParams, wilsonParams };
        } catch (err: any) {
            console.error(`Error fetching data for ${compoundName}:`, err.message);
            setError(`Data fetch failed for ${compoundName}: ${err.message}`);
            return null;
        }
    }

    const antoineBoilingPointSolverLocal = (antoineParams: AntoineParams | null, P_target: number): number | null => {
        if (!antoineParams) return null;
        try {
            let logP: number;
            const P_converted_to_antoine_units = P_target / (antoineParams.Units?.toLowerCase() === 'kpa' ? 1000 : 1);
            if (antoineParams.EquationNo === 1 || antoineParams.EquationNo === '1') {
                logP = Math.log10(P_converted_to_antoine_units);
            } else {
                logP = Math.log(P_converted_to_antoine_units);
            }
            if (antoineParams.A - logP === 0) return null;
            const T_K = antoineParams.B / (antoineParams.A - logP) - antoineParams.C;
            return (T_K > 0 && T_K < 1000) ? T_K : null;
        } catch { return null; }
    };
    
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
                const [fetched1, fetched2] = await Promise.all([fetchCompoundDataLocal(comp1Name), fetchCompoundDataLocal(comp2Name)]);
                if (!fetched1 || !fetched2) { if (showLoading) setLoading(false); return; }
                
                let fetchedParams: any;
                if (fluidPackage === 'unifac') {
                    if (!fetched1.unifacGroups || !fetched2.unifacGroups) throw new Error("UNIFAC groups missing.");
                    const allSubgroupIds = new Set<number>();
                    [fetched1, fetched2].forEach(comp => { if (comp.unifacGroups) Object.keys(comp.unifacGroups).forEach(id => allSubgroupIds.add(parseInt(id))); });
                    fetchedParams = await fetchUnifacInteractionParams(supabase, Array.from(allSubgroupIds));
                } else if (fluidPackage === 'nrtl') {
                    fetchedParams = await fetchNrtlParameters(supabase, fetched1.cas_number!, fetched2.cas_number!);
                } else if (fluidPackage === 'pr') {
                    fetchedParams = await fetchPrInteractionParams(supabase, fetched1.cas_number!, fetched2.cas_number!);
                } else if (fluidPackage === 'srk') {
                    fetchedParams = await fetchSrkInteractionParams(supabase, fetched1.cas_number!, fetched2.cas_number!);
                } else if (fluidPackage === 'uniquac') {
                    fetchedParams = await fetchUniquacInteractionParams(supabase, fetched1.cas_number!, fetched2.cas_number!);
                } else if (fluidPackage === 'wilson') {
                    fetchedParams = await fetchWilsonInteractionParams(supabase, fetched1.cas_number!, fetched2.cas_number!);
                }
                
                setComp1Data(fetched1); setComp2Data(fetched2); setInteractionParams(fetchedParams);
                d1 = fetched1; d2 = fetched2; params = fetchedParams;
            }
            if (!d1 || !d2 || !params) { throw new Error("Component data or parameters not available."); }
            const components: [CompoundData, CompoundData] = [d1, d2];
            
            const pointsCount = 51;
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
                        const fixedPressPa = useTemperatureForXY === false && diagramType === 'xy' ? pressureBar! * 1e5 : null;
                        
                        if (fixedTempK) { // Px-y or T-const x-y
                            const initialPressureGuess = (libCalculatePsat_Pa(d1.antoine!, fixedTempK) + libCalculatePsat_Pa(d2.antoine!, fixedTempK)) / 2 || 101325;
                            if(fluidPackage==='unifac') resultPoint=calculateBubblePressureUnifac(components,x1,fixedTempK,params);
                            else if(fluidPackage==='nrtl') resultPoint=calculateBubblePressureNrtl(components,x1,fixedTempK,params);
                            else if(fluidPackage==='pr') resultPoint=calculateBubblePressurePr(components,x1,fixedTempK,params,initialPressureGuess);
                            else if(fluidPackage==='srk') resultPoint=calculateBubblePressureSrk(components,x1,fixedTempK,params,initialPressureGuess);
                            else if(fluidPackage==='uniquac') resultPoint=calculateBubblePressureUniquac(components,x1,fixedTempK,params);
                            else if(fluidPackage==='wilson') resultPoint=calculateBubblePressureWilson(components,x1,fixedTempK,params);
                        } else if (fixedPressPa) { // P-const x-y
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
                
                if (resultPoint && resultPoint.error === undefined && typeof resultPoint.comp1_equilibrium === 'number') {
                    results.x.push(x1);
                    results.y.push(resultPoint.comp1_equilibrium);
                    if (diagramType === 'txy' && results.t && resultPoint.T_K) results.t.push(resultPoint.T_K);
                    else if (diagramType === 'pxy' && results.p && resultPoint.P_Pa) results.p.push(resultPoint.P_Pa);
                }
            }
            
            setChartData(results);
            const newParams = {
                comp1: d1.name, comp2: d2.name,
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
    }, [comp1Name, comp2Name, diagramType, fluidPackage, pressureBar, temperatureC, useTemperatureForXY, resolvedTheme, comp1Data, comp2Data, interactionParams, displayedParams.package]);

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
                                <Input type="text" value={pressureMax} onChange={(e) => setPressureMax(e.target.value)} className="w-16 h-8 text-xs"/>
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
                                <Input type="text" value={tempMax} onChange={(e) => setTempMax(e.target.value)} className="w-16 h-8 text-xs"/>
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
                                <Input type="text" value={useTemperatureForXY ? tempMax : pressureMax} onChange={e => useTemperatureForXY ? setTempMax(e.target.value) : setPressureMax(e.target.value)} className="w-16 h-8 text-xs"/>
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
            const { data, error } = await supabase.from('compounds').select('name').ilike('name', `${inputValue}%`).limit(5);
            if (error) { throw error; }
            const suggestions = data ? data.map(item => item.name) : [];
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
        setActiveSuggestionInput('comp1');
        fetchSuggestions(newValue, 'comp1');
    };

    const handleComp2NameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const newValue = e.target.value;
        setComp2Name(newValue);
        setActiveSuggestionInput('comp2');
        fetchSuggestions(newValue, 'comp2');
    };

    const handleSuggestionClick = (suggestion: string, inputTarget: 'comp1' | 'comp2') => {
        setAutoGenerateOnCompChange(true);
        if (inputTarget === 'comp1') {
            setComp1Name(suggestion);
            setShowComp1Suggestions(false);
            input1Ref.current?.focus();
        } else {
            setComp2Name(suggestion);
            setShowComp2Suggestions(false);
            input2Ref.current?.focus();
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
        return () => document.removeEventListener("mousedown", handleClickOutside);
    }, [activeSuggestionInput]);

    return (
        <div className="container mx-auto p-4 md:p-8">
            <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                <div className="lg:col-span-1 space-y-6">
                    <Card>
                        <CardHeader><CardTitle>VLE Diagram Generator</CardTitle></CardHeader>
                        <CardContent className="space-y-6">
                            <Tabs value={diagramType} onValueChange={(v) => setDiagramType(v as DiagramType)} className="w-full">
                                <TabsList className="grid w-full grid-cols-3"><TabsTrigger value="txy">T-x-y</TabsTrigger><TabsTrigger value="pxy">P-x-y</TabsTrigger><TabsTrigger value="xy">x-y</TabsTrigger></TabsList>
                            </Tabs>
                            <div className="space-y-4">
                                <div className="flex items-center gap-2">
                                    <div className="relative flex-1">
                                        <Input
                                            ref={input1Ref} id="comp1Name" value={comp1Name} onChange={handleComp1NameChange} onKeyDown={handleKeyDown}
                                            onFocus={() => { setActiveSuggestionInput('comp1'); fetchSuggestions(comp1Name, 'comp1'); }}
                                            placeholder="Methanol" required className="w-full" autoComplete="off"
                                        />
                                        {showComp1Suggestions && comp1Suggestions.length > 0 && (
                                            <div ref={suggestions1Ref} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1">
                                                {comp1Suggestions.map((suggestion, index) => (
                                                    <div key={index} onClick={() => handleSuggestionClick(suggestion, 'comp1')} className="px-3 py-2 hover:bg-accent cursor-pointer text-sm">
                                                        {suggestion}
                                                    </div>
                                                ))}
                                            </div>
                                        )}
                                    </div>
                                    <Button variant="ghost" size="icon" onClick={() => {
                                        const temp = comp1Name;
                                        setComp1Name(comp2Name);
                                        setComp2Name(temp);
                                        setAutoGenerateOnCompChange(true);
                                        generateDiagram();
                                    }} title="Swap Components"><ArrowLeftRight className="h-4 w-4" /></Button>
                                    <div className="relative flex-1">
                                        <Input
                                            ref={input2Ref} id="comp2Name" value={comp2Name} onChange={handleComp2NameChange} onKeyDown={handleKeyDown}
                                            onFocus={() => { setActiveSuggestionInput('comp2'); fetchSuggestions(comp2Name, 'comp2'); }}
                                            placeholder="Water" required className="w-full" autoComplete="off"
                                        />
                                        {showComp2Suggestions && comp2Suggestions.length > 0 && (
                                            <div ref={suggestions2Ref} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1">
                                                {comp2Suggestions.map((suggestion, index) => (
                                                    <div key={index} onClick={() => handleSuggestionClick(suggestion, 'comp2')} className="px-3 py-2 hover:bg-accent cursor-pointer text-sm">
                                                        {suggestion}
                                                    </div>
                                                ))}
                                            </div>
                                        )}
                                    </div>
                                </div>
                                {renderConditionalInput()}
                                <div className="flex items-center gap-2">
                                    <Label htmlFor="fluidPackage" className="text-sm font-medium whitespace-nowrap">Fluid Package:</Label>
                                    <Select value={fluidPackage} onValueChange={(v) => {
                                        setFluidPackage(v as FluidPackageType);
                                        // The generateDiagram() call is removed from here.
                                        // The useEffect hook will handle the regeneration
                                        // when the fluidPackage state changes.
                                    }}><SelectTrigger id="fluidPackage"><SelectValue /></SelectTrigger><SelectContent><SelectItem value="uniquac">UNIQUAC</SelectItem><SelectItem value="pr">Peng-Robinson</SelectItem><SelectItem value="wilson">Wilson</SelectItem><SelectItem value="nrtl">NRTL</SelectItem><SelectItem value="srk">SRK</SelectItem><SelectItem value="unifac">UNIFAC</SelectItem></SelectContent></Select>
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
                            <div className="p-4 bg-muted rounded-md">
                              <p className="text-sm font-medium">{comp1Name.charAt(0).toUpperCase() + comp1Name.slice(1)}</p>
                              <p className="text-2xl font-bold">{formatNumberToPrecision(chartData.t[chartData.x.indexOf(1)] - 273.15, 3)} °C</p>
                            </div>
                            <div className="p-4 bg-muted rounded-md">
                              <p className="text-sm font-medium">{comp2Name.charAt(0).toUpperCase() + comp2Name.slice(1)}</p>
                              <p className="text-2xl font-bold">{formatNumberToPrecision(chartData.t[chartData.x.indexOf(0)] - 273.15, 3)} °C</p>
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
                            <div className="p-4 bg-muted rounded-md">
                              <p className="text-sm font-medium">{comp1Name.charAt(0).toUpperCase() + comp1Name.slice(1)}</p>
                              <p className="text-2xl font-bold">{formatNumberToPrecision(chartData.p[chartData.x.indexOf(1)] / 1e5, 3)} bar</p>
                            </div>
                            <div className="p-4 bg-muted rounded-md">
                              <p className="text-sm font-medium">{comp2Name.charAt(0).toUpperCase() + comp2Name.slice(1)}</p>
                              <p className="text-2xl font-bold">{formatNumberToPrecision(chartData.p[chartData.x.indexOf(0)] / 1e5, 3)} bar</p>
                            </div>
                          </div>
                        </CardContent>
                      </Card>
                    )}
                </div>
                <div className="lg:col-span-2">
                    <Card className="h-full"><CardContent className="py-2 h-full">
                        <div className="relative aspect-square rounded-md h-full">
                           {loading && ( <div className="absolute inset-0 flex items-center justify-center text-muted-foreground"><div className="text-center"><div className="mb-2">Loading & Calculating...</div><div className="text-sm text-muted-foreground/70">Using { fluidPackage.toUpperCase()} model.</div></div></div> )}
                           {!loading && !chartData && !error && ( <div className="absolute inset-0 flex items-center justify-center text-muted-foreground">Provide inputs and generate a diagram.</div> )}
                           {error && !loading && ( <div className="absolute inset-0 flex items-center justify-center text-red-400">Error: {error}</div> )}
                           {!loading && chartData && Object.keys(echartsOptions).length > 0 && (
                            <ReactECharts echarts={echarts} option={echartsOptions} style={{ height: '100%', width: '100%' }} notMerge={false} lazyUpdate={true} />
                           )}
                        </div>
                    </CardContent></Card>
                </div>
            </div>
        </div>
    );
}