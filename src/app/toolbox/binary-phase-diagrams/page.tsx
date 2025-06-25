'use client';

import React, { useState, useEffect, useCallback, useRef } from 'react';
import { createClient, SupabaseClient } from '@supabase/supabase-js';
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

import {
    calculateVleDiagramData,
    fetchCompoundData,
    type VleChartData,
    type FluidPackageType,
    type DiagramType
} from '@/lib/phase-diagram-calculations';

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

export default function VleDiagramPage() {
    const { resolvedTheme } = useTheme();
    const [comp1Name, setComp1Name] = useState('methanol');
    const [comp2Name, setComp2Name] = useState('water');
    const [diagramType, setDiagramType] = useState<DiagramType>('txy');
    const [temperatureC, setTemperatureC] = useState<number | null>(60);
    const [pressureBar, setPressureBar] = useState<number | null>(1);
    const [fluidPackage, setFluidPackage] = useState<FluidPackageType>('unifac');
    const [temperatureInput, setTemperatureInput] = useState<string>(String(temperatureC));
    const [pressureInput, setPressureInput] = useState<string>(String(pressureBar));
    const [useTemperatureForXY, setUseTemperatureForXY] = useState(true);
    const [chartData, setChartData] = useState<VleChartData | null>(null);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
    const [displayedParams, setDisplayedParams] = useState({ comp1: '', comp2: '', temp: null as number | null, pressure: null as number | null, package: '', type: '' });
    const [displayedPressure, setDisplayedPressure] = useState<number | null>(null);
    const [displayedUseTemp, setDisplayedUseTemp] = useState(true);
    const [displayedFluidPackage, setDisplayedFluidPackage] = useState<FluidPackageType>('unifac');
    const [chartTitle, setChartTitle] = useState('Phase Diagram');

    const [comp1Suggestions, setComp1Suggestions] = useState<string[]>([]);
    const [comp2Suggestions, setComp2Suggestions] = useState<string[]>([]);
    const [showComp1Suggestions, setShowComp1Suggestions] = useState(false);
    const [showComp2Suggestions, setShowComp2Suggestions] = useState(false);
    const [activeSuggestionInput, setActiveSuggestionInput] = useState<'comp1' | 'comp2' | null>(null);
    const [autoGenerateOnCompChange, setAutoGenerateOnCompChange] = useState(false);

    const input1Ref = useRef<HTMLInputElement>(null);
    const input2Ref = useRef<HTMLInputElement>(null);
    const suggestions1Ref = useRef<HTMLDivElement>(null);
    const suggestions2Ref = useRef<HTMLDivElement>(null);

    const generateDiagram = useCallback(async () => {
        setLoading(true); setError(null); setChartData(null);
        if (!comp1Name || !comp2Name) {
            setError("Please enter both compound names.");
            setLoading(false);
            return;
        }
        try {
            const [data1, data2] = await Promise.all([fetchCompoundData(comp1Name), fetchCompoundData(comp2Name)]);
            
            let fixedConditionValue: number;
            let isTempFixed: boolean;

            if (diagramType === 'txy') {
                if (pressureBar === null) throw new Error("Pressure must be set.");
                fixedConditionValue = pressureBar; isTempFixed = false;
            } else if (diagramType === 'pxy') {
                if (temperatureC === null) throw new Error("Temperature must be set.");
                fixedConditionValue = temperatureC; isTempFixed = true;
            } else { // xy
                isTempFixed = useTemperatureForXY;
                const value = isTempFixed ? temperatureC : pressureBar;
                if (value === null) throw new Error("A fixed condition (T or P) must be set for x-y diagrams.");
                fixedConditionValue = value;
            }
            
            const resultData = await calculateVleDiagramData([data1, data2], diagramType, fluidPackage, fixedConditionValue, isTempFixed);
            setChartData(resultData);
            const newParams = {
                comp1: data1.name,
                comp2: data2.name,
                temp: isTempFixed ? fixedConditionValue : null,
                pressure: !isTempFixed ? fixedConditionValue : null,
                package: fluidPackage,
                type: diagramType,
                temperatureInput,
                pressureInput
            };
            setDisplayedParams(newParams);
            generateEchartsOptions(resultData, newParams, resolvedTheme);
        } catch (err: any) {
            setError(`Calculation failed: ${err.message}`);
        } finally {
            setLoading(false);
        }
    }, [comp1Name, comp2Name, diagramType, fluidPackage, pressureBar, temperatureC, useTemperatureForXY, resolvedTheme, temperatureInput, pressureInput]);

    useEffect(() => {
        if (supabase) {
            generateDiagram();
        }
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [supabase]);

    // Regenerate chart options on theme change
    useEffect(() => {
        if (chartData) {
            generateEchartsOptions(chartData, displayedParams, resolvedTheme);
        }
    }, [resolvedTheme]);

    // Auto-generate when diagram type or fluid package changes
    useEffect(() => {
        if (!loading) generateDiagram();
    }, [diagramType, fluidPackage]);

    const generateEchartsOptions = useCallback((data: VleChartData, params = displayedParams, themeOverride?: string) => {
        if (!data || data.x.length === 0) return;
        const series: SeriesOption[] = [];
        let yAxisName = "Mole Fraction", titleConditionText = "";
        const capitalize = (s: string) => s ? s.charAt(0).toUpperCase() + s.slice(1) : '';
        const comp1Label = capitalize(params.comp1);
        const comp2Label = capitalize(params.comp2);

        const interp = (xTarget:number, xArr:number[], yArr:number[]): number|null => {
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

        if (diagramType === 'txy' && data.t) {
            yAxisName = "Temperature (°C)";
            titleConditionText = `at ${(params.pressure ?? 0).toPrecision(4)} bar`;
            tC = data.t!.map(v => v - 273.15);
            dewTArray = data.x.map(xVal => interp(xVal, data.y, tC!) ?? NaN);
            series.push({ name: 'Bubble Point', type: 'line', data: data.x.map((x, i) => [x, tC![i]]), symbol: 'none', color: 'green', lineStyle: { width: 2.5 } });
            series.push({ name: 'Dew Point',    type: 'line', data: data.x.map((x, i) => [x, dewTArray![i]]).filter(p=>!isNaN(p[1])), symbol: 'none', color: '#3b82f6', lineStyle: { width: 2.5 } });
        } else if (diagramType === 'pxy' && data.p) {
            yAxisName = "Pressure (bar)";
            titleConditionText = `at ${(params.temp ?? 0).toFixed(1)} °C`;
            pressBar = data.p!.map(v => v/1e5);
            dewPArray = data.x.map(xVal => interp(xVal, data.y, pressBar!) ?? NaN);
            series.push({ name: 'Bubble Point', type: 'line', data: data.x.map((x,i)=>[x, pressBar![i]]), symbol: 'none', color: 'green', lineStyle: { width: 2.5 } });
            series.push({ name: 'Dew Point',    type: 'line', data: data.x.map((x,i)=>[x, dewPArray![i]]).filter(p=>!isNaN(p[1])), symbol: 'none', color: '#3b82f6', lineStyle: { width: 2.5 } });
        } else {
            yAxisName = `Vapor Mole Fraction ${comp1Label} (y)`;
            titleConditionText = useTemperatureForXY ? `at ${(params.temp ?? 0).toFixed(1)} °C` : `at ${(params.pressure ?? 0).toPrecision(4)} bar`;
            series.push({ name: 'Equilibrium', type: 'line', data: data.x.map((x: number, i: number) => [x, data.y[i]]), symbol: 'none', color: 'red', lineStyle: { width: 2.5 } });
            series.push({ name: 'y=x', type: 'line', data: [[0, 0], [1, 1]], symbol: 'none', lineStyle: { width: 1.5, type: 'dotted' }, color: 'cyan' });
        }
        
        const themeToUse = themeOverride ?? resolvedTheme;
        const isDark = themeToUse === 'dark';
        const textColor = isDark ? 'white' : '#000000';
        const tooltipBg = isDark ? '#08306b' : '#ffffff';
        const tooltipBorder = isDark ? '#55aaff' : '#333333';
        let inputTempString = (params as any).temperatureInput ?? temperatureInput;
        let inputPressureString = (params as any).pressureInput ?? pressureInput;
        let titleConditionTextRaw = '';
        if (diagramType === 'txy') {
            titleConditionTextRaw = pressureInput ? `at ${inputPressureString} bar` : '';
        } else if (diagramType === 'pxy') {
            titleConditionTextRaw = temperatureInput ? `at ${inputTempString} °C` : '';
        } else if (diagramType === 'xy') {
            if (useTemperatureForXY) {
                titleConditionTextRaw = temperatureInput ? `at ${inputTempString} °C` : '';
            } else {
                titleConditionTextRaw = pressureInput ? `at ${inputPressureString} bar` : '';
        }
        }
        const titleText = `${comp1Label}-${comp2Label} ${diagramType.toUpperCase()} Diagram ${titleConditionTextRaw}`;
        const xAxisName = diagramType === 'xy' ? `Liquid Mole Fraction ${comp1Label} (x)` : `Mole Fraction ${comp1Label} (x/y)`;

        setEchartsOptions({
            backgroundColor: 'transparent',
            title: { text: titleText, left: 'center', textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } },
            grid: {
                top: '5%',
                bottom: '10%',
                left: '5%',
                right: '5%',
                containLabel: true,
            },
            xAxis: {
                type: 'value', min: 0, max: 1, interval: 0.1,
                name: xAxisName,
                nameLocation: 'middle', nameGap: 40, nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
                axisLine: { lineStyle: { color: textColor } }, axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
                axisLabel: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans', formatter: '{value}' }, splitLine: { show: false },
            },
            yAxis: {
                type: 'value',
                name: yAxisName,
                nameLocation: 'middle', nameGap: 40, nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
                axisLine: { lineStyle: { color: textColor } }, axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
                axisLabel: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans', formatter: '{value}' }, splitLine: { show: false },
            },
            legend: {
                orient: 'horizontal',
                bottom: 10,
                left: 'center',
                data: series.map(s => s.name).filter((n): n is string => !!n),
                textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
                itemWidth: 25, itemHeight: 2,
            },
            tooltip: {
                trigger: 'axis',
                backgroundColor: tooltipBg,
                borderColor: tooltipBorder,
                borderWidth: 1,
                textStyle: { color: textColor, fontFamily: 'Merriweather Sans' },
                axisPointer: {
                    type: 'cross',
                    label: {
                        show: true,
                        backgroundColor: tooltipBg,
                        color: textColor,
                        borderColor: tooltipBorder,
                        borderWidth: 1,
                        shadowBlur: 0,
                        shadowColor: 'transparent',
                        fontFamily: 'Merriweather Sans',
                        formatter: (params: any) => {
                            if (params.axisDimension === 'x') {
                                return (diagramType === 'txy' || diagramType === 'pxy')
                                    ? `x/y: ${formatNumberToPrecision(params.value, 3)}`
                                    : `x: ${formatNumberToPrecision(params.value, 3)}`;
                            } else if (params.axisDimension === 'y') {
                                if (diagramType === 'txy') {
                                    return `T: ${formatNumberToPrecision(params.value, 3)} °C`;
                                } else if (diagramType === 'pxy') {
                                    return `P: ${formatNumberToPrecision(params.value, 3)} bar`;
                                } else {
                                    return `y: ${formatNumberToPrecision(params.value, 3)}`;
                                }
                            }
                            return `${formatNumberToPrecision(params.value, 3)}`;
                        }
                    },
                    crossStyle: { color: isDark ? '#999' : '#666' }
                },
                formatter: (params: any) => {
                    if (!params || !Array.isArray(params) || params.length === 0) return '';
                    const xVal = params[0].axisValue ?? params[0].value[0];
                    let html = (diagramType === 'txy' || diagramType === 'pxy')
                        ? `x/y: ${formatNumberToPrecision(xVal, 3)}<br/>`
                        : `x: ${formatNumberToPrecision(xVal, 3)}<br/>`;

                    const interpY = (xTarget:number, xArr:number[], yArr:number[]): number|null => {
                        if (xTarget<=xArr[0]) return yArr[0];
                        if (xTarget>=xArr[xArr.length-1]) return yArr[yArr.length-1];
                        for(let i=0;i<xArr.length-1;i++){
                            if((xArr[i]<=xTarget && xTarget<=xArr[i+1])||(xArr[i]>=xTarget && xTarget>=xArr[i+1])){
                                const x1=xArr[i],x2=xArr[i+1],y1=yArr[i],y2=yArr[i+1];
                                if(Math.abs(x2-x1)<1e-9) return y1;
                                return y1+(y2-y1)*(xTarget-x1)/(x2-x1);
                            }
                        }
                        return null;
                    };

                    if (diagramType === 'txy') {
                        const unit='°C';
                        const bubbleT = tC ? interp(xVal, data.x, tC) : null;
                        const dewT = tC && dewTArray ? interp(xVal, data.x, dewTArray) : null;
                        if(bubbleT!==null) html += `<span style=\"color: green;\">Bubble: ${formatNumberToPrecision(bubbleT,3)} ${unit}</span><br/>`;
                        if(dewT!==null)    html += `<span style=\"color: #3b82f6;\">Dew: ${formatNumberToPrecision(dewT,3)} ${unit}</span><br/>`;
                    } else if (diagramType === 'pxy') {
                        const unit='bar';
                        const bubbleP = pressBar ? interp(xVal, data.x, pressBar) : null;
                        const dewP = pressBar && dewPArray ? interp(xVal, data.x, dewPArray) : null;
                        if(bubbleP!==null) html += `<span style=\"color: green;\">Bubble: ${formatNumberToPrecision(bubbleP,3)} ${unit}</span><br/>`;
                        if(dewP!==null)    html += `<span style=\"color: #3b82f6;\">Dew: ${formatNumberToPrecision(dewP,3)} ${unit}</span><br/>`;
                    } else {
                        (params as any[]).forEach((p: any) => {
                            if (p.seriesName === 'Equilibrium') {
                                html += `<span style=\"color: ${p.color};\">y: ${formatNumberToPrecision(p.value[1],3)}</span><br/>`;
                            }
                        });
                    }
                    return html;
                }
            },
            series: series,
            animationDuration: 300, animationEasing: 'cubicInOut',
            toolbox: {
                show: true, orient: 'vertical', right: 0, top: 'bottom',
                feature: { saveAsImage: { show: true, title: 'Save as Image', name: `phase-diagram-${comp1Label}-${comp2Label}`, backgroundColor: isDark ? '#08306b' : '#ffffff', pixelRatio: 2 } },
                iconStyle: { borderColor: textColor }
            },
        });
    }, [diagramType, useTemperatureForXY, displayedParams, resolvedTheme]);

    const formatNumberToPrecision = (num: any, precision: number = 3): string => {
        if (typeof num === 'number') {
          if (num === 0) return '0';
          const fixed = num.toPrecision(precision);
          if (fixed.includes('.')) {
            return parseFloat(fixed).toString(); 
          }
          return fixed;
        }
        return String(num);
    };
    
    const handleKeyDown = (event: React.KeyboardEvent<HTMLInputElement>) => {
        if (event.key === 'Enter') {
            event.preventDefault();
            generateDiagram();
        }
    };

    const renderConditionalInput = () => {
        switch (diagramType) {
            case 'txy': return <div className="flex items-center gap-2"><Label htmlFor="pressure" className="whitespace-nowrap">Pressure:</Label><Input id="pressure" type="text" value={pressureInput} onChange={e => { setPressureInput(e.target.value); setPressureBar(parseFloat(e.target.value) || null); }} onKeyDown={handleKeyDown} placeholder="e.g., 1.0" className="flex-1" /><span className="text-sm text-muted-foreground w-8 text-left">bar</span></div>;
            case 'pxy': return <div className="flex items-center gap-2"><Label htmlFor="temperature" className="whitespace-nowrap">Temperature:</Label><Input id="temperature" type="text" value={temperatureInput} onChange={e => { setTemperatureInput(e.target.value); setTemperatureC(parseFloat(e.target.value) || null); }} onKeyDown={handleKeyDown} placeholder="e.g., 60" className="flex-1" /><span className="text-sm text-muted-foreground w-8 text-left">°C</span></div>;
            case 'xy': return <div className="flex items-center gap-2"><Tabs value={useTemperatureForXY ? "temperature" : "pressure"} onValueChange={(v) => setUseTemperatureForXY(v === "temperature")} className="flex-shrink-0"><TabsList><TabsTrigger value="temperature">Constant T</TabsTrigger><TabsTrigger value="pressure">Constant P</TabsTrigger></TabsList></Tabs>{useTemperatureForXY ? <Input id="temperature-xy" type="text" value={temperatureInput} onChange={e => { setTemperatureInput(e.target.value); setTemperatureC(parseFloat(e.target.value) || null); }} onKeyDown={handleKeyDown} placeholder="Temp" className="flex-1" /> : <Input id="pressure-xy" type="text" value={pressureInput} onChange={e => { setPressureInput(e.target.value); setPressureBar(parseFloat(e.target.value) || null); }} onKeyDown={handleKeyDown} placeholder="Pressure" className="flex-1" />}<span className="text-sm text-muted-foreground w-8 text-left">{useTemperatureForXY ? "°C" : "bar"}</span></div>;
        }
    };

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
                .from('compounds')
                .select('name')
                .ilike('name', `${inputValue}%`)
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
    }, []);

    const debouncedFetchSuggestions = useCallback(debounce(fetchSuggestions, 300), [fetchSuggestions]);

    const handleComp1NameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const newValue = e.target.value;
        setComp1Name(newValue);
        setActiveSuggestionInput('comp1');
        if (newValue.trim() === "") {
            setShowComp1Suggestions(false);
            setComp1Suggestions([]);
        } else {
            debouncedFetchSuggestions(newValue, 'comp1');
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
            debouncedFetchSuggestions(newValue, 'comp2');
        }
    };

    const handleSuggestionClick = (suggestion: string, inputTarget: 'comp1' | 'comp2') => {
        if (inputTarget === 'comp1') {
            setComp1Name(suggestion);
            setShowComp1Suggestions(false);
            setComp1Suggestions([]);
            input1Ref.current?.focus();
        } else {
            setComp2Name(suggestion);
            setShowComp2Suggestions(false);
            setComp2Suggestions([]);
            input2Ref.current?.focus();
        }
        setAutoGenerateOnCompChange(true);
    };

    useEffect(() => {
        function handleClickOutside(event: MouseEvent) {
            const target = event.target as Node;
            if (activeSuggestionInput === 'comp1') {
                if (suggestions1Ref.current && !suggestions1Ref.current.contains(target) && input1Ref.current && !input1Ref.current.contains(target)) {
                    setShowComp1Suggestions(false);
                }
            }
            if (activeSuggestionInput === 'comp2') {
                if (suggestions2Ref.current && !suggestions2Ref.current.contains(target) && input2Ref.current && !input2Ref.current.contains(target)) {
                    setShowComp2Suggestions(false);
                }
            }
        }
        document.addEventListener("mousedown", handleClickOutside);
        return () => {
            document.removeEventListener("mousedown", handleClickOutside);
        };
    }, [activeSuggestionInput]);

    // Auto-generate when comp names change due to suggestion or swap
    useEffect(() => {
        if (autoGenerateOnCompChange) {
            generateDiagram();
            setAutoGenerateOnCompChange(false);
        }
    }, [comp1Name, comp2Name]);

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
                                            ref={input1Ref}
                                            id="comp1Name"
                                            value={comp1Name}
                                            onChange={handleComp1NameChange}
                                            onKeyDown={handleKeyDown}
                                            onFocus={() => {
                                                setActiveSuggestionInput('comp1');
                                                if (comp1Name.trim() !== "" && comp1Suggestions.length > 0) setShowComp1Suggestions(true);
                                                else if (comp1Name.trim() !== "") debouncedFetchSuggestions(comp1Name, 'comp1');
                                            }}
                                            placeholder="Methanol"
                                            required
                                            className="w-full"
                                            autoComplete="off"
                                        />
                                        {showComp1Suggestions && comp1Suggestions.length > 0 && (
                                            <div ref={suggestions1Ref} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1 max-h-48 overflow-y-auto">
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
                                    <Button variant="ghost" size="icon" onClick={() => {
                                        const temp = comp1Name;
                                        setComp1Name(comp2Name);
                                        setComp2Name(temp);
                                        setAutoGenerateOnCompChange(true);
                                    }} title="Swap Components"><ArrowLeftRight className="h-4 w-4" /></Button>
                                    <div className="relative flex-1">
                                        <Input
                                            ref={input2Ref}
                                            id="comp2Name"
                                            value={comp2Name}
                                            onChange={handleComp2NameChange}
                                            onKeyDown={handleKeyDown}
                                            onFocus={() => {
                                                setActiveSuggestionInput('comp2');
                                                if (comp2Name.trim() !== "" && comp2Suggestions.length > 0) setShowComp2Suggestions(true);
                                                else if (comp2Name.trim() !== "") debouncedFetchSuggestions(comp2Name, 'comp2');
                                            }}
                                            placeholder="Water"
                                            required
                                            className="w-full"
                                            autoComplete="off"
                                        />
                                        {showComp2Suggestions && comp2Suggestions.length > 0 && (
                                            <div ref={suggestions2Ref} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1 max-h-48 overflow-y-auto">
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
                                {renderConditionalInput()}
                                <div className="flex items-center gap-2">
                                    <Label htmlFor="fluidPackage" className="text-sm font-medium whitespace-nowrap">Fluid Package:</Label>
                                    <Select value={fluidPackage} onValueChange={(v) => setFluidPackage(v as FluidPackageType)}><SelectTrigger id="fluidPackage"><SelectValue /></SelectTrigger><SelectContent><SelectItem value="uniquac">UNIQUAC</SelectItem><SelectItem value="pr">Peng-Robinson</SelectItem><SelectItem value="wilson">Wilson</SelectItem><SelectItem value="nrtl">NRTL</SelectItem><SelectItem value="srk">SRK</SelectItem><SelectItem value="unifac">UNIFAC</SelectItem></SelectContent></Select>
                                </div>
                            </div>
                            <Button onClick={generateDiagram} disabled={loading} className="w-full">{loading ? 'Calculating...' : 'Generate Diagram'}</Button>
                            {error && <p className="text-sm text-red-500 mt-2">{error}</p>}
                        </CardContent>
                    </Card>
                </div>
                <div className="lg:col-span-2">
                    <Card className="h-full"><CardContent className="py-2 h-full">
                        <div className="relative aspect-square rounded-md h-full">
                           {loading && ( <div className="absolute inset-0 flex items-center justify-center text-muted-foreground"><div className="text-center"><div className="mb-2">Loading & Calculating...</div><div className="text-sm text-muted-foreground/70">Using { fluidPackage.toUpperCase()} model.</div></div></div> )}
                           {!loading && !chartData && !error && ( <div className="absolute inset-0 flex items-center justify-center text-muted-foreground">Provide inputs and generate a diagram.</div> )}
                           {error && !loading && ( <div className="absolute inset-0 flex items-center justify-center text-red-400">Error: {error}</div> )}
                           {!loading && chartData && Object.keys(echartsOptions).length > 0 && (
                            <ReactECharts echarts={echarts} option={echartsOptions} style={{ height: '100%', width: '100%' }} notMerge={false} lazyUpdate={false} />
                           )}
                        </div>
                    </CardContent></Card>
                </div>
            </div>
        </div>
    );
}