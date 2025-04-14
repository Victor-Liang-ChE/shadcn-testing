"use client";
import React, { useState, useEffect, useCallback, useMemo, useRef } from 'react';

// ECharts imports
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsOption, SeriesOption } from 'echarts';
import { LineChart } from 'echarts/charts';
import { TitleComponent, TooltipComponent, GridComponent, LegendComponent, ToolboxComponent } from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

echarts.use([TitleComponent, TooltipComponent, GridComponent, LegendComponent, ToolboxComponent, LineChart, CanvasRenderer]);

// Shadcn UI imports
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Label } from "@/components/ui/label";
import { Slider } from "@/components/ui/slider";
import { Button } from "@/components/ui/button";
import { Checkbox } from "@/components/ui/checkbox";
import { Skeleton } from "@/components/ui/skeleton";

// --- Type Definitions ---
type OrderType = 'first' | 'second';
type FunctionType = 'step' | 'ramp';
interface ProcessParams {
    K: number; // Gain
    M: number; // Magnitude
    tau: number; // Time Constant
    zeta?: number; // Damping Ratio (only for second order)
}
interface Metrics {
    peakTime: number | null;
    overshoot: number | null;
    oscillationPeriod: number | null;
    decayRatio: number | null;
}
interface SimulationResult {
    t: number[];
    y: number[];
    y_input: number[];
    metrics: Metrics;
}

// --- Helper Functions ---
const formatNumber = (num: number | null | undefined, precision = 2): string => {
    if (num === null || num === undefined || isNaN(num)) return "---";
    return num.toFixed(precision);
};

// --- Main Component ---
export default function ProcessControlPage() {
    const [order, setOrder] = useState<OrderType>('first');
    const [functionType, setFunctionType] = useState<FunctionType>('step');
    const [params, setParams] = useState<ProcessParams>({ K: 1, M: 1, tau: 1, zeta: 1 });
    const [lockYAxis, setLockYAxis] = useState<boolean>(false);
    const [yAxisRange, setYAxisRange] = useState<[number, number]>([0, 1.1]); // Initialize with a default range
    const [simulationResult, setSimulationResult] = useState<SimulationResult | null>(null);
    const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
    const [isLoading, setIsLoading] = useState<boolean>(true); // Added loading state

    const echartsRef = useRef<ReactECharts | null>(null);

    // --- Calculation Logic ---
    const runSimulation = useCallback((): SimulationResult | null => {
        const { K, M, tau: rawTau, zeta: rawZeta } = params;
        const tau = Math.max(rawTau, 1e-9); // Avoid division by zero
        const zeta = rawZeta ?? 1; // Default zeta if undefined

        const t_final = 50;
        const num_points = 200; // Increased points for smoother curves
        const t = Array.from({ length: num_points }, (_, i) => i * t_final / (num_points - 1));

        let y: number[] = Array(t.length).fill(0);
        let y_input: number[] = Array(t.length).fill(0);
        let metrics: Metrics = { peakTime: null, overshoot: null, oscillationPeriod: null, decayRatio: null };

        try {
            if (order === 'first') {
                if (functionType === 'step') {
                    y = t.map(ti => K * M * (1 - Math.exp(-ti / tau)));
                    y_input = y_input.map(() => M);
                } else if (functionType === 'ramp') {
                    y = t.map(ti => K * M * (tau * (Math.exp(-ti / tau) - 1) + ti)); // Corrected ramp formula
                    y_input = t.map(ti => M * ti);
                }
            } else if (order === 'second') {
                if (functionType === 'step') {
                    if (zeta < 1) { // Underdamped
                        const wn = 1 / tau; // Natural frequency interpretation depends on definition, assuming tau = 1/wn here
                        const wd = wn * Math.sqrt(1 - zeta * zeta);
                        y = t.map(ti => K * M * (1 - Math.exp(-zeta * wn * ti) * (Math.cos(wd * ti) + (zeta / Math.sqrt(1 - zeta * zeta)) * Math.sin(wd * ti))));
                        metrics.peakTime = Math.PI / wd;
                        metrics.overshoot = Math.exp(-Math.PI * zeta / Math.sqrt(1 - zeta * zeta));
                        metrics.oscillationPeriod = 2 * Math.PI / wd;
                        metrics.decayRatio = metrics.overshoot * metrics.overshoot;
                    } else if (zeta === 1) { // Critically damped
                        const wn = 1 / tau;
                        y = t.map(ti => K * M * (1 - (1 + wn * ti) * Math.exp(-wn * ti)));
                    } else { // Overdamped
                        const wn = 1 / tau;
                        const r1 = wn * (-zeta + Math.sqrt(zeta * zeta - 1));
                        const r2 = wn * (-zeta - Math.sqrt(zeta * zeta - 1));
                        y = t.map(ti => K * M * (1 - (r2 * Math.exp(r1 * ti) - r1 * Math.exp(r2 * ti)) / (r2 - r1)));
                    }
                    y_input = y_input.map(() => M);
                }
                // Ramp for second order is not implemented based on the Python code
            }
            return { t, y, y_input, metrics };
        } catch (error) {
            console.error("Simulation error:", error);
            return null; // Return null on error
        }
    }, [order, functionType, params]);

    // --- ECharts Options Generation ---
    // Modified: Accepts finalYAxisRange, removed setYAxisRange call
    const generateEchartsOptions = useCallback((simResult: SimulationResult | null, finalYAxisRange: [number, number]): EChartsOption => {
        if (!simResult) {
            return {
                title: { text: 'Simulation Error or Not Run', left: 'center' },
                xAxis: { type: 'value' },
                yAxis: { type: 'value', min: finalYAxisRange[0], max: finalYAxisRange[1] }, // Use provided range
                series: [],
                backgroundColor: '#08306b'
            };
        }

        const { t, y, y_input } = simResult;

        const seriesData: SeriesOption[] = [
            {
                name: 'System Response', type: 'line', data: t.map((time, i) => [time, y[i]]),
                smooth: true, symbol: 'none', lineStyle: { color: '#ffff00', width: 2 }, // Yellow
                emphasis: { focus: 'series' }
            },
            {
                name: 'Input', type: 'line', data: t.map((time, i) => [time, y_input[i]]),
                smooth: false, symbol: 'none', lineStyle: { color: '#ff0000', type: 'dashed' }, // Red dashed
                emphasis: { focus: 'series' }
            }
        ];

        let titleText = '';
        if (order === 'first') titleText = `First Order ${functionType === 'step' ? 'Step' : 'Ramp'} Response`;
        else if (order === 'second') titleText = `Second Order ${functionType === 'step' ? 'Step' : 'Ramp'} Response`;


        return {
            backgroundColor: '#08306b',
            animation: false,
            title: {
                text: titleText, left: 'center',
                textStyle: { color: '#fff', fontSize: 18, fontFamily: 'sans-serif' }
            },
            tooltip: {
                trigger: 'axis',
                formatter: (params: any) => {
                    let tooltipText = `Time: ${params[0].axisValueLabel}<br/>`;
                    params.forEach((param: any) => {
                        tooltipText += `${param.marker}${param.seriesName}: ${param.value[1].toPrecision(4)}<br/>`;
                    });
                    return tooltipText;
                }
            },
            legend: {
                data: ['System Response', 'Input'],
                textStyle: { color: '#fff', fontSize: 12, fontFamily: 'sans-serif' },
                top: 'bottom',
                type: 'scroll'
            },
            grid: { left: '8%', right: '8%', bottom: '15%', top: '15%', containLabel: true },
            toolbox: {
                feature: { saveAsImage: { name: 'process-dynamics-plot', backgroundColor: '#08306b' } },
                iconStyle: { borderColor: '#fff' },
                orient: 'vertical', right: 10, bottom: 40
            },
            xAxis: {
                type: 'value', name: 'Time', nameLocation: 'middle', nameGap: 30,
                min: 0, max: 50, // Fixed X-axis range
                axisLabel: { color: '#fff', fontSize: 14, fontFamily: 'sans-serif', formatter: (v: number) => v.toFixed(1) },
                nameTextStyle: { color: '#fff', fontSize: 15, fontFamily: 'sans-serif' },
                axisLine: { lineStyle: { color: '#fff' } },
                axisTick: { lineStyle: { color: '#fff' } },
                splitLine: { show: false }
            },
            yAxis: {
                type: 'value', name: 'Response / Input', nameLocation: 'middle', nameGap: 50,
                min: finalYAxisRange[0], // Use the passed finalYAxisRange
                max: finalYAxisRange[1], // Use the passed finalYAxisRange
                axisLabel: { color: '#fff', fontSize: 14, fontFamily: 'sans-serif', formatter: (v: number) => v.toPrecision(3) },
                nameTextStyle: { color: '#fff', fontSize: 15, fontFamily: 'sans-serif' },
                axisLine: { lineStyle: { color: '#fff' } },
                axisTick: { lineStyle: { color: '#fff' } },
                splitLine: { show: false }
            },
            series: seriesData
        };
    }, [order, functionType]);

    // --- Effects ---
    useEffect(() => {
        setIsLoading(true);
        const result = runSimulation();
        setSimulationResult(result);

        let currentYMin = 0, currentYMax = 1.1; // Default range
        if (result) {
            const allY = [...result.y, ...result.y_input].filter(val => !isNaN(val) && isFinite(val));
            if (allY.length > 0) {
                currentYMin = Math.min(...allY);
                currentYMax = Math.max(...allY);
                const range = currentYMax - currentYMin;
                const buffer = Math.max(range * 0.1, 0.1);
                currentYMin -= buffer;
                currentYMax += buffer;
            }
        }

        // Determine the range to use for the chart
        const rangeForChart = lockYAxis ? yAxisRange : [currentYMin, currentYMax];

        // Update the stored range state *only* if not locked
        if (!lockYAxis) {
            // Check if the calculated range is different from the current state to avoid unnecessary updates
            if (yAxisRange[0] !== currentYMin || yAxisRange[1] !== currentYMax) {
                 setYAxisRange([currentYMin, currentYMax]);
            }
        }

        // Generate options using the determined rangeForChart
        // Add explicit type assertion here
        const options = generateEchartsOptions(result, rangeForChart as [number, number]); // <-- FIX: Explicit cast
        setEchartsOptions(options);
        setIsLoading(false);
    // yAxisRange is still needed here because lockYAxis depends on it.
    // generateEchartsOptions is needed as it might change based on order/functionType.
    }, [runSimulation, generateEchartsOptions, lockYAxis, yAxisRange]);

    // Reset function type if switching to second order and ramp is selected
    useEffect(() => {
        if (order === 'second' && functionType === 'ramp') {
            setFunctionType('step');
        }
    }, [order, functionType]);

    // --- Event Handlers ---
    const handleParamChange = (param: keyof ProcessParams, value: number) => {
        setParams(prev => ({ ...prev, [param]: value }));
    };

    const handleOrderChange = (newOrder: OrderType) => {
        setOrder(newOrder);
        // Reset zeta to default when switching order? Optional.
        // if (newOrder === 'first') setParams(prev => ({ ...prev, zeta: undefined }));
        // else if (newOrder === 'second' && params.zeta === undefined) setParams(prev => ({ ...prev, zeta: 1 }));
    };

    const handleFunctionChange = (newFunction: FunctionType) => {
        setFunctionType(newFunction);
    };

    // Modified: Simplified logic, relies on useEffect to update range if unlocking
    const handleLockChange = (checked: boolean | 'indeterminate') => {
        const isLocked = Boolean(checked);
        setLockYAxis(isLocked);
        // No immediate range update needed here. The main useEffect
        // will handle using/updating the range based on the new lockYAxis state.
        // If locking, it will use the existing yAxisRange state.
        // If unlocking, it will calculate and set the new yAxisRange state.
    };


    // --- Render ---
    const metrics = simulationResult?.metrics;

    return (
        <div className="container mx-auto p-4 md:p-8">
            <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                {/* Controls Column */}
                <div className="lg:col-span-1 space-y-6">
                    <Card>
                        <CardHeader><CardTitle>System Configuration</CardTitle></CardHeader>
                        <CardContent className="space-y-4">
                            <Label>Order:</Label>
                            <div className="flex gap-2">
                                <Button variant={order === 'first' ? 'default' : 'outline'} onClick={() => handleOrderChange('first')} className="flex-1">1st Order</Button>
                                <Button variant={order === 'second' ? 'default' : 'outline'} onClick={() => handleOrderChange('second')} className="flex-1">2nd Order</Button>
                            </div>
                            <Label>Input Function:</Label>
                            <div className="flex gap-2">
                                <Button variant={functionType === 'step' ? 'default' : 'outline'} onClick={() => handleFunctionChange('step')} className="flex-1">Step</Button>
                                <Button variant={functionType === 'ramp' ? 'default' : 'outline'} onClick={() => handleFunctionChange('ramp')} disabled={order === 'second'} className="flex-1">Ramp</Button>
                            </div>
                        </CardContent>
                    </Card>

                    <Card>
                        <CardHeader><CardTitle>Parameters</CardTitle></CardHeader>
                        <CardContent className="space-y-6">
                            {/* K Slider */}
                            <div className="space-y-2">
                                <Label htmlFor="k-slider" className="flex justify-between"><span>Gain (K):</span><span>{formatNumber(params.K)}</span></Label>
                                <Slider id="k-slider" min={0} max={10} step={0.1} value={[params.K]} onValueChange={([val]) => handleParamChange('K', val)} />
                            </div>
                            {/* M Slider */}
                            <div className="space-y-2">
                                <Label htmlFor="m-slider" className="flex justify-between"><span>Magnitude (M):</span><span>{formatNumber(params.M)}</span></Label>
                                <Slider id="m-slider" min={0} max={10} step={0.1} value={[params.M]} onValueChange={([val]) => handleParamChange('M', val)} />
                            </div>
                            {/* Tau Slider */}
                            <div className="space-y-2">
                                <Label htmlFor="tau-slider" className="flex justify-between"><span>Time Constant (&tau;):</span><span>{formatNumber(params.tau)}</span></Label>
                                <Slider id="tau-slider" min={0.1} max={10} step={0.1} value={[params.tau]} onValueChange={([val]) => handleParamChange('tau', val)} />
                            </div>
                            {/* Zeta Slider (Conditional) */}
                            {order === 'second' && (
                                <div className="space-y-2">
                                    <Label htmlFor="zeta-slider" className="flex justify-between"><span>Damping Ratio (&zeta;):</span><span>{formatNumber(params.zeta)}</span></Label>
                                    <Slider id="zeta-slider" min={0} max={2} step={0.01} value={[params.zeta ?? 1]} onValueChange={([val]) => handleParamChange('zeta', val)} />
                                </div>
                            )}
                        </CardContent>
                    </Card>

                    <Card>
                         <CardHeader><CardTitle>Display Options</CardTitle></CardHeader>
                         <CardContent>
                            <div className="flex items-center space-x-2">
                                <Checkbox id="lock-y-axis" checked={lockYAxis} onCheckedChange={handleLockChange} />
                                <Label htmlFor="lock-y-axis">Lock Y-axis</Label>
                            </div>
                         </CardContent>
                    </Card>

                    {/* Metrics Card (Conditional) */}
                    {order === 'second' && metrics && (
                        <Card>
                            <CardHeader><CardTitle>Performance Metrics</CardTitle></CardHeader>
                            <CardContent className="space-y-1 text-sm">
                                {isLoading ? (
                                    <> <Skeleton className="h-4 w-3/4 mb-1"/> <Skeleton className="h-4 w-2/3 mb-1"/> <Skeleton className="h-4 w-3/4 mb-1"/> <Skeleton className="h-4 w-2/3 mb-1"/> </>
                                ) : (
                                    <>
                                        <p>Peak Time: {formatNumber(metrics.peakTime)}</p>
                                        <p>Overshoot Ratio: {formatNumber(metrics.overshoot)}</p>
                                        <p>Oscillation Period: {formatNumber(metrics.oscillationPeriod)}</p>
                                        <p>Decay Ratio: {formatNumber(metrics.decayRatio)}</p>
                                    </>
                                )}
                            </CardContent>
                        </Card>
                    )}
                </div>

                {/* Plot Column */}
                <div className="lg:col-span-2">
                    <Card>
                        <CardContent className="pt-6">
                            <div className="relative h-[500px] md:h-[600px] rounded-md overflow-hidden border bg-card">
                                {isLoading && (
                                    <div className="absolute inset-0 flex items-center justify-center bg-background/60 dark:bg-background/70 z-10 backdrop-blur-sm">
                                        <Skeleton className="h-3/4 w-3/4 rounded-md"/>
                                    </div>
                                )}
                                {Object.keys(echartsOptions).length > 0 ? (
                                    <ReactECharts
                                        ref={echartsRef}
                                        echarts={echarts}
                                        option={echartsOptions}
                                        style={{ height: '100%', width: '100%' }}
                                        notMerge={true}
                                        lazyUpdate={false}
                                    />
                                ) : (
                                     <div className="absolute inset-0 flex items-center justify-center text-muted-foreground">
                                        {isLoading ? 'Loading...' : 'Configure parameters to view simulation.'}
                                    </div>
                                )}
                            </div>
                        </CardContent>
                    </Card>
                </div>
            </div>
        </div>
    );
}
