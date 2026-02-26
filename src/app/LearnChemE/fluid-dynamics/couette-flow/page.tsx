'use client';

import { useState, useEffect, useMemo } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
// Removed import of CouetteFlowSimulation as its content will be inlined

// Imports previously in CouetteFlowSimulation.tsx
import { Slider } from "@/components/ui/slider";
import { Label } from "@/components/ui/label";
// Input might not be used if SliderWithValue handles display, but keep if needed elsewhere
// import { Input } from "@/components/ui/input"; 
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import { LineChart } from 'echarts/charts';
import {
    TitleComponent,
    TooltipComponent,
    GridComponent,
    LegendComponent,
    MarkPointComponent,
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';
import { useTheme } from "next-themes";

echarts.use([
    TitleComponent, TooltipComponent, GridComponent, LegendComponent, LineChart, CanvasRenderer, MarkPointComponent
]);

// --- SliderWithValue component (previously in CouetteFlowSimulation.tsx) ---
const SliderWithValue = ({ id, label, min, max, step, value, unit, onValueChange, precision = 2, displayMultiplier = 1 }: {
    id: string, label: string, min: number, max: number, step: number, value: number, unit: string, onValueChange: (val: number) => void, precision?: number, displayMultiplier?: number
}) => (
    <div className="space-y-2">
        <div className="flex justify-between items-center">
            <Label htmlFor={id} className="text-sm font-medium">{label}</Label>
            <span className="text-sm text-muted-foreground w-28 text-right">{(value * displayMultiplier).toFixed(precision)} {unit}</span>
        </div>
        <Slider
            id={id}
            min={min}
            max={max}
            step={step}
            value={[value]}
            onValueChange={(vals) => {
                // console.log(`Slider ${id} event. New raw value:`, vals); // Keep if debugging
                if (vals && vals.length > 0 && typeof vals[0] === 'number') {
                    onValueChange(vals[0]);
                } else {
                    // console.warn(`Slider ${id} event. Invalid value received:`, vals); // Keep if debugging
                }
            }}
            className="w-full"
        />
    </div>
);

export default function CouetteFlowPage() {
    // --- State and Logic from CouetteFlowSimulation ---
    const { resolvedTheme } = useTheme();
    const [plateVelocityU, setPlateVelocityU] = useState<number>(0); // m/s
    const [plateDistanceH, setPlateDistanceH] = useState<number>(0.02); // m (2 cm)
    const [fluidViscosityMu, setFluidViscosityMu] = useState<number>(0.0001); // Pa.s
    const [pressureGradientDpDx, setPressureGradientDpDx] = useState<number>(-0.15); // Pa/m

    // DEBUG LOG: Log state values on re-render (optional, can be removed)
    // console.log(
    //     `CouetteFlowPage RENDER. U: ${plateVelocityU}, h: ${plateDistanceH}, mu: ${fluidViscosityMu}, dpdx: ${pressureGradientDpDx}`
    // );

    const [chartOptions, setChartOptions] = useState<echarts.EChartsCoreOption>({});
    const numPoints = 51; // Number of points for plotting the profile

    const velocityProfileData = useMemo(() => {
        const data: [number, number][] = [];
        for (let i = 0; i < numPoints; i++) {
            const y_m = (i / (numPoints - 1)) * plateDistanceH; // y is still in meters for calculation
            const term1 = (plateVelocityU * y_m) / plateDistanceH;
            const term2 = pressureGradientDpDx !== 0 && fluidViscosityMu !== 0 ?
                (1 / (2 * fluidViscosityMu)) * pressureGradientDpDx * (y_m * y_m - plateDistanceH * y_m) : 0;
            const u_y = term1 + term2;
            data.push([u_y, y_m * 100]); // Convert y to cm for plotting [velocity, y_position_cm]
        }
        return data;
    }, [plateVelocityU, plateDistanceH, fluidViscosityMu, pressureGradientDpDx, numPoints]);

    const shearStressBottom = useMemo(() => {
        const du_dy_0 = (plateVelocityU / plateDistanceH) -
                       (pressureGradientDpDx !== 0 && fluidViscosityMu !== 0 ? (plateDistanceH / (2 * fluidViscosityMu)) * pressureGradientDpDx : 0);
        return fluidViscosityMu * du_dy_0;
    }, [plateVelocityU, plateDistanceH, fluidViscosityMu, pressureGradientDpDx]);

    const markPointData = useMemo(() => {
        const points = [];
        // Only consider marking if there's a pressure gradient that could cause a parabolic peak
        if (pressureGradientDpDx < 0) {
            // y-coordinate (in meters) where shear stress is zero (du/dy = 0)
            // This formula is derived from setting du/dy = 0 for the combined Couette-Poiseuille flow equation
            const y_peak_m = (plateDistanceH / 2) - (fluidViscosityMu * plateVelocityU) / (plateDistanceH * pressureGradientDpDx);

            // Check if this peak is strictly within the plates
            if (y_peak_m > 0 && y_peak_m < plateDistanceH) {
                // Calculate velocity at this y_peak_m
                const term1 = (plateVelocityU * y_peak_m) / plateDistanceH;
                const term2 = (1 / (2 * fluidViscosityMu)) * pressureGradientDpDx * (y_peak_m * y_peak_m - plateDistanceH * y_peak_m);
                const u_at_peak = term1 + term2;

                points.push({
                    name: 'Parabolic Max', // Optional: for tooltip on the mark point
                    coord: [u_at_peak, y_peak_m * 100], // [velocity (m/s), y_position (cm)]
                    itemStyle: { color: 'red' },
                    symbol: 'circle',
                    symbolSize: 8
                });
            }
        }
        return points;
    }, [plateVelocityU, plateDistanceH, fluidViscosityMu, pressureGradientDpDx]);

    useEffect(() => {
        const isDark = resolvedTheme === 'dark';
        const textColor = isDark ? '#E5E7EB' : '#1F2937';
        const profileColor = isDark ? '#A7F3D0' : '#10B981';
        const chartBackgroundColor = isDark ? 'hsl(215, 40%, 30%)' : 'hsl(210, 40%, 96.1%)';
        const chartFontFamily = '"Merriweather Sans", sans-serif';

        setChartOptions({
            backgroundColor: chartBackgroundColor,
            title: {
                text: 'Velocity Profile u(y)',
                left: 'center',
                textStyle: { color: textColor, fontSize: 16, fontFamily: chartFontFamily }
            },
            tooltip: {
                trigger: 'axis',
                formatter: (params: any) => {
                    const point = params[0];
                    return `Velocity: ${parseFloat(point.value[0]).toFixed(4)} m/s<br />Height (y): ${parseFloat(point.value[1]).toFixed(4)} cm`; // y-axis is in cm
                },
                textStyle: { // Added font family for tooltip
                    fontFamily: chartFontFamily
                }
            },
            grid: { left: '15%', right: '10%', bottom: '15%', top: '15%' },
            xAxis: {
                type: 'value',
                name: 'Velocity u (m/s)',
                nameLocation: 'middle',
                nameGap: 30,
                axisLabel: { color: textColor, formatter: (val: number) => val.toFixed(4), fontFamily: chartFontFamily },
                nameTextStyle: { color: textColor, fontSize: 14, fontWeight: 'bold', fontFamily: chartFontFamily },
                axisLine: { lineStyle: { color: '#FFFFFF' } },
                splitLine: { show: false },
            },
            yAxis: {
                type: 'value',
                name: 'Distance y (cm)',
                nameLocation: 'middle',
                nameGap: 50,
                min: 0,
                max: plateDistanceH * 100,
                axisLabel: { color: textColor, formatter: (val: number) => val.toFixed(Math.max(0, Math.ceil(-Math.log10(plateDistanceH * 100))+1)), fontFamily: chartFontFamily },
                nameTextStyle: { color: textColor, fontSize: 14, fontWeight: 'bold', fontFamily: chartFontFamily },
                axisLine: { lineStyle: { color: '#FFFFFF' } },
                splitLine: { show: false },
            },
            series: [{
                name: 'Velocity Profile',
                type: 'line',
                smooth: true,
                symbol: 'none',
                lineStyle: { color: profileColor, width: 3 },
                data: velocityProfileData,
                animation: false,
                markPoint: { // Added markPoint configuration
                    data: markPointData
                }
            }]
        });
    }, [velocityProfileData, plateDistanceH, plateVelocityU, resolvedTheme, markPointData, fluidViscosityMu, pressureGradientDpDx]); // Added markPointData and its dependencies to ensure updates
    // --- End of State and Logic from CouetteFlowSimulation ---

    return (
        <div className="container mx-auto px-4 py-8 min-h-screen">
            <Card className="mt-6">
                <CardHeader>
                    {/* Title was removed as per previous request */}
                </CardHeader>
                <CardContent className="mt-4">
                    {/* --- JSX from CouetteFlowSimulation --- */}
                    <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
                        <div className="md:col-span-1 space-y-6">
                            <Card>
                                <CardHeader>
                                    <CardTitle>Parameters</CardTitle>
                                    <CardDescription>Adjust the flow conditions.</CardDescription>
                                </CardHeader>
                                <CardContent className="space-y-4">
                                    <SliderWithValue id="plateVelocityU" label="Top Plate Velocity (U)" min={0} max={5} step={0.1} value={plateVelocityU} unit="m/s" onValueChange={setPlateVelocityU} precision={1} />
                                    <SliderWithValue id="plateDistanceH" label="Plate Spacing (h)" min={0.001} max={0.1} step={0.001} value={plateDistanceH} unit="cm" onValueChange={setPlateDistanceH} precision={1} displayMultiplier={100} />
                                    <SliderWithValue id="fluidViscosityMu" label="Viscosity (µ)" min={0.0001} max={0.1} step={0.0001} value={fluidViscosityMu} unit="Pa·s" onValueChange={setFluidViscosityMu} precision={4} />
                                    <SliderWithValue id="pressureGradientDpDx" label="Pressure Gradient (dp/dx)" min={-0.3} max={0.1} step={0.01} value={pressureGradientDpDx} unit="Pa/m" onValueChange={setPressureGradientDpDx} precision={2} />
                                </CardContent>
                            </Card>
                            <Card>
                                <CardHeader>
                                    <CardTitle>Results</CardTitle>
                                </CardHeader>
                                <CardContent>
                                    <div className="flex justify-between text-sm">
                                        <span>Shear Stress (bottom wall, τ₀):</span>
                                        <span className="font-medium">{shearStressBottom.toPrecision(3)} Pa</span>
                                    </div>
                                </CardContent>
                            </Card>
                        </div>

                        <div className="md:col-span-2">
                            <Card className="h-[500px] md:h-[600px]">
                                <CardContent className="h-full p-2 md:p-4">
                                    {Object.keys(chartOptions).length > 0 ? (
                                        <ReactECharts
                                            echarts={echarts}
                                            option={chartOptions}
                                            style={{ height: '100%', width: '100%' }}
                                            notMerge={true}
                                            lazyUpdate={true}
                                            theme={resolvedTheme === 'dark' ? 'dark' : undefined}
                                        />
                                    ) : (
                                        <div className="flex items-center justify-center h-full text-muted-foreground">Loading chart...</div>
                                    )}
                                </CardContent>
                            </Card>
                        </div>
                    </div>
                    {/* --- End of JSX from CouetteFlowSimulation --- */}
                </CardContent>
            </Card>

            <footer className="py-6 mt-8 text-center text-sm text-muted-foreground">
                © {new Date().getFullYear()} Interactive Simulation.
            </footer>
        </div>
    );
}
