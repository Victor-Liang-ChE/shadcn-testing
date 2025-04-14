"use client";
import React, { useState, useEffect, useCallback, useMemo } from 'react';
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsOption, SeriesOption } from 'echarts';
import { LineChart } from 'echarts/charts';
import { TitleComponent, TooltipComponent, GridComponent, LegendComponent, ToolboxComponent } from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Label } from "@/components/ui/label";
import { Slider } from "@/components/ui/slider";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Input } from "@/components/ui/input";
import { Checkbox } from "@/components/ui/checkbox"; // Keep Checkbox for other sections
import { RadioGroup, RadioGroupItem } from "@/components/ui/radio-group";

// ECharts imports
echarts.use([TitleComponent, TooltipComponent, GridComponent, LegendComponent, ToolboxComponent, LineChart, CanvasRenderer]);

// --- Type Definitions ---
type CalculationType =
    | 'compounding_interest'
    | 'ear'
    | 'npv'
    | 'inflation'
    | 'annuity'
    | 'perpetuity'
    | 'depreciation'
    | 'cashflow'; // Added cashflow for completeness, though not implemented yet

// --- Helper Functions ---

// Abbreviate numbers (similar to Python version)
function abbreviateNumber(num: number): string {
    if (num >= 1e93) return '∞';
    const suffixes = [
        { value: 1e93, symbol: 'Tg' }, { value: 1e90, symbol: 'Nvg' }, { value: 1e87, symbol: 'Ovg' },
        { value: 1e84, symbol: 'Spvg' }, { value: 1e81, symbol: 'Sxvg' }, { value: 1e78, symbol: 'Qvg' },
        { value: 1e75, symbol: 'Qavg' }, { value: 1e72, symbol: 'Qav' }, { value: 1e69, symbol: 'Dvg' },
        { value: 1e66, symbol: 'Uvg' }, { value: 1e63, symbol: 'Vg' }, { value: 1e60, symbol: 'Nvg' }, // Corrected duplicate Vg
        { value: 1e57, symbol: 'Ocdc' }, { value: 1e54, symbol: 'Spdc' }, // Corrected Spvg
        { value: 1e51, symbol: 'Sxdc' }, { value: 1e48, symbol: 'Qidc' }, { value: 1e45, symbol: 'Qadc' },
        { value: 1e42, symbol: 'Tdc' }, { value: 1e39, symbol: 'Ddc' }, { value: 1e36, symbol: 'Udc' },
        { value: 1e33, symbol: 'Dc' }, { value: 1e30, symbol: 'Nm' }, { value: 1e27, symbol: 'Oc' },
        { value: 1e24, symbol: 'Sp' }, { value: 1e21, symbol: 'Sx' }, { value: 1e18, symbol: 'Qi' },
        { value: 1e15, symbol: 'Qa' }, { value: 1e12, symbol: 'T' }, { value: 1e9, symbol: 'B' },
        { value: 1e6, symbol: 'M' }, { value: 1e3, symbol: 'K' }
    ];

    for (let i = 0; i < suffixes.length; i++) {
        if (num >= suffixes[i].value && num < 1e93) {
            return (num / suffixes[i].value).toFixed(2) + suffixes[i].symbol;
        }
    }
    return num.toLocaleString('en-US', { minimumFractionDigits: 2, maximumFractionDigits: 2 });
}

// Format number to string with fixed precision
const formatNumber = (num: number | null | undefined, precision = 2): string => {
    if (num === null || num === undefined || isNaN(num)) return "---";
    return num.toFixed(precision);
};

// --- Main Component ---
export default function ChemEEconPage() {
    const [selectedCalculation, setSelectedCalculation] = useState<CalculationType>('compounding_interest');
    const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
    const [results, setResults] = useState<Record<string, string | number>>({});

    // --- State for Inputs (Grouped by Calculation Type) ---
    // Compounding Interest
    const [ciPV, setCiPV] = useState<number>(1000);
    const [ciYears, setCiYears] = useState<number>(1);
    const [ciCompoundsPerYear, setCiCompoundsPerYear] = useState<number>(1);
    const [ciInterestRate, setCiInterestRate] = useState<number>(5);
    const [ciCompareSimple, setCiCompareSimple] = useState<boolean>(false);

    // EAR - Removed earCompareMax state
    const [earCompounds, setEarCompounds] = useState<number>(1);
    const [earInterestRate, setEarInterestRate] = useState<number>(5);

    // NPV (Future Value Calculation)
    const [npvPV, setNpvPV] = useState<number>(1000);
    const [npvYears, setNpvYears] = useState<number>(1);
    const [npvDiscountRate, setNpvDiscountRate] = useState<number>(5);

    // Annuity
    const [annuityPV, setAnnuityPV] = useState<number>(1000);
    const [annuityYears, setAnnuityYears] = useState<number>(10);
    const [annuityInterestRate, setAnnuityInterestRate] = useState<number>(5);

    // Perpetuity
    const [perpCashFlow, setPerpCashFlow] = useState<number>(1000);
    const [perpInterestRate, setPerpInterestRate] = useState<number>(5);

    // Inflation
    const [infPV, setInfPV] = useState<number>(1000);
    const [infYears, setInfYears] = useState<number>(10);
    const [infInterestRate, setInfInterestRate] = useState<number>(5);
    const [infInflationRate, setInfInflationRate] = useState<number>(2);

    // Depreciation
    const [depFCI, setDepFCI] = useState<number>(10000);
    const [depLife, setDepLife] = useState<number>(10);
    const [depSalvage, setDepSalvage] = useState<number>(1000);
    const [depMethod, setDepMethod] = useState<'straight' | 'macrs'>('straight');
    const [depCompare, setDepCompare] = useState<boolean>(false);

    // --- Calculation and Plotting Logic ---
    const updateChartAndResults = useCallback(() => {
        let newOptions: EChartsOption = {};
        let newResults: Record<string, string | number> = {};
        const plotBgColor = '#08306b'; // Consistent dark blue background
        const textColor = '#fff'; // White text for dark background
        const axisColor = '#fff';
        const lineColors = ['#ffff00', '#ff0000', '#00ff00', '#00ffff']; // Yellow, Red, Green, Cyan

        try {
            switch (selectedCalculation) {
                case 'compounding_interest': {
                    const pv = ciPV;
                    const n = ciYears;
                    const m = ciCompoundsPerYear; // Compounds per year
                    const r = ciInterestRate / 100;
                    if (n <= 0 || m <= 0) throw new Error("Years and compounds must be positive.");

                    const totalPeriods = n * m;
                    const periodRate = r / m;
                    const fvCompound = pv * Math.pow(1 + periodRate, totalPeriods);
                    const fvSimple = pv * (1 + n * r);

                    newResults.fvCompoundText = `Future Value: $${abbreviateNumber(fvCompound)}`;
                    if (ciCompareSimple) {
                        newResults.fvSimpleText = `Future Value (Simple): $${abbreviateNumber(fvSimple)}`;
                        newResults.differenceText = `Difference: $${abbreviateNumber(fvCompound - fvSimple)}`;
                    }

                    const x_values = Array.from({ length: totalPeriods + 1 }, (_, x) => x / m); // Time in years
                    const y_compound = Array.from({ length: totalPeriods + 1 }, (_, x) => pv * Math.pow(1 + periodRate, x));
                    const y_simple = x_values.map(t => pv * (1 + r * t));

                    const series: SeriesOption[] = [{
                        name: 'Future Value (Compound)', type: 'line', data: x_values.map((x, i) => [x, y_compound[i]]),
                        lineStyle: { color: lineColors[0] }, symbol: 'none', smooth: true
                    }];
                    if (ciCompareSimple) {
                        series.push({
                            name: 'Future Value (Simple)', type: 'line', data: x_values.map((x, i) => [x, y_simple[i]]),
                            lineStyle: { color: lineColors[1] }, symbol: 'none', smooth: true
                        });
                    }

                    // Filter out undefined names before assigning to legend data
                    const legendData = series.map(s => s.name).filter((name): name is string => !!name);

                    newOptions = {
                        animation: false, // Disable global animations
                        backgroundColor: plotBgColor, title: { text: 'Future Value vs. Years', left: 'center', textStyle: { color: textColor } },
                        tooltip: { trigger: 'axis', valueFormatter: (value) => `$${formatNumber(value as number)}` },
                        legend: { data: legendData, top: 'bottom', textStyle: { color: textColor } },
                        grid: { containLabel: true, left: '10%', right: '10%' },
                        xAxis: { type: 'value', name: 'Years', min: 0, max: n, axisLabel: { color: textColor }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } }, // No grid lines
                        yAxis: { type: 'value', name: 'Future Value ($)', axisLabel: { color: textColor, formatter: (val: number) => abbreviateNumber(val) }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } }, // No grid lines
                        series: series
                    };
                    break;
                }
                case 'ear': {
                    const m = earCompounds;
                    const r = earInterestRate / 100;
                    if (m <= 0) throw new Error("Number of compounds must be positive.");

                    const ear = Math.pow(1 + r / m, m) - 1;
                    const maxEar = Math.exp(r) - 1; // Always calculate max EAR

                    // Always set results for max EAR
                    newResults.earText = `EAR: ${formatNumber(ear * 100)}%`;
                    newResults.maxEarText = `Max EAR: ${formatNumber(maxEar * 100)}%`;
                    newResults.differenceText = `Difference: ${formatNumber((maxEar - ear) * 100)}%`;

                    const compoundsRange = Array.from({ length: Math.max(m, 50) + 1 }, (_, i) => i).filter(i => i > 0);
                    const ears = compoundsRange.map(comp => Math.pow(1 + r / comp, comp) - 1);

                    const series: SeriesOption[] = [{
                        name: 'EAR', type: 'line', data: compoundsRange.map((comp, i) => [comp, ears[i] * 100]),
                        lineStyle: { color: lineColors[0] }, symbol: 'none', smooth: true
                    }];

                    // Always define shapes for the max EAR line
                    const shapes = [{
                        type: 'line' as const, x0: 1, y0: maxEar * 100, x1: Math.max(m, 50), y1: maxEar * 100,
                        lineStyle: { color: lineColors[1], type: 'dashed' as const }
                    }];

                    // Always define legend data including Max EAR
                    const legendData: string[] = ['EAR', 'Max EAR'];

                    newOptions = {
                        animation: false,
                        backgroundColor: plotBgColor, title: { text: 'EAR vs. Number of Compounds', left: 'center', textStyle: { color: textColor } },
                        tooltip: { trigger: 'axis', valueFormatter: (value) => `${formatNumber(value as number)}%` },
                        legend: { data: legendData, top: 'bottom', textStyle: { color: textColor } },
                        grid: { containLabel: true, left: '10%', right: '10%' },
                        xAxis: { type: 'value', name: 'Compounds per Year', min: 1, axisLabel: { color: textColor }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } },
                        yAxis: { type: 'value', name: 'EAR (%)', axisLabel: { color: textColor, formatter: '{value}%' }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } },
                        series: series,
                        graphic: { elements: shapes } // Always include graphic element
                    };
                    break;
                }
                case 'npv': { // Replicating Python's "Future Value" calculation
                    const pv = npvPV;
                    const n = npvYears;
                    const r = npvDiscountRate / 100;
                    if (n <= 0) throw new Error("Years must be positive.");

                    const fv = pv / Math.pow(1 + r, n); // This is actually PV of a future value, let's call it Discounted Value
                    newResults.fvText = `Discounted Value: $${formatNumber(fv)}`;

                    const yearsList = Array.from({ length: n + 1 }, (_, i) => i);
                    const discountedValues = yearsList.map(year => pv / Math.pow(1 + r, year));

                    const seriesName = 'Discounted Value';
                    newOptions = {
                        animation: false, // Disable global animations
                        backgroundColor: plotBgColor, title: { text: 'Discounted Value vs. Years', left: 'center', textStyle: { color: textColor } },
                        tooltip: { trigger: 'axis', valueFormatter: (value) => `$${formatNumber(value as number)}` },
                        legend: { data: [seriesName], top: 'bottom', textStyle: { color: textColor } },
                        grid: { containLabel: true, left: '10%', right: '10%' },
                        xAxis: { type: 'value', name: 'Years', min: 0, max: n, axisLabel: { color: textColor }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } }, // No grid lines
                        yAxis: { type: 'value', name: 'Discounted Value ($)', axisLabel: { color: textColor, formatter: (val: number) => abbreviateNumber(val) }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } }, // No grid lines
                        series: [{
                            name: seriesName, type: 'line', data: yearsList.map((yr, i) => [yr, discountedValues[i]]),
                            lineStyle: { color: lineColors[0] }, symbol: 'none', smooth: true
                        }]
                    };
                    break;
                }
                case 'annuity': {
                    const pv = annuityPV;
                    const n = annuityYears;
                    const r = annuityInterestRate / 100;
                    if (n <= 0) throw new Error("Years must be positive.");

                    const annuityPayment = (r === 0) ? pv / n : pv * (r * Math.pow(1 + r, n)) / (Math.pow(1 + r, n) - 1);
                    newResults.annuityPaymentText = `Required Annuity Payment: $${formatNumber(annuityPayment)}`;

                    const yearsRange = Array.from({ length: n }, (_, i) => i + 1);
                    const payments = yearsRange.map(yr => (r === 0) ? pv / yr : pv * (r * Math.pow(1 + r, yr)) / (Math.pow(1 + r, yr) - 1));

                    const seriesName = 'Annuity Payment';
                    newOptions = {
                        animation: false, // Disable global animations
                        backgroundColor: plotBgColor, title: { text: 'Required Annuity Payment vs. Years', left: 'center', textStyle: { color: textColor } },
                        tooltip: { trigger: 'axis', valueFormatter: (value) => `$${formatNumber(value as number)}` },
                        legend: { data: [seriesName], top: 'bottom', textStyle: { color: textColor } },
                        grid: { containLabel: true, left: '10%', right: '10%' },
                        xAxis: { type: 'value', name: 'Years', min: 1, max: n, axisLabel: { color: textColor }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } }, // No grid lines
                        yAxis: { type: 'value', name: 'Annuity Payment ($)', axisLabel: { color: textColor, formatter: (val: number) => abbreviateNumber(val) }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } }, // No grid lines
                        series: [{
                            name: seriesName, type: 'line', data: yearsRange.map((yr, i) => [yr, payments[i]]),
                            lineStyle: { color: lineColors[0] }, symbol: 'none', smooth: true
                        }]
                    };
                    break;
                }
                 case 'perpetuity': {
                    const a = perpCashFlow; // Uniform cash flow
                    const r = perpInterestRate / 100;
                    if (r <= 0) throw new Error("Interest rate must be positive.");

                    const pv = a / r;
                    newResults.pvText = `Present Value: $${formatNumber(pv)}`;

                    const cashFlowRange = Array.from({ length: 100 }, (_, i) => (i + 1) * (a / 10)); // Vary cash flow
                    const pvValues = cashFlowRange.map(cf => cf / r);

                    const seriesName = 'Present Value';
                    newOptions = {
                        animation: false, // Disable global animations
                        backgroundColor: plotBgColor, title: { text: 'Present Value vs. Uniform Cash Flow', left: 'center', textStyle: { color: textColor } },
                        tooltip: { trigger: 'axis', valueFormatter: (value) => `$${formatNumber(value as number)}` },
                        legend: { data: [seriesName], top: 'bottom', textStyle: { color: textColor } },
                        grid: { containLabel: true, left: '10%', right: '10%' },
                        xAxis: { type: 'value', name: 'Uniform Cash Flow ($)', axisLabel: { color: textColor, formatter: (val: number) => abbreviateNumber(val) }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } }, // No grid lines
                        yAxis: { type: 'value', name: 'Present Value ($)', axisLabel: { color: textColor, formatter: (val: number) => abbreviateNumber(val) }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } }, // No grid lines
                        series: [{
                            name: seriesName, type: 'line', data: cashFlowRange.map((cf, i) => [cf, pvValues[i]]),
                            lineStyle: { color: lineColors[0] }, symbol: 'none', smooth: true
                        }]
                    };
                    break;
                }
                case 'inflation': {
                    const pv = infPV;
                    const n = infYears;
                    const r = infInterestRate / 100;
                    const i = infInflationRate / 100; // Inflation rate
                    if (n <= 0) throw new Error("Years must be positive.");
                    if (1 + i <= 0) throw new Error("Inflation rate invalid.");

                    const realRateFactor = (1 + r) / (1 + i) - 1; // More stable calculation: (r-i)/(1+i)
                    const fpp = pv * Math.pow(1 + realRateFactor, n);
                    newResults.fppText = `Future Purchasing Power: $${formatNumber(fpp)}`;

                    const yearsList = Array.from({ length: n + 1 }, (_, y) => y);
                    const fppList = yearsList.map(yr => pv * Math.pow(1 + realRateFactor, yr));

                    const seriesName = 'Future Purchasing Power';
                    newOptions = {
                        animation: false, // Disable global animations
                        backgroundColor: plotBgColor, title: { text: 'Future Purchasing Power vs. Years', left: 'center', textStyle: { color: textColor } },
                        tooltip: { trigger: 'axis', valueFormatter: (value) => `$${formatNumber(value as number)}` },
                        legend: { data: [seriesName], top: 'bottom', textStyle: { color: textColor } },
                        grid: { containLabel: true, left: '10%', right: '10%' },
                        xAxis: { type: 'value', name: 'Years', min: 0, max: n, axisLabel: { color: textColor }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } }, // No grid lines
                        yAxis: { type: 'value', name: 'Future Purchasing Power ($)', axisLabel: { color: textColor, formatter: (val: number) => abbreviateNumber(val) }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } }, // No grid lines
                        series: [{
                            name: seriesName, type: 'line', data: yearsList.map((yr, i) => [yr, fppList[i]]),
                            lineStyle: { color: lineColors[0] }, symbol: 'none', smooth: true
                        }]
                    };
                    break;
                }
                case 'depreciation': {
                    const fci = depFCI;
                    const life = depLife;
                    const salvage = depSalvage;
                    if (life <= 0) throw new Error("Equipment life must be positive.");

                    const years = Array.from({ length: life + 1 }, (_, i) => i); // 0 to life

                    // Straight Line
                    const slDepPerYear = (fci - salvage) / life;
                    const slBookValue = years.map(yr => Math.max(salvage, fci - yr * slDepPerYear));
                    const slAnnualDep = years.map(yr => yr > 0 && yr <= life ? slDepPerYear : 0); // Annual dep amount

                    // MACRS (Simplified - using Python logic's structure)
                    // Note: Actual MACRS uses predefined tables based on asset class and recovery period.
                    // This JS implementation mimics the Python logic provided, which seems to be a custom DDB with switch.
                    const macrsAnnualDep: number[] = [0]; // Depreciation in year 0 is 0
                    const macrsBookValue: number[] = [fci];
                    let currentBookValue = fci;
                    const ddbRate = 2 / life;
                    let switchedToSL = false;

                    for (let yr = 1; yr <= life; yr++) {
                        let depreciationAmount = 0;
                        const remainingLifeForSL = life - yr + 1; // For SL switch check

                        // SL depreciation from current book value
                        const slFromCurrent = Math.max(0, (currentBookValue - salvage) / remainingLifeForSL);

                        // DDB depreciation
                        const ddbDep = currentBookValue * ddbRate;

                        // Apply half-year convention for first and potentially last year (simplified)
                        const factor = (yr === 1 || (yr === life && !switchedToSL)) ? 0.5 : 1.0; // Simplified half-year

                        if (!switchedToSL && ddbDep * factor > slFromCurrent * factor) {
                            depreciationAmount = ddbDep * factor;
                        } else {
                            switchedToSL = true;
                            depreciationAmount = slFromCurrent * factor; // Use SL amount
                        }

                        // Ensure book value doesn't go below salvage
                        depreciationAmount = Math.min(depreciationAmount, Math.max(0, currentBookValue - salvage));

                        macrsAnnualDep.push(depreciationAmount);
                        currentBookValue -= depreciationAmount;
                        macrsBookValue.push(currentBookValue);
                    }
                    // Adjust MACRS arrays to match length of years (0 to life) if needed
                    while (macrsAnnualDep.length < years.length) macrsAnnualDep.push(0);
                    while (macrsBookValue.length < years.length) macrsBookValue.push(salvage);


                    newResults.totalSLDep = `Total SL Dep: $${formatNumber(slAnnualDep.reduce((a, b) => a + b, 0))}`;
                    newResults.totalMACRSDep = `Total MACRS Dep: $${formatNumber(macrsAnnualDep.reduce((a, b) => a + b, 0))}`;

                    const series: SeriesOption[] = [];
                    const legendData: string[] = []; // Use string array for legend data

                    if (depMethod === 'straight' || depCompare) {
                        series.push({
                            name: 'Straight Line (Annual)', type: 'bar', data: years.map((yr, i) => [yr, slAnnualDep[i]]),
                            itemStyle: { color: lineColors[1] }, barGap: '-100%' // Overlay bars
                        });
                        legendData.push('Straight Line (Annual)');
                    }
                    if (depMethod === 'macrs' || depCompare) {
                        series.push({
                            name: 'MACRS (Annual)', type: 'bar', data: years.map((yr, i) => [yr, macrsAnnualDep[i]]),
                            itemStyle: { color: lineColors[0] }
                        });
                         legendData.push('MACRS (Annual)');
                    }
                     if (depCompare) {
                         series.push({
                             name: 'SL Book Value', type: 'line', data: years.map((yr, i) => [yr, slBookValue[i]]),
                             lineStyle: { color: lineColors[1], type: 'dashed' }, symbol: 'none', yAxisIndex: 1 // Use secondary axis
                         });
                         series.push({
                             name: 'MACRS Book Value', type: 'line', data: years.map((yr, i) => [yr, macrsBookValue[i]]),
                             lineStyle: { color: lineColors[0], type: 'dashed' }, symbol: 'none', yAxisIndex: 1 // Use secondary axis
                         });
                         legendData.push('SL Book Value', 'MACRS Book Value');
                     }


                    newOptions = {
                        animation: false, // Disable global animations
                        backgroundColor: plotBgColor, title: { text: 'Annual Depreciation vs. Year', left: 'center', textStyle: { color: textColor } },
                        tooltip: { trigger: 'axis', valueFormatter: (value) => `$${formatNumber(value as number)}` },
                        legend: { data: legendData, top: 'bottom', textStyle: { color: textColor } },
                        grid: { containLabel: true, left: '10%', right: '10%' },
                        xAxis: { type: 'value', name: 'Year', min: 0, max: life, axisLabel: { color: textColor }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } }, // No grid lines
                        yAxis: [
                            {
                                type: 'value', name: 'Annual Depreciation ($)', axisLabel: { color: textColor, formatter: (val: number) => abbreviateNumber(val) }, nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } // No grid lines
                            },
                             // Secondary axis for book value (optional)
                             ...(depCompare ? [{
                                 type: 'value' as const, name: 'Book Value ($)', position: 'right' as const,
                                 axisLabel: { color: textColor, formatter: (val: number) => abbreviateNumber(val) },
                                 nameTextStyle: { color: textColor }, axisLine: { lineStyle: { color: axisColor } }, splitLine: { show: false } // No grid lines
                             }] : [])
                        ],
                        series: series
                    };
                    break;
                }
                // Add cases for 'cashflow' etc. if needed
                default:
                    newOptions = { title: { text: 'Select Calculation Type', left: 'center' } };
            }
        } catch (error) {
            console.error("Calculation Error:", error);
            newOptions = { title: { text: `Error: ${error instanceof Error ? error.message : 'Unknown error'}`, left: 'center', textStyle: { color: 'red' } } };
            newResults.error = error instanceof Error ? error.message : 'Unknown error';
        }

        setEchartsOptions(newOptions);
        setResults(newResults);
    }, [selectedCalculation, ciPV, ciYears, ciCompoundsPerYear, ciInterestRate, ciCompareSimple, earCompounds, earInterestRate, npvPV, npvYears, npvDiscountRate, annuityPV, annuityYears, annuityInterestRate, perpCashFlow, perpInterestRate, infPV, infYears, infInterestRate, infInflationRate, depFCI, depLife, depSalvage, depMethod, depCompare]);

    // Update chart when inputs change
    useEffect(() => {
        updateChartAndResults();
    }, [updateChartAndResults]); // Dependency array includes the callback itself

    // --- Render Logic ---
    const renderInputs = () => {
        switch (selectedCalculation) {
            case 'compounding_interest':
                return (
                    <div className="space-y-4">
                        <div><Label>Present Value:</Label><Input type="number" value={ciPV} onChange={e => setCiPV(Number(e.target.value))} /></div>
                        <div><Label>Number of Years:</Label><Input type="number" value={ciYears} min={1} step={1} onChange={e => setCiYears(Number(e.target.value))} /></div>
                        <div><Label>Compounds per Year:</Label><Input type="number" value={ciCompoundsPerYear} min={1} step={1} onChange={e => setCiCompoundsPerYear(Number(e.target.value))} /></div>
                        <div className="mb-2"><Label>Interest Rate: {ciInterestRate.toFixed(1)}%</Label></div> {/* Added mb-2 */}
                        <Slider min={0} max={100} step={0.1} value={[ciInterestRate]} onValueChange={([val]) => setCiInterestRate(val)} />
                        <div className="flex items-center space-x-2"><Checkbox id="ciCompare" checked={ciCompareSimple} onCheckedChange={(checked) => setCiCompareSimple(Boolean(checked))} /><Label htmlFor="ciCompare">Compare with Simple Interest</Label></div>
                        <div className="text-sm text-muted-foreground space-y-1">
                            <p>{results.fvCompoundText}</p>
                            {ciCompareSimple && <p>{results.fvSimpleText}</p>}
                            {ciCompareSimple && <p>{results.differenceText}</p>}
                        </div>
                    </div>
                );
            case 'ear':
                 return (
                    <div className="space-y-4">
                        <div><Label>Compounds per Year:</Label><Input type="number" value={earCompounds} min={1} step={1} onChange={e => setEarCompounds(Number(e.target.value))} /></div>
                        <div className="mb-2"><Label>Annual Interest Rate: {earInterestRate.toFixed(1)}%</Label></div>
                        <Slider min={0} max={100} step={0.1} value={[earInterestRate]} onValueChange={([val]) => setEarInterestRate(val)} />
                        <div className="text-sm text-muted-foreground space-y-1">
                            {/* Always show max EAR results */}
                            <p>{results.earText}</p>
                            <p>{results.maxEarText}</p>
                            <p>{results.differenceText}</p>
                        </div>
                    </div>
                );
            case 'npv':
                 return (
                    <div className="space-y-4">
                        <div><Label>Future Value (at end of term):</Label><Input type="number" value={npvPV} onChange={e => setNpvPV(Number(e.target.value))} /></div>
                        <div><Label>Number of Years:</Label><Input type="number" value={npvYears} min={1} step={1} onChange={e => setNpvYears(Number(e.target.value))} /></div>
                        <div className="mb-2"><Label>Discount Rate: {npvDiscountRate.toFixed(1)}%</Label></div> {/* Added mb-2 */}
                        <Slider min={0} max={100} step={0.1} value={[npvDiscountRate]} onValueChange={([val]) => setNpvDiscountRate(val)} />
                        <div className="text-sm text-muted-foreground space-y-1">
                            <p>{results.fvText}</p> {/* Shows Discounted Value */}
                        </div>
                    </div>
                );
            case 'annuity':
                 return (
                    <div className="space-y-4">
                        <div><Label>Present Value (Loan Amount):</Label><Input type="number" value={annuityPV} onChange={e => setAnnuityPV(Number(e.target.value))} /></div>
                        <div><Label>Number of Years:</Label><Input type="number" value={annuityYears} min={1} step={1} onChange={e => setAnnuityYears(Number(e.target.value))} /></div>
                        <div className="mb-2"><Label>Annual Interest Rate: {annuityInterestRate.toFixed(1)}%</Label></div> {/* Added mb-2 */}
                        <Slider min={0} max={100} step={0.1} value={[annuityInterestRate]} onValueChange={([val]) => setAnnuityInterestRate(val)} />
                        <div className="text-sm text-muted-foreground space-y-1">
                            <p>{results.annuityPaymentText}</p>
                        </div>
                    </div>
                );
            case 'perpetuity':
                 return (
                    <div className="space-y-4">
                        <div><Label>Uniform Cash Flow (per period):</Label><Input type="number" value={perpCashFlow} onChange={e => setPerpCashFlow(Number(e.target.value))} /></div>
                        <div className="mb-2"><Label>Interest Rate (per period): {perpInterestRate.toFixed(1)}%</Label></div> {/* Added mb-2 */}
                        <Slider min={0.1} max={100} step={0.1} value={[perpInterestRate]} onValueChange={([val]) => setPerpInterestRate(val)} />
                        <div className="text-sm text-muted-foreground space-y-1">
                            <p>{results.pvText}</p>
                        </div>
                    </div>
                );
            case 'inflation':
                 return (
                    <div className="space-y-4">
                        <div><Label>Present Value:</Label><Input type="number" value={infPV} onChange={e => setInfPV(Number(e.target.value))} /></div>
                        <div><Label>Number of Years:</Label><Input type="number" value={infYears} min={1} step={1} onChange={e => setInfYears(Number(e.target.value))} /></div>
                        <div className="mb-2"><Label>Nominal Interest Rate: {infInterestRate.toFixed(1)}%</Label></div> {/* Added mb-2 */}
                        <Slider min={0} max={100} step={0.1} value={[infInterestRate]} onValueChange={([val]) => setInfInterestRate(val)} />
                        <div className="mb-2"><Label>Inflation Rate: {infInflationRate.toFixed(1)}%</Label></div> {/* Added mb-2 */}
                        <Slider min={0} max={100} step={0.1} value={[infInflationRate]} onValueChange={([val]) => setInfInflationRate(val)} />
                        <div className="text-sm text-muted-foreground space-y-1">
                            <p>{results.fppText}</p>
                        </div>
                    </div>
                );
            case 'depreciation':
                 return (
                    <div className="space-y-4">
                        <div><Label>Fixed Capital Investment (FCI):</Label><Input type="number" value={depFCI} onChange={e => setDepFCI(Number(e.target.value))} /></div>
                        <div><Label>Equipment Life (Years):</Label><Input type="number" value={depLife} min={1} step={1} onChange={e => setDepLife(Number(e.target.value))} /></div>
                        <div><Label>Salvage Value:</Label><Input type="number" value={depSalvage} min={0} onChange={e => setDepSalvage(Number(e.target.value))} /></div>
                        <div>
                            <Label className="mb-2 block">Depreciation Method:</Label>
                            <RadioGroup value={depMethod} onValueChange={(value: 'straight' | 'macrs') => setDepMethod(value)} className="flex space-x-4">
                                <div className="flex items-center space-x-2"><RadioGroupItem value="straight" id="depStraight" /><Label htmlFor="depStraight">Straight Line</Label></div>
                                <div className="flex items-center space-x-2"><RadioGroupItem value="macrs" id="depMACRS" /><Label htmlFor="depMACRS">MACRS (DDB w/ Switch)</Label></div>
                            </RadioGroup>
                        </div>
                        <div className="flex items-center space-x-2"><Checkbox id="depCompare" checked={depCompare} onCheckedChange={(checked) => setDepCompare(Boolean(checked))} /><Label htmlFor="depCompare">Compare Methods</Label></div>
                         <div className="text-sm text-muted-foreground space-y-1">
                            {(depMethod === 'straight' || depCompare) && <p>{results.totalSLDep}</p>}
                            {(depMethod === 'macrs' || depCompare) && <p>{results.totalMACRSDep}</p>}
                        </div>
                    </div>
                );
            default:
                return <p>Select a calculation type.</p>;
        }
    };

    return (
        <div className="container mx-auto p-4 md:p-8">
            <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                {/* Controls Column */}
                <div className="lg:col-span-1">
                    <Card>
                        <CardContent className="pt-6"> {/* Added padding-top */}
                            {/* Moved Select dropdown here */}
                            <div className="mb-6"> {/* Added margin-bottom */}
                                <Select value={selectedCalculation} onValueChange={(value) => setSelectedCalculation(value as CalculationType)}>
                                    <SelectTrigger className="w-full">
                                        <SelectValue placeholder="Select Calculation Type" />
                                    </SelectTrigger>
                                    <SelectContent>
                                        <SelectItem value="compounding_interest">Compounding Interest</SelectItem>
                                        <SelectItem value="ear">Effective Annual Rate (EAR)</SelectItem>
                                        <SelectItem value="npv">Discounted Value (PV of FV)</SelectItem>
                                        <SelectItem value="annuity">Annuity Payment (A from P)</SelectItem>
                                        <SelectItem value="perpetuity">Perpetuity (PV from A)</SelectItem>
                                        <SelectItem value="inflation">Inflationary Effects</SelectItem>
                                        <SelectItem value="depreciation">Depreciation</SelectItem>
                                        {/* <SelectItem value="cashflow">Cashflow Visualization</SelectItem> */}
                                    </SelectContent>
                                </Select>
                            </div>
                            {renderInputs()}
                            {results.error && <p className="text-red-600 mt-4">Error: {results.error}</p>}
                        </CardContent>
                    </Card>
                </div>

                {/* Plot Column */}
                <div className="lg:col-span-2">
                    <Card>
                        <CardContent className="pt-6">
                            <div className="relative h-[450px] md:h-[550px] rounded-md overflow-hidden border" style={{ backgroundColor: '#08306b' }}>
                                {Object.keys(echartsOptions).length > 0 ? (
                                    <ReactECharts
                                        echarts={echarts}
                                        option={echartsOptions}
                                        style={{ height: '100%', width: '100%' }}
                                        notMerge={true}
                                        lazyUpdate={false}
                                    />
                                ) : (
                                    <div className="absolute inset-0 flex items-center justify-center text-muted-foreground">
                                        Chart will appear here.
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
