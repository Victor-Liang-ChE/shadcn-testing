"use client";
import { useState, useEffect, useCallback } from 'react';
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsOption, SeriesOption } from 'echarts';
import { LineChart, BarChart } from 'echarts/charts';
import { TitleComponent, TooltipComponent, GridComponent, LegendComponent, ToolboxComponent, GraphicComponent } from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';
import { Card, CardContent } from "@/components/ui/card";
import { Label } from "@/components/ui/label";
import { Slider } from "@/components/ui/slider";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Input } from "@/components/ui/input";
import { Checkbox } from "@/components/ui/checkbox";
import { RadioGroup, RadioGroupItem } from "@/components/ui/radio-group";

// ECharts imports
echarts.use([
    TitleComponent,
    TooltipComponent,
    GridComponent,
    LegendComponent,
    ToolboxComponent,
    GraphicComponent,
    LineChart,
    BarChart,
    CanvasRenderer
]);

// --- Type Definitions ---
type CalculationType =
    | 'compounding_interest'
    | 'ear'
    | 'npv'
    | 'inflation'
    | 'annuity'
    | 'perpetuity'
    | 'depreciation'
    | 'cashflow';

// --- Helper Functions ---

// Abbreviate numbers
function abbreviateNumber(num: number): string {
    if (num >= 1e93) return '∞';
    const suffixes = [
        { value: 1e93, symbol: 'Tg' }, { value: 1e90, symbol: 'Nvg' }, { value: 1e87, symbol: 'Ovg' },
        { value: 1e84, symbol: 'Spvg' }, { value: 1e81, symbol: 'Sxvg' }, { value: 1e78, symbol: 'Qvg' },
        { value: 1e75, symbol: 'Qavg' }, { value: 1e72, symbol: 'Qav' }, { value: 1e69, symbol: 'Dvg' },
        { value: 1e66, symbol: 'Uvg' }, { value: 1e63, symbol: 'Vg' }, { value: 1e60, symbol: 'Nvg' },
        { value: 1e57, symbol: 'Ocdc' }, { value: 1e54, symbol: 'Spdc' },
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
    const [ciPVInput, setCiPVInput] = useState<string>('1000');
    const [ciYears, setCiYears] = useState<number>(10);
    const [ciYearsInput, setCiYearsInput] = useState<string>('10');
    const [ciCompoundsPerYear, setCiCompoundsPerYear] = useState<number>(12);
    const [ciCompoundsInput, setCiCompoundsInput] = useState<string>('12');
    const [ciInterestRate, setCiInterestRate] = useState<number>(5);
    const [ciCompareSimple, setCiCompareSimple] = useState<boolean>(true);

    // EAR
    const [earCompounds, setEarCompounds] = useState<number>(1);
    const [earCompoundsInput, setEarCompoundsInput] = useState<string>('1');
    const [earInterestRate, setEarInterestRate] = useState<number>(5);

    // PV of a Future Sum (formerly NPV)
    const [npvPV, setNpvPV] = useState<number>(1000);
    const [npvPVInput, setNpvPVInput] = useState<string>('1000');
    const [npvYears, setNpvYears] = useState<number>(1);
    const [npvYearsInput, setNpvYearsInput] = useState<string>('1');
    const [npvDiscountRate, setNpvDiscountRate] = useState<number>(5);

    // Annuity
    const [annuityPV, setAnnuityPV] = useState<number>(1000);
    const [annuityPVInput, setAnnuityPVInput] = useState<string>('1000');
    const [annuityYears, setAnnuityYears] = useState<number>(10);
    const [annuityYearsInput, setAnnuityYearsInput] = useState<string>('10');
    const [annuityInterestRate, setAnnuityInterestRate] = useState<number>(5);

    // Perpetuity
    const [perpCashFlow, setPerpCashFlow] = useState<number>(1000);
    const [perpCashFlowInput, setPerpCashFlowInput] = useState<string>('1000');
    const [perpInterestRate, setPerpInterestRate] = useState<number>(5);

    // Inflation
    const [infPV, setInfPV] = useState<number>(1000);
    const [infPVInput, setInfPVInput] = useState<string>('1000');
    const [infYears, setInfYears] = useState<number>(10);
    const [infYearsInput, setInfYearsInput] = useState<string>('10');
    const [infInterestRate, setInfInterestRate] = useState<number>(5);
    const [infInflationRate, setInfInflationRate] = useState<number>(2);

    // Depreciation
    const [depFCI, setDepFCI] = useState<number>(10000);
    const [depFCIInput, setDepFCIInput] = useState<string>('10000');
    const [depLife, setDepLife] = useState<number>(10);
    const [depLifeInput, setDepLifeInput] = useState<string>('10');
    const [depSalvage, setDepSalvage] = useState<number>(1000);
    const [depSalvageInput, setDepSalvageInput] = useState<string>('1000');
    const [depMethod, setDepMethod] = useState<'straight' | 'macrs'>('straight');
    const [depCompare, setDepCompare] = useState<boolean>(false);

    // --- Calculation and Plotting Logic ---
    
    // Helper function for consistent chart styling
    const createBaseChartOptions = (title: string, xAxisName: string, yAxisName: string, xMax?: number, yMax?: number): Partial<EChartsOption> => ({
        animation: false,
        backgroundColor: 'transparent',
        title: {
            text: title,
            left: 'center',
            textStyle: {
                color: 'white',
                fontSize: 18,
                fontFamily: 'Merriweather Sans'
            }
        },
        grid: { left: '5%', right: '5%', bottom: '12%', top: '5%', containLabel: true },
        xAxis: {
            type: 'value',
            name: xAxisName,
            min: 0,
            max: xMax,
            nameLocation: 'middle',
            nameGap: 50,
            nameTextStyle: {
                color: 'white',
                fontSize: 15,
                fontFamily: 'Merriweather Sans'
            },
            axisLine: { lineStyle: { color: 'white' } },
            axisTick: { lineStyle: { color: 'white' }, length: 5, inside: false },
            axisLabel: {
                color: 'white',
                fontSize: 16,
                fontFamily: 'Merriweather Sans'
            },
            splitLine: { show: false }
        },
        yAxis: {
            type: 'value',
            name: yAxisName,
            min: 0,
            max: yMax,
            nameLocation: 'middle',
            nameGap: 80,
            nameTextStyle: {
                color: 'white',
                fontSize: 15,
                fontFamily: 'Merriweather Sans'
            },
            axisLine: { lineStyle: { color: 'white' } },
            axisTick: { lineStyle: { color: 'white' }, length: 5, inside: false },
            axisLabel: {
                color: 'white',
                fontSize: 16,
                fontFamily: 'Merriweather Sans'
            },
            splitLine: { show: false }
        },
        legend: {
            bottom: '5%',
            textStyle: {
                color: 'white',
                fontSize: 12,
                fontFamily: 'Merriweather Sans'
            }
        },
        tooltip: { trigger: 'axis' }
    });

    const updateChartAndResults = useCallback(() => {
        let newOptions: EChartsOption = {};
        let newResults: Record<string, string | number> = {};
        const lineColors = ['#ffff00', '#ff0000', '#00ff00', '#00ffff'];

        // Validation function to check if inputs are valid for calculation
        const validateInputs = () => {
            switch (selectedCalculation) {
                case 'compounding_interest':
                    return ciPV > 0 && ciYears > 0 && ciCompoundsPerYear > 0;
                case 'ear':
                    return earCompounds > 0 && earInterestRate >= 0;
                case 'npv':
                    return npvPV > 0 && npvYears > 0 && npvDiscountRate >= 0;
                case 'annuity':
                    return annuityPV > 0 && annuityYears > 0 && annuityInterestRate >= 0;
                case 'perpetuity':
                    return perpCashFlow > 0 && perpInterestRate > 0;
                case 'inflation':
                    return infPV > 0 && infYears > 0 && infInterestRate >= 0 && infInflationRate >= 0;
                case 'depreciation':
                    return depFCI > 0 && depLife > 0 && depSalvage >= 0 && depFCI >= depSalvage;
                default:
                    return false;
            }
        };

        // If inputs are not valid, show a default empty chart
        if (!validateInputs()) {
            newOptions = {
                backgroundColor: 'transparent',
                title: { 
                    text: 'Enter valid inputs to see results', 
                    left: 'center', 
                    textStyle: { 
                        color: 'white', 
                        fontSize: 16, 
                        fontFamily: 'Merriweather Sans' 
                    } 
                }
            };
            setEchartsOptions(newOptions);
            setResults({});
            return;
        }

        try {
            switch (selectedCalculation) {
                case 'compounding_interest': {
                    const pv = ciPV || 0;
                    const n = ciYears || 1;
                    const m = ciCompoundsPerYear || 1;
                    const r = ciInterestRate / 100;
                    if (n <= 0 || m <= 0) {
                        // Use fallback values instead of throwing error
                        return;
                    }

                    const totalPeriods = n * m;
                    const periodRate = r / m;
                    const fvCompound = pv * Math.pow(1 + periodRate, totalPeriods);
                    const fvSimple = pv * (1 + n * r);

                    newResults.fvCompoundText = `Future Value: $${abbreviateNumber(fvCompound)}`;
                    if (ciCompareSimple) {
                        newResults.fvSimpleText = `Future Value (Simple): $${abbreviateNumber(fvSimple)}`;
                        newResults.differenceText = `Difference: $${abbreviateNumber(fvCompound - fvSimple)}`;
                    }

                    const x_values = Array.from({ length: totalPeriods + 1 }, (_, x) => x / m);
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

                    const baseOptions = createBaseChartOptions('Future Value vs. Years', 'Years', 'Future Value ($)', n);
                    const legendData = series.map(s => s.name).filter((name): name is string => !!name);

                    newOptions = {
                        ...baseOptions,
                        tooltip: { trigger: 'axis', valueFormatter: (value) => `$${formatNumber(value as number)}` },
                        legend: { 
                            ...baseOptions.legend, 
                            data: legendData 
                        },
                        yAxis: {
                            ...baseOptions.yAxis,
                            axisLabel: {
                                color: 'white',
                                fontSize: 16,
                                fontFamily: 'Merriweather Sans',
                                formatter: (val: number) => abbreviateNumber(val)
                            }
                        },
                        series: series
                    };
                    break;
                }
                case 'ear': {
                    const m = earCompounds;
                    const r = earInterestRate / 100;

                    const ear = Math.pow(1 + r / m, m) - 1;
                    const maxEar = Math.exp(r) - 1;

                    newResults.earText = `EAR: ${formatNumber(ear * 100)}%`;
                    newResults.maxEarText = `Max EAR: ${formatNumber(maxEar * 100)}%`;
                    newResults.differenceText = `Difference: ${formatNumber((maxEar - ear) * 100)}%`;

                    const compoundsRange = Array.from({ length: Math.max(m, 50) + 1 }, (_, i) => i).filter(i => i > 0);
                    const ears = compoundsRange.map(comp => Math.pow(1 + r / comp, comp) - 1);

                    const series: SeriesOption[] = [{
                        name: 'EAR', type: 'line', data: compoundsRange.map((comp, i) => [comp, ears[i] * 100]),
                        lineStyle: { color: lineColors[0] }, symbol: 'none', smooth: true
                    }];

                    const shapes = [{
                        type: 'line' as const, x0: 1, y0: maxEar * 100, x1: Math.max(m, 50), y1: maxEar * 100,
                        lineStyle: { color: lineColors[1], type: 'dashed' as const }
                    }];

                    const baseOptions = createBaseChartOptions('EAR vs. Number of Compounds', 'Compounds per Year', 'EAR (%)');
                    const legendData: string[] = ['EAR', 'Max EAR'];

                    newOptions = {
                        ...baseOptions,
                        tooltip: { trigger: 'axis', valueFormatter: (value) => `${formatNumber(value as number)}%` },
                        legend: { 
                            ...baseOptions.legend, 
                            data: legendData 
                        },
                        xAxis: {
                            ...baseOptions.xAxis,
                            min: 1
                        },
                        yAxis: {
                            ...baseOptions.yAxis,
                            nameGap: 60, // Keep original spacing for EAR chart
                            axisLabel: {
                                color: 'white',
                                fontSize: 16,
                                fontFamily: 'Merriweather Sans',
                                formatter: '{value}%'
                            }
                        },
                        series: series,
                        graphic: { elements: shapes }
                    };
                    break;
                }
                case 'npv': {
                    const fvAmount = npvPV; // Use more descriptive variable name
                    const n = npvYears;
                    const r = npvDiscountRate / 100;

                    // This formula calculates the Present Value (PV) of a Future Value (FV)
                    const pv = fvAmount / Math.pow(1 + r, n);
                    newResults.fvText = `Present Value: $${formatNumber(pv)}`; // Clarify result text

                    const yearsList = Array.from({ length: n + 1 }, (_, i) => i);
                    const discountedValues = yearsList.map(year => fvAmount / Math.pow(1 + r, year));

                    const baseOptions = createBaseChartOptions('Present Value of a Future Sum vs. Years', 'Years', 'Present Value ($)', n);
                    const seriesName = 'Present Value';
                    
                    newOptions = {
                        ...baseOptions,
                        tooltip: { trigger: 'axis', valueFormatter: (value) => `$${formatNumber(value as number)}` },
                        legend: { 
                            ...baseOptions.legend, 
                            data: [seriesName] 
                        },
                        yAxis: {
                            ...baseOptions.yAxis,
                            axisLabel: {
                                color: 'white',
                                fontSize: 16,
                                fontFamily: 'Merriweather Sans',
                                formatter: (val: number) => abbreviateNumber(val)
                            }
                        },
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

                    const annuityPayment = (r === 0) ? pv / n : pv * (r * Math.pow(1 + r, n)) / (Math.pow(1 + r, n) - 1);
                    newResults.annuityPaymentText = `Required Annuity Payment: $${formatNumber(annuityPayment)}`;

                    const yearsRange = Array.from({ length: n }, (_, i) => i + 1);
                    const payments = yearsRange.map(yr => (r === 0) ? pv / yr : pv * (r * Math.pow(1 + r, yr)) / (Math.pow(1 + r, yr) - 1));

                    const baseOptions = createBaseChartOptions('Required Annuity Payment vs. Years', 'Years', 'Annuity Payment ($)', n);
                    const seriesName = 'Annuity Payment';
                    
                    newOptions = {
                        ...baseOptions,
                        tooltip: { trigger: 'axis', valueFormatter: (value) => `$${formatNumber(value as number)}` },
                        legend: { 
                            ...baseOptions.legend, 
                            data: [seriesName] 
                        },
                        xAxis: {
                            ...baseOptions.xAxis,
                            min: 1
                        },
                        yAxis: {
                            ...baseOptions.yAxis,
                            axisLabel: {
                                color: 'white',
                                fontSize: 16,
                                fontFamily: 'Merriweather Sans',
                                formatter: (val: number) => abbreviateNumber(val)
                            }
                        },
                        series: [{
                            name: seriesName, type: 'line', data: yearsRange.map((yr, i) => [yr, payments[i]]),
                            lineStyle: { color: lineColors[0] }, symbol: 'none', smooth: true
                        }]
                    };
                    break;
                }
                case 'perpetuity': {
                    const a = perpCashFlow;
                    const r = perpInterestRate / 100;

                    const pv = a / r;
                    newResults.pvText = `Present Value: $${formatNumber(pv)}`;

                    const cashFlowRange = Array.from({ length: 100 }, (_, i) => (i + 1) * (a / 10));
                    const pvValues = cashFlowRange.map(cf => cf / r);

                    const baseOptions = createBaseChartOptions('Present Value vs. Uniform Cash Flow', 'Uniform Cash Flow ($)', 'Present Value ($)');
                    const seriesName = 'Present Value';
                    
                    newOptions = {
                        ...baseOptions,
                        tooltip: { trigger: 'axis', valueFormatter: (value) => `$${formatNumber(value as number)}` },
                        legend: { 
                            ...baseOptions.legend, 
                            data: [seriesName] 
                        },
                        xAxis: {
                            ...baseOptions.xAxis,
                            axisLabel: {
                                color: 'white',
                                fontSize: 16,
                                fontFamily: 'Merriweather Sans',
                                formatter: (val: number) => abbreviateNumber(val)
                            }
                        },
                        yAxis: {
                            ...baseOptions.yAxis,
                            axisLabel: {
                                color: 'white',
                                fontSize: 16,
                                fontFamily: 'Merriweather Sans',
                                formatter: (val: number) => abbreviateNumber(val)
                            }
                        },
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
                    const i = infInflationRate / 100;

                    const realRateFactor = (1 + r) / (1 + i) - 1;
                    const fpp = pv * Math.pow(1 + realRateFactor, n);
                    newResults.fppText = `Future Purchasing Power: $${formatNumber(fpp)}`;

                    const yearsList = Array.from({ length: n + 1 }, (_, y) => y);
                    const fppList = yearsList.map(yr => pv * Math.pow(1 + realRateFactor, yr));

                    const baseOptions = createBaseChartOptions('Future Purchasing Power vs. Years', 'Years', 'Future Purchasing Power ($)', n);
                    const seriesName = 'Future Purchasing Power';
                    
                    newOptions = {
                        ...baseOptions,
                        tooltip: { trigger: 'axis', valueFormatter: (value) => `$${formatNumber(value as number)}` },
                        legend: { 
                            ...baseOptions.legend, 
                            data: [seriesName] 
                        },
                        yAxis: {
                            ...baseOptions.yAxis,
                            axisLabel: {
                                color: 'white',
                                fontSize: 16,
                                fontFamily: 'Merriweather Sans',
                                formatter: (val: number) => abbreviateNumber(val)
                            }
                        },
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

                    const years = Array.from({ length: life + 1 }, (_, i) => i);

                    // --- Straight Line Calculation ---
                    const slDepPerYear = (fci - salvage) / life;
                    const slBookValue = years.map(yr => Math.max(salvage, fci - yr * slDepPerYear));

                    // --- Corrected MACRS (DDB with Switch to SL) Calculation ---
                    const macrsAnnualDep: number[] = [0];
                    const macrsBookValue: number[] = [fci];
                    let currentBookValue = fci;
                    const ddbRate = 2 / life;
                    let hasSwitchedToSL = false;

                    for (let yr = 1; yr <= life; yr++) {
                        const remainingLife = life - yr + 1;
                        const slDepFromCurrent = remainingLife > 0 ? (currentBookValue - salvage) / remainingLife : 0;
                        const ddbDepreciation = currentBookValue * ddbRate;
                        
                        let annualDepreciation = 0;

                        if (!hasSwitchedToSL && slDepFromCurrent > ddbDepreciation) {
                            hasSwitchedToSL = true;
                        }
                        
                        annualDepreciation = hasSwitchedToSL ? slDepFromCurrent : ddbDepreciation;
                        annualDepreciation = Math.max(0, Math.min(annualDepreciation, currentBookValue - salvage));
                        
                        macrsAnnualDep.push(annualDepreciation);
                        currentBookValue -= annualDepreciation;
                        macrsBookValue.push(currentBookValue);
                    }

                    newResults.totalSLDep = `Total SL Dep: $${formatNumber(fci - salvage)}`;
                    newResults.totalMACRSDep = `Total MACRS Dep: $${formatNumber(macrsAnnualDep.reduce((a, b) => a + b, 0))}`;

                    const series: SeriesOption[] = [];
                    const legendData: string[] = [];

                    if (depMethod === 'straight' || depCompare) {
                        series.push({
                            name: 'Straight Line (Annual)', type: 'bar', data: years.map((yr) => [yr, yr > 0 ? slDepPerYear : 0]),
                            itemStyle: { color: lineColors[1] }
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
                             lineStyle: { color: lineColors[1], type: 'dashed' }, symbol: 'none', yAxisIndex: 1
                         });
                         series.push({
                             name: 'MACRS Book Value', type: 'line', data: years.map((yr, i) => [yr, macrsBookValue[i]]),
                             lineStyle: { color: lineColors[0], type: 'dashed' }, symbol: 'none', yAxisIndex: 1
                         });
                         legendData.push('SL Book Value', 'MACRS Book Value');
                     }

                    const baseOptions = createBaseChartOptions('Annual Depreciation vs. Year', 'Year', 'Annual Depreciation ($)', life);

                    newOptions = {
                        ...baseOptions,
                        tooltip: { 
                            trigger: 'axis', 
                            valueFormatter: (value) => `$${formatNumber(value as number)}`,
                            axisPointer: { type: 'shadow' }
                        },
                        legend: { 
                            ...baseOptions.legend, 
                            data: legendData 
                        },
                        grid: { left: '5%', right: depCompare ? '15%' : '5%', bottom: '12%', top: '5%', containLabel: true },
                        yAxis: [
                            {
                                ...baseOptions.yAxis,
                                axisLabel: {
                                    color: 'white',
                                    fontSize: 16,
                                    fontFamily: 'Merriweather Sans',
                                    formatter: (val: number) => abbreviateNumber(val)
                                }
                            },
                             ...(depCompare ? [{
                                 type: 'value' as const, 
                                 name: 'Book Value ($)', 
                                 position: 'right' as const,
                                 nameLocation: 'middle' as const,
                                 nameGap: 60,
                                 nameTextStyle: {
                                     color: 'white',
                                     fontSize: 15,
                                     fontFamily: 'Merriweather Sans'
                                 },
                                 axisLine: { lineStyle: { color: 'white' } },
                                 axisTick: { lineStyle: { color: 'white' }, length: 5, inside: false },
                                 axisLabel: { 
                                     color: 'white', 
                                     fontSize: 16, 
                                     fontFamily: 'Merriweather Sans',
                                     formatter: (val: number) => abbreviateNumber(val) 
                                 },
                                 splitLine: { show: false }
                             }] : [])
                        ],
                        series: series
                    };
                    break;
                }
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

    useEffect(() => {
        updateChartAndResults();
    }, [updateChartAndResults]);

    const renderInputs = () => {
        switch (selectedCalculation) {
            case 'compounding_interest':
                return (
                    <div className="space-y-4">
                        <div><Label>Present Value:</Label><Input type="number" value={ciPVInput} onChange={e => {
                            setCiPVInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val > 0) setCiPV(val);
                        }} /></div>
                        <div><Label>Number of Years:</Label><Input type="number" value={ciYearsInput} min={1} step={1} onChange={e => {
                            setCiYearsInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val > 0) setCiYears(val);
                        }} /></div>
                        <div><Label>Compounds per Year:</Label><Input type="number" value={ciCompoundsInput} min={1} step={1} onChange={e => {
                            setCiCompoundsInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val > 0) setCiCompoundsPerYear(val);
                        }} /></div>
                        <div className="mb-2"><Label>Interest Rate: {ciInterestRate.toFixed(1)}%</Label></div>
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
                        <div><Label>Compounds per Year:</Label><Input type="number" value={earCompoundsInput} min={1} step={1} onChange={e => {
                            setEarCompoundsInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val > 0) setEarCompounds(val);
                        }} /></div>
                        <div className="mb-2"><Label>Annual Interest Rate: {earInterestRate.toFixed(1)}%</Label></div>
                        <Slider min={0} max={100} step={0.1} value={[earInterestRate]} onValueChange={([val]) => setEarInterestRate(val)} />
                        <div className="text-sm text-muted-foreground space-y-1">
                            <p>{results.earText}</p>
                            <p>{results.maxEarText}</p>
                            <p>{results.differenceText}</p>
                        </div>
                    </div>
                );
            case 'npv':
                 return (
                    <div className="space-y-4">
                        <div><Label>Future Sum (Value at end of term):</Label><Input type="number" value={npvPVInput} onChange={e => {
                            setNpvPVInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val > 0) setNpvPV(val);
                        }} /></div>
                        <div><Label>Number of Years:</Label><Input type="number" value={npvYearsInput} min={1} step={1} onChange={e => {
                            setNpvYearsInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val > 0) setNpvYears(val);
                        }} /></div>
                        <div className="mb-2"><Label>Discount Rate: {npvDiscountRate.toFixed(1)}%</Label></div>
                        <Slider min={0} max={100} step={0.1} value={[npvDiscountRate]} onValueChange={([val]) => setNpvDiscountRate(val)} />
                        <div className="text-sm text-muted-foreground space-y-1">
                            <p>{results.fvText}</p>
                        </div>
                    </div>
                );
            case 'annuity':
                 return (
                    <div className="space-y-4">
                        <div><Label>Present Value (Loan Amount):</Label><Input type="number" value={annuityPVInput} onChange={e => {
                            setAnnuityPVInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val > 0) setAnnuityPV(val);
                        }} /></div>
                        <div><Label>Number of Years:</Label><Input type="number" value={annuityYearsInput} min={1} step={1} onChange={e => {
                            setAnnuityYearsInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val > 0) setAnnuityYears(val);
                        }} /></div>
                        <div className="mb-2"><Label>Annual Interest Rate: {annuityInterestRate.toFixed(1)}%</Label></div>
                        <Slider min={0} max={100} step={0.1} value={[annuityInterestRate]} onValueChange={([val]) => setAnnuityInterestRate(val)} />
                        <div className="text-sm text-muted-foreground space-y-1">
                            <p>{results.annuityPaymentText}</p>
                        </div>
                    </div>
                );
            case 'perpetuity':
                 return (
                    <div className="space-y-4">
                        <div><Label>Uniform Cash Flow (per period):</Label><Input type="number" value={perpCashFlowInput} onChange={e => {
                            setPerpCashFlowInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val >= 0) setPerpCashFlow(val);
                        }} /></div>
                        <div className="mb-2"><Label>Interest Rate (per period): {perpInterestRate.toFixed(1)}%</Label></div>
                        <Slider min={0.1} max={100} step={0.1} value={[perpInterestRate]} onValueChange={([val]) => setPerpInterestRate(val)} />
                        <div className="text-sm text-muted-foreground space-y-1">
                            <p>{results.pvText}</p>
                        </div>
                    </div>
                );
            case 'inflation':
                 return (
                    <div className="space-y-4">
                        <div><Label>Present Value:</Label><Input type="number" value={infPVInput} onChange={e => {
                            setInfPVInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val >= 0) setInfPV(val);
                        }} /></div>
                        <div><Label>Number of Years:</Label><Input type="number" value={infYearsInput} min={1} step={1} onChange={e => {
                            setInfYearsInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val >= 1) setInfYears(val);
                        }} /></div>
                        <div className="mb-2"><Label>Nominal Interest Rate: {infInterestRate.toFixed(1)}%</Label></div>
                        <Slider min={0} max={100} step={0.1} value={[infInterestRate]} onValueChange={([val]) => setInfInterestRate(val)} />
                        <div className="mb-2"><Label>Inflation Rate: {infInflationRate.toFixed(1)}%</Label></div>
                        <Slider min={0} max={100} step={0.1} value={[infInflationRate]} onValueChange={([val]) => setInfInflationRate(val)} />
                        <div className="text-sm text-muted-foreground space-y-1">
                            <p>{results.fppText}</p>
                        </div>
                    </div>
                );
            case 'depreciation':
                 return (
                    <div className="space-y-4">
                        <div><Label>Fixed Capital Investment (FCI):</Label><Input type="number" value={depFCIInput} onChange={e => {
                            setDepFCIInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val >= 0) setDepFCI(val);
                        }} /></div>
                        <div><Label>Equipment Life (Years):</Label><Input type="number" value={depLifeInput} min={1} step={1} onChange={e => {
                            setDepLifeInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val >= 1) setDepLife(val);
                        }} /></div>
                        <div><Label>Salvage Value:</Label><Input type="number" value={depSalvageInput} min={0} onChange={e => {
                            setDepSalvageInput(e.target.value);
                            const val = Number(e.target.value);
                            if (!isNaN(val) && val >= 0) setDepSalvage(val);
                        }} /></div>
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
                <div className="lg:col-span-1">
                    <Card>
                        <CardContent className="pt-6">
                            <div className="mb-6">
                                <Select value={selectedCalculation} onValueChange={(value) => setSelectedCalculation(value as CalculationType)}>
                                    <SelectTrigger className="w-full">
                                        <SelectValue placeholder="Select Calculation Type" />
                                    </SelectTrigger>
                                    <SelectContent>
                                        <SelectItem value="compounding_interest">Compounding Interest</SelectItem>
                                        <SelectItem value="ear">Effective Annual Rate (EAR)</SelectItem>
                                        <SelectItem value="npv">Present Value of a Future Sum</SelectItem>
                                        <SelectItem value="annuity">Annuity Payment (A from P)</SelectItem>
                                        <SelectItem value="perpetuity">Perpetuity (PV from A)</SelectItem>
                                        <SelectItem value="inflation">Inflationary Effects</SelectItem>
                                        <SelectItem value="depreciation">Depreciation</SelectItem>
                                    </SelectContent>
                                </Select>
                            </div>
                            {renderInputs()}
                            {results.error && <p className="text-red-600 mt-4">Error: {results.error}</p>}
                        </CardContent>
                    </Card>
                </div>

                <div className="lg:col-span-2">
                    <Card>
                        <CardContent>
                            <div className="relative w-full aspect-square rounded-md" style={{ backgroundColor: '#08306b' }}>
                                {Object.keys(echartsOptions).length > 0 ? (
                                    <ReactECharts
                                        echarts={echarts}
                                        option={echartsOptions}
                                        style={{ height: '100%', width: '100%', borderRadius: '0.375rem', overflow: 'hidden' }}
                                        notMerge={true}
                                        lazyUpdate={false}
                                    />
                                ) : (
                                    <div className="absolute inset-0 flex items-center justify-center text-white">
                                        <div className="text-center">
                                            Chart will appear here.
                                        </div>
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