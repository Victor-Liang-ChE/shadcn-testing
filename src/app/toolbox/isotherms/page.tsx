'use client';

import React, { useState, useEffect, useCallback, useRef } from 'react';
import { useTheme } from 'next-themes';
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsOption, SeriesOption } from 'echarts';
import { LineChart, ScatterChart } from 'echarts/charts';
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  MarkLineComponent,
  MarkPointComponent,
  ToolboxComponent
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

echarts.use([
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  MarkLineComponent,
  MarkPointComponent,
  ToolboxComponent,
  LineChart,
  ScatterChart,
  CanvasRenderer
]);

import { Card, CardContent } from "@/components/ui/card";
import { Label } from "@/components/ui/label";
import { Slider } from "@/components/ui/slider";
import { TooltipProvider } from "@/components/ui/tooltip";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Button } from "@/components/ui/button";
// Function to format numbers to avoid long decimals - moved outside component
const formatNumberToPrecision = (num: any, precision: number = 3): string => {
  if (typeof num === 'number') {
    if (num === 0) return '0';
    const fixed = num.toPrecision(precision);
    // Remove trailing zeros after decimal point, but keep integer part if it ends in zero
    if (fixed.includes('.')) {
      return parseFloat(fixed).toString(); 
    }
    return fixed;
  }
  return String(num); // Return as string if not a number
};

export default function LangmuirIsothermPage() {
  const { resolvedTheme } = useTheme();
  
  // Isotherm Type State
  type IsothermType = 'langmuir' | 'freundlich' | 'temkin';
  const [selectedIsotherm, setSelectedIsotherm] = useState<IsothermType>('langmuir');

  // Input States
  // Langmuir
  const [langmuirK, setLangmuirK] = useState<number>(0.5);
  const [independentVar, setIndependentVar] = useState<'pressure' | 'concentration'>('pressure');

  // Freundlich
  const [freundlichKf, setFreundlichKf] = useState<number>(2);
  const [freundlichNf, setFreundlichNf] = useState<number>(2);

  // Temkin
  const [temkinAt, setTemkinAt] = useState<number>(1); // Temkin equilibrium binding constant
  const [temkinBt, setTemkinBt] = useState<number>(2); // Temkin constant related to heat of adsorption

  // Data & Control States
  const [isothermData, setIsothermData] = useState<{ c: number[], q: number[] } | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // State for ECharts options
  const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
  const echartsRef = useRef<ReactECharts | null>(null);

  // Function to calculate Isotherm
  const calculateIsotherm = useCallback(() => {
    setLoading(true);
    setError(null);
    try {
      const cValues: number[] = [];
      const qValues: number[] = [];
      let maxConcentration = 5; // Default max concentration for plotting
      const concentrationPoints = 100; // Fixed number of points

      for (let i = 0; i <= concentrationPoints; i++) {
        let c = (i / concentrationPoints) * maxConcentration;
        if (c === 0 && (selectedIsotherm === 'freundlich' || selectedIsotherm === 'temkin')) {
            c = 1e-12; // Use a very small value instead of zero for better precision
        }

        let q: number;
        switch (selectedIsotherm) {
          case 'langmuir':
            // θ_A = (K_eq * P_A) / (1 + K_eq * P_A)
            // Using c as the independent variable (P or C depending on user choice)
            q = (langmuirK * c) / (1 + langmuirK * c);
            break;
          case 'freundlich':
            q = freundlichKf * Math.pow(c, 1 / freundlichNf);
            break;
          case 'temkin':
            // q = B_temkin * ln(A_temkin * C)
            const temkinArg = temkinAt * c;
            if (temkinArg <= 0) {
                q = 0; // Handle logarithm domain issues
            } else {
                q = temkinBt * Math.log(temkinArg);
                q = Math.max(0, q); // Ensure q never negative
            }
            break;
          default:
            q = 0;
        }
        
        // Ensure q is never negative for all models
        q = Math.max(0, q);
        cValues.push(c);
        qValues.push(q);
      }
      setIsothermData({ c: cValues, q: qValues });
      setLoading(false);
    } catch (err: any) {
      setError(`Calculation error: ${err.message}`);
      setLoading(false);
    }
  }, [
    selectedIsotherm,
    langmuirK, independentVar,
    freundlichKf, freundlichNf,
    temkinAt, temkinBt
  ]);

  // useEffect for initial graph load and when parameters change
  useEffect(() => {
    calculateIsotherm();
  }, [calculateIsotherm]);

  // useEffect to update ECharts options when isothermData changes
  useEffect(() => {
    if (isothermData) {
      // Calculate yAxisMax based on actual data range
      const yAxisMax = Math.max(...isothermData.q) * 1.2;
      
      // Theme-aware colors
      const textColor = resolvedTheme === 'dark' ? '#ffffff' : '#000000';

      // Set chart title based on selected isotherm
      let chartTitle = 'Adsorption Isotherm';
      switch (selectedIsotherm) {
        case 'langmuir':
          chartTitle = 'Langmuir Adsorption Isotherm';
          break;
        case 'freundlich':
          chartTitle = 'Freundlich Adsorption Isotherm';
          break;
        case 'temkin':
          chartTitle = 'Temkin Adsorption Isotherm';
          break;
      }


      const newOptions: EChartsOption = {
        backgroundColor: 'transparent',
        animation: true,
        animationDuration: 300,
        animationEasing: 'cubicOut',
        title: {
          text: chartTitle,
          left: 'center',
          top: '20px',
          textStyle: {
            color: textColor,
            fontSize: 18,
            fontFamily: 'Merriweather Sans'
          }
        },
        grid: {
          left: '5%',
          right: '5%',
          bottom: '5%',
          top: '5%',
          containLabel: true
        },
        xAxis: {
          type: 'value',
          min: 0,
          max: Math.max(...isothermData.c),
          name: independentVar === 'pressure' ? 'Pressure, P (bar)' : 'Concentration, C (mol/L)',
          nameLocation: 'middle',
          nameGap: 30,
          nameTextStyle: {
            color: textColor,
            fontSize: 15,
            fontFamily: 'Merriweather Sans'
          },
          axisLine: {
            lineStyle: { color: textColor }
          },
          axisTick: {
            lineStyle: { color: textColor },
            length: 5,
            inside: false
          },
          axisLabel: {
            showMaxLabel: false, // Hide maximum value label on x-axis
            color: textColor,
            fontSize: 16,
            fontFamily: 'Merriweather Sans',
            formatter: (value: any) => formatNumberToPrecision(value, 3)
          },
          splitLine: { show: false }
        },
        yAxis: {
          type: 'value',
          min: 0,
          max: yAxisMax,
          name: selectedIsotherm === 'langmuir' 
            ? 'Surface Coverage (θ)'
            : 'Amount Adsorbed (q)',
          nameLocation: 'middle',
          nameGap: 60,
          nameTextStyle: {
            color: textColor,
            fontSize: 15,
            fontFamily: 'Merriweather Sans'
          },
          axisLine: {
            lineStyle: { color: textColor }
          },
          axisTick: {
            lineStyle: { color: textColor },
            length: 5,
            inside: false
          },
          axisLabel: {
            showMaxLabel: false, // Hide maximum value label on y-axis
            color: textColor,
            fontSize: 16,
            fontFamily: 'Merriweather Sans',
            formatter: (value: any) => formatNumberToPrecision(value, 3)
          },
          splitLine: { show: false }
        },
        tooltip: {
          show: true,
          trigger: 'axis',
          backgroundColor: resolvedTheme === 'dark' ? '#08306b' : '#ffffff',
          borderColor: resolvedTheme === 'dark' ? '#55aaff' : '#333333',
          borderWidth: 1,
          textStyle: {
            color: textColor,
            fontSize: 12,
            fontFamily: 'Merriweather Sans'
          },
          axisPointer: {
            type: 'cross',
            label: {
              show: true,
              backgroundColor: resolvedTheme === 'dark' ? '#08306b' : '#ffffff',
              color: textColor,
              borderColor: resolvedTheme === 'dark' ? '#55aaff' : '#333333',
              borderWidth: 1,
              fontFamily: 'Merriweather Sans',
              formatter: function (params: any) {
                if (params.axisDimension === 'x') {
                  const xLabel = independentVar === 'pressure' ? 'P' : 'C';
                  return `${xLabel}: ${params.value.toFixed(3)}`;
                } else {
                  const yLabel = selectedIsotherm === 'langmuir' ? 'θ' : 'q';
                  return `${yLabel}: ${params.value.toFixed(3)}`;
                }
              }
            }
          },
          formatter: function (params: any) {
            if (Array.isArray(params) && params.length > 0) {
              const point = params[0];
              const xLabel = independentVar === 'pressure' ? 'P' : 'C';
              const yLabel = selectedIsotherm === 'langmuir' ? 'θ' : 'q';
              return `${xLabel}: ${formatNumberToPrecision(point.axisValue, 3)}<br/><span style="color: ${point.color};">${yLabel}: ${formatNumberToPrecision(point.value[1], 3)}</span>`;
            }
            return '';
          }
        },
        series: [
          {
            name: selectedIsotherm.charAt(0).toUpperCase() + selectedIsotherm.slice(1) + ' Isotherm',
            type: 'line',
            smooth: true,
            data: isothermData.c.map((val, index) => [val, isothermData.q[index]]),
            showSymbol: false,
            color: 'green',
            lineStyle: { width: 3.5 },
            animation: false
          }
        ],
        toolbox: {
          show: true,
          orient: 'vertical',
          right: 0,
          top: 'bottom',
          feature: {
            saveAsImage: {
              show: true,
              title: 'Save as Image',
              name: 'langmuir-isotherm',
              backgroundColor: resolvedTheme === 'dark' ? '#08306b' : '#ffffff',
              pixelRatio: 2
            }
          },
          iconStyle: {
            borderColor: textColor
          }
        }
      };
      setEchartsOptions(newOptions);
    }
  }, [isothermData, selectedIsotherm, langmuirK, independentVar, freundlichKf, freundlichNf, temkinAt, temkinBt, resolvedTheme]);

  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8 px-8 md:px-32">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Controls Column */}
          <div className="lg:col-span-1 space-y-6">
            {/* Combined Slider Card with Isotherm Selector */}
            <Card>
              <CardContent className="space-y-8 pt-6 pb-6">
                {/* Isotherm Type Selector */}
                <div className="flex items-center gap-2">
                  <Label htmlFor="isothermType" className="text-sm font-medium whitespace-nowrap">Isotherm Model:</Label>
                  <Select
                    value={selectedIsotherm}
                    onValueChange={(value) => setSelectedIsotherm(value as IsothermType)}
                  >
                    <SelectTrigger id="isothermType" className="flex-1">
                      <SelectValue placeholder="Select isotherm model" />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="langmuir">Langmuir</SelectItem>
                      <SelectItem value="freundlich">Freundlich</SelectItem>
                      <SelectItem value="temkin">Temkin</SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                {selectedIsotherm === 'langmuir' && (
                  <>
                    {/* Langmuir K Slider */}
                    <div className="space-y-3">
                      <Label htmlFor="langmuirK">
                        Equilibrium Constant (K): {langmuirK.toFixed(3)}
                      </Label>
                      <Slider
                        id="langmuirK"
                        min={0.01}
                        max={5}
                        step={0.01}
                        value={[langmuirK]}
                        onValueChange={(value) => setLangmuirK(value[0])}
                        style={{ '--primary': 'hsl(262 84% 58%)' } as React.CSSProperties}
                      />
                    </div>

                  </>
                )}

                {selectedIsotherm === 'freundlich' && (
                  <>
                    <div className="space-y-3">
                      <Label htmlFor="freundlichKf">
                        <span>Freundlich Constant (K<sub style={{fontSize: '0.75em', lineHeight: '1'}}>f</sub>): {freundlichKf.toFixed(2)}</span>
                      </Label>
                      <Slider
                        id="freundlichKf"
                        min={0.1}
                        max={10}
                        step={0.1}
                        value={[freundlichKf]}
                        onValueChange={(value) => setFreundlichKf(value[0])}
                        style={{ '--primary': 'hsl(38 92% 50%)' } as React.CSSProperties}
                      />
                    </div>
                    <div className="space-y-3">
                      <Label htmlFor="freundlichNf">
                        <span>Freundlich Heterogeneity (n<sub style={{fontSize: '0.75em', lineHeight: '1'}}>f</sub>): {freundlichNf.toFixed(2)}</span>
                      </Label>
                      <Slider
                        id="freundlichNf"
                        min={0.5}
                        max={5}
                        step={0.1}
                        value={[freundlichNf]}
                        onValueChange={(value) => setFreundlichNf(value[0])}
                        style={{ '--primary': 'hsl(0 84% 60%)' } as React.CSSProperties}
                      />
                    </div>
                  </>
                )}

                {selectedIsotherm === 'temkin' && (
                  <>
                    <div className="space-y-3">
                      <Label htmlFor="temkinAt">
                        <span>Temkin Constant (A<sub style={{fontSize: '0.75em', lineHeight: '1'}}>T</sub>): {temkinAt.toFixed(2)}</span>
                      </Label>
                      <Slider
                        id="temkinAt"
                        min={0.1}
                        max={10}
                        step={0.1}
                        value={[temkinAt]}
                        onValueChange={(value) => setTemkinAt(value[0])}
                        style={{ '--primary': 'hsl(38 92% 50%)' } as React.CSSProperties}
                      />
                    </div>
                    <div className="space-y-3">
                      <Label htmlFor="temkinBt">
                        <span>Temkin Constant (B<sub style={{fontSize: '0.75em', lineHeight: '1'}}>T</sub>): {temkinBt.toFixed(2)}</span>
                      </Label>
                      <Slider
                        id="temkinBt"
                        min={0.1}
                        max={10}
                        step={0.1}
                        value={[temkinBt]}
                        onValueChange={(value) => setTemkinBt(value[0])}
                        style={{ '--primary': 'hsl(0 84% 60%)' } as React.CSSProperties}
                      />
                    </div>
                  </>
                )}
              </CardContent>
            </Card>
            {error && (
              <Card className="bg-red-100 border-red-500">
                <CardContent>
                  <p className="text-red-600">{error}</p>
                </CardContent>
              </Card>
            )}
          </div>

          {/* Column 2: Plot */}
          <div className="lg:col-span-2 space-y-6">
            <Card>
              <CardContent className="py-2">
                <div className="relative aspect-square rounded-md">
                  {/* Concentration/Pressure Switch - Top Right */}
                  <div className="absolute top-4 right-4 z-10">
                    <div className="flex items-center gap-2 bg-muted rounded-lg p-1">
                      <Button
                        variant={independentVar === 'concentration' ? 'default' : 'ghost'}
                        size="sm"
                        onClick={() => setIndependentVar('concentration')}
                        className="h-8 px-3"
                      >
                        Concentration
                      </Button>
                      <Button
                        variant={independentVar === 'pressure' ? 'default' : 'ghost'}
                        size="sm"
                        onClick={() => setIndependentVar('pressure')}
                        className="h-8 px-3"
                      >
                        Pressure
                      </Button>
                    </div>
                  </div>
                  {loading && (
                    <div className="absolute inset-0 flex items-center justify-center text-white">
                      <div className="text-center">
                        <div className="mb-2">Loading & Calculating Isotherm Data...</div>
                      </div>
                    </div>
                  )}
                  {!loading && !isothermData && !error && (
                    <div className="absolute inset-0 flex items-center justify-center text-white">
                      Please provide inputs and update graph.
                    </div>
                  )}
                  {error && !loading && (
                    <div className="absolute inset-0 flex items-center justify-center text-red-400">
                      Error: {error}
                    </div>
                  )}
                  {!loading && isothermData && Object.keys(echartsOptions).length > 0 && (
                    <ReactECharts
                      ref={echartsRef}
                      echarts={echarts}
                      option={echartsOptions}
                      style={{ height: '100%', width: '100%', borderRadius: '0.375rem', overflow: 'hidden' }}
                      notMerge={true}
                      lazyUpdate={true}
                    />
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
