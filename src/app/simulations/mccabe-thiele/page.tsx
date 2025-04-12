'use client';

import React, { useState, useEffect, useCallback, useRef } from 'react';
// Import ECharts components
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
// Import specific types from the main echarts package
import type { EChartsOption, SeriesOption } from 'echarts';
import { LineChart, ScatterChart } from 'echarts/charts';
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  MarkLineComponent,
  MarkPointComponent,
  ToolboxComponent // Import ToolboxComponent for saveAsImage feature
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

// Register necessary ECharts components
echarts.use([
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  MarkLineComponent,
  MarkPointComponent,
  ToolboxComponent, // Register ToolboxComponent
  LineChart,
  ScatterChart,
  CanvasRenderer
]);

import { Card, CardContent, CardHeader, CardTitle, CardDescription, CardFooter } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Slider } from "@/components/ui/slider";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Tooltip as ShadTooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { Download } from 'lucide-react'; // Import Download icon


// Define ECharts point type (can be number array or object)
type EChartsPoint = [number, number] | [number | null, number | null];

// Helper function for delay
const delay = (ms: number) => new Promise(resolve => setTimeout(resolve, ms));

export default function McCabeThielePage() {
  // Input States
  const [comp1, setComp1] = useState('methanol');
  const [comp2, setComp2] = useState('water');
  const [temperatureC, setTemperatureC] = useState(27);
  const [pressureBar, setPressureBar] = useState(1);
  const [useTemperature, setUseTemperature] = useState(true);

  // Parameter States
  const [xd, setXd] = useState(0.9);
  const [xb, setXb] = useState(0.1);
  const [xf, setXf] = useState(0.5);
  const [q, setQ] = useState(1.0);
  const [r, setR] = useState(1.5);

  const buffer = 0.01; // Define buffer for constraints

  // Data & Control States
  const [equilibriumData, setEquilibriumData] = useState<any>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [volatilityInfo, setVolatilityInfo] = useState<string | null>(null); // Keep for axis label logic

  // Result States
  const [stages, setStages] = useState<number | null>(null);
  const [feedStage, setFeedStage] = useState<number | null>(null);

  // State for ECharts options - Use the imported EChartsOption type
  const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
  const echartsRef = useRef<ReactECharts | null>(null); // Ref for ECharts instance
  // State for displayed parameters in the title
  const [displayedComp1, setDisplayedComp1] = useState('');
  const [displayedComp2, setDisplayedComp2] = useState('');
  const [displayedTemp, setDisplayedTemp] = useState<number | null>(null);
  const [displayedPressure, setDisplayedPressure] = useState<number | null>(null);
  const [displayedUseTemp, setDisplayedUseTemp] = useState(true);


  useEffect(() => {
    // Initial fetch without arguments
    fetchVLEData();
    setDisplayedComp1('methanol');
    setDisplayedComp2('water');
    setDisplayedTemp(27);
    setDisplayedPressure(null);
    setDisplayedUseTemp(true);
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []); // fetchVLEData should be stable due to useCallback

  // ... (keep navbar reset useEffect) ...

  const fetchVLEData = useCallback(async (retryCount = 0) => {
    const maxRetries = 3;
    const retryDelay = 500; // milliseconds

    console.log(`DEBUG: fetchVLEData called (Attempt ${retryCount + 1}) with components:`, comp1, comp2);
    console.log(`DEBUG: Using ${useTemperature ? 'temperature' : 'pressure'} mode with value ${useTemperature ? temperatureC : pressureBar}`);

    // Set loading only on the first attempt
    if (retryCount === 0) {
        setLoading(true);
        setError(null);
    }

    try {
      console.log("DEBUG: Preparing API call to McCabe-Thiele endpoint");

      const response = await fetch('/api/mccabe-thiele', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          comp1: comp1,
          comp2: comp2,
          temperature: useTemperature ? temperatureC : null,
          pressure: !useTemperature ? pressureBar : null,
        }),
      });
      console.log(`DEBUG: API response status: ${response.status}`);

      const data = await response.json();
      console.log(`DEBUG: API response received. Success: ${!data.error}`);

      if (data.error) {
        console.error("DEBUG ERROR: API returned error:", data.error);
        // Check for specific retryable error
        if (data.error.includes("'PiecewiseHeatCapacity' object has no attribute 'append'") && retryCount < maxRetries) {
            console.log(`DEBUG: Retrying API call due to specific error (Attempt ${retryCount + 2})...`);
            await delay(retryDelay);
            fetchVLEData(retryCount + 1); // Recursive call for retry
            return; // Exit current attempt
        }
        throw new Error(data.error); // Throw non-retryable or max-retry errors
      }

      // --- Process successful data ---
      console.log(`DEBUG: VLE data received. Points: ${data.x_values?.length || 0}`);
      console.log(`DEBUG: Temperature: ${data.temperature}K, Pressure: ${data.pressure ? data.pressure/1e5 + 'bar' : 'Not specified'}`);
      console.log("DEBUG: Volatility info:", data.volatility);

      setVolatilityInfo(data.volatility?.message || `${comp1} / ${comp2} VLE data`);

      const xValues = data.x_values || [];
      const yValues = data.y_values || [];

      console.log(`DEBUG: First few x values: ${xValues.slice(0, 5).join(', ')}`);
      console.log(`DEBUG: First few y values: ${yValues.slice(0, 5).join(', ')}`);

      setEquilibriumData({
        x: xValues,
        y: yValues,
        polyCoeffs: data.poly_coeffs || []
      });

      // Update displayed parameters *after* successful fetch
      setDisplayedComp1(comp1);
      setDisplayedComp2(comp2);
      setDisplayedUseTemp(useTemperature);
      if (useTemperature) {
        setDisplayedTemp(temperatureC);
        setDisplayedPressure(null);
      } else {
        setDisplayedTemp(null);
        setDisplayedPressure(pressureBar);
      }
      // --- End process successful data ---

    } catch (err) {
      // Handle errors that occurred during fetch or after retries failed
      console.error("DEBUG ERROR: Error fetching VLE data:", err);
      setError(`Failed to load equilibrium data. ${err instanceof Error ? err.message : 'Please try again.'}`);
      // Clear potentially stale data on final failure
      setEquilibriumData(null);
      setVolatilityInfo(null);
      setDisplayedComp1('');
      setDisplayedComp2('');
      setDisplayedTemp(null);
      setDisplayedPressure(null);

    } finally {
      // Set loading false only after the final attempt (success or failure)
      if (retryCount === 0 || error || equilibriumData) { // Ensure loading is set false on final state
          setLoading(false);
      }
      console.log("DEBUG: fetchVLEData attempt completed");
    }
  }, [comp1, comp2, useTemperature, temperatureC, pressureBar, error, equilibriumData]); // Dependencies remain the same

  // *** generateEChartsOptions function ***
  const generateEChartsOptions = useCallback((xValues: number[], yValues: number[]) => {
    if (!xValues || !yValues || xValues.length === 0) {
        setEchartsOptions({}); // Clear chart if no data
        return;
    }

    // Use the imported SeriesOption type
    const series: SeriesOption[] = [];

    // 1. Equilibrium Line
    series.push({
      name: 'Equilibrium Line',
      type: 'line',
      data: xValues.map((x, i) => [x, yValues[i]]),
      color: 'yellow',
      symbol: 'none', // No markers
      lineStyle: { width: 2.5 }, // Increased width
      z: 5, // Ensure it's drawn on top of stages
      animation: false,
    });

    // 2. y=x Line
    series.push({
      name: 'y = x Line',
      type: 'line',
      data: [[0, 0], [1, 1]],
      color: 'white',
      symbol: 'none',
      lineStyle: { width: 1.5, type: 'dotted' }, // Increased width
      animation: false,
    });

    // --- Operating Lines Calculations ---
    const rectifyingSlope = r / (r + 1);
    const rectifyingIntercept = xd / (r + 1);
    let feedSlope: number;
    if (Math.abs(q - 1) < 1e-6) { feedSlope = Infinity; } else { feedSlope = q / (q - 1); }
    const feedIntercept = (feedSlope === Infinity) ? NaN : xf - feedSlope * xf;
    let xIntersect: number;
    let yIntersect: number;
    if (feedSlope === Infinity) { xIntersect = xf; yIntersect = rectifyingSlope * xf + rectifyingIntercept; }
    else { xIntersect = (feedIntercept - rectifyingIntercept) / (rectifyingSlope - feedSlope); yIntersect = rectifyingSlope * xIntersect + rectifyingIntercept; }
    xIntersect = Math.max(0, Math.min(1, xIntersect));
    yIntersect = Math.max(0, Math.min(1, yIntersect));
    const strippingSlope = (xIntersect === xb) ? Infinity : (yIntersect - xb) / (xIntersect - xb);
    const strippingIntercept = (strippingSlope === Infinity) ? NaN : xb - strippingSlope * xb;


    // 3. Rectifying Line
    series.push({
        name: 'Rectifying Section',
        type: 'line',
        data: [[xd, xd], [xIntersect, yIntersect]],
        color: 'orange',
        symbol: 'none',
        lineStyle: { width: 2.5 }, // Increased width
        animation: false,
    });

    // 4. Feed Line
    const feedLineData: EChartsPoint[] = [[xf, xf]];
    if (feedSlope === Infinity) { feedLineData.push([xf, yIntersect]); }
    else { feedLineData.push([xIntersect, yIntersect]); }
    series.push({
        name: 'Feed Section',
        type: 'line',
        data: feedLineData,
        color: 'red',
        symbol: 'none',
        lineStyle: { width: 2.5 }, // Increased width
        animation: false,
    });

    // 5. Stripping Line
    series.push({
        name: 'Stripping Section',
        type: 'line',
        data: [[xIntersect, yIntersect], [xb, xb]],
        color: 'green',
        symbol: 'none',
        lineStyle: { width: 2.5 }, // Increased width
        animation: false,
    });

    // 6. Key Points (as a separate scatter series)
    series.push({
        name: 'Key Points',
        type: 'scatter',
        data: [
            { value: [xd, xd], itemStyle: { color: 'orange' }, name: 'Distillate' },
            { value: [xb, xb], itemStyle: { color: 'green' }, name: 'Bottoms' },
            { value: [xf, xf], itemStyle: { color: 'red' }, name: 'Feed' }
        ],
        symbolSize: 8,
        label: { show: false }, // Don't show labels on points
        tooltip: { // Custom tooltip for points if needed, otherwise disable globally
            formatter: (params: any) => `${params.name}: (${params.value[0].toFixed(3)}, ${params.value[1].toFixed(3)})`
        },
        animation: false,
    });

    // --- Stage Calculation ---
    let stageCount = 0;
    let feedStageCount = 0;
    let currentX = xd;
    let currentY = xd;
    const stageLineData: EChartsPoint[] = []; // Single array for all stage segments
    let previousSectionIsRectifying = true;

    while (currentX > xb + 0.005 && stageCount < 25) {
        // ... (find intersection logic - intersectX, intersectY) ...
        let intersectX = NaN;
        for (let i = 0; i < xValues.length - 1; i++) {
            if ((yValues[i] <= currentY && yValues[i+1] >= currentY) || (yValues[i] >= currentY && yValues[i+1] <= currentY)) {
                 if (Math.abs(yValues[i+1] - yValues[i]) < 1e-9) { intersectX = (yValues[i] === currentY) ? xValues[i] : xValues[i+1]; }
                 else { const fraction = (currentY - yValues[i]) / (yValues[i+1] - yValues[i]); intersectX = xValues[i] + fraction * (xValues[i+1] - xValues[i]); }
                 break;
            }
        }
        if (isNaN(intersectX)) {
            if (currentY <= yValues[0]) intersectX = xValues[0];
            else if (currentY >= yValues[yValues.length - 1]) intersectX = xValues[xValues.length - 1];
            else break;
        }
        intersectX = Math.max(0, Math.min(1, intersectX));

        // Add horizontal segment points
        stageLineData.push([currentX, currentY]);
        stageLineData.push([intersectX, currentY]);
        stageLineData.push([null, null]); // Use null point to break line

        // ... (determine nextY, feedStageCount logic) ...
        let nextY: number;
        const currentSectionIsRectifying = intersectX > xIntersect;
        if (currentSectionIsRectifying) { nextY = rectifyingSlope * intersectX + rectifyingIntercept; }
        else { nextY = (strippingSlope === Infinity) ? xb : strippingSlope * intersectX + strippingIntercept; }
        nextY = Math.max(0, Math.min(1, nextY));
        if (feedStageCount === 0 && currentSectionIsRectifying !== previousSectionIsRectifying) { feedStageCount = stageCount + 1; }
        previousSectionIsRectifying = currentSectionIsRectifying;

        const isLastStage = (intersectX <= xb + 0.005) || (stageCount >= 24);
        if (isLastStage) {
            // Add vertical segment stopping at y=x
            stageLineData.push([intersectX, currentY]);
            stageLineData.push([intersectX, intersectX]);
            stageLineData.push([null, null]); // Break line
            currentX = intersectX; // Exit loop condition
        } else {
            // Add normal vertical segment
            stageLineData.push([intersectX, currentY]);
            stageLineData.push([intersectX, nextY]);
            stageLineData.push([null, null]); // Break line
            currentX = intersectX;
            currentY = nextY;
        }
        stageCount++;
    }
    if (feedStageCount === 0 && stageCount > 0) { feedStageCount = stageCount; }

    // Add the consolidated stage lines series
    series.push({
        name: 'Stages',
        type: 'line',
        data: stageLineData,
        color: 'white',
        symbol: 'none',
        lineStyle: { width: 2 }, // Increased width
        connectNulls: false, // Important: Do not connect across null points
        legendHoverLink: false, // Disable legend interaction if desired
        animation: false,
    });

    setStages(stageCount);
    setFeedStage(feedStageCount);

    // --- Define ECharts Options ---
    const capitalizeFirst = (str: string) => str ? str.charAt(0).toUpperCase() + str.slice(1) : '';
    const dispComp1Cap = capitalizeFirst(displayedComp1);
    const dispComp2Cap = capitalizeFirst(displayedComp2);
    const titleCondition = displayedUseTemp ? `${displayedTemp} °C` : `${displayedPressure} bar`;
    const titleText = displayedComp1 && displayedComp2 ?
        `McCabe-Thiele Diagram: ${dispComp1Cap} & ${dispComp2Cap} at ${titleCondition}` :
        'McCabe-Thiele Diagram';
    const moreVolatile = volatilityInfo?.includes(`${displayedComp1} is more volatile`) ? dispComp1Cap : dispComp2Cap;

    setEchartsOptions({
        backgroundColor: '#08306b',
        title: {
            text: titleText,
            left: 'center',
            textStyle: { color: 'white', fontSize: 18, fontFamily: 'Merriweather Sans' }
        },
        grid: {
            left: '5%',   // Increased left margin to match bottom
            right: '5%',  
            bottom: '5%', // Keep bottom margin
            top: '5%',
            containLabel: true
        },
        xAxis: {
            type: 'value',
            min: 0,
            max: 1,
            interval: 0.1,
            name: `Liquid Mole Fraction ${moreVolatile}`,
            nameLocation: 'middle',
            nameGap: 30, // Slightly increased gap due to larger font
            nameTextStyle: { color: 'white', fontSize: 15, fontFamily: 'Merriweather Sans' }, // Increased font size
            axisLine: { lineStyle: { color: 'white' } },
            axisTick: { lineStyle: { color: 'white' }, length: 5, inside: false },
            axisLabel: { color: 'white', fontSize: 16, fontFamily: 'Merriweather Sans', formatter: '{value}' }, // Increased font size
            splitLine: { show: false },
        },
        yAxis: {
            type: 'value',
            min: 0,
            max: 1,
            interval: 0.1,
            name: `Vapor Mole Fraction ${moreVolatile}`,
            nameLocation: 'middle',
            nameGap: 40, // Slightly increased gap due to larger font
            nameTextStyle: { color: 'white', fontSize: 15, fontFamily: 'Merriweather Sans' }, // Increased font size
            axisLine: { lineStyle: { color: 'white' } },
            axisTick: { lineStyle: { color: 'white' }, length: 5, inside: false },
            axisLabel: { color: 'white', fontSize: 16, fontFamily: 'Merriweather Sans', formatter: '{value}' }, // Increased font size
            splitLine: { show: false },
        },
        legend: {
            orient: 'vertical',
            right: '2%',
            top: 'center',
            // Add type assertion to ensure the array contains only strings
            data: series.map(s => s.name).filter(name => name !== 'Key Points' && name !== 'Stages') as string[],
            textStyle: { color: 'white', fontSize: 12, fontFamily: 'Merriweather Sans' },
            itemWidth: 12,
            itemHeight: 12,
            icon: 'rect',
        },
        tooltip: {
            show: true, // Disable tooltips globally for now
            // trigger: 'axis', // Or 'item'
            // axisPointer: { type: 'cross' }
        },
        animationDuration: 300,
        animationEasing: 'cubicInOut',
        // Add Toolbox for saving image
        toolbox: {
            show: true,
            orient: 'vertical', // Position vertically
            right: 0, // Position from the right
            top: 'bottom', // Center vertically
            feature: {
                saveAsImage: {
                    show: true,
                    title: 'Save as Image',
                    name: `mccabe-thiele-${displayedComp1}-${displayedComp2}`, // Dynamic filename
                    backgroundColor: '#08306b', // Match chart background
                    pixelRatio: 2 // Increase resolution
                }
            },
            iconStyle: {
                borderColor: '#fff' // White icon border
            }
        },
        series: series,
    });

  }, [xd, xb, xf, q, r, equilibriumData, volatilityInfo, displayedComp1, displayedComp2, displayedTemp, displayedPressure, displayedUseTemp, buffer]);


  // --- Helper Functions ---

  const getFeedQualityState = () => {
    if (q === 1) return "Saturated Liquid";
    if (q === 0) return "Saturated Vapor";
    if (q > 1) return "Subcooled Liquid";
    if (q > 0 && q < 1) return "Partially Vaporized";
    return "Super Heated Vapor";
  };

  const handleKeyDown = (event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key === 'Enter') {
      event.preventDefault(); // Prevent default form submission behavior if any
      fetchVLEData();
    }
  };

  const updateCompositions = (type: 'xd' | 'xf' | 'xb', value: number) => {
    // Get the current state values to check against
    const currentXd = xd;
    const currentXf = xf;
    const currentXb = xb;

    let targetValue = value; // Start with the value from the slider event

    // Check constraints and adjust the targetValue if needed (snap to boundary)
    if (type === 'xd') {
        // Prevent xd from going below or equal to xf + buffer
        if (targetValue <= currentXf + buffer) {
            targetValue = currentXf + buffer; // Snap to the boundary
        }
        // Ensure targetValue doesn't exceed absolute max
        targetValue = Math.min(targetValue, 0.99);
        setXd(targetValue);
    } else if (type === 'xf') {
        // Prevent xf from going above or equal to xd - buffer
        if (targetValue >= currentXd - buffer) {
            targetValue = currentXd - buffer; // Snap to the boundary
        }
        // Prevent xf from going below or equal to xb + buffer
        if (targetValue <= currentXb + buffer) {
            targetValue = currentXb + buffer; // Snap to the boundary
        }
         // Ensure targetValue is within absolute bounds
        targetValue = Math.max(0.01, Math.min(targetValue, 0.99));
        setXf(targetValue);
    } else if (type === 'xb') {
        // Prevent xb from going above or equal to xf - buffer
        if (targetValue >= currentXf - buffer) {
            targetValue = currentXf - buffer; // Snap to the boundary
        }
        // Ensure targetValue doesn't go below absolute min
        targetValue = Math.max(targetValue, 0.01);
        setXb(targetValue);
    }

    // Note: We directly set the state for the changed slider.
    // The constraints prevent invalid values based on the *current* state of the other sliders.
    // This avoids complex cascading updates and relies on the user not being able to move sliders past the calculated boundaries.
  };

  // --- End Helper Functions ---

  // Update useEffect to call the new options generation function
  useEffect(() => {
    if (equilibriumData?.x && equilibriumData?.y) {
      generateEChartsOptions(equilibriumData.x, equilibriumData.y);
    } else {
        setEchartsOptions({});
        setStages(null);
        setFeedStage(null);
    }
  }, [equilibriumData, generateEChartsOptions]);


  // Wrapper function for the button click
  const handleUpdateGraphClick = () => {
    fetchVLEData(); // Call without arguments, retryCount defaults to 0
  };

  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8 px-8 md:px-32">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">

          {/* Column 1: Input and Slider Cards */}
          <div className="lg:col-span-1 space-y-6">
            {/* Input Card */}
            <Card>
              <CardContent className="space-y-6 p-4">
                {/* Combined Input Group */}
                <div className="space-y-4">
                  {/* Temperature/Pressure Toggle and Input FIRST */}
                  <Tabs value={useTemperature ? "temperature" : "pressure"} onValueChange={(value) => setUseTemperature(value === "temperature")}>
                      <TabsList className="grid w-full grid-cols-2 mb-4">
                          <TabsTrigger value="temperature">Temperature</TabsTrigger>
                          <TabsTrigger value="pressure">Pressure</TabsTrigger>
                      </TabsList>
                      <TabsContent value="temperature" className="space-y-2">
                           <div className="flex items-center gap-3">
                              <Label htmlFor="temperature" className="w-30 whitespace-nowrap">Temperature (°C):</Label>
                              <Input
                                 id="temperature" type="number" value={temperatureC}
                                 onChange={(e) => setTemperatureC(Number(e.target.value))}
                                 onKeyDown={handleKeyDown} // Usage is now after definition
                                 min="0" step="1"
                                 required={useTemperature} disabled={!useTemperature} className="flex-1"
                              />
                           </div>
                      </TabsContent>
                       <TabsContent value="pressure" className="space-y-2">
                           <div className="flex items-center gap-3">
                              <Label htmlFor="pressure" className="w-24 whitespace-nowrap">Pressure (bar):</Label>
                              <Input
                                 id="pressure" type="number" value={pressureBar}
                                 onChange={(e) => setPressureBar(Number(e.target.value))}
                                 onKeyDown={handleKeyDown} // Usage is now after definition
                                 min="0.1" step="0.1"
                                 required={!useTemperature} disabled={useTemperature} className="flex-1"
                              />
                           </div>
                       </TabsContent>
                   </Tabs>

                  {/* Component Inputs - Inline */}
                  <div className="space-y-3">
                     <div className="flex items-center gap-3">
                        <Label htmlFor="comp1" className="w-24 whitespace-nowrap">Component 1:</Label>
                        <Input
                          id="comp1" value={comp1} onChange={(e) => setComp1(e.target.value)}
                          onKeyDown={handleKeyDown} // Usage is now after definition
                          required className="flex-1"
                        />
                     </div>
                     <div className="flex items-center gap-3">
                        <Label htmlFor="comp2" className="w-24 whitespace-nowrap">Component 2:</Label>
                        <Input
                          id="comp2" value={comp2} onChange={(e) => setComp2(e.target.value)}
                          onKeyDown={handleKeyDown} // Usage is now after definition
                          required className="flex-1"
                        />
                     </div>
                  </div>
                </div>

                {/* Submit Button */}
                {/* Use the wrapper function here */}
                <Button onClick={handleUpdateGraphClick} disabled={loading} className="w-full">
                  {loading ? 'Calculating...' : 'Update Graph'}
                </Button>

                {error && <p className="text-sm text-red-500 mt-2">{error}</p>}
              </CardContent>
            </Card>

            {/* Slider Card */}
            <Card>
               <CardContent className="space-y-8 pt-6 pb-6">
                 {/* xd Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="xd">
                     <span dangerouslySetInnerHTML={{ __html: `Distillate Comp (x<sub>D</sub>): ${xd.toFixed(2)}` }} />
                   </Label>
                   <Slider
                     id="xd"
                     min={0.01}
                     max={0.99}
                     step={0.01}
                     value={[xd]}
                     onValueChange={(value) => updateCompositions('xd', value[0])} // Usage is now after definition
                     style={{ '--primary': 'hsl(38 92% 50%)' } as React.CSSProperties}
                   />
                 </div>
                 {/* xf Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="xf">
                     <span dangerouslySetInnerHTML={{ __html: `Feed Comp (x<sub>F</sub>): ${xf.toFixed(2)}` }} />
                   </Label>
                   <Slider
                     id="xf"
                     min={0.01}
                     max={0.99}
                     step={0.01}
                     value={[xf]}
                     onValueChange={(value) => updateCompositions('xf', value[0])} // Usage is now after definition
                     style={{ '--primary': 'hsl(0 84% 60%)' } as React.CSSProperties}
                   />
                 </div>
                 {/* xb Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="xb">
                     <span dangerouslySetInnerHTML={{ __html: `Bottoms Comp (x<sub>B</sub>): ${xb.toFixed(2)}` }} />
                   </Label>
                   <Slider
                     id="xb"
                     min={0.01}
                     max={0.99}
                     step={0.01}
                     value={[xb]}
                     onValueChange={(value) => updateCompositions('xb', value[0])} // Usage is now after definition
                     style={{ '--primary': 'hsl(142 71% 45%)' } as React.CSSProperties}
                   />
                 </div>
                 {/* q Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="q" className="flex items-center">
                     Feed Quality (q): {q.toFixed(2)}
                       <ShadTooltip>
                         <TooltipTrigger asChild>
                           <Button variant="ghost" size="icon" className="h-5 w-5 rounded-full">
                             <span className="text-xs">ⓘ</span>
                           </Button>
                         </TooltipTrigger>
                         <TooltipContent><p>{getFeedQualityState()}</p></TooltipContent> {/* Usage is now after definition */}
                       </ShadTooltip>
                   </Label>
                   <Slider
                     id="q" min={-1} max={2} step={0.05} value={[q]}
                     onValueChange={(value) => setQ(value[0])}
                     style={{ '--primary': 'hsl(0 0% 98%)' } as React.CSSProperties}
                   />
                 </div>
                 {/* r Slider */}
                 <div className="space-y-3">
                   <Label htmlFor="r">Reflux Ratio (R): {r.toFixed(2)}</Label>
                   <Slider
                     id="r" min={0.1} max={10} step={0.05} value={[r]}
                     onValueChange={(value) => setR(value[0])}
                     style={{ '--primary': 'hsl(262 84% 58%)' } as React.CSSProperties}
                   />
                 </div>
               </CardContent>
            </Card>
          </div>

          {/* Column 2: Plot and Results Cards */}
          <div className="lg:col-span-2 space-y-6">
            {/* Plot Card */}
            <Card>
              <CardContent className="py-2">
                {/* Chart Container */}
                <div className="relative h-[500px] md:h-[600px] rounded-md" style={{ backgroundColor: '#08306b' }}>
                   {loading && (
                    <div className="absolute inset-0 flex items-center justify-center text-white">
                      <div className="text-center">
                        <div className="mb-2">Loading data...</div>
                        <div className="text-sm text-gray-300">Calculating equilibrium curve</div>
                      </div>
                    </div>
                   )}
                   {!loading && !equilibriumData && !error && (
                    <div className="absolute inset-0 flex items-center justify-center text-white">Please provide inputs and update graph.</div>
                   )}
                   {error && !loading && (
                    <div className="absolute inset-0 flex items-center justify-center text-red-400">Error: {error}</div>
                   )}
                  {/* Render ReactECharts component */}
                  {!loading && equilibriumData && Object.keys(echartsOptions).length > 0 && (
                    <ReactECharts
                      ref={echartsRef} // Assign ref
                      echarts={echarts}
                      option={echartsOptions}
                      style={{
                        height: '100%',
                        width: '100%',
                        borderRadius: '0.375rem', // Use Tailwind's rounded-md value
                        overflow: 'hidden'
                      }}
                      notMerge={false}
                      lazyUpdate={false}
                    />
                  )}
                   {/* Custom Save Button (Optional - ECharts toolbox is usually preferred) */}
                   {/*
                   <Button
                       variant="outline"
                       size="icon"
                       onClick={saveChartAsImage}
                       className="absolute bottom-4 right-4 z-10 bg-background/70 hover:bg-background/90 backdrop-blur-sm"
                       title="Save Chart as PNG"
                   >
                       <Download className="h-4 w-4" />
                   </Button>
                   */}
                </div>

                {/* Results Display */}
                {stages !== null && feedStage !== null && (
                  <div className="grid grid-cols-2 gap-4 text-center mt-4 pt-2">
                    <div className="p-4 bg-muted rounded-md">
                      <p className="text-sm font-medium">Number of Stages</p>
                      <p className="text-2xl font-bold">{stages}</p>
                    </div>
                    <div className="p-4 bg-muted rounded-md">
                      <p className="text-sm font-medium">Feed Stage</p>
                      <p className="text-2xl font-bold">{feedStage}</p>
                    </div>
                  </div>
                )}
              </CardContent>
            </Card>

          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}
