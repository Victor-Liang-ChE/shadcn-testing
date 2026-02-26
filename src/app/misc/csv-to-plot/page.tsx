'use client';

import React, { useState, useCallback, useEffect, useRef } from 'react';
import Papa, { ParseResult } from 'papaparse';
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsOption, SeriesOption } from 'echarts';
import { LineChart, ScatterChart, BarChart } from 'echarts/charts'; // Import BarChart
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  ToolboxComponent,
  DataZoomComponent // Import DataZoomComponent for zooming/panning
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Slider } from "@/components/ui/slider";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Switch } from "@/components/ui/switch"; // Import Switch
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip"; // Import Tooltip components
import { UploadCloud, ArrowRightLeft } from 'lucide-react';
import { cn } from "@/lib/utils";
import { useTheme } from 'next-themes'; // Import useTheme

// Register necessary ECharts components
echarts.use([
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  ToolboxComponent,
  DataZoomComponent, // Register DataZoomComponent
  LineChart,
  ScatterChart,
  BarChart, // Register BarChart
  CanvasRenderer
]);

// Define type for parsed data row
type DataRow = { [key: string]: string | number | null };

// Helper function for linear regression
function calculateLinearRegression(data: [number, number][]): { slope: number; intercept: number; rSquared: number } | null {
  const n = data.length;
  if (n < 2) return null; // Need at least two points

  let sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
  for (const [x, y] of data) {
    sumX += x;
    sumY += y;
    sumXY += x * y;
    sumX2 += x * x;
    sumY2 += y * y;
  }

  const meanX = sumX / n;
  const meanY = sumY / n;

  const ssXX = sumX2 - n * meanX * meanX;
  const ssYY = sumY2 - n * meanY * meanY;
  const ssXY = sumXY - n * meanX * meanY;

  // Avoid division by zero if all X values are the same
  if (Math.abs(ssXX) < 1e-10) return null;

  const slope = ssXY / ssXX;
  const intercept = meanY - slope * meanX;

  // Calculate R-squared
  const predictedY = data.map(([x]) => slope * x + intercept);
  let ssRes = 0;
  for (let i = 0; i < n; i++) {
    ssRes += (data[i][1] - predictedY[i]) ** 2;
  }

  // Avoid division by zero if all Y values are the same
  const rSquared = Math.abs(ssYY) < 1e-10 ? 1 : 1 - (ssRes / ssYY);

  return { slope, intercept, rSquared };
}


export default function CsvToPlotPage() {
  // File & Data State
  const [csvData, setCsvData] = useState<DataRow[]>([]);
  const [headers, setHeaders] = useState<string[]>([]);
  const [fileName, setFileName] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [isDragging, setIsDragging] = useState(false);

  // Chart Configuration State
  const [xAxisColumn, setXAxisColumn] = useState<string | null>(null);
  const [yAxisColumn, setYAxisColumn] = useState<string | null>(null);
  const [chartType, setChartType] = useState<'scatter' | 'line' | 'bar'>('scatter'); // Add 'bar' type
  const [chartOptions, setChartOptions] = useState<EChartsOption>({});
  const [title, setTitle] = useState<string>('CSV Plot');
  const [xAxisLabel, setXAxisLabel] = useState<string>('X Axis');
  const [yAxisLabel, setYAxisLabel] = useState<string>('Y Axis');
  const [titleFontSize, setTitleFontSize] = useState<number>(16); // Default to 16px
  const [axisLabelFontSize, setAxisLabelFontSize] = useState<number>(14); // Default to 14px
  const [tickLabelFontSize, setTickLabelFontSize] = useState<number>(14); // Default to 14px
  const [showGrid, setShowGrid] = useState<boolean>(false); // Default grid lines OFF
  const [showBestFit, setShowBestFit] = useState<boolean>(false); // State for best fit line toggle

  // Regression State
  const [regressionSlope, setRegressionSlope] = useState<number | null>(null);
  const [regressionIntercept, setRegressionIntercept] = useState<number | null>(null);
  const [rSquared, setRSquared] = useState<number | null>(null);

  const echartsRef = useRef<ReactECharts | null>(null);
  const { resolvedTheme } = useTheme(); // Get the resolved theme ('light' or 'dark')

  // --- Drag and Drop Handlers ---
  const handleDragEnter = useCallback((e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(true);
  }, []);

  const handleDragLeave = useCallback((e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
    // Simplified logic: If leaving the specific drop zone, set dragging to false
    // This relies on dragEnter/dragOver on the same element to set it back to true if needed
    setIsDragging(false);
  }, []);

  const handleDragOver = useCallback((e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault(); // Necessary to allow drop
    e.stopPropagation();
    setIsDragging(true); // Keep true while dragging over the drop zone
  }, []);

  const handleDrop = useCallback((e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);
    setError(null); // Clear previous errors
    setCsvData([]); // Clear previous data
    setHeaders([]);
    setXAxisColumn(null);
    setYAxisColumn(null);
    setFileName(null);

    const files = e.dataTransfer.files;
    if (files && files.length > 0) {
      const file = files[0];
      setFileName(file.name);

      if (file.type === 'text/csv' || file.name.toLowerCase().endsWith('.csv')) {
        // Apply 'as any' to the config object
        Papa.parse<DataRow>(file, {
          header: true,
          skipEmptyLines: true,
          dynamicTyping: true, // Automatically convert numbers
          complete: (results: ParseResult<DataRow>) => { // Add type ParseResult<DataRow>
            if (results.errors.length > 0) {
              console.error("CSV Parsing Errors:", results.errors);
              setError(`Error parsing CSV: ${results.errors[0].message}`);
              return;
            }
            if (!results.data || results.data.length === 0) {
              setError("CSV file is empty or could not be parsed correctly.");
              return;
            }
            if (!results.meta.fields || results.meta.fields.length < 2) {
              setError("CSV must contain at least two columns with headers.");
              return;
            }

            console.log("Parsed Data:", results.data);
            console.log("Detected Headers:", results.meta.fields);

            setCsvData(results.data);
            setHeaders(results.meta.fields);

            // Attempt to set initial axes (find first two numeric columns)
            const numericHeaders = results.meta.fields.filter((header: string) => // Add type string
              results.data.some((row: DataRow) => typeof row[header] === 'number') // Add type DataRow
            );

            if (numericHeaders.length >= 2) {
              setXAxisColumn(numericHeaders[0]);
              setYAxisColumn(numericHeaders[1]);
              setXAxisLabel(numericHeaders[0]);
              setYAxisLabel(numericHeaders[1]);
              setTitle(`${numericHeaders[1]} vs ${numericHeaders[0]}`);
            } else if (results.meta.fields.length >= 2) {
              // Fallback to first two headers if not enough numeric found
              setXAxisColumn(results.meta.fields[0]);
              setYAxisColumn(results.meta.fields[1]);
              setXAxisLabel(results.meta.fields[0]);
              setYAxisLabel(results.meta.fields[1]);
              setTitle(`${results.meta.fields[1]} vs ${results.meta.fields[0]}`);
              setError("Warning: Could not automatically detect two numeric columns. Please select columns manually.");
            } else {
              setError("Could not determine initial X and Y axes. Please select columns.");
            }
          },
          error: (error: Error, _file?: File) => { // Keep Error type for basic info
            console.error("CSV Parsing Failed:", error);
            setError(`Failed to parse CSV file: ${error.message}`);
          }
        } as any); // Add 'as any' here
      } else {
        setError('Invalid file type. Please drop a .csv file.');
      }
    }
  }, []);

  // --- Chart Generation ---
  const generateChartOptions = useCallback(() => {
    if (!csvData.length || !xAxisColumn || !yAxisColumn) {
      setChartOptions({}); // Clear chart if no data or axes selected
      return;
    }

    // Filter data
    // For bar charts, we might want to handle categorical X-axis differently,
    // but ECharts often plots numeric X vs numeric Y acceptably for bars too.
    // Let's keep the filtering logic the same for now.
    const plotData = csvData
      .map(row => [row[xAxisColumn!], row[yAxisColumn!]])
      .filter(point => typeof point[0] === 'number' && typeof point[1] === 'number') as [number, number][];

    // ... (error check for plotData remains the same) ...
    if (plotData.length === 0) {
        setError(`No numeric data found for the selected columns: '${xAxisColumn}' and '${yAxisColumn}'.`);
        setChartOptions({});
        return;
    } else {
        setError(null); // Clear error if data is found
    }


    // Define colors based on theme
    const textColor = resolvedTheme === 'dark' ? '#E5E7EB' : '#1F2937';
    const mutedColor = resolvedTheme === 'dark' ? '#9CA3AF' : '#6B7280';
    // const borderColor = resolvedTheme === 'dark' ? 'hsl(var(--border))' : 'hsl(var(--border))'; // Use gridLineColor instead
    const chartFont = 'Merriweather Sans, sans-serif';
    const saveBackgroundColor = resolvedTheme === 'dark' ? 'hsl(215, 40%, 30%)' : '#FFFFFF';
    const dataColor = resolvedTheme === 'dark' ? '#FFFFFF' : '#000000'; // White for dark, Black for light
    const gridLineColor = resolvedTheme === 'dark' ? 'rgba(255, 255, 255, 0.2)' : 'rgba(0, 0, 0, 0.2)'; // Lighter/darker grid lines

    const series: SeriesOption[] = [{
      name: `${yAxisLabel} vs ${xAxisLabel}`,
      type: chartType,
      data: plotData,
      color: dataColor, // Use theme-based data color
      symbolSize: chartType === 'scatter' ? 6 : 4,
      smooth: chartType === 'line',
      showSymbol: chartType !== 'line', // Hide symbols only for line chart
      lineStyle: { width: chartType === 'line' ? 2 : 0 },
      itemStyle: {
          opacity: 0.8,
          // Ensure bar color is also white if needed, 'color' above should handle it
      },
      animation: false,
      // Bar chart specific styling (optional)
      barWidth: chartType === 'bar' ? '60%' : undefined, // Example: Adjust bar width
    }];

    // Add Best Fit Line Series if applicable
    if (chartType === 'line' && showBestFit && regressionSlope !== null && regressionIntercept !== null) {
      const xMin = Math.min(...plotData.map(p => p[0]));
      const xMax = Math.max(...plotData.map(p => p[0]));
      const yMin = regressionSlope * xMin + regressionIntercept;
      const yMax = regressionSlope * xMax + regressionIntercept;

      series.push({
        name: 'Best Fit Line',
        type: 'line',
        data: [[xMin, yMin], [xMax, yMax]],
        color: 'red', // Red color for best fit line
        symbol: 'none',
        lineStyle: {
          width: 2,
          type: 'dashed' // Dashed line style
        },
        tooltip: { // Optional: Custom tooltip for the line itself
          formatter: `y = ${regressionSlope.toFixed(3)}x + ${regressionIntercept.toFixed(3)}<br/>R² = ${rSquared?.toFixed(4)}`
        },
        animation: false,
        z: 10 // Ensure it's drawn on top
      });
    }

    // Determine axis type based on chart type (optional, ECharts is flexible)
    // For simplicity, keep both as 'value' for now. If strict categorical needed, adjust here.
    const xAxisType = 'value';
    const yAxisType = 'value';

    setChartOptions({
      backgroundColor: 'transparent',
      title: {
        text: title,
        left: 'center',
        // top: 10, // Example if explicit pixel positioning is needed
        textStyle: {
            color: textColor, // Use theme-based color
            fontSize: titleFontSize,
            fontFamily: chartFont // Use defined font
        }
      },
      grid: {
        left: '10%',
        right: '7%',
        bottom: '12%', // ADJUSTED: Reduced grid bottom margin to move axis up
        top: '10%',
        containLabel: true
      },
      tooltip: {
        trigger: chartType === 'bar' ? 'item' : 'axis', // Use 'item' trigger for better bar tooltips
        axisPointer: {
          type: chartType === 'bar' ? 'shadow' : 'cross' // Use 'shadow' for bar axis pointer
        },
        formatter: (params: any) => {
            if (!Array.isArray(params)) params = [params];
            const point = params[0];
            if (point && point.value) {
                 // Customize tooltip for bar if needed, otherwise default is fine
                 // Example: return `${point.seriesName}<br/>${xAxisLabel}: ${point.value[0]}<br/>${yAxisLabel}: ${point.value[1]}`;
                 return `<b>${point.seriesName}</b><br/>${xAxisLabel}: ${point.value[0]?.toFixed(3)}<br/>${yAxisLabel}: ${point.value[1]?.toFixed(3)}`;
            }
            return '';
        }
      },
      xAxis: {
        type: xAxisType,
        name: xAxisLabel,
        nameLocation: 'middle',
        nameGap: 35, // ADJUSTED: Increased gap to push label down relative to axis line
        nameTextStyle: { color: textColor, fontSize: axisLabelFontSize, fontFamily: chartFont }, // Use font
        axisLine: { lineStyle: { color: mutedColor } }, // Use theme-based muted color
        axisTick: { lineStyle: { color: mutedColor } }, // Use theme-based muted color
        axisLabel: { color: textColor, fontSize: tickLabelFontSize, fontFamily: chartFont }, // Use font
        splitLine: { show: showGrid, lineStyle: { type: 'dashed', color: gridLineColor } }, // Use theme-based grid color
        scale: true
      },
      yAxis: {
        type: yAxisType,
        name: yAxisLabel,
        nameLocation: 'middle',
        nameGap: 45,
        nameTextStyle: { color: textColor, fontSize: axisLabelFontSize, fontFamily: chartFont }, // Use font
        axisLine: { lineStyle: { color: mutedColor } }, // Use theme-based muted color
        axisTick: { lineStyle: { color: mutedColor } }, // Use theme-based muted color
        axisLabel: { color: textColor, fontSize: tickLabelFontSize, fontFamily: chartFont }, // Use font
        splitLine: { show: showGrid, lineStyle: { type: 'dashed', color: gridLineColor } }, // Use theme-based grid color
        scale: true
      },
      toolbox: {
        show: true,
        orient: 'horizontal', // Keep horizontal
        right: 20, // Keep right position
        top: 10, // ADJUSTED: Position from the top instead of bottom
        // bottom: 0, // REMOVED
        feature: {
          saveAsImage: {
            show: true,
            title: 'Save as Image',
            name: `plot-${fileName?.split('.')[0] ?? 'data'}`,
            backgroundColor: saveBackgroundColor, // Use theme-aware background color
            excludeComponents: ['dataZoom', 'toolbox'], // Exclude dataZoom AND toolbox
          },
          dataZoom: {
              show: true,
              title: {
                  zoom: 'Zoom',
                  back: 'Reset Zoom'
              }
          },
          restore: {
              show: true,
              title: 'Restore'
          },
        },
        iconStyle: {
            borderColor: textColor // Use theme-based color
        }
      },
      dataZoom: [
          { // X-axis Slider
              type: 'slider',
              xAxisIndex: 0,
              filterMode: 'filter',
              bottom: 10, // Keep slider positioned low
              height: 20,
              handleIcon: 'M10.7,11.9v-1.3H9.3v1.3c-4.9,0.3-8.8,4.4-8.8,9.4c0,5,3.9,9.1,8.8,9.4v1.3h1.3v-1.3c4.9-0.3,8.8-4.4,8.8-9.4C19.5,16.3,15.6,12.2,10.7,11.9z M13.3,24.4H6.7V23h6.6V24.4z M13.3,19.6H6.7v-1.4h6.6V19.6z',
              handleSize: '80%',
              handleStyle: {
                  color: '#fff',
                  shadowBlur: 3,
                  shadowColor: 'rgba(0, 0, 0, 0.6)',
                  shadowOffsetX: 2,
                  shadowOffsetY: 2
              },
              textStyle: {
                  color: textColor
              }
          },
          { // X-axis Inside (Scroll)
              type: 'inside',
              xAxisIndex: 0,
              filterMode: 'filter'
          },
          { // Y-axis Slider
              type: 'slider',
              yAxisIndex: 0,
              filterMode: 'filter',
              left: 0,
              width: 20,
              handleIcon: 'M10.7,11.9v-1.3H9.3v1.3c-4.9,0.3-8.8,4.4-8.8,9.4c0,5,3.9,9.1,8.8,9.4v1.3h1.3v-1.3c4.9-0.3,8.8-4.4,8.8-9.4C19.5,16.3,15.6,12.2,10.7,11.9z M13.3,24.4H6.7V23h6.6V24.4z M13.3,19.6H6.7v-1.4h6.6V19.6z',
              handleSize: '80%',
              handleStyle: {
                  color: '#fff',
                  shadowBlur: 3,
                  shadowColor: 'rgba(0, 0, 0, 0.6)',
                  shadowOffsetX: 2,
                  shadowOffsetY: 2
              },
               textStyle: {
                  color: textColor
              }
          },
          { // Y-axis Inside (Scroll)
              type: 'inside',
              yAxisIndex: 0,
              filterMode: 'filter'
          }
      ],
      series: series,
    });
  }, [csvData, xAxisColumn, yAxisColumn, chartType, title, xAxisLabel, yAxisLabel, titleFontSize, axisLabelFontSize, tickLabelFontSize, fileName, showGrid, resolvedTheme, showBestFit, regressionSlope, regressionIntercept, rSquared]);

  // --- Effects ---
  // Regenerate chart when data or config changes (now includes theme)
  useEffect(() => {
    generateChartOptions();
  }, [generateChartOptions]);

  // Calculate Regression when needed
  useEffect(() => {
    if (chartType === 'line' && showBestFit && csvData.length > 0 && xAxisColumn && yAxisColumn) {
      const plotData = csvData
        .map(row => [row[xAxisColumn], row[yAxisColumn]])
        .filter(point => typeof point[0] === 'number' && typeof point[1] === 'number') as [number, number][];

      if (plotData.length >= 2) {
        const regressionResult = calculateLinearRegression(plotData);
        if (regressionResult) {
          setRegressionSlope(regressionResult.slope);
          setRegressionIntercept(regressionResult.intercept);
          setRSquared(regressionResult.rSquared);
        } else {
          // Handle case where regression fails (e.g., vertical line)
          setRegressionSlope(null);
          setRegressionIntercept(null);
          setRSquared(null);
          // Optionally set an error message
        }
      } else {
        // Not enough data points
        setRegressionSlope(null);
        setRegressionIntercept(null);
        setRSquared(null);
      }
    } else {
      // Clear regression results if not applicable
      setRegressionSlope(null);
      setRegressionIntercept(null);
      setRSquared(null);
    }
  }, [chartType, showBestFit, csvData, xAxisColumn, yAxisColumn]); // Dependencies for calculation


  // --- UI Handlers ---
  const handleSwapAxes = () => {
    if (xAxisColumn && yAxisColumn) {
      const tempXCol = xAxisColumn;
      const tempXLabel = xAxisLabel;
      setXAxisColumn(yAxisColumn);
      setXAxisLabel(yAxisLabel);
      setYAxisColumn(tempXCol);
      setYAxisLabel(tempXLabel);
      // Update title automatically if it seems default
      if (title === `${yAxisLabel} vs ${xAxisLabel}`) {
          setTitle(`${tempXLabel} vs ${yAxisLabel}`);
      }
    }
  };

  const handleXAxisChange = (value: string) => {
    setXAxisColumn(value);
    setXAxisLabel(value); // Auto-update label
    if (yAxisColumn) setTitle(`${yAxisLabel} vs ${value}`);
  };

  const handleYAxisChange = (value: string) => {
    setYAxisColumn(value);
    setYAxisLabel(value); // Auto-update label
    if (xAxisColumn) setTitle(`${value} vs ${xAxisColumn}`);
  };


  // --- Render ---

  // Calculate plotData in the render scope for conditional checks
  const plotData = (csvData.length > 0 && xAxisColumn && yAxisColumn)
    ? csvData
        .map(row => [row[xAxisColumn], row[yAxisColumn]])
        .filter(point => typeof point[0] === 'number' && typeof point[1] === 'number') as [number, number][]
    : []; // Default to empty array if data/columns not ready

  return (
    <TooltipProvider>
      {/* Remove drag handlers from main element */}
      <main className="container mx-auto p-4 md:p-8 relative">

        {/* Main Content Grid */}
        {/* Grid logic: 1 column initially, 3 columns when data is loaded */}
        <div className={`grid grid-cols-1 ${csvData.length > 0 ? 'lg:grid-cols-3' : 'lg:grid-cols-1'} gap-6 transition-all duration-300 ease-in-out`}>

          {/* Controls Column - Render only when data is loaded */}
          {csvData.length > 0 && headers.length > 0 && (
            <div className="lg:col-span-1 space-y-6">
              <Card>
                <CardHeader>
                  <CardTitle>CSV Plotter</CardTitle>
                  <CardDescription>
                    {fileName ? `File: ${fileName}` : "Configure your plot."}
                  </CardDescription>
                </CardHeader>
                <CardContent className="space-y-6">
                  {/* Chart Type Selection - Moved UP */}
                   <div className="space-y-2">
                      <Label htmlFor="chart-type-select">Chart Type</Label>
                      <Select value={chartType} onValueChange={(v) => setChartType(v as 'scatter' | 'line' | 'bar')}>
                        <SelectTrigger id="chart-type-select">
                          <SelectValue placeholder="Select chart type" />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectItem value="scatter">Scatter Plot</SelectItem>
                          <SelectItem value="line">Line Plot</SelectItem>
                          <SelectItem value="bar">Bar Chart</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>

                  {/* Axis Selection - Adjusted Grid Layout */}
                  <div className="grid grid-cols-[1fr_auto_1fr] gap-2 items-end"> {/* 3 columns: Select, Button, Select */}
                    <div className="space-y-2">
                      <Label htmlFor="x-axis-select">X Axis</Label>
                      <Select value={xAxisColumn || ''} onValueChange={handleXAxisChange}>
                        {/* Added truncate classes */}
                        <SelectTrigger id="x-axis-select" className="w-full overflow-hidden">
                           <span className="truncate">
                             <SelectValue placeholder="Select X column" />
                           </span>
                        </SelectTrigger>
                        <SelectContent>
                          {headers.map(header => (
                            <SelectItem key={`x-${header}`} value={header}>{header}</SelectItem>
                          ))}
                        </SelectContent>
                      </Select>
                    </div>

                    {/* Swap Button - Icon only with Tooltip */}
                    <Tooltip>
                      <TooltipTrigger asChild>
                        <Button
                          variant="outline"
                          size="icon" // Use icon size
                          onClick={handleSwapAxes}
                          disabled={!xAxisColumn || !yAxisColumn}
                          className="mt-auto" // Align button vertically if needed
                        >
                          <ArrowRightLeft className="h-4 w-4" />
                          <span className="sr-only">Swap Axes</span> {/* Screen reader text */}
                        </Button>
                      </TooltipTrigger>
                      <TooltipContent>
                        <p>Swap Axes</p>
                      </TooltipContent>
                    </Tooltip>

                    <div className="space-y-2">
                      <Label htmlFor="y-axis-select">Y Axis</Label>
                      <Select value={yAxisColumn || ''} onValueChange={handleYAxisChange}>
                         {/* Added truncate classes */}
                        <SelectTrigger id="y-axis-select" className="w-full overflow-hidden">
                           <span className="truncate">
                             <SelectValue placeholder="Select Y column" />
                           </span>
                        </SelectTrigger>
                        <SelectContent>
                          {headers.map(header => (
                            <SelectItem key={`y-${header}`} value={header}>{header}</SelectItem>
                          ))}
                        </SelectContent>
                      </Select>
                    </div>
                  </div>
                  {/* Removed original Swap Button location */}

                  {/* Label Inputs */}
                  <div className="space-y-3">
                    <div className="space-y-2">
                      <Label htmlFor="title-input">Chart Title</Label>
                      <Input id="title-input" value={title} onChange={(e) => setTitle(e.target.value)} placeholder="Enter chart title" />
                    </div>
                    <div className="space-y-2">
                      <Label htmlFor="xlabel-input">X Axis Label</Label>
                      <Input id="xlabel-input" value={xAxisLabel} onChange={(e) => setXAxisLabel(e.target.value)} placeholder="Enter X axis label" />
                    </div>
                    <div className="space-y-2">
                      <Label htmlFor="ylabel-input">Y Axis Label</Label>
                      <Input id="ylabel-input" value={yAxisLabel} onChange={(e) => setYAxisLabel(e.target.value)} placeholder="Enter Y axis label" />
                    </div>
                  </div>

                  {/* Font Size Sliders & Grid Toggle */}
                  <div className="space-y-4">
                     <div className="space-y-2">
                        <Label htmlFor="title-fontsize">Title Font Size: {titleFontSize}px</Label>
                        <Slider id="title-fontsize" min={10} max={30} step={1} value={[titleFontSize]} onValueChange={(v) => setTitleFontSize(v[0])} />
                     </div>
                     <div className="space-y-2">
                        <Label htmlFor="axis-label-fontsize">Axis Label Font Size: {axisLabelFontSize}px</Label>
                        <Slider id="axis-label-fontsize" min={8} max={24} step={1} value={[axisLabelFontSize]} onValueChange={(v) => setAxisLabelFontSize(v[0])} />
                     </div>
                     <div className="space-y-2">
                        <Label htmlFor="tick-label-fontsize">Tick Label Font Size: {tickLabelFontSize}px</Label>
                        <Slider id="tick-label-fontsize" min={8} max={20} step={1} value={[tickLabelFontSize]} onValueChange={(v) => setTickLabelFontSize(v[0])} />
                     </div>

                     <div className="flex items-center justify-between space-x-2 pt-2">
                        <Label htmlFor="grid-toggle" className="cursor-pointer">Show Grid Lines</Label>
                        <Switch
                          id="grid-toggle"
                          checked={showGrid}
                          onCheckedChange={setShowGrid}
                        />
                     </div>

                     {/* Conditional Best Fit Line Toggle */}
                     {chartType === 'line' && (
                        <div className="space-y-2 pt-2">
                           <div className="flex items-center justify-between space-x-2">
                              <Label htmlFor="bestfit-toggle" className="cursor-pointer">Show Line of Best Fit</Label>
                              <Switch
                                id="bestfit-toggle"
                                checked={showBestFit}
                                onCheckedChange={setShowBestFit}
                              />
                           </div>
                           {/* Display Regression Results */}
                           {showBestFit && regressionSlope !== null && regressionIntercept !== null && rSquared !== null && (
                              <div className="text-xs text-muted-foreground pt-1 pl-1">
                                 <p>y = {regressionSlope.toFixed(4)}x + {regressionIntercept.toFixed(4)}</p>
                                 <p>R² = {rSquared.toFixed(4)}</p>
                              </div>
                           )}
                           {showBestFit && (regressionSlope === null || rSquared === null) && plotData.length >= 2 && (
                               <p className="text-xs text-destructive pt-1 pl-1">Could not calculate regression (check data).</p>
                           )}
                        </div>
                     )}
                  </div>

                </CardContent>
                {error && !plotData.length && ( // Show plot-related errors here
                  (<CardContent>
                    <p className="text-sm text-red-500">{error}</p>
                  </CardContent>)
                )}
              </Card>
            </div>
          )}

          {/* Chart Column - Render only when data is loaded */}
          {csvData.length > 0 && (
            <div className="lg:col-span-2">
              <Card>
                <CardContent className="pt-6">
                  <div className="relative h-[500px] md:h-[600px] rounded-md border bg-card">
                    {Object.keys(chartOptions).length > 0 ? (
                      <ReactECharts
                        ref={echartsRef}
                        echarts={echarts}
                        option={chartOptions}
                        style={{ height: '100%', width: '100%' }}
                        notMerge={true} // Force full redraw on option change
                        lazyUpdate={true} // Update lazily for performance
                      />
                    ) : (
                      <div className="absolute inset-0 flex items-center justify-center text-muted-foreground">
                        {error ? `Error: ${error}` : "Select X and Y axes to plot data."}
                      </div>
                    )}
                  </div>
                </CardContent>
              </Card>
            </div>
          )}

          {/* Placeholder and Drop Zone - Render only when NO data is loaded */}
          {csvData.length === 0 && (
            // This div spans the full width (lg:col-span-1 in a 1-column grid)
            // Add drag handlers HERE
            (<div
              className="lg:col-span-1 relative flex items-center justify-center h-[400px] border-2 border-dashed border-muted-foreground/50 rounded-lg bg-card"
              onDragEnter={handleDragEnter}
              onDragLeave={handleDragLeave}
              onDragOver={handleDragOver}
              onDrop={handleDrop}
            >
              {/* Text Content */}
              {!isDragging && (
                <div className="text-center text-muted-foreground">
                  <UploadCloud className="w-12 h-12 mx-auto mb-4" />
                  <p className="text-lg font-medium">Drag and drop your CSV file here</p>
                  {error && <p className="text-sm text-red-500 mt-4">{error}</p>}
                </div>
              )}
              {/* Drag Overlay - Show only when dragging over THIS specific element */}
              {isDragging && (
                <div className={cn(
                  "absolute bg-primary/20 backdrop-blur-sm flex flex-col items-center justify-center z-10 rounded-lg border-2 border-dashed border-primary pointer-events-none",
                  // Make it cover the placeholder box exactly
                  "inset-0"
                )}>
                  <UploadCloud className="w-16 h-16 text-primary mb-4 animate-bounce" />
                  <p className="text-xl font-semibold text-primary">Drop your .csv file here</p>
                </div>
              )}
            </div>)
          )}

        </div>
      </main>
    </TooltipProvider>
  );
}
