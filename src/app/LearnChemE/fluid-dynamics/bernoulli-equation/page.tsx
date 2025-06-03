// src/app/LearnChemE/fluid-dynamics/bernoulli-equation/page.tsx
'use client';

import React, { useState, useEffect, useMemo } from 'react';
import { Slider } from "@/components/ui/slider";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
// Import TeX as the default export from react-katex
import TeX from '@matejmazur/react-katex';
import 'katex/dist/katex.min.css';

const ATM_PRESSURE_PA = 101325; // Pa
const GRAVITY = 9.81; // m/s^2

const FIXED_FLUID_DENSITY = 1000; // kg/m^3 (Water)

export default function BernoulliPage() {
  // --- Inputs (Matching Image Inputs) ---
  const [d1Input, setD1Input] = useState<number>(0.05); // m, Input for Inlet Diameter
  const [d2, setD2] = useState<number>(0.029); // m, Input for Outlet Diameter, label as cm
  const [z1] = useState<number>(0); // m, Inlet Elevation (fixed at 0 for simplicity, Δz is the input)
  const [deltaZInput, setDeltaZInput] = useState<number>(4.5); // m, Input for Change in Height
  const [p1AbsKPaInput, setP1AbsKPaInput] = useState<number>(308); // kPa, Input for Inlet Absolute Pressure
  const [v1Input, setV1Input] = useState<number>(4.3); // m/s, Input for Inlet Velocity

  // --- Calculated Outputs ---
  const [z2, setZ2] = useState<number>(z1 + deltaZInput); // m, Calculated Outlet Elevation
  const [v1, setV1] = useState<number>(v1Input); // m/s, Set from input v1Input
  const [v2, setV2] = useState<number>(0); // m/s
  const [p1AbsKPa, setP1AbsKPa] = useState<number>(p1AbsKPaInput); // kPa, Set from input
  const [p2AbsKPa, setP2AbsKPa] = useState<number>(0); // kPa
  const [p1GaugeKPa, setP1GaugeKPa] = useState<number>(p1AbsKPaInput - ATM_PRESSURE_PA / 1000); // kPa
  const [p2GaugeKPa, setP2GaugeKPa] = useState<number>(0); // kPa
  const [massFlowRate, setMassFlowRate] = useState<number>(0); // kg/s
  const [volFlowRate, setVolFlowRate] = useState<number>(0); // m^3/s

  // Area calculations now use the dynamic d1Input and state d2
  const area1 = useMemo(() => Math.PI * (d1Input / 2) ** 2, [d1Input]);
  const area2 = useMemo(() => Math.PI * (d2 / 2) ** 2, [d2]);

  // --- Recalculation Logic ---
  useEffect(() => {
    // Update calculated states based on inputs
    setZ2(z1 + deltaZInput); // z1 is fixed at 0
    setV1(v1Input);
    setP1AbsKPa(p1AbsKPaInput);
    setP1GaugeKPa(p1AbsKPaInput - ATM_PRESSURE_PA / 1000);

    // Ensure diameters are positive for calculations
    const safeD1 = d1Input > 1e-9 ? d1Input : 1e-9;
    const safeD2 = d2 > 1e-9 ? d2 : 1e-9;

    // Calculate V2 using continuity equation (A1*v1 = A2*v2)
    const calculatedV2 = v1Input * (Math.PI * (safeD1 / 2) ** 2) / (Math.PI * (safeD2 / 2) ** 2); // v2 = v1 * (D1/D2)^2
    setV2(Number.isFinite(calculatedV2) ? calculatedV2 : 0);

    // Calculate P2 using Bernoulli Equation (Absolute Pressures)
    const p1AbsPa = p1AbsKPaInput * 1000; // Convert kPa to Pa
    const calculatedP2AbsPa =
      p1AbsPa +
      0.5 * FIXED_FLUID_DENSITY * (v1Input ** 2 - calculatedV2 ** 2) +
      FIXED_FLUID_DENSITY * GRAVITY * (z1 - (z1 + deltaZInput));

    // Ensure P2 is not below absolute zero (though negative gauge is fine)
    const safeP2AbsPa = Math.max(0, calculatedP2AbsPa);
    setP2AbsKPa(safeP2AbsPa / 1000); // Convert Pa back to kPa
    setP2GaugeKPa((safeP2AbsPa - ATM_PRESSURE_PA) / 1000); // Convert Pa to gauge kPa

    // Calculate Flow Rates
    const calculatedVolFlowRate = Math.PI * (safeD1 / 2) ** 2 * v1Input;
    setVolFlowRate(Number.isFinite(calculatedVolFlowRate) ? calculatedVolFlowRate : 0);
    setMassFlowRate(FIXED_FLUID_DENSITY * (Number.isFinite(calculatedVolFlowRate) ? calculatedVolFlowRate : 0));
  }, [d1Input, d2, deltaZInput, p1AbsKPaInput, v1Input, z1]);

  // --- Output Display Component ---
  const OutputDisplay = ({ label, unit, value }: { label: string; unit: string; value: number }) => (
    <div className="flex items-center justify-between space-x-2">
      <Label htmlFor={label} className="text-sm whitespace-nowrap">
        {label}
      </Label>
      <div className="flex items-center space-x-1">
        <Input
          id={label}
          readOnly
          value={
            Math.abs(value) > 1e5 || (Math.abs(value) < 1e-4 && value !== 0)
              ? value.toExponential(3)
              : value.toPrecision(3)
          }
          className="h-8 text-right text-sm font-mono tabular-nums w-32" // Added w-32 for consistent width
        />
        <span className="text-sm text-muted-foreground whitespace-nowrap min-w-[40px] text-left">{unit}</span>
      </div>
    </div>
  );

  // --- Bernoulli Equation String ---
  const bernoulliEq = `P_{in} + \\frac{1}{2}\\rho v_{in}^2 + \\rho g z_{in} = P_{out} + \\frac{1}{2}\\rho v_{out}^2 + \\rho g z_{out}`;

  // --- Calculate dynamic SVG coordinates based on elevation and diameters ---
  const scaleY = 10; // Pixel per meter (increased scale slightly)
  const svgHeight = 370; // Previously 300 - increased for more top space
  const svgWidth = 500;
  const y_svg_baseline = svgHeight - 50; // SVG y for a reference level (e.g., z=0)

  // Inlet pipe positioning (uses dynamic diameter)
  const inlet_pipe_x = 50;
  const inlet_pipe_width = 100;
  const inlet_diameter_px = d1Input * 100 * 2; // Convert diameter to pixels for visualization (arbitrary scale factor 2)
  const y_inlet_center_svg = y_svg_baseline - z1 * scaleY; // Center of pipe at z1 (which is 0)
  const y_inlet_top_svg = y_inlet_center_svg - inlet_diameter_px / 2;
  const y_inlet_bottom_svg = y_inlet_center_svg + inlet_diameter_px / 2;

  // Outlet pipe positioning (uses state diameter)
  const outlet_pipe_width = 100;
  const outlet_diameter_px = d2 * 100 * 2; // Convert diameter to pixels for visualization
  const y_outlet_center_svg = y_svg_baseline - z2 * scaleY; // Center of pipe at z2
  const y_outlet_top_svg = y_outlet_center_svg - outlet_diameter_px / 2;
  const y_outlet_bottom_svg = y_outlet_center_svg + outlet_diameter_px / 2;

  // Vertical section starts immediately after the inlet pipe ends horizontally
  const vertical_section_start_x = inlet_pipe_x + inlet_pipe_width;
  const vertical_section_width = inlet_diameter_px;
  const vertical_section_end_x = vertical_section_start_x + vertical_section_width;

  // Reducer starts at the end of the vertical pipe section and connects to the start of the outlet pipe
  const reducer_start_x = vertical_section_end_x;
  const outlet_pipe_x = reducer_start_x + 50;

  // Delta Z dimension line coordinates
  const delta_z_line_x = Math.max(outlet_pipe_x + outlet_pipe_width + 50, vertical_section_end_x + 50);
  const delta_z_line_y1 = y_inlet_center_svg;
  const delta_z_line_y2 = y_outlet_center_svg;

  // Coordinates for the P, u, D text
  const text_vertical_spacing = 20;
  const text_y_offset_from_pipe_top = 10;

  // Inlet text positioning
  const p1_text_x = inlet_pipe_x + inlet_pipe_width / 2 - 5;
  const p1_text_y_start = y_inlet_top_svg - text_y_offset_from_pipe_top - 2 * text_vertical_spacing;

  // Outlet text positioning
  const p2_text_x = outlet_pipe_x + outlet_pipe_width / 2;
  const p2_text_y_start = y_outlet_top_svg - text_y_offset_from_pipe_top - 2 * text_vertical_spacing;

  return (
    <div className="container mx-auto p-4 md:p-8">
      <style jsx global>{`
        .diagram-text {
          fill: black;
        }
        html.dark .diagram-text {
          fill: white;
        }
        .responsive-svg-container {
          width: 100%;
          padding-bottom: 74%;
          position: relative;
          overflow: hidden;
          box-sizing: border-box;
        }
        .responsive-svg-container svg {
          position: absolute;
          top: 0;
          left: 0;
          width: 100%;
          height: 100%;
        }
      `}</style>

      <div className="grid grid-cols-1 md:grid-cols-3 gap-6 lg:gap-8">
        {/* --- Input Controls Column --- */}
        <Card className="md:col-span-1">
          <CardHeader>
            <CardTitle>Inputs</CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            {/* Inlet Pipe Diameter */}
            <div className="space-y-3 mb-6">
              <Label htmlFor="d1Input" className="font-medium flex justify-between items-center">
                <span>Inlet pipe diameter (D<sub>in</sub>):</span>
                <span className="font-mono w-[80px] text-right tabular-nums">
                  {(d1Input * 100).toFixed(1)} cm
                </span>
              </Label>
              <Slider
                id="d1Input"
                min={0.1}
                max={30}
                step={0.1}
                value={[d1Input * 100]}
                onValueChange={(vals) => setD1Input(vals[0] / 100)}
                className="my-2"
              />
            </div>

            {/* Outlet Pipe Diameter */}
            <div className="space-y-3 mb-6">
              <Label htmlFor="d2" className="font-medium flex justify-between items-center">
                <span>Outlet pipe diameter (D<sub>out</sub>):</span>
                <span className="font-mono w-[80px] text-right tabular-nums">
                  {(d2 * 100).toFixed(1)} cm
                </span>
              </Label>
              <Slider
                id="d2"
                min={0.1}
                max={30}
                step={0.1}
                value={[d2 * 100]}
                onValueChange={(vals) => setD2(vals[0] / 100)}
                className="my-2"
              />
            </div>

            {/* Change in Height (Delta Z) */}
            <div className="space-y-3 mb-6">
              <Label htmlFor="deltaZInput" className="font-medium flex justify-between items-center">
                <span>Change in height (Δz):</span>
                <span className="font-mono w-[80px] text-right tabular-nums">
                  {deltaZInput.toFixed(1)} m
                </span>
              </Label>
              <Slider
                id="deltaZInput"
                min={0}
                max={20}
                step={0.1}
                value={[deltaZInput]}
                onValueChange={(vals) => setDeltaZInput(vals[0])}
                className="my-2"
              />
            </div>

            {/* Inlet Absolute Pressure */}
            <div className="space-y-3 mb-6">
              <Label htmlFor="p1AbsKPaInput" className="font-medium flex justify-between items-center">
                <span>Inlet pressure (P<sub>in</sub>):</span>
                <span className="font-mono w-[80px] text-right tabular-nums">
                  {p1AbsKPaInput.toFixed(1)} kPa
                </span>
              </Label>
              <Slider
                id="p1AbsKPaInput"
                min={1}
                max={1000}
                step={1}
                value={[p1AbsKPaInput]}
                onValueChange={(vals) => setP1AbsKPaInput(vals[0])}
                className="my-2"
              />
            </div>

            {/* Inlet Velocity */}
            <div className="space-y-3 mb-6">
              <Label htmlFor="v1Input" className="font-medium flex justify-between items-center">
                <span>Inlet velocity (u<sub>in</sub>):</span>
                <span className="font-mono w-[80px] text-right tabular-nums">
                  {v1Input.toFixed(1)} m/s
                </span>
              </Label>
              <Slider
                id="v1Input"
                min={0}
                max={20}
                step={0.1}
                value={[v1Input]}
                onValueChange={(vals) => setV1Input(vals[0])}
                className="my-2"
              />
            </div>
          </CardContent>
        </Card>

        {/* --- Diagram and Outputs Column --- */}
        <div className="md:col-span-2 space-y-6">
          {/* --- System Diagram --- */}
          <Card>
            <CardContent className="flex justify-center items-center min-h-[350px] p-2">
              <div className="responsive-svg-container">
                <svg
                  width={svgWidth}
                  height={svgHeight}
                  viewBox={`0 0 ${svgWidth} ${svgHeight}`}
                  xmlns="http://www.w3.org/2000/svg"
                  className="rounded-md"
                >
                  {/* Legend */}
                  <text x="20" y="30" fontSize="12" className="diagram-text">
                    D = pipe diameter
                  </text>
                  <text x="20" y="45" fontSize="12" className="diagram-text">
                    P = pressure (Abs)
                  </text>
                  <text x="20" y="60" fontSize="12" className="diagram-text">
                    u = fluid velocity
                  </text>

                  {/* Pipe sections */}
                  <rect
                    x={inlet_pipe_x}
                    y={y_inlet_top_svg}
                    width={inlet_pipe_width}
                    height={inlet_diameter_px}
                    fill="#cccccc"
                    stroke="#555555"
                    strokeWidth="1"
                  />
                  <rect
                    x={vertical_section_start_x}
                    y={Math.min(y_inlet_center_svg, y_outlet_center_svg) - vertical_section_width / 2}
                    width={vertical_section_width}
                    height={Math.abs(y_outlet_center_svg - y_inlet_center_svg) + vertical_section_width}
                    fill="#cccccc"
                    stroke="#555555"
                    strokeWidth="1"
                  />
                  <polygon
                    points={`
                      ${vertical_section_start_x + vertical_section_width},${y_outlet_center_svg - vertical_section_width / 2}
                      ${outlet_pipe_x},${y_outlet_top_svg}
                      ${outlet_pipe_x},${y_outlet_bottom_svg}
                      ${vertical_section_start_x + vertical_section_width},${y_outlet_center_svg + vertical_section_width / 2}
                    `}
                    fill="#cccccc"
                    stroke="#555555"
                    strokeWidth="1"
                  />
                  <rect
                    x={outlet_pipe_x}
                    y={y_outlet_top_svg}
                    width={outlet_pipe_width}
                    height={outlet_diameter_px}
                    fill="#cccccc"
                    stroke="#555555"
                    strokeWidth="1"
                  />

                  {/* Point 1 Text Outputs */}
                  <text
                    x={p1_text_x}
                    y={p1_text_y_start}
                    fontFamily="Merriweather Sans, monospace"
                    fontSize="12"
                    className="diagram-text"
                    textAnchor="middle"
                  >
                    D<tspan baselineShift="sub">in</tspan> = {(d1Input * 100).toFixed(1)} cm
                  </text>
                  <text
                    x={p1_text_x}
                    y={p1_text_y_start + text_vertical_spacing}
                    fontFamily="Merriweather Sans, monospace"
                    fontSize="12"
                    className="diagram-text"
                    textAnchor="middle"
                  >
                    P<tspan baselineShift="sub">in</tspan> = {p1AbsKPa.toFixed(1)} kPa
                  </text>
                  <text
                    x={p1_text_x}
                    y={p1_text_y_start + 2 * text_vertical_spacing}
                    fontFamily="Merriweather Sans, monospace"
                    fontSize="12"
                    className="diagram-text"
                    textAnchor="middle"
                  >
                    u<tspan baselineShift="sub">in</tspan> = {v1.toFixed(1)} m/s
                  </text>

                  {/* Point 2 Text Outputs */}
                  <text
                    x={p2_text_x}
                    y={p2_text_y_start}
                    fontFamily="Merriweather Sans, monospace"
                    fontSize="12"
                    className="diagram-text"
                    textAnchor="middle"
                  >
                    D<tspan baselineShift="sub">out</tspan> = {(d2 * 100).toFixed(1)} cm
                  </text>
                  <text
                    x={p2_text_x}
                    y={p2_text_y_start + text_vertical_spacing}
                    fontFamily="Merriweather Sans, monospace"
                    fontSize="12"
                    className="diagram-text"
                    textAnchor="middle"
                  >
                    P<tspan baselineShift="sub">out</tspan> = {p2AbsKPa.toFixed(1)} kPa
                  </text>
                  <text
                    x={p2_text_x}
                    y={p2_text_y_start + 2 * text_vertical_spacing}
                    fontFamily="Merriweather Sans, monospace"
                    fontSize="12"
                    className="diagram-text"
                    textAnchor="middle"
                  >
                    u<tspan baselineShift="sub">out</tspan> = {v2.toFixed(1)} m/s
                  </text>

                  {/* Delta Z dimension */}
                  <line
                    x1={delta_z_line_x}
                    y1={delta_z_line_y1}
                    x2={delta_z_line_x}
                    y2={delta_z_line_y2}
                    stroke="black"
                    strokeWidth="1"
                  />
                  <line
                    x1={inlet_pipe_x + inlet_pipe_width}
                    y1={delta_z_line_y1}
                    x2={delta_z_line_x}
                    y2={delta_z_line_y1}
                    stroke="black"
                    strokeDasharray="2,2"
                  />
                  <line
                    x1={outlet_pipe_x + outlet_pipe_width}
                    y1={delta_z_line_y2}
                    x2={delta_z_line_x}
                    y2={delta_z_line_y2}
                    stroke="black"
                    strokeDasharray="2,2"
                  />
                  <polygon
                    points={`${delta_z_line_x - 3},${delta_z_line_y2 + (delta_z_line_y1 > delta_z_line_y2 ? -5 : 5)} ${
                      delta_z_line_x + 3
                    },${delta_z_line_y2 + (delta_z_line_y1 > delta_z_line_y2 ? -5 : 5)} ${delta_z_line_x},${delta_z_line_y2}`}
                    fill="black"
                  />
                  <polygon
                    points={`${delta_z_line_x - 3},${delta_z_line_y1 + (delta_z_line_y1 > delta_z_line_y2 ? 5 : -5)} ${
                      delta_z_line_x + 3
                    },${delta_z_line_y1 + (delta_z_line_y2 ? 5 : -5)} ${delta_z_line_x},${delta_z_line_y1}`}
                    fill="black"
                  />
                  <text
                    x={delta_z_line_x + 5}
                    y={(delta_z_line_y1 + delta_z_line_y2) / 2}
                    fontSize="12"
                    fontWeight="bold"
                    className="diagram-text"
                    dominantBaseline="middle"
                    textAnchor="start"
                  >
                    Δz = {deltaZInput.toFixed(1)} m
                  </text>
                </svg>
              </div>
            </CardContent>
          </Card>

          {/* --- Main Calculated Outputs --- */}
          <Card>
            <CardHeader>
              <CardTitle>Calculated Values</CardTitle>
            </CardHeader>
            <CardContent className="space-y-3">
              <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
                <OutputDisplay label="Change in Pressure (ΔP)" unit="kPa" value={p2AbsKPa - p1AbsKPa} />
                <OutputDisplay label="Change in Velocity (Δu)" unit="m/s" value={v2 - v1} />
                <OutputDisplay label="Difference in Diameter (ΔD)" unit="cm" value={(d2 * 100) - (d1Input * 100)} />
                <OutputDisplay label="Change in Height (Δz)" unit="m" value={deltaZInput} />
                <OutputDisplay label="Mass Flow Rate (ṁ)" unit="kg/s" value={massFlowRate} />
                <OutputDisplay label="Volumetric Flow Rate (V̇)" unit="m³/s" value={volFlowRate} />
              </div>
            </CardContent>
          </Card>

          {/* --- Equation Display --- */}
          <Card>
            <CardHeader>
              <CardTitle>Governing Equation (Simplified Bernoulli)</CardTitle>
            </CardHeader>
            <CardContent>
              {/* Use TeX component with block prop for block math */}
              <TeX math={bernoulliEq} block />
              <p className="text-sm text-muted-foreground mt-2">
                {/* Use TeX component (default inline) for inline math */}
                Where <TeX math="P" /> is absolute pressure,{" "}
                <TeX math="\rho" /> is fluid density, <TeX math="v" /> is
                fluid velocity, <TeX math="g" /> is gravitational acceleration, and{" "}
                <TeX math="z" /> is elevation. This simplified form neglects friction and pump
                work.
              </p>
            </CardContent>
          </Card>
        </div>
      </div>
    </div>
  );
}