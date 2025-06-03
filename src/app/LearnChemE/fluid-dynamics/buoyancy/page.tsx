'use client';

import React, { useState, useEffect, useRef } from 'react';
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Slider } from "@/components/ui/slider";
import { Label } from "@/components/ui/label";
import { Input } from "@/components/ui/input";
import TeX from '@matejmazur/react-katex';
import 'katex/dist/katex.min.css';

const GRAVITY = 9.81; // m/s^2
const DENSITY_WATER = 1000; // kg/m³
const L_FIXED = 0.1; // Fixed actual side length of the cube in meters

const OutputDisplay = ({ label, unit, value, texLabel, className }: { label: string; unit: string; value: string | number; texLabel?: string; className?: string }) => {
  const spanRef = useRef<HTMLSpanElement>(null);
  const [inputWidth, setInputWidth] = useState<number>(0);

  useEffect(() => {
    if (spanRef.current) {
      setInputWidth(spanRef.current.offsetWidth + 25); // Add padding to prevent cutting off the last number
    }
  }, [value]);

  return (
    <div className="flex items-center justify-between space-x-2 p-2 border-b last:border-b-0">
      <Label htmlFor={label.replace(/\s+/g, '-')} className="text-sm whitespace-nowrap">
        {texLabel ? <TeX math={texLabel} /> : label}
      </Label>
      <div className="flex items-center space-x-1">
        {/* Hidden span to measure text width */}
        <span ref={spanRef} className="absolute invisible whitespace-nowrap">
          {value}
        </span>
        <Input
          id={label.replace(/\s+/g, '-')}
          readOnly
          value={typeof value === 'number' ? value.toPrecision(3) : value}
          className={`h-8 text-right text-sm font-mono tabular-nums`}
          style={{ width: `${inputWidth}px` }} // Dynamically set width
        />
        <span className="text-sm text-muted-foreground whitespace-nowrap min-w-[40px] text-left">{unit}</span>
      </div>
    </div>
  );
};

export default function BuoyancyPage() {
  const [objectSG, setObjectSG] = useState<number>(0.5); // Specific Gravity of Cube
  const [liquidSG, setLiquidSG] = useState<number>(1.0); // Specific Gravity of Fluid

  const [submergedDepth, setSubmergedDepth] = useState<number>(0); // m
  const [cubeWeight, setCubeWeight] = useState<number>(0); // N
  const [buoyantForce, setBuoyantForce] = useState<number>(0); // N
  const [status, setStatus] = useState<string>("Floats");

  useEffect(() => {
    const L = L_FIXED; // Use fixed side length
    const rho_cube = objectSG * DENSITY_WATER;
    const rho_fluid = liquidSG * DENSITY_WATER;

    const W = rho_cube * Math.pow(L, 3) * GRAVITY;
    setCubeWeight(W);

    let h_submerged: number;
    let currentStatus: string;

    if (rho_cube < rho_fluid) {
      h_submerged = (rho_cube / rho_fluid) * L;
      currentStatus = "Floats";
    } else if (rho_cube === rho_fluid) {
      h_submerged = L;
      currentStatus = "Neutrally Buoyant";
    } else { // rho_cube > rho_fluid
      h_submerged = L; // Fully submerged
      currentStatus = "Sinks";
    }
    h_submerged = Math.min(h_submerged, L);
    setSubmergedDepth(h_submerged);
    setStatus(currentStatus);

    const V_submerged = Math.pow(L, 2) * h_submerged;
    const F_B = rho_fluid * V_submerged * GRAVITY;
    setBuoyantForce(F_B);

  }, [objectSG, liquidSG]);

  const instructionalVideoUrl = "https://learncheme.github.io/website_resources/youtube_iframe.html?id=7qQimrfB5Fk";
  const openInstructionalVideo = () => {
    window.open(instructionalVideoUrl, 'Instructional Video', 'location=no,toolbar=no,menubar=no,resizable=yes,height=420,width=580,scrollbars=no,directories=no,status=no');
  };

  const svgWidth = 300;
  const svgHeight = 350;
  const fluidHeightRatio = 0.7;
  const fluidLevelY = svgHeight * (1 - fluidHeightRatio);
  const fluidVisualHeight = svgHeight * fluidHeightRatio;

  const CUBE_VISUAL_SIDE_SVG = 80; // Fixed visual size of the cube in SVG pixels

  const submergedDisplayHeight = (submergedDepth / L_FIXED) * CUBE_VISUAL_SIDE_SVG;

  let cubeYPosition: number;

  if (status === "Sinks") {
    cubeYPosition = svgHeight - CUBE_VISUAL_SIDE_SVG; // Bottom of cube at SVG bottom
  } else {
    cubeYPosition = fluidLevelY - (CUBE_VISUAL_SIDE_SVG - submergedDisplayHeight);
  }

  cubeYPosition = Math.max(0, cubeYPosition);
  if (cubeYPosition + CUBE_VISUAL_SIDE_SVG > svgHeight) {
    cubeYPosition = svgHeight - CUBE_VISUAL_SIDE_SVG;
  }

  const equations = [
    { name: "Density Calculation", formula: `\\rho_{obj} = SG_{obj} \\cdot \\rho_{water}` },
    { name: "Density Calculation", formula: `\\rho_{fluid} = SG_{fluid} \\cdot \\rho_{water}` },
    { name: "Weight of Cube", formula: `W = \\rho_{obj} \\cdot L^3 \\cdot g` },
    { name: "Buoyant Force", formula: `F_B = \\rho_{fluid} \\cdot V_{submerged} \\cdot g` },
    { name: "Submerged Volume", formula: `V_{submerged} = L^2 \\cdot h` },
    { name: "Submerged Depth (h)", formula: `\\text{If } \\rho_{obj} \\le \\rho_{fluid}, h = \\frac{\\rho_{obj}}{\\rho_{fluid}} L \\text{, else } h = L` },
  ];

  // Initial static path for the water body
  const initialWaterPathD = `
    M 0 ${fluidLevelY}
    Q ${svgWidth / 2} ${fluidLevelY} ${svgWidth} ${fluidLevelY}
    L ${svgWidth} ${svgHeight}
    L 0 ${svgHeight}
    Z
  `;

  // Condition for showing Delta H indicator
  const showDeltaHIndicator = status === "Floats" && cubeYPosition < fluidLevelY;
  const deltaHSvg = fluidLevelY - cubeYPosition; // Visual height difference in SVG units
  // Calculate real-world delta H in centimeters
  const deltaHRealCm = deltaHSvg * (L_FIXED / CUBE_VISUAL_SIDE_SVG) * 100;

  return (
    <div className="container mx-auto p-4 md:p-8">
      <header className="mb-8">
      </header>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        <div className="lg:col-span-1 space-y-6">
          <Card>
            <CardHeader>
              <CardTitle>Inputs</CardTitle>
            </CardHeader>
            <CardContent className="space-y-6 pt-2">
              <div className="space-y-3">
                <Label htmlFor="objectSG" className="font-medium flex justify-between">
                  <span>Object Specific Gravity (<TeX math="SG_{obj}" />)</span> <span>{objectSG.toFixed(2)}</span>
                </Label>
                <Slider id="objectSG" min={0.1} max={2.0} step={0.01} value={[objectSG]} onValueChange={(vals) => setObjectSG(vals[0])} />
              </div>
              <div className="space-y-3">
                <Label htmlFor="liquidSG" className="font-medium flex justify-between">
                  <span>Liquid Specific Gravity (<TeX math="SG_{fluid}" />)</span> <span>{liquidSG.toFixed(2)}</span>
                </Label>
                <Slider id="liquidSG" min={0.5} max={1.5} step={0.01} value={[liquidSG]} onValueChange={(vals) => setLiquidSG(vals[0])} />
              </div>
            </CardContent>
          </Card>

          <Card>
            <CardHeader>
              <CardTitle>Governing Equations</CardTitle>
            </CardHeader>
            <CardContent className="space-y-3 text-sm">
              {equations.map(eq => (
                <div key={eq.formula + eq.name}>
                  <p className="font-medium">{eq.name}:</p>
                  <TeX math={eq.formula} block className="p-2 bg-muted/30 rounded-md text-center overflow-x-auto"/>
                </div>
              ))}
            </CardContent>
          </Card>
        </div>

        <div className="lg:col-span-2 space-y-6">
          <Card>
            <CardHeader>
              <CardTitle>Visualization</CardTitle>
            </CardHeader>
            <CardContent className="flex justify-center items-center min-h-[380px] p-2">
              <svg width={svgWidth} height={svgHeight} viewBox={`0 0 ${svgWidth} ${svgHeight}`} className="border rounded-md bg-background">
                {/* Fluid */}
                <path d={initialWaterPathD} fill="#a7d7f9" />
                
                {/* Cube */}
                <rect
                  x={(svgWidth - CUBE_VISUAL_SIDE_SVG) / 2}
                  y={cubeYPosition}
                  width={CUBE_VISUAL_SIDE_SVG}
                  height={CUBE_VISUAL_SIDE_SVG}
                  fill="#facc15"
                  stroke="#ca8a04"
                  strokeWidth="2"
                />
                {/* Water Level Line */}
                <line x1="0" y1={fluidLevelY} x2={svgWidth} y2={fluidLevelY} stroke="#60a5fa" strokeWidth="1" strokeDasharray="4 2"/>

                {/* Delta H Indicator */}
                {showDeltaHIndicator && deltaHSvg > 1 && (
                  <g stroke="orange" strokeWidth="1.5">
                    {/* Vertical line for Delta H */}
                    <line
                      x1={(svgWidth + CUBE_VISUAL_SIDE_SVG) / 2 + 15}
                      y1={cubeYPosition}
                      x2={(svgWidth + CUBE_VISUAL_SIDE_SVG) / 2 + 15}
                      y2={fluidLevelY}
                    />
                    {/* Top tick */}
                    <line
                      x1={(svgWidth + CUBE_VISUAL_SIDE_SVG) / 2 + 10}
                      y1={cubeYPosition}
                      x2={(svgWidth + CUBE_VISUAL_SIDE_SVG) / 2 + 20}
                      y2={cubeYPosition}
                    />
                    {/* Bottom tick */}
                    <line
                      x1={(svgWidth + CUBE_VISUAL_SIDE_SVG) / 2 + 10}
                      y1={fluidLevelY}
                      x2={(svgWidth + CUBE_VISUAL_SIDE_SVG) / 2 + 20}
                      y2={fluidLevelY}
                    />
                    <text
                      x={(svgWidth + CUBE_VISUAL_SIDE_SVG) / 2 + 25}
                      y={cubeYPosition + deltaHSvg / 2}
                      fill="orange"
                      fontSize="10"
                      dominantBaseline="middle"
                    >
                      Δh = {deltaHRealCm.toFixed(1)} cm
                    </text>
                  </g>
                )}
              </svg>
            </CardContent>
          </Card>

          <Card>
            <CardHeader>
              <CardTitle>Calculated Values</CardTitle>
            </CardHeader>
            <CardContent className="pt-0">
              <OutputDisplay label="Cube Weight" texLabel="W" unit="N" value={cubeWeight.toFixed(2)} />
              <OutputDisplay label="Buoyant Force" texLabel="F_B" unit="N" value={buoyantForce.toFixed(2)} />
              <OutputDisplay label="Submerged Depth" texLabel="h" unit="cm" value={(submergedDepth * 100).toFixed(1)} />
              <OutputDisplay label="Object Density" texLabel="\rho_{obj}" unit="kg/m³" value={(objectSG * DENSITY_WATER).toFixed(0)} />
              <OutputDisplay label="Fluid Density" texLabel="\rho_{fluid}" unit="kg/m³" value={(liquidSG * DENSITY_WATER).toFixed(0)} />
              <OutputDisplay label="Status" unit="" value={status} />
            </CardContent>
          </Card>
        </div>
      </div>
    </div>
  );
}
