// src/app/LearnChemE/fluid-dynamics/bernoulli-equation/page.tsx
'use client';

import React, { useState, useEffect, useMemo } from 'react';
import { Slider } from "@/components/ui/slider";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Checkbox } from "@/components/ui/checkbox";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import BlockMath from '@matejmazur/react-katex';
import 'katex/dist/katex.min.css';

const ATM_PRESSURE_PA = 101325; // Pa
const GRAVITY = 9.81; // m/s^2

const parseFloatSafe = (value: string | number | readonly string[] | undefined, defaultValue = 0): number => {
  if (value === undefined || value === null) return defaultValue;
  const num = parseFloat(value.toString());
  return isNaN(num) ? defaultValue : num;
};


export default function BernoulliPage() {
  // --- Inputs ---
  const [d1, setD1] = useState<number>(0.1); // m
  const [d2, setD2] = useState<number>(0.05); // m
  const [z1, setZ1] = useState<number>(0); // m
  const [z2, setZ2] = useState<number>(2.9); // m - Default based on new diagram example
  const [p1GaugeKPa, setP1GaugeKPa] = useState<number>(100); // kPa - Default based on P_kjm3 = 100
  const [rhoMan, setRhoMan] = useState<number>(13600); // kg/m^3 (Mercury) - Less relevant now?
  const [hMan, setHMan] = useState<number>(0.1); // m - Less relevant now?
  const [rhoFluid, setRhoFluid] = useState<number>(1000); // kg/m^3 (Water)
  const [pumpWorkKJkg, setPumpWorkKJkg] = useState<number>(0); // kJ/kg - No pump in new diagram
  const [includeFriction, setIncludeFriction] = useState<boolean>(true); // Friction often included
  const [frictionLossKJkg, setFrictionLossKJkg] = useState<number>(10); // kJ/kg

  // --- Outputs ---
  const [v1, setV1] = useState<number>(0); // m/s
  const [v2, setV2] = useState<number>(0); // m/s
  const [p2GaugeKPa, setP2GaugeKPa] = useState<number>(0); // kPa
  const [p1AbsKPa, setP1AbsKPa] = useState<number>(0); // kPa
  const [p2AbsKPa, setP2AbsKPa] = useState<number>(0); // kPa
  const [massFlowRate, setMassFlowRate] = useState<number>(0); // kg/s
  const [volFlowRate, setVolFlowRate] = useState<number>(0); // m^3/s
  // --- NEW Outputs for Diagram (kJ/m^3) ---
  const [ke1, setKe1] = useState<number>(0);
  const [pe1, setPe1] = useState<number>(0);
  const [p1kjm3, setP1kjm3] = useState<number>(0);
  const [ke2, setKe2] = useState<number>(0);
  const [pe2, setPe2] = useState<number>(0);
  const [p2kjm3, setP2kjm3] = useState<number>(0);

  // --- Add a new state variable for the test slider ---
  const [testSliderVal, setTestSliderVal] = useState(50);

  const area1 = useMemo(() => Math.PI * (d1 / 2) ** 2, [d1]);
  const area2 = useMemo(() => Math.PI * (d2 / 2) ** 2, [d2]);

  // --- Recalculation Logic ---
  useEffect(() => {
    // Convert inputs to base SI units for calculation
    const p1GaugePa = p1GaugeKPa * 1000;
    const pumpWorkJkg = pumpWorkKJkg * 1000; // Note: Pump work input might be removed if not needed
    const frictionLossJkg = includeFriction ? frictionLossKJkg * 1000 : 0;

    // Calculate absolute pressure at 1 (needed for P display and Bernoulli)
    const p1AbsPa = p1GaugePa + ATM_PRESSURE_PA;
    setP1AbsKPa(p1AbsPa / 1000);
    setP1kjm3(p1AbsPa / 1000); // P in kJ/m3

    // Calculate PE at 1
    setPe1((rhoFluid * GRAVITY * z1) / 1000);

    // --- Solve Bernoulli ---
    // Using FIXED v1 = 1.0 m/s for now based on example diagram
    const EXAMPLE_V1 = 1.0;
    const calculatedV1 = EXAMPLE_V1;
    setV1(calculatedV1);
    setKe1(0.5 * rhoFluid * calculatedV1**2 / 1000); // KE1 calculation

    const diameterRatio = (d1 > 1e-9 && d2 > 1e-9) ? (d1 / d2) ** 2 : 1;
    const calculatedV2 = calculatedV1 * diameterRatio;
    setV2(Number.isFinite(calculatedV2) ? calculatedV2 : 0);
    setKe2(0.5 * rhoFluid * calculatedV2**2 / 1000); // KE2 calculation

    // Now calculate P2 using Bernoulli with known v1, v2
    let p2AbsPa = p1AbsPa + rhoFluid * (calculatedV1**2 - calculatedV2**2) / 2 + rhoFluid * GRAVITY * (z1 - z2) + rhoFluid*pumpWorkJkg - rhoFluid*frictionLossJkg;

    if (!Number.isFinite(p2AbsPa)) p2AbsPa = 0;

    setP2AbsKPa(p2AbsPa / 1000);
    setP2kjm3(p2AbsPa / 1000); // P2 in kJ/m3
    setP2GaugeKPa((p2AbsPa - ATM_PRESSURE_PA) / 1000);

    // Calculate PE at 2
    setPe2((rhoFluid * GRAVITY * z2) / 1000);

    // --- Calculate Flow Rates ---
    setMassFlowRate(rhoFluid * area1 * calculatedV1);
    setVolFlowRate(area1 * calculatedV1);

  }, [
      // Ensure dependencies are correct
      d1, d2, z1, z2, p1GaugeKPa, rhoFluid,
      pumpWorkKJkg, includeFriction, frictionLossKJkg, area1, area2
  ]);

  // --- Output Display Component ---
 const OutputDisplay = ({ label, unit, value }: { label: string; unit: string; value: number }) => (
    <div className="flex items-center justify-between space-x-2">
        <Label htmlFor={label} className="text-sm whitespace-nowrap">{label} ({unit})</Label>
        <Input
            id={label}
            readOnly
            value={Math.abs(value) > 1e6 || (Math.abs(value) < 1e-3 && value !== 0) ? value.toExponential(3) : value.toPrecision(4)}
            className="h-8 text-right text-sm font-mono tabular-nums"
        />
    </div>
);

  // --- Bernoulli Equation String ---
  const bernoulliEq = `\\frac{P_1}{\\rho g} + \\frac{v_1^2}{2g} + z_1 + h_{pump} = \\frac{P_2}{\\rho g} + \\frac{v_2^2}{2g} + z_2 + h_{friction}`;

  return (
    <div className="container mx-auto p-4 md:p-8">
      <h1 className="text-3xl font-bold mb-6 text-center">Bernoulli Equation Simulation</h1>

      <div className="grid grid-cols-1 md:grid-cols-3 gap-6 lg:gap-8">
        {/* --- Input Controls Column --- */}
        <Card className="md:col-span-1">
          <CardHeader>
            <CardTitle>Inputs</CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            {/* --- REPLACED InputSlider call with direct structure --- */}
            {/* Diameter 1 */}
            <div className="space-y-3">
              <Label htmlFor="d1" className="font-medium">
                Diameter 1 (D₁) (m): <span className="font-mono inline-block w-[80px] text-right tabular-nums">{d1.toFixed(3)}</span>
              </Label>
              <Slider
                id="d1"
                min={0.01}
                max={0.5}
                step={0.001}
                value={[d1]}
                onValueChange={(vals) => setD1(vals[0])}
                // onValueCommit={(vals) => setD1(vals[0])} // OR use commit if preferred
                className="my-2"
              />
            </div>

            {/* Diameter 2 */}
            <div className="space-y-3">
              <Label htmlFor="d2" className="font-medium">
                Diameter 2 (D₂) (m): <span className="font-mono inline-block w-[80px] text-right tabular-nums">{d2.toFixed(3)}</span>
              </Label>
              <Slider
                id="d2"
                min={0.01}
                max={0.5}
                step={0.001}
                value={[d2]}
                onValueChange={(vals) => setD2(vals[0])}
                // onValueCommit={(vals) => setD2(vals[0])} // OR use commit if preferred
                className="my-2"
              />
            </div>

            {/* Elevation 1 */}
            <div className="space-y-3">
              <Label htmlFor="z1" className="font-medium">
                Elevation 1 (z₁) (m): <span className="font-mono inline-block w-[80px] text-right tabular-nums">{z1.toFixed(1)}</span> {/* Adjusted precision */}
              </Label>
              <Slider
                id="z1"
                min={-10} max={20} step={0.1}
                value={[z1]}
                onValueChange={(vals) => setZ1(vals[0])}
                className="my-2"
              />
            </div>

            {/* Elevation 2 */}
            <div className="space-y-3">
              <Label htmlFor="z2" className="font-medium">
                Elevation 2 (z₂) (m): <span className="font-mono inline-block w-[80px] text-right tabular-nums">{z2.toFixed(1)}</span> {/* Adjusted precision */}
              </Label>
              <Slider
                id="z2"
                min={-10} max={20} step={0.1}
                value={[z2]}
                onValueChange={(vals) => setZ2(vals[0])}
                className="my-2"
              />
            </div>

            {/* Gauge Pressure 1 */}
            <div className="space-y-3">
              <Label htmlFor="p1GaugeKPa" className="font-medium">
                Gauge Pressure 1 (P₁) (kPa): <span className="font-mono inline-block w-[80px] text-right tabular-nums">{p1GaugeKPa.toFixed(1)}</span> {/* Adjusted precision */}
              </Label>
              <Slider
                id="p1GaugeKPa"
                min={0} max={500} step={1}
                value={[p1GaugeKPa]}
                onValueChange={(vals) => setP1GaugeKPa(vals[0])}
                className="my-2"
              />
            </div>

            {/* Fluid Density */}
            <div className="space-y-3">
              <Label htmlFor="rhoFluid" className="font-medium">
                Fluid Density (ρ) (kg/m³): <span className="font-mono inline-block w-[80px] text-right tabular-nums">{rhoFluid.toFixed(0)}</span> {/* Adjusted precision */}
              </Label>
              <Slider
                id="rhoFluid"
                min={500} max={2000} step={5}
                value={[rhoFluid]}
                onValueChange={(vals) => setRhoFluid(vals[0])}
                className="my-2"
              />
            </div>

            {/* Friction Checkbox and Slider */}
            <div className="flex items-center space-x-2 pt-2">
              <Checkbox
                id="includeFriction"
                checked={includeFriction}
                onCheckedChange={(checked) => setIncludeFriction(Boolean(checked))}
              />
              <Label htmlFor="includeFriction">Include Friction?</Label>
            </div>
            {includeFriction && (
               <div className="space-y-3">
                 <Label htmlFor="frictionLoss" className="font-medium">
                   Friction Loss (F) (kJ/kg): <span className="font-mono inline-block w-[80px] text-right tabular-nums">{frictionLossKJkg.toFixed(1)}</span> {/* Adjusted precision */}
                 </Label>
                 <Slider
                   id="frictionLoss"
                   min={0} max={50} step={0.1}
                   value={[frictionLossKJkg]}
                   onValueChange={(vals) => setFrictionLossKJkg(vals[0])}
                   className="my-2"
                 />
               </div>
            )}
          </CardContent>
        </Card>

        {/* --- Diagram and Outputs Column --- */}
        <div className="md:col-span-2 space-y-6">
            {/* --- NEW Diagram --- */}
             <Card>
                <CardHeader>
                    <CardTitle>System Diagram & Results</CardTitle>
                </CardHeader>
                <CardContent className="flex justify-center items-center min-h-[350px]">
                     <svg width="500" height="300" viewBox="0 0 500 300" xmlns="http://www.w3.org/2000/svg" className="rounded-md">
                        {/* Pipe sections */}
                        <rect x="50" y="180" width="100" height="20" fill="#cccccc" stroke="#555555" strokeWidth="1"/> {/* Inlet pipe */}
                        <rect x="150" y="100" width="20" height="80" fill="#cccccc" stroke="#555555" strokeWidth="1"/> {/* Vertical pipe */}
                        <rect x="170" y="100" width="150" height="12" fill="#cccccc" stroke="#555555" strokeWidth="1"/> {/* Reducer visually */}
                        <rect x="320" y="100" width="100" height="12" fill="#cccccc" stroke="#555555" strokeWidth="1"/> {/* Outlet pipe */}

                        {/* Elbows/Connectors */}
                        <path d="M 150 180 Q 150 100 170 100" stroke="#555555" strokeWidth="1" fill="none"/> {/* Bottom elbow */}
                         {/* Simplified reducer connection */}
                        <line x1="170" y1="112" x2="170" y2="100" stroke="#555555" strokeWidth="1"/>

                        {/* Flow Arrows */}
                        <path d="M 30 190 L 70 190 L 70 185 L 80 190 L 70 195 L 70 190 Z" fill="blue" opacity="0.6"/>
                        <path d="M 420 106 L 460 106 L 460 101 L 470 106 L 460 111 L 460 106 Z" fill="blue" opacity="0.6"/>

                        {/* Point 1 Text Outputs */}
                        <text x="50" y="225" fontFamily="monospace" fontSize="12">K.E. = {ke1.toFixed(1)} kJ/m³</text>
                        <text x="50" y="240" fontFamily="monospace" fontSize="12">P.E. = {pe1.toFixed(1)} kJ/m³</text>
                        <text x="50" y="255" fontFamily="monospace" fontSize="12">P    = {p1kjm3.toFixed(1)} kJ/m³</text>

                         {/* Point 2 Text Outputs */}
                        <text x="350" y="130" fontFamily="monospace" fontSize="12">K.E. = {ke2.toFixed(1)} kJ/m³</text>
                        <text x="350" y="145" fontFamily="monospace" fontSize="12">P.E. = {pe2.toFixed(1)} kJ/m³</text>
                        <text x="350" y="160" fontFamily="monospace" fontSize="12">P    = {p2kjm3.toFixed(1)} kJ/m³</text>

                        {/* Delta Z dimension */}
                         <line x1="180" y1="180" x2="200" y2="180" stroke="black" strokeDasharray="2,2"/> {/* Line from bottom pipe */}
                        <line x1="180" y1="106" x2="200" y2="106" stroke="black" strokeDasharray="2,2"/> {/* Line from top pipe */}
                         <line x1="190" y1="180" x2="190" y2="106" stroke="black"/> {/* Vertical line */}
                         <polygon points="187,116 193,116 190,106" fill="black"/> {/* Arrow top */}
                         <polygon points="187,170 193,170 190,180" fill="black"/> {/* Arrow bottom */}
                         <text x="205" y="145" fontSize="12" fontWeight="bold">Δz = {(z2 - z1).toFixed(1)} m</text>

                        {/* Legend */}
                         <text x="50" y="20" fontSize="12">K.E. = Kinetic energy</text>
                         <text x="50" y="35" fontSize="12">P.E. = Potential energy</text>
                         <text x="50" y="50" fontSize="12">P    = Pressure</text>
                     </svg>
                </CardContent>
            </Card>

            {/* --- Additional Outputs (Optional Now) --- */}
             {/* You might still want to display v1, v2, flow rates etc. here */}
             <Card>
                <CardHeader>
                    <CardTitle>Other Calculated Values</CardTitle>
                </CardHeader>
                <CardContent className="space-y-3">
                   <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
                        <OutputDisplay label="Velocity 1 (v₁)" unit="m/s" value={v1} />
                        <OutputDisplay label="Velocity 2 (v₂)" unit="m/s" value={v2} />
                        <OutputDisplay label="Gauge Pressure 2 (P₂)" unit="kPa" value={p2GaugeKPa} />
                        {/* <OutputDisplay label="Absolute Pressure 1 (P₁)" unit="kPa" value={p1AbsKPa} /> */}
                        {/* <OutputDisplay label="Absolute Pressure 2 (P₂)" unit="kPa" value={p2AbsKPa} /> */}
                        <OutputDisplay label="Mass Flow Rate (ṁ)" unit="kg/s" value={massFlowRate} />
                        <OutputDisplay label="Volumetric Flow Rate (V̇)" unit="m³/s" value={volFlowRate} />
                   </div>
                </CardContent>
            </Card>

            {/* --- Equation Display --- */}
            <Card>
                <CardHeader>
                    <CardTitle>Governing Equation</CardTitle>
                </CardHeader>
                <CardContent>
                     <BlockMath math={bernoulliEq} />
                     <p className="text-sm text-muted-foreground mt-2">
                         Where <code className="font-mono">$h_{"{pump}"} = W_s / g$</code> and <code className="font-mono">$h_{"{friction}"} = F / g$</code>. Pressures <code className="font-mono">$P_1, P_2$</code> are absolute.
                     </p>
                </CardContent>
            </Card>
        </div>
      </div>
    </div>
  );
}