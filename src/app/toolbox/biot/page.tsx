'use client';

import { useState, useEffect, useCallback, useRef, useMemo } from 'react';
import { useTheme } from "next-themes";
import * as THREE from 'three';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, Center } from '@react-three/drei';

// Import ECharts components
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsOption } from 'echarts';
import { LineChart } from 'echarts/charts';
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  MarkLineComponent
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

echarts.use([
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  MarkLineComponent,
  LineChart,
  CanvasRenderer
]);

import { Card, CardContent } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Slider } from "@/components/ui/slider";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { TooltipProvider } from "@/components/ui/tooltip";
import { Play, Pause, RotateCcw } from 'lucide-react';

// --- Types ---
type GeometryType = 'plane' | 'cylinder' | 'sphere' | 'fin' | 'hollow_cylinder';

type SimulationState = {
  temps: number[];
  time: number;
  isSteady: boolean;
  steadyTime?: number | null;
};

// --- Constants ---
const NODES = 21;
const DR_REAL = 0.01;

const formatDuration = (secondsRaw: number): string => {
    const total = Math.round(secondsRaw);
    const hrs = Math.floor(total / 3600);
    const mins = Math.floor((total % 3600) / 60);
    const secs = total % 60;
    if (hrs > 0) return `${hrs}h ${mins}m ${secs}s`;
    if (mins > 0) return `${mins}m ${secs}s`;
    return `${secs}s`;
};

// --- SHADER CODE ---
const vertexShader = `
varying vec3 vPos;
varying vec3 vNormal;

void main() {
  vPos = position; 
  vNormal = normalize(normalMatrix * normal);
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}
`;

const fragmentShader = `
uniform float uTemps[21];
uniform float uMin;
uniform float uMax;
uniform float uType; // 0=Plane, 1=Cyl, 2=Sphere, 3=Fin, 4=HollowCylinder

varying vec3 vPos;
varying vec3 vNormal;

vec3 heatMap(float t) {
    t = clamp(t, 0.0, 1.0);
    vec3 blue   = vec3(0.0, 0.0, 1.0);
    vec3 cyan   = vec3(0.0, 1.0, 1.0);
    vec3 green  = vec3(0.0, 1.0, 0.0);
    vec3 yellow = vec3(1.0, 1.0, 0.0);
    vec3 red    = vec3(1.0, 0.0, 0.0);

    if (t < 0.25) return mix(blue, cyan, t * 4.0);
    if (t < 0.50) return mix(cyan, green, (t - 0.25) * 4.0);
    if (t < 0.75) return mix(green, yellow, (t - 0.50) * 4.0);
    return mix(yellow, red, (t - 0.75) * 4.0);
}

void main() {
    float r = 0.0;
    
    // --- GEOMETRY MAPPING ---
    if (uType < 0.5) { 
        // PLANE (Type 0)
        r = abs(vPos.x); 
    } 
    else if (uType < 1.5) { 
        // CYLINDER (Type 1)
        r = length(vPos.xz);
    } 
    else if (uType < 2.5) { 
        // SPHERE (Type 2)
        r = length(vPos);
        if (r > 1.0) discard;
    }
    else if (uType < 3.5) {
        // PIN FIN (Type 3)
        // Map X axis (-1 to 1) to (0 to 1) for the gradient
        // Fin base is at -1 (Hot), Tip is at +1 (Cold)
        r = (vPos.x + 1.0) * 0.5;
    }
    else {
        // HOLLOW CYLINDER (Type 4)
        r = length(vPos.xz);
        // Discard the center to make it hollow (Inner Radius = 0.5)
        if (r < 0.5) discard;
        
        // Re-normalize r so the color gradient (0 to 1) fits the wall thickness (0.5 to 1.0)
        r = (r - 0.5) * 2.0;
    }

    r = clamp(r, 0.0, 1.0);
    float index = r * 20.0;
    int i = int(index);
    float f = fract(index);
    if (i >= 20) i = 19;
    
    float t1 = uTemps[i];
    float t2 = uTemps[i+1];
    float tempVal = mix(t1, t2, f);
    
    float nT = (tempVal - uMin) / (uMax - uMin);
    vec3 baseColor = heatMap(nT);

    vec3 lightDir = normalize(vec3(1.0, 1.0, 1.0)); 
    float diff = max(dot(vNormal, lightDir), 0.2); 
    
    vec3 finalColor = baseColor * (diff * 0.8 + 0.4); 

    gl_FragColor = vec4(finalColor, 1.0);
}
`;

const ActiveHeatMaterial = ({ temps, min, max, type }: { temps: number[], min: number, max: number, type: number }) => {
    const matRef = useRef<THREE.ShaderMaterial>(null);
    const uniforms = useMemo(() => ({
        uTemps: { value: temps },
        uMin: { value: min },
        uMax: { value: max },
        uType: { value: type }
    }), []);

    useFrame(() => {
        if (matRef.current) {
            matRef.current.uniforms.uTemps.value = temps;
            matRef.current.uniforms.uMin.value = min;
            matRef.current.uniforms.uMax.value = max;
            matRef.current.uniforms.uType.value = type;
        }
    });

    return (
        <shaderMaterial
            ref={matRef}
            vertexShader={vertexShader}
            fragmentShader={fragmentShader}
            uniforms={uniforms}
            side={THREE.DoubleSide}
        />
    );
};

const HeatMesh = ({ geometryType, temps, minTemp, maxTemp }: { geometryType: GeometryType, temps: number[], minTemp: number, maxTemp: number }) => {
    
    const typeId = geometryType === 'plane' ? 0.0 : (geometryType === 'cylinder' ? 1.0 : (geometryType === 'sphere' ? 2.0 : (geometryType === 'fin' ? 3.0 : 4.0)) );
    const rMin = Math.min(minTemp, maxTemp);
    const rMax = Math.max(minTemp, maxTemp);

    const { cylFill1, cylFill2, sphereFill1, sphereFill2 } = useMemo(() => {
        const matrix = new THREE.Matrix4();
        
        // --- 1. BASE FLAP (Rectangle) ---
        // We use this SAME base for both Cylinder and Sphere to ensure perfect alignment.
        // For Sphere, the shader simply discards the corners to make it round.
        const baseFlap = new THREE.PlaneGeometry(1, 2, 20, 20); // Using segments > 1 helps lighting
        matrix.makeTranslation(0.5, 0, 0); // Extend 0 to 1 on X
        baseFlap.applyMatrix4(matrix);

        // --- 2. CYLINDER POSITIONING ---
        const c1 = baseFlap.clone();
        matrix.makeRotationY(-Math.PI / 2); // -90 deg
        c1.applyMatrix4(matrix);

        const c2 = baseFlap.clone();
        matrix.makeRotationY(Math.PI); // 180 deg — close the 270° cut
        c2.applyMatrix4(matrix);

        // --- 3. SPHERE POSITIONING ---
        // Use CircleGeometry semi-circle on +X side so no vertex translation needed
        const sphereBase = new THREE.CircleGeometry(1, 32, -Math.PI / 2, Math.PI);

        const s1 = sphereBase.clone();
        matrix.makeRotationY(Math.PI / 2); // 90 deg
        s1.applyMatrix4(matrix);

        const s2 = sphereBase.clone();
        matrix.makeRotationY(Math.PI); // 180 deg
        s2.applyMatrix4(matrix);

        return {
            cylFill1: c1, cylFill2: c2,
            sphereFill1: s1, sphereFill2: s2
        };
    }, []);

    return (
        <group rotation={[0, Math.PI / 4, 0]}>
            {geometryType === 'plane' && (
                <mesh>
                    <boxGeometry args={[2, 1, 1]} />
                    <ActiveHeatMaterial temps={temps} min={rMin} max={rMax} type={typeId} />
                </mesh>
            )}

            {geometryType === 'cylinder' && (
                <>
                    <mesh>
                        <cylinderGeometry args={[1, 1, 2, 64, 1, false, 0, Math.PI * 1.5]} />
                        <ActiveHeatMaterial temps={temps} min={rMin} max={rMax} type={typeId} />
                    </mesh>
                    <mesh geometry={cylFill1}>
                        <ActiveHeatMaterial temps={temps} min={rMin} max={rMax} type={typeId} />
                    </mesh>
                    <mesh geometry={cylFill2}>
                        <ActiveHeatMaterial temps={temps} min={rMin} max={rMax} type={typeId} />
                    </mesh>
                </>
            )}

            {geometryType === 'sphere' && (
                <>
                    <mesh>
                        <sphereGeometry args={[1, 64, 32, 0, Math.PI * 1.5, 0, Math.PI]} />
                        <ActiveHeatMaterial temps={temps} min={rMin} max={rMax} type={typeId} />
                    </mesh>
                    <mesh geometry={sphereFill1}>
                        <ActiveHeatMaterial temps={temps} min={rMin} max={rMax} type={typeId} />
                    </mesh>
                    <mesh geometry={sphereFill2}>
                        <ActiveHeatMaterial temps={temps} min={rMin} max={rMax} type={typeId} />
                    </mesh>
                </>
            )}

            {/* PIN FIN (Rod) - maps along X axis */}
            {geometryType === 'fin' && (
                <mesh rotation={[0, 0, -Math.PI / 2]}>
                    <cylinderGeometry args={[0.5, 0.5, 2, 64]} />
                    <ActiveHeatMaterial temps={temps} min={rMin} max={rMax} type={typeId} />
                </mesh>
            )}

            {/* HOLLOW CYLINDER (Pipe) - shader makes the inner hole */}
            {geometryType === 'hollow_cylinder' && (
                <>
                    <mesh>
                        <cylinderGeometry args={[1, 1, 2, 64, 1, false, 0, Math.PI * 1.5]} />
                        <ActiveHeatMaterial temps={temps} min={rMin} max={rMax} type={typeId} />
                    </mesh>
                    <mesh geometry={cylFill1}>
                        <ActiveHeatMaterial temps={temps} min={rMin} max={rMax} type={typeId} />
                    </mesh>
                    <mesh geometry={cylFill2}>
                        <ActiveHeatMaterial temps={temps} min={rMin} max={rMax} type={typeId} />
                    </mesh>
                </>
            )}
        </group>
    );
};


export default function HeatDiffusionPage() {
    const { resolvedTheme } = useTheme();

    const [geometry, setGeometry] = useState<GeometryType>('sphere');
    const [k, setK] = useState(50);
    const [h, setH] = useState(200);
    const [alphaScaled, setAlphaScaled] = useState(10);
    const [tempInit, setTempInit] = useState(25);
    const [tempInf, setTempInf] = useState(100);

    // Local input strings to allow empty input without overwriting numeric state with 0
    const [tempInitInput, setTempInitInput] = useState<string>(tempInit.toString());
    const [tempInfInput, setTempInfInput] = useState<string>(tempInf.toString());

    useEffect(() => { setTempInitInput(tempInit.toString()); }, [tempInit]);
    useEffect(() => { setTempInfInput(tempInf.toString()); }, [tempInf]);

    const [isRunning, setIsRunning] = useState(true);
    const [simState, setSimState] = useState<SimulationState>({
        temps: Array(NODES).fill(25),
        time: 0,
        isSteady: false,
        steadyTime: null
    });

    const requestRef = useRef<number>(0);
    const simDataRef = useRef<number[]>(Array(NODES).fill(25));
    const timeRef = useRef<number>(0);
    const steadyTimeRef = useRef<number | null>(null);

    // Derived Calculations
    const Lc = geometry === 'plane' ? DR_REAL * NODES : (geometry === 'cylinder' ? (DR_REAL * NODES) / 2 : (DR_REAL * NODES) / 3);
    const biotNumber = (h * Lc) / k;
    const alpha = alphaScaled * 1e-6;
    
    // --- ROBUST STABILITY CONTROLLER ---
    
    // 1. Calculate Local Biot (Surface Constraint)
    const localBi = (h * DR_REAL) / k;

    // 2. Determine Geometric Stability Limit (Center Constraint)
    // Sphere is the most fragile (1/6), then Cylinder (1/4), Plane is toughest (1/2)
    let maxFo_Geometry = 0.5;
    if (geometry === 'cylinder' || geometry === 'fin' || geometry === 'hollow_cylinder') maxFo_Geometry = 0.25;
    if (geometry === 'sphere') maxFo_Geometry = 1.0 / 6.0;

    // 3. Calculate Critical Time Steps
    // Limit A: Surface Stability (prevents crash at high h / low k)
    const limitSurface = (DR_REAL * DR_REAL) / (2 * alpha * (1 + localBi));
    
    // Limit B: Center/Interior Stability (prevents "Dot" artifact and crash at low h)
    const limitCenter = (maxFo_Geometry * DR_REAL * DR_REAL) / alpha;

    // 4. Pick the "Weakest Link" (Smallest allowable dt)
    // We use 0.9 safety factor to stay safely under the theoretical limit
    const dt = Math.min(limitSurface, limitCenter) * 0.9;

    const updateSimulation = useCallback(() => {
        if (!isRunning) return;

        const currentTemps = [...simDataRef.current];
        const newTemps = [...currentTemps];
        const Fo = (alpha * dt) / (DR_REAL * DR_REAL);
        const Bi = (h * DR_REAL) / k;

        // FIN Solver (1D fin losing to ambient along its side)
        if (geometry === 'fin') {
            const D = 1.0; // assumed diameter
            const sinkConst = (alpha * dt) * (4 * h) / (k * D);

            // Node 0: Base (fixed to wall temp)
            newTemps[0] = tempInit;

            // Interior Nodes
            for (let i = 1; i < NODES - 1; i++) {
                const conduction = Fo * (currentTemps[i + 1] + currentTemps[i - 1]) + (1 - 2 * Fo) * currentTemps[i];
                const convectionLoss = sinkConst * (currentTemps[i] - tempInf);
                newTemps[i] = conduction - convectionLoss;
            }

            // Tip (insulated)
            newTemps[NODES - 1] = newTemps[NODES - 2];
        }
        // HOLLOW CYLINDER Solver (inner surface at node 0)
        else if (geometry === 'hollow_cylinder') {
            const r_inner = 0.5;
            // Node 0: inner surface fixed to tempInit
            newTemps[0] = tempInit;

            for (let i = 1; i < NODES - 1; i++) {
                const r = r_inner + (i * DR_REAL);
                const termMinus = 1 - (DR_REAL) / (2 * r);
                const termPlus = 1 + (DR_REAL) / (2 * r);
                newTemps[i] = Fo * (termMinus * currentTemps[i - 1] + termPlus * currentTemps[i + 1]) + (1 - 2 * Fo) * currentTemps[i];
            }

            // Outer Surface (convection)
            const T_inner = currentTemps[NODES - 2];
            const T_surf = currentTemps[NODES - 1];
            newTemps[NODES - 1] = (2 * Fo * T_inner + 2 * Fo * Bi * tempInf + (1 - 2 * Fo - 2 * Fo * Bi) * T_surf);
        }
        // Standard Solvers (plane/cylinder/sphere)
        else {
            const n = geometry === 'plane' ? 0 : (geometry === 'cylinder' ? 1 : 2);

            // 1. Center Node (Symmetry)
            newTemps[0] = currentTemps[0] + 2 * (n + 1) * Fo * (currentTemps[1] - currentTemps[0]);

            // 2. Interior Nodes
            for (let i = 1; i < NODES - 1; i++) {
                if (geometry === 'plane') {
                    newTemps[i] = Fo * (currentTemps[i + 1] + currentTemps[i - 1]) + (1 - 2 * Fo) * currentTemps[i];
                } else {
                    const r = i * DR_REAL;
                    const termMinus = 1 - (n * DR_REAL) / (2 * r);
                    const termPlus = 1 + (n * DR_REAL) / (2 * r);
                    newTemps[i] = Fo * (termMinus * currentTemps[i - 1] + termPlus * currentTemps[i + 1]) + (1 - 2 * Fo) * currentTemps[i];
                }
            }

            // 3. Surface Node (Convection)
            const T_inner = currentTemps[NODES - 2];
            const T_surf = currentTemps[NODES - 1];
            newTemps[NODES - 1] = (2 * Fo * T_inner + 2 * Fo * Bi * tempInf + (1 - 2 * Fo - 2 * Fo * Bi) * T_surf);
        }

        if (isNaN(newTemps[0]) || !isFinite(newTemps[0])) {
            setIsRunning(false);
            console.warn("Simulation unstable - halted.");
            return;
        }

        simDataRef.current = newTemps;
        timeRef.current += dt;

        const tolerance = 0.1;
        const isSteady = Math.abs(newTemps[0] - tempInf) < tolerance && Math.abs(newTemps[NODES - 1] - tempInf) < tolerance;

            if (isSteady) {
            // Record the time to steady state the first time we detect it
            if (steadyTimeRef.current === null) steadyTimeRef.current = timeRef.current;

            simDataRef.current = Array(NODES).fill(tempInit);
            timeRef.current = 0;
        }

        setSimState({
            temps: newTemps,
            time: timeRef.current,
            isSteady,
            steadyTime: steadyTimeRef.current
        });

        requestRef.current = requestAnimationFrame(updateSimulation);
    }, [isRunning, alpha, dt, h, k, tempInf, tempInit, geometry]);

    useEffect(() => {
        if (isRunning) requestRef.current = requestAnimationFrame(updateSimulation);
        return () => { if (requestRef.current) cancelAnimationFrame(requestRef.current); };
    }, [updateSimulation, isRunning]);

    useEffect(() => {
        simDataRef.current = Array(NODES).fill(tempInit);
        timeRef.current = 0;
        steadyTimeRef.current = null;
        setSimState({ temps: Array(NODES).fill(tempInit), time: 0, isSteady: false, steadyTime: null });
    }, [tempInit, geometry]);

    // --- ECharts Config ---
    const generateLineChartOption = (): EChartsOption => {
        const xMax = (NODES - 1) * DR_REAL;
        const isDark = resolvedTheme === 'dark';
        const textColor = isDark ? '#eee' : '#333';

        return {
            backgroundColor: 'transparent',
            animation: false,
            title: { text: `Temperature Profile`, left: 'center', textStyle: { color: textColor, fontSize: 14, fontFamily: "'Merriweather Sans', sans-serif" } },
            tooltip: { 
                trigger: 'axis', triggerOn: 'mousemove', showDelay: 0, hideDelay: 300, enterable: true, textStyle: { fontFamily: "'Merriweather Sans', sans-serif" },
                // Use a permissive any type to avoid strict CallbackDataParams typing issues
                formatter: (params: any) => {
                    const p = Array.isArray(params) ? params[0] : params;
                    const raw = (p && (p.data !== undefined)) ? p.data : (p && p.value !== undefined ? p.value : null);
                    let temp: number | null = null;
                    let axisVal: number | null = null;
                    if (Array.isArray(raw)) {
                        axisVal = Number(raw[0]);
                        temp = Number(raw[1]);
                    } else if (typeof raw === 'number') {
                        temp = raw;
                        axisVal = (p && p.name !== undefined) ? Number(p.name) : null;
                    }
                    const tempStr = temp !== null ? temp.toFixed(1) : '';
                    const axisStr = axisVal !== null ? axisVal.toFixed(2) : '';
                    // Removed p.marker to avoid the colored dot next to the temperature
                    return `${p.seriesName}: ${tempStr} °C<br/>Radius: ${axisStr} m`;
                }
            },
            grid: { top: 40, right: 30, bottom: 40, left: 50 },
            visualMap: {
                type: 'continuous',
                dimension: 1, // Map based on y-value (temperature)
                min: Math.min(tempInit, tempInf),
                max: Math.max(tempInit, tempInf),
                inRange: {
                    color: ['#2563eb', '#f59e0b', '#dc2626']
                },
                show: false
            },
            xAxis: {
                type: 'value', min: 0, max: xMax, interval: 0.05, name: 'Radius, r (m)', nameLocation: 'middle', nameGap: 25,
                axisLine: { show: true, lineStyle: { color: textColor } }, axisTick: { show: true }, axisLabel: { color: textColor, fontFamily: "'Merriweather Sans', sans-serif", formatter: (val: any) => Number(val).toFixed(2) },
                nameTextStyle: { fontFamily: "'Merriweather Sans', sans-serif", color: textColor },
                splitLine: { show: false }
            },
            yAxis: {
                type: 'value', name: 'Temperature, T (°C)',
                min: Math.min(tempInit, tempInf) - 5, 
                max: Math.max(tempInit, tempInf) + 5,
                axisLine: { show: true, lineStyle: { color: textColor } }, axisLabel: { color: textColor, fontFamily: "'Merriweather Sans', sans-serif" },
                nameTextStyle: { fontFamily: "'Merriweather Sans', sans-serif", color: textColor },
                nameLocation: 'middle', nameRotate: 90, nameGap: 30,
                splitLine: { show: false }
            },

            series: [{
                name: 'Temperature', type: 'line', smooth: true, showSymbol: false,
                data: simState.temps.map((t, i) => [(i * DR_REAL), t]),
                lineStyle: { width: 3 },
                areaStyle: { opacity: 0.3 },
                labelLine: { lineStyle: { color: textColor } }
            }]
        };
    };

    // Memoize chart option to reduce object churn (helps tooltip stability)
    const chartOption = useMemo(() => {
        return generateLineChartOption();
    }, [simState.time, tempInit, tempInf, resolvedTheme]);

    // Characteristic Length (Radius or Half-thickness) for Fourier Calc
    const R_char = DR_REAL * NODES; 
    
    // Process Fourier Number (Dimensionless Time)
    // Fo = (alpha * time) / L^2
    const fourierNumber = (alpha * simState.time) / (R_char * R_char);

    return (
        <TooltipProvider>
            <div className="container mx-auto p-4 md:p-8 px-4 md:px-16">
                <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                    {/* LEFT COLUMN */}
                    <div className="lg:col-span-1 space-y-6">
                        <Card>
                            <CardContent className="space-y-6">

                                <div className="flex gap-2">
                                    <Button variant={isRunning ? "secondary" : "default"} className="flex-1" onClick={() => setIsRunning(!isRunning)}>
                                        {isRunning ? <><Pause className="mr-2 h-4 w-4" /> Pause</> : <><Play className="mr-2 h-4 w-4" /> Resume</>}
                                    </Button>
                                    <Button variant="outline" size="icon" onClick={() => {
                                        simDataRef.current = Array(NODES).fill(tempInit);
                                        timeRef.current = 0;
                                        steadyTimeRef.current = null;
                                        setSimState({ temps: Array(NODES).fill(tempInit), time: 0, isSteady: false, steadyTime: null });
                                        setIsRunning(true);
                                    }}>
                                        <RotateCcw className="h-4 w-4" />
                                    </Button>
                                </div>

                                <div className="space-y-2">
                                    <Label>Geometry</Label>
                                    <Select value={geometry} onValueChange={(val: GeometryType) => setGeometry(val)}>
                                        <SelectTrigger className="w-full"><SelectValue /></SelectTrigger>
                                        <SelectContent>
                                            <SelectItem value="plane">Plane Wall</SelectItem>
                                            <SelectItem value="cylinder">Cylinder</SelectItem>
                                            <SelectItem value="sphere">Sphere</SelectItem>
                                            <SelectItem value="fin">Pin Fin</SelectItem>
                                            <SelectItem value="hollow_cylinder">Hollow Cylinder</SelectItem>
                                        </SelectContent>
                                    </Select>
                                </div>
                                <div className="space-y-4">
                                    <div className="space-y-2">
                                        <div className="flex justify-between text-sm"><Label>Thermal Conductivity (k)</Label><span>{k} W/(m·K)</span></div>
                                        <Slider min={1} max={200} value={[k]} onValueChange={(v) => setK(v[0])} />
                                    </div>
                                    <div className="space-y-2">
                                        <div className="flex justify-between text-sm"><Label>Heat Transfer Coefficient (h)</Label><span>{h} W/(m²·K)</span></div>
                                        <Slider min={5} max={500} value={[h]} onValueChange={(v) => setH(v[0])} />
                                    </div>
                                    <div className="space-y-2">
                                        <div className="flex justify-between text-sm"><Label>Diffusivity (α scale)</Label><span>{alphaScaled}</span></div>
                                        <Slider min={1} max={50} value={[alphaScaled]} onValueChange={(v) => setAlphaScaled(v[0])} />
                                    </div>
                                </div>
                                <div className="grid grid-cols-2 gap-4">
                                    <div className="space-y-1">
                                        <Label><span className="tight-sub">T<sub>init</sub></span></Label>
                                        <div className="relative">
                                            <Input
                                                type="number"
                                                value={tempInitInput}
                                                onChange={(e) => {
                                                    const v = e.target.value;
                                                    setTempInitInput(v);
                                                    if (v === '') return; // preserve last numeric value when input cleared
                                                    const n = Number(v);
                                                    if (!Number.isNaN(n)) setTempInit(n);
                                                }}
                                                onBlur={() => { if (tempInitInput === '') setTempInitInput(tempInit.toString()); }}
                                                className={`pr-12 min-w-[88px] overflow-hidden whitespace-nowrap ${tempInitInput && tempInitInput.length > 7 ? 'text-transparent caret-transparent' : 'text-center font-mono'}`}
                                            />
                                            <span className="absolute right-3 inset-y-0 flex items-center text-xs text-muted-foreground pointer-events-none">°C</span>
                                        </div>
                                    </div>

                                    <div className="space-y-1">
                                        <Label><span className="tight-sub">T<sub>fluid</sub></span></Label>
                                        <div className="relative">
                                            <Input
                                                type="number"
                                                value={tempInfInput}
                                                onChange={(e) => {
                                                    const v = e.target.value;
                                                    setTempInfInput(v);
                                                    if (v === '') return; // preserve last numeric value when input cleared
                                                    const n = Number(v);
                                                    if (!Number.isNaN(n)) setTempInf(n);
                                                }}
                                                onBlur={() => { if (tempInfInput === '') setTempInfInput(tempInf.toString()); }}
                                                className={`pr-12 min-w-[88px] overflow-hidden whitespace-nowrap ${tempInfInput && tempInfInput.length > 7 ? 'text-transparent caret-transparent' : 'text-center font-mono'}`}
                                            />
                                            <span className="absolute right-3 inset-y-0 flex items-center text-xs text-muted-foreground pointer-events-none">°C</span>
                                        </div>
                                    </div>
                                </div>
                                <div className="p-4 bg-muted/50 rounded-lg text-center space-y-4">
                                    <div className="grid grid-cols-2 gap-4">
                                        <div>
                                            <span className="text-xs uppercase text-muted-foreground font-bold tracking-wider">Biot Number</span>
                                            <div className="text-2xl font-bold">{biotNumber.toFixed(2)}</div>
                                        </div>
                                        <div>
                                            <span className="text-xs uppercase text-muted-foreground font-bold tracking-wider">Fourier No.</span>
                                            <div className="text-2xl font-bold">
                                                {fourierNumber.toFixed(2)}
                                            </div>
                                        </div>
                                    </div>
                                    
                                    <div className="pt-2 border-t border-muted-foreground/20">
                                        <span className="text-xs uppercase text-muted-foreground font-bold tracking-wider">Simulation Time</span>
                                        <div className="font-mono text-lg">{formatDuration(simState.time)}</div>
                                        {simState.isSteady && (simState.steadyTime !== null && simState.steadyTime !== undefined) && (
                                            <div className="text-xs text-muted-foreground mt-1">Steady state reached in {formatDuration(simState.steadyTime)}</div>
                                        )}
                                    </div>
                                </div>
                            </CardContent>
                        </Card>
                    </div>

                    {/* RIGHT COLUMN */}
                    <div className="lg:col-span-2 space-y-6">
                        <Card className="h-[400px] flex flex-col">

                            <CardContent className="flex-1 p-0 relative" style={{ background: 'var(--color-card)' }}>
                                <div className="absolute bottom-2 right-2 z-10 flex items-center gap-3 p-2 rounded">
                                    <span className="tight-sub text-xs">T<sub>init</sub></span>
                                    <div className="w-36 h-2 rounded-full" style={{ background: 'linear-gradient(90deg, #2563eb 0%, #f59e0b 50%, #dc2626 100%)' }} />
                                    <span className="tight-sub text-xs">T<sub>fluid</sub></span>
                                </div>

                                <Canvas camera={{ position: [2, 2, 4], fov: 45 }} style={{ background: 'transparent' }} gl={{ alpha: true, antialias: true }} onCreated={({ gl }) => { gl.setClearColor(new THREE.Color(0, 0, 0), 0); }}>
                                    <ambientLight intensity={0.5} />
                                    <pointLight position={[10, 10, 10]} intensity={1} />
                                    <Center>
                                        <HeatMesh
                                            geometryType={geometry}
                                            temps={simState.temps}
                                            minTemp={tempInit}
                                            maxTemp={tempInf}
                                        />
                                    </Center>
                                    <OrbitControls makeDefault />
                                </Canvas>
                            </CardContent>
                        </Card>
                        <Card className="h-[300px]">
                            <CardContent className="h-full p-2">
                                <ReactECharts option={chartOption} style={{ height: '100%', width: '100%' }} notMerge={false} lazyUpdate={true} />
                            </CardContent>
                        </Card>
                    </div>
                </div>
            </div>
        </TooltipProvider>
    );
}