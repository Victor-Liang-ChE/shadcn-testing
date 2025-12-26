'use client'

import React, { useState, useEffect, useCallback, useRef } from 'react'
import { useTheme } from "next-themes"

// ECharts imports
import ReactECharts from 'echarts-for-react'
import * as echarts from 'echarts/core'
import type { EChartsOption } from 'echarts'
import { LineChart, ScatterChart, CustomChart } from 'echarts/charts'
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  VisualMapComponent
} from 'echarts/components'
import { CanvasRenderer } from 'echarts/renderers'

// React Three Fiber imports for 3D visualization
import { Canvas } from "@react-three/fiber"
import { OrbitControls, PerspectiveCamera } from "@react-three/drei"
import * as THREE from "three"

echarts.use([
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  VisualMapComponent,
  LineChart,
  ScatterChart,
  CustomChart,
  CanvasRenderer
])

import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Label } from "@/components/ui/label"
import { Button } from "@/components/ui/button"
import { Slider } from "@/components/ui/slider"
import 'katex/dist/katex.min.css'
import { BlockMath } from 'react-katex'
import { 
  Play, 
  Pause, 
  RotateCcw, 
  Move, 
  ArrowRight, 
  MousePointer2, 
  Activity, 
  List, 
  Zap,
  Info
} from 'lucide-react'

// Type declaration for react-katex
declare module 'react-katex' {
  export function BlockMath(props: { math: string }): React.ReactElement
}

// Types for our simulation
type Particle = {
  id: number
  x: number
  y: number
  z: number
  vx: number
  vy: number
  vz: number
  fx: number
  fy: number
  fz: number
}

// Simulation Constants
const BOX_SIZE = 100
const DT_DEFAULT = 0.002 // Strict: Lower timestep for stability without velocity caps
const BOLTZMANN_K = 1.0 // Normalized units

// RDF Configuration
const RDF_BINS = 100
const RDF_CUTOFF = BOX_SIZE / 2.0
const RDF_BIN_WIDTH = RDF_CUTOFF / RDF_BINS

export default function MolecularDynamicsPage() {
  const { resolvedTheme } = useTheme()
    const isDarkTheme = resolvedTheme === 'dark'
    const velColor = isDarkTheme ? '#4ade80' : '#16a34a'
    const forceColor = isDarkTheme ? '#facc15' : '#d97706'
    const [mounted, setMounted] = useState(false)
    useEffect(() => { setMounted(true) }, [])
  
  // --- Simulation State ---
  const [running, setRunning] = useState(false)
  const [particles, setParticles] = useState<Particle[]>([])
  const [selectedId, setSelectedId] = useState<number | null>(null)
  const [isThreeD, setIsThreeD] = useState(false)
  
  // --- Control States ---
  const [numParticles, setNumParticles] = useState<number>(50)
  const [temperature, setTemperature] = useState<number>(1.0)
  const [timeStep, setTimeStep] = useState<number>(DT_DEFAULT)
  
  // New L-J Parameters
  const [epsilon, setEpsilon] = useState<number>(50.0)
  const [sigma, setSigma] = useState<number>(5.0)
  
  // Potential Type and Generic Parameter C
    type PotentialType = 'LJ' | 'WCA' | 'MORSE' | 'SOFT' | 'YUKAWA' | 'GAUSS' | 'BUCK'
  const [potentialType, setPotentialType] = useState<PotentialType>('LJ')
  const [paramC, setParamC] = useState<number>(1.5) // Width for Morse, Exponent for Soft

  // --- Visualization States ---
  const [showForce, setShowForce] = useState(false)
  const [showVelocity, setShowVelocity] = useState(false)
  
  // --- Analysis State ---
  const [kineticEnergy, setKineticEnergy] = useState<number>(0)
  const [thermostatZeta, setThermostatZeta] = useState<number>(0) // Visualizing the friction
  const [pressure, setPressure] = useState<number>(0)
  const [rdfData, setRdfData] = useState<{r: number, g: number}[]>([])
  const [energyHistory, setEnergyHistory] = useState<{time: number, ke: number, pe: number, tot: number}[]>([])
    const energyHistoryRef = useRef<{time: number, ke: number, pe: number, tot: number}[]>([])
  const [analysisView, setAnalysisView] = useState<'rdf' | 'energy' | 'msd'>('rdf')
  const [energyMode, setEnergyMode] = useState<'sliding' | 'expanding'>('sliding')
  
  // Add MSD state and refs
  const [msdData, setMsdData] = useState<{time: number, msd: number}[]>([])
  const unwrappedRef = useRef<{x: number, y: number, z: number, x0: number, y0: number, z0: number}[]>([])
  
  // Structure Metrics State
  const [structureMetrics, setStructureMetrics] = useState({
      peakDist: 0,
      peakHeight: 0,
      coordination: 0,
      phaseEstimate: 'Unknown'
  })
  const msdHistoryRef = useRef<{time: number, msd: number}[]>([])

    // Keep tooltips visible while data updates
    const analysisChartRef = useRef<any>(null)
    const analysisHoverRef = useRef(false)
    const analysisPointerRef = useRef<{ x: number; y: number } | null>(null)
    const analysisTooltipTimerRef = useRef<number | null>(null)
  
  // --- Selected Particle View State ---
  const particleHistoryRef = useRef<{time: number, vx: number, vy: number, f: number}[]>([])
  const [historyTrigger, setHistoryTrigger] = useState(0) 

  // Refs for Simulation Loop
  const requestRef = useRef<number | null>(null)
  const particlesRef = useRef<Particle[]>([])
  const isThreeDRef = useRef(isThreeD)
  const initial3DViewRef = useRef(true)
  
  // Refs for current parameter values
  const paramsRef = useRef({
      dt: DT_DEFAULT,
      epsilon: 50.0,
      sigma: 5.0,
      paramC: 1.5,
      potentialType: 'LJ' as PotentialType,
      box: BOX_SIZE,
      targetTemp: 1.0
  })
  
  // Nosé-Hoover State Variable
  // Zeta is the "friction" coefficient that evolves over time
  const nhStateRef = useRef({
      zeta: 0.0,
      Q: 10.0 // Thermal mass (N * k * T * tau^2)
  })
  
  const timeRef = useRef(0)
  
  // RDF Accumulation State
  const rdfHistogramRef = useRef<number[]>(new Array(RDF_BINS).fill(0))
  const rdfFramesRef = useRef<number>(0)

  // Update refs when state changes
  useEffect(() => {
      paramsRef.current.dt = timeStep
      paramsRef.current.epsilon = epsilon
      paramsRef.current.sigma = sigma
      paramsRef.current.paramC = paramC
      paramsRef.current.potentialType = potentialType
      paramsRef.current.targetTemp = temperature
      
      // Update NH Thermal Mass Q roughly based on system size
      // Q approx N * T * tau^2. Let tau = 0.5ps
      nhStateRef.current.Q = numParticles * temperature * 0.25
  }, [timeStep, epsilon, sigma, paramC, potentialType, temperature, numParticles])

  // --- PHYSICS KERNEL: Force Calculation (3D O(N^2)) ---
  const calculateForces = (currentParticles: Particle[]): { pe: number, virial: number } => {
      const { epsilon: eps, sigma: sig, box } = paramsRef.current
      
      let potentialEnergy = 0
      let virialSum = 0
      
      // Reset Forces
      for (let p of currentParticles) {
          p.fx = 0; p.fy = 0; p.fz = 0
      }

      const N = currentParticles.length
      const cutoffSq = 6.25 * sig * sig
      
      // O(N^2) for simplicity in 3D transition
      for (let i = 0; i < N; i++) {
          const p1 = currentParticles[i]
          for (let j = i + 1; j < N; j++) {
              const p2 = currentParticles[j]

              let dx = p2.x - p1.x
              let dy = p2.y - p1.y
              let dz = p2.z - p1.z

              // Minimum Image Convention (3D)
              dx -= box * Math.round(dx / box)
              dy -= box * Math.round(dy / box)
              dz -= box * Math.round(dz / box)

              const distSq = dx*dx + dy*dy + dz*dz
              
              if (distSq > cutoffSq) continue

              const r = Math.sqrt(distSq)
              const r_inv = 1.0 / r
              const r2_inv = r_inv * r_inv
              
              let forceScalar = 0 
              let pairPotential = 0

              switch (paramsRef.current.potentialType) {
                  case 'LJ': {
                      const s2 = sig * sig
                      const s2_r2 = s2 * r2_inv
                      const s6_r6 = s2_r2 * s2_r2 * s2_r2
                      const s12_r12 = s6_r6 * s6_r6
                      
                      forceScalar = (24 * eps * r2_inv) * (2 * s12_r12 - s6_r6)
                      pairPotential = 4 * eps * (s12_r12 - s6_r6)
                      break
                  }
                  case 'WCA': {
                      const cutoffSq = 1.25992 * sig * sig 
                      if (distSq < cutoffSq) {
                          const s2 = sig * sig
                          const s2_r2 = s2 * r2_inv
                          const s6_r6 = s2_r2 * s2_r2 * s2_r2
                          const s12_r12 = s6_r6 * s6_r6
                          
                          forceScalar = (24 * eps * r2_inv) * (2 * s12_r12 - s6_r6)
                          pairPotential = 4 * eps * (s12_r12 - s6_r6) + eps
                      } else {
                          forceScalar = 0
                          pairPotential = 0
                      }
                      break
                  }
                  case 'MORSE': {
                      const width = paramsRef.current.paramC
                      const displacement = r - sig
                      const expTerm = Math.exp(-width * displacement)
                      const term = 1 - expTerm
                      
                      pairPotential = eps * term * term
                      
                      const f_r = 2 * eps * width * term * expTerm
                      forceScalar = f_r * r_inv
                      break
                  }
                  case 'YUKAWA': {
                      const screened = Math.exp(-r / sig)
                      pairPotential = (eps * screened) * r_inv
                      const f_r_yuk = pairPotential * (r_inv + (1.0 / sig))
                      forceScalar = f_r_yuk * r_inv
                      break
                  }
                  case 'GAUSS': {
                      const r2_s2 = distSq / (sig * sig)
                      const gauss = Math.exp(-0.5 * r2_s2)
                      pairPotential = eps * gauss
                      forceScalar = (eps * gauss) / (sig * sig)
                      break
                  }
                  case 'BUCK': {
                      const r_safe = Math.max(r, 0.5 * sig)
                      const A = eps
                      const rho = sig
                      const C = paramsRef.current.paramC
                      const expTerm = Math.exp(-r_safe / rho)
                      const r2 = r_safe * r_safe
                      const r6 = r2 * r2 * r2
                      const r7 = r6 * r_safe
                      pairPotential = A * expTerm - (C / r6)
                      const f_r_buck = (A / rho) * expTerm - (6 * C / r7)
                      forceScalar = f_r_buck / r_safe
                      break
                  }
                  case 'SOFT': {
                      const n = Math.max(1, paramsRef.current.paramC)
                      const ratio = sig * r_inv
                      const ratio_n = Math.pow(ratio, n)
                      pairPotential = eps * ratio_n
                      const f_r = n * eps * ratio_n * r_inv
                      forceScalar = f_r * r_inv
                      break
                  }
              }

              // Accumulate Energy & Virial
              potentialEnergy += pairPotential
              virialSum += forceScalar * distSq 

              // Force Components
              const fx = dx * forceScalar
              const fy = dy * forceScalar
              const fz = dz * forceScalar

              p1.fx -= fx; p1.fy -= fy; p1.fz -= fz
              p2.fx += fx; p2.fy += fy; p2.fz += fz
          }
      }
      return { pe: potentialEnergy, virial: virialSum }
  }

  // --- ALGORITHM 1: Steepest Descent Minimization ---
  // Strictly moves particles "downhill" to fix overlaps
  const runMinimization = (pList: Particle[]) => {
      const MAX_STEPS = 200
      const GAMMA = 0.01 
      const box = paramsRef.current.box

      for (let step = 0; step < MAX_STEPS; step++) {
          const { pe: _ } = calculateForces(pList)
          
          let maxForce = 0
          for (let p of pList) {
              // 3D Magnitude
              const fMag = Math.sqrt(p.fx**2 + p.fy**2 + p.fz**2)
              if (fMag > maxForce) maxForce = fMag
              
              const scale = Math.min(GAMMA, 0.5 / (fMag + 1e-6))
              
              p.x += p.fx * scale
              p.y += p.fy * scale
              p.z += p.fz * scale
              
              // Enforce PBC 3D
              if (p.x < 0) p.x += box; if (p.x >= box) p.x -= box
              if (p.y < 0) p.y += box; if (p.y >= box) p.y -= box
              if (p.z < 0) p.z += box; if (p.z >= box) p.z -= box
              
              // Zero velocities
              p.vx = 0; p.vy = 0; p.vz = 0
          }
          if (maxForce < 10.0) break 
      }
  }

  // --- ALGORITHM 3: Accumulate Radial Distribution Function ---
  const sampleRDF = (currentParticles: Particle[]) => {
      const { box } = paramsRef.current
      const hist = rdfHistogramRef.current
      const N = currentParticles.length

      // O(N^2) Loop required for long-range RDF (up to L/2)
      // Cell lists are too short-range for full structure analysis
      for (let i = 0; i < N; i++) {
          for (let j = i + 1; j < N; j++) {
              let dx = currentParticles[j].x - currentParticles[i].x
              let dy = currentParticles[j].y - currentParticles[i].y

              // MIC (Minimum Image Convention)
              dx -= box * Math.round(dx / box)
              dy -= box * Math.round(dy / box)

              const r = Math.sqrt(dx*dx + dy*dy)

              if (r < RDF_CUTOFF) {
                  const bin = Math.floor(r / RDF_BIN_WIDTH)
                  if (bin >= 0 && bin < RDF_BINS) {
                      hist[bin] += 2 // Add 2 because we found pair (i,j) AND (j,i)
                  }
              }
          }
      }
      rdfFramesRef.current += 1
  }

  // STEP 2: Normalize to g(r)
  const computeRDFGraph = () => {
      const hist = rdfHistogramRef.current
      const frames = rdfFramesRef.current
      const N = particlesRef.current.length
      const { box } = paramsRef.current
      if (frames === 0 || N === 0) return []

      // Global Density (rho)
      const density = N / (box * box)
      
      return hist.map((count, i) => {
          const r = (i + 1) * RDF_BIN_WIDTH
          // Area of annulus shell in 2D: dA = 2 * pi * r * dr
          const area = 2 * Math.PI * r * RDF_BIN_WIDTH
          
          // Expected count for ideal gas = rho * area
          const idealCount = density * area
          
          // g(r) = Observed / (Ideal * N * Frames)
          // We divide by N because the histogram sums over ALL particles as centers
          const g = count / (idealCount * N * frames)
          
          return { r, g }
      })
  }

  // Analyze RDF to find Peak and Coordination Number
  const analyzeStructure = (data: {r: number, g: number}[]) => {
      if (data.length < 5) return

      // 1. Find First Peak (Location and Height)
      let maxVal = 0
      let maxIdx = 0
      // Only look in the first half of the box to avoid PBC noise
      const searchLimit = data.length / 2 
      
      for (let i = 0; i < searchLimit; i++) {
          if (data[i].g > maxVal) {
              maxVal = data[i].g
              maxIdx = i
          }
      }
      
      const r_peak = data[maxIdx].r

      // 2. Find First Minimum after the peak (The "First Shell" Cutoff)
      // We look for the point where slope turns positive or g drops to ~1
      let minIdx = maxIdx
      for (let i = maxIdx + 1; i < data.length - 1; i++) {
          // Simple local minimum detection
          if (data[i].g < data[i+1].g && data[i].g < 1.5) { 
              minIdx = i
              break
          }
      }
      // Fallback: if no minimum found (e.g. gas), integrate up to 1.5 * peak
      if (minIdx === maxIdx) minIdx = Math.min(data.length - 1, Math.floor(maxIdx * 1.5))

      // 3. Compute Coordination Number (Integration)
      // 3D: Integral of 4 * pi * r^2 * rho * g(r) dr
      // 2D: Integral of 2 * pi * r * rho * g(r) dr
      const { box } = paramsRef.current
      const N = numParticles
      const volume = isThreeD ? Math.pow(box, 3) : Math.pow(box, 2)
      const rho = N / volume
      
      let integral = 0
      for (let i = 0; i < minIdx; i++) {
          const r = data[i].r
          const g = data[i].g
          const dr = RDF_BIN_WIDTH
          
          // Geometry factor depends on 2D vs 3D
          const shellVolume = isThreeD 
              ? 4 * Math.PI * r * r * dr 
              : 2 * Math.PI * r * dr
              
          integral += shellVolume * rho * g
      }

      // 4. Estimate Phase (Updated Logic)
      let phase = 'Fluid'
      
      // Heuristic: Check the peak height. 
      // Liquids usually have a first peak height < 2.5. 
      // Solids usually have sharp peaks > 3.0.
      const isStructurallyOrdered = maxVal > 2.8 
      
      if (isThreeD) {
          if (integral > 10.0) phase = 'Solid (FCC/HCP)'     // Perfect Crystal
          else if (integral > 7.0) phase = 'Solid (BCC)'     // Less dense Crystal
          else if (isStructurallyOrdered) phase = 'Solid (Amorphous)' // Sharp peak, low coord (Glass)
          else phase = 'Liquid/Gas'
      } else {
          // 2D Logic
          if (integral > 5.0) phase = 'Solid (Hexagonal)'
          else if (integral > 3.5) phase = 'Solid (Square)'
          else if (isStructurallyOrdered) phase = 'Solid (Glassy)'
          else phase = 'Liquid'
      }
      
      // Override: If density is very low, it's likely a gas regardless of local structure
      const globalRho = numParticles / (isThreeD ? Math.pow(box, 3) : Math.pow(box, 2))
      if (globalRho < 0.05) phase = 'Gas'

      setStructureMetrics({
          peakDist: r_peak,
          peakHeight: maxVal,
          coordination: integral,
          phaseEstimate: phase
      })
  }

  // Detect System State (Liquid vs Crystalline vs Amorphous Solid)
  const detectSystemState = () => {
      // 1. DYNAMICS CHECK (Is it Solid or Liquid?)
      // Calculate Slope of MSD (Diffusion Coefficient)
      let isLiquid = false
      let diffusionCoeff = 0
      
      if (msdData.length > 20) {
          // Look at last 50% of data to ignore initial ballistic regime
          const window = msdData.slice(Math.floor(msdData.length * 0.5))
          const dt = window[window.length - 1].time - window[0].time
          const dMSD = window[window.length - 1].msd - window[0].msd
          
          // Slope = MSD / t. 
          // In 3D, D = Slope / 6. In 2D, D = Slope / 4. 
          // We just check the raw slope for a threshold.
          const slope = dMSD / dt
          diffusionCoeff = slope
          
          // Threshold: If slope > 0.1 sigma^2/ps, it's diffusing significantly
          isLiquid = slope > 0.1
      }

      // 2. STRUCTURE CHECK (Is it Crystalline or Amorphous?)
      // We use the "Peak Height" from your existing structureMetrics
      // Crystalline solids usually have g(r) peaks > 3.0 or 4.0
      // Amorphous solids/Liquids usually have g(r) peaks < 2.8
      const isOrdered = structureMetrics.peakHeight > 2.8

      // 3. COMBINE
      let detectedState = 'Unknown'
      let color = 'text-gray-500'

      if (isLiquid) {
          detectedState = 'Liquid'
          color = 'text-blue-500'
      } else {
          if (isOrdered) {
              detectedState = 'Crystalline Solid'
              color = 'text-green-500'
          } else {
              detectedState = 'Amorphous Solid (Glass)'
              color = 'text-amber-500'
          }
      }
      
      return { state: detectedState, color, diffusion: diffusionCoeff }
  }

  // --- Initialization Logic ---
  const initializeParticles = useCallback(() => {
    const newParticles: Particle[] = []
    
    // 3D Logic: Cube Root vs Square Root
    // If 125 particles -> 5x5x5 grid
    const dim = isThreeD 
      ? Math.ceil(Math.pow(numParticles, 1/3)) 
      : Math.ceil(Math.sqrt(numParticles))
      
    const spacing = BOX_SIZE / dim

    // 1. Grid Placement (3D or 2D)
    for (let i = 0; i < numParticles; i++) {
      let x, y, z

      if (isThreeD) {
          // 3D Grid Indices
          const ix = i % dim
          const iy = Math.floor(i / dim) % dim
          const iz = Math.floor(i / (dim * dim))
          
          x = ix * spacing + spacing / 2
          y = iy * spacing + spacing / 2
          z = iz * spacing + spacing / 2
      } else {
          // 2D Grid Indices
          const ix = i % dim
          const iy = Math.floor(i / dim)
          
          x = ix * spacing + spacing / 2
          y = iy * spacing + spacing / 2
          z = 0
      }

      // Add jitter to avoid perfect lattice (helps minimizer)
      x += (Math.random() - 0.5) * (spacing * 0.5)
      y += (Math.random() - 0.5) * (spacing * 0.5)
      z += isThreeD ? (Math.random() - 0.5) * (spacing * 0.5) : 0

      // Safety Clamp to prevent initial wrapping
      if (x < 0) x += BOX_SIZE; if (x >= BOX_SIZE) x -= BOX_SIZE;
      if (y < 0) y += BOX_SIZE; if (y >= BOX_SIZE) y -= BOX_SIZE;
      if (isThreeD) {
          if (z < 0) z += BOX_SIZE; if (z >= BOX_SIZE) z -= BOX_SIZE;
      } else {
          z = 0; // Force flat 0 for 2D
      }

      newParticles.push({
        id: i,
        x, y, z,
        vx: (Math.random() - 0.5) * temperature,
        vy: (Math.random() - 0.5) * temperature,
        vz: isThreeD ? (Math.random() - 0.5) * temperature : 0,
        fx: 0, fy: 0, fz: 0
      })
    }
    // 3. Assign Maxwell-Boltzmann Velocities
    for (let p of newParticles) {
        p.vx = (Math.random() - 0.5) * temperature
        p.vy = (Math.random() - 0.5) * temperature
        p.vz = isThreeD ? (Math.random() - 0.5) * temperature : 0
    }

    // 4. Remove Center of Mass Motion (Optional but good practice)
    let vcmX = 0, vcmY = 0, vcmZ = 0
    for (let p of newParticles) { vcmX += p.vx; vcmY += p.vy; vcmZ += p.vz }
    vcmX /= numParticles; vcmY /= numParticles; vcmZ /= numParticles
    for (let p of newParticles) { p.vx -= vcmX; p.vy -= vcmY; p.vz -= vcmZ }

    // Seed energy history at t=0 so the chart doesn't start with an empty left section.
    const { pe: pe0 } = calculateForces(newParticles) // <--- Destructure just PE
    let ke0 = 0
    for (const p of newParticles) ke0 += 0.5 * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz)
    
    setParticles(newParticles)
    particlesRef.current = newParticles
    setSelectedId(null)
    particleHistoryRef.current = []
    rdfHistogramRef.current = new Array(RDF_BINS).fill(0) // Clear RDF bins
    rdfFramesRef.current = 0 // Reset RDF frame counter
    setRdfData([]) // Clear RDF display
    
    // Initialize Unwrapped Coordinates for MSD
    unwrappedRef.current = newParticles.map(p => ({
        x: p.x, y: p.y, z: p.z,
        x0: p.x, y0: p.y, z0: p.z
    }))
    
    // Reset MSD History
    msdHistoryRef.current = [{ time: 0, msd: 0 }]
    setMsdData(msdHistoryRef.current)
        energyHistoryRef.current = [{ time: 0, ke: ke0, pe: pe0, tot: ke0 + pe0 }]
        setEnergyHistory(energyHistoryRef.current)
        setKineticEnergy(ke0)
        setThermostatZeta(0)
        setPressure(0)
        setHistoryTrigger(0)
    timeRef.current = 0
    nhStateRef.current.zeta = 0 // Reset Thermostat
  }, [numParticles, temperature, isThreeD])

  useEffect(() => {
    initializeParticles()
  }, [])

  // --- ALGORITHM 2: Velocity Verlet with Nosé-Hoover Thermostat ---
  const updatePhysics = () => {
    const currentParticles = [...particlesRef.current]
    const { dt, box, targetTemp } = paramsRef.current
    const { Q } = nhStateRef.current
    let zeta = nhStateRef.current.zeta
    
    timeRef.current += dt

    const N = currentParticles.length
    // Degrees of Freedom: 3N - 3 for 3D, 2N - 2 for 2D
    const dof = isThreeDRef.current ? (3 * N - 3) : (2 * N - 2)
    
    // 1. Calculate Forces at t
    calculateForces(currentParticles)

    // 2. First Half-step Update
    for (let p of currentParticles) {
        const u = unwrappedRef.current[p.id]
        
        // Update Velocity (Half Step) with thermostat
        p.vx += 0.5 * dt * (p.fx - zeta * p.vx)
        p.vy += 0.5 * dt * (p.fy - zeta * p.vy)
        p.vz += 0.5 * dt * (p.fz - zeta * p.vz)

        // Update Position (Full Step)
        p.x += p.vx * dt
        p.y += p.vy * dt
        p.z += p.vz * dt

        // Update UNWRAPPED Position (Full Step) - FOR MSD
        u.x += p.vx * dt
        u.y += p.vy * dt
        u.z += p.vz * dt

        // PBC
        if (p.x < 0) p.x += box
        if (p.x >= box) p.x -= box
        if (p.y < 0) p.y += box
        if (p.y >= box) p.y -= box
        if (p.z < 0) p.z += box
        if (p.z >= box) p.z -= box
    }

    // 3. Calculate Forces at t + dt
    const { pe, virial } = calculateForces(currentParticles)

    // 4. Update Thermostat Variable (Zeta)
    let keSum = 0
    for (let p of currentParticles) {
        keSum += 0.5 * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz)
    }
    const currentTemp = (2.0 * keSum) / dof
    
    zeta += dt * (dof * BOLTZMANN_K * (currentTemp - targetTemp)) / Q

    // 5. Second Half-step Update
    let kineticEnergyFull = 0
    for (let p of currentParticles) {
        const vxh = p.vx
        const vyh = p.vy
        const vzh = p.vz
        
        p.vx = (vxh + 0.5 * dt * p.fx) / (1.0 + 0.5 * dt * zeta)
        p.vy = (vyh + 0.5 * dt * p.fy) / (1.0 + 0.5 * dt * zeta)
        p.vz = (vzh + 0.5 * dt * p.fz) / (1.0 + 0.5 * dt * zeta)
        
        kineticEnergyFull += 0.5 * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz)
    }

    nhStateRef.current.zeta = zeta

    // --- Calculate Pressure ---
    // 2D: Volume = L^2, Divisor = 2
    // 3D: Volume = L^3, Divisor = 3
    const volume = isThreeD ? Math.pow(box, 3) : Math.pow(box, 2)
    const virialDivisor = isThreeD ? 3.0 : 2.0
    
    const finalTemp = (2.0 * kineticEnergyFull) / dof
    const pressureVal = (N * BOLTZMANN_K * finalTemp + (virial / virialDivisor)) / volume

    setPressure(pressureVal)

    // --- Calculate MSD ---
    let sqDispSum = 0
    for (let i = 0; i < N; i++) {
        const u = unwrappedRef.current[i]
        const dx = u.x - u.x0
        const dy = u.y - u.y0
        const dz = u.z - u.z0
        // Use 3D distance if isThreeD, otherwise 2D
        sqDispSum += (dx*dx + dy*dy + (isThreeDRef.current ? dz*dz : 0))
    }
    const currentMSD = sqDispSum / N
    
    // Update MSD History
    msdHistoryRef.current.push({ time: timeRef.current, msd: currentMSD })

    // Update State
    particlesRef.current = currentParticles
    setKineticEnergy(kineticEnergyFull)
    setThermostatZeta(zeta)

    // Track Energy History every physics step (real samples, no synthetic padding)
    const tot = kineticEnergyFull + pe
    energyHistoryRef.current.push({ time: timeRef.current, ke: kineticEnergyFull, pe, tot })
    
    // Sample RDF statistics (every step)
    sampleRDF(currentParticles)

    setFrameCount(prev => {
        const next = prev + 1
        if (next % 4 === 0) { // Render every 4th physics frame
            setParticles(currentParticles)
            
            // Compute and Update RDF Graph Data (every 4th frame)
            const newRdfData = computeRDFGraph()
            setRdfData(newRdfData)
            
            // Analyze structure metrics from RDF
            analyzeStructure(newRdfData)
            
            // Update chart data (batched to UI render frames)
            setEnergyHistory([...energyHistoryRef.current])
            setMsdData([...msdHistoryRef.current])
            
            if (selectedId !== null) setHistoryTrigger(prev => prev + 1)
        }
        return next
    })
  }

  // Helper for re-render
  const [_, setFrameCount] = useState(0)

  // --- Animation Loop ---
  const animate = () => {
    updatePhysics()
    requestRef.current = requestAnimationFrame(animate)
  }

  useEffect(() => {
    if (running) {
      requestRef.current = requestAnimationFrame(animate)
    } else {
      if (requestRef.current) cancelAnimationFrame(requestRef.current)
    }
    return () => {
      if (requestRef.current) cancelAnimationFrame(requestRef.current)
    }
  }, [running])

  // --- Handlers ---
  const handleParticleCountChange = (newCount: number) => {
    // For scientific consistency, changing particle count necessitates a reset/minimization
    // Otherwise inserting random atoms creates massive energy spikes.
    setNumParticles(newCount)
    // We delay the re-init to the effect hook or user can press reset.
  }

  const handleStartStop = () => setRunning(!running)
  const handleReset = () => {
      setRunning(false)
      initializeParticles()
  }
  
  // Reset effect when numParticles changes
  useEffect(() => {
      initializeParticles()
  }, [initializeParticles])

  // Keep isThreeDRef in sync with isThreeD state
  useEffect(() => {
      isThreeDRef.current = isThreeD
  }, [isThreeD])

  // Reset camera view flag when switching to 3D
  useEffect(() => {
      if (isThreeD) {
          initial3DViewRef.current = true
          // Force a re-init to generate Z coordinates
          initializeParticles()
      }
  }, [isThreeD, initializeParticles])

  const refreshAnalysisTooltip = useCallback(() => {
      if (!analysisHoverRef.current) return
      const pt = analysisPointerRef.current
      if (!pt) return
      const inst = analysisChartRef.current?.getEchartsInstance?.()
      if (!inst) return
      inst.dispatchAction({ type: 'showTip', x: pt.x, y: pt.y })
  }, [])

  const stopAnalysisTooltipLoop = useCallback(() => {
      if (analysisTooltipTimerRef.current !== null) {
          window.clearInterval(analysisTooltipTimerRef.current)
          analysisTooltipTimerRef.current = null
      }
  }, [])

  const startAnalysisTooltipLoop = useCallback(() => {
      if (analysisTooltipTimerRef.current !== null) return
      analysisTooltipTimerRef.current = window.setInterval(() => {
          refreshAnalysisTooltip()
      }, 50)
  }, [refreshAnalysisTooltip])

  // After each data refresh, re-show the tooltip at the last cursor position.
  useEffect(() => {
      // Defer to the next frame so it runs after ECharts applies setOption.
      requestAnimationFrame(() => refreshAnalysisTooltip())
  }, [rdfData, energyHistory, analysisView, refreshAnalysisTooltip])

  useEffect(() => {
      return () => stopAnalysisTooltipLoop()
  }, [stopAnalysisTooltipLoop])

  const onChartClick = (params: any) => {
      const dataValues = params.data.value || params.data;
      
      if (Array.isArray(dataValues) && dataValues.length >= 6) {
          const clickedId = dataValues[5];
          setSelectedId(prevId => {
              if (prevId === clickedId) {
                  particleHistoryRef.current = []
                  return null
              }
              particleHistoryRef.current = []
              return clickedId
          });
      }
  }

  // --- Chart Configuration ---
  const getParticleOption = (): EChartsOption => {
    const isDark = resolvedTheme === 'dark'

    // --- 2D MODE CONFIGURATION ---
    const renderArrow = (params: any, api: any) => {
        const x = api.value(0)
        const y = api.value(1)
        const dxVal = api.value(2)
        const dyVal = api.value(3)
        const type = api.value(4)
        const id = api.value(5)
        
        const isSelected = id === selectedId
        
        let isVisible = false
        if (selectedId !== null) {
            isVisible = isSelected
        } else {
            isVisible = (type === 0 && showVelocity) || (type === 1 && showForce)
        }

        if (!isVisible) return

        if (Math.abs(dxVal) < 0.1 && Math.abs(dyVal) < 0.1) return

        const scale = type === 0 ? 0.5 : 0.05
        const start = api.coord([x, y])
        const end = api.coord([x + dxVal * scale, y + dyVal * scale])
        
        const dx = end[0] - start[0]
        const dy = end[1] - start[1]
        const angle = Math.atan2(dy, dx)
        const headLen = 8
        const headAngle = Math.PI / 6

        const head1X = end[0] - headLen * Math.cos(angle - headAngle)
        const head1Y = end[1] - headLen * Math.sin(angle - headAngle)
        const head2X = end[0] - headLen * Math.cos(angle + headAngle)
        const head2Y = end[1] - headLen * Math.sin(angle + headAngle)

        const color = type === 0 
          ? (isDark ? '#4ade80' : '#16a34a') 
          : (isDark ? '#facc15' : '#d97706')

        return {
            type: 'group',
            children: [
                {
                    type: 'line',
                    shape: { x1: start[0], y1: start[1], x2: end[0], y2: end[1] },
                    style: { stroke: color, lineWidth: 2 }
                },
                {
                    type: 'polygon',
                    shape: { points: [[end[0], end[1]], [head1X, head1Y], [head2X, head2Y]] },
                    style: { fill: color }
                }
            ]
        } as any
    }

    const forceData = particles.map(p => [p.x, p.y, p.fx, p.fy, 1, p.id])
    const velData = particles.map(p => [p.x, p.y, p.vx, p.vy, 0, p.id])
    
    const scatterData = particles.map(p => {
        const isSelected = p.id === selectedId
        return {
            value: [p.x, p.y, Math.sqrt(p.vx**2 + p.vy**2), p.id, 0, p.id],
            itemStyle: {
                color: new echarts.graphic.LinearGradient(0, 0, 1, 0, [
                    { offset: 0, color: '#3b82f6' },
                    { offset: 1, color: '#ef4444' } 
                ]),
                borderColor: isSelected ? '#facc15' : 'transparent',
                borderWidth: isSelected ? 3 : 0,
                shadowBlur: isSelected ? 10 : 0,
                shadowColor: '#facc15'
            }
        }
    })

    return {
      backgroundColor: 'transparent',
      grid: { left: 0, right: 0, top: 0, bottom: 0 },
      xAxis: { type: 'value', min: 0, max: BOX_SIZE, show: false },
      yAxis: { type: 'value', min: 0, max: BOX_SIZE, show: false },
      series: [
        {
            type: 'custom',
            renderItem: renderArrow,
            data: forceData,
            z: 4,
            silent: true,
            animation: false
        },
        {
            type: 'custom',
            renderItem: renderArrow,
            data: velData,
            z: 4,
            silent: true,
            animation: false
        },
        {
          type: 'scatter',
          symbolSize: (data: any) => {
             return Math.max(5, Math.min(30, sigma * 2.5))
          },
          data: scatterData,
          emphasis: { scale: false, focus: 'none' },
          z: 10,
          animation: false
        } as any
      ] as any
    }
  }

  const selectedParticleData = particles.find(p => p.id === selectedId)

  // --- RDF Visualization ---
  const getRdfOption = (): EChartsOption => {
      const isDark = resolvedTheme === 'dark'
      const textColor = isDark ? '#ffffff' : '#000000'
      
      return {
          backgroundColor: 'transparent',
          title: {
              text: 'Radial Distribution Function g(r)',
              left: 'center',
              textStyle: { fontSize: 12, color: textColor, fontFamily: 'Merriweather Sans' }
          },
          tooltip: { 
              trigger: 'axis',
              triggerOn: 'mousemove|click',
              alwaysShowContent: false,
              hideDelay: 100,
              confine: true,
              enterable: true,
              backgroundColor: isDark ? 'rgba(31, 41, 55, 0.95)' : 'rgba(255, 255, 255, 0.95)',
              borderColor: isDark ? '#374151' : '#e5e7eb',
              textStyle: { fontFamily: 'Merriweather Sans', color: textColor },
              axisPointer: { type: 'line' },
              formatter: (params: any) => {
                  if (!Array.isArray(params)) return ''
                  const p = params[0]
                  if (!p || p.seriesName === 'Ideal Gas') return ''
                  // Remove trailing zeros from distance
                  const r = parseFloat(p.value[0]).toFixed(3).replace(/\.?0+$/, '')
                  const g = parseFloat(p.value[1]).toPrecision(3)
                  return `r: ${r} σ<br/>g(r): ${g}`
              }
          },
          grid: { left: 40, right: 10, top: 30, bottom: 40 },
          xAxis: {
              type: 'value',
              name: 'r (Distance) [σ]',
              nameLocation: 'middle',
              nameGap: 30,
              nameTextStyle: { color: textColor, fontSize: 10, fontFamily: 'Merriweather Sans' },
              axisLabel: { color: textColor, fontSize: 10, fontFamily: 'Merriweather Sans' },
              axisLine: { show: true, lineStyle: { color: textColor } },
              splitLine: { show: false }
          },
          yAxis: {
              type: 'value',
              name: 'g(r)',
              min: 0,
              nameTextStyle: { color: textColor, fontSize: 10, fontFamily: 'Merriweather Sans' },
              axisLabel: { color: textColor, fontSize: 10, fontFamily: 'Merriweather Sans', formatter: (v: number) => v.toPrecision(3) },
              axisLine: { show: true, lineStyle: { color: textColor } },
              splitLine: { show: false }
          },
          series: [
              {
                  type: 'line',
                  showSymbol: false,
                  data: rdfData.map(d => [d.r, d.g]),
                  lineStyle: { color: '#8b5cf6', width: 2 },
                  areaStyle: { opacity: 0.1, color: '#8b5cf6' },
                  markLine: rdfData.length > 0 ? {
                      symbol: 'none',
                      data: [
                          { 
                              xAxis: structureMetrics.peakDist, 
                              lineStyle: { type: 'dashed', color: '#f59e0b' },
                              label: { formatter: 'Peak', position: 'end', color: textColor, fontFamily: 'Merriweather Sans', fontSize: 10 }
                          }
                      ]
                  } : undefined
              },
              {
                  type: 'line',
                  name: 'Ideal Gas',
                  data: [[0,1], [RDF_CUTOFF, 1]],
                  symbol: 'none',
                  lineStyle: { type: 'dashed', color: isDark ? '#555' : '#999', width: 1 }
              }
          ],
          animation: false
      }
  }

  // --- Energy Chart Configuration ---
  const getEnergyOption = (): EChartsOption => {
      const isDark = resolvedTheme === 'dark'
      const textColor = isDark ? '#ffffff' : '#000000'
        // Choose window mode based on energyMode state
        let minTime: number, maxTime: number
        if (energyMode === 'sliding') {
            const windowSize = 0.5 // ps
            maxTime = Math.max(timeRef.current, windowSize)
            minTime = Math.max(0, maxTime - windowSize)
        } else {
            // expanding mode: min at 0, max at current time
            minTime = 0
            maxTime = timeRef.current
        }

        // No filtering/padding here: we always plot real samples only.
        // ECharts will clip to [minTime, maxTime] via the xAxis min/max.
        const totalData = energyHistory.map((h): [number, number] => [h.time, h.tot])
        const peData = energyHistory.map((h): [number, number] => [h.time, h.pe])
        const keData = energyHistory.map((h): [number, number] => [h.time, h.ke])

                        return {
                    // Explicit palette ensures legend markers match series line colors
                    color: ['#10b981', '#3b82f6', '#ef4444'],
                    title: { text: 'Energy', left: 'center', textStyle: { fontSize: 12, color: textColor, fontFamily: 'Merriweather Sans' } },
          tooltip: {
              trigger: 'axis',
              triggerOn: 'mousemove|click',
              alwaysShowContent: false,
              hideDelay: 100,
              confine: true,
              enterable: true,
              backgroundColor: isDark ? 'rgba(31, 41, 55, 0.95)' : 'rgba(255, 255, 255, 0.95)',
              borderColor: isDark ? '#374151' : '#e5e7eb',
              textStyle: { fontFamily: 'Merriweather Sans', color: textColor },
              formatter: (params: any) => {
                  if (!Array.isArray(params) || params.length === 0) return ''
                  // Show the time once, then unique series entries (dedupe by seriesName)
                  const first = params[0]
                  const time = first && Array.isArray(first.value) ? Number(first.value[0]) : NaN
                  const timeLine = `Time: ${Number.isFinite(time) ? time.toFixed(3) + ' ps' : '—'}`
                  const seen = new Set<string>()
                  const lines: string[] = []
                  for (const p of params) {
                      const name = String(p.seriesName ?? p.seriesIndex)
                      if (seen.has(name)) continue
                      seen.add(name)
                      const val = Array.isArray(p.value) ? p.value[1] : p.value
                      lines.push(`${name}: ${Number(val).toPrecision(4)} ε`)
                  }
                  return [timeLine, ...lines].join('<br/>')
              }
          },
          legend: { bottom: -5, symbolKeepAspect: true, itemWidth: 8, itemHeight: 8, textStyle: { color: textColor, fontFamily: 'Merriweather Sans' } },
          grid: { left: 70, right: 20, top: 50, bottom: 60, containLabel: false },
          xAxis: {
              type: 'value',
              min: minTime,
              max: maxTime,
              boundaryGap: [0, 0],
              show: true,
              position: 'bottom',
              name: 'Time (ps)',
              nameLocation: 'middle',
              nameGap: 25,
              nameTextStyle: { color: isDark ? '#ffffff' : '#000', fontSize: 10, fontFamily: 'Merriweather Sans' },
              axisLine: { show: true, lineStyle: { color: isDark ? '#ffffff' : '#000' }, onZero: false },
              axisTick: { show: true },
              splitLine: { show: false },
              axisLabel: { showMinLabel: false, showMaxLabel: false, color: isDark ? '#ffffff' : '#000', fontFamily: 'Merriweather Sans', margin: 5, formatter: (v: number) => v.toFixed(1) }
          },
          yAxis: { type: 'value', name: 'Energy (ε)', nameLocation: 'middle', nameGap: 55, nameTextStyle: { color: textColor, fontSize: 10, fontFamily: 'Merriweather Sans' }, splitLine: { show: false }, axisLine: { show: true, lineStyle: { color: textColor }, onZero: false }, axisTick: { show: true }, axisLabel: { color: textColor, fontFamily: 'Merriweather Sans', margin: 5 } },
          series: [
              { name: 'Total E', type: 'line', data: totalData, showSymbol: false, lineStyle: { width: 2, color: '#10b981' } },
              { name: 'Potential (PE)', type: 'line', data: peData, showSymbol: false, lineStyle: { width: 1.5, color: '#3b82f6' } },
              { name: 'Kinetic (KE)', type: 'line', data: keData, showSymbol: false, lineStyle: { width: 1.5, color: '#ef4444' } }
          ],
          animation: false
      }
  }

  const getMsdOption = (): EChartsOption => {
      const isDark = resolvedTheme === 'dark'
      const textColor = isDark ? '#ffffff' : '#000000'

      // Expanding window: min stays at 0, max continuously grows
      const minTime = 0
      const maxTime = timeRef.current
      
      // --- PLATEAU DETECTION LOGIC ---
      let markLineOpt: any = undefined
      
      // Only check if we have enough data points (e.g., > 50 frames) and time has passed
      if (msdData.length > 50 && timeRef.current > 2.0) {
          // Look at the last 20% of the data
          const windowSize = Math.floor(msdData.length * 0.2)
          const recentData = msdData.slice(-windowSize)
          
          // Calculate Average MSD in this window
          const sum = recentData.reduce((acc, cur) => acc + cur.msd, 0)
          const avg = sum / recentData.length
          
          // Calculate Slope (End - Start) / Time
          const start = recentData[0]
          const end = recentData[recentData.length - 1]
          const slope = (end.msd - start.msd) / (end.time - start.time)
          
          // Heuristic: If slope is very small (< 0.1), we consider it a "Plateau" (Solid)
          // If slope is high, it's diffusing (Liquid), so we don't show a plateau line.
          if (Math.abs(slope) < 0.1) {
              markLineOpt = {
                  symbol: ['none', 'none'],
                  label: {
                      position: 'insideEndTop',
                      formatter: '{b}: {c}',
                      color: textColor,
                      fontFamily: 'Merriweather Sans',
                      fontSize: 11
                  },
                  data: [
                      {
                          yAxis: parseFloat(avg.toPrecision(3)),
                          name: 'Vibration Limit (Solid)',
                          lineStyle: { type: 'dashed', color: '#ef4444', width: 2 },
                          label: {
                              show: true,
                              position: 'insideEndTop',
                              formatter: (params: any) => {
                                  const val = Number(params.value)
                                  return `${params.name}: ${val.toPrecision(3)}`
                              },
                              color: textColor,
                              fontFamily: 'Merriweather Sans',
                              fontSize: 11
                          }
                      }
                  ]
              }
          }
      }

      return {
          title: { 
              text: 'Mean Squared Displacement', 
              left: 'center', 
              textStyle: { fontSize: 12, color: textColor, fontFamily: 'Merriweather Sans' } 
          },
          tooltip: {
              trigger: 'axis',
              triggerOn: 'mousemove|click',
              alwaysShowContent: false,
              hideDelay: 100,
              confine: true,
              enterable: true,
              backgroundColor: isDark ? 'rgba(31, 41, 55, 0.95)' : 'rgba(255, 255, 255, 0.95)',
              borderColor: isDark ? '#374151' : '#e5e7eb',
              textStyle: { fontFamily: 'Merriweather Sans', color: textColor },
              formatter: (params: any) => {
                  if (!params[0]) return ''
                  const t = params[0].value[0].toFixed(2)
                  const msd = params[0].value[1].toFixed(3)
                  return `Time: ${t} ps<br/>MSD: ${msd} σ²`
              }
          },
          grid: { left: 50, right: 20, top: 40, bottom: 40 },
          xAxis: {
              type: 'value',
              min: minTime,
              max: maxTime,
              boundaryGap: [0, 0],
              show: true,
              position: 'bottom',
              name: 'Time (ps)',
              nameLocation: 'middle',
              nameGap: 25,
              nameTextStyle: { color: textColor, fontSize: 10, fontFamily: 'Merriweather Sans' },
              axisLine: { show: true, lineStyle: { color: textColor }, onZero: false },
              axisTick: { show: true },
              splitLine: { show: false },
              axisLabel: { showMinLabel: false, showMaxLabel: false, color: textColor, fontFamily: 'Merriweather Sans', margin: 5, formatter: (v: number) => v.toFixed(1) }
          },
          yAxis: {
              type: 'value',
              name: 'MSD (σ²)',
              nameLocation: 'middle',
              nameGap: 35,
              nameTextStyle: { color: textColor, fontSize: 10, fontFamily: 'Merriweather Sans' },
              axisLine: { show: true, lineStyle: { color: textColor }, onZero: false },
              axisTick: { show: true },
              splitLine: { show: false },
              axisLabel: { color: textColor, fontFamily: 'Merriweather Sans', margin: 5 }
          },
          series: [{
              type: 'line',
              showSymbol: false,
              data: msdData.map(d => [d.time, d.msd]),
              lineStyle: { width: 2, color: '#f59e0b' },
              areaStyle: { opacity: 0.1, color: '#f59e0b' },
              markLine: markLineOpt
          }],
          animation: false
      }
  }

  // Helper to get labels based on current potential type
  const getLabels = () => {
      switch (potentialType) {
          case 'MORSE': return {
              eps: 'Well Depth (D)',
              sig: <span className="font-mono">Equilibrium (r<sub className="font-mono">e</sub>)</span>
          }
          case 'SOFT': return {
              eps: 'Strength (ε)',
              sig: 'Scale (σ)'
          }
          case 'WCA': return {
              eps: 'Repulsion Strength (ε)',
              sig: 'Core Radius (σ)'
          }
          case 'YUKAWA': return {
              eps: 'Coupling (ε)',
              sig: 'Screening Length (λ)'
          }
          case 'BUCK': return {
              eps: 'Repulsion (A)',
              sig: 'Decay Length (ρ)'
          }
          case 'GAUSS': return {
              eps: 'Energy Barrier (ε)',
              sig: 'Gaussian Width (σ)'
          }
          default: return {
              eps: 'L-J Epsilon (ε)',
              sig: 'L-J Sigma (σ)'
          }
      }
  }

  const labels = getLabels()

  const getEquationInfo = (type: string) => {
    switch (type) {
      case 'LJ':
        return (
          <>
            <p className="font-bold mb-1 text-center">Lennard-Jones</p>
            <div className="text-xs">
              <BlockMath math={'V(r) = 4\\epsilon \\left[ \\left(\\frac{\\sigma}{r}\\right)^{12} - \\left(\\frac{\\sigma}{r}\\right)^{6} \\right]'} />
            </div>
            <p className="text-[10px] mt-2 text-muted-foreground text-center">Standard model for noble gases.</p>
          </>
        )
      case 'WCA':
        return (
          <>
            <p className="font-bold mb-1 text-center">WCA Potential</p>
            <div className="text-xs">
               <BlockMath math={'V(r) = \\begin{cases} V_{LJ}(r) + \\epsilon & r < 2^{1/6}\\sigma \\\\ 0 & r \\ge 2^{1/6}\\sigma \\end{cases}'} />
            </div>
            <p className="text-[10px] mt-2 text-muted-foreground text-center">Purely repulsive hard-sphere approximation.</p>
          </>
        )
            case 'MORSE':
                return (
                    <>
                        <p className="font-bold mb-1 text-center">Morse Potential</p>
                        <div className="text-xs">
                            <BlockMath math={'V(r) = D \\left[ 1 - e^{-a(r-r_e)} \\right]^2'} />
                        </div>
                        <p className="text-[10px] mt-2 text-muted-foreground text-center">Models chemical bonds and vibrations.</p>
                    </>
                )
            case 'YUKAWA':
                return (
                    <>
                        <p className="font-bold mb-1 text-center">Yukawa (Screened Coulomb)</p>
                        <div className="text-xs">
                            <BlockMath math={'V(r) = \\epsilon \\frac{e^{-r/\\lambda}}{r}'} />
                        </div>
                        <p className="text-[10px] mt-2 text-muted-foreground text-center">Models charged particles in plasma or salt solution.</p>
                    </>
                )
            case 'GAUSS':
                return (
                    <>
                        <p className="font-bold mb-1 text-center">Gaussian Core</p>
                        <div className="text-xs">
                            <BlockMath math={'V(r) = \\epsilon e^{-r^2 / 2\\sigma^2}'} />
                        </div>
                        <p className="text-[10px] mt-2 text-muted-foreground text-center">Finite repulsion at overlap. Used for polymers.</p>
                    </>
                )
            case 'BUCK':
                return (
                    <>
                        <p className="font-bold mb-1 text-center">Buckingham</p>
                        <div className="text-xs">
                            <BlockMath math={'V(r) = A e^{-r/\\rho} - \\frac{C}{r^6}'} />
                        </div>
                        <p className="text-[10px] mt-2 text-muted-foreground text-center">More accurate repulsion for ionic crystals.</p>
                    </>
                )
            case 'SOFT':
                return (
                    <>
                        <p className="font-bold mb-1 text-center">Soft Sphere</p>
                        <div className="text-xs">
                            <BlockMath math={'V(r) = \\epsilon \\left( \\frac{\\sigma}{r} \\right)^n'} />
                        </div>
                        <p className="text-[10px] mt-2 text-muted-foreground text-center">Generic inverse-power repulsion.</p>
                    </>
                )
            default:
                return null
        }
    }

  return (
    <div className="container mx-auto p-4 md:p-8 px-4 md:px-16">
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        
        {/* Left Column: Controls */}
        <div className="lg:col-span-1 space-y-6">
          <Card>
            <CardHeader>
              <CardTitle>System Controls</CardTitle>
            </CardHeader>
            <CardContent className="space-y-6">
              
              <div className="flex gap-2">
                <Button 
                    onClick={handleStartStop} 
                    className={`flex-1 ${running ? 'bg-amber-600 hover:bg-amber-700' : 'bg-green-600 hover:bg-green-700'}`}
                >
                    {running ? <><Pause className="mr-2 h-4 w-4"/> Pause</> : <><Play className="mr-2 h-4 w-4"/> Start</>}
                </Button>
                <Button variant="outline" size="icon" onClick={handleReset} title="Reset">
                    <RotateCcw className="h-4 w-4" />
                </Button>
                <Button variant="outline" size="icon" onClick={() => initializeParticles()} title="Zap (Re-Minimize)">
                    <Zap className="h-4 w-4 text-yellow-500" />
                </Button>
              </div>

              <div className="space-y-4">
                {/* PARTICLE COUNT: Unsafe during run */}
                <div className={`space-y-2 ${running ? 'opacity-50 pointer-events-none' : ''}`}>
                  <Label>Particle Count: {numParticles}</Label>
                  <Slider 
                    disabled={running}
                    value={[numParticles]} min={10} max={200} step={10}
                    onValueChange={(v) => handleParticleCountChange(v[0])}
                  />
                </div>

                {/* TEMPERATURE: SAFE to change during run */}
                <div className="space-y-2">
                  <Label>Target Temp (T): {temperature.toFixed(1)} K</Label>
                  <Slider 
                    value={[temperature]} 
                    min={0.1} 
                    max={5.0} 
                    step={0.1}
                    onValueChange={(v) => setTemperature(v[0])}
                  />
                </div>

                {/* TIME STEP: Safe to change during run */}
                <div className="space-y-2">
                    <Label>Time Step (dt): {timeStep.toFixed(4)} ps</Label>
                    <Slider
                        value={[timeStep]} min={0.0005} max={0.05} step={0.0005}
                        onValueChange={(v) => setTimeStep(v[0])}
                    />
                </div>

                {/* 1. POTENTIAL SELECTOR WITH TOOLTIP */}
                <div className={`space-y-2 pt-2 border-t ${running ? 'opacity-50 pointer-events-none' : ''}`}>
                    <div className="flex items-center justify-between">
                        <Label>Interatomic Potential</Label>
                    </div>

                    <div className="flex items-center gap-2">
                        <select 
                            className="flex-1 h-9 items-center rounded-md border border-input bg-background px-3 py-2 text-sm shadow-sm ring-offset-background placeholder:text-muted-foreground focus:outline-none focus:ring-1 focus:ring-ring disabled:cursor-not-allowed disabled:opacity-50"
                            value={potentialType}
                            onChange={(e) => setPotentialType(e.target.value as any)}
                            disabled={running}
                        >
                            <option value="LJ">Lennard-Jones (Standard)</option>
                            <option value="WCA">WCA (Pure Repulsion)</option>
                            <option value="MORSE">Morse (Bonding)</option>
                            <option value="SOFT">Soft Sphere (Power Law)</option>
                            <option value="YUKAWA">Yukawa (Screened Coulomb)</option>
                            <option value="GAUSS">Gaussian Core (Polymers)</option>
                            <option value="BUCK">Buckingham (Ionic)</option>
                        </select>

                        {/* Info icon placed to the right of the select */}
                        <div className="relative group flex items-center justify-center">
                            <Info className="h-4 w-4 text-muted-foreground cursor-help hover:text-foreground transition-colors" />
                            <div className="absolute left-1/2 -translate-x-1/2 top-full mt-2 z-50 hidden w-80 p-3 text-sm bg-popover text-popover-foreground border rounded-md shadow-md group-hover:block animate-in fade-in zoom-in-95 duration-200">
                                {getEquationInfo(potentialType)}
                            </div>
                        </div>
                    </div>
                </div>

                {/* EPSILON / DEPTH / STRENGTH */}
                <div className={`space-y-2 ${running ? 'opacity-50 pointer-events-none' : ''}`}>
                    <Label className="text-white">
                        {labels.eps}: {epsilon.toFixed(1)}
                    </Label>
                    <Slider
                        disabled={running}
                        value={[epsilon]} min={1.0} max={200.0} step={1.0}
                        onValueChange={(v) => setEpsilon(v[0])}
                    />
                </div>

                {/* SIGMA / EQUILIBRIUM / RADIUS */}
                <div className={`space-y-2 ${running ? 'opacity-50 pointer-events-none' : ''}`}>
                    <Label className="text-white">
                        {labels.sig}: {sigma.toFixed(1)}
                    </Label>
                    <Slider
                        disabled={running}
                        value={[sigma]} min={1.0} max={15.0} step={0.5}
                        onValueChange={(v) => setSigma(v[0])}
                    />
                </div>

                {/* DYNAMIC "PARAMETER C" SLIDER */}
                {(potentialType === 'MORSE' || potentialType === 'SOFT' || potentialType === 'BUCK') && (
                    <div className={`space-y-2 ${running ? 'opacity-50 pointer-events-none' : ''}`}>
                        <Label className="text-white">
                            {potentialType === 'MORSE' ? `Width (a): ${paramC.toFixed(2)}` : 
                             potentialType === 'BUCK' ? `Attraction (C): ${paramC.toFixed(1)}` :
                             `Exponent (n): ${paramC.toFixed(1)}`}
                        </Label>
                        <Slider
                            disabled={running}
                            value={[paramC]} 
                            min={potentialType === 'MORSE' ? 0.5 : (potentialType === 'BUCK' ? 0.0 : 1.0)} 
                            max={potentialType === 'MORSE' ? 3.0 : (potentialType === 'BUCK' ? 500.0 : 18.0)} 
                            step={potentialType === 'BUCK' ? 10.0 : 0.1}
                            onValueChange={(v) => setParamC(v[0])}
                        />
                    </div>
                )}
              </div>

              <div className="pt-4 border-t space-y-4">
                 <div className="grid grid-cols-2 gap-2 text-center">
                    <div className="p-2 bg-muted rounded">
                        <p className="text-[10px] text-muted-foreground">Thermostat (ζ)</p>
                        <p className={`font-mono font-bold text-xs ${thermostatZeta > 0 ? 'text-red-500' : 'text-blue-500'}`}>
                            {thermostatZeta.toFixed(3)}
                        </p>
                    </div>
                    <div className="p-2 bg-muted rounded">
                        <p className="text-[10px] text-muted-foreground">Pressure (P)</p>
                        <p className="font-mono font-bold text-xs text-purple-600">
                            {pressure.toFixed(3)} ε/σ²
                        </p>
                    </div>
                </div>

                <div className="grid grid-cols-2 gap-2">
                    <Button
                        variant={showVelocity ? "default" : "outline"}
                        onClick={() => setShowVelocity(!showVelocity)}
                        className={`text-xs cursor-pointer ${showVelocity ? 'bg-green-600 text-white' : 'text-green-600'}`}
                        style={mounted ? (showVelocity ? { backgroundColor: velColor, color: '#ffffff' } : { color: velColor }) : undefined}
                    >
                        Velocity
                    </Button>
                    <Button
                        variant={showForce ? "default" : "outline"}
                        onClick={() => setShowForce(!showForce)}
                        className={`text-xs cursor-pointer ${showForce ? 'bg-amber-500 text-black' : 'text-amber-500'}`}
                        style={mounted ? (showForce ? { backgroundColor: forceColor, color: '#111827' } : { color: forceColor }) : undefined}
                    >
                        Force
                    </Button>
                </div>
              </div>
            </CardContent>
          </Card>

          {/* Selected Particle Info */}
          <Card>
              <CardHeader className="py-3">
                  <CardTitle className="text-sm">Selected Particle</CardTitle>
              </CardHeader>
              
              <CardContent className="pt-0">
                {selectedParticleData ? (
                    <div className="grid grid-cols-2 gap-3 text-sm">
                        {/* Magnitudes */}
                        <div className="p-2 bg-muted rounded col-span-1">
                            <p className="text-muted-foreground text-xs mb-1">Speed</p>
                            <p className="font-mono font-bold">
                                {Math.sqrt(selectedParticleData.vx**2 + selectedParticleData.vy**2 + selectedParticleData.vz**2).toFixed(2)} σ/ps
                            </p>
                        </div>
                        <div className="p-2 bg-muted rounded col-span-1">
                            <p className="text-muted-foreground text-xs mb-1">Force</p>
                            <p className="font-mono font-bold">
                                {Math.sqrt(selectedParticleData.fx**2 + selectedParticleData.fy**2 + selectedParticleData.fz**2).toFixed(2)} ε/σ
                            </p>
                        </div>

                        {/* Velocity Components */}
                        <div className={`col-span-2 grid ${isThreeD ? 'grid-cols-3' : 'grid-cols-2'} gap-2`}>
                             <div className="p-2 bg-muted rounded text-center">
                                <p className="text-muted-foreground text-[10px] mb-1">Vx</p>
                                <p className="font-mono font-bold text-xs">{selectedParticleData.vx.toFixed(2)}</p>
                             </div>
                             <div className="p-2 bg-muted rounded text-center">
                                <p className="text-muted-foreground text-[10px] mb-1">Vy</p>
                                <p className="font-mono font-bold text-xs">{selectedParticleData.vy.toFixed(2)}</p>
                             </div>
                             {isThreeD && (
                                 <div className="p-2 bg-muted rounded text-center">
                                    <p className="text-muted-foreground text-[10px] mb-1">Vz</p>
                                    <p className="font-mono font-bold text-xs">{selectedParticleData.vz.toFixed(2)}</p>
                                 </div>
                             )}
                        </div>

                        {/* Force Components */}
                        <div className={`col-span-2 grid ${isThreeD ? 'grid-cols-3' : 'grid-cols-2'} gap-2`}>
                             <div className="p-2 bg-muted rounded text-center">
                                <p className="text-muted-foreground text-[10px] mb-1">Fx</p>
                                <p className="font-mono font-bold text-xs">{selectedParticleData.fx.toFixed(2)}</p>
                             </div>
                             <div className="p-2 bg-muted rounded text-center">
                                <p className="text-muted-foreground text-[10px] mb-1">Fy</p>
                                <p className="font-mono font-bold text-xs">{selectedParticleData.fy.toFixed(2)}</p>
                             </div>
                             {isThreeD && (
                                 <div className="p-2 bg-muted rounded text-center">
                                    <p className="text-muted-foreground text-[10px] mb-1">Fz</p>
                                    <p className="font-mono font-bold text-xs">{selectedParticleData.fz.toFixed(2)}</p>
                                 </div>
                             )}
                        </div>
                    </div>
                ) : (
                    <p className="text-xs text-muted-foreground">
                        {running 
                            ? "Pause simulation to select a particle." 
                            : "Click on a particle to view its properties."}
                    </p>
                )}
              </CardContent>
          </Card>
        </div>

        {/* Right Column: Simulation & RDF */}
        <div className="lg:col-span-2 flex flex-col gap-6">
           <Card className="w-full p-0 overflow-hidden bg-muted/20 border-border relative" style={{ height: '500px' }}>
              {/* NEW: 2D/3D Switch Overlay (Top Right) */}
              <div className="absolute top-2 right-2 z-10 flex gap-2 items-center">
                  <Button
                      size="sm"
                      variant={!isThreeD ? 'default' : 'outline'}
                      onClick={() => {
                          setIsThreeD(false)
                          setRunning(false) // Stop simulation on switch for safety
                      }}
                      className={`h-8 text-xs shadow-none hover:bg-background/80 active:bg-background/80 focus-visible:ring-0 focus-visible:outline-none font-[Merriweather_Sans] ${
                          mounted ? (isDarkTheme
                              ? 'bg-background/80 text-white border border-white/20'
                              : 'bg-background/80 text-black border border-black/20')
                          : ''
                      }`}
                  >
                      2D View
                  </Button>
                  <Button
                      size="sm"
                      variant={isThreeD ? 'default' : 'outline'}
                      onClick={() => {
                          setIsThreeD(true)
                          setRunning(false)
                      }}
                      className={`h-8 text-xs shadow-none hover:bg-background/80 active:bg-background/80 focus-visible:ring-0 focus-visible:outline-none font-[Merriweather_Sans] ${
                          mounted ? (isDarkTheme
                              ? 'bg-background/80 text-white border border-white/20'
                              : 'bg-background/80 text-black border border-black/20')
                          : ''
                      }`}
                  >
                      3D View
                  </Button>
              </div>

              {/* Chart Component Swap */}
              {mounted && (
                isThreeD ? (
                    <div className="h-full w-full">
                        <Canvas>
                            <PerspectiveCamera makeDefault position={[BOX_SIZE * 2.5, BOX_SIZE * 2, BOX_SIZE * 2.5]} fov={45} />
                            <OrbitControls target={[BOX_SIZE/2, BOX_SIZE/2, BOX_SIZE/2]} />
                            <Molecule3D 
                                particles={particles} 
                                boxSize={BOX_SIZE} 
                                sigma={sigma} 
                                selectedId={selectedId}
                                isDark={isDarkTheme}
                                showVelocity={showVelocity}
                                showForce={showForce}
                                onParticleClick={(id) => {
                                    if (!running) {
                                        setSelectedId(prev => prev === id ? null : id)
                                    }
                                }}
                            />
                        </Canvas>
                    </div>
                ) : (
                    <ReactECharts
                        key='2d-mode'
                        option={getParticleOption()}
                        style={{ height: '100%', width: '100%' }}
                        onEvents={{ click: onChartClick }}
                        notMerge={false}
                        lazyUpdate={true}
                    />
                )
              )}
           </Card>

           <Card className="relative">
              <CardContent className="p-1 h-[320px]">
                  {/* Energy Mode Switcher (Top Left) - Only show for Energy graph */}
                  {analysisView === 'energy' && (
                      <div className="absolute top-2 left-2 z-10 flex gap-2 items-center">
                          <Button
                              size="sm"
                              variant={energyMode === 'sliding' ? 'default' : 'outline'}
                              onClick={() => setEnergyMode('sliding')}
                              className={`h-8 text-xs shadow-none hover:bg-background/80 active:bg-background/80 focus-visible:ring-0 focus-visible:outline-none font-[Merriweather_Sans] ${
                                  mounted ? (isDarkTheme
                                      ? 'bg-background/80 text-white border border-white/20'
                                      : 'bg-background/80 text-black border border-black/20')
                                  : ''
                              }`}
                          >
                              Sliding
                          </Button>
                          <Button
                              size="sm"
                              variant={energyMode === 'expanding' ? 'default' : 'outline'}
                              onClick={() => setEnergyMode('expanding')}
                              className={`h-8 text-xs shadow-none hover:bg-background/80 active:bg-background/80 focus-visible:ring-0 focus-visible:outline-none font-[Merriweather_Sans] ${
                                  mounted ? (isDarkTheme
                                      ? 'bg-background/80 text-white border border-white/20'
                                      : 'bg-background/80 text-black border border-black/20')
                                  : ''
                              }`}
                          >
                              Expanding
                          </Button>
                      </div>
                  )}

                  {/* Structure Metrics Info Card (Center Right) - Only show for RDF */}
                  {mounted && analysisView === 'rdf' && (() => {
                      const { state, color } = detectSystemState()
                      return (
                          <div className={`absolute top-1/4 -translate-y-1/2 right-2 z-10 backdrop-blur-sm p-2 rounded border text-[10px] space-y-1 shadow-sm w-40 ${
                              isDarkTheme
                                  ? 'bg-background/80 border-white/20'
                                  : 'bg-background/80 border-black/20'
                          }`}>
                              <div className="flex justify-between">
                                  <span className="text-muted-foreground">1st Peak (r):</span>
                                  <span className="font-mono font-bold">{structureMetrics.peakDist.toFixed(2)} σ</span>
                              </div>
                              <div className="flex justify-between">
                                  <span className="text-muted-foreground">Coordination Number:</span>
                                  <span className="font-mono font-bold text-blue-500">{structureMetrics.coordination.toFixed(2)}</span>
                              </div>
                              <div className="flex justify-between border-t pt-1 mt-1">
                                  <span className="text-muted-foreground">System State:</span>
                                  <span className={`font-bold ${color}`}>{state}</span>
                              </div>
                          </div>
                      )
                  })()}

                  <div className="absolute top-2 right-2 z-10 flex gap-2 items-center">
                      <Button
                          size="sm"
                          variant={analysisView === 'rdf' ? 'default' : 'outline'}
                          onClick={() => setAnalysisView('rdf')}
                          className={`h-8 text-xs shadow-none hover:bg-background/80 active:bg-background/80 focus-visible:ring-0 focus-visible:outline-none font-[Merriweather_Sans] ${
                              mounted ? (isDarkTheme
                                  ? 'bg-background/80 text-white border border-white/20'
                                  : 'bg-background/80 text-black border border-black/20')
                              : ''
                          }`}
                      >
                          RDF
                      </Button>
                      <Button
                          size="sm"
                          variant={analysisView === 'energy' ? 'default' : 'outline'}
                          onClick={() => setAnalysisView('energy')}
                          className={`h-8 text-xs shadow-none hover:bg-background/80 active:bg-background/80 focus-visible:ring-0 focus-visible:outline-none font-[Merriweather_Sans] ${
                              mounted ? (isDarkTheme
                                  ? 'bg-background/80 text-white border border-white/20'
                                  : 'bg-background/80 text-black border border-black/20')
                              : ''
                          }`}
                      >
                          Energy
                      </Button>
                      <Button
                          size="sm"
                          variant={analysisView === 'msd' ? 'default' : 'outline'}
                          onClick={() => setAnalysisView('msd')}
                          className={`h-8 text-xs shadow-none hover:bg-background/80 active:bg-background/80 focus-visible:ring-0 focus-visible:outline-none font-[Merriweather_Sans] ${
                              mounted ? (isDarkTheme
                                  ? 'bg-background/80 text-white border border-white/20'
                                  : 'bg-background/80 text-black border border-black/20')
                              : ''
                          }`}
                      >
                          MSD
                      </Button>
                  </div>
                  <div
                      className="h-full w-full"
                      onMouseEnter={() => {
                          analysisHoverRef.current = true
                          startAnalysisTooltipLoop()
                          refreshAnalysisTooltip()
                      }}
                      onMouseMove={(e) => {
                          analysisHoverRef.current = true
                          const rect = (e.currentTarget as HTMLDivElement).getBoundingClientRect()
                          analysisPointerRef.current = { x: e.clientX - rect.left, y: e.clientY - rect.top }
                          startAnalysisTooltipLoop()
                          refreshAnalysisTooltip()
                      }}
                      onMouseLeave={() => {
                          analysisHoverRef.current = false
                          analysisPointerRef.current = null
                          stopAnalysisTooltipLoop()
                          const inst = analysisChartRef.current?.getEchartsInstance?.()
                          inst?.dispatchAction({ type: 'hideTip' })
                      }}
                  >
                      <ReactECharts
                          key={analysisView}
                          ref={analysisChartRef}
                          option={
                              analysisView === 'rdf' ? getRdfOption() : 
                              analysisView === 'energy' ? getEnergyOption() : 
                              getMsdOption()
                          }
                          style={{ height: '100%', width: '100%' }}
                          notMerge={false}
                          lazyUpdate={true}
                      />
                  </div>
              </CardContent>
           </Card>
        </div>

      </div>
    </div>
  )

}

// Helper component for rendering 3D arrows
const Arrow3D = ({ 
  start, 
  vector, 
  color, 
  scale 
}: { 
  start: [number, number, number], 
  vector: [number, number, number], 
  color: string, 
  scale: number 
}) => {
  const dir = new THREE.Vector3(...vector)
  const length = dir.length()
  
  // Don't render tiny vectors
  if (length < 1e-6) return null

  dir.normalize()
  const arrowLength = length * scale
  
  // Create the ArrowHelper once per render cycle
  // args: [dir, origin, length, color, headLength, headWidth]
  return (
    <primitive 
      object={new THREE.ArrowHelper(
        dir, 
        new THREE.Vector3(...start), 
        arrowLength, 
        color, 
        Math.min(arrowLength * 0.2, 2), // Cap head size
        Math.min(arrowLength * 0.08, 1)
      )} 
    />
  )
}

const Molecule3D = ({ 
  particles, 
  boxSize, 
  sigma, 
  selectedId, 
  isDark,
  showVelocity,
  showForce,
  onParticleClick
}: { 
  particles: Particle[], 
  boxSize: number, 
  sigma: number, 
  selectedId: number | null, 
  isDark: boolean,
  showVelocity: boolean,
  showForce: boolean,
  onParticleClick: (id: number) => void
}) => {
  const particleColor = isDark ? "#ffffff" : "#3b82f6"
  const wireframeColor = isDark ? "#ffffff" : "#000000"
  
  // Colors for arrows
  const velColor = isDark ? '#4ade80' : '#16a34a'
  const forceColor = isDark ? '#facc15' : '#d97706'

  // Arrow Scaling Factors (matching 2D logic)
  const VEL_SCALE = 0.5
  const FORCE_SCALE = 0.05

  return (
    <group>
      {/* Simulation Bounding Box */}
      <mesh position={[boxSize/2, boxSize/2, boxSize/2]}>
        <boxGeometry args={[boxSize, boxSize, boxSize]} />
        <meshBasicMaterial 
          color={wireframeColor}
          wireframe 
          transparent 
          opacity={0.15} 
        />
      </mesh>
      
      {/* Particles */}
      {particles.map((p) => {
        const isSelected = p.id === selectedId
        
        // Visibility Logic:
        // 1. If a particle is selected, ONLY show arrows for that particle.
        // 2. If nothing is selected, show arrows based on global toggles.
        const showV = selectedId !== null ? isSelected : showVelocity
        const showF = selectedId !== null ? isSelected : showForce

        return (
          <group key={p.id}>
             {/* The Particle Sphere */}
             <mesh 
                position={[p.x, p.y, p.z]} 
                onClick={(e) => {
                  e.stopPropagation()
                  onParticleClick(p.id)
                }}
             >
              <sphereGeometry args={[sigma * 0.5, 16, 16]} />
              <meshStandardMaterial 
                color={isSelected ? "#facc15" : particleColor}
                emissive={isSelected ? "#facc15" : "#000000"}
                emissiveIntensity={isSelected ? 0.5 : 0}
                roughness={0.2}
                metalness={0.5}
              />
            </mesh>

            {/* Velocity Arrow */}
            {showV && (
              <Arrow3D 
                start={[p.x, p.y, p.z]} 
                vector={[p.vx, p.vy, p.vz]} 
                color={velColor} 
                scale={VEL_SCALE} 
              />
            )}

            {/* Force Arrow */}
            {showF && (
              <Arrow3D 
                start={[p.x, p.y, p.z]} 
                vector={[p.fx, p.fy, p.fz]} 
                color={forceColor} 
                scale={FORCE_SCALE} 
              />
            )}
          </group>
        )
      })}
      
      <ambientLight intensity={0.6} />
      <pointLight position={[boxSize * 1.5, boxSize * 1.5, boxSize * 1.5]} intensity={1} />
    </group>
  )
}