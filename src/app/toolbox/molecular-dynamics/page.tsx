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
import { Play, Pause, RotateCcw, Move, ArrowRight, MousePointer2, Activity, List, Zap } from 'lucide-react'

// Types for our simulation
type Particle = {
  id: number
  x: number
  y: number
  vx: number
  vy: number
  fx: number
  fy: number
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
  
  // --- Control States ---
  const [numParticles, setNumParticles] = useState<number>(50)
  const [temperature, setTemperature] = useState<number>(1.0)
  const [timeStep, setTimeStep] = useState<number>(DT_DEFAULT)
  
  // New L-J Parameters
  const [epsilon, setEpsilon] = useState<number>(50.0)
  const [sigma, setSigma] = useState<number>(5.0)

  // --- Visualization States ---
  const [showForce, setShowForce] = useState(false)
  const [showVelocity, setShowVelocity] = useState(false)
  
  // --- Analysis State ---
  const [kineticEnergy, setKineticEnergy] = useState<number>(0)
  const [thermostatZeta, setThermostatZeta] = useState<number>(0) // Visualizing the friction
  const [rdfData, setRdfData] = useState<{r: number, g: number}[]>([])
  const [energyHistory, setEnergyHistory] = useState<{time: number, ke: number, pe: number, tot: number}[]>([])
  const [analysisView, setAnalysisView] = useState<'rdf' | 'energy'>('rdf')
  
  // --- Selected Particle View State ---
  const particleHistoryRef = useRef<{time: number, vx: number, vy: number, f: number}[]>([])
  const [historyTrigger, setHistoryTrigger] = useState(0) 

  // Refs for Simulation Loop
  const requestRef = useRef<number | null>(null)
  const particlesRef = useRef<Particle[]>([])
  
  // Refs for current parameter values
  const paramsRef = useRef({
      dt: DT_DEFAULT,
      epsilon: 50.0,
      sigma: 5.0,
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
      paramsRef.current.targetTemp = temperature
      
      // Update NH Thermal Mass Q roughly based on system size
      // Q approx N * T * tau^2. Let tau = 0.5ps
      nhStateRef.current.Q = numParticles * temperature * 0.25
  }, [timeStep, epsilon, sigma, temperature, numParticles])

  // --- PHYSICS KERNEL: Force Calculation ---
  // Abstracted so it can be used by both MD Loop and Minimizer
  const calculateForces = (currentParticles: Particle[]): number => {
      const { epsilon: eps, sigma: sig, box } = paramsRef.current
      
      let potentialEnergy = 0
      
      // 1. Reset Forces
      for (let p of currentParticles) {
          p.fx = 0
          p.fy = 0
      }

      // 2. Cell Lists Setup
      const cutoff = 2.5 * sig
      const cellSize = Math.max(cutoff, 2.0)
      const gridCols = Math.floor(box / cellSize)
      const gridRows = Math.floor(box / cellSize)
      
      if (gridCols <= 0 || gridRows <= 0) return 0

      const getGridIndex = (x: number, y: number) => {
          const cx = Math.max(0, Math.min(gridCols - 1, Math.floor(x / cellSize)))
          const cy = Math.max(0, Math.min(gridRows - 1, Math.floor(y / cellSize)))
          return cy * gridCols + cx
      }

      const grid: number[][] = Array.from({ length: gridCols * gridRows }, () => [])
      for (let p of currentParticles) {
          const idx = getGridIndex(p.x, p.y)
          grid[idx].push(p.id)
      }

      // 3. Pairwise Interactions
      for (let p1 of currentParticles) {
          const cx = Math.max(0, Math.min(gridCols - 1, Math.floor(p1.x / cellSize)))
          const cy = Math.max(0, Math.min(gridRows - 1, Math.floor(p1.y / cellSize)))

          for (let i = -1; i <= 1; i++) {
              for (let j = -1; j <= 1; j++) {
                  let nx = cx + i
                  let ny = cy + j
                  
                  // Periodic wrapping
                  if (nx < 0) nx += gridCols
                  if (nx >= gridCols) nx -= gridCols
                  if (ny < 0) ny += gridRows
                  if (ny >= gridRows) ny -= gridRows

                  const cellParticles = grid[ny * gridCols + nx]
                  if (!cellParticles) continue

                  for (let otherId of cellParticles) {
                      if (otherId <= p1.id) continue
                      const p2 = currentParticles[otherId]
                      if (!p2) continue

                      let dx = p2.x - p1.x
                      let dy = p2.y - p1.y

                      // Minimum Image Convention
                      dx -= box * Math.round(dx / box)
                      dy -= box * Math.round(dy / box)

                      const distSq = dx * dx + dy * dy

                      if (distSq > cutoff * cutoff) continue

                      // STRICT: No force clamping needed due to minimization
                      const r2 = Math.max(1e-6, distSq)
                      const invR2 = 1.0 / r2
                      const s2 = sig * sig
                      const s2_r2 = s2 * invR2
                      const s6_r6 = s2_r2 * s2_r2 * s2_r2
                      const s12_r12 = s6_r6 * s6_r6

                      const forceScalar = (24 * eps * invR2) * (2 * s12_r12 - s6_r6)
                      
                      // Lennard-Jones Potential Energy: V = 4*eps*(s12 - s6)
                      potentialEnergy += 4 * eps * (s12_r12 - s6_r6)

                      const fx = dx * forceScalar
                      const fy = dy * forceScalar

                      p1.fx -= fx
                      p1.fy -= fy
                      p2.fx += fx
                      p2.fy += fy
                  }
              }
          }
      }
      
      return potentialEnergy
  }

  // --- ALGORITHM 1: Steepest Descent Minimization ---
  // Strictly moves particles "downhill" to fix overlaps
  const runMinimization = (pList: Particle[]) => {
      const MAX_STEPS = 200
      const GAMMA = 0.01 // Descent rate
      const box = paramsRef.current.box

      for (let step = 0; step < MAX_STEPS; step++) {
          calculateForces(pList)
          
          let maxForce = 0
          for (let p of pList) {
              const fMag = Math.sqrt(p.fx * p.fx + p.fy * p.fy)
              if (fMag > maxForce) maxForce = fMag
              
              // Move particle along force vector (Gradient Descent)
              const scale = Math.min(GAMMA, 0.5 / (fMag + 1e-6))
              
              p.x += p.fx * scale
              p.y += p.fy * scale
              
              // Enforce PBC
              if (p.x < 0) p.x += box
              if (p.x >= box) p.x -= box
              if (p.y < 0) p.y += box
              if (p.y >= box) p.y -= box
              
              // Zero out velocity (Descent extracts energy)
              p.vx = 0
              p.vy = 0
          }
          if (maxForce < 10.0) break // Converged
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
          const r = (i + 0.5) * RDF_BIN_WIDTH
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

  // --- Initialization Logic ---
  const initializeParticles = useCallback(() => {
    const newParticles: Particle[] = []
    const gridDim = Math.ceil(Math.sqrt(numParticles))
    const spacing = BOX_SIZE / gridDim

    // 1. Grid Placement
    for (let i = 0; i < numParticles; i++) {
      newParticles.push({
        id: i,
        x: (i % gridDim) * spacing + spacing / 2 + (Math.random() - 0.5),
        y: Math.floor(i / gridDim) * spacing + spacing / 2 + (Math.random() - 0.5),
        vx: 0,
        vy: 0,
        fx: 0,
        fy: 0
      })
    }
    
    // 2. Run Energy Minimization (Steepest Descent)
    // This fixes the "Exploding Atom" problem scientifically
    runMinimization(newParticles)

    // 3. Assign Maxwell-Boltzmann Velocities
    for (let p of newParticles) {
        p.vx = (Math.random() - 0.5) * temperature
        p.vy = (Math.random() - 0.5) * temperature
    }

    // 4. Remove Center of Mass Motion (Optional but good practice)
    let vcmX = 0, vcmY = 0
    for (let p of newParticles) { vcmX += p.vx; vcmY += p.vy }
    vcmX /= numParticles; vcmY /= numParticles
    for (let p of newParticles) { p.vx -= vcmX; p.vy -= vcmY }
    
    setParticles(newParticles)
    particlesRef.current = newParticles
    setSelectedId(null)
    particleHistoryRef.current = []
    rdfHistogramRef.current = new Array(RDF_BINS).fill(0) // Clear RDF bins
    rdfFramesRef.current = 0 // Reset RDF frame counter
    setRdfData([]) // Clear RDF display
    timeRef.current = 0
    nhStateRef.current.zeta = 0 // Reset Thermostat
  }, [numParticles, temperature])

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
    const dof = 2 * N - 2 // Degrees of freedom (minus constraints)
    
    // 1. Calculate Forces at t
    calculateForces(currentParticles)

    // 2. First Half-step Update
    // v(t + dt/2) and zeta(t + dt/2)
    // We use a simplified explicit update for stability in JS
    
    for (let p of currentParticles) {
        // Update Velocity (Half Step) with thermostat
        // v <- v + 0.5 * dt * (F - zeta * v)
        p.vx += 0.5 * dt * (p.fx - zeta * p.vx)
        p.vy += 0.5 * dt * (p.fy - zeta * p.vy)
        
        // Update Position (Full Step)
        p.x += p.vx * dt
        p.y += p.vy * dt

        // PBC
        if (p.x < 0) p.x += box
        if (p.x >= box) p.x -= box
        if (p.y < 0) p.y += box
        if (p.y >= box) p.y -= box
    }

    // 3. Calculate Forces at t + dt
    const pe = calculateForces(currentParticles)

    // 4. Update Thermostat Variable (Zeta)
    // Based on kinetic energy at half-step (approx using current v)
    let keSum = 0
    for (let p of currentParticles) {
        keSum += 0.5 * (p.vx * p.vx + p.vy * p.vy)
    }
    const currentTemp = (2.0 * keSum) / dof
    
    // Zeta integration: dZeta/dt = (KE - TargetKE) / Q
    zeta += dt * (dof * BOLTZMANN_K * (currentTemp - targetTemp)) / Q

    // 5. Second Half-step Update
    let kineticEnergyFull = 0
    for (let p of currentParticles) {
        // Explicit approximation for v(t+dt)
        // Note: Strict NH requires iterative solution here
        const vxh = p.vx
        const vyh = p.vy
        
        p.vx = (vxh + 0.5 * dt * p.fx) / (1.0 + 0.5 * dt * zeta)
        p.vy = (vyh + 0.5 * dt * p.fy) / (1.0 + 0.5 * dt * zeta)
        
        // STRICT: Removed MAX_VELOCITY
        // If simulation explodes, it's because dt is too high, which is scientifically accurate.
        
        kineticEnergyFull += 0.5 * (p.vx * p.vx + p.vy * p.vy)
    }

    nhStateRef.current.zeta = zeta

    // Update State
    particlesRef.current = currentParticles
    setKineticEnergy(kineticEnergyFull)
    setThermostatZeta(zeta)
    
    // Sample RDF statistics (every step)
    sampleRDF(currentParticles)

    setFrameCount(prev => {
        const next = prev + 1
        if (next % 4 === 0) { // Render every 4th physics frame
            setParticles(currentParticles)
            
            // Compute and Update RDF Graph Data (every 4th frame)
            const newRdfData = computeRDFGraph()
            setRdfData(newRdfData)
            
            // Track Energy History
            setEnergyHistory(prev => {
                const newHistory = [...prev, { 
                    time: timeRef.current, 
                    ke: kineticEnergyFull, 
                    pe: pe, 
                    tot: kineticEnergyFull + pe 
                }]
                if (newHistory.length > 100) newHistory.shift()
                return newHistory
            })
            
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
  }, [numParticles])

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

    const renderArrow = (params: any, api: any) => {
        const x = api.value(0)
        const y = api.value(1)
        const dxVal = api.value(2)
        const dyVal = api.value(3)
        const type = api.value(4)
        const id = api.value(5)
        
        const isSelected = id === selectedId
        
        // If ANY particle is selected, only show arrows for THAT particle.
        // If NO particle is selected, show arrows based on global toggles.
        let isVisible = false
        if (selectedId !== null) {
            isVisible = isSelected
        } else {
            isVisible = (type === 0 && showVelocity) || (type === 1 && showForce)
        }

        if (!isVisible) return

        if (Math.abs(dxVal) < 0.1 && Math.abs(dyVal) < 0.1) return

        const scale = type === 0 ? 0.5 : 0.05 // Reduced scale for Force vectors as L-J forces can be large
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
             // Visual size proportional to Sigma, but clamped so it doesn't get too crazy
             // Sigma 5.0 -> approx 12px
             return Math.max(5, Math.min(30, sigma * 2.5))
          },
          data: scatterData,
          emphasis: {
             scale: false,
             focus: 'none'
          },
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
      
      return {
          backgroundColor: 'transparent',
          title: {
              text: 'Radial Distribution Function g(r)',
              left: 'center',
              textStyle: { fontSize: 12, color: isDark ? '#fff' : '#000', fontFamily: 'Merriweather Sans' }
          },
          tooltip: { 
              trigger: 'axis',
              confine: true,
              enterable: true,
              backgroundColor: isDark ? 'rgba(31, 41, 55, 0.95)' : 'rgba(255, 255, 255, 0.95)',
              borderColor: isDark ? '#374151' : '#e5e7eb',
              textStyle: { fontFamily: 'Merriweather Sans', color: '#ffffff' },
              formatter: (params: any) => {
                  if (!Array.isArray(params)) return ''
                  const p = params[0]
                  if (!p || p.seriesName === 'Ideal Gas') return ''
                  const r = parseFloat(p.value[0]).toFixed(3)
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
              nameTextStyle: { color: '#ffffff', fontSize: 10, fontFamily: 'Merriweather Sans' },
              axisLabel: { color: '#ffffff', fontSize: 10, fontFamily: 'Merriweather Sans' },
              axisLine: { show: true, lineStyle: { color: '#ffffff' } },
              splitLine: { show: false }
          },
          yAxis: {
              type: 'value',
              name: 'g(r)',
              min: 0,
              nameTextStyle: { color: '#ffffff', fontSize: 10, fontFamily: 'Merriweather Sans' },
              axisLabel: { color: '#ffffff', fontSize: 10, fontFamily: 'Merriweather Sans', formatter: (v: number) => v.toPrecision(3) },
              axisLine: { show: true, lineStyle: { color: '#ffffff' } },
              splitLine: { show: false }
          },
          series: [
              {
                  type: 'line',
                  showSymbol: false,
                  data: rdfData.map(d => [d.r, d.g]),
                  lineStyle: { color: '#8b5cf6', width: 2 },
                  areaStyle: { opacity: 0.1, color: '#8b5cf6' }
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
      const times = energyHistory.map(h => h.time.toFixed(1))
      const displayTimes = times.map((t, i) => (i === 0 || t !== times[i - 1]) ? t : '')
      
      return {
          title: { text: 'Energy', left: 'center', textStyle: { fontSize: 12, color: isDark ? '#fff' : '#000', fontFamily: 'Merriweather Sans' } },
          tooltip: { 
              trigger: 'axis',
              confine: true,
              enterable: true,
              backgroundColor: isDark ? 'rgba(31, 41, 55, 0.95)' : 'rgba(255, 255, 255, 0.95)',
              borderColor: isDark ? '#374151' : '#e5e7eb',
              textStyle: { fontFamily: 'Merriweather Sans', color: '#ffffff' },
              formatter: (params: any) => {
                  if (!Array.isArray(params)) return ''
                  return params.map((p: any) => `${p.seriesName}: ${Number(p.value).toPrecision(3)}`).join('<br/>')
              }
          },
          legend: { bottom: -5, symbolKeepAspect: true, itemWidth: 8, itemHeight: 8, textStyle: { color: isDark ? '#fff' : '#000', fontFamily: 'Merriweather Sans' } },
          grid: { left: 70, right: 20, top: 50, bottom: 60, containLabel: false },
          xAxis: { type: 'category', data: times, boundaryGap: false, show: true, position: 'bottom', name: 'Time (ps)', nameLocation: 'middle', nameGap: 25, nameTextStyle: { color: '#ffffff', fontSize: 10, fontFamily: 'Merriweather Sans' }, axisLine: { show: true, lineStyle: { color: '#ffffff' }, onZero: false }, axisTick: { show: true, alignWithLabel: true, interval: 'auto' }, axisLabel: { show: times.length > 0, showMinLabel: false, color: '#ffffff', margin: 12, interval: 'auto', formatter: (_: string, index: number) => displayTimes[index] } },
          yAxis: { type: 'value', name: 'Energy (ε)', nameLocation: 'middle', nameGap: 55, nameTextStyle: { color: '#ffffff', fontSize: 10, fontFamily: 'Merriweather Sans' }, splitLine: { show: false }, axisLine: { show: true, lineStyle: { color: '#ffffff' }, onZero: false }, axisTick: { show: true }, axisLabel: { color: isDark ? '#fff' : '#000', fontFamily: 'Merriweather Sans', margin: 5 } },
          series: [
              { name: 'Total E', type: 'line', data: energyHistory.map(h => h.tot), showSymbol: false, lineStyle: { width: 2, color: '#10b981' } },
              { name: 'Potential (PE)', type: 'line', data: energyHistory.map(h => h.pe), showSymbol: false, lineStyle: { width: 1.5, color: '#3b82f6' } },
              { name: 'Kinetic (KE)', type: 'line', data: energyHistory.map(h => h.ke), showSymbol: false, lineStyle: { width: 1.5, color: '#ef4444' } }
          ],
          animation: false
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

                {/* TIME STEP: Unsafe during run */}
                <div className={`space-y-2 ${running ? 'opacity-50 pointer-events-none' : ''}`}>
                    <Label>Time Step (dt): {timeStep.toFixed(4)} ps</Label>
                    <Slider
                        disabled={running}
                        value={[timeStep]} min={0.0005} max={0.005} step={0.0005}
                        onValueChange={(v) => setTimeStep(v[0])}
                    />
                </div>

                {/* EPSILON: Unsafe during run */}
                <div className={`space-y-2 pt-2 border-t ${running ? 'opacity-50 pointer-events-none' : ''}`}>
                    <Label className="text-amber-600">L-J Epsilon: {epsilon.toFixed(1)} ε</Label>
                    <Slider
                        disabled={running}
                        value={[epsilon]} min={1.0} max={200.0} step={1.0}
                        onValueChange={(v) => setEpsilon(v[0])}
                    />
                </div>

                {/* SIGMA: HIGHLY Unsafe during run */}
                <div className={`space-y-2 ${running ? 'opacity-50 pointer-events-none' : ''}`}>
                    <Label className="text-amber-600">L-J Sigma: {sigma.toFixed(1)} σ</Label>
                    <Slider
                        disabled={running}
                        value={[sigma]} min={1.0} max={15.0} step={0.5}
                        onValueChange={(v) => setSigma(v[0])}
                    />
                </div>
              </div>

              <div className="pt-4 border-t space-y-4">
                 <div className="grid grid-cols-3 gap-2 text-center">
                    <div className="p-2 bg-muted rounded">
                        <p className="text-[10px] text-muted-foreground">Total KE</p>
                        <p className="font-mono font-bold text-xs">{kineticEnergy.toFixed(1)} ε</p>
                    </div>
                    <div className="p-2 bg-muted rounded">
                        <p className="text-[10px] text-muted-foreground">Thermostat (ζ)</p>
                        <p className={`font-mono font-bold text-xs ${thermostatZeta > 0 ? 'text-red-500' : 'text-blue-500'}`}>
                            {thermostatZeta.toFixed(3)}
                        </p>
                    </div>
                    <div className="p-2 bg-muted rounded">
                        <p className="text-[10px] text-muted-foreground">Particles</p>
                        <p className="font-mono font-bold text-xs">{particles.length}</p>
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
                        <div className="p-2 bg-muted rounded">
                            <p className="text-muted-foreground text-xs mb-1">Velocity</p>
                            <p className="font-mono font-bold">{Math.sqrt(selectedParticleData.vx**2 + selectedParticleData.vy**2).toFixed(2)} σ/ps</p>
                        </div>
                        <div className="p-2 bg-muted rounded">
                            <p className="text-muted-foreground text-xs mb-1">Force</p>
                            <p className="font-mono font-bold">{Math.sqrt(selectedParticleData.fx**2 + selectedParticleData.fy**2).toFixed(2)} ε/σ</p>
                        </div>
                        <div className="p-2 bg-muted rounded">
                            <p className="text-muted-foreground text-xs mb-1">Vx</p>
                            <p className="font-mono font-bold">{selectedParticleData.vx.toFixed(2)} σ/ps</p>
                        </div>
                        <div className="p-2 bg-muted rounded">
                            <p className="text-muted-foreground text-xs mb-1">Vy</p>
                            <p className="font-mono font-bold">{selectedParticleData.vy.toFixed(2)} σ/ps</p>
                        </div>
                    </div>
                ) : (
                    <p className="text-xs text-muted-foreground">Click on a particle while paused to show its values.</p>
                )}
              </CardContent>
          </Card>
        </div>

        {/* Right Column: Simulation & RDF */}
        <div className="lg:col-span-2 flex flex-col gap-6">
           <Card className="w-full p-0 overflow-hidden bg-muted/20 border-border relative" style={{ height: '500px' }}>
              <ReactECharts
                  option={getParticleOption()}
                  style={{ height: '100%', width: '100%' }}
                  onEvents={{
                      click: onChartClick
                  }}
                  notMerge={true} 
              />
           </Card>

           <Card className="relative">
              <CardContent className="p-1 h-[320px]">
                  <div className="absolute top-2 right-2 z-10 flex gap-2 items-center">
                      <Button
                          size="sm"
                          variant={analysisView === 'rdf' ? 'default' : 'outline'}
                          onClick={() => setAnalysisView('rdf')}
                          className="h-8 text-xs bg-background/80 text-white border border-white/20 shadow-none hover:bg-background/80 active:bg-background/80 focus-visible:ring-0 focus-visible:outline-none font-[Merriweather_Sans]"
                      >
                          RDF
                      </Button>
                      <Button
                          size="sm"
                          variant={analysisView === 'energy' ? 'default' : 'outline'}
                          onClick={() => setAnalysisView('energy')}
                          className="h-8 text-xs bg-background/80 text-white border border-white/20 shadow-none hover:bg-background/80 active:bg-background/80 focus-visible:ring-0 focus-visible:outline-none font-[Merriweather_Sans]"
                      >
                          Energy
                      </Button>
                  </div>
                  <ReactECharts
                      option={analysisView === 'rdf' ? getRdfOption() : getEnergyOption()}
                      style={{ height: '100%', width: '100%' }}
                      notMerge={true}
                  />
              </CardContent>
           </Card>
        </div>

      </div>
    </div>
  )
}