'use client'

import { useState, useEffect, useCallback, useRef } from 'react'
import { useTheme } from "next-themes"

// Import ECharts components
import ReactECharts from 'echarts-for-react'
import * as echarts from 'echarts/core'
import type { EChartsOption } from 'echarts'
import { ScatterChart } from 'echarts/charts'
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  ToolboxComponent
} from 'echarts/components'
import { CanvasRenderer } from 'echarts/renderers'

// Register necessary ECharts components
echarts.use([
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  ToolboxComponent,
  ScatterChart,
  CanvasRenderer
])

import { Card, CardContent } from "@/components/ui/card"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import { Button } from "@/components/ui/button"
import { Slider } from "@/components/ui/slider"
import { Badge } from "@/components/ui/badge"
import { Tabs, TabsList, TabsTrigger } from "@/components/ui/tabs"

// --- TYPES ---
type Particle = {
  id: number
  x: number
  y: number
}

type SimulationStats = {
  step: number
  energy: number
  acceptanceRate: number
  time: number
}

type SimMode = 'metropolis' | 'kinetic'

export default function MonteCarloPage() {
  const { resolvedTheme } = useTheme()
  const isDark = resolvedTheme === 'dark'
  
  // --- STATE: Inputs ---
  const [simMode, setSimMode] = useState<SimMode>('metropolis')

  // Common Params
  const [temperature, setTemperature] = useState<number>(1.0) // Reduced Temp T*
  
  // Metropolis Specific
  const [density, setDensity] = useState<number>(0.1) // Number density rho*
  const [particleCount, setParticleCount] = useState<number>(100)

  // Kinetic Specific
  const [latticeSize, setLatticeSize] = useState<number>(20) // 20x20 grid
  const [vacancyPercent, setVacancyPercent] = useState<number>(0.1) // 10% empty spots
  const [activationEnergy, setActivationEnergy] = useState<number>(0.5) // Barrier Ea (eV)
  const [attemptFreq] = useState<number>(1e12) // nu (Hz)
  
  // --- STATE: Simulation Control ---
  const [isRunning, setIsRunning] = useState(false)
  const [stats, setStats] = useState<SimulationStats>({ step: 0, energy: 0, acceptanceRate: 0, time: 0.0 })
  
  // --- REFS: Physics Engine ---
  const particlesRef = useRef<Particle[]>([]) // Used for Metropolis
  const latticeRef = useRef<number[][]>([]) // Used for Kinetic (2D grid: 0=empty, 1=particle)
  const boxSizeRef = useRef<number>(10) 
  
  const animationFrameRef = useRef<number | null>(null)
  const stepsRef = useRef(0)
  const acceptedMovesRef = useRef(0)
  const totalMovesRef = useRef(0)
  const timeRef = useRef(0.0)
  const echartsRef = useRef<ReactECharts | null>(null)

  // ==========================================
  // PHYSICS ENGINE 1: METROPOLIS (Off-Lattice)
  // ==========================================

  const calculateLennardJonesEnergy = (p1: Particle, p2: Particle, boxSize: number, newX?: number, newY?: number): number => {
    const x1 = newX ?? p1.x
    const y1 = newY ?? p1.y
    let dx = x1 - p2.x
    let dy = y1 - p2.y
    
    // Minimum Image Convention
    dx = dx - boxSize * Math.round(dx / boxSize)
    dy = dy - boxSize * Math.round(dy / boxSize)
    
    const r2 = dx * dx + dy * dy
    // Cutoff at 2.5 sigma or if particles are too close (avoid division issues)
    if (r2 > 6.25 || r2 < 1e-6) return 0
    
    const r6_inv = 1 / (r2 * r2 * r2)
    const r12_inv = r6_inv * r6_inv
    return 4 * (r12_inv - r6_inv)
  }

  const calculateTotalFluidEnergy = (): number => {
    let totalE = 0
    const ps = particlesRef.current
    const box = boxSizeRef.current
    for (let i = 0; i < ps.length; i++) {
      for (let j = i + 1; j < ps.length; j++) {
        totalE += calculateLennardJonesEnergy(ps[i], ps[j], box)
      }
    }
    return totalE
  }

  const performMetropolisSteps = () => {
    const ps = particlesRef.current
    const L = boxSizeRef.current
    const T = temperature
    const N = ps.length
    const delta = 0.5 
    
    for (let k = 0; k < N; k++) {
      const idx = Math.floor(Math.random() * N)
      const p = ps[idx]
      
      const dx = (Math.random() - 0.5) * delta
      const dy = (Math.random() - 0.5) * delta
      let newX = p.x + dx
      let newY = p.y + dy
      
      if (newX < 0) newX += L
      if (newX >= L) newX -= L
      if (newY < 0) newY += L
      if (newY >= L) newY -= L

      let currentE = 0
      let newE = 0
      for (let j = 0; j < N; j++) {
        if (idx === j) continue
        currentE += calculateLennardJonesEnergy(p, ps[j], L)
        newE += calculateLennardJonesEnergy(p, ps[j], L, newX, newY)
      }
      
      const dE = newE - currentE
      if (dE < 0 || Math.random() < Math.exp(-dE / T)) {
        p.x = newX
        p.y = newY
        acceptedMovesRef.current++
      }
      totalMovesRef.current++
    }
    stepsRef.current++
  }

  // ==========================================
  // PHYSICS ENGINE 2: KINETIC (On-Lattice)
  // ==========================================
  // Model: Vacancy mediated diffusion with NN repulsion (Ising-like)

  const performKMCSteps = () => {
    const grid = latticeRef.current
    const N = latticeSize
    const T = temperature
    const Ea = activationEnergy
    
    // 1. Calculate Rate of a single jump
    // Rate = freq * exp(-Barrier / T)
    const baseRate = attemptFreq * Math.exp(-Ea / T)
    
    // 2. Find all possible events (Vacancy jumps)
    // An event is: "Particle at (nx, ny) moves into Vacancy at (x, y)"
    const events: {vx: number, vy: number, px: number, py: number, rate: number}[] = []
    
    for(let x=0; x<N; x++) {
        for(let y=0; y<N; y++) {
            if (grid[x][y] === 0) { // Found a Vacancy
                // Check 4 neighbors
                const dirs = [[0,1], [0,-1], [1,0], [-1,0]]
                for(let d=0; d<dirs.length; d++) {
                    let nx = x + dirs[d][0]
                    let ny = y + dirs[d][1]
                    // Periodic Wrap
                    if (nx < 0) nx += N; if (nx >= N) nx -= N;
                    if (ny < 0) ny += N; if (ny >= N) ny -= N;
                    
                    if (grid[nx][ny] !== 0) {
                        // Found a particle that can jump into this vacancy
                        events.push({
                            vx: x, vy: y,   // Vacancy pos
                            px: nx, py: ny, // Particle pos
                            rate: baseRate
                        })
                    }
                }
            }
        }
    }
    
    if (events.length === 0) return // Nothing can move
    
    // 3. Calculate Total Rate (R)
    // Since we simplified to uniform rates: R = count * baseRate
    const R_total = events.length * baseRate
    
    // 4. Advance Time (The "Clock")
    // dt = -ln(random) / R_total
    const u1 = Math.random()
    const dt = -Math.log(u1) / R_total
    timeRef.current += dt
    
    // 5. Pick an Event
    // We pick a random event from the list weighted by its rate.
    // Since all our rates are equal, just pick uniform random index.
    const eventIdx = Math.floor(Math.random() * events.length)
    const event = events[eventIdx]
    
    // 6. Execute Move
    grid[event.vx][event.vy] = grid[event.px][event.py] // Move particle to vacancy
    grid[event.px][event.py] = 0 // Old spot becomes vacancy
    
    acceptedMovesRef.current++
    totalMovesRef.current++
    stepsRef.current++
  }

  // --- INITIALIZATION ---
  const initializeSystem = useCallback(() => {
    setIsRunning(false)
    stepsRef.current = 0
    acceptedMovesRef.current = 0
    totalMovesRef.current = 0
    timeRef.current = 0.0
    
    if (simMode === 'metropolis') {
        // FLUID INIT
        const safeDensity = Math.max(0.01, density)
        const L = Math.sqrt(particleCount / safeDensity)
        boxSizeRef.current = L

        const newParticles: Particle[] = []
        const gridSize = Math.ceil(Math.sqrt(particleCount))
        const spacing = L / gridSize
        
        let count = 0
        for (let i = 0; i < gridSize; i++) {
          for (let j = 0; j < gridSize; j++) {
            if (count >= particleCount) break
            newParticles.push({
              id: count,
              x: (i + 0.5) * spacing,
              y: (j + 0.5) * spacing
            })
            count++
          }
        }
        particlesRef.current = newParticles
        const initialE = calculateTotalFluidEnergy()
        setStats({ step: 0, energy: initialE, acceptanceRate: 0, time: 0.0 })

    } else {
        // LATTICE INIT (Binary Alloy)
        boxSizeRef.current = latticeSize
        const N = latticeSize
        const totalSites = N * N
        // 10% Vacancies, 45% Type A, 45% Type B
        const numVacancies = Math.floor(totalSites * vacancyPercent)
        
        // Create grid
        const newGrid = Array(N).fill(0).map(() => Array(N).fill(0))
        
        // Linear array of types to shuffle
        const types = []
        for(let i=0; i<numVacancies; i++) types.push(0) // Vacancies
        
        // Fill remaining spots split 50/50 between Type 1 and Type 2
        const remaining = totalSites - numVacancies
        const numType1 = Math.floor(remaining / 2)
        const numType2 = remaining - numType1
        
        for(let i=0; i<numType1; i++) types.push(1) // Type A (Blue)
        for(let i=0; i<numType2; i++) types.push(2) // Type B (Orange/Red)
        
        // Fisher-Yates Shuffle
        for (let i = types.length - 1; i > 0; i--) {
            const j = Math.floor(Math.random() * (i + 1));
            [types[i], types[j]] = [types[j], types[i]];
        }
        
        // Fill Grid
        let idx = 0
        for(let i=0; i<N; i++) {
            for(let j=0; j<N; j++) {
                newGrid[i][j] = types[idx++]
            }
        }
        
        latticeRef.current = newGrid
        setStats({ step: 0, energy: 0, acceptanceRate: 0, time: 0.0 })
    }
    
    // Clean update
    setTimeout(() => updateChart(), 0)
  }, [simMode, density, particleCount, latticeSize, vacancyPercent])

  // --- ANIMATION LOOP ---
  const animate = () => {
    // Speed control: KMC needs fewer loops per frame to look natural, MMC needs more for convergence
    const loops = simMode === 'metropolis' ? 5 : 1
    
    for(let i=0; i<loops; i++) {
        if (simMode === 'metropolis') performMetropolisSteps()
        else performKMCSteps()
    }

    if (stepsRef.current % 5 === 0) {
        // Recalc energy occasionally
        let currentEnergy = 0
        if (simMode === 'metropolis') currentEnergy = calculateTotalFluidEnergy()
        // For KMC we skip expensive full energy sum every frame, just show 0 or delta
        
        const accRate = totalMovesRef.current > 0 ? acceptedMovesRef.current / totalMovesRef.current : 0
        setStats({
            step: stepsRef.current,
            energy: currentEnergy,
            acceptanceRate: accRate,
            time: timeRef.current
        })
        updateChart()
    }

    if (isRunning) {
      animationFrameRef.current = requestAnimationFrame(animate)
    }
  }

  const updateChart = () => {
    if (echartsRef.current) {
        const instance = echartsRef.current.getEchartsInstance()
        
        // COMMON SETTINGS
        // We must always provide the full series configuration to avoid "undefined series" or clearing data
        
        if (simMode === 'metropolis') {
             const data = particlesRef.current.map(p => [p.x, p.y])
             
             instance.setOption({
                xAxis: { max: boxSizeRef.current },
                yAxis: { max: boxSizeRef.current },
                series: [
                    {
                        name: 'Fluid',
                        type: 'scatter',
                        data: data,
                        symbol: 'circle',
                        symbolSize: 8,
                        itemStyle: { color: isDark ? '#4ade80' : '#16a34a' }
                    }, 
                    // We must explicitely clear the other series if we switch modes
                    { name: 'Type B', type: 'scatter', data: [] } 
                ]
             })
        } else {
             // KINETIC MODE
             const grid = latticeRef.current
             const N = latticeSize
             const dataType1: number[][] = []
             const dataType2: number[][] = []
             
             // Extract positions for both types
             for(let i=0; i<N; i++) {
                 for(let j=0; j<N; j++) {
                     if (grid[i][j] === 1) dataType1.push([i, j])
                     else if (grid[i][j] === 2) dataType2.push([i, j])
                 }
             }

             const size = (500 / boxSizeRef.current) * 0.8

             instance.setOption({
                xAxis: { max: N },
                yAxis: { max: N },
                series: [
                    {
                        name: 'Type A',
                        type: 'scatter',
                        data: dataType1,
                        symbol: 'rect',
                        symbolSize: size,
                        itemStyle: { color: '#3b82f6' } // Blue
                    },
                    {
                        name: 'Type B',
                        type: 'scatter',
                        data: dataType2, // <--- Ensure this is being populated
                        symbol: 'rect',
                        symbolSize: size,
                        itemStyle: { color: '#f97316' } // Orange
                    }
                ]
            })
        }
    }
  }

  // --- HELPER: Format time ---
  const formatTime = (t: number) => {
      if (t < 1e-9) return `${(t*1e12).toFixed(2)} ps`
      if (t < 1e-6) return `${(t*1e9).toFixed(2)} ns`
      if (t < 1e-3) return `${(t*1e6).toFixed(2)} μs`
      if (t < 1.0) return `${(t*1e3).toFixed(2)} ms`
      return `${t.toFixed(4)} s`
  }

  // --- HANDLERS ---
  const handleToggleSimulation = () => {
    if (isRunning) {
      setIsRunning(false)
      if (animationFrameRef.current) cancelAnimationFrame(animationFrameRef.current)
    } else {
      setIsRunning(true)
    }
  }

  useEffect(() => {
    if (isRunning) {
      animationFrameRef.current = requestAnimationFrame(animate)
    }
    return () => {
      if (animationFrameRef.current) cancelAnimationFrame(animationFrameRef.current)
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [isRunning])

  useEffect(() => {
    initializeSystem()
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [simMode]) // Only auto-init on mode switch

  // --- CHART OPTIONS ---
  const getChartOption = (): EChartsOption => {
    const textColor = isDark ? '#eee' : '#333'
    return {
      backgroundColor: 'transparent',
      title: {
        text: simMode === 'metropolis' ? '2D Fluid (Off-Lattice)' : 'Crystal Diffusion (On-Lattice)',
        left: 'center',
        textStyle: { color: textColor, fontSize: 16 }
      },
      grid: { left: '5%', right: '5%', top: '10%', bottom: '10%', containLabel: true },
      tooltip: { show: false },
      xAxis: {
        type: 'value',
        min: 0,
        splitLine: { show: false },
        axisLabel: { show: false },
        axisTick: { show: false }
      },
      yAxis: {
        type: 'value',
        min: 0,
        splitLine: { show: false },
        axisLabel: { show: false },
        axisTick: { show: false }
      },
      animation: false,
      series: [
        {
          type: 'scatter',
          data: []
        }
      ]
    }
  }

  return (
    <div className="container mx-auto p-4 md:p-8 px-4 md:px-16">
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        
        {/* LEFT COLUMN: Controls */}
        <div className="lg:col-span-1 space-y-6">
          <Card>
            <CardContent className="space-y-6 p-4">
              
              {/* MODE SWITCHER */}
              <div className="w-full">
                  <Tabs value={simMode} onValueChange={(v) => {
                      setSimMode(v as SimMode)
                      setIsRunning(false)
                  }}>
                    <TabsList className="grid w-full grid-cols-2">
                        <TabsTrigger value="metropolis">Metropolis (Fluid)</TabsTrigger>
                        <TabsTrigger value="kinetic">Kinetic (Lattice)</TabsTrigger>
                    </TabsList>
                  </Tabs>
              </div>

              <div className="space-y-4">
                 <div className="flex items-center justify-between">
                    <h3 className="font-bold text-lg">Simulation Controls</h3>
                    <Badge variant={isRunning ? "default" : "secondary"}>
                        {isRunning ? "Running" : "Paused"}
                    </Badge>
                 </div>
                 
                 <Button 
                    onClick={handleToggleSimulation} 
                    className={`w-full ${isRunning ? 'bg-amber-600 hover:bg-amber-700' : ''}`}
                 >
                    {isRunning ? 'Pause' : 'Start Simulation'}
                 </Button>

                 <Button variant="outline" onClick={initializeSystem} className="w-full">
                    Reset / Apply Settings
                 </Button>
              </div>

              <div className="space-y-6 pt-4 border-t">
                {/* Temperature Slider (Common) */}
                <div className="space-y-3">
                   <div className="flex justify-between">
                     <Label>Temperature (T*)</Label>
                     <span className="text-sm font-mono">{temperature.toFixed(2)}</span>
                   </div>
                   <Slider 
                     min={0.1} max={5.0} step={0.1} 
                     value={[temperature]} 
                     onValueChange={([v]) => setTemperature(v)}
                   />
                   <p className="text-xs text-muted-foreground">
                     {simMode === 'metropolis' ? 'High T = Gas, Low T = Liquid.' : 'High T = Fast Diffusion.'}
                   </p>
                </div>

                {/* Conditional Inputs based on Mode */}
                {simMode === 'metropolis' ? (
                    <>
                        <div className="space-y-3">
                           <div className="flex justify-between">
                             <Label>Density (&rho;*)</Label>
                             <span className="text-sm font-mono">{density.toFixed(2)}</span>
                           </div>
                           <Slider 
                             min={0.05} max={0.8} step={0.05} 
                             value={[density]} 
                             onValueChange={([v]) => setDensity(v)} 
                             disabled={isRunning}
                           />
                        </div>
                        <div className="space-y-3">
                           <Label>Number of Particles (N)</Label>
                           <Input 
                                type="number" 
                                value={particleCount}
                                onChange={(e) => setParticleCount(Number(e.target.value))}
                                disabled={isRunning}
                           />
                        </div>
                    </>
                ) : (
                    <>
                        <div className="space-y-3">
                           <div className="flex justify-between">
                             <Label>Activation Barrier (E<sub>a</sub>)</Label>
                             <span className="text-sm font-mono">{activationEnergy.toFixed(2)} eV</span>
                           </div>
                           <Slider 
                             min={0.1} max={1.0} step={0.05} 
                             value={[activationEnergy]} 
                             onValueChange={([v]) => setActivationEnergy(v)} 
                           />
                           <p className="text-xs text-muted-foreground">
                             Higher barrier = Slower diffusion (Time increases faster).
                           </p>
                        </div>

                        <div className="space-y-3">
                           <div className="flex justify-between">
                             <Label>Vacancy Percent</Label>
                             <span className="text-sm font-mono">{(vacancyPercent * 100).toFixed(0)}%</span>
                           </div>
                           <Slider 
                             min={0.01} max={0.9} step={0.01} 
                             value={[vacancyPercent]} 
                             onValueChange={([v]) => setVacancyPercent(v)} 
                             disabled={isRunning}
                           />
                           <p className="text-xs text-muted-foreground">
                               Empty spots required for hopping.
                           </p>
                        </div>
                        <div className="space-y-3">
                           <Label>Lattice Size (NxN)</Label>
                           <Input 
                                type="number" 
                                value={latticeSize}
                                onChange={(e) => setLatticeSize(Number(e.target.value))}
                                disabled={isRunning}
                           />
                        </div>
                    </>
                )}
              </div>
            </CardContent>
          </Card>

          {/* Statistics Card */}
          <Card>
            <CardContent className="p-4 space-y-4">
               <h3 className="font-bold text-md">Real-time Statistics</h3>
               
               <div className="grid grid-cols-2 gap-4">
                  <div className="bg-muted p-3 rounded-md text-center">
                    <div className="text-xs text-muted-foreground uppercase">Steps</div>
                    <div className="text-xl font-bold font-mono">{stats.step}</div>
                  </div>
                  <div className="bg-muted p-3 rounded-md text-center">
                    <div className="text-xs text-muted-foreground uppercase">Acceptance</div>
                    <div className="text-xl font-bold font-mono">{(stats.acceptanceRate * 100).toFixed(1)}%</div>
                  </div>
                  {simMode === 'kinetic' && (
                    <div className="col-span-2 bg-muted p-3 rounded-md text-center">
                      <div className="text-xs text-muted-foreground uppercase">Sim Time</div>
                      <div className="text-xl font-bold font-mono">{formatTime(stats.time)}</div>
                    </div>
                  )}
               </div>
            </CardContent>
          </Card>
        </div>

        {/* RIGHT COLUMN: Visualization */}
        <div className="lg:col-span-2">
            <Card className="h-full min-h-[500px]">
              <CardContent className="p-2 h-full flex flex-col">
                <div className="relative flex-1 w-full h-full min-h-[500px] rounded-md overflow-hidden bg-background">
                   <ReactECharts 
                      ref={echartsRef}
                      echarts={echarts}
                      option={getChartOption()}
                      style={{ height: '100%', width: '100%' }}
                      notMerge={true} 
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