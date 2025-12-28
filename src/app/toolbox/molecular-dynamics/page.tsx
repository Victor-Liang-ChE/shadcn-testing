'use client'

import React, { useState, useEffect, useCallback, useRef, useMemo } from 'react'
import { useTheme } from "next-themes"

// ECharts imports
import ReactECharts from 'echarts-for-react'
import * as echarts from 'echarts/core'
import type { EChartsOption } from 'echarts'
import { LineChart, ScatterChart } from 'echarts/charts'
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  VisualMapComponent
} from 'echarts/components'
import { CanvasRenderer } from 'echarts/renderers'

// React Three Fiber imports for 3D visualization
import { Canvas, useFrame } from "@react-three/fiber"
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
  MousePointer2, 
  Activity, 
  List, 
  Zap,
  Info,
  AlertCircle
} from 'lucide-react'

// Type declaration for react-katex
declare module 'react-katex' {
  export function BlockMath(props: { math: string }): React.ReactElement
}

// Types for our simulation
type Particle = {
  id: number
  type: 0 | 1
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

// --- UNIT CONSTANTS & CONVERSIONS ---
// System: Energy=Kelvin, Dist=Angstrom, Mass=Dalton, Time=ps
// Boltzmann constant is 1.0 in internal units (E/kB).
// The real conversion factors are needed for derived properties.
//
// 1. Acceleration Factor (a = F/m conversion)
// Acceleration = Force / Mass
// Units: Å² / (ps²·K)
// Conversion: 1 K/(Å·Da) ≈ 0.831446 Å²/(ps²·K)
const KB_OVER_MU = 0.831446

// 2. Pressure Factor
// P = Energy/Volume. Units: K/Å³.
// 1 K/Å³ = 138.06 bar
const PRESSURE_CONV_BAR = 138.06

const MASS_ARGON = 39.948 // Daltons

// Simulation Constants
const BOX_SIZE = 40.0 // Angstroms (typical nano-box for Argon)
const DT_DEFAULT = 0.005 // 5 femtoseconds (0.005 ps)

// RDF Configuration
const RDF_BINS = 100
const RDF_CUTOFF = BOX_SIZE / 2.0
const RDF_BIN_WIDTH = RDF_CUTOFF / RDF_BINS

// History caps (prevents unbounded memory growth)
const MAX_HISTORY_POINTS = 12000

type LJReference = {
    formula: string
    substance: string
    sigmaAngstrom: number
    epsilonOverKbK: number
}

// Reference Lennard-Jones parameters (σ in Å, ε/kB in K).
// Used only for UI hints to help users recognize common fluids.
const LJ_REFERENCES: LJReference[] = [
    { formula: 'Ar', substance: 'Argon', sigmaAngstrom: 3.542, epsilonOverKbK: 93.3 },
    { formula: 'He', substance: 'Helium', sigmaAngstrom: 2.551, epsilonOverKbK: 10.22 },
    { formula: 'Kr', substance: 'Krypton', sigmaAngstrom: 3.655, epsilonOverKbK: 178.9 },
    { formula: 'Ne', substance: 'Neon', sigmaAngstrom: 2.82, epsilonOverKbK: 32.8 },
    { formula: 'Xe', substance: 'Xenon', sigmaAngstrom: 4.047, epsilonOverKbK: 231 },
    { formula: 'Air', substance: 'Air', sigmaAngstrom: 3.711, epsilonOverKbK: 78.6 },
    { formula: 'AsH₃', substance: 'Arsine', sigmaAngstrom: 4.145, epsilonOverKbK: 259.8 },
    { formula: 'BCl₃', substance: 'Boron chloride', sigmaAngstrom: 5.127, epsilonOverKbK: 337.7 },
    { formula: 'BF₃', substance: 'Boron floride', sigmaAngstrom: 4.198, epsilonOverKbK: 186.3 },
    { formula: 'B(OCH₃)₃', substance: 'Methyl borate', sigmaAngstrom: 5.503, epsilonOverKbK: 396.7 },
    { formula: 'Br₂', substance: 'Bromine', sigmaAngstrom: 4.296, epsilonOverKbK: 507.9 },
    { formula: 'CCl₄', substance: 'Carbon tetrachloride', sigmaAngstrom: 5.947, epsilonOverKbK: 322.7 },
    { formula: 'CF₄', substance: 'Carbon tetrafluoride', sigmaAngstrom: 4.662, epsilonOverKbK: 134 },
    { formula: 'CHCl₃', substance: 'Chloroform', sigmaAngstrom: 5.389, epsilonOverKbK: 340.2 },
    { formula: 'CH₂Cl₂', substance: 'Methylene chloride', sigmaAngstrom: 4.898, epsilonOverKbK: 356.3 },
    { formula: 'CH₃Br', substance: 'Methyl bromide', sigmaAngstrom: 4.118, epsilonOverKbK: 449.2 },
    { formula: 'CH₃Cl', substance: 'Methyl chloride', sigmaAngstrom: 4.182, epsilonOverKbK: 350 },
    { formula: 'CH₃OH', substance: 'Methanol', sigmaAngstrom: 3.626, epsilonOverKbK: 481.8 },
    { formula: 'CH₄', substance: 'Methane', sigmaAngstrom: 3.758, epsilonOverKbK: 148.6 },
    { formula: 'CO', substance: 'Carbon monoxide', sigmaAngstrom: 3.69, epsilonOverKbK: 91.7 },
    { formula: 'COS', substance: 'Carbonyl sulfide', sigmaAngstrom: 4.13, epsilonOverKbK: 336 },
    { formula: 'CO₂', substance: 'Carbon dioxide', sigmaAngstrom: 3.941, epsilonOverKbK: 195.2 },
    { formula: 'CS₂', substance: 'Carbon disulfide', sigmaAngstrom: 4.483, epsilonOverKbK: 467 },
    { formula: 'C₂H₂', substance: 'Acetylene', sigmaAngstrom: 4.033, epsilonOverKbK: 231.8 },
    { formula: 'C₂H₄', substance: 'Ethylene', sigmaAngstrom: 4.163, epsilonOverKbK: 224.7 },
    { formula: 'C₂H₆', substance: 'Ethane', sigmaAngstrom: 4.443, epsilonOverKbK: 215.7 },
    { formula: 'C₂H₅Cl', substance: 'Ethyl chloride', sigmaAngstrom: 4.898, epsilonOverKbK: 300 },
    { formula: 'C₂H₅OH', substance: 'Ethanol', sigmaAngstrom: 4.53, epsilonOverKbK: 362.6 },
    { formula: 'C₂N₂', substance: 'Cyanogen', sigmaAngstrom: 4.361, epsilonOverKbK: 348.6 },
    { formula: 'CH₃OCH₃', substance: 'Methyl ether', sigmaAngstrom: 4.307, epsilonOverKbK: 395 },
    { formula: 'CH₂CHCH₃', substance: 'Propylene', sigmaAngstrom: 4.678, epsilonOverKbK: 298.9 },
    { formula: 'CH₃CCH', substance: 'Methylacetylene', sigmaAngstrom: 4.761, epsilonOverKbK: 251.8 },
    { formula: 'C₃H₆', substance: 'Cyclopropane', sigmaAngstrom: 4.807, epsilonOverKbK: 248.9 },
    { formula: 'C₃H₈', substance: 'Propane', sigmaAngstrom: 5.118, epsilonOverKbK: 237.1 },
    { formula: 'n-C₃H₇OH', substance: 'n-Propyl alcohol', sigmaAngstrom: 4.549, epsilonOverKbK: 576.7 },
    { formula: 'CH₃COCH₃', substance: 'Acetone', sigmaAngstrom: 4.6, epsilonOverKbK: 560.2 },
    { formula: 'CH₃COOCH₃', substance: 'Methyl acetate', sigmaAngstrom: 4.936, epsilonOverKbK: 469.8 },
    { formula: 'n-C₄H₁₀', substance: 'n-Butane', sigmaAngstrom: 4.687, epsilonOverKbK: 531.4 },
    { formula: 'iso-C₄H₁₀', substance: 'Isobutane', sigmaAngstrom: 5.278, epsilonOverKbK: 330.1 },
    { formula: 'C₂H₅OC₂H₅', substance: 'Ethyl ether', sigmaAngstrom: 5.678, epsilonOverKbK: 313.8 },
    { formula: 'CH₃COOC₂H₅', substance: 'Ethyl acetate', sigmaAngstrom: 5.205, epsilonOverKbK: 521.3 },
    { formula: 'n-C₅H₁₂', substance: 'n-Pentane', sigmaAngstrom: 5.784, epsilonOverKbK: 341.1 },
    { formula: 'C(CH₃)₄', substance: '2,2-Dimethylpropane', sigmaAngstrom: 6.464, epsilonOverKbK: 193.4 },
    { formula: 'C₆H₆', substance: 'Benzene', sigmaAngstrom: 5.349, epsilonOverKbK: 412.3 },
    { formula: 'C₆H₁₂', substance: 'Cyclohexane', sigmaAngstrom: 6.182, epsilonOverKbK: 297.1 },
    { formula: 'n-C₆H₁₄', substance: 'n-Hexane', sigmaAngstrom: 5.949, epsilonOverKbK: 399.3 },
    { formula: 'Cl₂', substance: 'Chlorine', sigmaAngstrom: 4.217, epsilonOverKbK: 316 },
    { formula: 'F₂', substance: 'Fluorine', sigmaAngstrom: 3.357, epsilonOverKbK: 112.6 },
    { formula: 'HBr', substance: 'Hydrogen bromide', sigmaAngstrom: 3.353, epsilonOverKbK: 449 },
    { formula: 'HCN', substance: 'Hydrogen cyanide', sigmaAngstrom: 3.63, epsilonOverKbK: 569.1 },
    { formula: 'HCl', substance: 'Hydrogen chloride', sigmaAngstrom: 3.339, epsilonOverKbK: 344.7 },
    { formula: 'HF', substance: 'Hydrogen fluoride', sigmaAngstrom: 3.148, epsilonOverKbK: 330 },
    { formula: 'HI', substance: 'Hydrogen iodide', sigmaAngstrom: 4.211, epsilonOverKbK: 288.7 },
    { formula: 'H₂', substance: 'Hydrogen', sigmaAngstrom: 2.827, epsilonOverKbK: 59.7 },
    { formula: 'H₂O', substance: 'Water', sigmaAngstrom: 2.641, epsilonOverKbK: 809.1 },
    { formula: 'H₂O₂', substance: 'Hydrogen peroxide', sigmaAngstrom: 4.196, epsilonOverKbK: 289.3 },
    { formula: 'H₂S', substance: 'Hydrogen sulfide', sigmaAngstrom: 3.623, epsilonOverKbK: 301.1 },
    { formula: 'Hg', substance: 'Mercury', sigmaAngstrom: 2.969, epsilonOverKbK: 750 },
    { formula: 'HgBr₂', substance: 'Mercuric bromide', sigmaAngstrom: 5.08, epsilonOverKbK: 686.2 },
    { formula: 'HgCl₂', substance: 'Mercuric chloride', sigmaAngstrom: 4.55, epsilonOverKbK: 750 },
    { formula: 'HgI₂', substance: 'Mercuric iodide', sigmaAngstrom: 5.625, epsilonOverKbK: 695.6 },
    { formula: 'I₂', substance: 'Iodine', sigmaAngstrom: 5.16, epsilonOverKbK: 474.2 },
    { formula: 'NH₃', substance: 'Ammonia', sigmaAngstrom: 2.9, epsilonOverKbK: 558.3 },
    { formula: 'NO', substance: 'Nitric oxide', sigmaAngstrom: 3.492, epsilonOverKbK: 116.7 },
    { formula: 'NOCl', substance: 'Nitrosyl chloride', sigmaAngstrom: 4.112, epsilonOverKbK: 395.3 },
    { formula: 'N₂', substance: 'Nitrogen', sigmaAngstrom: 3.798, epsilonOverKbK: 71.4 },
    { formula: 'N₂O', substance: 'Nitrous oxide', sigmaAngstrom: 3.828, epsilonOverKbK: 232.4 },
    { formula: 'O₂', substance: 'Oxygen', sigmaAngstrom: 3.467, epsilonOverKbK: 106.7 },
    { formula: 'PH₃', substance: 'Phosphine', sigmaAngstrom: 3.981, epsilonOverKbK: 251.5 },
    { formula: 'SF₆', substance: 'Sulfur hexafluoride', sigmaAngstrom: 5.128, epsilonOverKbK: 222.1 },
    { formula: 'SO₂', substance: 'Sulfur dioxide', sigmaAngstrom: 4.112, epsilonOverKbK: 335.4 },
    { formula: 'SiF₄', substance: 'Silicon tetrafluoride', sigmaAngstrom: 4.88, epsilonOverKbK: 171.9 },
    { formula: 'SiH₄', substance: 'Silicon hydride', sigmaAngstrom: 4.084, epsilonOverKbK: 207.6 },
    { formula: 'SnBr₄', substance: 'Stannic bromide', sigmaAngstrom: 6.388, epsilonOverKbK: 563.7 },
    { formula: 'UF₆', substance: 'Uranium hexafluoride', sigmaAngstrom: 5.967, epsilonOverKbK: 236.8 }
]

const LJ_EPSILON_STEP_K = 5
const LJ_SIGMA_STEP_A = 0.01

const roundToStep = (value: number, step: number) => Math.round(value / step) * step
const nearlyEqual = (a: number, b: number, eps = 1e-9) => Math.abs(a - b) <= eps

const getLJMatchesForPair = (epsilonK: number, sigmaA: number): LJReference[] => {
    if (!Number.isFinite(epsilonK) || !Number.isFinite(sigmaA)) return []

    const epsRounded = roundToStep(epsilonK, LJ_EPSILON_STEP_K)
    const sigRounded = roundToStep(sigmaA, LJ_SIGMA_STEP_A)

    return LJ_REFERENCES.filter((ref) => {
        const refEpsRounded = roundToStep(ref.epsilonOverKbK, LJ_EPSILON_STEP_K)
        const refSigRounded = roundToStep(ref.sigmaAngstrom, LJ_SIGMA_STEP_A)
        return nearlyEqual(refEpsRounded, epsRounded) && nearlyEqual(refSigRounded, sigRounded)
    })
}

const getLJMatchesForValue = (
    value: number,
    key: 'sigmaAngstrom' | 'epsilonOverKbK'
): LJReference[] => {
    if (!Number.isFinite(value)) return []
    // Keep this slightly looser than the slider step so users can still see nearby references.
    const tolerance = key === 'epsilonOverKbK' ? 5 : 0.15 // ±5 K for epsilon, ±0.15 Å for sigma
    const matches = LJ_REFERENCES.filter((ref) => Math.abs(ref[key] - value) <= tolerance)
    matches.sort((a, b) => Math.abs(a[key] - value) - Math.abs(b[key] - value))
    return matches
}

const POLAR_LJ_FORMULAS = new Set(['H₂O', 'NH₃', 'HF'])
const isPolarLJReference = (ref: LJReference): boolean => {
    if (POLAR_LJ_FORMULAS.has(ref.formula)) return true
    const substance = ref.substance.toLowerCase()
    // Alcohols (and a couple of common explicit names) are hydrogen-bonding / polar.
    return substance.includes('alcohol') || substance === 'methanol' || substance === 'ethanol'
}

export default function MolecularDynamicsPage() {
  const { resolvedTheme } = useTheme()
    const isDarkTheme = resolvedTheme === 'dark'
    const velColor = isDarkTheme ? '#4ade80' : '#16a34a'
    const forceColor = isDarkTheme ? '#facc15' : '#d97706'
    const [mounted, setMounted] = useState(false)
    useEffect(() => { setMounted(true) }, [])
  
  // --- Simulation State ---
  const [running, setRunning] = useState(false)
  const [isCrashed, setIsCrashed] = useState(false)
  const [particles, setParticles] = useState<Particle[]>([])
  const [selectedId, setSelectedId] = useState<number | null>(null)
  const [isThreeD, setIsThreeD] = useState(false)
  
  // --- Binary Mixture State ---
  const [isBinary, setIsBinary] = useState(false)
  
  // Particle Counts
  const [numParticles, setNumParticles] = useState<number>(60)
  const [numParticlesA, setNumParticlesA] = useState(40)
  const [numParticlesB, setNumParticlesB] = useState(40)
  
  // --- Control States (Real Units for Argon) ---
  const [temperature, setTemperature] = useState<number>(120.0) // Kelvin (Liquid Argon region)
  const [temperatureUnit, setTemperatureUnit] = useState<'K' | 'C' | 'F'>('K')
  const [timeStep, setTimeStep] = useState<number>(DT_DEFAULT)
  
  // Parameters for Species A (Default Argon)
  const [epsilonA, setEpsilonA] = useState<number>(120.0) // Kelvin (eps/kB)
  const [sigmaA, setSigmaA] = useState<number>(3.4) // Angstroms

  // Parameters for Species B (Default Krypton-ish)
  const [epsilonB, setEpsilonB] = useState<number>(160.0)
  const [sigmaB, setSigmaB] = useState<number>(4.0)
  
  // Keep backward-compat alias for existing code
  const epsilon = epsilonA
  const setEpsilon = setEpsilonA
  const sigma = sigmaA
  const setSigma = setSigmaA

    // Use reduced temperature T* = T/ε as the primary UI control.
    // Note: since ε is stored in Kelvin (i.e., ε/kB), reduced temperature is simply T/ε.
    const REDUCED_STEP = 0.01
    const snapToStep = (v: number, step: number) => Math.round(v / step) * step
    const [reducedTemp, setReducedTemp] = useState<number>(() => snapToStep(temperature / epsilon, REDUCED_STEP))

    // Compute slider-friendly min/max for reduced temperature so endpoints land on "nice" multiples
    const reducedMinRaw = 10 / epsilon
    const reducedMaxRaw = 500 / epsilon
    let reducedMin = Math.ceil(reducedMinRaw / REDUCED_STEP) * REDUCED_STEP
    let reducedMax = Math.floor(reducedMaxRaw / REDUCED_STEP) * REDUCED_STEP
    if (reducedMin > reducedMax) {
        // Fallback: if epsilon is such that no multiple exists within range, use raw values
        reducedMin = reducedMinRaw
        reducedMax = reducedMaxRaw
    }

    // Keep reducedTemp in-range and snapped when epsilon changes
    useEffect(() => {
        setReducedTemp((prev) => {
            const snapped = snapToStep(prev, REDUCED_STEP)
            const clamped = Math.min(reducedMax, Math.max(reducedMin, snapped))
            return clamped
        })
    }, [epsilon])

    // Keep the actual target temperature consistent with the reduced temperature.
    useEffect(() => {
        setTemperature(Number((reducedTemp * epsilon).toFixed(3)))
    }, [reducedTemp, epsilon])
  
  // Potential Type and Generic Parameter C
    type PotentialType = 'LJ' | 'WCA' | 'MORSE' | 'SOFT' | 'YUKAWA' | 'GAUSS' | 'BUCK'
  const [potentialType, setPotentialType] = useState<PotentialType>('LJ')
  const [paramC, setParamC] = useState<number>(1.5) // Width for Morse, Exponent for Soft
  const [useGravity, setUseGravity] = useState(false)
  const [gravityStrength, setGravityStrength] = useState<number>(5.0)

    const [showGravityNotice, setShowGravityNotice] = useState(false)
    const [gravityNoticeFading, setGravityNoticeFading] = useState(false)
    const gravityNoticeTimersRef = useRef<{ fade: number | null; hide: number | null }>({ fade: null, hide: null })

    const showLJReferenceHints = potentialType === 'LJ'

    const ljEpsilonMatchesA = useMemo(
        () => (showLJReferenceHints ? getLJMatchesForValue(epsilonA, 'epsilonOverKbK') : []),
        [showLJReferenceHints, epsilonA]
    )
    const ljSigmaMatchesA = useMemo(
        () => (showLJReferenceHints ? getLJMatchesForValue(sigmaA, 'sigmaAngstrom') : []),
        [showLJReferenceHints, sigmaA]
    )
    const ljEpsilonMatchesB = useMemo(
        () => (showLJReferenceHints ? getLJMatchesForValue(epsilonB, 'epsilonOverKbK') : []),
        [showLJReferenceHints, epsilonB]
    )
    const ljSigmaMatchesB = useMemo(
        () => (showLJReferenceHints ? getLJMatchesForValue(sigmaB, 'sigmaAngstrom') : []),
        [showLJReferenceHints, sigmaB]
    )

    const ljRefId = (ref: LJReference) => `${ref.substance}::${ref.formula}`

    // Intersection-based “pair match”: whatever appears in BOTH the epsilon list and the sigma list.
    const ljPairMatchesA = useMemo(() => {
        if (!showLJReferenceHints) return []
        const epsilonIds = new Set(ljEpsilonMatchesA.map(ljRefId))
        return ljSigmaMatchesA.filter((ref) => epsilonIds.has(ljRefId(ref)))
    }, [showLJReferenceHints, ljEpsilonMatchesA, ljSigmaMatchesA])

    const ljPairMatchesB = useMemo(() => {
        if (!showLJReferenceHints) return []
        const epsilonIds = new Set(ljEpsilonMatchesB.map(ljRefId))
        return ljSigmaMatchesB.filter((ref) => epsilonIds.has(ljRefId(ref)))
    }, [showLJReferenceHints, ljEpsilonMatchesB, ljSigmaMatchesB])

    const polarLJRefsShown = useMemo(() => {
        if (!showLJReferenceHints) return [] as LJReference[]
        const refs: LJReference[] = [
            ...ljEpsilonMatchesA,
            ...ljSigmaMatchesA,
            ...ljPairMatchesA,
            ...(isBinary ? [...ljEpsilonMatchesB, ...ljSigmaMatchesB, ...ljPairMatchesB] : [])
        ]

        const seen = new Set<string>()
        const unique: LJReference[] = []
        for (const ref of refs) {
            const id = ljRefId(ref)
            if (seen.has(id)) continue
            seen.add(id)
            unique.push(ref)
        }

        return unique.filter(isPolarLJReference)
    }, [
        showLJReferenceHints,
        ljEpsilonMatchesA,
        ljSigmaMatchesA,
        ljPairMatchesA,
        isBinary,
        ljEpsilonMatchesB,
        ljSigmaMatchesB,
        ljPairMatchesB
    ])

    const showPolarLJWarning = showLJReferenceHints && polarLJRefsShown.length > 0
    const polarLJNames = useMemo(() => {
        if (!showPolarLJWarning) return ''
        const names = polarLJRefsShown.map((r) => r.substance)
        // Keep it short so it doesn't dominate the UI.
        const uniqueNames = Array.from(new Set(names))
        return uniqueNames.slice(0, 6).join(', ') + (uniqueNames.length > 6 ? ', …' : '')
    }, [showPolarLJWarning, polarLJRefsShown])

    useEffect(() => {
        // Show a short-lived notice when gravity is enabled (because it changes Y boundary conditions).
        if (!useGravity) return

        const timers = gravityNoticeTimersRef.current
        if (timers.fade !== null) window.clearTimeout(timers.fade)
        if (timers.hide !== null) window.clearTimeout(timers.hide)

        setGravityNoticeFading(false)
        setShowGravityNotice(true)

        timers.fade = window.setTimeout(() => setGravityNoticeFading(true), 2500)
        timers.hide = window.setTimeout(() => {
            setShowGravityNotice(false)
            setGravityNoticeFading(false)
        }, 3200)

        return () => {
            const t = gravityNoticeTimersRef.current
            if (t.fade !== null) window.clearTimeout(t.fade)
            if (t.hide !== null) window.clearTimeout(t.hide)
            t.fade = null
            t.hide = null
        }
    }, [useGravity])

    // Recalculate forces when gravity settings change (even when paused/stopped) so force arrows update
    useEffect(() => {
        if (!running && particlesRef.current.length > 0) {
            calculateForces(particlesRef.current)
            setParticles([...particlesRef.current])
        }
    }, [useGravity, gravityStrength, running])

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
  
  // Refs for current parameter values
  const paramsRef = useRef({
      dt: DT_DEFAULT,
      epsilonA: 120.0,
      sigmaA: 3.4,
      epsilonB: 160.0,
      sigmaB: 4.0,
      paramC: 1.5,
      potentialType: 'LJ' as PotentialType,
      box: BOX_SIZE,
      targetTemp: 120.0,
      isBinary: false,
      useGravity: false,
      gravityStrength: 5.0
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
      paramsRef.current.epsilonA = epsilonA
      paramsRef.current.sigmaA = sigmaA
      paramsRef.current.epsilonB = epsilonB
      paramsRef.current.sigmaB = sigmaB
      paramsRef.current.paramC = paramC
      paramsRef.current.potentialType = potentialType
      paramsRef.current.targetTemp = temperature
      paramsRef.current.isBinary = isBinary
      paramsRef.current.useGravity = useGravity
      paramsRef.current.gravityStrength = gravityStrength
      
      // Update Mass/Q if needed (assuming same mass for now to keep simple)
      const totalN = isBinary ? (numParticlesA + numParticlesB) : numParticles
      nhStateRef.current.Q = totalN * temperature * 0.25
  }, [timeStep, epsilonA, sigmaA, epsilonB, sigmaB, paramC, potentialType, temperature, numParticles, numParticlesA, numParticlesB, isBinary, useGravity, gravityStrength])

  // --- PHYSICS KERNEL: Force Calculation (3D O(N^2)) with Mixing Rules ---
  const calculateForces = (currentParticles: Particle[]): { pe: number, virial: number } => {
      const { epsilonA, sigmaA, epsilonB, sigmaB, box, isBinary, useGravity } = paramsRef.current
      
      let potentialEnergy = 0
      let virialSum = 0
      
      // Reset Forces
      for (let p of currentParticles) {
          p.fx = 0; p.fy = 0; p.fz = 0
          
          // Add Gravity Force (External Field)
          if (useGravity) {
              // Apply a downward force (negative Y direction)
              p.fy -= paramsRef.current.gravityStrength
          }
      }

      const N = currentParticles.length
      
      // Determine max sigma for cutoff (use larger of the two in binary mode)
      const maxSigma = isBinary ? Math.max(sigmaA, sigmaB) : sigmaA
      const cutoffSq = 6.25 * maxSigma * maxSigma
      
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
              
              // --- SELECT PARAMETERS via Lorentz-Berthelot Mixing Rules ---
              let eps, sig
              
              if (!isBinary) {
                  // Simple Case: Pure substance
                  eps = epsilonA
                  sig = sigmaA
              } else {
                  // Binary Case: Apply mixing rules for cross interactions
                  if (p1.type === 0 && p2.type === 0) {
                      eps = epsilonA
                      sig = sigmaA
                  } else if (p1.type === 1 && p2.type === 1) {
                      eps = epsilonB
                      sig = sigmaB
                  } else {
                      // Cross interaction: σ_mix = (σ_A + σ_B)/2, ε_mix = sqrt(ε_A * ε_B)
                      sig = 0.5 * (sigmaA + sigmaB)
                      eps = Math.sqrt(epsilonA * epsilonB)
                  }
              }
              
              // Pre-compute LJ shift at cutoff for this pair's parameters
              const rc = 2.5 * sig
              const rc_inv = 1.0 / rc
              const rc2_inv = rc_inv * rc_inv
              const s2 = sig * sig
              const s2_rc2 = s2 * rc2_inv
              const s6_rc6 = s2_rc2 * s2_rc2 * s2_rc2
              const s12_rc12 = s6_rc6 * s6_rc6
              const ljShiftAtCutoff = 4 * eps * (s12_rc12 - s6_rc6)
              
              let forceScalar = 0 
              let pairPotential = 0

              switch (paramsRef.current.potentialType) {
                  case 'LJ': {
                      const s2_r2 = s2 * r2_inv
                      const s6_r6 = s2_r2 * s2_r2 * s2_r2
                      const s12_r12 = s6_r6 * s6_r6
                      
                      forceScalar = (24 * eps * r2_inv) * (2 * s12_r12 - s6_r6)
                      pairPotential = 4 * eps * (s12_r12 - s6_r6) - ljShiftAtCutoff
                      break
                  }
                  case 'WCA': {
                      const wcaCutoffSq = 1.25992 * sig * sig 
                      if (distSq < wcaCutoffSq) {
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
                      // V(r) = eps * exp(-r/sig) / r
                      // F(r) = -dV/dr
                      // Using product rule: d/dr (u*v) where u = exp(-r/sig) and v = 1/r
                      // d/dr[exp(-r/sig)] = -(1/sig) * exp(-r/sig)
                      // d/dr[1/r] = -1/r²
                      // Result: F(r) = V(r) * (1/r + 1/sig)
                      const screened = Math.exp(-r / sig)
                      pairPotential = (eps * screened) * r_inv
                      const f_r_yuk = pairPotential * (r_inv + (1.0 / sig))
                      forceScalar = f_r_yuk * r_inv
                      break
                  }
                  case 'GAUSS': {
                      // V(r) = eps * exp(-r²/2σ²)
                      // F(r) = -dV/dr = (eps * r / σ²) * exp(...)
                      // Force vector form: F_vec = F(r) * (r_vec / r) = (eps / σ²) * exp(...) * r_vec
                      // code uses: fx = dx * forceScalar, so forceScalar = (eps / σ²) * exp(...)
                      // No extra r_inv is needed here.
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
              let dz = currentParticles[j].z - currentParticles[i].z

              // MIC (Minimum Image Convention)
              dx -= box * Math.round(dx / box)
              dy -= box * Math.round(dy / box)

              if (isThreeDRef.current) {
                  dz -= box * Math.round(dz / box)
              } else {
                  dz = 0
              }

              const r = Math.sqrt(dx*dx + dy*dy + dz*dz)

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
      const volume = isThreeDRef.current ? Math.pow(box, 3) : Math.pow(box, 2)
      const density = N / volume
      
      return hist.map((count, i) => {
          const r = (i + 1) * RDF_BIN_WIDTH
          // Shell measure depends on dimensionality
          // 2D: dA = 2π r dr
          // 3D: dV = 4π r^2 dr
          const shell = isThreeDRef.current
              ? 4 * Math.PI * r * r * RDF_BIN_WIDTH
              : 2 * Math.PI * r * RDF_BIN_WIDTH
          
          // Expected count for ideal gas = rho * shell
          const idealCount = density * shell
          
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
    const randomNormal = useCallback((): number => {
        // Box–Muller transform
        let u = 0
        let v = 0
        while (u === 0) u = Math.random()
        while (v === 0) v = Math.random()
        return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v)
    }, [])

  const initializeParticles = useCallback(() => {
    const newParticles: Particle[] = []
    
    // Determine Total Count
    const totalParticles = isBinary ? (numParticlesA + numParticlesB) : numParticles
    
    // 3D Logic: Cube Root vs Square Root
    // If 125 particles -> 5x5x5 grid
    const dim = isThreeD 
      ? Math.ceil(Math.pow(totalParticles, 1/3)) 
      : Math.ceil(Math.sqrt(totalParticles))
      
    const spacing = BOX_SIZE / dim

    // 1. Grid Placement (3D or 2D)
    for (let i = 0; i < totalParticles; i++) {
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

      // Determine Type: First N_A are type 0, rest are type 1
      let pType: 0 | 1 = 0
      if (isBinary) {
          pType = i < numParticlesA ? 0 : 1
      }

      newParticles.push({
        id: i,
        type: pType,
        x, y, z,
                vx: 0,
                vy: 0,
                vz: 0,
        fx: 0, fy: 0, fz: 0
      })
    }

// 3. Assign Maxwell–Boltzmann Velocities (Gaussian with std = sqrt(kT/m))
    // v_std = sqrt(T * KB_OVER_MU / m) where T is in Kelvin, m in Daltons
    // This gives velocities in Å/ps
    const vStd = Math.sqrt((temperature * KB_OVER_MU) / MASS_ARGON)
        for (const p of newParticles) {
            p.vx = randomNormal() * vStd
            p.vy = randomNormal() * vStd
            p.vz = isThreeD ? randomNormal() * vStd : 0
        }

    // 4. Remove Center of Mass Motion (Optional but good practice)
    let vcmX = 0, vcmY = 0, vcmZ = 0
    for (let p of newParticles) { vcmX += p.vx; vcmY += p.vy; vcmZ += p.vz }
    vcmX /= totalParticles; vcmY /= totalParticles; vcmZ /= totalParticles
    for (let p of newParticles) { p.vx -= vcmX; p.vy -= vcmY; p.vz -= vcmZ }

    // Seed energy history at t=0 so the chart doesn't start with an empty left section.
    const { pe: pe0 } = calculateForces(newParticles)
    let sumV2 = 0
    for (const p of newParticles) sumV2 += (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz)
    // KE in Kelvin = 0.5 * m * sum(v^2) / KB_OVER_MU
    const ke0 = 0.5 * MASS_ARGON * sumV2 / KB_OVER_MU
    
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
    }, [numParticles, numParticlesA, numParticlesB, isBinary, temperature, isThreeD, randomNormal])

  useEffect(() => {
    initializeParticles()
  }, [])

  const convertTemperature = (tempK: number, unit: 'K' | 'C' | 'F'): number => {
    if (unit === 'K') return tempK
    if (unit === 'C') return tempK - 273.15
    if (unit === 'F') return (tempK - 273.15) * 9 / 5 + 32
    return tempK
  }

  const handleTemperatureUnitClick = () => {
    const nextUnit: 'K' | 'C' | 'F' = temperatureUnit === 'K' ? 'C' : temperatureUnit === 'C' ? 'F' : 'K'
    setTemperatureUnit(nextUnit)
  }

  // --- ALGORITHM 2: Velocity Verlet with Nosé-Hoover Thermostat (Real Units) ---
  const updatePhysics = (): boolean => {
    const currentParticles = [...particlesRef.current]
    const { dt, box, targetTemp } = paramsRef.current
    let zeta = nhStateRef.current.zeta
    
    timeRef.current += dt

    const N = currentParticles.length
    // Degrees of Freedom: 3N - 3 for 3D, 2N - 2 for 2D
    const dof = isThreeDRef.current ? (3 * N - 3) : (2 * N - 2)
    
    // Pre-calculate acceleration factor: (0.5 * dt * KB_OVER_MU) / MASS_ARGON
    const accFactor = (0.5 * dt * KB_OVER_MU) / MASS_ARGON
    
    // 1. Calculate Forces at t
    calculateForces(currentParticles)

    // 2. First Half-step Update
    for (let p of currentParticles) {
        const u = unwrappedRef.current[p.id]
        
        // Update Velocity (Half Step) with thermostat
        // F[K/Å] * accFactor -> a[Å/ps²]
        p.vx += (p.fx * accFactor) - (0.5 * dt * zeta * p.vx)
        p.vy += (p.fy * accFactor) - (0.5 * dt * zeta * p.vy)
        p.vz += (p.fz * accFactor) - (0.5 * dt * zeta * p.vz)

        // Update Position (Full Step)
        p.x += p.vx * dt
        p.y += p.vy * dt
        p.z += p.vz * dt

        // Update UNWRAPPED Position (Full Step) - FOR MSD
        u.x += p.vx * dt
        u.y += p.vy * dt
        u.z += p.vz * dt

        // --- BOUNDARY CONDITIONS ---
        const { useGravity, box } = paramsRef.current

        // X-Axis: Always Periodic (Infinite sides)
        if (p.x < 0) p.x += box
        if (p.x >= box) p.x -= box

        // Z-Axis: Always Periodic (Infinite depth)
        if (isThreeDRef.current) {
            if (p.z < 0) p.z += box
            if (p.z >= box) p.z -= box
        }

        // Y-Axis: Conditional (Periodic vs Reflective)
        if (useGravity) {
            // FLOOR (Reflective Hard Wall)
            if (p.y < 0) {
                p.y = -p.y // Reflect position
                p.vy *= -1 // Reverse velocity (Bounce)
            }
            // CEILING (Reflective Hard Wall)
            if (p.y >= box) {
                p.y = (2 * box) - p.y // Reflect position
                p.vy *= -1 // Reverse velocity (Bounce)
            }
        } else {
            // Standard Periodic Wrapping
            if (p.y < 0) p.y += box
            if (p.y >= box) p.y -= box
        }
    }

    // 3. Calculate Forces at t + dt
    const { pe, virial } = calculateForces(currentParticles)

    // 4. Update Thermostat Variable (Zeta)
    // T_inst (K) = (m * sum(v^2)) / (N_dof * KB_OVER_MU)
    let sumV2 = 0
    for (let p of currentParticles) {
        sumV2 += (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz)
    }
    const tempKelvin = (MASS_ARGON * sumV2) / (dof * KB_OVER_MU)

    // Update Thermostat (Zeta)
    // Q should have units of [Energy * Time^2]. In internal units: [K * ps^2].
    // Reference: Frenkel & Smit (common choice): Q = dof * T_target * (tau/2π)^2
    // Note: dof already scales with N, so we do NOT multiply by N again.
    const tau = 0.5
    const Q = (dof * targetTemp * tau * tau) / (4 * Math.PI * Math.PI)
    zeta += dt * (dof * (tempKelvin - targetTemp)) / Q

    // 5. Second Half-step Update
    let sumV2Full = 0
    for (let p of currentParticles) {
        const vxh = p.vx
        const vyh = p.vy
        const vzh = p.vz
        
        // Standard Nosé-Hoover update
        const divisor = 1.0 + 0.5 * dt * zeta
        p.vx = (vxh + p.fx * accFactor) / divisor
        p.vy = (vyh + p.fy * accFactor) / divisor
        p.vz = (vzh + p.fz * accFactor) / divisor
        
        sumV2Full += (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz)
    }

    nhStateRef.current.zeta = zeta

    // Final KE in Kelvin: KE = 0.5 * m * sum(v^2) / KB_OVER_MU
    const kineticEnergyFull = 0.5 * MASS_ARGON * sumV2Full / KB_OVER_MU

    // --- Calculate Pressure ---
    // 2D vs 3D Pressure Logic
    const volume = isThreeDRef.current ? Math.pow(box, 3) : Math.pow(box, 2)
    const virialDivisor = isThreeDRef.current ? 3.0 : 2.0

    // Temperature from full-step velocities (more consistent with pressure sampling)
    const tempKelvinFull = (MASS_ARGON * sumV2Full) / (dof * KB_OVER_MU)

    let pressureBar = NaN
    if (isThreeDRef.current) {
        // Standard 3D Pressure: (N k T + W/d) / V
        // internalPressure is in K/Å³
        const internalPressure = (N * tempKelvinFull + virial / virialDivisor) / volume
        pressureBar = internalPressure * PRESSURE_CONV_BAR
    } else {
        // 2D Pressure is ill-defined in "bar" without a thickness.
        pressureBar = NaN
    }
    
    // --- CRASH DETECTION ---
    // Check 1: Is the Temperature absurdly high? (e.g. > 10,000 K)
    // Check 2: Are the numbers broken? (NaN or Infinity)
    if (Number.isNaN(tempKelvinFull) || !Number.isFinite(tempKelvinFull) || tempKelvinFull > 10000) {
        return false // FAILURE: Simulation Exploded
    }
    
    setPressure(pressureBar)

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

        if (msdHistoryRef.current.length > MAX_HISTORY_POINTS) {
            msdHistoryRef.current.splice(0, msdHistoryRef.current.length - MAX_HISTORY_POINTS)
        }

    // Update State
    particlesRef.current = currentParticles
    setKineticEnergy(kineticEnergyFull)
    setThermostatZeta(zeta)

    // Track Energy History every physics step (real samples, no synthetic padding)
    // Both KE and PE are in Kelvin (kB-normalized)
    const tot = kineticEnergyFull + pe
    energyHistoryRef.current.push({ time: timeRef.current, ke: kineticEnergyFull, pe, tot })

        if (energyHistoryRef.current.length > MAX_HISTORY_POINTS) {
            energyHistoryRef.current.splice(0, energyHistoryRef.current.length - MAX_HISTORY_POINTS)
        }
    
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
    
    return true // SUCCESS: Physics step valid
  }

  // Helper for re-render
  const [_, setFrameCount] = useState(0)

  // --- Animation Loop ---
  const animate = () => {
    const isValid = updatePhysics()
    
    if (!isValid) {
        setRunning(false)
        setIsCrashed(true)
        return // Stop the loop immediately
    }
    
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
      setIsCrashed(false)
      initializeParticles()
  }

  const handleExportData = () => {
      // 1. Create CSV Header
      let csvContent = "data:text/csv;charset=utf-8,"
      csvContent += "Time(ps),TotalEnergy(K),Potential(K),Kinetic(K),MSD(A^2)\n"

      // 2. Merge Energy and MSD history
      // Assuming arrays are roughly synchronized by frame index
      energyHistoryRef.current.forEach((row, index) => {
          const msdVal = msdHistoryRef.current[index] ? msdHistoryRef.current[index].msd : 0
          const rowString = `${row.time.toFixed(4)},${row.tot.toFixed(2)},${row.pe.toFixed(2)},${row.ke.toFixed(2)},${msdVal.toFixed(4)}`
          csvContent += rowString + "\n"
      })

      // 3. Trigger Download
      const encodedUri = encodeURI(csvContent)
      const link = document.createElement("a")
      link.setAttribute("href", encodedUri)
      link.setAttribute("download", "simulation_data.csv")
      document.body.appendChild(link)
      link.click()
      document.body.removeChild(link)
  }
  
  // Reset effect when numParticles changes
  useEffect(() => {
      initializeParticles()
  }, [initializeParticles])

  // Keep isThreeDRef in sync with isThreeD state
  useEffect(() => {
      isThreeDRef.current = isThreeD
  }, [isThreeD])

  // Re-initialize particles with Z coordinates when switching to 3D
  useEffect(() => {
      if (isThreeD) {
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

        const getPercentile = (values: number[], p: number) => {
                if (values.length === 0) return 0
                const sorted = [...values].sort((a, b) => a - b)
                const idx = Math.min(sorted.length - 1, Math.max(0, Math.floor(p * (sorted.length - 1))))
                return sorted[idx]
        }

        // Keep arrow lengths visually stable across parameter/temperature changes.
        // We normalize by the 95th percentile magnitude (avoids single outliers dominating).
        const maxArrowLenData = BOX_SIZE * 0.12
        const velRef2D = getPercentile(particles.map(p => Math.hypot(p.vx, p.vy)), 0.95)
        const forceRef2D = getPercentile(particles.map(p => Math.hypot(p.fx, p.fy)), 0.95)
        const velScale2D = velRef2D > 1e-12 ? (maxArrowLenData / velRef2D) : 0
        const forceScale2D = forceRef2D > 1e-12 ? (maxArrowLenData / forceRef2D) : 0

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

        const scale = type === 0 ? velScale2D : forceScale2D
    // Hard clamp per-arrow length so a single extreme outlier cannot draw across the box.
    const dxScaled = dxVal * scale
    const dyScaled = dyVal * scale
    const lenScaled = Math.hypot(dxScaled, dyScaled)
    if (!Number.isFinite(lenScaled) || lenScaled < 1e-9) return
    const clampFactor = lenScaled > maxArrowLenData ? (maxArrowLenData / lenScaled) : 1

    const start = api.coord([x, y])
    const end = api.coord([x + dxScaled * clampFactor, y + dyScaled * clampFactor])
        
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
        // Color Logic: Type 0 = Blue, Type 1 = Red
        const baseColor = p.type === 0 ? '#3b82f6' : '#ef4444'
        
        return {
            // Keep id at index 5 for click handler; store type at index 6 for symbol sizing.
            value: [p.x, p.y, Math.sqrt(p.vx**2 + p.vy**2), p.id, 0, p.id, p.type],
            itemStyle: {
                color: baseColor,
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
                         // data is the raw `value` array; we stash particle type at index 6.
                         const pType = (Array.isArray(data) && (data[6] === 1 || data[6] === 0)) ? (data[6] as 0 | 1) : 0
                         const sigLocal = pType === 0 ? sigmaA : sigmaB
                         return Math.max(5, Math.min(30, sigLocal * 2.5))
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
                  return `r: ${r} Å<br/>g(r): ${g}`
              }
          },
          grid: { left: 40, right: 10, top: 30, bottom: 40 },
          xAxis: {
              type: 'value',
              name: 'r (Distance) [Å]',
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
                  smooth: true,
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
                      lines.push(`${name}: ${Number(val).toPrecision(4)} K`)
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
          yAxis: { type: 'value', name: 'Energy (K)', nameLocation: 'middle', nameGap: 55, nameTextStyle: { color: textColor, fontSize: 10, fontFamily: 'Merriweather Sans' }, splitLine: { show: false }, axisLine: { show: true, lineStyle: { color: textColor }, onZero: false }, axisTick: { show: true }, axisLabel: { color: textColor, fontFamily: 'Merriweather Sans', margin: 5 } },
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
                  return `Time: ${t} ps<br/>MSD: ${msd} Å²`
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
              name: 'MSD (Å²)',
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
              eps: <span>L-J Epsilon (ε/k<sub>B</sub>)</span>,
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
                <div className="relative group">
                    <Button variant="outline" size="icon" onClick={handleReset} title="Reset">
                        <RotateCcw className="h-4 w-4" />
                    </Button>
                    <div className="absolute left-1/2 -translate-x-1/2 bottom-full mb-2 z-50 hidden w-32 p-2 text-xs bg-popover text-popover-foreground border rounded-md shadow-md group-hover:block animate-in fade-in zoom-in-95 duration-200 text-center whitespace-normal">
                        Reset Simulation
                    </div>
                </div>
                <div className="relative group">
                    <Button variant="outline" size="icon" onClick={() => initializeParticles()} title="Zap (Re-Minimize)">
                        <Zap className="h-4 w-4 text-yellow-500" />
                    </Button>
                    <div className="absolute left-1/2 -translate-x-1/2 bottom-full mb-2 z-50 hidden w-32 p-2 text-xs bg-popover text-popover-foreground border rounded-md shadow-md group-hover:block animate-in fade-in zoom-in-95 duration-200 text-center whitespace-normal">
                        Minimize Energy
                    </div>
                </div>
              </div>

              <div className="space-y-4 border-b pb-4">
                {/* PARTICLE COUNTS */}
                <div className={`space-y-2 ${running ? 'opacity-50 pointer-events-none' : ''}`}>
                    {!isBinary ? (
                        <>
                            <Label>Count: {numParticles}</Label>
                            <Slider 
                                disabled={running}
                                value={[numParticles]} min={10} max={200} step={10}
                                onValueChange={(v) => handleParticleCountChange(v[0])}
                            />
                        </>
                    ) : (
                        <div className="grid grid-cols-2 gap-4">
                            <div className="space-y-2">
                                <Label className="text-blue-500">Count A: {numParticlesA}</Label>
                                <Slider 
                                    disabled={running}
                                    value={[numParticlesA]} min={5} max={100} step={5}
                                    onValueChange={(v) => setNumParticlesA(v[0])}
                                />
                            </div>
                            <div className="space-y-2">
                                <Label className="text-red-500">Count B: {numParticlesB}</Label>
                                <Slider 
                                    disabled={running}
                                    value={[numParticlesB]} min={5} max={100} step={5}
                                    onValueChange={(v) => setNumParticlesB(v[0])}
                                />
                            </div>
                        </div>
                    )}
                </div>

                {/* TEMPERATURE: SAFE to change during run */}
                <div className="space-y-2">
                    <div className="flex justify-between items-center">
                        <Label><span>Reduced Temp (T/ε): {reducedTemp.toFixed(2)}</span></Label>
                        <Label 
                            className="text-muted-foreground text-xs cursor-pointer hover:text-foreground transition-colors"
                            onClick={handleTemperatureUnitClick}
                        >
                            <span>(T = {convertTemperature(temperature, temperatureUnit).toFixed(1)} {temperatureUnit})</span>
                        </Label>
                    </div>
                    <Slider 
                        value={[reducedTemp]} 
                        min={reducedMin} 
                        max={reducedMax} 
                        step={REDUCED_STEP}
                        onValueChange={(v) => setReducedTemp(snapToStep(v[0], REDUCED_STEP))}
                    />
                </div>

                {/* TIME STEP: Safe to change during run */}
                <div className="space-y-2">
                    <Label>Time Step (dt): {timeStep.toFixed(4)} ps</Label>
                    <Slider
                        value={[timeStep]} min={0.0005} max={0.03} step={0.0005}
                        onValueChange={(v) => setTimeStep(v[0])}
                    />
                </div>

                {useGravity && (
                    <div className={`space-y-2 ${running ? 'opacity-50 pointer-events-none' : ''} animate-in fade-in`}>
                        <Label>
                            Gravity Strength (G-Force): {gravityStrength.toFixed(1)} K/Å
                        </Label>
                        <Slider
                            disabled={running}
                            value={[gravityStrength]} min={0.1} max={20.0} step={0.5}
                            onValueChange={(v) => setGravityStrength(v[0])}
                        />
                    </div>
                )}
              </div>

              <div className="space-y-4">
                {/* 1. POTENTIAL SELECTOR WITH TOOLTIP */}
                <div className={`space-y-2 ${running ? 'opacity-50 pointer-events-none' : ''}`}>
                    <div className="flex items-center justify-between">
                        <div className="flex items-center gap-2">
                            <Label>
                                {!isBinary && showLJReferenceHints && ljPairMatchesA.length === 1
                                    ? `${ljPairMatchesA[0].substance}'s Interatomic Potential`
                                    : 'Interatomic Potential'}
                            </Label>
                            <div className={`relative group flex items-center justify-center ${showPolarLJWarning ? 'visible' : 'invisible'}`}>
                                <AlertCircle className="h-4 w-4 text-amber-500 cursor-help hover:text-amber-400 transition-colors" />
                                <div className="absolute left-1/2 -translate-x-1/2 top-full mt-2 z-50 hidden w-80 p-3 text-xs bg-popover text-popover-foreground border rounded-md shadow-md group-hover:block animate-in fade-in zoom-in-95 duration-200 space-y-1">
                                    <div className="font-semibold">Polar Molecules</div>
                                    <div>Polar / hydrogen-bonding molecules shown use <span className="font-semibold">effective LJ parameters only</span> (no electrostatics / hydrogen bonding).</div>
                                    <div>Results are <span className="font-semibold">qualitative</span>, not predictive.{polarLJNames ? ` Detected: ${polarLJNames}.` : ''}</div>
                                </div>
                            </div>
                        </div>
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

                {/* PARAMETERS (Epsilon & Sigma) */}
                <div className={`space-y-2 ${running ? 'opacity-50 pointer-events-none' : ''}`}>
                     {!isBinary ? (
                         /* OLD SINGLE SLIDERS */
                         <>
                             <Label className="text-white"><span>{labels.eps}: {epsilonA.toFixed(0)} K</span></Label>
                             <Slider
                                 disabled={running}
                                 value={[epsilonA]} min={5} max={1000} step={5}
                                 onValueChange={(v) => setEpsilonA(v[0])}
                             />
                             {showLJReferenceHints && ljPairMatchesA.length !== 1 && (
                                 <div className="text-[10px] leading-tight text-muted-foreground">
                                     {ljEpsilonMatchesA.map((m) => (
                                         <div key={`epsA-${m.formula}`}>
                                             {m.substance} ({m.formula}) — {m.epsilonOverKbK.toFixed(1)} K
                                         </div>
                                     ))}
                                 </div>
                             )}
                             
                             <Label className="text-white"><span>{labels.sig}: {sigmaA.toFixed(2)} Å</span></Label>
                             <Slider
                                 disabled={running}
                                 value={[sigmaA]} min={2.0} max={6.5} step={0.01}
                                 onValueChange={(v) => setSigmaA(v[0])}
                             />
                             {showLJReferenceHints && ljPairMatchesA.length !== 1 && (
                                 <div className="text-[10px] leading-tight text-muted-foreground">
                                     {ljSigmaMatchesA.map((m) => (
                                         <div key={`sigA-${m.formula}`}>
                                             {m.substance} ({m.formula}) — {m.sigmaAngstrom.toFixed(3)} Å
                                         </div>
                                     ))}
                                 </div>
                             )}
                         </>
                     ) : (
                         /* NEW SPLIT SLIDERS */
                         <div className="grid grid-cols-2 gap-4">
                             {/* Left Column: Species A (Blue) */}
                             <div className="space-y-2 border-r pr-2">
                                 <Label className="text-blue-500 text-xs font-bold">
                                     {showLJReferenceHints && ljPairMatchesA.length === 1
                                         ? ljPairMatchesA[0].substance
                                         : 'Species A (Blue)'}
                                 </Label>
                                 
                                 <Label className="text-[10px]"><span>{labels.eps}: {epsilonA.toFixed(0)} K</span></Label>
                                 <Slider 
                                     disabled={running}
                                     value={[epsilonA]} min={5} max={1000} step={5}
                                     onValueChange={(v) => setEpsilonA(v[0])} 
                                 />
                                 {showLJReferenceHints && ljPairMatchesA.length !== 1 && (
                                     <div className="text-[10px] leading-tight text-muted-foreground">
                                         {ljEpsilonMatchesA.map((m) => (
                                             <div key={`epsA2-${m.formula}`}>
                                                 {m.substance} ({m.formula})
                                             </div>
                                         ))}
                                     </div>
                                 )}
                                 
                                 <Label className="text-[10px]"><span>{labels.sig}: {sigmaA.toFixed(2)} Å</span></Label>
                                 <Slider 
                                     disabled={running}
                                     value={[sigmaA]} min={2} max={6} step={0.01}
                                     onValueChange={(v) => setSigmaA(v[0])} 
                                 />
                                 {showLJReferenceHints && ljPairMatchesA.length !== 1 && (
                                     <div className="text-[10px] leading-tight text-muted-foreground">
                                         {ljSigmaMatchesA.map((m) => (
                                             <div key={`sigA2-${m.formula}`}>
                                                 {m.substance} ({m.formula})
                                             </div>
                                         ))}
                                     </div>
                                 )}
                             </div>

                             {/* Right Column: Species B (Red) */}
                             <div className="space-y-2">
                                 <Label className="text-red-500 text-xs font-bold">
                                     {showLJReferenceHints && ljPairMatchesB.length === 1
                                         ? ljPairMatchesB[0].substance
                                         : 'Species B (Red)'}
                                 </Label>
                                 
                                 <Label className="text-[10px]"><span>{labels.eps}: {epsilonB.toFixed(0)} K</span></Label>
                                 <Slider 
                                     disabled={running}
                                     value={[epsilonB]} min={5} max={1000} step={5}
                                     onValueChange={(v) => setEpsilonB(v[0])} 
                                 />
                                 {showLJReferenceHints && ljPairMatchesB.length !== 1 && (
                                     <div className="text-[10px] leading-tight text-muted-foreground">
                                         {ljEpsilonMatchesB.map((m) => (
                                             <div key={`epsB-${m.formula}`}>
                                                 {m.substance} ({m.formula})
                                             </div>
                                         ))}
                                     </div>
                                 )}
                                 
                                 <Label className="text-[10px]"><span>{labels.sig}: {sigmaB.toFixed(2)} Å</span></Label>
                                 <Slider 
                                     disabled={running}
                                     value={[sigmaB]} min={2} max={6} step={0.01}
                                     onValueChange={(v) => setSigmaB(v[0])} 
                                 />
                                 {showLJReferenceHints && ljPairMatchesB.length !== 1 && (
                                     <div className="text-[10px] leading-tight text-muted-foreground">
                                         {ljSigmaMatchesB.map((m) => (
                                             <div key={`sigB-${m.formula}`}>
                                                 {m.substance} ({m.formula})
                                             </div>
                                         ))}
                                     </div>
                                 )}
                             </div>
                         </div>
                     )}
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
                    <div className={`p-2 bg-muted rounded ${!isThreeD ? 'col-span-2' : ''}`}>
                        <p className="text-[10px] text-muted-foreground">Thermostat (ζ)</p>
                        <p className={`font-mono font-bold text-xs ${thermostatZeta > 0 ? 'text-red-500' : 'text-blue-500'}`}>
                            {thermostatZeta.toFixed(3)}
                        </p>
                    </div>
                    {isThreeD && (
                        <div className="p-2 bg-muted rounded">
                            <p className="text-[10px] text-muted-foreground">Pressure (P)</p>
                            <p className="font-mono font-bold text-xs text-purple-600">
                                {pressure.toFixed(1)} bar
                            </p>
                        </div>
                    )}
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
              
              <CardContent className="pt-0 -mt-5">
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
              <Card className={`w-full p-0 overflow-hidden bg-muted/20 border border-border relative transition-colors duration-700 ${!isThreeD && showGravityNotice && !gravityNoticeFading ? 'border-b-amber-500' : 'border-b-border'}`} style={{ height: '500px' }}>
              {/* Gravity Toggle - Top Left */}
              <div className="absolute top-2 left-2 z-10">
                  <Button
                      size="sm"
                      variant={useGravity ? 'default' : 'outline'}
                      onClick={() => setUseGravity(!useGravity)}
                      className={`h-8 text-xs shadow-none hover:bg-background/80 active:bg-background/80 focus-visible:ring-0 focus-visible:outline-none font-[Merriweather_Sans] ${
                          mounted ? (isDarkTheme
                              ? 'bg-background/80 text-white border border-white/20'
                              : 'bg-background/80 text-black border border-black/20')
                          : ''
                      }`}
                  >
                      {useGravity ? '⬇ Gravity ON' : '⬇ Gravity OFF'}
                  </Button>
              </div>

              {/* Mode Switch - Bottom Left */}
              <div className="absolute bottom-2 left-2 z-10">
                  <Button
                      size="sm"
                      variant="outline"
                      disabled={running}
                      onClick={() => {
                          setRunning(false)
                          setIsCrashed(false)
                          setSelectedId(null)
                          setIsBinary(prev => !prev)
                      }}
                      className={`h-8 text-xs shadow-none hover:bg-background/80 active:bg-background/80 focus-visible:ring-0 focus-visible:outline-none font-[Merriweather_Sans] ${
                          mounted ? (isDarkTheme
                              ? 'bg-background/80 text-white border border-white/20'
                              : 'bg-background/80 text-black border border-black/20')
                          : ''
                      }`}
                  >
                      {isBinary ? 'Switch to Pure' : 'Switch to Binary'}
                  </Button>
              </div>

              {showGravityNotice && (
                  <div className="absolute top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2 z-10 pointer-events-none">
                      <div
                          className={`rounded-md border border-amber-500/40 bg-amber-500/10 px-3 py-2 text-[11px] leading-snug text-amber-900 dark:text-amber-100 backdrop-blur-sm shadow-sm transition-opacity duration-700 whitespace-nowrap ${
                              gravityNoticeFading ? 'opacity-0' : 'opacity-100'
                          }`}
                      >
                          <span className="font-semibold">Gravity ON:</span> Y becomes a floor/ceiling (reflective walls; no periodic wrap).
                      </div>
                  </div>
              )}

              {/* CRASH WARNING ALERT - Centered Overlay */}
              {isCrashed && (
                  <div className="absolute inset-0 z-20 flex items-center justify-center pointer-events-none">
                      <div className="bg-red-500/15 border border-red-500/50 rounded-lg p-6 flex items-start gap-4 animate-in fade-in zoom-in-95 w-96 shadow-2xl pointer-events-auto">
                          <div className="p-1 bg-red-500 rounded-full mt-1 flex-shrink-0">
                              <Activity className="h-5 w-5 text-white" />
                          </div>
                          <div className="space-y-2">
                              <h4 className="font-bold text-base text-red-500">Simulation Unstable!</h4>
                              <div className="text-sm text-muted-foreground space-y-1">
                                  <div>The system exploded.</div>
                                  <div>This could be due to <strong>Time Step (dt)</strong> being too high, model parameters (ε, σ) being too extreme, or too many particles.</div>
                                  <div>Try reducing dt, lowering particle count, or adjusting potential parameters.</div>
                              </div>
                              <div className="pt-3">
                                  <Button 
                                      size="sm" 
                                      variant="destructive" 
                                      onClick={() => {
                                          handleReset()
                                      }}
                                      className="h-8 text-sm"
                                  >
                                      Reset
                                  </Button>
                              </div>
                          </div>
                      </div>
                  </div>
              )}
              
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
                                sigmaA={sigmaA}
                                sigmaB={sigmaB}
                                selectedId={selectedId}
                                isDark={isDarkTheme}
                                showVelocity={showVelocity}
                                showForce={showForce}
                                showGravityFloorHighlight={showGravityNotice}
                                gravityFloorFading={gravityNoticeFading}
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

                  {/* Structure Metrics Info Card (Center Right) - Only show for RDF when simulation has data */}
                  {mounted && analysisView === 'rdf' && rdfData.length > 0 && (() => {
                      const { state, color } = detectSystemState()
                      return (
                          <div className={`absolute top-1/4 -translate-y-1/2 right-2 z-10 backdrop-blur-sm p-2 rounded border text-[10px] space-y-1 shadow-sm w-40 ${
                              isDarkTheme
                                  ? 'bg-background/80 border-white/20'
                                  : 'bg-background/80 border-black/20'
                          }`}>
                              <div className="flex justify-between">
                                  <span className="text-muted-foreground">1st Peak (r):</span>
                                  <span className="font-mono font-bold">{structureMetrics.peakDist.toFixed(2)} Å</span>
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

                  {/* MSD Diffusion Info Box (Placed at top-3/4 as requested) */}
                  {mounted && analysisView === 'msd' && msdData.length > 20 && (() => {
                      // 1. Get raw slope from your existing detector
                      const { state, color, diffusion } = detectSystemState()
                      
                      // 2. Calculate D based on dimensionality
                      // D = Slope / (2 * dimensions)
                      const dimDivisor = isThreeD ? 6.0 : 4.0
                      const diffCoeff = diffusion / dimDivisor

                      return (
                          <div className={`absolute top-3/5 -translate-y-1/2 right-2 z-10 backdrop-blur-sm p-2 rounded border text-[10px] space-y-1 shadow-sm w-44 ${
                              isDarkTheme
                                  ? 'bg-background/80 border-white/20'
                                  : 'bg-background/80 border-black/20'
                          }`}>
                              <div className="flex justify-between">
                                  <span className="text-muted-foreground">Slope (MSD/t):</span>
                                  <span className="font-mono font-bold">{diffusion.toFixed(3)}</span>
                              </div>
                              <div className="flex justify-between items-center">
                                  <span className="text-muted-foreground">Diff. Coeff (D):</span>
                                  <span className="font-mono font-bold text-amber-500 text-xs">
                                    {diffCoeff.toFixed(4)} Å²/ps
                                  </span>
                              </div>
                              <div className="flex justify-between border-t pt-1 mt-1">
                                  <span className="text-muted-foreground">Est. Phase:</span>
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

                  {/* Export Data Button - Bottom Right */}
                  <div className="absolute bottom-2 right-2 z-10">
                      <Button
                          size="sm"
                          variant="outline"
                          onClick={handleExportData}
                          className={`h-8 text-xs shadow-none hover:bg-background/80 active:bg-background/80 focus-visible:ring-0 focus-visible:outline-none font-[Merriweather_Sans] ${
                              mounted ? (isDarkTheme
                                  ? 'bg-background/80 text-white border border-white/20'
                                  : 'bg-background/80 text-black border border-black/20')
                              : ''
                          }`}
                      >
                          ⬇ Export
                      </Button>
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
    scale,
    maxLength
}: { 
  start: [number, number, number], 
  vector: [number, number, number], 
  color: string, 
    scale: number,
    maxLength: number
}) => {
  const dir = new THREE.Vector3(...vector)
  const length = dir.length()
  
  // Don't render tiny vectors
  if (length < 1e-6) return null

  dir.normalize()
    // Scale by percentile reference, then hard clamp so outliers never draw absurdly long arrows.
    const arrowLengthRaw = length * scale
    const arrowLength = Math.min(arrowLengthRaw, Math.max(0, maxLength))
    if (!Number.isFinite(arrowLength) || arrowLength < 1e-6) return null

    // Ensure geometry stays valid even when clamped very small.
    const headLength = Math.min(arrowLength * 0.2, 2, arrowLength)
    const headWidth = Math.min(arrowLength * 0.08, 1)
    const shaftLength = Math.max(0, arrowLength - headLength)
  const shaftRadius = headWidth * 0.25

  // Calculate rotation to align with direction vector
  const quaternion = new THREE.Quaternion()
  quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir)

  // Memoize material to avoid recreating on every render
  const material = useMemo(() => new THREE.MeshBasicMaterial({ color }), [color])

  return (
    <group position={start} quaternion={[quaternion.x, quaternion.y, quaternion.z, quaternion.w]}>
      {/* Shaft (cylinder) */}
      <mesh position={[0, shaftLength / 2, 0]} material={material}>
        <cylinderGeometry args={[shaftRadius, shaftRadius, shaftLength, 8]} />
      </mesh>
      
      {/* Head (cone) */}
      <mesh position={[0, shaftLength + headLength / 2, 0]} material={material}>
        <coneGeometry args={[headWidth, headLength, 8]} />
      </mesh>
    </group>
  )
}

const Molecule3D = ({ 
  particles, 
  boxSize, 
    sigmaA,
    sigmaB,
  selectedId, 
  isDark,
  showVelocity,
  showForce,
    showGravityFloorHighlight,
    gravityFloorFading,
  onParticleClick
}: { 
  particles: Particle[], 
  boxSize: number, 
    sigmaA: number,
    sigmaB: number,
  selectedId: number | null, 
  isDark: boolean,
  showVelocity: boolean,
  showForce: boolean,
    showGravityFloorHighlight: boolean,
    gravityFloorFading: boolean,
  onParticleClick: (id: number) => void
}) => {
  const particleColor = isDark ? "#ffffff" : "#3b82f6"
  const wireframeColor = isDark ? "#ffffff" : "#000000"
  
  // Colors for arrows
  const velColor = isDark ? '#4ade80' : '#16a34a'
  const forceColor = isDark ? '#facc15' : '#d97706'

    const getPercentile = (values: number[], p: number) => {
        if (values.length === 0) return 0
        const sorted = [...values].sort((a, b) => a - b)
        const idx = Math.min(sorted.length - 1, Math.max(0, Math.floor(p * (sorted.length - 1))))
        return sorted[idx]
    }

    // Normalize arrow sizes by current magnitudes so they don't explode with temperature/params.
    const maxArrowLen3D = boxSize * 0.12
    const velRef3D = useMemo(
        () => getPercentile(particles.map(p => Math.hypot(p.vx, p.vy, p.vz)), 0.95),
        [particles]
    )
    const forceRef3D = useMemo(
        () => getPercentile(particles.map(p => Math.hypot(p.fx, p.fy, p.fz)), 0.95),
        [particles]
    )
    const velScale3D = velRef3D > 1e-12 ? (maxArrowLen3D / velRef3D) : 0
    const forceScale3D = forceRef3D > 1e-12 ? (maxArrowLen3D / forceRef3D) : 0

    // Floor highlight (used when gravity notice is shown in 3D)
    const FloorHighlight = ({ boxSize, isDark, show, fading }: { boxSize: number; isDark: boolean; show: boolean; fading: boolean }) => {
        const materialRef = React.useRef<THREE.MeshBasicMaterial | null>(null)
        const targetOpacity = show ? (fading ? 0 : 0.22) : 0

        useFrame((_, delta) => {
            if (!materialRef.current) return
            const current = materialRef.current.opacity ?? 0
            // Exponential smoothing; tuned to roughly match the 700ms UI fade.
            const t = 1 - Math.exp(-delta * 10)
            materialRef.current.opacity = THREE.MathUtils.lerp(current, targetOpacity, t)
        })

        if (!show) return null

        return (
            <mesh position={[boxSize / 2, 0.01, boxSize / 2]} rotation={[-Math.PI / 2, 0, 0]}>
                <planeGeometry args={[boxSize, boxSize]} />
                <meshBasicMaterial
                    ref={materialRef}
                    color={isDark ? '#facc15' : '#f59e0b'}
                    transparent
                    opacity={0}
                    depthWrite={false}
                    side={THREE.DoubleSide}
                />
            </mesh>
        )
    }

  return (
    <group>
            <FloorHighlight boxSize={boxSize} isDark={isDark} show={showGravityFloorHighlight} fading={gravityFloorFading} />
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
        
        // Color by type: Type 0 = Blue, Type 1 = Red
        const baseColor = p.type === 0 
            ? (isDark ? "#60a5fa" : "#3b82f6") // Blue for A
            : (isDark ? "#f87171" : "#ef4444") // Red for B

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
                            <sphereGeometry args={[(p.type === 0 ? sigmaA : sigmaB) * 0.5, 16, 16]} />
              <meshStandardMaterial 
                color={isSelected ? "#facc15" : baseColor}
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
                                                                scale={velScale3D}
                                                                maxLength={maxArrowLen3D}
              />
            )}

            {/* Force Arrow */}
            {showF && (
              <Arrow3D 
                start={[p.x, p.y, p.z]} 
                vector={[p.fx, p.fy, p.fz]} 
                color={forceColor} 
                                                                scale={forceScale3D}
                                                                maxLength={maxArrowLen3D}
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