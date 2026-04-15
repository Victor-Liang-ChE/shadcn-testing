// ─── Canopy Zustand Store ────────────────────────────────────────────────────
// Manages all simulation state: compounds, streams, units, and solver results.
// ─────────────────────────────────────────────────────────────────────────────

import { create } from 'zustand';
import type {
  CanopyCompound,
  FluidPackage,
  MaterialStream,
  UnitOperation,
  UnitOpType,
  UnitPort,
} from './types';
import { createDefaultStream } from './types';
import { flashPT_Ideal } from './thermo';
import { solveMixer, solveHeater, solveFlashDrum, solveValve } from './solver';

// ────────────────────────────────────────────────────────────────
// Store Interface
// ────────────────────────────────────────────────────────────────

export type SimScreen = 'setup' | 'flowsheet';

/** Unit system for display (all internal calcs use SI) */
export interface UnitSystem {
  temperature: 'C' | 'K' | 'F';
  pressure: 'bar' | 'atm' | 'kPa' | 'Pa' | 'psi';
  flow: 'mol/s' | 'kmol/h' | 'kg/s' | 'kg/h';
  energy: 'kW' | 'W' | 'BTU/h' | 'kcal/h';
}

export const DEFAULT_UNIT_SYSTEM: UnitSystem = {
  temperature: 'C',
  pressure: 'bar',
  flow: 'mol/s',
  energy: 'kW',
};

/** Tab can be a unit or a stream */
export type TabItem = { type: 'unit'; id: string } | { type: 'stream'; id: string };

interface CanopyState {
  // Navigation
  currentScreen: SimScreen;
  setScreen: (screen: SimScreen) => void;

  // Setup
  compounds: CanopyCompound[];
  fluidPackage: FluidPackage;
  unitSystem: UnitSystem;
  setCompounds: (compounds: CanopyCompound[]) => void;
  addCompound: (compound: CanopyCompound) => void;
  removeCompound: (name: string) => void;
  setFluidPackage: (pkg: FluidPackage) => void;
  setUnitSystem: (patch: Partial<UnitSystem>) => void;

  // Flowsheet
  streams: Record<string, MaterialStream>;
  units: Record<string, UnitOperation>;
  addStream: (stream: MaterialStream) => void;
  updateStream: (id: string, patch: Partial<MaterialStream>) => void;
  removeStream: (id: string) => void;
  addUnit: (unit: UnitOperation) => void;
  updateUnit: (id: string, patch: Partial<UnitOperation>) => void;
  removeUnit: (id: string) => void;
  /** Add a new unit operation to the flowsheet from the object palette */
  addUnitFromPalette: (type: UnitOpType, position: { x: number; y: number }) => void;
  /** Add a new feed stream to the flowsheet */
  addFeedStream: () => void;
  /** Get the next available ID for a stream or unit */
  nextStreamId: () => string;
  nextUnitId: () => string;

  // Tabs (open unit/stream details)
  openTabs: TabItem[];
  activeTab: TabItem | null;
  openTab: (item: TabItem) => void;
  closeTab: (item: TabItem) => void;
  setActiveTab: (item: TabItem | null) => void;

  // Solver
  solved: boolean;
  solverErrors: string[];
  solveAll: () => void;

  // Preload
  loadPreset: (preset: 'benzene-toluene') => void;

  // Node positions (for React Flow drag persistence)
  nodePositions: Record<string, { x: number; y: number }>;
  setNodePosition: (nodeId: string, pos: { x: number; y: number }) => void;
}

// ────────────────────────────────────────────────────────────────
// Preloaded Data — Benzene + Toluene Flash Separation
// ────────────────────────────────────────────────────────────────

// DIPPR coefficients from the Aspen/NIST compound_properties table.
// These are hardcoded for the preloaded example so no DB fetch is needed on first load.
const BENZENE: CanopyCompound = {
  name: 'BENZENE',
  displayName: 'Benzene',
  molecularWeight: 78.112,
  Tc_K: 562.1,
  Pc_Pa: 4893997.5,
  omega: 0.212,
  dipprCoeffs: {
    // eqno 16: A + exp(B/T + C + D·T + E·T²)  [J/kmol/K]
    IdealGasHeatCapacityCp: { A: 34010.24, B: -588.0978, C: 12.81777, D: -0.000197306, E: 5.142899e-8 },
    // eqno 106: A·(1-Tr)^(B + C·Tr + D·Tr² + E·Tr³)  [J/kmol]
    HeatOfVaporization: { A: 4.881e7, B: 0.61066, C: -0.25882, D: 0.032238, E: 0.022475 },
    // eqno 16: A + exp(B/T + C + D·T + E·T²)  [J/kmol/K]
    LiquidHeatCapacityCp: { A: 111460, B: -1854.3, C: 22.399, D: -0.028936, E: 0.000028991 },
  },
  antoine: {
    C1: 73.862398, C2: -5970.4385, C3: 0, C4: 0.0055376032,
    C5: -8.0797644, C6: 6.6129752e-18, C7: 6,
    Tmin_K: 293.36, Tmax_K: 562.1,
  },
};

const TOLUENE: CanopyCompound = {
  name: 'TOLUENE',
  displayName: 'Toluene',
  molecularWeight: 92.138,
  Tc_K: 591.7,
  Pc_Pa: 4113795,
  omega: 0.257,
  dipprCoeffs: {
    IdealGasHeatCapacityCp: { A: 47225, B: -565.85, C: 12.856, D: 0.000005535, E: -1.998e-8 },
    HeatOfVaporization: { A: 5.3752e7, B: 0.50341, C: 0.24755, D: -0.72898, E: 0.37794 },
    LiquidHeatCapacityCp: { A: 28291, B: 48.171, C: 10.912, D: 0.0020542, E: 8.7875e-7 },
  },
  antoine: {
    C1: 71.277464, C2: -6413.286, C3: 0, C4: 0.0041663022,
    C5: -7.5053519, C6: 5.419977e-18, C7: 6,
    Tmin_K: 318.72, Tmax_K: 591.7,
  },
};

function createPresetFlowsheet(): {
  streams: Record<string, MaterialStream>;
  units: Record<string, UnitOperation>;
} {
  const n = 2; // benzene + toluene

  // ── Feed Stream ──
  const feed: MaterialStream = {
    ...createDefaultStream('S1', 'Feed', n),
    T_K: 303.15,          // 30°C
    P_Pa: 101325,         // 1 atm
    totalFlow_molps: 100, // 100 mol/s
    moleFractions: [0.5, 0.5], // 50/50 benzene/toluene
    targetUnitId: 'U1',
  };

  // ── Heater Outlet ──
  const heaterOut: MaterialStream = {
    ...createDefaultStream('S2', 'Heater Out', n),
    T_K: 373.15,
    P_Pa: 101325,
    totalFlow_molps: 100,
    moleFractions: [0.5, 0.5],
    sourceUnitId: 'U1',
    targetUnitId: 'U2',
  };

  // ── Flash Vapor ──
  const flashVapor: MaterialStream = {
    ...createDefaultStream('S3', 'Flash Vapor', n),
    T_K: 373.15,
    P_Pa: 101325,
    sourceUnitId: 'U2',
  };

  // ── Flash Liquid ──
  const flashLiquid: MaterialStream = {
    ...createDefaultStream('S4', 'Flash Liquid', n),
    T_K: 373.15,
    P_Pa: 101325,
    sourceUnitId: 'U2',
  };

  // ── Units ──
  const heater: UnitOperation = {
    id: 'U1',
    name: 'Heater',
    type: 'Heater',
    ports: [
      { id: 'U1-in', label: 'Inlet', type: 'inlet', position: 'left', streamId: 'S1' },
      { id: 'U1-out', label: 'Outlet', type: 'outlet', position: 'right', streamId: 'S2' },
    ],
    params: { targetT_K: 368.15, outletP_Pa: 101325 },
    solved: false,
    duty_W: 0,
    errors: [],
  };

  const flash: UnitOperation = {
    id: 'U2',
    name: 'Flash Drum',
    type: 'Flash',
    ports: [
      { id: 'U2-in', label: 'Inlet', type: 'inlet', position: 'left', streamId: 'S2' },
      { id: 'U2-vap', label: 'Vapor', type: 'outlet', position: 'top', streamId: 'S3' },
      { id: 'U2-liq', label: 'Liquid', type: 'outlet', position: 'bottom', streamId: 'S4' },
    ],
    params: { T_K: 368.15, P_Pa: 101325 },
    solved: false,
    duty_W: 0,
    errors: [],
  };

  return {
    streams: { S1: feed, S2: heaterOut, S3: flashVapor, S4: flashLiquid },
    units: { U1: heater, U2: flash },
  };
}

// ────────────────────────────────────────────────────────────────
// Store
// ────────────────────────────────────────────────────────────────

export const useCanopyStore = create<CanopyState>((set, get) => ({
  // Navigation
  currentScreen: 'setup',
  setScreen: (screen) => set({ currentScreen: screen }),

  // Setup
  compounds: [],
  fluidPackage: 'Ideal',
  unitSystem: { ...DEFAULT_UNIT_SYSTEM },
  setCompounds: (compounds) => set({ compounds }),
  addCompound: (compound) => set(s => ({
    compounds: s.compounds.some(c => c.name === compound.name)
      ? s.compounds
      : [...s.compounds, compound],
  })),
  removeCompound: (name) => set(s => ({
    compounds: s.compounds.filter(c => c.name !== name),
  })),
  setFluidPackage: (pkg) => set({ fluidPackage: pkg }),
  setUnitSystem: (patch) => set(s => ({ unitSystem: { ...s.unitSystem, ...patch } })),

  // Flowsheet
  streams: {},
  units: {},
  addStream: (stream) => set(s => ({ streams: { ...s.streams, [stream.id]: stream } })),
  updateStream: (id, patch) => set(s => ({
    streams: { ...s.streams, [id]: { ...s.streams[id], ...patch } },
  })),
  removeStream: (id) => set(s => {
    const { [id]: _, ...rest } = s.streams;
    return { streams: rest };
  }),
  addUnit: (unit) => set(s => ({ units: { ...s.units, [unit.id]: unit } })),
  updateUnit: (id, patch) => set(s => ({
    units: { ...s.units, [id]: { ...s.units[id], ...patch } },
  })),
  removeUnit: (id) => set(s => {
    const { [id]: _, ...rest } = s.units;
    return { units: rest };
  }),

  // Add unit from palette
  addUnitFromPalette: (type, position) => {
    const state = get();
    const unitId = state.nextUnitId();
    const portDefs: Record<UnitOpType, () => UnitPort[]> = {
      Heater: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Flash: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-vap`, label: 'Vapor', type: 'outlet', position: 'top', streamId: null },
        { id: `${unitId}-liq`, label: 'Liquid', type: 'outlet', position: 'bottom', streamId: null },
      ],
      Mixer: () => [
        { id: `${unitId}-in1`, label: 'Inlet 1', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-in2`, label: 'Inlet 2', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Valve: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out`, label: 'Outlet', type: 'outlet', position: 'right', streamId: null },
      ],
      Splitter: () => [
        { id: `${unitId}-in`, label: 'Inlet', type: 'inlet', position: 'left', streamId: null },
        { id: `${unitId}-out1`, label: 'Outlet 1', type: 'outlet', position: 'right', streamId: null },
        { id: `${unitId}-out2`, label: 'Outlet 2', type: 'outlet', position: 'right', streamId: null },
      ],
    };
    const defaultParams: Record<UnitOpType, Record<string, number>> = {
      Heater: { targetT_K: 373.15, outletP_Pa: 101325 },
      Flash: { T_K: 373.15, P_Pa: 101325 },
      Mixer: { outletP_Pa: 101325 },
      Valve: { outletP_Pa: 50000 },
      Splitter: { splitRatio: 0.5, outletP_Pa: 101325 },
    };
    const newUnit: UnitOperation = {
      id: unitId,
      name: `${type} ${unitId.slice(1)}`,
      type,
      ports: portDefs[type](),
      params: defaultParams[type],
      solved: false,
      duty_W: 0,
      errors: [],
    };
    set(s => ({
      units: { ...s.units, [unitId]: newUnit },
      nodePositions: { ...s.nodePositions, [unitId]: position },
      solved: false,
    }));
  },

  // Add feed stream
  addFeedStream: () => {
    const state = get();
    const streamId = state.nextStreamId();
    const n = state.compounds.length;
    const newStream: MaterialStream = {
      ...createDefaultStream(streamId, `Feed ${streamId.slice(1)}`, n),
      T_K: 298.15,
      P_Pa: 101325,
      totalFlow_molps: 100,
      moleFractions: n > 0 ? new Array(n).fill(1 / n) : [],
    };
    set(s => ({
      streams: { ...s.streams, [streamId]: newStream },
      nodePositions: { ...s.nodePositions, [`label-${streamId}`]: { x: 50, y: 200 } },
      solved: false,
    }));
  },

  nextStreamId: () => {
    const ids = Object.keys(get().streams).map(id => parseInt(id.slice(1))).filter(n => !isNaN(n));
    return `S${Math.max(0, ...ids) + 1}`;
  },
  nextUnitId: () => {
    const ids = Object.keys(get().units).map(id => parseInt(id.slice(1))).filter(n => !isNaN(n));
    return `U${Math.max(0, ...ids) + 1}`;
  },

  // Tabs (supports both units and streams)
  openTabs: [],
  activeTab: null,
  openTab: (item) => set(s => ({
    openTabs: s.openTabs.some(t => t.type === item.type && t.id === item.id)
      ? s.openTabs
      : [...s.openTabs, item],
    activeTab: item,
  })),
  closeTab: (item) => set(s => {
    const newTabs = s.openTabs.filter(t => !(t.type === item.type && t.id === item.id));
    const isActive = s.activeTab?.type === item.type && s.activeTab?.id === item.id;
    return {
      openTabs: newTabs,
      activeTab: isActive ? (newTabs[newTabs.length - 1] ?? null) : s.activeTab,
    };
  }),
  setActiveTab: (item) => set({ activeTab: item }),

  // Solver
  solved: false,
  solverErrors: [],
  solveAll: () => {
    const { compounds, streams, units } = get();
    if (compounds.length === 0) {
      set({ solverErrors: ['No compounds selected.'], solved: false });
      return;
    }

    const updatedStreams = { ...streams };
    const updatedUnits = { ...units };
    const errors: string[] = [];

    // First: flash the feed stream(s) — any stream with no source unit
    for (const [id, stream] of Object.entries(updatedStreams)) {
      if (!stream.sourceUnitId) {
        const flash = flashPT_Ideal(compounds, stream.moleFractions, stream.T_K, stream.P_Pa);
        updatedStreams[id] = {
          ...stream,
          phase: flash.phase,
          vaporFraction: flash.vaporFraction,
          H_Jpmol: flash.H_Jpmol,
          x_liquid: flash.x,
          y_vapor: flash.y,
          solved: true,
        };
      }
    }

    // Topological sort: solve units in order
    // Simple approach: iterate units, solve if all inlets are solved
    const solvedUnits = new Set<string>();
    let changed = true;
    let maxPasses = 10;

    while (changed && maxPasses-- > 0) {
      changed = false;
      for (const [unitId, unit] of Object.entries(updatedUnits)) {
        if (solvedUnits.has(unitId)) continue;

        // Check if all inlet streams are solved
        const inletPorts = unit.ports.filter(p => p.type === 'inlet');
        const inletStreams = inletPorts
          .map(p => p.streamId ? updatedStreams[p.streamId] : null)
          .filter((s): s is MaterialStream => s !== null && s.solved);

        if (inletStreams.length < inletPorts.length) continue; // Not all inlets ready

        try {
          switch (unit.type) {
            case 'Heater': {
              const result = solveHeater(compounds, inletStreams[0], {
                targetT_K: unit.params.targetT_K as number,
                outletP_Pa: unit.params.outletP_Pa as number,
              });
              const outPort = unit.ports.find(p => p.type === 'outlet');
              if (outPort?.streamId) {
                updatedStreams[outPort.streamId] = {
                  ...updatedStreams[outPort.streamId],
                  ...result.outlet,
                  id: outPort.streamId,
                  name: updatedStreams[outPort.streamId].name,
                  sourceUnitId: unitId,
                  sourcePortId: outPort.id,
                  targetUnitId: updatedStreams[outPort.streamId].targetUnitId,
                  targetPortId: updatedStreams[outPort.streamId].targetPortId,
                } as MaterialStream;
              }
              updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
              break;
            }

            case 'Flash': {
              const result = solveFlashDrum(compounds, inletStreams[0], {
                T_K: unit.params.T_K as number,
                P_Pa: unit.params.P_Pa as number,
              });
              const vapPort = unit.ports.find(p => p.id.includes('vap') || p.label === 'Vapor');
              const liqPort = unit.ports.find(p => p.id.includes('liq') || p.label === 'Liquid');
              if (vapPort?.streamId) {
                updatedStreams[vapPort.streamId] = {
                  ...updatedStreams[vapPort.streamId],
                  ...result.vapor,
                  id: vapPort.streamId,
                  name: updatedStreams[vapPort.streamId].name,
                  sourceUnitId: unitId,
                  sourcePortId: vapPort.id,
                  targetUnitId: updatedStreams[vapPort.streamId].targetUnitId,
                  targetPortId: updatedStreams[vapPort.streamId].targetPortId,
                } as MaterialStream;
              }
              if (liqPort?.streamId) {
                updatedStreams[liqPort.streamId] = {
                  ...updatedStreams[liqPort.streamId],
                  ...result.liquid,
                  id: liqPort.streamId,
                  name: updatedStreams[liqPort.streamId].name,
                  sourceUnitId: unitId,
                  sourcePortId: liqPort.id,
                  targetUnitId: updatedStreams[liqPort.streamId].targetUnitId,
                  targetPortId: updatedStreams[liqPort.streamId].targetPortId,
                } as MaterialStream;
              }
              updatedUnits[unitId] = { ...unit, duty_W: result.duty_W, solved: true, errors: [] };
              break;
            }

            case 'Mixer': {
              const result = solveMixer(compounds, inletStreams,
                (unit.params.outletP_Pa as number) ?? inletStreams[0].P_Pa);
              const outPort = unit.ports.find(p => p.type === 'outlet');
              if (outPort?.streamId) {
                updatedStreams[outPort.streamId] = {
                  ...updatedStreams[outPort.streamId],
                  ...result.outlet,
                  id: outPort.streamId,
                  name: updatedStreams[outPort.streamId].name,
                  sourceUnitId: unitId,
                  sourcePortId: outPort.id,
                  targetUnitId: updatedStreams[outPort.streamId].targetUnitId,
                  targetPortId: updatedStreams[outPort.streamId].targetPortId,
                } as MaterialStream;
              }
              updatedUnits[unitId] = { ...unit, duty_W: 0, solved: true, errors: [] };
              break;
            }

            case 'Valve': {
              const result = solveValve(compounds, inletStreams[0], {
                outletP_Pa: unit.params.outletP_Pa as number,
              });
              const outPort = unit.ports.find(p => p.type === 'outlet');
              if (outPort?.streamId) {
                updatedStreams[outPort.streamId] = {
                  ...updatedStreams[outPort.streamId],
                  ...result.outlet,
                  id: outPort.streamId,
                  name: updatedStreams[outPort.streamId].name,
                  sourceUnitId: unitId,
                  sourcePortId: outPort.id,
                  targetUnitId: updatedStreams[outPort.streamId].targetUnitId,
                  targetPortId: updatedStreams[outPort.streamId].targetPortId,
                } as MaterialStream;
              }
              updatedUnits[unitId] = { ...unit, duty_W: 0, solved: true, errors: [] };
              break;
            }

            default:
              errors.push(`Unknown unit type: ${unit.type}`);
              continue;
          }
          solvedUnits.add(unitId);
          changed = true;
        } catch (e: any) {
          errors.push(`Error solving ${unit.name}: ${e.message}`);
          updatedUnits[unitId] = { ...unit, solved: false, errors: [e.message] };
        }
      }
    }

    set({
      streams: updatedStreams,
      units: updatedUnits,
      solved: errors.length === 0 && Object.values(updatedUnits).every(u => u.solved),
      solverErrors: errors,
    });
  },

  // Preload
  loadPreset: (preset) => {
    if (preset === 'benzene-toluene') {
      const { streams, units } = createPresetFlowsheet();
      set({
        compounds: [BENZENE, TOLUENE],
        fluidPackage: 'Ideal',
        unitSystem: { ...DEFAULT_UNIT_SYSTEM },
        streams,
        units,
        currentScreen: 'setup',
        solved: false,
        solverErrors: [],
        openTabs: [],
        activeTab: null,
        nodePositions: {
          'label-S1': { x: 50, y: 200 },
          'U1': { x: 250, y: 200 },
          'U2': { x: 550, y: 200 },
          'label-S3': { x: 750, y: 50 },
          'label-S4': { x: 750, y: 380 },
        },
      });
    }
  },

  // Node positions
  nodePositions: {},
  setNodePosition: (nodeId, pos) => set(s => ({
    nodePositions: { ...s.nodePositions, [nodeId]: pos },
  })),
}));
