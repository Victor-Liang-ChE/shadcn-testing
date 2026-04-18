'use client';

import { useCallback, useMemo, useState, useRef, useEffect } from 'react';
import {
  ReactFlow,
  Background,
  Controls,
  MiniMap,
  ConnectionLineType,
  type Edge,
  type Node,
  type NodeChange,
  type EdgeChange,
  type Connection,
  MarkerType,
  Panel,
  applyNodeChanges,
  applyEdgeChanges,
  useReactFlow,
  ReactFlowProvider,
} from '@xyflow/react';
import '@xyflow/react/dist/style.css';

/* Selected-edge highlight — orange glow */
const selectedEdgeCSS = `
  .react-flow__edge.selected path,
  .react-flow__edge:focus path,
  .react-flow__edge:focus-visible path {
    stroke: #f97316 !important;
    stroke-width: 3px !important;
    filter: drop-shadow(0 0 4px rgba(249,115,22,0.5));
  }
`;

import { useCanopyStore, type TabItem } from '@/lib/canopy/store';
import { canopyNodeTypes } from './FlowNodes';
import {
  HeaterIcon, FlashDrumIcon, MixerIcon, SplitterIcon, ValveIcon,
  PumpIcon, CompressorIcon, HeatExchangerIcon, ColumnIcon, CSTRIcon, PFRIcon,
  PipeIcon, AbsorberIcon, SeparatorIcon, DecanterIcon, CrystallizerIcon,
  CrusherIcon, DryerIcon, ReactorIcon, RGibbsIcon, RYieldIcon, REquilIcon, RBatchIcon,
} from './UnitIcons';
import { X, ArrowRight, Trash2, Unlink, Undo2, Redo2, Download, Upload, Pause, Play, FileText } from 'lucide-react';
import { Button } from '@/components/ui/button';
import { Accordion, AccordionContent, AccordionItem, AccordionTrigger } from '@/components/ui/accordion';
import {
  fromSI_T, toSI_T, fromSI_P, toSI_P, fromSI_F, toSI_F, fromSI_Q,
  tempLabel, pressLabel, flowLabel, energyLabel,
} from '@/lib/canopy/units';
import type { UnitOpType, MaterialStream, CanopyCompound } from '@/lib/canopy/types';
import { calculateEconomics, calculateEquipmentCost, DEFAULT_UTILITY_RATES } from '@/lib/canopy/economics';
import type { EconomicsSummary, UtilityRates } from '@/lib/canopy/economics';
import { generateTxy, generatePxy, generateXY } from '@/lib/canopy/analysis';
import type { TxyData, PxyData, XYData } from '@/lib/canopy/analysis';
import { getRigorousColumnSpecDiagnostics, type RigorousColumnSpecMode } from '@/lib/canopy/solver';

function fmt(val: number, decimals = 2): string {
  if (!isFinite(val)) return '—';
  return val.toFixed(decimals);
}

/** Extract stream ID from edge id patterns: e-feed-{id}, e-prod-{id}, e-{id}-a, e-{id}-b, e-{id} */
function extractStreamIdFromEdge(edgeId: string): string | null {
  if (edgeId.startsWith('e-feed-')) return edgeId.slice(7);
  if (edgeId.startsWith('e-prod-')) return edgeId.slice(7);
  // e-{streamId}-a or e-{streamId}-b (intermediate split edges)
  const m = edgeId.match(/^e-(.+)-[ab]$/);
  if (m) return m[1];
  // Legacy: e-{streamId}
  if (edgeId.startsWith('e-')) return edgeId.slice(2);
  return null;
}

function streamSummary(s: MaterialStream, us: ReturnType<typeof useCanopyStore.getState>['unitSystem']): string {
  const lines = [
    `T: ${fmt(fromSI_T(s.T_K, us.temperature), 1)}${tempLabel(us.temperature)}`,
    `P: ${fmt(fromSI_P(s.P_Pa, us.pressure), 2)} ${pressLabel(us.pressure)}`,
    `Flow: ${fmt(fromSI_F(s.totalFlow_molps, us.flow), 1)} ${flowLabel(us.flow)}`,
    `Phase: ${s.phase}`,
  ];
  if (s.vaporFraction > 0 && s.vaporFraction < 1) {
    lines.push(`VF: ${fmt(s.vaporFraction * 100, 1)}%`);
  }
  return lines.join('\n');
}

// ─── Object Palette ─────────────────────────────────────────────
// Each item is draggable onto the canvas (HTML5 drag-and-drop).
// Also supports click-to-add at a default position.

type PaletteItem = { type: UnitOpType; label: string; Icon: React.ComponentType<{ size?: number; className?: string }> };
type PaletteGroup = { id: string; label: string; items: PaletteItem[] };

const PALETTE_GROUPS: PaletteGroup[] = [
  {
    id: 'columns',
    label: 'Columns',
    items: [
      { type: 'Column', label: 'Column', Icon: ColumnIcon },
      { type: 'RadFrac', label: 'RadFrac', Icon: ColumnIcon },
      { type: 'RateFrac', label: 'RateFrac', Icon: ColumnIcon },
      { type: 'Absorber', label: 'Absorber', Icon: AbsorberIcon },
      { type: 'Extractor', label: 'Extractor', Icon: AbsorberIcon },
    ],
  },
  {
    id: 'reactors',
    label: 'Reactors',
    items: [
      { type: 'RStoic', label: 'RStoic', Icon: ReactorIcon },
      { type: 'CSTR', label: 'CSTR', Icon: CSTRIcon },
      { type: 'PFR', label: 'PFR', Icon: PFRIcon },
      { type: 'RYield', label: 'RYield', Icon: RYieldIcon },
      { type: 'REquil', label: 'REquil', Icon: REquilIcon },
      { type: 'RGibbs', label: 'RGibbs', Icon: RGibbsIcon },
      { type: 'RBatch', label: 'RBatch', Icon: RBatchIcon },
      { type: 'ElectrolyteEquilibrium', label: 'Elec Equil', Icon: REquilIcon },
    ],
  },
  {
    id: 'heating',
    label: 'Heat Exchangers',
    items: [
      { type: 'Heater', label: 'Heater', Icon: HeaterIcon },
      { type: 'HeatX', label: 'Heat Exchanger', Icon: HeatExchangerIcon },
      { type: 'MHeatX', label: 'Multi-HeatX', Icon: HeatExchangerIcon },
      { type: 'Flash', label: 'Flash Drum', Icon: FlashDrumIcon },
      { type: 'ThreePhaseFlash', label: 'Flash3', Icon: FlashDrumIcon },
      { type: 'Decanter', label: 'Decanter', Icon: DecanterIcon },
    ],
  },
  {
    id: 'pressure',
    label: 'Pressure Changers',
    items: [
      { type: 'Valve', label: 'Valve', Icon: ValveIcon },
      { type: 'Pump', label: 'Pump', Icon: PumpIcon },
      { type: 'Compressor', label: 'Compressor', Icon: CompressorIcon },
      { type: 'Pipe', label: 'Pipe', Icon: PipeIcon },
    ],
  },
  {
    id: 'mixsep',
    label: 'Mixers And Separators',
    items: [
      { type: 'Mixer', label: 'Mixer', Icon: MixerIcon },
      { type: 'Splitter', label: 'Splitter', Icon: SplitterIcon },
      { type: 'ComponentSeparator', label: 'Sep', Icon: SeparatorIcon },
      { type: 'Membrane', label: 'Membrane', Icon: SeparatorIcon },
      { type: 'Cyclone', label: 'Cyclone', Icon: SeparatorIcon },
      { type: 'Filter', label: 'Filter', Icon: SeparatorIcon },
      { type: 'Screen', label: 'Screen', Icon: SeparatorIcon },
      { type: 'Centrifuge', label: 'Centrifuge', Icon: SeparatorIcon },
    ],
  },
  {
    id: 'solids',
    label: 'Solids',
    items: [
      { type: 'Crystallizer', label: 'Crystallizer', Icon: CrystallizerIcon },
      { type: 'Crusher', label: 'Crusher', Icon: CrusherIcon },
      { type: 'Dryer', label: 'Dryer', Icon: DryerIcon },
    ],
  },
  {
    id: 'utilities',
    label: 'Utilities',
    items: [
      { type: 'SteamDrum', label: 'Steam Drum', Icon: FlashDrumIcon },
      { type: 'SteamHeater', label: 'Steam Heater', Icon: HeaterIcon },
      { type: 'SteamTurbine', label: 'Steam Turbine', Icon: CompressorIcon },
      { type: 'SteamValve', label: 'Steam Valve', Icon: ValveIcon },
      { type: 'SteamHeader', label: 'Steam Header', Icon: MixerIcon },
      { type: 'SteamTrap', label: 'Steam Trap', Icon: FlashDrumIcon },
    ],
  },
  {
    id: 'control',
    label: 'Control And Analysis',
    items: [
      { type: 'PIDController', label: 'PID Ctrl', Icon: ValveIcon },
      { type: 'LeadLag', label: 'Lead-Lag', Icon: ValveIcon },
      { type: 'DeadTime', label: 'Dead Time', Icon: ValveIcon },
      { type: 'SignalSelector', label: 'Selector', Icon: ValveIcon },
      { type: 'Analyzer', label: 'Analyzer', Icon: FileText },
      { type: 'ChargeBalance', label: 'Charge Bal', Icon: SeparatorIcon },
    ],
  },
];

function ObjectPalette() {
  const addUnitFromPalette = useCanopyStore(s => s.addUnitFromPalette);
  const addFeedStream = useCanopyStore(s => s.addFeedStream);
  const [openGroup, setOpenGroup] = useState<string>('columns');

  const onDragStart = (e: React.DragEvent, type: UnitOpType | 'feed') => {
    e.dataTransfer.setData('canopy/object-type', type);
    e.dataTransfer.effectAllowed = 'copy';
  };

  return (
    <div
      className="pointer-events-auto flex h-[calc(100dvh-8rem)] max-h-[calc(100dvh-8rem)] w-44 flex-col overflow-hidden rounded-lg border border-border bg-card/95 p-2 shadow-lg backdrop-blur"
      style={{ touchAction: 'pan-y' }}
      onWheelCapture={(event) => event.stopPropagation()}
    >
      <div className="text-[10px] font-semibold text-muted-foreground uppercase tracking-wider px-1 mb-1">
        Objects
      </div>
      <div className="text-[10px] text-muted-foreground px-1 mb-2">
        Drag or click to add
      </div>
      <div
        className="flex min-h-0 flex-1 flex-col gap-1 overflow-y-auto overscroll-contain pr-1"
        onWheel={(event) => event.stopPropagation()}
        style={{ touchAction: 'pan-y' }}
      >
        <div
          draggable
          onDragStart={(e) => onDragStart(e, 'feed')}
          onClick={() => addFeedStream({ x: 50, y: 200 })}
          className="mb-2 flex items-center gap-2 rounded border border-border/60 px-2 py-1.5 text-xs transition-colors hover:bg-accent cursor-grab active:cursor-grabbing select-none"
          title="Drag or click to add Stream"
        >
          <ArrowRight className="h-3 w-3 text-blue-500" />
          <span>Material Stream</span>
        </div>
        <Accordion type="single" collapsible value={openGroup} onValueChange={setOpenGroup} className="w-full">
          {PALETTE_GROUPS.map(group => (
            <AccordionItem key={group.id} value={group.id} className="border-b border-border/60">
              <AccordionTrigger className="py-2 text-[11px] font-semibold text-muted-foreground hover:no-underline">
                {group.label}
              </AccordionTrigger>
              <AccordionContent className="pb-2">
                <div className="flex flex-col gap-1">
                  {group.items.map(item => (
                    <div
                      key={item.type}
                      draggable
                      onDragStart={(e) => onDragStart(e, item.type)}
                      onClick={() => addUnitFromPalette(item.type, { x: 350, y: 200 })}
                      className="flex items-center gap-2 rounded px-2 py-1.5 text-xs transition-colors hover:bg-accent cursor-grab active:cursor-grabbing select-none"
                      title={`Drag or click to add ${item.label}`}
                    >
                      <item.Icon size={14} className="shrink-0 opacity-70" />
                      <span>{item.label}</span>
                    </div>
                  ))}
                </div>
              </AccordionContent>
            </AccordionItem>
          ))}
        </Accordion>
      </div>
    </div>
  );
}

// ─── Main Flowsheet Screen ──────────────────────────────────────
// Wrapped in ReactFlowProvider to allow useReactFlow() inside.

export default function FlowsheetScreen() {
  return (
    <ReactFlowProvider>
      <FlowsheetInner />
    </ReactFlowProvider>
  );
}

function FlowsheetInner() {
  const {
    compounds, streams, units, solved, unitSystem,
    openTabs, activeTab, openTab, closeTab, setActiveTab,
    nodePositions, setNodePosition, addUnitFromPalette, addFeedStream, removeUnit, removeStream,
    disconnectStream, connectPorts, updateStream,
    solverPaused, setSolverPaused,
    solverErrors, dbWarnings,
  } = useCanopyStore();

  const { screenToFlowPosition } = useReactFlow();

  // ── Auto-solve: re-run solver when streams/units/compounds change (debounced) ──
  const solveTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  const streamsRef = useRef(streams);
  const unitsRef = useRef(units);
  const compoundsRef = useRef(compounds);
  useEffect(() => {
    // Only auto-solve if something actually changed and we have data
    const changed = streamsRef.current !== streams || unitsRef.current !== units || compoundsRef.current !== compounds;
    streamsRef.current = streams;
    unitsRef.current = units;
    compoundsRef.current = compounds;
    if (!changed || compounds.length === 0) return;

    if (solveTimerRef.current) clearTimeout(solveTimerRef.current);
    solveTimerRef.current = setTimeout(() => {
      useCanopyStore.getState().solveAll();
    }, 400);
    return () => { if (solveTimerRef.current) clearTimeout(solveTimerRef.current); };
  }, [streams, units, compounds]);

  // ── Derive the "source of truth" nodes from store (re-computed when store changes) ──
  const storeNodes: Node[] = useMemo(() => {
    const result: Node[] = [];
    const feedStreams = Object.values(streams).filter(s => !s.sourceUnitId);
    feedStreams.forEach((s, i) => {
      const pos = nodePositions[`label-${s.id}`] ?? { x: 50, y: 200 + i * 120 };
      result.push({
        id: `label-${s.id}`,
        type: 'streamLabel',
        position: pos,
        data: {
          label: s.name,
          isSource: true,
          streamId: s.id,
          streamInfo: s.solved ? streamSummary(s, unitSystem) : '',
        },
        draggable: true,
      });
    });
    const unitEntries = Object.values(units);
    unitEntries.forEach((u, i) => {
      const nodeType = u.type.toLowerCase() as string;
      const pos = nodePositions[u.id] ?? { x: 250 + i * 300, y: 200 };
      const inletPorts = u.ports.filter(p => p.type === 'inlet');
      const outletPorts = u.ports.filter(p => p.type === 'outlet');
      result.push({
        id: u.id,
        type: nodeType === 'flash' ? 'flash' : nodeType,
        position: pos,
        data: {
          label: u.name,
          duty: u.solved ? u.duty_W : null,
          unitSystem,
          inletCount: inletPorts.length,
          outletCount: outletPorts.length,
          propertyCode: u.params.propertyCode,
          propertyLabel: u.params.propertyLabel,
          signalValue: u.params.signalValue,
        },
        draggable: true,
      });
    });
    const productStreams = Object.values(streams).filter(s => !s.targetUnitId && s.sourceUnitId);
    productStreams.forEach((s, i) => {
      const pos = nodePositions[`label-${s.id}`] ?? { x: 700, y: 50 + i * 160 };
      result.push({
        id: `label-${s.id}`,
        type: 'streamLabel',
        position: pos,
        data: {
          label: s.name,
          isSource: false,
          streamId: s.id,
          streamInfo: s.solved ? streamSummary(s, unitSystem) : '',
        },
        draggable: true,
      });
    });

    // Intermediate stream labels (connecting two units) — draggable nodes
    const intermediateStreams = Object.values(streams).filter(s => s.sourceUnitId && s.targetUnitId);
    intermediateStreams.forEach((s) => {
      // Use stored position; only compute midpoint as initial fallback
      const pos = nodePositions[`label-${s.id}`] ?? (() => {
        const srcPos = nodePositions[s.sourceUnitId!] ?? { x: 300, y: 200 };
        const tgtPos = nodePositions[s.targetUnitId!] ?? { x: 600, y: 200 };
        return {
          x: (srcPos.x + tgtPos.x) / 2,
          y: (srcPos.y + tgtPos.y) / 2 - 30,
        };
      })();
      result.push({
        id: `label-${s.id}`,
        type: 'streamLabel',
        position: pos,
        data: {
          label: s.name,
          isSource: false,
          isIntermediate: true,
          streamId: s.id,
          streamInfo: s.solved ? streamSummary(s, unitSystem) : '',
        },
        draggable: true,
      });
    });

    return result;
  }, [streams, units, compounds, solved, nodePositions, unitSystem]);

  // Persist initial positions for intermediate labels that lack a stored position
  // (deferred to useEffect to avoid setState-during-render)
  useEffect(() => {
    const intermediateStreams = Object.values(streams).filter(s => s.sourceUnitId && s.targetUnitId);
    for (const s of intermediateStreams) {
      const key = `label-${s.id}`;
      if (!nodePositions[key]) {
        const srcPos = nodePositions[s.sourceUnitId!] ?? { x: 300, y: 200 };
        const tgtPos = nodePositions[s.targetUnitId!] ?? { x: 600, y: 200 };
        setNodePosition(key, {
          x: (srcPos.x + tgtPos.x) / 2,
          y: (srcPos.y + tgtPos.y) / 2 - 30,
        });
      }
    }
  }, [streams, nodePositions, setNodePosition]);

  // ── Local React Flow node state for fluid dragging ──
  // We keep a local copy that React Flow mutates during drag.
  // We sync from storeNodes whenever the store changes (excluding drag).
  const [rfNodes, setRfNodes] = useState<Node[]>(storeNodes);
  const isDragging = useRef(false);

  // Sync store → local when store changes (but not while dragging)
  useEffect(() => {
    if (!isDragging.current) {
      setRfNodes(storeNodes);
    }
  }, [storeNodes]);

  // Handle React Flow node changes (position, selection, etc.)
  const onNodesChange = useCallback((changes: NodeChange[]) => {
    // Apply changes to local state immediately (fluid drag animation)
    setRfNodes(prev => applyNodeChanges(changes, prev));

    for (const change of changes) {
      if (change.type === 'position') {
        if (change.dragging) {
          isDragging.current = true;
        } else if (isDragging.current && change.position) {
          // Drag ended — persist final position to store
          isDragging.current = false;
          setNodePosition(change.id, change.position);
        }
      }
    }
  }, [setNodePosition]);

  // ── Edges ──
  const edges: Edge[] = useMemo(() => {
    const result: Edge[] = [];
    for (const u of Object.values(units)) {
      for (const port of u.ports) {
        if (port.type === 'inlet' && port.streamId) {
          const stream = streams[port.streamId];
          if (stream && !stream.sourceUnitId) {
            result.push({
              id: `e-feed-${port.streamId}`,
              source: `label-${port.streamId}`,
              sourceHandle: 'out',
              target: u.id,
              targetHandle: port.id.slice(u.id.length + 1) || 'in',
              type: 'smoothstep',
              style: { stroke: '#3b82f6', strokeWidth: 2 },
              markerEnd: { type: MarkerType.ArrowClosed, color: '#3b82f6' },
              label: stream.name,
              labelStyle: { fontSize: 10 },
            });
          } else if (stream?.sourceUnitId) {
            const srcUnit = units[stream.sourceUnitId];
            const srcPort = srcUnit?.ports.find(p => p.streamId === port.streamId && p.type === 'outlet');
            if (srcUnit && srcPort) {
              // Split into two edges through draggable label node
              result.push({
                id: `e-${port.streamId}-a`,
                source: srcUnit.id,
                sourceHandle: srcPort.id.slice(srcUnit.id.length + 1) || 'out',
                target: `label-${port.streamId}`,
                targetHandle: 'in',
                type: 'smoothstep',
                style: { stroke: '#3b82f6', strokeWidth: 2 },
              });
              result.push({
                id: `e-${port.streamId}-b`,
                source: `label-${port.streamId}`,
                sourceHandle: 'out',
                target: u.id,
                targetHandle: port.id.slice(u.id.length + 1) || 'in',
                type: 'smoothstep',
                style: { stroke: '#3b82f6', strokeWidth: 2 },
                markerEnd: { type: MarkerType.ArrowClosed, color: '#3b82f6' },
              });
            }
          }
        }
      }
    }
    for (const u of Object.values(units)) {
      for (const port of u.ports) {
        if (port.type === 'outlet' && port.streamId) {
          const stream = streams[port.streamId];
          if (stream && !stream.targetUnitId) {
            const handleId = port.id.slice(u.id.length + 1) || 'out';
            result.push({
              id: `e-prod-${port.streamId}`,
              source: u.id,
              sourceHandle: handleId,
              target: `label-${port.streamId}`,
              targetHandle: 'in',
              type: 'smoothstep',
              style: { stroke: '#3b82f6', strokeWidth: 2 },
              markerEnd: {
                type: MarkerType.ArrowClosed,
                color: '#3b82f6',
              },
              label: stream.name,
              labelStyle: { fontSize: 10 },
            });
          }
        }
      }
    }
    return result;
  }, [streams, units, solved]);

  // ── Drag-and-drop from palette onto canvas ──
  const onDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.dataTransfer.dropEffect = 'copy';
  }, []);

  const onDrop = useCallback((e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    const objectType = e.dataTransfer.getData('canopy/object-type') as UnitOpType | 'feed' | '';
    if (!objectType) return;
    const pos = screenToFlowPosition({ x: e.clientX, y: e.clientY });
    if (objectType === 'feed') {
      addFeedStream(pos);
    } else {
      addUnitFromPalette(objectType, pos);
    }
  }, [screenToFlowPosition, addUnitFromPalette, addFeedStream]);

  // ── Double click: open unit OR stream tab ──
  const onNodeDoubleClick = useCallback((_event: React.MouseEvent, node: Node) => {
    if (units[node.id]) {
      openTab({ type: 'unit', id: node.id });
    } else if (node.id.startsWith('label-')) {
      const streamId = node.id.replace('label-', '');
      if (streams[streamId]) {
        openTab({ type: 'stream', id: streamId });
      }
    }
  }, [units, streams, openTab]);

  // ── Double click edge: open intermediate stream tab ──
  const onEdgeDoubleClick = useCallback((_event: React.MouseEvent, edge: Edge) => {
    const streamId = extractStreamIdFromEdge(edge.id);
    if (streamId && streams[streamId]) {
      openTab({ type: 'stream', id: streamId });
    }
  }, [streams, openTab]);

  // ── Delete selected edge (disconnect stream) ──
  const onEdgesDelete = useCallback((deletedEdges: Edge[]) => {
    for (const edge of deletedEdges) {
      const streamId = extractStreamIdFromEdge(edge.id);
      if (streamId) disconnectStream(streamId);
    }
  }, [disconnectStream]);

  // ── Connect nodes: create stream when user drags handle to handle ──
  const onConnect = useCallback((connection: Connection) => {
    const { source, target, sourceHandle, targetHandle } = connection;
    if (!source || !target) return;

    // Feed stream → unit connection
    if (source.startsWith('label-') && units[target]) {
      const streamId = source.replace('label-', '');
      const stream = streams[streamId];
      if (stream && !stream.sourceUnitId) {
        // Connect this feed stream to the target unit inlet
        const tgtUnit = units[target];
        const tgtPort = tgtUnit.ports.find(p => {
          const hId = p.id.slice(target!.length + 1) || '';
          return p.type === 'inlet' && (hId === targetHandle || p.id === targetHandle);
        });
        if (tgtPort) {
          updateStream(streamId, { targetUnitId: target, targetPortId: tgtPort.id });
          // Update port on unit
          const updatedPorts = tgtUnit.ports.map(p =>
            p.id === tgtPort.id ? { ...p, streamId } : p
          );
          useCanopyStore.getState().updateUnit(target, { ports: updatedPorts });
        }
      }
      return;
    }

    // Unit → unit connection
    if (units[source] && units[target] && sourceHandle && targetHandle) {
      connectPorts(source, sourceHandle, target, targetHandle);
    }
  }, [streams, units, connectPorts, updateStream]);

  // ── Edge context menu (right-click) ──
  const [edgeMenu, setEdgeMenu] = useState<{
    x: number; y: number; streamId: string;
  } | null>(null);

  const onEdgeContextMenu = useCallback((event: React.MouseEvent, edge: Edge) => {
    event.preventDefault();
    const streamId = extractStreamIdFromEdge(edge.id);
    if (streamId && streams[streamId]) {
      setEdgeMenu({ x: event.clientX, y: event.clientY, streamId });
    }
  }, [streams]);

  // Close edge menu on click anywhere
  useEffect(() => {
    const close = () => setEdgeMenu(null);
    window.addEventListener('click', close);
    return () => window.removeEventListener('click', close);
  }, []);

  // ── Edge changes (selection) ──
  const [rfEdges, setRfEdges] = useState<Edge[]>(edges);
  useEffect(() => { setRfEdges(edges); }, [edges]);

  // ── Bottom panel visibility (reset when solver re-runs) ──
  const [bottomPanelHidden, setBottomPanelHidden] = useState(false);
  const prevSolvedRef = useRef(solved);
  useEffect(() => {
    if (solved && !prevSolvedRef.current) setBottomPanelHidden(false);
    prevSolvedRef.current = solved;
  }, [solved]);

  const onEdgesChange = useCallback((changes: EdgeChange[]) => {
    setRfEdges(prev => applyEdgeChanges(changes, prev));
  }, []);

  // ── Delete selected node on Delete/Backspace key (units, feed streams, or edges) ──
  const onKeyDown = useCallback((e: React.KeyboardEvent) => {
    // Undo / Redo
    if ((e.ctrlKey || e.metaKey) && e.key === 'z' && !e.shiftKey) {
      e.preventDefault();
      useCanopyStore.getState().undo();
      return;
    }
    if ((e.ctrlKey || e.metaKey) && (e.key === 'y' || (e.key === 'z' && e.shiftKey))) {
      e.preventDefault();
      useCanopyStore.getState().redo();
      return;
    }
    if (e.key === 'Delete' || e.key === 'Backspace') {
      // Check for selected edge first
      const selectedEdge = rfEdges.find(edge => edge.selected);
      if (selectedEdge) {
        const streamId = extractStreamIdFromEdge(selectedEdge.id);
        if (streamId) disconnectStream(streamId);
        return;
      }
      const selected = rfNodes.find(n => n.selected);
      if (!selected) return;
      if (units[selected.id]) {
        removeUnit(selected.id);
      } else if (selected.id.startsWith('label-')) {
        const streamId = selected.id.replace('label-', '');
        const stream = streams[streamId];
        if (stream && !stream.sourceUnitId) {
          removeStream(streamId);
        }
      }
    }
  }, [rfNodes, rfEdges, units, streams, removeUnit, removeStream, disconnectStream]);

  // ── Selected edge for toolbar delete button ──
  const selectedEdgeStreamId = (() => {
    const sel = rfEdges.find(edge => edge.selected);
    if (!sel) return undefined;
    return extractStreamIdFromEdge(sel.id) ?? undefined;
  })();

  // ── Selected node for toolbar delete button ──
  const selectedUnitId = rfNodes.find(n => n.selected && units[n.id])?.id;
  const selectedFeedId = (() => {
    const sel = rfNodes.find(n => n.selected && n.id.startsWith('label-'));
    if (!sel) return undefined;
    const sid = sel.id.replace('label-', '');
    const s = streams[sid];
    return s && !s.sourceUnitId ? sid : undefined;
  })();

  const tabLabel = (tab: TabItem) => {
    if (tab.type === 'unit') return units[tab.id]?.name ?? tab.id;
    return streams[tab.id]?.name ?? tab.id;
  };

  const isActiveTab = (tab: TabItem) =>
    activeTab?.type === tab.type && activeTab?.id === tab.id;

  return (
    <div className="flex-1 flex flex-col overflow-hidden">
      {/* Tab bar */}
      {openTabs.length > 0 && (
        <div className="flex items-center border-b border-border bg-muted/30 px-2 overflow-x-auto">
          <button
            onClick={() => setActiveTab(null)}
            className={`px-3 py-1.5 text-xs font-medium border-b-2 transition-colors mr-1
              ${activeTab === null ? 'border-blue-500 text-foreground' : 'border-transparent text-muted-foreground hover:text-foreground'}`}
          >
            PFD
          </button>
          {openTabs.map(tab => (
            <div key={`${tab.type}-${tab.id}`} className="flex items-center">
              <button
                onClick={() => setActiveTab(tab)}
                className={`px-3 py-1.5 text-xs font-medium border-b-2 transition-colors
                  ${isActiveTab(tab) ? 'border-blue-500 text-foreground' : 'border-transparent text-muted-foreground hover:text-foreground'}`}
              >
                {tab.type === 'stream' ? '📊 ' : ''}{tabLabel(tab)}
              </button>
              <button
                onClick={() => closeTab(tab)}
                className="p-0.5 text-muted-foreground hover:text-foreground transition-colors"
              >
                <X className="w-3 h-3" />
              </button>
            </div>
          ))}
        </div>
      )}

      {/* Content area */}
      {activeTab === null ? (
        <div
          className="flex-1 relative"
          onDragOver={onDragOver}
          onDrop={onDrop}
          onKeyDown={onKeyDown}
          tabIndex={0}
          style={{ outline: 'none' }}
        >
          <style dangerouslySetInnerHTML={{ __html: selectedEdgeCSS }} />
          <ReactFlow
            nodes={rfNodes}
            edges={rfEdges}
            nodeTypes={canopyNodeTypes}
            onNodeDoubleClick={onNodeDoubleClick}
            onEdgeDoubleClick={onEdgeDoubleClick}
            onEdgesDelete={onEdgesDelete}
            onEdgesChange={onEdgesChange}
            onEdgeContextMenu={onEdgeContextMenu}
            onNodesChange={onNodesChange}
            onConnect={onConnect}
            connectionLineType={ConnectionLineType.SmoothStep}
            defaultEdgeOptions={{
              type: 'smoothstep',
              focusable: true,
              style: { strokeWidth: 2, cursor: 'pointer' },
            }}
            fitView
            proOptions={{ hideAttribution: true }}
            minZoom={0.3}
            maxZoom={2}
            nodesDraggable={true}
            edgesReconnectable={false}
            deleteKeyCode="Delete"
            elevateNodesOnSelect={true}
          >
            <Background gap={20} size={1} />
            <Controls className="!bg-card !border-border !shadow-md [&>button]:!bg-card [&>button]:!border-border [&>button]:!text-foreground [&>button:hover]:!bg-accent" />
            <MiniMap
              className="!bg-card !border !border-border !shadow-md"
              maskColor="rgba(128,128,128,0.15)"
              nodeColor={(n) => {
                if (n.type === 'heater') return '#f97316';
                if (n.type === 'flash') return '#a855f7';
                if (n.type === 'mixer') return '#22c55e';
                return '#6b7280';
              }}
            />
            <Panel position="top-left">
              <ObjectPalette />
            </Panel>
            <Panel position="top-right">
              <div className="flex flex-col gap-2 items-end">
                {/* Undo / Redo / Export / Import */}
                <div className="flex items-center gap-1">
                  <button
                    onClick={() => useCanopyStore.getState().undo()}
                    disabled={!useCanopyStore.getState().canUndo}
                    className="p-1.5 text-xs bg-card/90 border border-border rounded shadow hover:bg-accent disabled:opacity-30 disabled:cursor-not-allowed"
                    title="Undo (Ctrl+Z)"
                  >
                    <Undo2 className="w-3.5 h-3.5" />
                  </button>
                  <button
                    onClick={() => useCanopyStore.getState().redo()}
                    disabled={!useCanopyStore.getState().canRedo}
                    className="p-1.5 text-xs bg-card/90 border border-border rounded shadow hover:bg-accent disabled:opacity-30 disabled:cursor-not-allowed"
                    title="Redo (Ctrl+Y)"
                  >
                    <Redo2 className="w-3.5 h-3.5" />
                  </button>
                  <button
                    onClick={() => {
                      const json = useCanopyStore.getState().exportFlowsheet();
                      const blob = new Blob([json], { type: 'application/json' });
                      const url = URL.createObjectURL(blob);
                      const a = document.createElement('a');
                      a.href = url;
                      a.download = 'canopy-flowsheet.json';
                      a.click();
                      URL.revokeObjectURL(url);
                    }}
                    className="p-1.5 text-xs bg-card/90 border border-border rounded shadow hover:bg-accent"
                    title="Export Flowsheet"
                  >
                    <Download className="w-3.5 h-3.5" />
                  </button>
                  <button
                    onClick={() => {
                      const input = document.createElement('input');
                      input.type = 'file';
                      input.accept = '.json';
                      input.onchange = (e) => {
                        const file = (e.target as HTMLInputElement).files?.[0];
                        if (!file) return;
                        const reader = new FileReader();
                        reader.onload = (ev) => {
                          const ok = useCanopyStore.getState().importFlowsheet(ev.target?.result as string);
                          if (!ok) alert('Failed to import flowsheet — invalid file format.');
                        };
                        reader.readAsText(file);
                      };
                      input.click();
                    }}
                    className="p-1.5 text-xs bg-card/90 border border-border rounded shadow hover:bg-accent"
                    title="Import Flowsheet"
                  >
                    <Upload className="w-3.5 h-3.5" />
                  </button>
                  <button
                    onClick={() => {
                      const state = useCanopyStore.getState();
                      const { streams: ss, units: uu, compounds: cc, unitSystem: us2, fluidPackage: fp, solved: sv } = state;
                      const eco = (() => { try { return calculateEconomics(uu, ss); } catch { return null; } })();
                      const lines: string[] = [];
                      lines.push('CANOPY SIMULATION REPORT');
                      lines.push('='.repeat(60));
                      lines.push(`Generated: ${new Date().toISOString()}`);
                      lines.push(`Fluid Package: ${fp}`);
                      lines.push(`Solved: ${sv ? 'Yes' : 'No'}`);
                      lines.push(`Units: ${us2.temperature} / ${us2.pressure} / ${us2.flow}`);
                      lines.push('');
                      lines.push('COMPOUNDS');
                      lines.push('-'.repeat(40));
                      cc.forEach((c: CanopyCompound, i: number) => lines.push(`  ${i + 1}. ${c.displayName} (${c.name})`));
                      lines.push('');
                      lines.push('STREAM TABLE');
                      lines.push('-'.repeat(40));
                      lines.push(`${'Stream'.padEnd(16)} ${'T'.padStart(10)} ${'P'.padStart(10)} ${'Flow'.padStart(12)} ${'VF'.padStart(6)}`);
                      Object.values(ss).forEach((s: MaterialStream) => {
                        lines.push(`${s.name.padEnd(16)} ${fromSI_T(s.T_K, us2.temperature).toFixed(1).padStart(10)} ${fromSI_P(s.P_Pa, us2.pressure).toFixed(2).padStart(10)} ${fromSI_F(s.totalFlow_molps, us2.flow).toFixed(2).padStart(12)} ${(s.vaporFraction ?? 0).toFixed(3).padStart(6)}`);
                        if (s.moleFractions.length > 0 && cc.length > 0) {
                          const fracs = s.moleFractions.map((z: number, ci: number) => `${cc[ci]?.displayName ?? ci}: ${z.toFixed(4)}`).join(', ');
                          lines.push(`  Composition: ${fracs}`);
                        }
                      });
                      lines.push('');
                      lines.push('UNIT OPERATIONS');
                      lines.push('-'.repeat(40));
                      Object.entries(uu).forEach(([uid, u]) => {
                        lines.push(`${u.name} (${u.type})`);
                        if (u.params) Object.entries(u.params).forEach(([k, v]) => lines.push(`  ${k}: ${v}`));
                        if (u.duty_W !== undefined) lines.push(`  Duty: ${u.duty_W.toFixed(0)} W`);
                      });
                      if (eco) {
                        lines.push('');
                        lines.push('ECONOMICS SUMMARY');
                        lines.push('-'.repeat(40));
                        lines.push(`  Total Equipment Cost: $${eco.totalPurchaseCost.toFixed(0)}`);
                        lines.push(`  Total Bare Module Cost: $${eco.totalBareModuleCost.toFixed(0)}`);
                        lines.push(`  Grass Roots Cost: $${eco.totalGrassRootsCost.toFixed(0)}`);
                        lines.push(`  Annual Utility Cost: $${eco.totalAnnualUtilityCost.toFixed(0)}/yr`);
                      }
                      lines.push('');
                      lines.push('END OF REPORT');
                      const blob = new Blob([lines.join('\n')], { type: 'text/plain' });
                      const url = URL.createObjectURL(blob);
                      const a = document.createElement('a'); a.href = url; a.download = 'canopy_report.txt'; a.click();
                      URL.revokeObjectURL(url);
                    }}
                    className="p-1.5 text-xs bg-card/90 border border-border rounded shadow hover:bg-accent"
                    title="Generate Report"
                  >
                    <FileText className="w-3.5 h-3.5" />
                  </button>
                  {solved && (
                    <button
                      onClick={() => setBottomPanelHidden(hidden => !hidden)}
                      className="px-2 py-1.5 text-xs bg-card/90 border border-border rounded shadow hover:bg-accent"
                      title={bottomPanelHidden ? 'Open bottom results panel' : 'Hide bottom results panel'}
                    >
                      {bottomPanelHidden ? 'Results' : 'Hide Results'}
                    </button>
                  )}
                  <div className="w-px h-4 bg-border" />
                  <button
                    onClick={() => setSolverPaused(!solverPaused)}
                    className={`p-1.5 text-xs border rounded shadow ${
                      solverPaused
                        ? 'bg-yellow-500/90 border-yellow-600 text-white hover:bg-yellow-600'
                        : 'bg-card/90 border-border hover:bg-accent'
                    }`}
                    title={solverPaused ? 'Resume solver' : 'Pause solver (freeze results while editing)'}
                  >
                    {solverPaused ? <Play className="w-3.5 h-3.5" /> : <Pause className="w-3.5 h-3.5" />}
                  </button>
                </div>
                {solverPaused && (
                  <div className="text-xs text-yellow-500 bg-card/80 px-2 py-1 rounded border border-yellow-500/50">
                    Solver paused — changes won&apos;t auto-solve
                  </div>
                )}
                {selectedUnitId && (
                  <button
                    onClick={() => removeUnit(selectedUnitId)}
                    className="flex items-center gap-1.5 px-2 py-1 text-xs bg-red-500/90 hover:bg-red-600 text-white rounded shadow border border-red-600"
                    title="Delete selected unit (or press Delete)"
                  >
                    <Trash2 className="w-3 h-3" />
                    Delete {units[selectedUnitId]?.name}
                  </button>
                )}
                {selectedFeedId && !selectedUnitId && (
                  <button
                    onClick={() => removeStream(selectedFeedId)}
                    className="flex items-center gap-1.5 px-2 py-1 text-xs bg-red-500/90 hover:bg-red-600 text-white rounded shadow border border-red-600"
                    title="Delete selected feed stream (or press Delete)"
                  >
                    <Trash2 className="w-3 h-3" />
                    Delete {streams[selectedFeedId]?.name}
                  </button>
                )}
                {selectedEdgeStreamId && !selectedUnitId && !selectedFeedId && (
                  <button
                    onClick={() => disconnectStream(selectedEdgeStreamId)}
                    className="flex items-center gap-1.5 px-2 py-1 text-xs bg-red-500/90 hover:bg-red-600 text-white rounded shadow border border-red-600"
                    title="Delete selected connection (or press Delete)"
                  >
                    <Trash2 className="w-3 h-3" />
                    Delete connection {streams[selectedEdgeStreamId]?.name}
                  </button>
                )}
                <div className="text-xs text-muted-foreground bg-card/80 px-2 py-1 rounded border border-border">
                  Double-click to edit · Select + Delete to remove/disconnect
                </div>
              </div>
            </Panel>
          </ReactFlow>

          {/* Edge right-click context menu */}
          {edgeMenu && (
            <div
              className="fixed z-50 bg-card border border-border rounded-lg shadow-xl py-1 min-w-[160px]"
              style={{ left: edgeMenu.x, top: edgeMenu.y }}
              onClick={(e) => e.stopPropagation()}
            >
              <button
                className="w-full text-left px-3 py-1.5 text-xs hover:bg-accent transition-colors flex items-center gap-2"
                onClick={() => {
                  openTab({ type: 'stream', id: edgeMenu.streamId });
                  setEdgeMenu(null);
                }}
              >
                📊 View Stream Details
              </button>
              <button
                className="w-full text-left px-3 py-1.5 text-xs hover:bg-accent transition-colors flex items-center gap-2 text-red-500"
                onClick={() => {
                  disconnectStream(edgeMenu.streamId);
                  setEdgeMenu(null);
                }}
              >
                <Unlink className="w-3 h-3" /> Disconnect Stream
              </button>
            </div>
          )}
        </div>
      ) : activeTab.type === 'unit' ? (
        <UnitDetailPanel unitId={activeTab.id} />
      ) : (
        <StreamDetailPanel streamId={activeTab.id} />
      )}

      {/* Stream Results Table + Utility Summary */}
      {solved && activeTab === null && !bottomPanelHidden && (
        <div className="border-t border-border bg-card max-h-[32vh] overflow-auto">
          <BottomPanelTabs onClose={() => setBottomPanelHidden(true)} />
        </div>
      )}

      {/* Status Bar */}
      <div className="flex items-center justify-between px-3 py-1 border-t border-border bg-muted/30 text-xs">
        <div className="flex items-center gap-3">
          <span className={solved ? 'text-green-600' : solverErrors.length > 0 ? 'text-red-500' : 'text-muted-foreground'}>
            {solved ? '● Converged' : solverErrors.length > 0 ? '● Errors' : '○ Not solved'}
          </span>
          <span className="text-muted-foreground">
            {Object.values(units).filter(u => u.solved).length}/{Object.values(units).length} units solved
          </span>
          <span className="text-muted-foreground">
            {Object.values(streams).filter(s => s.solved).length}/{Object.values(streams).length} streams
          </span>
          {solverPaused && <span className="text-yellow-500 font-medium">⏸ Paused</span>}
        </div>
        <div className="flex items-center gap-2">
          {solved && (
            <button
              onClick={() => setBottomPanelHidden(hidden => !hidden)}
              className="rounded border border-border px-2 py-0.5 text-xs text-foreground hover:bg-accent"
              title={bottomPanelHidden ? 'Open bottom results panel' : 'Hide bottom results panel'}
            >
              {bottomPanelHidden ? 'Open Results' : 'Hide Results'}
            </button>
          )}
          {solverErrors.length > 0 && (
            <span className="text-red-500 truncate max-w-[300px]" title={solverErrors.join('\n')}>
              {solverErrors[0]}
            </span>
          )}
          {dbWarnings.filter(w => w.startsWith('Validation:')).length > 0 && (
            <span className="text-amber-500" title={dbWarnings.filter(w => w.startsWith('Validation:')).join('\n')}>
              ⚠ {dbWarnings.filter(w => w.startsWith('Validation:')).length} warning(s)
            </span>
          )}
          <span className="text-muted-foreground">Use Results or double-click any stream to inspect it.</span>
          <span className="text-muted-foreground">{compounds.length} compounds · {unitSystem.temperature}/{unitSystem.pressure}</span>
        </div>
      </div>
    </div>
  );
}

// ─── Stream Detail / Editor Panel ───────────────────────────────

function StreamDetailPanel({ streamId }: { streamId: string }) {
  const { streams, compounds, updateStream, units, unitSystem: us, renameStream } = useCanopyStore();
  const stream = streams[streamId];
  const [editingName, setEditingName] = useState(false);
  const [nameInput, setNameInput] = useState('');

  if (!stream) return <div className="p-4 text-muted-foreground">Stream not found.</div>;

  const isUpstreamDetermined = !!stream.sourceUnitId;
  const sourceUnit = stream.sourceUnitId ? units[stream.sourceUnitId] : null;

  const handleTempChange = (val: number) => {
    updateStream(streamId, { T_K: toSI_T(val, us.temperature), solved: false });
  };
  const handlePressChange = (val: number) => {
    updateStream(streamId, { P_Pa: toSI_P(val, us.pressure), solved: false });
  };
  const handleFlowChange = (val: number) => {
    updateStream(streamId, { totalFlow_molps: toSI_F(val, us.flow), solved: false });
  };
  const handleMoleFracChange = (idx: number, val: number) => {
    const newFracs = [...stream.moleFractions];
    newFracs[idx] = val;
    updateStream(streamId, { moleFractions: newFracs, solved: false });
  };

  const normalizeFracs = () => {
    const sum = stream.moleFractions.reduce((a, b) => a + b, 0);
    if (sum > 0) {
      updateStream(streamId, {
        moleFractions: stream.moleFractions.map(x => x / sum),
        solved: false,
      });
    }
  };

  const fracSum = stream.moleFractions.reduce((a, b) => a + b, 0);

  const startRename = () => {
    setNameInput(stream.name);
    setEditingName(true);
  };

  const commitRename = () => {
    if (nameInput.trim()) renameStream(streamId, nameInput.trim());
    setEditingName(false);
  };

  return (
    <div className="flex-1 overflow-auto p-6 space-y-6">
      <div className="flex items-start justify-between gap-3">
        <div className="space-y-2">
          <div className="flex flex-wrap items-center gap-2">
            <span className="text-xs px-2 py-0.5 rounded bg-muted text-muted-foreground">
              {stream.id}
            </span>
            {isUpstreamDetermined && (
              <span className="text-xs px-2 py-0.5 rounded bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300">
                Set by {sourceUnit?.name ?? 'upstream unit'}
              </span>
            )}
            {stream.solved && (
              <span className="text-xs px-2 py-0.5 rounded bg-green-100 text-green-800 dark:bg-green-900/40 dark:text-green-300">
                Solved
              </span>
            )}
          </div>
          {editingName ? (
            <input
              value={nameInput}
              onChange={e => setNameInput(e.target.value)}
              onBlur={commitRename}
              onKeyDown={e => { if (e.key === 'Enter') commitRename(); if (e.key === 'Escape') setEditingName(false); }}
              className="text-lg font-bold border-b-2 border-blue-500 outline-none bg-transparent px-0 py-0.5"
              autoFocus
            />
          ) : (
            <button
              type="button"
              className="text-left text-lg font-bold text-foreground hover:text-blue-500 transition-colors"
              onClick={startRename}
              title="Click to rename"
            >
              {stream.name}
            </button>
          )}
        </div>
      </div>

      {/* State Variables */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
        <FieldInput
          label={<>Temperature ({tempLabel(us.temperature)})</>}
          value={fromSI_T(stream.T_K, us.temperature)}
          onChange={handleTempChange}
          locked={isUpstreamDetermined}
          decimals={2}
        />
        <FieldInput
          label={<>Pressure ({pressLabel(us.pressure)})</>}
          value={fromSI_P(stream.P_Pa, us.pressure)}
          onChange={handlePressChange}
          locked={isUpstreamDetermined}
          decimals={3}
        />
        <FieldInput
          label={<>Molar Flow ({flowLabel(us.flow)})</>}
          value={fromSI_F(stream.totalFlow_molps, us.flow)}
          onChange={handleFlowChange}
          locked={isUpstreamDetermined}
          decimals={2}
        />
      </div>

      {/* Composition */}
      <div className="space-y-2">
        <div className="flex items-center justify-between">
          <h3 className="text-sm font-semibold">
            Composition (mole fractions)
          </h3>
          {!isUpstreamDetermined && (
            <div className="flex items-center gap-2">
              <span className={`text-xs ${Math.abs(fracSum - 1) < 0.001 ? 'text-green-600' : 'text-amber-600'}`}>
                Σz = {fracSum.toFixed(4)}
              </span>
              <Button size="sm" variant="outline" onClick={normalizeFracs} className="text-xs h-6 px-2">
                Normalize
              </Button>
            </div>
          )}
        </div>
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-2">
          {compounds.map((c, i) => (
            <FieldInput
              key={c.name}
              label={<>z<sub>{c.displayName}</sub></>}
              value={stream.moleFractions[i] ?? 0}
              onChange={(v) => handleMoleFracChange(i, v)}
              locked={isUpstreamDetermined}
              decimals={4}
              min={0}
              max={1}
              step={0.01}
            />
          ))}
        </div>

        {/* Composition Pie Chart */}
        {compounds.length > 0 && (() => {
          const PIE_COLORS = ['#3b82f6', '#ef4444', '#22c55e', '#f59e0b', '#8b5cf6', '#ec4899', '#06b6d4', '#84cc16', '#f97316', '#6366f1'];
          const fracs = stream.moleFractions;
          const total = fracs.reduce((a, b) => a + b, 0) || 1;
          let cumAngle = 0;
          const R = 50, CX = 60, CY = 60;
          return (
            <div className="flex items-center gap-4 mt-2">
              <svg width="120" height="120" viewBox="0 0 120 120">
                {fracs.map((f, i) => {
                  const frac = f / total;
                  if (frac < 0.001) return null;
                  const startAngle = cumAngle;
                  cumAngle += frac * 2 * Math.PI;
                  const endAngle = cumAngle;
                  const largeArc = frac > 0.5 ? 1 : 0;
                  const x1 = CX + R * Math.cos(startAngle - Math.PI / 2);
                  const y1 = CY + R * Math.sin(startAngle - Math.PI / 2);
                  const x2 = CX + R * Math.cos(endAngle - Math.PI / 2);
                  const y2 = CY + R * Math.sin(endAngle - Math.PI / 2);
                  if (frac > 0.999) {
                    return <circle key={i} cx={CX} cy={CY} r={R} fill={PIE_COLORS[i % PIE_COLORS.length]} />;
                  }
                  return (
                    <path key={i}
                      d={`M ${CX} ${CY} L ${x1} ${y1} A ${R} ${R} 0 ${largeArc} 1 ${x2} ${y2} Z`}
                      fill={PIE_COLORS[i % PIE_COLORS.length]}
                      stroke="var(--card)" strokeWidth="1"
                    />
                  );
                })}
              </svg>
              <div className="flex flex-col gap-0.5 text-xs">
                {compounds.map((c, i) => (
                  <div key={c.name} className="flex items-center gap-1.5">
                    <span className="w-2.5 h-2.5 rounded-sm inline-block" style={{ backgroundColor: PIE_COLORS[i % PIE_COLORS.length] }} />
                    <span>{c.displayName}: {((fracs[i] / total) * 100).toFixed(1)}%</span>
                  </div>
                ))}
              </div>
            </div>
          );
        })()}
      </div>

      {/* Results (if solved) */}
      {stream.solved && (
        <div className="space-y-3 border-t border-border pt-4">
          <div className="grid grid-cols-2 md:grid-cols-4 gap-3 text-sm">
            <ResultField label="Phase" value={stream.phase} />
            <ResultField label="Vapor Fraction" value={fmt(stream.vaporFraction, 4)} />
            <ResultField label="Enthalpy (J/mol)" value={fmt(stream.H_Jpmol, 1)} />
            <ResultField label={<>T ({tempLabel(us.temperature)})</>} value={fmt(fromSI_T(stream.T_K, us.temperature), 2)} />
          </div>

          {stream.vaporFraction > 0 && stream.vaporFraction < 1 && (
            <div className="space-y-2">
              <h4 className="text-xs font-semibold text-muted-foreground">Phase Compositions</h4>
              <table className="text-xs w-full">
                <thead>
                  <tr className="border-b border-border">
                    <th className="text-left p-1">Compound</th>
                    <th className="text-right p-1">z<sub>i</sub> (feed)</th>
                    <th className="text-right p-1">x<sub>i</sub> (liquid)</th>
                    <th className="text-right p-1">y<sub>i</sub> (vapor)</th>
                  </tr>
                </thead>
                <tbody>
                  {compounds.map((c, i) => (
                    <tr key={c.name} className="border-b border-border/50">
                      <td className="p-1">{c.displayName}</td>
                      <td className="p-1 text-right">{fmt(stream.moleFractions[i], 4)}</td>
                      <td className="p-1 text-right">{fmt(stream.x_liquid[i], 4)}</td>
                      <td className="p-1 text-right">{fmt(stream.y_vapor[i], 4)}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

// ─── Parameter Label Formatting ─────────────────────────────────

type UsType = ReturnType<typeof useCanopyStore.getState>['unitSystem'];

function formatParamLabel(key: string, us: UsType): React.ReactNode {
  const map: Record<string, (us: UsType) => React.ReactNode> = {
    targetT_K: (us) => <>Target T ({tempLabel(us.temperature)})</>,
    outletT_K: (us) => <>Outlet T ({tempLabel(us.temperature)})</>,
    T_K: (us) => <>Temperature ({tempLabel(us.temperature)})</>,
    outletP_Pa: (us) => <>Outlet P ({pressLabel(us.pressure)})</>,
    P_Pa: (us) => <>Pressure ({pressLabel(us.pressure)})</>,
    condenserP_Pa: (us) => <>Condenser P ({pressLabel(us.pressure)})</>,
    reboilerP_Pa: (us) => <>Reboiler P ({pressLabel(us.pressure)})</>,
    pressure_Pa: (us) => <>Pressure ({pressLabel(us.pressure)})</>,
    gasT_in_K: (us) => <>Gas T<sub>in</sub> ({tempLabel(us.temperature)})</>,
    T_ref_K: (us) => <>T<sub>ref</sub> ({tempLabel(us.temperature)})</>,
    splitRatio: () => 'Split Ratio',
    efficiency: () => 'Efficiency',
    model: () => 'Model',
    spec: () => 'Specification',
    specValue: () => 'Spec Value',
    flowArrangement: () => 'Flow Arrangement',
    hotOutletP_Pa: (us) => <>Hot Outlet P ({pressLabel(us.pressure)})</>,
    coldOutletP_Pa: (us) => <>Cold Outlet P ({pressLabel(us.pressure)})</>,
    driverEfficiency: () => 'Driver Efficiency',
    length_m: () => 'Length (m)',
    diameter_m: () => 'Diameter (m)',
    roughness_m: () => 'Roughness (m)',
    elevation_m: () => 'Elevation (m)',
    K_fittings: () => <>K<sub>fittings</sub></>,
    lightKeyIndex: () => 'Light Key Index',
    heavyKeyIndex: () => 'Heavy Key Index',
    lightKeyRecovery: () => 'Light Key Recovery',
    heavyKeyRecovery: () => 'Heavy Key Recovery',
    refluxRatioMultiplier: () => <>R / R<sub>min</sub></>,
    condenserType: () => 'Condenser Type',
    rigorousMode: () => 'Rigorous Mode',
    nStages: () => 'Number of Stages',
    feedStage: () => 'Feed Stage',
    refluxRatio: () => 'Reflux Ratio',
    distillateRate_molps: () => 'Distillate Rate (mol/s)',
    volume_m3: () => 'Volume (m³)',
    batchTime_s: () => 'Batch Time (s)',
    residenceTime_s: () => 'Residence Time (s)',
    soluteIndex: () => 'Solute Index',
    C_sat_molpm3: () => <>C<sub>sat</sub> (mol/m³)</>,
    k_g: () => <>k<sub>g</sub> (growth)</>,
    g_exponent: () => 'Growth Exponent',
    k_b: () => <>k<sub>b</sub> (nucleation)</>,
    b_exponent: () => 'Nucleation Exponent',
    workIndex_kWhpton: () => 'Bond Work Index (kWh/ton)',
    feedSize_um: () => <>d<sub>80,feed</sub> (μm)</>,
    productSize_um: () => <>d<sub>80,product</sub> (μm)</>,
    solidFlow_kgps: () => 'Solid Flow (kg/s)',
    X_in_kgpkg: () => <>X<sub>in</sub> (kg/kg)</>,
    X_out_kgpkg: () => <>X<sub>out</sub> (kg/kg)</>,
    X_c_kgpkg: () => <>X<sub>c</sub> (critical, kg/kg)</>,
    X_eq_kgpkg: () => <>X<sub>eq</sub> (equilibrium, kg/kg)</>,
    dryFlow_kgps: () => 'Dry Flow (kg/s)',
    gasFlow_kgps: () => 'Gas Flow (kg/s)',
    Lv_Jpkg: () => <>L<sub>v</sub> (J/kg)</>,
    deltaT_min_K: () => 'ΔT min (K)',
    UA_WpK: () => 'UA (W/K)',
    isothermal: () => 'Isothermal',
    nSteps: () => 'Integration Steps',
    duty_W: () => 'Duty (W)',
    Ea_Jpmol: () => <>E<sub>a</sub> (J/mol)</>,
    k_v: () => <>k<sub>v</sub> (shape factor)</>,
    rho_crystal_kgpm3: () => <>ρ<sub>crystal</sub> (kg/m³)</>,
    elementMatrix: () => 'Element Matrix',
    measuredObjectType: () => 'Measured Object Type',
    measuredObjectId: () => 'Measured Object',
    propertyCode: () => 'Property Code',
    componentIndex: () => 'Component Index',
    signalValue: () => 'Signal Value',
    propertyLabel: () => 'Property Label',
    adjustComponentIndex: () => 'Adjustment Species',
    targetCharge: () => 'Target Charge',
    residualCharge: () => 'Residual Charge',
    adjustedComponentFlow_molps: () => 'Adjusted Species Flow (mol/s)',
    dirtyWaterMode: () => 'Dirty Water Mode',
    freeWaterSplitFraction: () => 'Free Water Split Fraction',
    waterCarryoverFraction: () => 'Water Carryover Fraction',
    hydrocarbonToAqueousFraction: () => 'HC to Aqueous Fraction',
    freeWaterSeparated_molps: () => 'Free Water Separated (mol/s)',
    waterSaturationFraction: () => 'Water Saturation Fraction',
    aqueousPhaseIndex: () => 'Aqueous Phase Outlet',
    waterSaturatedFeed: () => 'Water-Saturated Feed',
    waterAddedForSaturation_molps: () => 'Saturation Water Added (mol/s)',
    excessFreeWater_molps: () => 'Excess Free Water (mol/s)',
    apparentSpeciesIndex: () => 'Apparent Species',
    cationIndex: () => 'Cation Species',
    anionIndex: () => 'Anion Species',
    cationStoich: () => 'Cation Stoich',
    anionStoich: () => 'Anion Stoich',
    log10K: () => 'log10 K',
    maxDissociationFraction: () => 'Max Dissociation Fraction',
    reactionSetJson: () => 'Reaction Set JSON',
    dissociationFraction: () => 'Dissociation Fraction',
    ionicStrength: () => 'Ionic Strength',
    equilibriumConstant: () => 'Equilibrium Constant',
    activeReactionCount: () => 'Active Reaction Count',
    equilibriumPasses: () => 'Equilibrium Passes',
    internalsMode: () => 'Internals Mode',
    pressureRelaxation: () => 'Pressure Relaxation',
    pressureCouplingResidual_Pa: () => 'Pressure Coupling Residual (Pa)',
  };
  const fn = map[key];
  if (fn) return fn(us);
  // Fallback: convert camelCase/snake_case to readable form
  return key
    .replace(/_/g, ' ')
    .replace(/([A-Z])/g, ' $1')
    .replace(/^\s+/, '')
    .split(' ')
    .filter(Boolean)
    .map(w => w.charAt(0).toUpperCase() + w.slice(1).toLowerCase())
    .join(' ');
}

function isTemperatureParam(key: string): boolean {
  return key.endsWith('T_K') || key === 'T_K' || key === 'gasT_in_K' || key === 'T_ref_K';
}

function isPressureParam(key: string): boolean {
  return key.endsWith('P_Pa') || key === 'P_Pa' || key === 'pressure_Pa' || key === 'hotOutletP_Pa' || key === 'coldOutletP_Pa';
}

const COLUMN_SPEC_OPTIONS: { value: RigorousColumnSpecMode; label: string }[] = [
  { value: 'None', label: 'None' },
  { value: 'DistillateRate', label: 'Distillate rate' },
  { value: 'DistillateFraction', label: 'Distillate fraction' },
  { value: 'BottomsRate', label: 'Bottoms rate' },
  { value: 'RefluxRatio', label: 'Reflux ratio' },
  { value: 'BoilupRatio', label: 'Boilup ratio' },
  { value: 'DistillateComposition', label: 'Distillate composition' },
];

const COLUMN_SPEC_PARAM_KEYS = new Set([
  'spec1Mode', 'spec1Value', 'spec1ComponentIndex',
  'spec2Mode', 'spec2Value', 'spec2ComponentIndex',
]);

const ANALYZER_PROPERTY_OPTIONS = [
  { value: 'TEMP', label: 'Temperature' },
  { value: 'PRES', label: 'Pressure' },
  { value: 'ENTH', label: 'Molar Enthalpy' },
  { value: 'MOLE_FLOW', label: 'Molar Flow' },
  { value: 'MASS_FLOW', label: 'Mass Flow' },
  { value: 'VFRA', label: 'Vapor Fraction' },
  { value: 'LFRA', label: 'Liquid Fraction' },
  { value: 'MOLE_FRAC', label: 'Component Mole Fraction' },
  { value: 'MASS_FRAC', label: 'Component Mass Fraction' },
  { value: 'TBPT', label: 'Average TBP' },
  { value: 'VABP', label: 'Volume Avg Boiling Point' },
  { value: 'CETA', label: 'Cetane Number' },
  { value: 'AIT', label: 'Autoignition Temp' },
  { value: 'D86T', label: 'ASTM D86 Equivalent Temp' },
  { value: 'D1160T', label: 'ASTM D1160 Equivalent Temp' },
  { value: 'D2887T', label: 'ASTM D2887 Equivalent Temp' },
  { value: 'BUBPT', label: 'Bubble Point Temp' },
  { value: 'DEWPT', label: 'Dew Point Temp' },
  { value: 'WATBUB', label: 'Water-Sat Bubble Point' },
  { value: 'WATDEW', label: 'Water-Sat Dew Point' },
  { value: 'SG', label: 'Specific Gravity' },
  { value: 'API', label: 'API Gravity' },
  { value: 'REFINDEX', label: 'Refractive Index' },
  { value: 'SMX', label: 'Mixture Entropy' },
  { value: 'WAT', label: 'Water Content' },
  { value: 'WATSAT', label: 'Water Saturation' },
  { value: 'GRS', label: 'Gross Heating Value' },
  { value: 'NET', label: 'Net Heating Value' },
  { value: 'FLASHPOINT', label: 'Flash Point' },
  { value: 'ANILPT', label: 'Aniline Point' },
  { value: 'REIDVP', label: 'Reid Vapor Pressure' },
] as const;

const SPECIAL_UNIT_PARAM_KEYS = new Set([
  'measuredObjectType',
  'measuredObjectId',
  'propertyCode',
  'componentIndex',
  'signalValue',
  'propertyLabel',
  'adjustComponentIndex',
  'targetCharge',
  'residualCharge',
  'adjustedComponentFlow_molps',
  'freeWaterSeparated_molps',
  'waterSaturationFraction',
  'aqueousPhaseIndex',
  'waterAddedForSaturation_molps',
  'excessFreeWater_molps',
  'apparentSpeciesIndex',
  'cationIndex',
  'anionIndex',
  'cationStoich',
  'anionStoich',
  'log10K',
  'maxDissociationFraction',
  'reactionSetJson',
  'dissociationFraction',
  'ionicStrength',
  'equilibriumConstant',
  'activeReactionCount',
  'equilibriumPasses',
  'pressureCouplingResidual_Pa',
]);

// ─── Unit Detail Panel ──────────────────────────────────────

function UnitDetailPanel({ unitId }: { unitId: string }) {
  const { units, streams, compounds, updateUnit, openTab, unitSystem: us, renameUnit } = useCanopyStore();
  const unit = units[unitId];
  const [editingName, setEditingName] = useState(false);
  const [nameInput, setNameInput] = useState('');

  if (!unit) return <div className="p-4 text-muted-foreground">Unit not found.</div>;

  const displayValue = (key: string, val: number): number => {
    if (isTemperatureParam(key)) return fromSI_T(val, us.temperature);
    if (isPressureParam(key)) return fromSI_P(val, us.pressure);
    return val;
  };

  const fromDisplay = (key: string, display: number): number => {
    if (isTemperatureParam(key)) return toSI_T(display, us.temperature);
    if (isPressureParam(key)) return toSI_P(display, us.pressure);
    return display;
  };

  const isRigorousColumn = unit.type === 'RadFrac' || unit.type === 'RateFrac' || (unit.type === 'Column' && Boolean(unit.params.rigorousMode));
  const isAnalyzer = unit.type === 'Analyzer';
  const isChargeBalance = unit.type === 'ChargeBalance';
  const isDecanter = unit.type === 'Decanter';
  const isElectrolyteEquilibrium = unit.type === 'ElectrolyteEquilibrium';
  const isWaterHandledUnit = unit.type === 'Flash' || unit.type === 'Heater' || unit.type === 'Valve';
  const feedPort = unit.ports.find(port => port.type === 'inlet' && port.streamId);
  const feedStream = feedPort?.streamId ? streams[feedPort.streamId] : null;
  const specDiagnostics = isRigorousColumn
    ? getRigorousColumnSpecDiagnostics({
        nStages: (unit.params.nStages as number) ?? 20,
        feedStage: (unit.params.feedStage as number) ?? 10,
        refluxRatio: (unit.params.refluxRatio as number) ?? 2,
        distillateRate_molps: (unit.params.distillateRate_molps as number) ?? feedStream?.totalFlow_molps ?? 0,
        condenserP_Pa: (unit.params.condenserP_Pa as number) ?? 101325,
        reboilerP_Pa: (unit.params.reboilerP_Pa as number) ?? 101325,
        condenserType: (unit.params.condenserType as 'Total' | 'Partial') ?? 'Total',
        spec1Mode: unit.params.spec1Mode as RigorousColumnSpecMode | undefined,
        spec1Value: unit.params.spec1Value as number | undefined,
        spec1ComponentIndex: unit.params.spec1ComponentIndex as number | undefined,
        spec2Mode: unit.params.spec2Mode as RigorousColumnSpecMode | undefined,
        spec2Value: unit.params.spec2Value as number | undefined,
        spec2ComponentIndex: unit.params.spec2ComponentIndex as number | undefined,
      })
    : null;

  // Filter out non-editable params (arrays, objects)
  const editableParams = Object.entries(unit.params).filter(
    ([key, value]) =>
      !COLUMN_SPEC_PARAM_KEYS.has(key)
      && !SPECIAL_UNIT_PARAM_KEYS.has(key)
      && (typeof value === 'number' || typeof value === 'string' || typeof value === 'boolean')
  );

  const startRename = () => {
    setNameInput(unit.name);
    setEditingName(true);
  };

  const commitRename = () => {
    if (nameInput.trim()) renameUnit(unitId, nameInput.trim());
    setEditingName(false);
  };

  return (
    <div className="flex-1 overflow-auto p-6 space-y-4">
      <div className="flex items-start justify-between gap-3">
        <div className="space-y-2">
          <div className="flex flex-wrap items-center gap-2">
            <span className="text-xs px-2 py-0.5 rounded bg-muted text-muted-foreground">{unit.type}</span>
            {unit.solved && (
              <span className="text-xs px-2 py-0.5 rounded bg-green-100 text-green-800 dark:bg-green-900/40 dark:text-green-300">
                Solved
              </span>
            )}
          </div>
          {editingName ? (
            <input
              value={nameInput}
              onChange={e => setNameInput(e.target.value)}
              onBlur={commitRename}
              onKeyDown={e => { if (e.key === 'Enter') commitRename(); if (e.key === 'Escape') setEditingName(false); }}
              className="text-lg font-bold border-b-2 border-blue-500 outline-none bg-transparent px-0 py-0.5"
              autoFocus
            />
          ) : (
            <button
              type="button"
              className="text-left text-lg font-bold text-foreground hover:text-blue-500 transition-colors"
              onClick={startRename}
              title="Click to rename"
            >
              {unit.name}
            </button>
          )}
        </div>
      </div>

      {/* Parameters */}
      {isAnalyzer && (
        <div className="space-y-3 rounded border border-border/60 bg-muted/20 p-3">
          <div className="flex items-center justify-between gap-3">
            <div>
              <h3 className="text-sm font-semibold">Analyzer</h3>
              <p className="text-xs text-muted-foreground">
                Aspen-style property monitor for streams and signal-producing units.
              </p>
            </div>
            <div className="rounded bg-muted px-2 py-1 text-xs text-muted-foreground">
              {typeof unit.params.signalValue === 'number' ? Number(unit.params.signalValue).toFixed(4) : 'No signal'}
            </div>
          </div>
          <div className="grid grid-cols-1 gap-3 md:grid-cols-2">
            <div className="space-y-1">
              <label className="text-xs font-medium text-muted-foreground">Source Type</label>
              <select
                value={String(unit.params.measuredObjectType ?? 'stream')}
                onChange={(e) => updateUnit(unitId, { params: { ...unit.params, measuredObjectType: e.target.value } })}
                className="w-full px-2 py-1.5 text-sm border rounded bg-background border-input focus:outline-none focus:ring-2 focus:ring-ring"
              >
                <option value="stream">Stream</option>
                <option value="unit">Unit</option>
              </select>
            </div>
            <div className="space-y-1">
              <label className="text-xs font-medium text-muted-foreground">Measured Object</label>
              <select
                value={String(unit.params.measuredObjectId ?? '')}
                onChange={(e) => updateUnit(unitId, { params: { ...unit.params, measuredObjectId: e.target.value } })}
                className="w-full px-2 py-1.5 text-sm border rounded bg-background border-input focus:outline-none focus:ring-2 focus:ring-ring"
              >
                <option value="">Connected stream / none</option>
                {String(unit.params.measuredObjectType ?? 'stream') === 'stream'
                  ? Object.values(streams).map(stream => <option key={stream.id} value={stream.id}>{stream.name}</option>)
                  : Object.values(units).filter(candidate => candidate.id !== unitId).map(candidate => (
                      <option key={candidate.id} value={candidate.id}>{candidate.name}</option>
                    ))}
              </select>
            </div>
            <div className="space-y-1">
              <label className="text-xs font-medium text-muted-foreground">Property</label>
              <select
                value={String(unit.params.propertyCode ?? 'TEMP')}
                onChange={(e) => updateUnit(unitId, { params: { ...unit.params, propertyCode: e.target.value } })}
                className="w-full px-2 py-1.5 text-sm border rounded bg-background border-input focus:outline-none focus:ring-2 focus:ring-ring"
              >
                {ANALYZER_PROPERTY_OPTIONS.map(option => (
                  <option key={option.value} value={option.value}>{option.label}</option>
                ))}
              </select>
            </div>
            {(unit.params.propertyCode === 'MOLE_FRAC' || unit.params.propertyCode === 'MASS_FRAC') && (
              <div className="space-y-1">
                <label className="text-xs font-medium text-muted-foreground">Component</label>
                <select
                  value={Number(unit.params.componentIndex ?? 0)}
                  onChange={(e) => updateUnit(unitId, { params: { ...unit.params, componentIndex: parseInt(e.target.value, 10) } })}
                  className="w-full px-2 py-1.5 text-sm border rounded bg-background border-input focus:outline-none focus:ring-2 focus:ring-ring"
                >
                  {compounds.map((compound, index) => (
                    <option key={compound.name} value={index}>{compound.displayName}</option>
                  ))}
                </select>
              </div>
            )}
          </div>
          {typeof unit.params.propertyLabel === 'string' && (
            <div className="text-xs text-muted-foreground">Output: {String(unit.params.propertyLabel)}</div>
          )}
        </div>
      )}

      {isChargeBalance && (
        <div className="space-y-3 rounded border border-border/60 bg-muted/20 p-3">
          <div>
            <h3 className="text-sm font-semibold">Charge Balance</h3>
            <p className="text-xs text-muted-foreground">
              Adjust one ionic species to satisfy the requested net charge on the feed stream.
            </p>
          </div>
          <div className="grid grid-cols-1 gap-3 md:grid-cols-2">
            <div className="space-y-1">
              <label className="text-xs font-medium text-muted-foreground">Adjustment Species</label>
              <select
                value={Number(unit.params.adjustComponentIndex ?? 0)}
                onChange={(e) => updateUnit(unitId, { params: { ...unit.params, adjustComponentIndex: parseInt(e.target.value, 10) } })}
                className="w-full px-2 py-1.5 text-sm border rounded bg-background border-input focus:outline-none focus:ring-2 focus:ring-ring"
              >
                {compounds.map((compound, index) => (
                  <option key={compound.name} value={index}>
                    {compound.displayName}{compound.chargeNumber ? ` (${compound.chargeNumber > 0 ? '+' : ''}${compound.chargeNumber})` : ''}
                  </option>
                ))}
              </select>
            </div>
            <FieldInput
              label="Target Charge"
              value={Number(unit.params.targetCharge ?? 0)}
              onChange={(v) => updateUnit(unitId, { params: { ...unit.params, targetCharge: v } })}
              locked={false}
              decimals={4}
            />
          </div>
          <div className="grid grid-cols-1 gap-3 text-xs text-muted-foreground md:grid-cols-2">
            <div>Residual charge: {Number(unit.params.residualCharge ?? 0).toExponential(3)}</div>
            <div>Adjusted species flow: {Number(unit.params.adjustedComponentFlow_molps ?? 0).toFixed(4)} mol/s</div>
          </div>
        </div>
      )}

      {isDecanter && (
        <div className="space-y-3 rounded border border-border/60 bg-muted/20 p-3">
          <div>
            <h3 className="text-sm font-semibold">Decanter Results</h3>
            <p className="text-xs text-muted-foreground">
              Dirty-water mode applies APV140-style water and hydrocarbon solubility carryover estimates before the generic VLLE fallback.
            </p>
          </div>
          <div className="grid grid-cols-1 gap-3 text-xs text-muted-foreground md:grid-cols-3">
            <div>Free water separated: {Number(unit.params.freeWaterSeparated_molps ?? 0).toFixed(4)} mol/s</div>
            <div>Aqueous outlet: {unit.params.aqueousPhaseIndex == null ? 'Auto / none' : `Liquid ${String(unit.params.aqueousPhaseIndex)}`}</div>
            <div>Water saturation: {(Number(unit.params.waterSaturationFraction ?? 0) * 100).toFixed(1)}%</div>
          </div>
        </div>
      )}

      {isWaterHandledUnit && (
        <div className="space-y-3 rounded border border-border/60 bg-muted/20 p-3">
          <div>
            <h3 className="text-sm font-semibold">Water Handling</h3>
            <p className="text-xs text-muted-foreground">
              Enable water-saturated feed when the unit should assume contact with free water before equilibrium.
            </p>
          </div>
          <div className="grid grid-cols-1 gap-3 text-xs text-muted-foreground md:grid-cols-3">
            <div>Water added for saturation: {Number(unit.params.waterAddedForSaturation_molps ?? 0).toFixed(4)} mol/s</div>
            <div>Excess free water: {Number(unit.params.excessFreeWater_molps ?? 0).toFixed(4)} mol/s</div>
            <div>Water-saturated feed: {Boolean(unit.params.waterSaturatedFeed) ? 'Enabled' : 'Disabled'}</div>
          </div>
        </div>
      )}

      {isElectrolyteEquilibrium && (
        <div className="space-y-3 rounded border border-border/60 bg-muted/20 p-3">
          <div>
            <h3 className="text-sm font-semibold">Electrolyte Equilibrium</h3>
            <p className="text-xs text-muted-foreground">
              Apparent-species dissociation to ionic products with activity-based equilibrium closure. Add multiple reactions with the JSON set below when one apparent electrolyte is not enough.
            </p>
          </div>
          <div className="grid grid-cols-1 gap-3 md:grid-cols-3">
            <div className="space-y-1">
              <label className="text-xs font-medium text-muted-foreground">Apparent Species</label>
              <select
                value={Number(unit.params.apparentSpeciesIndex ?? 0)}
                onChange={(e) => updateUnit(unitId, { params: { ...unit.params, apparentSpeciesIndex: parseInt(e.target.value, 10) } })}
                className="w-full px-2 py-1.5 text-sm border rounded bg-background border-input focus:outline-none focus:ring-2 focus:ring-ring"
              >
                {compounds.map((compound, index) => (
                  <option key={compound.name} value={index}>{compound.displayName}</option>
                ))}
              </select>
            </div>
            <div className="space-y-1">
              <label className="text-xs font-medium text-muted-foreground">Cation</label>
              <select
                value={Number(unit.params.cationIndex ?? 0)}
                onChange={(e) => updateUnit(unitId, { params: { ...unit.params, cationIndex: parseInt(e.target.value, 10) } })}
                className="w-full px-2 py-1.5 text-sm border rounded bg-background border-input focus:outline-none focus:ring-2 focus:ring-ring"
              >
                {compounds.map((compound, index) => (
                  <option key={compound.name} value={index}>{compound.displayName}</option>
                ))}
              </select>
            </div>
            <div className="space-y-1">
              <label className="text-xs font-medium text-muted-foreground">Anion</label>
              <select
                value={Number(unit.params.anionIndex ?? 0)}
                onChange={(e) => updateUnit(unitId, { params: { ...unit.params, anionIndex: parseInt(e.target.value, 10) } })}
                className="w-full px-2 py-1.5 text-sm border rounded bg-background border-input focus:outline-none focus:ring-2 focus:ring-ring"
              >
                {compounds.map((compound, index) => (
                  <option key={compound.name} value={index}>{compound.displayName}</option>
                ))}
              </select>
            </div>
            <FieldInput
              label="Cation Stoich"
              value={Number(unit.params.cationStoich ?? 1)}
              onChange={(v) => updateUnit(unitId, { params: { ...unit.params, cationStoich: Math.max(1, v) } })}
              locked={false}
              decimals={0}
              min={1}
            />
            <FieldInput
              label="Anion Stoich"
              value={Number(unit.params.anionStoich ?? 1)}
              onChange={(v) => updateUnit(unitId, { params: { ...unit.params, anionStoich: Math.max(1, v) } })}
              locked={false}
              decimals={0}
              min={1}
            />
            <FieldInput
              label="log10 K"
              value={Number(unit.params.log10K ?? 2)}
              onChange={(v) => updateUnit(unitId, { params: { ...unit.params, log10K: v } })}
              locked={false}
              decimals={3}
            />
            <FieldInput
              label="Max Dissociation"
              value={Number(unit.params.maxDissociationFraction ?? 0.999)}
              onChange={(v) => updateUnit(unitId, { params: { ...unit.params, maxDissociationFraction: Math.max(0, Math.min(0.999999, v)) } })}
              locked={false}
              decimals={4}
              min={0}
              max={1}
              step={0.01}
            />
          </div>
          <div className="space-y-1">
            <label className="text-xs font-medium text-muted-foreground">Reaction Set JSON</label>
            <textarea
              value={String(unit.params.reactionSetJson ?? '[]')}
              onChange={(e) => updateUnit(unitId, { params: { ...unit.params, reactionSetJson: e.target.value } })}
              className="min-h-28 w-full rounded border border-input bg-background px-2 py-1.5 font-mono text-xs focus:outline-none focus:ring-2 focus:ring-ring"
              spellCheck={false}
            />
            <p className="text-[11px] text-muted-foreground">
              Legacy example: [{`{"name":"NaCl","apparentSpeciesIndex":1,"cationIndex":2,"anionIndex":3,"log10K":4,"active":true}`}]
              <br />
              Generic stoichiometric example: [{`{"name":"A2 dimerization","stoichiometry":[-2,1,0],"log10K":1.2,"referenceT_K":298.15,"deltaH_Jmol":-12000,"active":true}`}]
            </p>
          </div>
          <div className="grid grid-cols-1 gap-3 text-xs text-muted-foreground md:grid-cols-3">
            <div>Dissociation: {(Number(unit.params.dissociationFraction ?? 0) * 100).toFixed(2)}%</div>
            <div>Ionic strength: {Number(unit.params.ionicStrength ?? 0).toFixed(4)}</div>
            <div>Equilibrium K: {Number(unit.params.equilibriumConstant ?? 0).toExponential(3)}</div>
            <div>Active reactions: {Number(unit.params.activeReactionCount ?? 0)}</div>
            <div>Equilibrium passes: {Number(unit.params.equilibriumPasses ?? 0)}</div>
          </div>
        </div>
      )}

      {isRigorousColumn && (
        <div className="space-y-3 rounded border border-border/60 bg-muted/20 p-3">
          <div className="flex items-center justify-between gap-3">
            <div>
              <h3 className="text-sm font-semibold">Column Specifications</h3>
              <p className="text-xs text-muted-foreground">
                Aspen-style operating specs. Two active specs are required for a rigorous column.
              </p>
            </div>
            <div className={`rounded px-2 py-1 text-xs ${
              specDiagnostics?.status === 'ok'
                ? 'bg-green-100 text-green-800 dark:bg-green-900/40 dark:text-green-300'
                : specDiagnostics?.status === 'underspecified'
                  ? 'bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300'
                  : 'bg-red-100 text-red-800 dark:bg-red-900/40 dark:text-red-300'
            }`}>
              {specDiagnostics?.status ?? 'ok'}
            </div>
          </div>

          {(['spec1', 'spec2'] as const).map(slot => {
            const modeKey = `${slot}Mode` as const;
            const valueKey = `${slot}Value` as const;
            const componentKey = `${slot}ComponentIndex` as const;
            const mode = (unit.params[modeKey] as RigorousColumnSpecMode | undefined) ?? 'None';
            return (
              <div key={slot} className="grid grid-cols-1 gap-3 rounded border border-border/60 bg-background/60 p-3 md:grid-cols-3">
                <div className="space-y-1">
                  <label className="text-xs font-medium text-muted-foreground">{slot === 'spec1' ? 'Spec 1' : 'Spec 2'}</label>
                  <select
                    value={mode}
                    onChange={(e) => updateUnit(unitId, { params: { ...unit.params, [modeKey]: e.target.value as RigorousColumnSpecMode } })}
                    className="w-full px-2 py-1.5 text-sm border rounded bg-background border-input focus:outline-none focus:ring-2 focus:ring-ring"
                  >
                    {COLUMN_SPEC_OPTIONS.map(option => (
                      <option key={option.value} value={option.value}>{option.label}</option>
                    ))}
                  </select>
                </div>
                <FieldInput
                  label="Value"
                  value={Number(unit.params[valueKey] ?? 0)}
                  onChange={(v) => updateUnit(unitId, { params: { ...unit.params, [valueKey]: v } })}
                  locked={mode === 'None'}
                  decimals={mode === 'DistillateComposition' || mode === 'DistillateFraction' ? 4 : 2}
                />
                <div className="space-y-1">
                  <label className="text-xs font-medium text-muted-foreground">Component</label>
                  <select
                    value={Number(unit.params[componentKey] ?? 0)}
                    onChange={(e) => updateUnit(unitId, { params: { ...unit.params, [componentKey]: parseInt(e.target.value, 10) } })}
                    disabled={mode !== 'DistillateComposition'}
                    className="w-full px-2 py-1.5 text-sm border rounded bg-background border-input focus:outline-none focus:ring-2 focus:ring-ring disabled:cursor-not-allowed disabled:bg-muted/30"
                  >
                    {compounds.map((compound, index) => (
                      <option key={compound.name} value={index}>{compound.displayName}</option>
                    ))}
                  </select>
                </div>
              </div>
            );
          })}

          <div className="space-y-1 text-xs">
            {specDiagnostics?.messages.map(message => (
              <div key={message} className="text-muted-foreground">{message}</div>
            ))}
            {!feedStream && (
              <div className="text-muted-foreground">Connect a feed stream to let fraction- and bottoms-based specs resolve against feed flow.</div>
            )}
          </div>
        </div>
      )}

      <div className="space-y-2">
        <h3 className="text-sm font-semibold">Parameters</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
          {editableParams.map(([key, value]) => {
            if (typeof value === 'boolean') {
              return (
                <div key={key} className="flex items-center gap-2">
                  <label className="text-xs font-medium text-muted-foreground">{formatParamLabel(key, us)}</label>
                  <input
                    type="checkbox"
                    checked={value}
                    onChange={(e) => updateUnit(unitId, { params: { ...unit.params, [key]: e.target.checked } })}
                    className="rounded border-input"
                  />
                </div>
              );
            }
            if (typeof value === 'string') {
              return (
                <div key={key} className="space-y-1">
                  <label className="text-xs font-medium text-muted-foreground">{formatParamLabel(key, us)}</label>
                  <input
                    type="text"
                    value={value}
                    onChange={(e) => updateUnit(unitId, { params: { ...unit.params, [key]: e.target.value } })}
                    className="w-full px-2 py-1.5 text-sm border rounded bg-background border-input focus:outline-none focus:ring-2 focus:ring-ring"
                  />
                </div>
              );
            }
            return (
              <FieldInput
                key={key}
                label={formatParamLabel(key, us)}
                value={displayValue(key, value as number)}
                onChange={(v) => {
                  updateUnit(unitId, {
                    params: { ...unit.params, [key]: fromDisplay(key, v) },
                  });
                }}
                locked={false}
                decimals={key === 'splitRatio' || key === 'efficiency' || key.includes('Recovery') ? 3 : 2}
              />
            );
          })}
        </div>
      </div>

      {/* Reaction Editor — for RStoic, CSTR, PFR, REquil, RBatch */}
      {['RStoic', 'CSTR', 'PFR', 'REquil', 'RBatch'].includes(unit.type) && (() => {
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        const rxns = ((unit.params as any).reactions ?? []) as any[];
        const setRxns = (next: any[]) =>
          updateUnit(unitId, { params: { ...unit.params, reactions: next as any } });
        return (
        <div className="space-y-2">
          <h3 className="text-sm font-semibold">Reactions</h3>
          {rxns.length === 0 && (
            <div className="text-xs text-muted-foreground italic">No reactions defined.</div>
          )}
          {rxns.map((rxn: any, ri: number) => (
            <div key={ri} className="border border-border rounded p-2 space-y-1.5 bg-muted/20">
              <div className="flex items-center justify-between">
                <span className="text-xs font-semibold">Reaction {ri + 1}</span>
                <button
                  onClick={() => {
                    const next = [...rxns];
                    next.splice(ri, 1);
                    setRxns(next);
                  }}
                  className="text-red-500 hover:text-red-600 text-xs"
                >
                  <Trash2 className="w-3 h-3" />
                </button>
              </div>
              <div className="text-[10px] text-muted-foreground">
                {(rxn.stoichCoeffs ?? []).map((c: number, ci: number) => {
                  const name = compounds[ci]?.displayName ?? `Comp ${ci}`;
                  if (c === 0) return null;
                  return <span key={ci} className="mr-1">{c > 0 ? '+' : ''}{c} {name}</span>;
                })}
              </div>
              <div className="grid grid-cols-2 gap-2">
                {compounds.map((comp, ci) => (
                  <div key={ci} className="flex items-center gap-1">
                    <label className="text-[10px] text-muted-foreground w-16 truncate" title={comp.displayName}>{comp.displayName}</label>
                    <input
                      type="number"
                      step="0.1"
                      value={rxn.stoichCoeffs?.[ci] ?? 0}
                      onChange={(e) => {
                        const next = [...rxns];
                        const coeffs = [...(next[ri].stoichCoeffs ?? new Array(compounds.length).fill(0))];
                        while (coeffs.length < compounds.length) coeffs.push(0);
                        coeffs[ci] = parseFloat(e.target.value) || 0;
                        next[ri] = { ...next[ri], stoichCoeffs: coeffs };
                        setRxns(next);
                      }}
                      className="w-16 px-1 py-0.5 text-xs border rounded bg-background border-input"
                    />
                  </div>
                ))}
              </div>
              <div className="grid grid-cols-2 gap-2">
                <div className="space-y-0.5">
                  <label className="text-[10px] text-muted-foreground">Key Component</label>
                  <select
                    value={rxn.keyComponentIndex ?? 0}
                    onChange={(e) => {
                      const next = [...rxns];
                      next[ri] = { ...next[ri], keyComponentIndex: parseInt(e.target.value) };
                      setRxns(next);
                    }}
                    className="w-full px-1 py-0.5 text-xs border rounded bg-background border-input"
                  >
                    {compounds.map((comp, ci) => (
                      <option key={ci} value={ci}>{comp.displayName}</option>
                    ))}
                  </select>
                </div>
                <FieldInput
                  label="Fractional Conversion"
                  value={rxn.fractionalConversion ?? 0}
                  onChange={(v) => {
                    const next = [...rxns];
                    next[ri] = { ...next[ri], fractionalConversion: Math.max(0, Math.min(1, v)) };
                    setRxns(next);
                  }}
                  locked={false}
                  decimals={3}
                />
              </div>
            </div>
          ))}
          <Button
            variant="outline"
            size="sm"
            onClick={() => {
              setRxns([...rxns, {
                stoichCoeffs: new Array(compounds.length).fill(0),
                keyComponentIndex: 0,
                fractionalConversion: 0.5,
              }]);
            }}
            className="text-xs"
          >
            + Add Reaction
          </Button>
          {compounds.length === 0 && (
            <div className="text-xs text-amber-600">⚠ No compounds selected. Add compounds in the Setup screen first.</div>
          )}
        </div>
        );
      })()}

      {/* Connected Streams — clickable */}
      <div className="space-y-2">
        <h3 className="text-sm font-semibold">Connected Streams</h3>
        {unit.ports.map(port => {
          const stream = port.streamId ? streams[port.streamId] : null;
          return (
            <button
              key={port.id}
              onClick={() => stream && openTab({ type: 'stream', id: stream.id })}
              className="w-full text-left text-xs border border-border rounded p-2.5 bg-muted/30 hover:bg-accent/50 transition-colors"
            >
              <div className="font-medium flex items-center gap-1.5">
                <span className={`w-2 h-2 rounded-full ${port.type === 'inlet' ? 'bg-blue-500' : 'bg-green-500'}`} />
                {port.label} ({port.type}): {stream?.name ?? '— not connected'}
              </div>
              {stream?.solved && (
                <div className="mt-1 text-muted-foreground grid grid-cols-2 gap-x-4">
                  <span>T: {fmt(fromSI_T(stream.T_K, us.temperature), 1)} {tempLabel(us.temperature)}</span>
                  <span>P: {fmt(fromSI_P(stream.P_Pa, us.pressure), 2)} {pressLabel(us.pressure)}</span>
                  <span>Flow: {fmt(fromSI_F(stream.totalFlow_molps, us.flow), 2)} {flowLabel(us.flow)}</span>
                  <span>Phase: {stream.phase}</span>
                </div>
              )}
            </button>
          );
        })}
      </div>

      {/* Duty */}
      {unit.solved && (
        <div className="space-y-3 border-t border-border pt-4">
          <h3 className="text-sm font-semibold">Results</h3>
          <div className="grid grid-cols-2 md:grid-cols-3 gap-3 text-sm">
            <ResultField label="Heat Duty" value={`${fmt(fromSI_Q(unit.duty_W, us.energy), 2)} ${energyLabel(us.energy)}`} />
            <ResultField label="Status" value={unit.solved ? '✓ Converged' : '✗ Not solved'} />
            {/* Mass balance */}
            {(() => {
              const inPorts = unit.ports.filter(p => p.type === 'inlet' && p.streamId);
              const outPorts = unit.ports.filter(p => p.type === 'outlet' && p.streamId);
              const inFlow = inPorts.reduce((s, p) => {
                const st = streams[p.streamId!];
                return s + (st?.totalFlow_molps ?? 0);
              }, 0);
              const outFlow = outPorts.reduce((s, p) => {
                const st = streams[p.streamId!];
                return s + (st?.totalFlow_molps ?? 0);
              }, 0);
              const inH = inPorts.reduce((s, p) => {
                const st = streams[p.streamId!];
                return s + (st?.totalFlow_molps ?? 0) * (st?.H_Jpmol ?? 0);
              }, 0);
              const outH = outPorts.reduce((s, p) => {
                const st = streams[p.streamId!];
                return s + (st?.totalFlow_molps ?? 0) * (st?.H_Jpmol ?? 0);
              }, 0);
              // Molecular weight average for mass flow
              const avgMW = compounds.reduce((s, c, i) => {
                const inSt = inPorts[0]?.streamId ? streams[inPorts[0].streamId] : null;
                return s + c.molecularWeight * (inSt?.moleFractions[i] ?? 0);
              }, 0) || 1;
              const massFlowIn = inFlow * avgMW / 1000; // kg/s
              return (
                <>
                  <ResultField label={`Inlet Flow (${flowLabel(us.flow)})`} value={fmt(fromSI_F(inFlow, us.flow), 2)} />
                  <ResultField label={`Outlet Flow (${flowLabel(us.flow)})`} value={fmt(fromSI_F(outFlow, us.flow), 2)} />
                  <ResultField label="Mass Flow (kg/s)" value={fmt(massFlowIn, 3)} />
                  <ResultField label="Q_in (J/s)" value={fmt(inH, 0)} />
                  <ResultField label="Q_out (J/s)" value={fmt(outH, 0)} />
                  <ResultField label="Energy Balance" value={`${fmt(outH - inH - unit.duty_W, 1)} W`} />
                </>
              );
            })()}
          </div>

          {/* Equipment cost estimate */}
          {(() => {
            const ec = calculateEquipmentCost(unit, streams);
            const fmtUSD = (v: number) => v >= 1e6 ? `$${(v/1e6).toFixed(2)}M` : v >= 1e3 ? `$${(v/1e3).toFixed(1)}k` : `$${v.toFixed(0)}`;
            return (
              <div className="text-xs text-muted-foreground mt-1 border border-border/50 rounded p-2 bg-muted/20">
                <span className="font-semibold text-foreground">Equipment Cost Estimate: </span>
                Purchase {fmtUSD(ec.purchaseCost)} · Installed {fmtUSD(ec.bareModuleCost)}
                <span className="ml-2">({ec.sizingParam.toFixed(1)} {ec.sizingUnit})</span>
              </div>
            );
          })()}
        </div>
      )}

      {/* Errors */}
      {unit.errors.length > 0 && (
        <div className="text-xs text-red-500 space-y-1">
          {unit.errors.map((e, i) => <div key={i}>⚠ {e}</div>)}
        </div>
      )}
    </div>
  );
}

// ─── Shared Field Input ─────────────────────────────────────────

function FieldInput({
  label,
  value,
  onChange,
  locked,
  decimals = 2,
  min,
  max,
  step,
}: {
  label: React.ReactNode;
  value: number;
  onChange: (v: number) => void;
  locked: boolean;
  decimals?: number;
  min?: number;
  max?: number;
  step?: number;
}) {
  return (
    <div className="space-y-1">
      <label className="text-xs font-medium text-muted-foreground flex items-center gap-1">
        {locked && <span title="Determined by upstream unit">🔒</span>}
        {label}
      </label>
      <input
        type="number"
        value={parseFloat(value.toFixed(decimals))}
        onChange={(e) => onChange(parseFloat(e.target.value) || 0)}
        disabled={locked}
        min={min}
        max={max}
        step={step ?? Math.pow(10, -decimals)}
        className={`w-full px-2 py-1.5 text-sm border rounded bg-background
          ${locked
            ? 'border-muted text-muted-foreground cursor-not-allowed bg-muted/30'
            : 'border-input focus:outline-none focus:ring-2 focus:ring-ring'
          }`}
      />
    </div>
  );
}

function ResultField({ label, value }: { label: React.ReactNode; value: string }) {
  return (
    <div>
      <div className="text-xs text-muted-foreground">{label}</div>
      <div className="text-sm font-medium">{value}</div>
    </div>
  );
}

// ─── Bottom Panel with Tabs ───────────────────────────────────

function BottomPanelTabs({ onClose }: { onClose?: () => void }) {
  const [activeBottomTab, setActiveBottomTab] = useState<'streams' | 'utilities' | 'economics' | 'designSpecs' | 'sensitivity' | 'analysis' | 'calculators'>('streams');

  const tabClass = (tab: string) =>
    `px-3 py-1.5 text-xs font-medium transition-colors border-b-2 ${
      activeBottomTab === tab
        ? 'border-blue-500 text-blue-600 dark:text-blue-400'
        : 'border-transparent text-muted-foreground hover:text-foreground'
    }`;

  return (
    <div>
      <div className="flex items-center gap-0 border-b border-border bg-muted/30 px-2">
        <button onClick={() => setActiveBottomTab('streams')} className={tabClass('streams')}>
          Stream Results
        </button>
        <button onClick={() => setActiveBottomTab('utilities')} className={tabClass('utilities')}>
          Utility Summary
        </button>
        <button onClick={() => setActiveBottomTab('economics')} className={tabClass('economics')}>
          Economics
        </button>
        <button onClick={() => setActiveBottomTab('designSpecs')} className={tabClass('designSpecs')}>
          Design Specs
        </button>
        <button onClick={() => setActiveBottomTab('sensitivity')} className={tabClass('sensitivity')}>
          Sensitivity
        </button>
        <button onClick={() => setActiveBottomTab('analysis')} className={tabClass('analysis')}>
          Property Analysis
        </button>
        <button onClick={() => setActiveBottomTab('calculators')} className={tabClass('calculators')}>
          Calculators
        </button>
        {onClose && (
          <button
            onClick={onClose}
            className="ml-auto p-1 text-muted-foreground hover:text-foreground hover:bg-accent rounded"
            title="Close panel"
          >
            <X className="w-3.5 h-3.5" />
          </button>
        )}
      </div>
      {activeBottomTab === 'streams' && <StreamResultsTable />}
      {activeBottomTab === 'utilities' && <UtilitySummary />}
      {activeBottomTab === 'economics' && <CostingPanel />}
      {activeBottomTab === 'designSpecs' && <DesignSpecsPanel />}
      {activeBottomTab === 'sensitivity' && <SensitivityPanel />}
      {activeBottomTab === 'analysis' && <AnalysisPanel />}
      {activeBottomTab === 'calculators' && <CalculatorBlockPanel />}
    </div>
  );
}

// ─── Utility Summary ────────────────────────────────────────────

function UtilitySummary() {
  const { units, unitSystem: us, openTab } = useCanopyStore();
  const unitList = Object.values(units).filter(u => u.solved);

  const heaters = unitList.filter(u => u.duty_W > 0);
  const coolers = unitList.filter(u => u.duty_W < 0);
  const totalHeating = heaters.reduce((sum, u) => sum + u.duty_W, 0);
  const totalCooling = coolers.reduce((sum, u) => sum + u.duty_W, 0);
  const netDuty = totalHeating + totalCooling;

  return (
    <div className="p-3">
      {/* Summary cards */}
      <div className="grid grid-cols-3 gap-3 mb-4">
        <div className="rounded-md border border-red-500/30 bg-red-500/5 px-3 py-2">
          <div className="text-[10px] text-muted-foreground uppercase tracking-wider mb-1">Total Heating</div>
          <div className="text-sm font-semibold text-red-600 dark:text-red-400">
            {fmt(fromSI_Q(totalHeating, us.energy), 2)} {energyLabel(us.energy)}
          </div>
        </div>
        <div className="rounded-md border border-blue-500/30 bg-blue-500/5 px-3 py-2">
          <div className="text-[10px] text-muted-foreground uppercase tracking-wider mb-1">Total Cooling</div>
          <div className="text-sm font-semibold text-blue-600 dark:text-blue-400">
            {fmt(fromSI_Q(totalCooling, us.energy), 2)} {energyLabel(us.energy)}
          </div>
        </div>
        <div className="rounded-md border border-border px-3 py-2">
          <div className="text-[10px] text-muted-foreground uppercase tracking-wider mb-1">Net Duty</div>
          <div className={`text-sm font-semibold ${netDuty >= 0 ? 'text-red-600 dark:text-red-400' : 'text-blue-600 dark:text-blue-400'}`}>
            {fmt(fromSI_Q(netDuty, us.energy), 2)} {energyLabel(us.energy)}
          </div>
        </div>
      </div>

      {/* Per-unit breakdown table */}
      <div className="overflow-x-auto">
        <table className="text-xs w-full border-collapse">
          <thead>
            <tr className="border-b border-border">
              <th className="text-left p-1.5 font-semibold text-muted-foreground">Unit</th>
              <th className="text-left p-1.5 font-semibold text-muted-foreground">Type</th>
              <th className="text-right p-1.5 font-semibold text-muted-foreground">Duty ({energyLabel(us.energy)})</th>
              <th className="text-left p-1.5 font-semibold text-muted-foreground">Classification</th>
            </tr>
          </thead>
          <tbody>
            {unitList.map(u => (
              <tr key={u.id} className="border-b border-border/50">
                <td className="p-1.5">
                  <button
                    onClick={() => openTab({ type: 'unit', id: u.id })}
                    className="hover:text-blue-500 transition-colors underline decoration-dotted"
                  >
                    {u.name}
                  </button>
                </td>
                <td className="p-1.5 text-muted-foreground">{u.type}</td>
                <td className={`p-1.5 text-right font-mono ${u.duty_W > 0 ? 'text-red-500' : u.duty_W < 0 ? 'text-blue-500' : ''}`}>
                  {fmt(fromSI_Q(u.duty_W, us.energy), 2)}
                </td>
                <td className="p-1.5 text-muted-foreground">
                  {u.duty_W > 100 ? 'Heating' : u.duty_W < -100 ? 'Cooling' : 'Adiabatic'}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

// ─── Economics / Costing Panel ──────────────────────────────────

function CostingPanel() {
  const { units, streams } = useCanopyStore();
  const [rates, setRates] = useState<UtilityRates>({ ...DEFAULT_UTILITY_RATES });

  const econ: EconomicsSummary = useMemo(
    () => calculateEconomics(units, streams, rates),
    [units, streams, rates],
  );

  const fmtUSD = (v: number) => {
    if (v >= 1e6) return `$${(v / 1e6).toFixed(2)}M`;
    if (v >= 1e3) return `$${(v / 1e3).toFixed(1)}k`;
    return `$${v.toFixed(0)}`;
  };

  return (
    <div className="p-3 space-y-4">
      {/* Summary Cards */}
      <div className="grid grid-cols-4 gap-3 text-xs">
        <div className="bg-blue-500/10 border border-blue-500/20 rounded p-2 text-center">
          <div className="text-muted-foreground mb-0.5">Total Equipment Purchase</div>
          <div className="text-sm font-bold text-blue-500">{fmtUSD(econ.totalPurchaseCost)}</div>
        </div>
        <div className="bg-green-500/10 border border-green-500/20 rounded p-2 text-center">
          <div className="text-muted-foreground mb-0.5">Total Bare Module (Installed)</div>
          <div className="text-sm font-bold text-green-500">{fmtUSD(econ.totalBareModuleCost)}</div>
        </div>
        <div className="bg-purple-500/10 border border-purple-500/20 rounded p-2 text-center">
          <div className="text-muted-foreground mb-0.5">Grass Roots Capital</div>
          <div className="text-sm font-bold text-purple-500">{fmtUSD(econ.totalGrassRootsCost)}</div>
        </div>
        <div className="bg-orange-500/10 border border-orange-500/20 rounded p-2 text-center">
          <div className="text-muted-foreground mb-0.5">Annual Utility Cost</div>
          <div className="text-sm font-bold text-orange-500">{fmtUSD(econ.totalAnnualUtilityCost)}/yr</div>
        </div>
      </div>

      {/* Equipment Cost Table */}
      <div>
        <h4 className="text-xs font-semibold mb-1">Equipment Costs</h4>
        <table className="w-full text-xs border border-border">
          <thead className="bg-muted/50">
            <tr>
              <th className="text-left px-2 py-1 border-b border-border">Unit</th>
              <th className="text-left px-2 py-1 border-b border-border">Type</th>
              <th className="text-right px-2 py-1 border-b border-border">Sizing</th>
              <th className="text-right px-2 py-1 border-b border-border">Purchase Cost</th>
              <th className="text-right px-2 py-1 border-b border-border">Installed Cost</th>
            </tr>
          </thead>
          <tbody>
            {econ.equipmentCosts.map(ec => (
              <tr key={ec.unitId} className="hover:bg-muted/30">
                <td className="px-2 py-0.5 border-b border-border">{ec.unitName}</td>
                <td className="px-2 py-0.5 border-b border-border text-muted-foreground">{ec.unitType}</td>
                <td className="px-2 py-0.5 border-b border-border text-right font-mono">
                  {ec.sizingParam.toFixed(2)} {ec.sizingUnit}
                </td>
                <td className="px-2 py-0.5 border-b border-border text-right font-mono">{fmtUSD(ec.purchaseCost)}</td>
                <td className="px-2 py-0.5 border-b border-border text-right font-mono font-semibold">{fmtUSD(ec.bareModuleCost)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Utility Cost Table */}
      <div>
        <h4 className="text-xs font-semibold mb-1">Operating / Utility Costs ({econ.operatingHours} hr/yr)</h4>
        <table className="w-full text-xs border border-border">
          <thead className="bg-muted/50">
            <tr>
              <th className="text-left px-2 py-1 border-b border-border">Unit</th>
              <th className="text-right px-2 py-1 border-b border-border">Duty</th>
              <th className="text-left px-2 py-1 border-b border-border">Utility</th>
              <th className="text-right px-2 py-1 border-b border-border">Rate</th>
              <th className="text-right px-2 py-1 border-b border-border">Annual Cost</th>
            </tr>
          </thead>
          <tbody>
            {econ.utilityCosts.map(uc => (
              <tr key={uc.unitId} className="hover:bg-muted/30">
                <td className="px-2 py-0.5 border-b border-border">{uc.unitName}</td>
                <td className="px-2 py-0.5 border-b border-border text-right font-mono">
                  {(uc.duty_W / 1000).toFixed(1)} kW
                </td>
                <td className="px-2 py-0.5 border-b border-border">
                  <span className={uc.duty_W > 0 ? 'text-orange-500' : 'text-blue-500'}>
                    {uc.utilityType === 'electricity' ? '⚡ Electricity'
                      : uc.utilityType === 'coolingWater' ? '💧 Cooling Water'
                      : uc.utilityType === 'lpSteam' ? '♨️ LP Steam'
                      : uc.utilityType === 'mpSteam' ? '♨️ MP Steam'
                      : uc.utilityType === 'hpSteam' ? '♨️ HP Steam'
                      : uc.utilityType === 'refrigeration' ? '❄️ Refrigeration'
                      : '🔥 Fuel'}
                  </span>
                </td>
                <td className="px-2 py-0.5 border-b border-border text-right font-mono text-muted-foreground">
                  ${uc.rateUsed.toFixed(3)}/{uc.utilityType === 'electricity' ? 'kWh' : 'GJ'}
                </td>
                <td className="px-2 py-0.5 border-b border-border text-right font-mono font-semibold">{fmtUSD(uc.annualCost)}/yr</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Editable utility rates */}
      <details className="text-xs">
        <summary className="cursor-pointer text-muted-foreground hover:text-foreground">
          Edit Utility Rates
        </summary>
        <div className="grid grid-cols-4 gap-2 mt-2">
          {(Object.keys(rates) as Array<keyof UtilityRates>).map(key => (
            <label key={key} className="flex flex-col gap-0.5">
              <span className="text-muted-foreground capitalize">{key.replace(/([A-Z])/g, ' $1')}</span>
              <input
                type="number"
                step="0.01"
                value={rates[key]}
                onChange={e => setRates(prev => ({ ...prev, [key]: parseFloat(e.target.value) || 0 }))}
                className="w-full px-1.5 py-0.5 border border-border rounded bg-background text-xs font-mono"
              />
            </label>
          ))}
        </div>
      </details>
    </div>
  );
}

// ─── Design Specs Panel ─────────────────────────────────────────

function DesignSpecsPanel() {
  const { designSpecs, addDesignSpec, updateDesignSpec, removeDesignSpec, streams, units } = useCanopyStore();
  const [showAdd, setShowAdd] = useState(false);
  const [newSpec, setNewSpec] = useState({
    name: '', targetType: 'stream' as 'stream' | 'unit', targetId: '', targetProp: '', targetVal: '',
    manipUnitId: '', manipParam: '', lower: '', upper: '', tolerance: '0.001',
  });

  const streamList = Object.values(streams);
  const unitList = Object.values(units);

  const handleAdd = () => {
    if (!newSpec.name || !newSpec.targetId || !newSpec.targetProp || !newSpec.targetVal || !newSpec.manipUnitId || !newSpec.manipParam) return;
    addDesignSpec({
      id: `DS-${Date.now()}`,
      name: newSpec.name,
      target: {
        objectType: newSpec.targetType,
        objectId: newSpec.targetId,
        property: newSpec.targetProp,
        value: parseFloat(newSpec.targetVal),
      },
      manipulated: {
        unitId: newSpec.manipUnitId,
        paramName: newSpec.manipParam,
        lowerBound: parseFloat(newSpec.lower) || 0,
        upperBound: parseFloat(newSpec.upper) || 1e6,
      },
      tolerance: parseFloat(newSpec.tolerance) || 0.001,
      active: true,
    });
    setShowAdd(false);
    setNewSpec({ name: '', targetType: 'stream', targetId: '', targetProp: '', targetVal: '', manipUnitId: '', manipParam: '', lower: '', upper: '', tolerance: '0.001' });
  };

  const streamProps = ['T_K', 'P_Pa', 'totalFlow_molps', 'vaporFraction', 'H_Jpmol'];
  const unitProps = ['duty_W'];

  return (
    <div className="p-3">
      <div className="flex items-center justify-between mb-3">
        <button onClick={() => setShowAdd(!showAdd)} className="text-xs px-2 py-1 rounded bg-blue-600 text-white hover:bg-blue-700">
          {showAdd ? 'Cancel' : '+ Add Spec'}
        </button>
      </div>

      {showAdd && (
        <div className="grid grid-cols-2 md:grid-cols-4 gap-2 mb-3 p-3 border border-border rounded-md bg-muted/30">
          <input placeholder="Spec Name" value={newSpec.name} onChange={e => setNewSpec({ ...newSpec, name: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background" />
          <select value={newSpec.targetType} onChange={e => setNewSpec({ ...newSpec, targetType: e.target.value as any })} className="px-2 py-1 text-xs border border-input rounded bg-background">
            <option value="stream">Target: Stream</option>
            <option value="unit">Target: Unit</option>
          </select>
          <select value={newSpec.targetId} onChange={e => setNewSpec({ ...newSpec, targetId: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background">
            <option value="">Select target...</option>
            {(newSpec.targetType === 'stream' ? streamList : unitList).map(o => (
              <option key={o.id} value={o.id}>{o.name}</option>
            ))}
          </select>
          <select value={newSpec.targetProp} onChange={e => setNewSpec({ ...newSpec, targetProp: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background">
            <option value="">Property...</option>
            {(newSpec.targetType === 'stream' ? streamProps : unitProps).map(p => (
              <option key={p} value={p}>{p}</option>
            ))}
          </select>
          <input placeholder="Target value" value={newSpec.targetVal} onChange={e => setNewSpec({ ...newSpec, targetVal: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background" type="number" />
          <select value={newSpec.manipUnitId} onChange={e => setNewSpec({ ...newSpec, manipUnitId: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background">
            <option value="">Manipulated unit...</option>
            {unitList.map(u => (
              <option key={u.id} value={u.id}>{u.name}</option>
            ))}
          </select>
          <input placeholder="Param name (e.g. targetT_K)" value={newSpec.manipParam} onChange={e => setNewSpec({ ...newSpec, manipParam: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background" />
          <div className="flex gap-1">
            <input placeholder="Lower" value={newSpec.lower} onChange={e => setNewSpec({ ...newSpec, lower: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background w-1/2" type="number" />
            <input placeholder="Upper" value={newSpec.upper} onChange={e => setNewSpec({ ...newSpec, upper: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background w-1/2" type="number" />
          </div>
          <button onClick={handleAdd} className="col-span-2 md:col-span-4 text-xs px-2 py-1.5 rounded bg-green-600 text-white hover:bg-green-700">
            Add Design Spec
          </button>
        </div>
      )}

      {designSpecs.length === 0 ? (
        <p className="text-xs text-muted-foreground italic">No design specifications defined. Design specs let you specify a target value (e.g., stream temperature) and the solver will adjust a unit parameter to meet it.</p>
      ) : (
        <table className="text-xs w-full border-collapse">
          <thead>
            <tr className="border-b border-border">
              <th className="text-left p-1.5 font-semibold text-muted-foreground">Name</th>
              <th className="text-left p-1.5 font-semibold text-muted-foreground">Target</th>
              <th className="text-right p-1.5 font-semibold text-muted-foreground">Value</th>
              <th className="text-left p-1.5 font-semibold text-muted-foreground">Manipulated</th>
              <th className="text-center p-1.5 font-semibold text-muted-foreground">Active</th>
              <th className="p-1.5"></th>
            </tr>
          </thead>
          <tbody>
            {designSpecs.map(ds => (
              <tr key={ds.id} className="border-b border-border/50">
                <td className="p-1.5 font-medium">{ds.name}</td>
                <td className="p-1.5 text-muted-foreground">{ds.target.objectId}.{ds.target.property}</td>
                <td className="p-1.5 text-right font-mono">{ds.target.value}</td>
                <td className="p-1.5 text-muted-foreground">{ds.manipulated.unitId}.{ds.manipulated.paramName}</td>
                <td className="p-1.5 text-center">
                  <input type="checkbox" checked={ds.active} onChange={e => updateDesignSpec(ds.id, { active: e.target.checked })} />
                </td>
                <td className="p-1.5">
                  <button onClick={() => removeDesignSpec(ds.id)} className="text-red-400 hover:text-red-600 text-xs">
                    <Trash2 className="w-3 h-3" />
                  </button>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      )}
    </div>
  );
}

// ─── Sensitivity Analysis Panel ─────────────────────────────────

function SensitivityPanel() {
  const { sensitivities, addSensitivity, removeSensitivity, runSensitivity, units, streams } = useCanopyStore();
  const [showAdd, setShowAdd] = useState(false);
  const [newSA, setNewSA] = useState({
    name: '', unitId: '', paramName: '', from: '', to: '', steps: '10',
    outType: 'stream' as 'stream' | 'unit', outId: '', outProp: '', outLabel: '',
  });

  const unitList = Object.values(units);
  const streamList = Object.values(streams);

  const handleAdd = () => {
    if (!newSA.name || !newSA.unitId || !newSA.paramName || !newSA.from || !newSA.to || !newSA.outId || !newSA.outProp) return;
    const from = parseFloat(newSA.from);
    const to = parseFloat(newSA.to);
    const steps = parseInt(newSA.steps) || 10;
    const values = Array.from({ length: steps + 1 }, (_, i) => from + (to - from) * i / steps);

    addSensitivity({
      id: `SA-${Date.now()}`,
      name: newSA.name,
      variable: { unitId: newSA.unitId, paramName: newSA.paramName, values },
      outputs: [{
        objectType: newSA.outType,
        objectId: newSA.outId,
        property: newSA.outProp,
        label: newSA.outLabel || `${newSA.outId}.${newSA.outProp}`,
      }],
      results: [],
      active: true,
    });
    setShowAdd(false);
  };

  const outputProps = ['T_K', 'P_Pa', 'totalFlow_molps', 'vaporFraction', 'H_Jpmol', 'duty_W'];

  return (
    <div className="p-3">
      <div className="flex items-center justify-between mb-3">
        <button onClick={() => setShowAdd(!showAdd)} className="text-xs px-2 py-1 rounded bg-blue-600 text-white hover:bg-blue-700">
          {showAdd ? 'Cancel' : '+ Add Study'}
        </button>
      </div>

      {showAdd && (
        <div className="grid grid-cols-2 md:grid-cols-4 gap-2 mb-3 p-3 border border-border rounded-md bg-muted/30">
          <input placeholder="Study name" value={newSA.name} onChange={e => setNewSA({ ...newSA, name: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background" />
          <select value={newSA.unitId} onChange={e => setNewSA({ ...newSA, unitId: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background">
            <option value="">Vary unit...</option>
            {unitList.map(u => <option key={u.id} value={u.id}>{u.name}</option>)}
          </select>
          <input placeholder="Parameter (e.g. targetT_K)" value={newSA.paramName} onChange={e => setNewSA({ ...newSA, paramName: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background" />
          <div className="flex gap-1">
            <input placeholder="From" value={newSA.from} onChange={e => setNewSA({ ...newSA, from: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background w-1/3" type="number" />
            <input placeholder="To" value={newSA.to} onChange={e => setNewSA({ ...newSA, to: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background w-1/3" type="number" />
            <input placeholder="Steps" value={newSA.steps} onChange={e => setNewSA({ ...newSA, steps: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background w-1/3" type="number" />
          </div>
          <select value={newSA.outType} onChange={e => setNewSA({ ...newSA, outType: e.target.value as any })} className="px-2 py-1 text-xs border border-input rounded bg-background">
            <option value="stream">Output: Stream</option>
            <option value="unit">Output: Unit</option>
          </select>
          <select value={newSA.outId} onChange={e => setNewSA({ ...newSA, outId: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background">
            <option value="">Select output...</option>
            {(newSA.outType === 'stream' ? streamList : unitList).map(o => (
              <option key={o.id} value={o.id}>{o.name}</option>
            ))}
          </select>
          <select value={newSA.outProp} onChange={e => setNewSA({ ...newSA, outProp: e.target.value })} className="px-2 py-1 text-xs border border-input rounded bg-background">
            <option value="">Property...</option>
            {outputProps.map(p => <option key={p} value={p}>{p}</option>)}
          </select>
          <button onClick={handleAdd} className="text-xs px-2 py-1.5 rounded bg-green-600 text-white hover:bg-green-700">
            Add Study
          </button>
        </div>
      )}

      {sensitivities.length === 0 ? (
        <p className="text-xs text-muted-foreground italic">No sensitivity studies defined. Sweep a unit parameter across a range and observe how outputs change.</p>
      ) : (
        <div className="space-y-4">
          {sensitivities.map(sa => (
            <div key={sa.id} className="border border-border rounded-md p-3">
              <div className="flex items-center justify-between mb-2">
                <span className="text-sm font-semibold">{sa.name}</span>
                <div className="flex items-center gap-2">
                  <button
                    onClick={() => runSensitivity(sa.id)}
                    className="text-xs px-2 py-1 rounded bg-green-600 text-white hover:bg-green-700"
                  >
                    ▶ Run
                  </button>
                  <button onClick={() => removeSensitivity(sa.id)} className="text-red-400 hover:text-red-600">
                    <Trash2 className="w-3 h-3" />
                  </button>
                </div>
              </div>
              <p className="text-xs text-muted-foreground mb-2">
                Vary {sa.variable.unitId}.{sa.variable.paramName} from {sa.variable.values[0]} to {sa.variable.values[sa.variable.values.length - 1]} ({sa.variable.values.length} points)
                → {sa.outputs.map(o => o.label).join(', ')}
              </p>
              {sa.results.length > 0 && (
                <div className="overflow-x-auto">
                  <table className="text-xs w-full border-collapse">
                    <thead>
                      <tr className="border-b border-border">
                        <th className="text-right p-1 font-semibold text-muted-foreground">{sa.variable.paramName}</th>
                        {sa.outputs.map((o, i) => (
                          <th key={i} className="text-right p-1 font-semibold text-muted-foreground">{o.label}</th>
                        ))}
                      </tr>
                    </thead>
                    <tbody>
                      {sa.results.map((row, ri) => (
                        <tr key={ri} className="border-b border-border/50">
                          {row.map((val, ci) => (
                            <td key={ci} className="p-1 text-right font-mono">{isNaN(val) ? '—' : fmt(val, 3)}</td>
                          ))}
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}
            </div>
          ))}
        </div>
      )}
    </div>
  );
}

// ─── Calculator Block UI ────────────────────────────────────────

function CalculatorBlockPanel() {
  const { calculatorBlocks, addCalculatorBlock, removeCalculatorBlock, updateCalculatorBlock,
    streams, units } = useCanopyStore();
  const [editIdx, setEditIdx] = useState<number | null>(null);

  const handleAdd = () => {
    const id = `calc-${Date.now()}`;
    addCalculatorBlock({
      id,
      name: `Calc-${calculatorBlocks.length + 1}`,
      imports: [],
      exports: [],
      expression: '// result = imported_var * 2',
      executionOrder: 'before',
      active: true,
    });
    setEditIdx(calculatorBlocks.length);
  };

  const editing = editIdx !== null ? calculatorBlocks[editIdx] : null;

  return (
    <div className="p-3 space-y-3">
      <div className="flex items-center justify-between">
        <Button variant="outline" size="sm" onClick={handleAdd} className="text-xs h-7">+ Add Block</Button>
      </div>
      <div className="flex gap-4">
        {/* Block list */}
        <div className="w-48 space-y-1">
          {calculatorBlocks.length === 0 && (
            <div className="text-xs text-muted-foreground">No calculator blocks yet.</div>
          )}
          {calculatorBlocks.map((cb, i) => (
            <div key={cb.id}
              onClick={() => setEditIdx(i)}
              className={`flex items-center justify-between px-2 py-1 text-xs border rounded cursor-pointer ${
                editIdx === i ? 'border-primary bg-primary/10' : 'border-border hover:bg-muted'
              }`}>
              <span className="truncate">{cb.name}</span>
              <div className="flex items-center gap-1">
                <button
                  onClick={(e) => { e.stopPropagation(); updateCalculatorBlock(cb.id, { active: !cb.active }); }}
                  className={`w-4 h-4 rounded text-[9px] ${cb.active ? 'bg-green-500 text-white' : 'bg-muted text-muted-foreground'}`}>
                  {cb.active ? '✓' : '✗'}
                </button>
                <button
                  onClick={(e) => { e.stopPropagation(); removeCalculatorBlock(cb.id); if (editIdx === i) setEditIdx(null); }}
                  className="w-4 h-4 rounded bg-destructive/20 text-destructive text-[9px] hover:bg-destructive/40">✕</button>
              </div>
            </div>
          ))}
        </div>
        {/* Editor */}
        {editing && editIdx !== null && (
          <div className="flex-1 space-y-2 min-w-0">
            <div className="flex items-center gap-2">
              <label className="text-xs text-muted-foreground">Name</label>
              <input value={editing.name}
                onChange={e => updateCalculatorBlock(editing.id, { name: e.target.value })}
                className="flex-1 px-1.5 py-1 border border-border rounded bg-background text-xs" />
              <label className="text-xs text-muted-foreground">Order</label>
              <select value={editing.executionOrder}
                onChange={e => updateCalculatorBlock(editing.id, { executionOrder: e.target.value as 'before' | 'after' })}
                className="px-1.5 py-1 border border-border rounded bg-background text-xs">
                <option value="before">Before Solve</option>
                <option value="after">After Solve</option>
              </select>
            </div>
            <div className="grid grid-cols-2 gap-2">
              <div>
                <div className="text-xs text-muted-foreground mb-1">Imports (stream/unit → variable)</div>
                {editing.imports.map((imp, ii) => (
                  <div key={ii} className="flex items-center gap-1 mb-1">
                    <select value={imp.objectId}
                      onChange={e => {
                        const newImps = [...editing.imports];
                        const isStream = Object.values(streams).some((s: MaterialStream) => s.id === e.target.value);
                        newImps[ii] = { ...imp, objectId: e.target.value, objectType: isStream ? 'stream' : 'unit' };
                        updateCalculatorBlock(editing.id, { imports: newImps });
                      }}
                      className="w-24 px-1 py-0.5 border border-border rounded bg-background text-[10px]">
                      <option value="">Source</option>
                      {Object.values(streams).map((s: MaterialStream) => <option key={s.id} value={s.id}>S: {s.name}</option>)}
                      {Object.entries(units).map(([uid, u]) => <option key={uid} value={uid}>U: {u.name}</option>)}
                    </select>
                    <input value={imp.property} placeholder="property"
                      onChange={e => { const newImps = [...editing.imports]; newImps[ii] = { ...imp, property: e.target.value }; updateCalculatorBlock(editing.id, { imports: newImps }); }}
                      className="w-20 px-1 py-0.5 border border-border rounded bg-background text-[10px]" />
                    <span className="text-[10px]">→</span>
                    <input value={imp.varName} placeholder="var"
                      onChange={e => { const newImps = [...editing.imports]; newImps[ii] = { ...imp, varName: e.target.value }; updateCalculatorBlock(editing.id, { imports: newImps }); }}
                      className="w-16 px-1 py-0.5 border border-border rounded bg-background text-[10px]" />
                    <button onClick={() => { const newImps = editing.imports.filter((_: unknown, j: number) => j !== ii); updateCalculatorBlock(editing.id, { imports: newImps }); }}
                      className="text-[10px] text-destructive">✕</button>
                  </div>
                ))}
                <button onClick={() => updateCalculatorBlock(editing.id, { imports: [...editing.imports, { varName: 'T_in', objectType: 'stream' as const, objectId: '', property: 'T' }] })}
                  className="text-[10px] text-primary hover:underline">+ Import</button>
              </div>
              <div>
                <div className="text-xs text-muted-foreground mb-1">Exports (variable → unit param)</div>
                {editing.exports.map((exp, ei) => (
                  <div key={ei} className="flex items-center gap-1 mb-1">
                    <input value={exp.varName} placeholder="var"
                      onChange={e => { const newExps = [...editing.exports]; newExps[ei] = { ...exp, varName: e.target.value }; updateCalculatorBlock(editing.id, { exports: newExps }); }}
                      className="w-16 px-1 py-0.5 border border-border rounded bg-background text-[10px]" />
                    <span className="text-[10px]">→</span>
                    <select value={exp.unitId}
                      onChange={e => { const newExps = [...editing.exports]; newExps[ei] = { ...exp, unitId: e.target.value }; updateCalculatorBlock(editing.id, { exports: newExps }); }}
                      className="w-24 px-1 py-0.5 border border-border rounded bg-background text-[10px]">
                      <option value="">Target Unit</option>
                      {Object.entries(units).map(([uid, u]) => <option key={uid} value={uid}>{u.name}</option>)}
                    </select>
                    <input value={exp.paramName} placeholder="param"
                      onChange={e => { const newExps = [...editing.exports]; newExps[ei] = { ...exp, paramName: e.target.value }; updateCalculatorBlock(editing.id, { exports: newExps }); }}
                      className="w-20 px-1 py-0.5 border border-border rounded bg-background text-[10px]" />
                    <button onClick={() => { const newExps = editing.exports.filter((_: unknown, j: number) => j !== ei); updateCalculatorBlock(editing.id, { exports: newExps }); }}
                      className="text-[10px] text-destructive">✕</button>
                  </div>
                ))}
                <button onClick={() => updateCalculatorBlock(editing.id, { exports: [...editing.exports, { varName: 'result', unitId: '', paramName: '' }] })}
                  className="text-[10px] text-primary hover:underline">+ Export</button>
              </div>
            </div>
            <div>
              <div className="text-xs text-muted-foreground mb-1">Expression (JavaScript)</div>
              <textarea value={editing.expression}
                onChange={e => updateCalculatorBlock(editing.id, { expression: e.target.value })}
                rows={4}
                className="w-full px-2 py-1 border border-border rounded bg-background text-xs font-mono resize-y"
                placeholder="// Imported variables are available by name&#10;// Set result = ... to export" />
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

// ─── Property Analysis Panel ────────────────────────────────────

function AnalysisPanel() {
  const { compounds, fluidPackage, interactionParams, unitSystem: us } = useCanopyStore();
  const [comp1, setComp1] = useState(0);
  const [comp2, setComp2] = useState(Math.min(1, compounds.length - 1));
  const [chartType, setChartType] = useState<'Txy' | 'Pxy' | 'xy'>('Txy');
  const [pressure, setPressure] = useState(101325);
  const [temperature, setTemperature] = useState(350);
  const [data, setData] = useState<TxyData | PxyData | XYData | null>(null);

  const generate = () => {
    if (compounds.length < 2) return;
    try {
      if (chartType === 'Txy') {
        setData(generateTxy(compounds, comp1, comp2, pressure, fluidPackage, interactionParams));
      } else if (chartType === 'Pxy') {
        setData(generatePxy(compounds, comp1, comp2, temperature, fluidPackage, interactionParams));
      } else {
        setData(generateXY(compounds, comp1, comp2, pressure, fluidPackage, interactionParams));
      }
    } catch {
      setData(null);
    }
  };

  if (compounds.length < 2) {
    return <div className="p-3 text-xs text-muted-foreground">Need at least 2 compounds for property analysis.</div>;
  }

  // SVG chart rendering
  const renderChart = () => {
    if (!data) return null;
    const W = 400, H = 280, ML = 55, MR = 20, MT = 20, MB = 35;
    const pw = W - ML - MR, ph = H - MT - MB;

    if (chartType === 'Txy' && 'bubbleCurve' in data) {
      const txy = data as TxyData;
      const allT = [...txy.bubbleCurve.map(p => p.T_K), ...txy.dewCurve.map(p => p.T_K)];
      const Tmin = Math.min(...allT), Tmax = Math.max(...allT);
      const Trange = Tmax - Tmin || 10;
      const toX = (x: number) => ML + x * pw;
      const toY = (T: number) => MT + ph - ((T - Tmin) / Trange) * ph;

      const bubblePath = txy.bubbleCurve.map((p, i) => `${i === 0 ? 'M' : 'L'} ${toX(p.x1).toFixed(1)} ${toY(p.T_K).toFixed(1)}`).join(' ');
      const dewPath = txy.dewCurve.map((p, i) => `${i === 0 ? 'M' : 'L'} ${toX(p.x1).toFixed(1)} ${toY(p.T_K).toFixed(1)}`).join(' ');

      return (
        <svg width={W} height={H} className="bg-card border border-border rounded">
          {/* Grid */}
          {[0, 0.25, 0.5, 0.75, 1].map(v => (
            <line key={`gx-${v}`} x1={toX(v)} y1={MT} x2={toX(v)} y2={MT + ph} stroke="var(--border)" strokeWidth="0.5" />
          ))}
          {Array.from({ length: 5 }, (_, i) => Tmin + (i / 4) * Trange).map(T => (
            <g key={`gy-${T}`}>
              <line x1={ML} y1={toY(T)} x2={ML + pw} y2={toY(T)} stroke="var(--border)" strokeWidth="0.5" />
              <text x={ML - 4} y={toY(T) + 3} textAnchor="end" className="fill-muted-foreground" fontSize="9">
                {fromSI_T(T, us.temperature).toFixed(1)}
              </text>
            </g>
          ))}
          {/* Axes labels */}
          <text x={ML + pw / 2} y={H - 4} textAnchor="middle" className="fill-foreground" fontSize="10">
            x, y ({txy.comp1})
          </text>
          <text x={12} y={MT + ph / 2} textAnchor="middle" className="fill-foreground" fontSize="10"
            transform={`rotate(-90, 12, ${MT + ph / 2})`}>
            T ({tempLabel(us.temperature)})
          </text>
          {/* Curves */}
          <path d={bubblePath} fill="none" stroke="#3b82f6" strokeWidth="2" />
          <path d={dewPath} fill="none" stroke="#ef4444" strokeWidth="2" />
          {/* Legend */}
          <line x1={ML + 10} y1={MT + 8} x2={ML + 25} y2={MT + 8} stroke="#3b82f6" strokeWidth="2" />
          <text x={ML + 28} y={MT + 11} className="fill-foreground" fontSize="9">Bubble</text>
          <line x1={ML + 70} y1={MT + 8} x2={ML + 85} y2={MT + 8} stroke="#ef4444" strokeWidth="2" />
          <text x={ML + 88} y={MT + 11} className="fill-foreground" fontSize="9">Dew</text>
        </svg>
      );
    }

    if (chartType === 'Pxy' && 'bubbleCurve' in data) {
      const pxy = data as PxyData;
      const allP = [...pxy.bubbleCurve.map(p => p.P_Pa), ...pxy.dewCurve.map(p => p.P_Pa)];
      const Pmin = Math.min(...allP), Pmax = Math.max(...allP);
      const Prange = Pmax - Pmin || 1000;
      const toX = (x: number) => ML + x * pw;
      const toY = (P: number) => MT + ph - ((P - Pmin) / Prange) * ph;

      const bubblePath = pxy.bubbleCurve.map((p, i) => `${i === 0 ? 'M' : 'L'} ${toX(p.x1).toFixed(1)} ${toY(p.P_Pa).toFixed(1)}`).join(' ');
      const dewPath = pxy.dewCurve.map((p, i) => `${i === 0 ? 'M' : 'L'} ${toX(p.x1).toFixed(1)} ${toY(p.P_Pa).toFixed(1)}`).join(' ');

      return (
        <svg width={W} height={H} className="bg-card border border-border rounded">
          {[0, 0.25, 0.5, 0.75, 1].map(v => (
            <line key={`gx-${v}`} x1={toX(v)} y1={MT} x2={toX(v)} y2={MT + ph} stroke="var(--border)" strokeWidth="0.5" />
          ))}
          {Array.from({ length: 5 }, (_, i) => Pmin + (i / 4) * Prange).map(P => (
            <g key={`gy-${P}`}>
              <line x1={ML} y1={toY(P)} x2={ML + pw} y2={toY(P)} stroke="var(--border)" strokeWidth="0.5" />
              <text x={ML - 4} y={toY(P) + 3} textAnchor="end" className="fill-muted-foreground" fontSize="9">
                {fromSI_P(P, us.pressure).toFixed(2)}
              </text>
            </g>
          ))}
          <text x={ML + pw / 2} y={H - 4} textAnchor="middle" className="fill-foreground" fontSize="10">
            x, y ({pxy.comp1})
          </text>
          <text x={12} y={MT + ph / 2} textAnchor="middle" className="fill-foreground" fontSize="10"
            transform={`rotate(-90, 12, ${MT + ph / 2})`}>
            P ({pressLabel(us.pressure)})
          </text>
          <path d={bubblePath} fill="none" stroke="#3b82f6" strokeWidth="2" />
          <path d={dewPath} fill="none" stroke="#ef4444" strokeWidth="2" />
          <line x1={ML + 10} y1={MT + 8} x2={ML + 25} y2={MT + 8} stroke="#3b82f6" strokeWidth="2" />
          <text x={ML + 28} y={MT + 11} className="fill-foreground" fontSize="9">Bubble</text>
          <line x1={ML + 70} y1={MT + 8} x2={ML + 85} y2={MT + 8} stroke="#ef4444" strokeWidth="2" />
          <text x={ML + 88} y={MT + 11} className="fill-foreground" fontSize="9">Dew</text>
        </svg>
      );
    }

    if (chartType === 'xy' && 'points' in data) {
      const xy = data as XYData;
      const toX = (x: number) => ML + x * pw;
      const toY = (y: number) => MT + ph - y * ph;
      const curvePath = xy.points.map((p, i) => `${i === 0 ? 'M' : 'L'} ${toX(p.x1).toFixed(1)} ${toY(p.y1).toFixed(1)}`).join(' ');
      const diagPath = `M ${toX(0)} ${toY(0)} L ${toX(1)} ${toY(1)}`;

      return (
        <svg width={W} height={H} className="bg-card border border-border rounded">
          {[0, 0.25, 0.5, 0.75, 1].map(v => (
            <g key={`g-${v}`}>
              <line x1={toX(v)} y1={MT} x2={toX(v)} y2={MT + ph} stroke="var(--border)" strokeWidth="0.5" />
              <line x1={ML} y1={toY(v)} x2={ML + pw} y2={toY(v)} stroke="var(--border)" strokeWidth="0.5" />
            </g>
          ))}
          <text x={ML + pw / 2} y={H - 4} textAnchor="middle" className="fill-foreground" fontSize="10">
            x ({xy.comp1})
          </text>
          <text x={12} y={MT + ph / 2} textAnchor="middle" className="fill-foreground" fontSize="10"
            transform={`rotate(-90, 12, ${MT + ph / 2})`}>
            y ({xy.comp1})
          </text>
          <path d={diagPath} fill="none" stroke="var(--border)" strokeWidth="1" strokeDasharray="4 2" />
          <path d={curvePath} fill="none" stroke="#22c55e" strokeWidth="2" />
        </svg>
      );
    }

    return null;
  };

  return (
    <div className="p-3 space-y-3">
      <div className="flex items-end gap-3 flex-wrap">
        <label className="space-y-0.5 text-xs">
          <span className="text-muted-foreground">Component 1</span>
          <select value={comp1} onChange={e => setComp1(Number(e.target.value))}
            className="block w-32 px-1.5 py-1 border border-border rounded bg-background text-xs">
            {compounds.map((c, i) => <option key={i} value={i}>{c.displayName}</option>)}
          </select>
        </label>
        <label className="space-y-0.5 text-xs">
          <span className="text-muted-foreground">Component 2</span>
          <select value={comp2} onChange={e => setComp2(Number(e.target.value))}
            className="block w-32 px-1.5 py-1 border border-border rounded bg-background text-xs">
            {compounds.map((c, i) => <option key={i} value={i}>{c.displayName}</option>)}
          </select>
        </label>
        <label className="space-y-0.5 text-xs">
          <span className="text-muted-foreground">Chart Type</span>
          <select value={chartType} onChange={e => setChartType(e.target.value as any)}
            className="block w-24 px-1.5 py-1 border border-border rounded bg-background text-xs">
            <option value="Txy">T-xy</option>
            <option value="Pxy">P-xy</option>
            <option value="xy">x-y</option>
          </select>
        </label>
        {(chartType === 'Txy' || chartType === 'xy') && (
          <label className="space-y-0.5 text-xs">
            <span className="text-muted-foreground">P ({pressLabel(us.pressure)})</span>
            <input type="number" value={fromSI_P(pressure, us.pressure)} step="0.1"
              onChange={e => setPressure(toSI_P(parseFloat(e.target.value) || 101325, us.pressure))}
              className="block w-24 px-1.5 py-1 border border-border rounded bg-background text-xs font-mono" />
          </label>
        )}
        {chartType === 'Pxy' && (
          <label className="space-y-0.5 text-xs">
            <span className="text-muted-foreground">T ({tempLabel(us.temperature)})</span>
            <input type="number" value={fromSI_T(temperature, us.temperature)} step="1"
              onChange={e => setTemperature(toSI_T(parseFloat(e.target.value) || 350, us.temperature))}
              className="block w-24 px-1.5 py-1 border border-border rounded bg-background text-xs font-mono" />
          </label>
        )}
        <Button variant="outline" size="sm" onClick={generate} className="text-xs h-7">
          Generate
        </Button>
      </div>
      {data && (
        <div className="flex gap-4">
          {renderChart()}
          <div className="text-xs text-muted-foreground">
            <div className="font-semibold text-foreground mb-1">
              {chartType === 'Txy' ? 'T-xy Diagram' : chartType === 'Pxy' ? 'P-xy Diagram' : 'x-y Diagram'}
            </div>
            <div>{(data as any).comp1} / {(data as any).comp2}</div>
            {'P_Pa' in data && <div>P = {fromSI_P((data as any).P_Pa, us.pressure).toFixed(2)} {pressLabel(us.pressure)}</div>}
            {'T_K' in data && (data as PxyData).T_K && <div>T = {fromSI_T((data as PxyData).T_K, us.temperature).toFixed(1)} {tempLabel(us.temperature)}</div>}
            <div>Fluid Package: {fluidPackage}</div>
            <div className="mt-2">
              {'bubbleCurve' in data && <div>{(data as any).bubbleCurve.length} points</div>}
              {'points' in data && <div>{(data as XYData).points.length} points</div>}
            </div>
          </div>
        </div>
      )}
    </div>
  );
}

// ─── Stream Results Table ───────────────────────────────────────

function StreamResultsTable() {
  const { streams, compounds, unitSystem: us, openTab } = useCanopyStore();
  const sortedStreams = Object.values(streams).sort((a, b) => a.id.localeCompare(b.id));

  // Compute MW and mass flows for each stream
  const streamMW = sortedStreams.map(s => {
    let mw = 0;
    for (let i = 0; i < compounds.length; i++) mw += (s.moleFractions[i] ?? 0) * (compounds[i]?.molecularWeight ?? 0);
    return mw;
  });
  const streamMassFlow = sortedStreams.map((s, si) => s.totalFlow_molps * (streamMW[si] / 1000)); // kg/s

  // Mass fractions per stream
  const massFracs = sortedStreams.map((s, si) => {
    const mw = streamMW[si];
    if (mw === 0) return s.moleFractions.map(() => 0);
    return s.moleFractions.map((z, ci) => z * (compounds[ci]?.molecularWeight ?? 0) / mw);
  });
  const xFracs = sortedStreams.map(s => s.x_liquid ?? new Array(compounds.length).fill(0));
  const yFracs = sortedStreams.map(s => s.y_vapor ?? new Array(compounds.length).fill(0));

  return (
    <div className="p-3">
      <div className="flex items-center justify-between mb-2">
        <div>
          <p className="text-[11px] text-muted-foreground">Click any stream name below to open its detailed tab.</p>
        </div>
        <div className="flex items-center gap-1">
          <button
            onClick={() => {
              // Build TSV for clipboard
              const headers = ['Property', ...sortedStreams.map(s => s.name)];
              const rows: string[][] = [headers];
              rows.push([`T (${tempLabel(us.temperature)})`, ...sortedStreams.map(s => fmt(fromSI_T(s.T_K, us.temperature), 2))]);
              rows.push([`P (${pressLabel(us.pressure)})`, ...sortedStreams.map(s => fmt(fromSI_P(s.P_Pa, us.pressure), 3))]);
              rows.push([`Flow (${flowLabel(us.flow)})`, ...sortedStreams.map(s => fmt(fromSI_F(s.totalFlow_molps, us.flow), 3))]);
              rows.push(['Mass Flow (kg/s)', ...streamMassFlow.map(v => fmt(v, 4))]);
              rows.push(['MW (g/mol)', ...streamMW.map(v => fmt(v, 2))]);
              rows.push(['Phase', ...sortedStreams.map(s => s.phase)]);
              rows.push(['Vapor Frac', ...sortedStreams.map(s => fmt(s.vaporFraction, 4))]);
              rows.push(['H (J/mol)', ...sortedStreams.map(s => fmt(s.H_Jpmol, 1))]);
              for (let ci = 0; ci < compounds.length; ci++) {
                rows.push([`z_${compounds[ci].displayName}`, ...sortedStreams.map(s => fmt(s.moleFractions[ci], 4))]);
              }
              for (let ci = 0; ci < compounds.length; ci++) {
                rows.push([`w_${compounds[ci].displayName}`, ...massFracs.map(mf => fmt(mf[ci], 4))]);
              }
              const tsv = rows.map(r => r.join('\t')).join('\n');
              navigator.clipboard.writeText(tsv);
            }}
            className="px-2 py-0.5 text-xs bg-card border border-border rounded hover:bg-accent"
            title="Copy table to clipboard (TSV)"
          >
            📋 Copy
          </button>
          <button
            onClick={() => {
              // Build CSV download
              const headers = ['Property', ...sortedStreams.map(s => s.name)];
              const rows: string[][] = [headers];
              rows.push([`T (K)`, ...sortedStreams.map(s => s.T_K.toFixed(2))]);
              rows.push([`P (Pa)`, ...sortedStreams.map(s => s.P_Pa.toFixed(0))]);
              rows.push([`Flow (mol/s)`, ...sortedStreams.map(s => s.totalFlow_molps.toFixed(4))]);
              rows.push(['Phase', ...sortedStreams.map(s => s.phase)]);
              rows.push(['VapFrac', ...sortedStreams.map(s => s.vaporFraction.toFixed(4))]);
              rows.push(['H (J/mol)', ...sortedStreams.map(s => s.H_Jpmol.toFixed(2))]);
              for (let ci = 0; ci < compounds.length; ci++) {
                rows.push([`z_${compounds[ci].name}`, ...sortedStreams.map(s => (s.moleFractions[ci] ?? 0).toFixed(6))]);
              }
              const csv = rows.map(r => r.join(',')).join('\n');
              const blob = new Blob([csv], { type: 'text/csv' });
              const url = URL.createObjectURL(blob);
              const a = document.createElement('a');
              a.href = url; a.download = 'canopy-streams.csv'; a.click();
              URL.revokeObjectURL(url);
            }}
            className="px-2 py-0.5 text-xs bg-card border border-border rounded hover:bg-accent"
            title="Download as CSV"
          >
            ⬇ CSV
          </button>
        </div>
      </div>
      <div className="overflow-x-auto">
        <table className="text-xs w-full border-collapse">
          <thead>
            <tr className="border-b border-border">
              <th className="text-left p-1.5 font-semibold text-muted-foreground">Property</th>
              {sortedStreams.map(s => (
                <th key={s.id} className="text-right p-1.5 font-semibold min-w-[100px]">
                  <button
                    onClick={() => openTab({ type: 'stream', id: s.id })}
                    className="hover:text-blue-500 transition-colors underline decoration-dotted"
                  >
                    {s.name}
                  </button>
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            <tr className="border-b border-border/50">
              <td className="p-1.5 text-muted-foreground">Temperature ({tempLabel(us.temperature)})</td>
              {sortedStreams.map(s => (
                <td key={s.id} className="p-1.5 text-right">{fmt(fromSI_T(s.T_K, us.temperature), 1)}</td>
              ))}
            </tr>
            <tr className="border-b border-border/50">
              <td className="p-1.5 text-muted-foreground">Pressure ({pressLabel(us.pressure)})</td>
              {sortedStreams.map(s => (
                <td key={s.id} className="p-1.5 text-right">{fmt(fromSI_P(s.P_Pa, us.pressure), 2)}</td>
              ))}
            </tr>
            <tr className="border-b border-border/50">
              <td className="p-1.5 text-muted-foreground">Molar Flow ({flowLabel(us.flow)})</td>
              {sortedStreams.map(s => (
                <td key={s.id} className="p-1.5 text-right">{fmt(fromSI_F(s.totalFlow_molps, us.flow), 2)}</td>
              ))}
            </tr>
            <tr className="border-b border-border/50">
              <td className="p-1.5 text-muted-foreground">Mass Flow (kg/s)</td>
              {sortedStreams.map((_, si) => (
                <td key={sortedStreams[si].id} className="p-1.5 text-right">{fmt(streamMassFlow[si], 4)}</td>
              ))}
            </tr>
            <tr className="border-b border-border/50">
              <td className="p-1.5 text-muted-foreground">MW (g/mol)</td>
              {sortedStreams.map((_, si) => (
                <td key={sortedStreams[si].id} className="p-1.5 text-right">{fmt(streamMW[si], 2)}</td>
              ))}
            </tr>
            <tr className="border-b border-border/50">
              <td className="p-1.5 text-muted-foreground">Phase</td>
              {sortedStreams.map(s => <td key={s.id} className="p-1.5 text-right">{s.phase}</td>)}
            </tr>
            <tr className="border-b border-border/50">
              <td className="p-1.5 text-muted-foreground">Vapor Fraction</td>
              {sortedStreams.map(s => <td key={s.id} className="p-1.5 text-right">{fmt(s.vaporFraction, 4)}</td>)}
            </tr>
            <tr className="border-b border-border/50">
              <td className="p-1.5 text-muted-foreground">Enthalpy (J/mol)</td>
              {sortedStreams.map(s => <td key={s.id} className="p-1.5 text-right">{fmt(s.H_Jpmol, 1)}</td>)}
            </tr>
          </tbody>
        </table>
      </div>

      <Accordion type="multiple" className="mt-2 rounded border border-border/60 px-3">
        <AccordionItem value="z">
          <AccordionTrigger className="py-2 text-xs">Mole Fractions z</AccordionTrigger>
          <AccordionContent>
            <div className="overflow-x-auto">
              <table className="text-xs w-full border-collapse">
                <thead>
                  <tr className="border-b border-border">
                    <th className="text-left p-1.5 font-semibold text-muted-foreground">Component</th>
                    {sortedStreams.map(s => (
                      <th key={s.id} className="text-right p-1.5 font-semibold min-w-[100px]">{s.name}</th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {compounds.map((c, ci) => (
                    <tr key={`z-${c.name}`} className="border-b border-border/50">
                      <td className="p-1.5 text-muted-foreground">{c.displayName}</td>
                      {sortedStreams.map(s => (
                        <td key={s.id} className="p-1.5 text-right">{fmt(s.moleFractions[ci], 4)}</td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </AccordionContent>
        </AccordionItem>
        <AccordionItem value="w">
          <AccordionTrigger className="py-2 text-xs">Mass Fractions w</AccordionTrigger>
          <AccordionContent>
            <div className="overflow-x-auto">
              <table className="text-xs w-full border-collapse">
                <thead>
                  <tr className="border-b border-border">
                    <th className="text-left p-1.5 font-semibold text-muted-foreground">Component</th>
                    {sortedStreams.map(s => (
                      <th key={s.id} className="text-right p-1.5 font-semibold min-w-[100px]">{s.name}</th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {compounds.map((c, ci) => (
                    <tr key={`w-${c.name}`} className="border-b border-border/50">
                      <td className="p-1.5 text-muted-foreground">{c.displayName}</td>
                      {sortedStreams.map((_, si) => (
                        <td key={sortedStreams[si].id} className="p-1.5 text-right">{fmt(massFracs[si][ci], 4)}</td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </AccordionContent>
        </AccordionItem>
        <AccordionItem value="phase">
          <AccordionTrigger className="py-2 text-xs">Phase Compositions x / y</AccordionTrigger>
          <AccordionContent>
            <div className="overflow-x-auto">
              <table className="text-xs w-full border-collapse">
                <thead>
                  <tr className="border-b border-border">
                    <th className="text-left p-1.5 font-semibold text-muted-foreground">Component</th>
                    {sortedStreams.map(s => (
                      <th key={`${s.id}-x`} className="text-right p-1.5 font-semibold min-w-[140px]">{s.name} x / y</th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {compounds.map((c, ci) => (
                    <tr key={`phase-${c.name}`} className="border-b border-border/50">
                      <td className="p-1.5 text-muted-foreground">{c.displayName}</td>
                      {sortedStreams.map((_, si) => (
                        <td key={`${sortedStreams[si].id}-${c.name}`} className="p-1.5 text-right">
                          {fmt(xFracs[si][ci], 4)} / {fmt(yFracs[si][ci], 4)}
                        </td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </AccordionContent>
        </AccordionItem>
      </Accordion>
    </div>
  );
}
