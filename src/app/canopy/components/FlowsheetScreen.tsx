'use client';

import { useCallback, useMemo, useState, useRef, useEffect } from 'react';
import {
  ReactFlow,
  Background,
  Controls,
  MiniMap,
  type Edge,
  type Node,
  type NodeChange,
  MarkerType,
  Panel,
  applyNodeChanges,
  useReactFlow,
  ReactFlowProvider,
} from '@xyflow/react';
import '@xyflow/react/dist/style.css';
import { useCanopyStore, type TabItem } from '@/lib/canopy/store';
import { canopyNodeTypes } from './FlowNodes';
import { X, ArrowRight, Trash2 } from 'lucide-react';
import { Button } from '@/components/ui/button';
import {
  fromSI_T, toSI_T, fromSI_P, toSI_P, fromSI_F, toSI_F, fromSI_Q,
  tempLabel, pressLabel, flowLabel, energyLabel,
} from '@/lib/canopy/units';
import type { UnitOpType, MaterialStream } from '@/lib/canopy/types';

function fmt(val: number, decimals = 2): string {
  if (!isFinite(val)) return '—';
  return val.toFixed(decimals);
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

const PALETTE_ITEMS: { type: UnitOpType; label: string; icon: string }[] = [
  { type: 'Heater', label: 'Heater', icon: '🔥' },
  { type: 'Flash', label: 'Flash Drum', icon: '⚗️' },
  { type: 'Mixer', label: 'Mixer', icon: '🔀' },
  { type: 'Valve', label: 'Valve', icon: '🔧' },
  { type: 'Splitter', label: 'Splitter', icon: '↗️' },
];

function ObjectPalette() {
  const addUnitFromPalette = useCanopyStore(s => s.addUnitFromPalette);
  const addFeedStream = useCanopyStore(s => s.addFeedStream);

  const onDragStart = (e: React.DragEvent, type: UnitOpType | 'feed') => {
    e.dataTransfer.setData('canopy/object-type', type);
    e.dataTransfer.effectAllowed = 'copy';
  };

  return (
    <div className="flex flex-col gap-1 bg-card/95 backdrop-blur border border-border rounded-lg p-2 shadow-lg w-36">
      <div className="text-[10px] font-semibold text-muted-foreground uppercase tracking-wider px-1 mb-1">
        Objects — drag or click
      </div>
      {PALETTE_ITEMS.map(item => (
        <div
          key={item.type}
          draggable
          onDragStart={(e) => onDragStart(e, item.type)}
          onClick={() => addUnitFromPalette(item.type, { x: 350, y: 200 })}
          className="flex items-center gap-2 px-2 py-1.5 text-xs rounded hover:bg-accent transition-colors cursor-grab active:cursor-grabbing select-none"
          title={`Drag or click to add ${item.label}`}
        >
          <span>{item.icon}</span>
          <span>{item.label}</span>
        </div>
      ))}
      <hr className="border-border my-1" />
      <div
        draggable
        onDragStart={(e) => onDragStart(e, 'feed')}
        onClick={addFeedStream}
        className="flex items-center gap-2 px-2 py-1.5 text-xs rounded hover:bg-accent transition-colors cursor-grab active:cursor-grabbing select-none"
        title="Drag or click to add Feed Stream"
      >
        <ArrowRight className="w-3 h-3 text-blue-500" />
        <span>Feed Stream</span>
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
    nodePositions, setNodePosition, addUnitFromPalette, addFeedStream, removeUnit,
  } = useCanopyStore();

  const { screenToFlowPosition } = useReactFlow();

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
      result.push({
        id: u.id,
        type: nodeType === 'flash' ? 'flash' : nodeType,
        position: pos,
        data: {
          label: u.name,
          duty: u.solved ? u.duty_W : null,
          unitSystem,
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
    return result;
  }, [streams, units, compounds, solved, nodePositions, unitSystem]);

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
              targetHandle: port.id.split('-').pop() || 'in',
              animated: true,
              style: { stroke: '#3b82f6', strokeWidth: 2 },
              markerEnd: { type: MarkerType.ArrowClosed, color: '#3b82f6' },
              label: stream.name,
              labelStyle: { fontSize: 10 },
            });
          } else if (stream?.sourceUnitId) {
            const srcUnit = units[stream.sourceUnitId];
            const srcPort = srcUnit?.ports.find(p => p.streamId === port.streamId && p.type === 'outlet');
            if (srcUnit && srcPort) {
              result.push({
                id: `e-${port.streamId}`,
                source: srcUnit.id,
                sourceHandle: srcPort.id.split('-').pop() || 'out',
                target: u.id,
                targetHandle: port.id.split('-').pop() || 'in',
                animated: solved,
                style: { stroke: '#3b82f6', strokeWidth: 2 },
                markerEnd: { type: MarkerType.ArrowClosed, color: '#3b82f6' },
                label: stream.name,
                labelStyle: { fontSize: 10 },
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
            const handleId = port.id.split('-').pop() || 'out';
            result.push({
              id: `e-prod-${port.streamId}`,
              source: u.id,
              sourceHandle: handleId,
              target: `label-${port.streamId}`,
              targetHandle: 'in',
              animated: solved,
              style: {
                stroke: handleId === 'vap' ? '#ef4444' : '#3b82f6',
                strokeWidth: 2,
              },
              markerEnd: {
                type: MarkerType.ArrowClosed,
                color: handleId === 'vap' ? '#ef4444' : '#3b82f6',
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
      addFeedStream();
      // The store adds it but doesn't know position yet — we'll let it default and user can drag
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

  // ── Delete selected unit on Delete/Backspace key ──
  const onKeyDown = useCallback((e: React.KeyboardEvent) => {
    if (e.key === 'Delete' || e.key === 'Backspace') {
      const selected = rfNodes.find(n => n.selected && units[n.id]);
      if (selected) {
        removeUnit(selected.id);
      }
    }
  }, [rfNodes, units, removeUnit]);

  // ── Selected unit for toolbar delete button ──
  const selectedUnitId = rfNodes.find(n => n.selected && units[n.id])?.id;

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
          <ReactFlow
            nodes={rfNodes}
            edges={edges}
            nodeTypes={canopyNodeTypes}
            onNodeDoubleClick={onNodeDoubleClick}
            onNodesChange={onNodesChange}
            fitView
            proOptions={{ hideAttribution: true }}
            minZoom={0.3}
            maxZoom={2}
            nodesDraggable={true}
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
                <div className="text-xs text-muted-foreground bg-card/80 px-2 py-1 rounded border border-border">
                  Double-click to edit · Select + Delete to remove
                </div>
              </div>
            </Panel>
          </ReactFlow>
        </div>
      ) : activeTab.type === 'unit' ? (
        <UnitDetailPanel unitId={activeTab.id} />
      ) : (
        <StreamDetailPanel streamId={activeTab.id} />
      )}

      {/* Stream Results Table */}
      {solved && activeTab === null && (
        <div className="border-t border-border bg-card max-h-[35vh] overflow-auto">
          <StreamResultsTable />
        </div>
      )}
    </div>
  );
}

// ─── Stream Detail / Editor Panel ───────────────────────────────

function StreamDetailPanel({ streamId }: { streamId: string }) {
  const { streams, compounds, updateStream, units, unitSystem: us } = useCanopyStore();
  const stream = streams[streamId];
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

  return (
    <div className="flex-1 overflow-auto p-6 space-y-6">
      <div className="flex items-center gap-3">
        <h2 className="text-lg font-bold">{stream.name}</h2>
        <span className="text-xs px-2 py-0.5 rounded bg-muted text-muted-foreground">
          {stream.id}
        </span>
        {isUpstreamDetermined && (
          <span className="text-xs px-2 py-0.5 rounded bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300">
            🔒 Set by {sourceUnit?.name ?? 'upstream unit'}
          </span>
        )}
        {stream.solved && (
          <span className="text-xs px-2 py-0.5 rounded bg-green-100 text-green-800 dark:bg-green-900/40 dark:text-green-300">
            ✓ Solved
          </span>
        )}
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
      </div>

      {/* Results (if solved) */}
      {stream.solved && (
        <div className="space-y-3 border-t border-border pt-4">
          <h3 className="text-sm font-semibold">Flash Results</h3>
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

// ─── Unit Detail Panel ──────────────────────────────────────

function UnitDetailPanel({ unitId }: { unitId: string }) {
  const { units, streams, compounds, updateUnit, openTab, unitSystem: us } = useCanopyStore();
  const unit = units[unitId];
  if (!unit) return <div className="p-4 text-muted-foreground">Unit not found.</div>;

  const paramLabels: Record<string, React.ReactNode> = {
    targetT_K: <>Target T ({tempLabel(us.temperature)})</>,
    outletP_Pa: <>Outlet P ({pressLabel(us.pressure)})</>,
    T_K: <>T ({tempLabel(us.temperature)})</>,
    P_Pa: <>P ({pressLabel(us.pressure)})</>,
    splitRatio: <>Split Ratio</>,
  };

  const displayValue = (key: string, val: number): number => {
    if (key === 'targetT_K' || key === 'T_K') return fromSI_T(val, us.temperature);
    if (key === 'outletP_Pa' || key === 'P_Pa') return fromSI_P(val, us.pressure);
    return val;
  };

  const fromDisplay = (key: string, display: number): number => {
    if (key === 'targetT_K' || key === 'T_K') return toSI_T(display, us.temperature);
    if (key === 'outletP_Pa' || key === 'P_Pa') return toSI_P(display, us.pressure);
    return display;
  };

  return (
    <div className="flex-1 overflow-auto p-6 space-y-4">
      <div className="flex items-center gap-3">
        <h2 className="text-lg font-bold">{unit.name}</h2>
        <span className="text-xs px-2 py-0.5 rounded bg-muted text-muted-foreground">{unit.type}</span>
        {unit.solved && (
          <span className="text-xs px-2 py-0.5 rounded bg-green-100 text-green-800 dark:bg-green-900/40 dark:text-green-300">
            ✓ Solved
          </span>
        )}
      </div>

      {/* Parameters */}
      <div className="space-y-2">
        <h3 className="text-sm font-semibold">Parameters</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
          {Object.entries(unit.params).map(([key, value]) => (
            <FieldInput
              key={key}
              label={paramLabels[key] ?? key}
              value={displayValue(key, typeof value === 'number' ? value : 0)}
              onChange={(v) => {
                updateUnit(unitId, {
                  params: { ...unit.params, [key]: fromDisplay(key, v) },
                });
              }}
              locked={false}
              decimals={key === 'splitRatio' ? 3 : 2}
            />
          ))}
        </div>
      </div>

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
        <div className="text-sm">
          <span className="font-semibold">Heat Duty:</span>{' '}
          <span className={unit.duty_W >= 0 ? 'text-red-500' : 'text-blue-500'}>
            {fmt(fromSI_Q(unit.duty_W, us.energy), 2)} {energyLabel(us.energy)}
          </span>
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

// ─── Stream Results Table ───────────────────────────────────────

function StreamResultsTable() {
  const { streams, compounds, unitSystem: us, openTab } = useCanopyStore();
  const sortedStreams = Object.values(streams).sort((a, b) => a.id.localeCompare(b.id));

  return (
    <div className="p-3">
      <h3 className="text-sm font-bold mb-2">Stream Results</h3>
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
            {compounds.map((c, ci) => (
              <tr key={c.name} className="border-b border-border/50">
                <td className="p-1.5 text-muted-foreground">z<sub>{c.displayName}</sub></td>
                {sortedStreams.map(s => (
                  <td key={s.id} className="p-1.5 text-right">{fmt(s.moleFractions[ci], 4)}</td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
