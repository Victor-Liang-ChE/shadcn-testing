'use client';

import { memo } from 'react';
import { Handle, Position, type NodeProps } from '@xyflow/react';
import { fromSI_Q, energyLabel } from '@/lib/canopy/units';

// ─── Heater Node ──────────────────────────────────────────
export const HeaterNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div
      className={`relative flex flex-col items-center justify-center w-[100px] h-[80px] rounded-lg border-2 transition-colors
        ${selected ? 'border-blue-500 shadow-lg' : 'border-orange-400'}
        bg-orange-50 dark:bg-orange-950/40`}
    >
      <Handle type="target" position={Position.Left} id="in"
        className="!w-3 !h-3 !bg-blue-500 !border-2 !border-white" />
      <div className="text-xs font-bold text-orange-700 dark:text-orange-300">🔥</div>
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      {d.duty != null && d.unitSystem && (
        <div className="text-[9px] text-muted-foreground mt-0.5">
          Q = {fromSI_Q(d.duty, d.unitSystem.energy).toFixed(1)} {energyLabel(d.unitSystem.energy)}
        </div>
      )}
      <Handle type="source" position={Position.Right} id="out"
        className="!w-3 !h-3 !bg-blue-500 !border-2 !border-white" />
    </div>
  );
});
HeaterNode.displayName = 'HeaterNode';

// ─── Flash Drum Node ──────────────────────────────────────
export const FlashNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div
      className={`relative flex flex-col items-center justify-center w-[100px] h-[120px] rounded-xl border-2 transition-colors
        ${selected ? 'border-blue-500 shadow-lg' : 'border-purple-400'}
        bg-purple-50 dark:bg-purple-950/40`}
    >
      <Handle type="target" position={Position.Left} id="in"
        className="!w-3 !h-3 !bg-blue-500 !border-2 !border-white" />
      <div className="text-lg">⚗️</div>
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      {d.duty != null && d.unitSystem && (
        <div className="text-[9px] text-muted-foreground mt-0.5">
          Q = {fromSI_Q(d.duty, d.unitSystem.energy).toFixed(1)} {energyLabel(d.unitSystem.energy)}
        </div>
      )}
      <Handle type="source" position={Position.Top} id="vap"
        className="!w-3 !h-3 !bg-red-400 !border-2 !border-white" />
      <Handle type="source" position={Position.Bottom} id="liq"
        className="!w-3 !h-3 !bg-blue-400 !border-2 !border-white" />
    </div>
  );
});
FlashNode.displayName = 'FlashNode';

// ─── Mixer Node ───────────────────────────────────────────
export const MixerNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div
      className={`relative flex flex-col items-center justify-center w-[80px] h-[80px] rounded-lg border-2 transition-colors
        ${selected ? 'border-blue-500 shadow-lg' : 'border-green-400'}
        bg-green-50 dark:bg-green-950/40`}
    >
      <Handle type="target" position={Position.Left} id="in1"
        className="!w-3 !h-3 !bg-blue-500 !border-2 !border-white" style={{ top: '30%' }} />
      <Handle type="target" position={Position.Left} id="in2"
        className="!w-3 !h-3 !bg-blue-500 !border-2 !border-white" style={{ top: '70%' }} />
      <div className="text-xs font-bold">🔀</div>
      <div className="text-[10px] font-semibold text-center">{d.label}</div>
      <Handle type="source" position={Position.Right} id="out"
        className="!w-3 !h-3 !bg-blue-500 !border-2 !border-white" />
    </div>
  );
});
MixerNode.displayName = 'MixerNode';

// ─── Valve Node ────────────────────────────────────────────
export const ValveNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div
      className={`relative flex flex-col items-center justify-center w-[80px] h-[60px] border-2 transition-colors
        ${selected ? 'border-blue-500 shadow-lg' : 'border-gray-400'}
        bg-gray-50 dark:bg-gray-800/40`}
      style={{ clipPath: 'polygon(0 0, 100% 30%, 100% 70%, 0 100%)' }}
    >
      <Handle type="target" position={Position.Left} id="in"
        className="!w-3 !h-3 !bg-blue-500 !border-2 !border-white" />
      <div className="text-[10px] font-semibold">{d.label}</div>
      <Handle type="source" position={Position.Right} id="out"
        className="!w-3 !h-3 !bg-blue-500 !border-2 !border-white" />
    </div>
  );
});
ValveNode.displayName = 'ValveNode';

// ─── Splitter Node ─────────────────────────────────────────
export const SplitterNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div
      className={`relative flex flex-col items-center justify-center w-[80px] h-[80px] rounded-lg border-2 transition-colors
        ${selected ? 'border-blue-500 shadow-lg' : 'border-yellow-400'}
        bg-yellow-50 dark:bg-yellow-950/40`}
    >
      <Handle type="target" position={Position.Left} id="in"
        className="!w-3 !h-3 !bg-blue-500 !border-2 !border-white" />
      <div className="text-xs font-bold">↗️</div>
      <div className="text-[10px] font-semibold text-center">{d.label}</div>
      <Handle type="source" position={Position.Right} id="out1"
        className="!w-3 !h-3 !bg-blue-500 !border-2 !border-white" style={{ top: '30%' }} />
      <Handle type="source" position={Position.Right} id="out2"
        className="!w-3 !h-3 !bg-blue-500 !border-2 !border-white" style={{ top: '70%' }} />
    </div>
  );
});
SplitterNode.displayName = 'SplitterNode';

// ─── Feed/Product Arrow Node (stream label) ───────────────
export const StreamLabelNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  const isSource = d.isSource; // true = feed arrow (has outlet handle), false = product (has inlet handle)
  return (
    <div
      className={`relative flex flex-col items-center justify-center px-3 py-1 rounded border transition-colors
        ${selected ? 'border-blue-500' : 'border-dashed border-gray-400 dark:border-gray-600'}
        bg-card`}
    >
      {!isSource && (
        <Handle type="target" position={Position.Left} id="in"
          className="!w-3 !h-3 !bg-blue-500 !border-2 !border-white" />
      )}
      <div className="text-[10px] font-semibold">{d.label}</div>
      {d.streamInfo && (
        <div className="text-[8px] text-muted-foreground whitespace-pre">{d.streamInfo}</div>
      )}
      {isSource && (
        <Handle type="source" position={Position.Right} id="out"
          className="!w-3 !h-3 !bg-blue-500 !border-2 !border-white" />
      )}
    </div>
  );
});
StreamLabelNode.displayName = 'StreamLabelNode';

// Export node types map for React Flow
export const canopyNodeTypes = {
  heater: HeaterNode,
  flash: FlashNode,
  mixer: MixerNode,
  valve: ValveNode,
  splitter: SplitterNode,
  streamLabel: StreamLabelNode,
};
