'use client';

import { memo } from 'react';
import { Handle, Position, type NodeProps } from '@xyflow/react';
import { FileText } from 'lucide-react';
import { fromSI_Q, energyLabel } from '@/lib/canopy/units';
import {
  HeaterIcon, CoolerIcon, FlashDrumIcon, MixerIcon, SplitterIcon, ValveIcon,
  PumpIcon, CompressorIcon, HeatExchangerIcon, ColumnIcon, CSTRIcon, PFRIcon,
  PipeIcon, AbsorberIcon, SeparatorIcon, DecanterIcon, CrystallizerIcon,
  CrusherIcon, DryerIcon, ReactorIcon, RGibbsIcon, RYieldIcon, REquilIcon, RBatchIcon,
} from './UnitIcons';

const H = "!w-3 !h-3 !border-2 !border-white dark:!border-gray-900";

function DutyBadge({ duty, unitSystem }: { duty: number | null; unitSystem: any }) {
  if (duty == null || !unitSystem) return null;
  return (
    <div className="text-[9px] text-muted-foreground mt-0.5">
      Q = {fromSI_Q(duty, unitSystem.energy).toFixed(1)} {energyLabel(unitSystem.energy)}
    </div>
  );
}

function WorkBadge({ duty, unitSystem }: { duty: number | null; unitSystem: any }) {
  if (duty == null || !unitSystem) return null;
  return (
    <div className="text-[9px] text-muted-foreground mt-0.5">
      W = {fromSI_Q(duty, unitSystem.energy).toFixed(1)} {energyLabel(unitSystem.energy)}
    </div>
  );
}

// ─── Heater Node ──────────────────────────────────────────
export const HeaterNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[100px] h-[80px] rounded-lg border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-orange-400'} bg-orange-50 dark:bg-orange-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <HeaterIcon size={24} className="text-orange-600 dark:text-orange-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <DutyBadge duty={d.duty} unitSystem={d.unitSystem} />
      <Handle type="source" position={Position.Right} id="out" className={`${H} !bg-blue-500`} />
    </div>
  );
});
HeaterNode.displayName = 'HeaterNode';

// ─── Flash Drum Node ──────────────────────────────────────
export const FlashNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[100px] h-[120px] rounded-xl border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-purple-400'} bg-purple-50 dark:bg-purple-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <FlashDrumIcon size={22} className="text-purple-600 dark:text-purple-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <DutyBadge duty={d.duty} unitSystem={d.unitSystem} />
      <Handle type="source" position={Position.Top} id="vap" className={`${H} !bg-red-400`} />
      <Handle type="source" position={Position.Bottom} id="liq" className={`${H} !bg-blue-400`} />
    </div>
  );
});
FlashNode.displayName = 'FlashNode';

// ─── Mixer Node ───────────────────────────────────────────
export const MixerNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  const inletCount = Math.max(2, d.inletCount ?? 2);
  const inletHandles = [];
  for (let i = 1; i <= inletCount; i++) {
    const pct = ((i) / (inletCount + 1)) * 100;
    inletHandles.push(
      <Handle key={`in${i}`} type="target" position={Position.Left} id={`in${i}`}
        className={`${H} !bg-blue-500`} style={{ top: `${pct}%` }} />
    );
  }
  return (
    <div className={`relative flex flex-col items-center justify-center w-[80px] h-[80px] rounded-lg border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-green-400'} bg-green-50 dark:bg-green-950/40`}>
      {inletHandles}
      <MixerIcon size={22} className="text-green-600 dark:text-green-400" />
      <div className="text-[10px] font-semibold text-center">{d.label}</div>
      <Handle type="source" position={Position.Right} id="out" className={`${H} !bg-blue-500`} />
    </div>
  );
});
MixerNode.displayName = 'MixerNode';

// ─── Valve Node ────────────────────────────────────────────
export const ValveNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[80px] h-[60px] rounded-md border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-gray-400'} bg-gray-50 dark:bg-gray-800/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <ValveIcon size={20} className="text-gray-600 dark:text-gray-400" />
      <div className="text-[10px] font-semibold">{d.label}</div>
      <Handle type="source" position={Position.Right} id="out" className={`${H} !bg-blue-500`} />
    </div>
  );
});
ValveNode.displayName = 'ValveNode';

// ─── Splitter Node ─────────────────────────────────────────
export const SplitterNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  const outletCount = Math.max(2, d.outletCount ?? 2);
  const outletHandles = [];
  for (let i = 1; i <= outletCount; i++) {
    const pct = ((i) / (outletCount + 1)) * 100;
    outletHandles.push(
      <Handle key={`out${i}`} type="source" position={Position.Right} id={`out${i}`}
        className={`${H} !bg-blue-500`} style={{ top: `${pct}%` }} />
    );
  }
  return (
    <div className={`relative flex flex-col items-center justify-center w-[80px] h-[80px] rounded-lg border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-yellow-400'} bg-yellow-50 dark:bg-yellow-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <SplitterIcon size={22} className="text-yellow-600 dark:text-yellow-500" />
      <div className="text-[10px] font-semibold text-center">{d.label}</div>
      {outletHandles}
    </div>
  );
});
SplitterNode.displayName = 'SplitterNode';

// ─── Pump Node ─────────────────────────────────────────────
export const PumpNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[90px] h-[70px] rounded-lg border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-cyan-400'} bg-cyan-50 dark:bg-cyan-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <PumpIcon size={22} className="text-cyan-600 dark:text-cyan-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <WorkBadge duty={d.duty} unitSystem={d.unitSystem} />
      <Handle type="source" position={Position.Right} id="out" className={`${H} !bg-blue-500`} />
    </div>
  );
});
PumpNode.displayName = 'PumpNode';

// ─── Compressor Node ─────────────────────────────────────────
export const CompressorNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[100px] h-[80px] rounded-lg border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-sky-400'} bg-sky-50 dark:bg-sky-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <CompressorIcon size={22} className="text-sky-600 dark:text-sky-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <WorkBadge duty={d.duty} unitSystem={d.unitSystem} />
      <Handle type="source" position={Position.Right} id="out" className={`${H} !bg-blue-500`} />
    </div>
  );
});
CompressorNode.displayName = 'CompressorNode';

// ─── Heat Exchanger Node ─────────────────────────────────────
export const HeatXNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[100px] h-[100px] rounded-full border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-rose-400'} bg-rose-50 dark:bg-rose-950/40`}>
      <Handle type="target" position={Position.Left} id="hot-in" className={`${H} !bg-red-500`} style={{ top: '30%' }} />
      <Handle type="target" position={Position.Left} id="cold-in" className={`${H} !bg-blue-400`} style={{ top: '70%' }} />
      <HeatExchangerIcon size={24} className="text-rose-600 dark:text-rose-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <DutyBadge duty={d.duty} unitSystem={d.unitSystem} />
      <Handle type="source" position={Position.Right} id="hot-out" className={`${H} !bg-red-500`} style={{ top: '30%' }} />
      <Handle type="source" position={Position.Right} id="cold-out" className={`${H} !bg-blue-400`} style={{ top: '70%' }} />
    </div>
  );
});
HeatXNode.displayName = 'HeatXNode';

// ─── Pipe Node ────────────────────────────────────────────
export const PipeNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[120px] h-[40px] rounded border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-stone-400'} bg-stone-50 dark:bg-stone-800/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <PipeIcon size={18} className="text-stone-600 dark:text-stone-400" />
      <div className="text-[10px] font-semibold text-center">{d.label}</div>
      <Handle type="source" position={Position.Right} id="out" className={`${H} !bg-blue-500`} />
    </div>
  );
});
PipeNode.displayName = 'PipeNode';

// ─── RStoic (Stoichiometric Reactor) Node ─────────────────
export const RStoicNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[100px] h-[80px] rounded-lg border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-emerald-400'} bg-emerald-50 dark:bg-emerald-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <ReactorIcon size={22} className="text-emerald-600 dark:text-emerald-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <DutyBadge duty={d.duty} unitSystem={d.unitSystem} />
      <Handle type="source" position={Position.Right} id="out" className={`${H} !bg-blue-500`} />
    </div>
  );
});
RStoicNode.displayName = 'RStoicNode';

// ─── Column (Distillation) Node ───────────────────────────
export const ColumnNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  const hasSideVapor = (d.outletCount ?? 2) > 2;
  const hasSideLiquid = (d.outletCount ?? 2) > 3;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[80px] h-[140px] rounded-xl border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-indigo-400'} bg-indigo-50 dark:bg-indigo-950/40`}>
      <Handle type="target" position={Position.Left} id="feed" className={`${H} !bg-blue-500`} />
      <ColumnIcon size={18} className="text-indigo-600 dark:text-indigo-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <DutyBadge duty={d.duty != null ? Math.abs(d.duty) : null} unitSystem={d.unitSystem} />
      <Handle type="source" position={Position.Top} id="dist" className={`${H} !bg-red-400`} />
      {hasSideVapor && <Handle type="source" position={Position.Right} id="side-vap" className={`${H} !bg-red-300`} style={{ top: '34%' }} />}
      {hasSideLiquid && <Handle type="source" position={Position.Right} id="side-liq" className={`${H} !bg-blue-300`} style={{ top: '68%' }} />}
      <Handle type="source" position={Position.Bottom} id="bot" className={`${H} !bg-blue-400`} />
    </div>
  );
});
ColumnNode.displayName = 'ColumnNode';

// ─── Three-Phase Flash Node ───────────────────────────────
export const ThreePhaseFlashNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[100px] h-[130px] rounded-xl border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-violet-400'} bg-violet-50 dark:bg-violet-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <FlashDrumIcon size={20} className="text-violet-600 dark:text-violet-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <div className="text-[8px] text-muted-foreground">3-Phase</div>
      <Handle type="source" position={Position.Top} id="vap" className={`${H} !bg-red-400`} />
      <Handle type="source" position={Position.Right} id="liq1" className={`${H} !bg-blue-400`} />
      <Handle type="source" position={Position.Bottom} id="liq2" className={`${H} !bg-teal-400`} />
    </div>
  );
});
ThreePhaseFlashNode.displayName = 'ThreePhaseFlashNode';

// ─── CSTR Node ────────────────────────────────────────────
export const CSTRNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[100px] h-[80px] rounded-lg border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-emerald-400'} bg-emerald-50 dark:bg-emerald-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <CSTRIcon size={22} className="text-emerald-600 dark:text-emerald-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <DutyBadge duty={d.duty} unitSystem={d.unitSystem} />
      <Handle type="source" position={Position.Right} id="out" className={`${H} !bg-blue-500`} />
    </div>
  );
});
CSTRNode.displayName = 'CSTRNode';

// ─── PFR Node ─────────────────────────────────────────────
export const PFRNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[120px] h-[50px] rounded-full border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-emerald-400'} bg-emerald-50 dark:bg-emerald-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <PFRIcon size={18} className="text-emerald-600 dark:text-emerald-400" />
      <div className="text-[10px] font-semibold text-center">{d.label}</div>
      <Handle type="source" position={Position.Right} id="out" className={`${H} !bg-blue-500`} />
    </div>
  );
});
PFRNode.displayName = 'PFRNode';

// ─── Absorber Node ────────────────────────────────────────
export const AbsorberNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[80px] h-[140px] rounded-xl border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-teal-400'} bg-teal-50 dark:bg-teal-950/40`}>
      <Handle type="target" position={Position.Bottom} id="gas-in" className={`${H} !bg-red-400`} />
      <Handle type="target" position={Position.Top} id="liq-in" className={`${H} !bg-blue-400`} style={{ left: '25%' }} />
      <AbsorberIcon size={18} className="text-teal-600 dark:text-teal-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <Handle type="source" position={Position.Top} id="gas-out" className={`${H} !bg-red-400`} style={{ left: '75%' }} />
      <Handle type="source" position={Position.Bottom} id="liq-out" className={`${H} !bg-blue-400`} style={{ left: '75%' }} />
    </div>
  );
});
AbsorberNode.displayName = 'AbsorberNode';

export const ExtractorNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[90px] h-[140px] rounded-xl border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-amber-400'} bg-amber-50 dark:bg-amber-950/40`}>
      <Handle type="target" position={Position.Left} id="feed-in" className={`${H} !bg-blue-500`} />
      <Handle type="target" position={Position.Top} id="solvent-in" className={`${H} !bg-cyan-400`} style={{ left: '28%' }} />
      <AbsorberIcon size={18} className="text-amber-600 dark:text-amber-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <Handle type="source" position={Position.Top} id="extract-out" className={`${H} !bg-cyan-400`} style={{ left: '72%' }} />
      <Handle type="source" position={Position.Right} id="raffinate-out" className={`${H} !bg-blue-500`} />
    </div>
  );
});
ExtractorNode.displayName = 'ExtractorNode';

// ─── Component Separator Node ─────────────────────────────
export const ComponentSeparatorNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[100px] h-[90px] rounded-lg border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-amber-400'} bg-amber-50 dark:bg-amber-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <SeparatorIcon size={20} className="text-amber-600 dark:text-amber-500" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <Handle type="source" position={Position.Top} id="top" className={`${H} !bg-red-400`} />
      <Handle type="source" position={Position.Bottom} id="bot" className={`${H} !bg-blue-400`} />
    </div>
  );
});
ComponentSeparatorNode.displayName = 'ComponentSeparatorNode';

export const MembraneNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[100px] h-[100px] rounded-lg border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-fuchsia-400'} bg-fuchsia-50 dark:bg-fuchsia-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <SeparatorIcon size={20} className="text-fuchsia-600 dark:text-fuchsia-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <div className="text-[8px] text-muted-foreground">Selective</div>
      <Handle type="source" position={Position.Top} id="perm" className={`${H} !bg-cyan-400`} />
      <Handle type="source" position={Position.Right} id="ret" className={`${H} !bg-blue-500`} />
    </div>
  );
});
MembraneNode.displayName = 'MembraneNode';

export const CycloneNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[96px] h-[110px] rounded-lg border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-orange-400'} bg-orange-50 dark:bg-orange-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <SeparatorIcon size={18} className="text-orange-600 dark:text-orange-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <div className="text-[8px] text-muted-foreground">d50 cut</div>
      <Handle type="source" position={Position.Top} id="ovhd" className={`${H} !bg-blue-400`} />
      <Handle type="source" position={Position.Bottom} id="uflow" className={`${H} !bg-stone-500`} />
    </div>
  );
});
CycloneNode.displayName = 'CycloneNode';

export const FilterNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[100px] h-[96px] rounded-lg border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-slate-400'} bg-slate-50 dark:bg-slate-900/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <SeparatorIcon size={18} className="text-slate-600 dark:text-slate-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <div className="text-[8px] text-muted-foreground">Retained / clear</div>
      <Handle type="source" position={Position.Right} id="filtrate" className={`${H} !bg-blue-400`} />
      <Handle type="source" position={Position.Bottom} id="cake" className={`${H} !bg-stone-500`} />
    </div>
  );
});
FilterNode.displayName = 'FilterNode';

export const SteamDrumNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[100px] h-[120px] rounded-xl border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-sky-400'} bg-sky-50 dark:bg-sky-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <FlashDrumIcon size={20} className="text-sky-600 dark:text-sky-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <div className="text-[8px] text-muted-foreground">Steam drum</div>
      <DutyBadge duty={d.duty} unitSystem={d.unitSystem} />
      <Handle type="source" position={Position.Top} id="vap" className={`${H} !bg-red-400`} />
      <Handle type="source" position={Position.Bottom} id="liq" className={`${H} !bg-blue-400`} />
    </div>
  );
});
SteamDrumNode.displayName = 'SteamDrumNode';

export const AnalyzerNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[108px] h-[88px] rounded-lg border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-violet-400'} bg-violet-50 dark:bg-violet-950/40`}>
      <FileText size={18} className="text-violet-600 dark:text-violet-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-1">{d.label}</div>
      <div className="text-[8px] text-muted-foreground text-center px-2">
        {d.propertyLabel || d.propertyCode || 'Analyzer'}
      </div>
      {typeof d.signalValue === 'number' && (
        <div className="text-[9px] text-muted-foreground mt-0.5">{d.signalValue.toFixed(3)}</div>
      )}
    </div>
  );
});
AnalyzerNode.displayName = 'AnalyzerNode';

// ─── Decanter Node ────────────────────────────────────────
export const DecanterNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  return (
    <div className={`relative flex flex-col items-center justify-center w-[100px] h-[90px] rounded-lg border-2 transition-colors
      ${selected ? 'border-blue-500 shadow-lg' : 'border-sky-400'} bg-sky-50 dark:bg-sky-950/40`}>
      <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
      <DecanterIcon size={20} className="text-sky-600 dark:text-sky-400" />
      <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
      <Handle type="source" position={Position.Top} id="liq1" className={`${H} !bg-orange-400`} />
      <Handle type="source" position={Position.Bottom} id="liq2" className={`${H} !bg-blue-400`} />
    </div>
  );
});
DecanterNode.displayName = 'DecanterNode';

// ─── Generic Unit Node (for RYield, REquil, RGibbs, RBatch, Crystallizer, Crusher, Dryer) ──
const GENERIC_CFG: Record<string, {
  border: string; bg: string;
  Icon: React.ComponentType<{ size?: number; className?: string }>;
  showWork?: boolean;
}> = {
  ryield:       { border: 'border-lime-400', bg: 'bg-lime-50 dark:bg-lime-950/40', Icon: RYieldIcon },
  requil:       { border: 'border-emerald-400', bg: 'bg-emerald-50 dark:bg-emerald-950/40', Icon: REquilIcon },
  rgibbs:       { border: 'border-amber-400', bg: 'bg-amber-50 dark:bg-amber-950/40', Icon: RGibbsIcon },
  rbatch:       { border: 'border-emerald-400', bg: 'bg-emerald-50 dark:bg-emerald-950/40', Icon: RBatchIcon },
  crystallizer: { border: 'border-blue-400', bg: 'bg-blue-50 dark:bg-blue-950/40', Icon: CrystallizerIcon },
  crusher:      { border: 'border-stone-400', bg: 'bg-stone-50 dark:bg-stone-800/40', Icon: CrusherIcon, showWork: true },
  dryer:        { border: 'border-orange-400', bg: 'bg-orange-50 dark:bg-orange-950/40', Icon: DryerIcon },
  screen:       { border: 'border-zinc-400', bg: 'bg-zinc-50 dark:bg-zinc-900/40', Icon: SeparatorIcon },
  centrifuge:   { border: 'border-cyan-400', bg: 'bg-cyan-50 dark:bg-cyan-950/40', Icon: SeparatorIcon },
  steamheater:  { border: 'border-rose-400', bg: 'bg-rose-50 dark:bg-rose-950/40', Icon: HeaterIcon },
  steamturbine: { border: 'border-sky-400', bg: 'bg-sky-50 dark:bg-sky-950/40', Icon: CompressorIcon, showWork: true },
  steamvalve:   { border: 'border-slate-400', bg: 'bg-slate-50 dark:bg-slate-900/40', Icon: ValveIcon },
  steamheader:  { border: 'border-indigo-400', bg: 'bg-indigo-50 dark:bg-indigo-950/40', Icon: MixerIcon },
  pidcontroller:{ border: 'border-slate-400', bg: 'bg-slate-50 dark:bg-slate-900/40', Icon: ValveIcon },
  leadlag:      { border: 'border-violet-400', bg: 'bg-violet-50 dark:bg-violet-950/40', Icon: ValveIcon },
  deadtime:     { border: 'border-fuchsia-400', bg: 'bg-fuchsia-50 dark:bg-fuchsia-950/40', Icon: ValveIcon },
  signalselector:{ border: 'border-teal-400', bg: 'bg-teal-50 dark:bg-teal-950/40', Icon: ValveIcon },
  chargebalance:{ border: 'border-cyan-400', bg: 'bg-cyan-50 dark:bg-cyan-950/40', Icon: SeparatorIcon },
  electrolyteequilibrium:{ border: 'border-emerald-400', bg: 'bg-emerald-50 dark:bg-emerald-950/40', Icon: REquilIcon },
};

function makeGenericNode(nodeType: string) {
  const cfg = GENERIC_CFG[nodeType];
  const Node = memo(({ data, selected }: NodeProps) => {
    const d = data as any;
    const { border, bg, Icon, showWork } = cfg ?? {
      border: 'border-gray-400', bg: 'bg-gray-50 dark:bg-gray-800/40',
      Icon: ReactorIcon, showWork: false,
    };
    return (
      <div className={`relative flex flex-col items-center justify-center w-[100px] h-[80px] rounded-lg border-2 transition-colors
        ${selected ? 'border-blue-500 shadow-lg' : border} ${bg}`}>
        <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />
        <Icon size={22} className="text-current opacity-70" />
        <div className="text-[10px] font-semibold text-center leading-tight mt-0.5">{d.label}</div>
        {showWork ? <WorkBadge duty={d.duty} unitSystem={d.unitSystem} /> : <DutyBadge duty={d.duty} unitSystem={d.unitSystem} />}
        <Handle type="source" position={Position.Right} id="out" className={`${H} !bg-blue-500`} />
      </div>
    );
  });
  Node.displayName = `${nodeType}Node`;
  return Node;
}

// ─── Feed/Product/Intermediate Arrow Node (stream label) ───────────────
export const StreamLabelNode = memo(({ data, selected }: NodeProps) => {
  const d = data as any;
  const isSource = d.isSource;
  const isIntermediate = d.isIntermediate;
  return (
    <div className={`relative flex flex-col items-center justify-center px-3 py-1 rounded border transition-colors
      ${selected ? 'border-blue-500' : 'border-dashed border-gray-400 dark:border-gray-600'} bg-card`}>
      {(!isSource || isIntermediate) && <Handle type="target" position={Position.Left} id="in" className={`${H} !bg-blue-500`} />}
      <div className="text-[10px] font-semibold">{d.label}</div>
      {d.streamInfo && <div className="text-[8px] text-muted-foreground whitespace-pre">{d.streamInfo}</div>}
      {(isSource || isIntermediate) && <Handle type="source" position={Position.Right} id="out" className={`${H} !bg-blue-500`} />}
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
  pump: PumpNode,
  compressor: CompressorNode,
  heatx: HeatXNode,
  pipe: PipeNode,
  rstoic: RStoicNode,
  column: ColumnNode,
  radfrac: ColumnNode,
  ratefrac: ColumnNode,
  threephaseflash: ThreePhaseFlashNode,
  cstr: CSTRNode,
  pfr: PFRNode,
  absorber: AbsorberNode,
  extractor: ExtractorNode,
  componentseparator: ComponentSeparatorNode,
  decanter: DecanterNode,
  membrane: MembraneNode,
  cyclone: CycloneNode,
  filter: FilterNode,
  screen: makeGenericNode('screen'),
  centrifuge: makeGenericNode('centrifuge'),
  steamdrum: SteamDrumNode,
  ryield: makeGenericNode('ryield'),
  requil: makeGenericNode('requil'),
  rgibbs: makeGenericNode('rgibbs'),
  mheatx: HeatXNode,
  rbatch: makeGenericNode('rbatch'),
  crystallizer: makeGenericNode('crystallizer'),
  crusher: makeGenericNode('crusher'),
  dryer: makeGenericNode('dryer'),
  steamheater: makeGenericNode('steamheater'),
  steamturbine: makeGenericNode('steamturbine'),
  steamvalve: makeGenericNode('steamvalve'),
  steamheader: makeGenericNode('steamheader'),
  steamtrap: makeGenericNode('steamtrap'),
  pidcontroller: makeGenericNode('pidcontroller'),
  leadlag: makeGenericNode('leadlag'),
  deadtime: makeGenericNode('deadtime'),
  signalselector: makeGenericNode('signalselector'),
  analyzer: AnalyzerNode,
  chargebalance: makeGenericNode('chargebalance'),
  electrolyteequilibrium: makeGenericNode('electrolyteequilibrium'),
  streamLabel: StreamLabelNode,
};
