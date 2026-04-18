'use client';

// Process engineering SVG icons for unit operations
// These replace the emoji-based icons with proper P&ID-style symbols

import { memo } from 'react';

interface IconProps {
  size?: number;
  className?: string;
}

/** Heater / Cooler — circle with flame or snowflake */
export const HeaterIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <circle cx="24" cy="24" r="20" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* flame shape */}
    <path d="M24 12c0 0-8 8-8 16a8 8 0 0016 0c0-8-8-16-8-16z" fill="currentColor" fillOpacity="0.15" stroke="currentColor" strokeWidth="1.5" />
    <path d="M24 22c0 0-3 3-3 6a3 3 0 006 0c0-3-3-6-3-6z" fill="currentColor" fillOpacity="0.3" />
  </svg>
));
HeaterIcon.displayName = 'HeaterIcon';

export const CoolerIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <circle cx="24" cy="24" r="20" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* snowflake */}
    <line x1="24" y1="12" x2="24" y2="36" stroke="currentColor" strokeWidth="2" />
    <line x1="12" y1="24" x2="36" y2="24" stroke="currentColor" strokeWidth="2" />
    <line x1="15.5" y1="15.5" x2="32.5" y2="32.5" stroke="currentColor" strokeWidth="1.5" />
    <line x1="32.5" y1="15.5" x2="15.5" y2="32.5" stroke="currentColor" strokeWidth="1.5" />
    {/* crystal tips */}
    <line x1="24" y1="12" x2="21" y2="15" stroke="currentColor" strokeWidth="1.5" />
    <line x1="24" y1="12" x2="27" y2="15" stroke="currentColor" strokeWidth="1.5" />
    <line x1="24" y1="36" x2="21" y2="33" stroke="currentColor" strokeWidth="1.5" />
    <line x1="24" y1="36" x2="27" y2="33" stroke="currentColor" strokeWidth="1.5" />
  </svg>
));
CoolerIcon.displayName = 'CoolerIcon';

/** Flash Drum — vertical cylinder with level line */
export const FlashDrumIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 64" fill="none" className={className}>
    {/* vessel body */}
    <rect x="10" y="8" width="28" height="48" rx="6" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* liquid level (wavy) */}
    <path d="M10 38 Q17 34, 24 38 Q31 42, 38 38" stroke="currentColor" strokeWidth="1.5" strokeDasharray="3 2" />
    {/* liquid fill */}
    <path d="M10 38 Q17 34, 24 38 Q31 42, 38 38 L38 50 Q38 56 32 56 L16 56 Q10 56 10 50 Z" fill="currentColor" fillOpacity="0.1" />
    {/* vapor arrow up */}
    <line x1="24" y1="20" x2="24" y2="14" stroke="currentColor" strokeWidth="1.5" />
    <polyline points="21,17 24,14 27,17" stroke="currentColor" strokeWidth="1.5" fill="none" />
  </svg>
));
FlashDrumIcon.displayName = 'FlashDrumIcon';

/** Mixer — converging triangle/funnel */
export const MixerIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <polygon points="6,8 42,24 6,40" stroke="currentColor" strokeWidth="2.5" fill="currentColor" fillOpacity="0.08" strokeLinejoin="round" />
  </svg>
));
MixerIcon.displayName = 'MixerIcon';

/** Splitter — diverging triangle */
export const SplitterIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <polygon points="6,24 42,8 42,40" stroke="currentColor" strokeWidth="2.5" fill="currentColor" fillOpacity="0.08" strokeLinejoin="round" />
  </svg>
));
SplitterIcon.displayName = 'SplitterIcon';

/** Valve — bowtie shape */
export const ValveIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <polygon points="4,8 24,24 4,40" stroke="currentColor" strokeWidth="2.5" fill="currentColor" fillOpacity="0.08" strokeLinejoin="round" />
    <polygon points="44,8 24,24 44,40" stroke="currentColor" strokeWidth="2.5" fill="currentColor" fillOpacity="0.08" strokeLinejoin="round" />
    {/* stem */}
    <line x1="24" y1="24" x2="24" y2="6" stroke="currentColor" strokeWidth="2" />
    <line x1="18" y1="6" x2="30" y2="6" stroke="currentColor" strokeWidth="2.5" />
  </svg>
));
ValveIcon.displayName = 'ValveIcon';

/** Pump — circle with discharge triangle */
export const PumpIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <circle cx="22" cy="26" r="16" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* discharge nozzle */}
    <line x1="38" y1="18" x2="46" y2="12" stroke="currentColor" strokeWidth="2.5" />
    {/* arrow inside showing rotation */}
    <path d="M14 26 L22 18 L30 26" stroke="currentColor" strokeWidth="2" fill="none" />
  </svg>
));
PumpIcon.displayName = 'PumpIcon';

/** Compressor — trapezoid with reduced outlet */
export const CompressorIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <polygon points="6,6 42,16 42,32 6,42" stroke="currentColor" strokeWidth="2.5" fill="currentColor" fillOpacity="0.08" strokeLinejoin="round" />
    {/* rotation arrow */}
    <path d="M18 20 A6 6 0 1 1 18 28" stroke="currentColor" strokeWidth="1.5" fill="none" />
    <polyline points="15,27 18,28 18,25" stroke="currentColor" strokeWidth="1.5" fill="none" />
  </svg>
));
CompressorIcon.displayName = 'CompressorIcon';

/** Heat Exchanger — circle with S-curve (shell & tube) */
export const HeatExchangerIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <circle cx="24" cy="24" r="20" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* S-curve tubes */}
    <path d="M8 18 Q16 18, 24 24 Q32 30, 40 30" stroke="#ef4444" strokeWidth="2" fill="none" />
    <path d="M8 30 Q16 30, 24 24 Q32 18, 40 18" stroke="#3b82f6" strokeWidth="2" fill="none" />
  </svg>
));
HeatExchangerIcon.displayName = 'HeatExchangerIcon';

/** Column / Distillation — tall vertical vessel with trays */
export const ColumnIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size * 1.8} viewBox="0 0 40 72" fill="none" className={className}>
    {/* vessel */}
    <rect x="6" y="4" width="28" height="64" rx="5" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* trays */}
    {[16, 24, 32, 40, 48, 56].map(y => (
      <line key={y} x1="10" y1={y} x2="30" y2={y} stroke="currentColor" strokeWidth="1" strokeOpacity="0.5" />
    ))}
    {/* condenser (top) */}
    <circle cx="20" cy="4" r="4" stroke="currentColor" strokeWidth="1.5" fill="currentColor" fillOpacity="0.1" />
    {/* reboiler (bottom) */}
    <circle cx="20" cy="68" r="4" stroke="currentColor" strokeWidth="1.5" fill="currentColor" fillOpacity="0.1" />
  </svg>
));
ColumnIcon.displayName = 'ColumnIcon';

/** Reactor (CSTR) — cylinder with agitator */
export const CSTRIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <rect x="8" y="6" width="32" height="36" rx="4" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* agitator shaft */}
    <line x1="24" y1="2" x2="24" y2="28" stroke="currentColor" strokeWidth="2" />
    {/* impeller blades */}
    <line x1="16" y1="28" x2="32" y2="28" stroke="currentColor" strokeWidth="2.5" />
    {/* motor */}
    <rect x="20" y="0" width="8" height="4" rx="1" stroke="currentColor" strokeWidth="1.5" fill="currentColor" fillOpacity="0.2" />
  </svg>
));
CSTRIcon.displayName = 'CSTRIcon';

/** PFR — horizontal tube with arrow */
export const PFRIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size * 1.5} height={size} viewBox="0 0 72 36" fill="none" className={className}>
    <rect x="4" y="6" width="64" height="24" rx="12" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* flow arrow */}
    <path d="M20 18 L48 18" stroke="currentColor" strokeWidth="1.5" />
    <polyline points="44,14 48,18 44,22" stroke="currentColor" strokeWidth="1.5" fill="none" />
  </svg>
));
PFRIcon.displayName = 'PFRIcon';

/** Pipe — simple horizontal line with flanges */
export const PipeIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size * 1.5} height={size * 0.5} viewBox="0 0 72 18" fill="none" className={className}>
    <line x1="4" y1="9" x2="68" y2="9" stroke="currentColor" strokeWidth="3" />
    <line x1="4" y1="3" x2="4" y2="15" stroke="currentColor" strokeWidth="2" />
    <line x1="68" y1="3" x2="68" y2="15" stroke="currentColor" strokeWidth="2" />
  </svg>
));
PipeIcon.displayName = 'PipeIcon';

/** Absorber — tall packed column */
export const AbsorberIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size * 1.6} viewBox="0 0 40 64" fill="none" className={className}>
    <rect x="6" y="4" width="28" height="56" rx="5" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* packing representation */}
    {[14, 22, 30, 38, 46].map(y => (
      <g key={y}>
        <line x1="10" y1={y} x2="18" y2={y + 4} stroke="currentColor" strokeWidth="1" strokeOpacity="0.4" />
        <line x1="22" y1={y} x2="30" y2={y + 4} stroke="currentColor" strokeWidth="1" strokeOpacity="0.4" />
        <line x1="14" y1={y + 4} x2="26" y2={y} stroke="currentColor" strokeWidth="1" strokeOpacity="0.4" />
      </g>
    ))}
  </svg>
));
AbsorberIcon.displayName = 'AbsorberIcon';

/** Component Separator — box with split arrow */
export const SeparatorIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <rect x="6" y="6" width="36" height="36" rx="4" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* dividing line */}
    <line x1="6" y1="24" x2="42" y2="24" stroke="currentColor" strokeWidth="1.5" strokeDasharray="4 2" />
    {/* up arrow */}
    <path d="M24 20 L24 10" stroke="currentColor" strokeWidth="1.5" />
    <polyline points="21,13 24,10 27,13" stroke="currentColor" strokeWidth="1.5" fill="none" />
    {/* down arrow */}
    <path d="M24 28 L24 38" stroke="currentColor" strokeWidth="1.5" />
    <polyline points="21,35 24,38 27,35" stroke="currentColor" strokeWidth="1.5" fill="none" />
  </svg>
));
SeparatorIcon.displayName = 'SeparatorIcon';

/** Decanter — horizontal vessel with phase line */
export const DecanterIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size * 1.3} height={size} viewBox="0 0 56 40" fill="none" className={className}>
    <rect x="4" y="4" width="48" height="32" rx="8" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* phase boundary */}
    <line x1="4" y1="20" x2="52" y2="20" stroke="currentColor" strokeWidth="1.5" strokeDasharray="4 2" />
    {/* light phase label area */}
    <text x="28" y="14" textAnchor="middle" fontSize="7" fill="currentColor" fillOpacity="0.5">L</text>
    {/* heavy phase label area */}
    <text x="28" y="30" textAnchor="middle" fontSize="7" fill="currentColor" fillOpacity="0.5">H</text>
  </svg>
));
DecanterIcon.displayName = 'DecanterIcon';

/** Crystallizer — vessel with crystal shapes */
export const CrystallizerIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <rect x="8" y="8" width="32" height="32" rx="4" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* agitator */}
    <line x1="24" y1="4" x2="24" y2="20" stroke="currentColor" strokeWidth="1.5" />
    <line x1="18" y1="20" x2="30" y2="20" stroke="currentColor" strokeWidth="2" />
    {/* crystal shapes */}
    <polygon points="16,30 20,26 24,30 20,34" stroke="currentColor" strokeWidth="1" fill="currentColor" fillOpacity="0.2" />
    <polygon points="28,32 32,28 36,32 32,36" stroke="currentColor" strokeWidth="1" fill="currentColor" fillOpacity="0.2" />
  </svg>
));
CrystallizerIcon.displayName = 'CrystallizerIcon';

/** Crusher — jaw crusher profile */
export const CrusherIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    {/* crusher jaws */}
    <polygon points="8,4 24,24 8,44" stroke="currentColor" strokeWidth="2.5" fill="currentColor" fillOpacity="0.08" />
    <polygon points="40,4 24,24 40,44" stroke="currentColor" strokeWidth="2.5" fill="currentColor" fillOpacity="0.08" />
    {/* discharge */}
    <line x1="24" y1="28" x2="24" y2="46" stroke="currentColor" strokeWidth="2" />
    <polyline points="20,42 24,46 28,42" stroke="currentColor" strokeWidth="1.5" fill="none" />
  </svg>
));
CrusherIcon.displayName = 'CrusherIcon';

/** Dryer — horizontal cylinder with heat wavy lines */
export const DryerIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size * 1.3} height={size} viewBox="0 0 56 40" fill="none" className={className}>
    <rect x="4" y="6" width="48" height="28" rx="6" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* heat lines (wavy) */}
    <path d="M14 16 Q17 13, 20 16 Q23 19, 26 16" stroke="#ef4444" strokeWidth="1.5" fill="none" />
    <path d="M30 16 Q33 13, 36 16 Q39 19, 42 16" stroke="#ef4444" strokeWidth="1.5" fill="none" />
    {/* moisture drops */}
    <circle cx="18" cy="28" r="2" fill="#3b82f6" fillOpacity="0.4" />
    <circle cx="28" cy="30" r="1.5" fill="#3b82f6" fillOpacity="0.3" />
    <circle cx="38" cy="28" r="2" fill="#3b82f6" fillOpacity="0.4" />
  </svg>
));
DryerIcon.displayName = 'DryerIcon';

/** RStoic / General reactor — circle with Rx */
export const ReactorIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <rect x="8" y="6" width="32" height="36" rx="4" stroke="currentColor" strokeWidth="2.5" fill="none" />
    <text x="24" y="28" textAnchor="middle" fontSize="12" fontWeight="bold" fill="currentColor" fillOpacity="0.7">Rx</text>
  </svg>
));
ReactorIcon.displayName = 'ReactorIcon';

/** RGibbs — reactor with G */
export const RGibbsIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <rect x="8" y="6" width="32" height="36" rx="4" stroke="currentColor" strokeWidth="2.5" fill="none" />
    <text x="24" y="28" textAnchor="middle" fontSize="14" fontWeight="bold" fill="currentColor" fillOpacity="0.7">G</text>
  </svg>
));
RGibbsIcon.displayName = 'RGibbsIcon';

/** RYield — reactor with Y */
export const RYieldIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <rect x="8" y="6" width="32" height="36" rx="4" stroke="currentColor" strokeWidth="2.5" fill="none" />
    <text x="24" y="28" textAnchor="middle" fontSize="14" fontWeight="bold" fill="currentColor" fillOpacity="0.7">Y</text>
  </svg>
));
RYieldIcon.displayName = 'RYieldIcon';

/** REquil — reactor with balanced scale */
export const REquilIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <rect x="8" y="6" width="32" height="36" rx="4" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* equilibrium arrows */}
    <path d="M16 20 L32 20" stroke="currentColor" strokeWidth="1.5" />
    <polyline points="28,17 32,20 28,23" stroke="currentColor" strokeWidth="1.5" fill="none" />
    <path d="M32 28 L16 28" stroke="currentColor" strokeWidth="1.5" />
    <polyline points="20,25 16,28 20,31" stroke="currentColor" strokeWidth="1.5" fill="none" />
  </svg>
));
REquilIcon.displayName = 'REquilIcon';

/** RBatch — vessel with timer */
export const RBatchIcon = memo(({ size = 28, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 48 48" fill="none" className={className}>
    <rect x="8" y="10" width="32" height="32" rx="4" stroke="currentColor" strokeWidth="2.5" fill="none" />
    {/* agitator */}
    <line x1="24" y1="6" x2="24" y2="22" stroke="currentColor" strokeWidth="1.5" />
    <line x1="18" y1="22" x2="30" y2="22" stroke="currentColor" strokeWidth="2" />
    {/* timer clock */}
    <circle cx="38" cy="10" r="6" stroke="currentColor" strokeWidth="1.5" fill="white" />
    <line x1="38" y1="7" x2="38" y2="10" stroke="currentColor" strokeWidth="1.5" />
    <line x1="38" y1="10" x2="40" y2="12" stroke="currentColor" strokeWidth="1.5" />
  </svg>
));
RBatchIcon.displayName = 'RBatchIcon';

/** Stream arrow — simple arrow for feed/product labels */
export const StreamArrowIcon = memo(({ size = 16, className }: IconProps) => (
  <svg width={size} height={size} viewBox="0 0 24 24" fill="none" className={className}>
    <path d="M4 12 L18 12" stroke="currentColor" strokeWidth="2.5" />
    <polyline points="14,8 18,12 14,16" stroke="currentColor" strokeWidth="2.5" fill="none" />
  </svg>
));
StreamArrowIcon.displayName = 'StreamArrowIcon';
