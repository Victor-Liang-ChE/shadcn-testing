'use client'

import React, { useState, useEffect, useRef, useCallback } from 'react'
import { Skeleton } from '@/components/ui/skeleton'
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from '@/components/ui/tooltip'
import {
  ChevronDown,
  ChevronRight,
  AlertTriangle,
  DollarSign,
  Info,
  User,
} from 'lucide-react'
import { AnimatePresence } from 'framer-motion'
import { cn } from '@/lib/utils'
import { wcBadgeStyle } from '@/lib/wc-colors'
import { supabase } from '@/lib/supabaseClient'
import {
  FighterDetailsModal,
  type ChampionEntry,
} from '@/components/FighterDetailsModal'

// ─── Types ───

type MarketMode = 'sportsbook' | 'exchange'
type BetAdvice = 'advised' | 'caution' | 'limited'
interface OddsPoint {
  timestamp: number
  decimal: number
  american: number
}

interface OddsEvolution {
  mean: OddsPoint[]
  per_book: Record<string, OddsPoint[]>
}

interface FightPrediction {
  rFighter: string
  bFighter: string
  rOdds: number | null
  bOdds: number | null
  rImplied: number
  bImplied: number
  rImpliedNorm: number
  bImpliedNorm: number
  rStreak: number
  bStreak: number
  cbProbRed: number
  lgProbRed: number
  avgProbRed: number
  modelPick: string
  confidence: number
  betSignal: 'bet' | 'none'
  betFighter: string
  betOdds: number
  betProbModel: number
  breakEvenSB: number
  breakEvenEX: number
  evSB: number
  evEX: number
  betAdviceSB: BetAdvice
  betAdviceEX: BetAdvice
  rPerBook: Record<string, number>
  bPerBook: Record<string, number>
  rEvolution?: OddsEvolution
  bEvolution?: OddsEvolution
  rDebut: boolean
  bDebut: boolean
  rRecord: string
  bRecord: string
  weightClass?: string
  modelProb: number
  matchupId?: number
}

interface EventPrediction {
  name: string
  url?: string
  date?: string
  fights: FightPrediction[]
  activeBets: number
  props: PropBet[]
}

interface PropBet {
  rName: string
  bName: string
  rOdds: number | null
  bOdds: number | null
  rPerBook: Record<string, number>
  bPerBook: Record<string, number>
}

interface ModelInfo {
  name: string
  features: string[]
  hp: string
  trainingFights: number
  source?: 'storage' | 'bundle'
}

interface PredictionsData {
  generatedAt: string
  model: ModelInfo
  events: EventPrediction[]
  totalBets: number
  totalFights: number
  oddsScrapedAt: string | null
}

type TottEntry = {
  height?: string | null
  weight?: string | null
  reach?: string | null
  stance?: string | null
  dob?: string | null
}

// ─── Helpers ───

const WC_ABBREV: Record<string, string> = {
  heavyweight: 'HW',
  lightheavyweight: 'LHW',
  middleweight: 'MW',
  welterweight: 'WW',
  lightweight: 'LW',
  featherweight: 'FEW',
  bantamweight: 'BW',
  flyweight: 'FLW',
  women_bantamweight: 'W-BW',
  women_featherweight: 'W-FEW',
  women_flyweight: 'W-FLW',
  women_strawweight: 'W-SW',
}

const WC_FULL_NAMES: Record<string, string> = {
  heavyweight: 'Heavyweight (265 lbs)',
  lightheavyweight: 'Light Heavyweight (205 lbs)',
  middleweight: 'Middleweight (185 lbs)',
  welterweight: 'Welterweight (170 lbs)',
  lightweight: 'Lightweight (155 lbs)',
  featherweight: 'Featherweight (145 lbs)',
  bantamweight: 'Bantamweight (135 lbs)',
  flyweight: 'Flyweight (125 lbs)',
  women_bantamweight: "Women's Bantamweight (135 lbs)",
  women_featherweight: "Women's Featherweight (145 lbs)",
  women_flyweight: "Women's Flyweight (125 lbs)",
  women_strawweight: "Women's Strawweight (115 lbs)",
}

const BOOK_URLS: Record<string, string> = {
  fanduel: 'https://www.fanduel.com',
  draftkings: 'https://www.draftkings.com',
  betmgm: 'https://www.betmgm.com',
  caesars: 'https://www.caesars.com/sportsbook-and-casino',
  'espn bet': 'https://espnbet.com',
  espnbet: 'https://espnbet.com',
  bet365: 'https://www.bet365.com',
  pointsbet: 'https://www.pointsbet.com',
  superbook: 'https://www.superbook.com',
  wynn: 'https://www.wynnbet.com',
  unibet: 'https://www.unibet.com',
  betrivers: 'https://www.betrivers.com',
  barstool: 'https://www.barstoolsportsbook.com',
  kalshi: 'https://www.kalshi.com',
  polymarket: 'https://polymarket.com',
  betfair: 'https://www.betfair.com',
  robinhood: 'https://robinhood.com',
  betway: 'https://betway.com/',
}

function formatOdds(odds: number | null): string {
  if (odds == null) return '—'
  return odds > 0 ? `+${odds}` : `${odds}`
}

function formatTimestamp(ts: number): string {
  const d = new Date(ts)
  return d.toLocaleDateString('en-US', { month: 'short', day: 'numeric' })
}

function formatTimestampFull(ts: number): string {
  const d = new Date(ts)
  const date = d.toLocaleDateString('en-US', { month: 'short', day: 'numeric' })
  const time = d.toLocaleTimeString('en-US', { hour: 'numeric', minute: '2-digit', hour12: true })
  return `${date}, ${time}`
}

function normKey(name: string): string {
  return name
    .normalize('NFD')
    .replace(/[\u0300-\u036f]/g, '')
    .replace(/[.'\u2011]/g, '')
    .toLowerCase()
    .trim()
}

function americanToImplied(odds: number): number {
  return odds < 0 ? Math.abs(odds) / (Math.abs(odds) + 100) * 100 : 100 / (odds + 100) * 100
}

// Odds display that visually centers around the number only (sign floats outside)
function OddsDisplay({ odds, isFav, implied, showAsImplied, className }: { odds: number | null; isFav: boolean; implied?: number; showAsImplied?: boolean; className?: string }) {
  if (odds == null) return <span className={className}>—</span>
  if (showAsImplied && implied != null) {
    return <span className={cn('tabular-nums', isFav ? 'text-yellow-300' : 'text-foreground', className)}>{implied.toFixed(1)}%</span>
  }
  const sign = odds >= 0 ? '+' : '−'
  const abs = Math.abs(odds)
  const oddsEl = (
    <span className={cn('relative inline-flex items-baseline justify-center', isFav ? 'text-yellow-300' : 'text-foreground', className)}>
      <span className="absolute right-full top-1/2 -translate-y-1/2 font-bold opacity-90 text-[0.9em] pr-0.5 pointer-events-none select-none">
        {sign}
      </span>
      <span className="tabular-nums">{abs}</span>
    </span>
  )
  if (implied == null) return oddsEl
  return (
    <TooltipProvider>
      <Tooltip>
        <TooltipTrigger asChild>{oddsEl}</TooltipTrigger>
        <TooltipContent side="bottom">
          <span className="text-xs">{implied.toFixed(1)}% implied odds</span>
        </TooltipContent>
      </Tooltip>
    </TooltipProvider>
  )
}

function InlineOdds({ odds, showAsImplied, className }: { odds: number | null; showAsImplied?: boolean; className?: string }) {
  if (odds == null) return <span className={className}>—</span>
  if (showAsImplied) {
    return <span className={cn('tabular-nums', className)}>{americanToImplied(odds).toFixed(1)}%</span>
  }

  const sign = odds >= 0 ? '+' : '−'
  const abs = Math.abs(odds)
  return (
    <span className={cn('relative inline-flex items-baseline justify-center', className)}>
      <span className="absolute right-full top-1/2 -translate-y-1/2 font-bold opacity-90 text-[0.9em] pr-0.5 pointer-events-none select-none">
        {sign}
      </span>
      <span className="tabular-nums">{abs}</span>
    </span>
  )
}

function formatDob(dob: string): string {
  const d = new Date(dob)
  if (isNaN(d.getTime())) return dob
  const today = new Date()
  let age = today.getFullYear() - d.getUTCFullYear()
  if (today.getMonth() < d.getUTCMonth() || (today.getMonth() === d.getUTCMonth() && today.getDate() < d.getUTCDate())) age--
  return `${d.getUTCMonth() + 1}/${d.getUTCDate()}/${d.getUTCFullYear()} (Age ${age})`
}

// ─── Sparkline Component (interactive + theme-aware) ───

function Sparkline({
  points,
  height = 52,
  color = '#6366f1',
  showAxes = false,
  onHoverChange,
  formatY,
}: {
  points: OddsPoint[]
  height?: number
  color?: string
  showAxes?: boolean
  onHoverChange?: (pt: OddsPoint | null) => void
  formatY?: (value: number) => string
}) {
  const [hoverState, setHoverState] = useState<{ idx: number; x: number } | null>(null)
  const svgRef = useRef<SVGSVGElement>(null)

  const W = 460
  if (!points || points.length < 2) return null
  const americans = points.map((p) => p.american)
  const minVal = Math.min(...americans)
  const maxVal = Math.max(...americans)
  const range = maxVal - minVal || 1

  const lp = showAxes ? 36 : 2
  const rp = 2
  const tp = showAxes ? 8 : 2
  const bp = showAxes ? 18 : 2
  const cW = W - lp - rp
  const cH = height - tp - bp

  const pts = points.map((p, i) => ({
    x: lp + (i / (points.length - 1)) * cW,
    y: tp + (1 - (p.american - minVal) / range) * cH,
    ...p,
  }))

  const pathD = pts.map((pt, i) => `${i === 0 ? 'M' : 'L'}${pt.x.toFixed(1)},${pt.y.toFixed(1)}`).join(' ')
  const first = americans[0]
  const last = americans[americans.length - 1]
  const endColor = color
  const yFormatter = formatY ?? formatOdds

  const handleMouseMove = useCallback((e: React.MouseEvent<SVGSVGElement>) => {
    if (!svgRef.current) return
    const ctm = svgRef.current.getScreenCTM()
    if (!ctm) return
    const svgX = (e.clientX - ctm.e) / ctm.a
    const clampedX = Math.max(lp, Math.min(lp + cW, svgX))
    const frac = Math.max(0, Math.min(1, (clampedX - lp) / cW))
    const idx = Math.round(frac * (pts.length - 1))
    setHoverState({ idx, x: clampedX })
    onHoverChange?.(points[idx])
  }, [pts.length, cW, lp, onHoverChange, points])

  const hPt = hoverState !== null ? pts[hoverState.idx] : null

  return (
    <svg ref={svgRef} width="100%" height={height} viewBox={`0 0 ${W} ${height}`} className="overflow-visible" style={{ cursor: 'crosshair' }} onMouseMove={handleMouseMove} onMouseLeave={() => { setHoverState(null); onHoverChange?.(null) }}>
      {showAxes && (
        <>
          <text x={lp - 3} y={tp + 5} textAnchor="end" style={{ fontSize: 8, fill: 'var(--muted-foreground)' }}>{yFormatter(maxVal)}</text>
          <text x={lp - 3} y={tp + cH} textAnchor="end" style={{ fontSize: 8, fill: 'var(--muted-foreground)' }}>{yFormatter(minVal)}</text>
          <text x={lp} y={height - 1} textAnchor="start" style={{ fontSize: 8, fill: 'var(--muted-foreground)' }}>{formatTimestamp(points[0].timestamp)}</text>
          <text x={lp + cW} y={height - 1} textAnchor="end" style={{ fontSize: 8, fill: 'var(--muted-foreground)' }}>{formatTimestamp(points[points.length - 1].timestamp)}</text>
          <line x1={lp} y1={tp + cH} x2={lp + cW} y2={tp + cH} style={{ stroke: 'var(--border)', strokeWidth: 0.5 }} />
          <line x1={lp} y1={tp} x2={lp} y2={tp + cH} style={{ stroke: 'var(--border)', strokeWidth: 0.5 }} />
        </>
      )}
      <path d={pathD} fill="none" stroke={color} strokeWidth={1.5} />
      <circle cx={pts[0].x} cy={pts[0].y} r={2} fill={color} />
      <circle cx={pts[pts.length - 1].x} cy={pts[pts.length - 1].y} r={2.5} fill={endColor} />
      {/* Invisible full-area hit zone for hover */}
      <rect x={0} y={0} width={W} height={height} fill="transparent" />
      {/* Hover crosshair */}
      {hPt && hoverState && (
        <>
          <circle cx={hPt.x} cy={hPt.y} r={3.5} fill={color} style={{ stroke: 'var(--background)', strokeWidth: 1.5 }} />
        </>
      )}
    </svg>
  )
}

// ─── Per-Book Odds Table ───

function PerBookOdds({
  rPerBook,
  bPerBook,
  rFighter,
  bFighter,
  marketMode,
  rEvolution,
  bEvolution,
}: {
  rPerBook: Record<string, number>
  bPerBook: Record<string, number>
  rFighter: string
  bFighter: string
  marketMode: MarketMode
  rEvolution?: OddsEvolution
  bEvolution?: OddsEvolution
}) {
  const [selectedBook, setSelectedBook] = useState<string | null>(null)
  const [rHoverPt, setRHoverPt] = useState<OddsPoint | null>(null)
  const [bHoverPt, setBHoverPt] = useState<OddsPoint | null>(null)
  const books = [...new Set([...Object.keys(rPerBook), ...Object.keys(bPerBook)])].sort()
  if (books.length === 0) return null

  const isExchange = marketMode === 'exchange'

  const rVals = books.map(b => rPerBook[b]).filter((v): v is number => v != null)
  const bVals = books.map(b => bPerBook[b]).filter((v): v is number => v != null)
  const rAvgOdds = rVals.length > 0 ? Math.round(rVals.reduce((a, b) => a + b) / rVals.length) : null
  const bAvgOdds = bVals.length > 0 ? Math.round(bVals.reduce((a, b) => a + b) / bVals.length) : null

  const getPoints = (book: string, side: 'r' | 'b'): OddsPoint[] | null => {
    const evo = side === 'r' ? rEvolution : bEvolution
    if (!evo) return null
    if (book === '__avg__') return (evo.mean?.length ?? 0) > 0 ? evo.mean : null
    return (evo.per_book[book]?.length ?? 0) > 0 ? evo.per_book[book] : null
  }

  const rLast = rFighter.split(' ').pop() ?? rFighter
  const bLast = bFighter.split(' ').pop() ?? bFighter

  const handleValueClick = (e: React.MouseEvent, book: string) => {
    e.stopPropagation()
    setSelectedBook(prev => prev === book ? null : book)
    setRHoverPt(null)
    setBHoverPt(null)
  }

  const rSelPts = selectedBook ? getPoints(selectedBook, 'r') : null
  const bSelPts = selectedBook ? getPoints(selectedBook, 'b') : null

  return (
    <div className="space-y-2">
      <div className="space-y-0">
        {/* Book rows */}
        {books.map((book) => {
          const rVal = rPerBook[book]
          const bVal = bPerBook[book]
          const rBest = rVal != null && Object.values(rPerBook).every((v) => v == null || v <= rVal)
          const bBest = bVal != null && Object.values(bPerBook).every((v) => v == null || v <= bVal)
          const isSelected = selectedBook === book
          return (
            <div
              key={book}
              className={cn('grid grid-cols-[1fr_auto_1fr] gap-x-3 border-b border-border/30 last:border-0 cursor-pointer', isSelected && 'bg-muted/20')}
              onClick={e => handleValueClick(e, book)}
            >
              <div className={cn('py-0.5 tabular-nums text-xs justify-self-center', rBest ? 'text-foreground' : 'text-muted-foreground')}>
                <InlineOdds odds={rVal} showAsImplied={isExchange} className="text-xs" />
              </div>
              <div className="min-w-[52px] py-0.5 text-center text-xs">
                {(() => {
                  const url = BOOK_URLS[book.toLowerCase()]
                  return url ? (
                    <a href={url} target="_blank" rel="noopener noreferrer" className="text-foreground hover:underline" onClick={e => e.stopPropagation()}>{book}</a>
                  ) : (
                    <span className="text-foreground">{book}</span>
                  )
                })()}
              </div>
              <div className={cn('py-0.5 tabular-nums text-xs justify-self-center', bBest ? 'text-foreground' : 'text-muted-foreground')}>
                <InlineOdds odds={bVal} showAsImplied={isExchange} className="text-xs" />
              </div>
            </div>
          )
        })}
        {/* Average row */}
        {(rAvgOdds != null || bAvgOdds != null) && (
          <div
            className={cn('grid grid-cols-[1fr_auto_1fr] gap-x-3 border-t-2 border-border/60 font-semibold cursor-pointer', selectedBook === '__avg__' && 'bg-muted/20')}
            onClick={e => handleValueClick(e, '__avg__')}
          >
            <div className="py-1 tabular-nums text-xs text-foreground justify-self-center"><InlineOdds odds={rAvgOdds != null ? rAvgOdds : null} showAsImplied={isExchange} className="text-xs" /></div>
            <div className="min-w-[52px] py-1 text-center text-muted-foreground text-[11px] uppercase tracking-wide">Average</div>
            <div className="py-1 tabular-nums text-xs text-foreground justify-self-center"><InlineOdds odds={bAvgOdds != null ? bAvgOdds : null} showAsImplied={isExchange} className="text-xs" /></div>
          </div>
        )}
      </div>
      {selectedBook && (rSelPts || bSelPts) && (
        <div className="pt-2 border-t border-border/30">
          <div className="flex justify-end mb-1">
            <button onClick={e => { e.stopPropagation(); setSelectedBook(null) }} className="text-[10px] text-muted-foreground/50 hover:text-muted-foreground px-1">✕</button>
          </div>
          <div className="grid grid-cols-[1fr_auto_1fr] gap-x-3">
            <div className="flex flex-col items-center">
              {rSelPts ? (
                <>
                  <div className="mb-1 text-center">
                    <div className="text-sm font-semibold tabular-nums text-red-400">
                      <InlineOdds odds={(rHoverPt ?? rSelPts[rSelPts.length - 1]).american} showAsImplied={isExchange} className="text-sm" />
                    </div>
                    <div className="text-[10px] text-muted-foreground">
                      {formatTimestampFull((rHoverPt ?? rSelPts[rSelPts.length - 1]).timestamp)}
                    </div>
                  </div>
                  <Sparkline
                    points={rSelPts}
                    height={52}
                    color="#f87171"
                    showAxes
                    onHoverChange={setRHoverPt}
                    formatY={isExchange ? (v) => `${americanToImplied(v).toFixed(1)}%` : formatOdds}
                  />
                </>
              ) : (
                <div className="text-xs text-muted-foreground/50 text-center py-4">No data</div>
              )}
            </div>
            <div className="min-w-[52px]" />
            <div className="flex flex-col items-center">
              {bSelPts ? (
                <>
                  <div className="mb-1 text-center">
                    <div className="text-sm font-semibold tabular-nums text-blue-400">
                      <InlineOdds odds={(bHoverPt ?? bSelPts[bSelPts.length - 1]).american} showAsImplied={isExchange} className="text-sm" />
                    </div>
                    <div className="text-[10px] text-muted-foreground">
                      {formatTimestampFull((bHoverPt ?? bSelPts[bSelPts.length - 1]).timestamp)}
                    </div>
                  </div>
                  <Sparkline
                    points={bSelPts}
                    height={52}
                    color="#60a5fa"
                    showAxes
                    onHoverChange={setBHoverPt}
                    formatY={isExchange ? (v) => `${americanToImplied(v).toFixed(1)}%` : formatOdds}
                  />
                </>
              ) : (
                <div className="text-xs text-muted-foreground/50 text-center py-4">No data</div>
              )}
            </div>
          </div>
        </div>
      )}
    </div>
  )
}

// ─── Evolution Timeline ───

function EvolutionTimeline({
  evolution,
  fighterName,
  color = '#6366f1',
  labelColor,
}: {
  evolution?: OddsEvolution
  fighterName: string
  color?: string
  labelColor?: string
}) {
  const [showPerBook, setShowPerBook] = useState(false)
  if (!evolution || evolution.mean.length === 0) return null

  const books = Object.keys(evolution.per_book).filter((k) => evolution.per_book[k].length > 0)
  const last = evolution.mean[evolution.mean.length - 1].american

  return (
    <div className="space-y-1">
      <div className="flex items-center justify-between mb-1">
        <span className="text-xs">
          <span className={labelColor ?? 'text-muted-foreground'}>{fighterName}</span>
          <span className="text-muted-foreground"> avg odds evolution</span>
        </span>
        <span className="text-[11px] tabular-nums font-semibold text-foreground">{formatOdds(last)}</span>
      </div>
      <div className="w-full">
        <Sparkline points={evolution.mean} height={52} color={color} showAxes />
      </div>
      {books.length > 0 && (
        <>
          <button
            onClick={() => setShowPerBook(!showPerBook)}
            className="text-[10px] text-muted-foreground/70 hover:text-muted-foreground flex items-center gap-0.5 mt-1"
          >
            {showPerBook ? <ChevronDown className="w-3 h-3" /> : <ChevronRight className="w-3 h-3" />}
            {books.length} books
          </button>
          {showPerBook && (
            <div className="pl-2 space-y-2 mt-1">
              {books.map((book) => (
                <div key={book} className="space-y-0.5">
                  <span className="text-[10px] text-muted-foreground/60">{book}</span>
                  <div className="w-full">
                    <Sparkline points={evolution.per_book[book]} height={52} color="#94a3b8" showAxes />
                  </div>
                </div>
              ))}
            </div>
          )}
        </>
      )}
    </div>
  )
}

// ─── Probability Bar ───

function ProbBar({
  prob,
  label,
  isRed,
}: {
  prob: number
  label: string
  isRed: boolean
}) {
  const pct = Math.round(prob * 100)
  return (
    <div className="flex items-center gap-2 text-xs">
      <span className="w-8 text-right text-muted-foreground tabular-nums">
        {pct}%
      </span>
      <div className="flex-1 h-2 bg-muted rounded-full overflow-hidden">
        <div
          className={cn(
            'h-full rounded-full transition-all',
            isRed ? 'bg-red-500/70' : 'bg-blue-500/70'
          )}
          style={{ width: `${pct}%` }}
        />
      </div>
      <span className="text-[10px] text-muted-foreground">{label}</span>
    </div>
  )
}

// ─── Fighter Hover Tooltip ───

function FighterTooltip({
  name,
  record,
  streak,
  weightClass,
  headshot,
  tott,
  children,
}: {
  name: string
  record: string
  streak: number
  weightClass?: string
  headshot?: string | null
  tott?: TottEntry | null
  children: React.ReactNode
}) {
  const wcAbbrev = weightClass ? WC_ABBREV[weightClass] : undefined
  return (
    <TooltipProvider>
      <Tooltip>
        <TooltipTrigger asChild>{children}</TooltipTrigger>
        <TooltipContent
          side="right"
          className="p-0 overflow-hidden bg-card text-card-foreground border border-border shadow-lg"
        >
          <div className="flex gap-3 p-3 min-w-[200px] max-w-[240px]">
            <div className="w-14 h-14 rounded-full overflow-hidden bg-muted flex items-center justify-center flex-shrink-0 border border-border">
              {headshot ? (
                <img src={headshot} alt={name} className="w-full h-full object-cover object-top" />
              ) : (
                <User className="w-6 h-6 text-muted-foreground" />
              )}
            </div>
            <div className="flex flex-col gap-0.5 text-xs min-w-0 justify-center">
              <div className="font-bold text-sm truncate leading-tight">{name}</div>
              {wcAbbrev && (
                <span
                  className="text-[9px] font-semibold uppercase px-1.5 py-0.5 rounded border w-fit"
                  style={wcBadgeStyle(weightClass!)}
                >
                  {wcAbbrev}
                </span>
              )}
              <div className="text-muted-foreground">UFC {record}</div>
              {streak !== 0 && (
                <div className={cn('font-medium', streak > 0 ? 'text-emerald-500' : 'text-red-400')}>
                  {streak > 0 ? `W${streak}` : `L${Math.abs(streak)}`} streak
                </div>
              )}
              {tott?.height && <div><span className="font-medium">Height:</span> {tott.height}</div>}
              {tott?.reach && <div><span className="font-medium">Reach:</span> {tott.reach}</div>}
              {tott?.stance && <div><span className="font-medium">Stance:</span> {tott.stance}</div>}
              {tott?.dob && <div><span className="font-medium">Born:</span> {formatDob(tott.dob)}</div>}
            </div>
          </div>
        </TooltipContent>
      </Tooltip>
    </TooltipProvider>
  )
}

// ─── Market Toggle ───

function MarketToggle({
  mode,
  onChange,
}: {
  mode: MarketMode
  onChange: (m: MarketMode) => void
}) {
  const [hoveredMode, setHoveredMode] = useState<MarketMode | null>(null)
  const [anchorRect, setAnchorRect] = useState<DOMRect | null>(null)

  const bookHints: Record<MarketMode, string> = {
    sportsbook: 'FanDuel · DraftKings · BetMGM · Caesars',
    exchange: 'Kalshi · Polymarket · Betfair · Robinhood',
  }
  const labels: Record<MarketMode, string> = {
    sportsbook: 'Sportsbook',
    exchange: 'Exchange / No-Vig',
  }

  return (
    <div className="w-full relative">
      <div className="flex w-full rounded-lg overflow-hidden border border-border bg-muted p-0.5 gap-0.5">
        {(['sportsbook', 'exchange'] as MarketMode[]).map((m) => {
          const active = mode === m
          return (
            <button
              key={m}
              onClick={() => onChange(m)}
              onMouseEnter={(e) => { setHoveredMode(m); setAnchorRect(e.currentTarget.getBoundingClientRect()) }}
              onMouseLeave={() => { setHoveredMode(null); setAnchorRect(null) }}
              className={cn(
                'flex-1 py-2 px-3 text-sm font-medium rounded-md transition-all',
                active
                  ? 'bg-background text-foreground shadow-sm'
                  : 'text-muted-foreground hover:text-foreground'
              )}
            >
              {labels[m]}
            </button>
          )
        })}
      </div>
      {/* Speech-bubble tooltip — absolutely positioned, no layout shift */}
      {hoveredMode && anchorRect && (
        <div
          className="fixed z-50 pointer-events-none"
          style={{
            left: anchorRect.left + anchorRect.width / 2,
            top: anchorRect.bottom + 6,
            transform: 'translateX(-50%)',
          }}
        >
          <div className="relative bg-popover border border-border rounded-md px-3 py-1.5 text-xs text-muted-foreground shadow-md whitespace-nowrap">
            {/* Arrow up */}
            <span className="absolute -top-[5px] left-1/2 -translate-x-1/2 w-0 h-0"
              style={{ borderLeft: '5px solid transparent', borderRight: '5px solid transparent', borderBottom: '5px solid var(--border)' }} />
            <span className="absolute -top-[4px] left-1/2 -translate-x-1/2 w-0 h-0"
              style={{ borderLeft: '5px solid transparent', borderRight: '5px solid transparent', borderBottom: '5px solid var(--popover)' }} />
            {bookHints[hoveredMode]}
          </div>
        </div>
      )}
    </div>
  )
}

// ─── Bet Advice Banner ───

function BetAdviceBanner({
  fight,
  mode,
}: {
  fight: FightPrediction
  mode: MarketMode
}) {
  const advice = mode === 'sportsbook' ? fight.betAdviceSB : fight.betAdviceEX
  const ev = mode === 'sportsbook' ? fight.evSB : fight.evEX
  const breakEven = mode === 'sportsbook' ? fight.breakEvenSB : fight.breakEvenEX
  const modelPct = Math.round(fight.betProbModel * 100)
  const edge = modelPct - breakEven

  const styles: Record<
    BetAdvice,
    { bg: string; border: string; icon: string; label: string; labelColor: string; sublabel: string }
  > = {
    advised: {
      bg: 'bg-emerald-950/30',
      border: 'border-emerald-800/50',
      icon: '✓',
      label: 'ADVISED',
      labelColor: 'text-emerald-400',
      sublabel:
        mode === 'sportsbook'
          ? 'Historically +EV range on FanDuel / DraftKings'
          : 'Historically +EV range on exchange / Kalshi / Polymarket',
    },
    caution: {
      bg: 'bg-amber-950/30',
      border: 'border-amber-800/50',
      icon: '⚠',
      label: 'CAUTION',
      labelColor: 'text-amber-400',
      sublabel:
        mode === 'sportsbook'
          ? '+130–175 range was historically -EV on sportsbook'
          : '+130–175 range was slightly -EV; better on exchange than sportsbook',
    },
    limited: {
      bg: 'bg-muted/30',
      border: 'border-border',
      icon: '~',
      label: 'LIMITED DATA',
      labelColor: 'text-muted-foreground',
      sublabel: 'Fewer than 5 historical bets at these odds — insufficient sample',
    },
  }

  const s = styles[advice]

  return (
    <div className={cn('mt-3 rounded-lg overflow-hidden border', s.bg, s.border)}>
      <div className="flex items-center gap-2.5 p-3">
        <DollarSign className={cn('w-4 h-4 flex-shrink-0', s.labelColor)} />
        <div className="flex-1 min-w-0">
          <span className={cn('text-sm font-bold', s.labelColor)}>
            {s.icon} {s.label} · BET {fight.betFighter} ({formatOdds(fight.betOdds)})
          </span>
          <div className="text-xs text-muted-foreground mt-0.5">{s.sublabel}</div>
        </div>
      </div>
      <div className="grid grid-cols-3 border-t border-border/50">
        {[
          { label: 'Est. EV/$100', value: `${ev >= 0 ? '+' : ''}$${ev}`, color: ev >= 0 ? 'text-emerald-400' : 'text-red-400', tip: 'Expected return per $100' },
          { label: 'Model gives', value: `${modelPct}%`, color: 'text-violet-400', tip: `Model probability for ${fight.betFighter}` },
          { label: 'Break-even', value: `${breakEven}%`, color: 'text-muted-foreground', tip: 'Break-even win % at these odds' },
        ].map(({ label, value, color, tip }) => (
          <div key={label} title={tip} className="p-2.5 border-r border-border/50 last:border-r-0 cursor-help">
            <div className="text-[10px] text-muted-foreground uppercase tracking-wide">{label}</div>
            <div className={cn('text-base font-bold mt-0.5', color)}>{value}</div>
            {label !== 'Break-even' && (
              <div className={cn('text-[10px]', edge >= 0 ? 'text-emerald-400' : 'text-red-400')}>
                {edge >= 0 ? '+' : ''}{edge.toFixed(1)}pp edge
              </div>
            )}
          </div>
        ))}
      </div>
      <div className="px-3 py-1.5 border-t border-border/50 text-[10px] text-muted-foreground flex items-center gap-1">
        <Info className="w-3 h-3 flex-shrink-0" />
        {mode === 'sportsbook'
          ? 'Profitable range: dog odds ≤ +130 · Avoid: +130–175 (historically -EV)'
          : 'Exchange pays fair odds. Same profitable range, ~$5/bet extra vs sportsbook.'}
      </div>
    </div>
  )
}

// ─── Fight Card (versus layout) ───

function FightCard({
  fight,
  marketMode,
  onNameClick,
  headshots,
  totts,
}: {
  fight: FightPrediction
  marketMode: MarketMode
  onNameClick: (name: string) => void
  headshots: Map<string, string>
  totts: Map<string, TottEntry>
}) {
  const [expanded, setExpanded] = useState(false)
  const isBet = fight.betSignal === 'bet'

  // rFighter always LEFT (red corner), bFighter always RIGHT (blue corner)
  const rIsFav = fight.rOdds != null && fight.bOdds != null ? fight.rOdds <= fight.bOdds : true

  const rKey = normKey(fight.rFighter)
  const bKey = normKey(fight.bFighter)
  const rWords = fight.rFighter.trim().split(/\s+/)
  const bWords = fight.bFighter.trim().split(/\s+/)
  const rLastKey = normKey(fight.rFighter.split(' ').pop() ?? fight.rFighter)
  const bLastKey = normKey(fight.bFighter.split(' ').pop() ?? fight.bFighter)
  const rHeadshot = headshots.get(rKey) ?? (rWords.length === 1 ? headshots.get(rLastKey) : undefined)
  const bHeadshot = headshots.get(bKey) ?? (bWords.length === 1 ? headshots.get(bLastKey) : undefined)
  const rTott = totts.get(rKey) ?? (rWords.length === 1 ? totts.get(rLastKey) : undefined)
  const bTott = totts.get(bKey) ?? (bWords.length === 1 ? totts.get(bLastKey) : undefined)
  const wcAbbrev = fight.weightClass ? WC_ABBREV[fight.weightClass] : undefined
  const wcLabel = fight.weightClass ? (WC_FULL_NAMES[fight.weightClass] ?? undefined) : undefined

  const rPick = fight.modelPick === fight.rFighter
  const bPick = fight.modelPick === fight.bFighter
  // Triangle size: based on picked fighter's model probability
  const probPct = fight.avgProbRed * 100
  const pickedProb = Math.max(50, Math.min(100, rPick ? probPct : (100 - probPct)))
  // Triangle size scales with model confidence (50–100%).
  // Increase max size so very confident picks are more visually prominent.
  const triH = 20 + Math.round(((pickedProb - 50) / 50) * 56)
  const triW = Math.round(triH * 0.52)

  const rIsDebut = fight.rDebut || fight.rRecord === '—'
  const bIsDebut = fight.bDebut || fight.bRecord === '—'

  return (
    <div className={cn('border-b border-border/60 last:border-b-0')}>
      {/* Versus card row */}
      <button onClick={() => setExpanded(!expanded)} className="w-full px-4 py-3 text-left">
        <div className="grid grid-cols-[1fr_auto_1fr] gap-x-3 items-center">
          {/* RED FIGHTER (left) */}
          <div className="flex flex-col items-center gap-1">
            <div className="w-12 h-12 rounded-full overflow-hidden flex-shrink-0">
              {rHeadshot ? (
                <img src={rHeadshot} alt={fight.rFighter} className="w-full h-full object-cover object-top" />
              ) : (
                <div className="w-full h-full bg-muted flex items-center justify-center rounded-full">
                  <User className="w-5 h-5 text-muted-foreground" />
                </div>
              )}
            </div>
            {rIsDebut ? (
              <span
                onClick={(e) => { e.stopPropagation(); onNameClick(fight.rFighter) }}
                className="font-semibold text-xl leading-tight cursor-pointer hover:underline text-center text-red-400"
              >
                {fight.rFighter}
              </span>
            ) : (
              <FighterTooltip name={fight.rFighter} record={fight.rRecord} streak={fight.rStreak} weightClass={fight.weightClass} headshot={rHeadshot} tott={rTott}>
                <span
                  onClick={(e) => { e.stopPropagation(); onNameClick(fight.rFighter) }}
                  className="font-semibold text-xl leading-tight cursor-pointer hover:underline text-center text-red-400"
                >
                  {fight.rFighter}
                </span>
              </FighterTooltip>
            )}
            {rIsDebut ? (
              <span className="text-xs text-muted-foreground">DEBUT</span>
            ) : (
              <span className="text-xs text-muted-foreground tabular-nums">
                {fight.rRecord}
                {fight.rStreak != null && fight.rStreak !== 0 && (
                  <span className={cn(
                    'px-1.5 py-0.5 rounded ml-1',
                    fight.rStreak > 0 
                      ? 'bg-emerald-500/10 text-emerald-500' 
                      : 'bg-red-500/10 text-red-500'
                  )}>
                    {fight.rStreak > 0 ? `W${fight.rStreak}` : `L${Math.abs(fight.rStreak)}`}
                  </span>
                )}
              </span>
            )}
            <OddsDisplay odds={fight.rOdds} isFav={rIsFav} implied={fight.rImpliedNorm} showAsImplied={marketMode === 'exchange'} className="text-xl font-bold" />
          </div>

          {/* MIDDLE — vs + WC badge, centered; triangles colored by fighter + tooltip */}
          <div className="relative flex flex-col items-center justify-center gap-1.5 min-w-[52px] self-stretch">
            {/* Triangle for rFighter pick: points LEFT, red */}
            {rPick && (
              <TooltipProvider>
                <Tooltip>
                  <TooltipTrigger asChild>
                    <span
                      className="absolute cursor-help"
                      style={{ right: '100%', top: '50%', transform: 'translateY(-50%)', paddingRight: '6px' }}
                    >
                      <svg width={triW} height={triH} viewBox={`0 0 ${triW} ${triH}`}>
                        <polygon points={`0,${triH / 2} ${triW},0 ${triW},${triH}`} fill="#f87171" opacity="0.85" />
                      </svg>
                    </span>
                  </TooltipTrigger>
                  <TooltipContent side="left" className="p-2 min-w-[150px]">
                    <div className="space-y-1 text-xs">
                      <div className="font-medium mb-1.5">
                        Model favors <span className="text-red-400">{fight.rFighter.split(' ').pop()}</span>
                      </div>
                      <div className="flex justify-between gap-3">
                        <span className="text-muted-foreground">CatBoost</span>
                        <span className="tabular-nums">{Math.round(fight.cbProbRed * 100)}%</span>
                      </div>
                      <div className="flex justify-between gap-3">
                        <span className="text-muted-foreground">LightGBM</span>
                        <span className="tabular-nums">{Math.round(fight.lgProbRed * 100)}%</span>
                      </div>
                      <div className="flex justify-between gap-3 font-semibold border-t border-border/50 pt-1 mt-1">
                        <span className="text-muted-foreground">Average</span>
                        <span className="tabular-nums">{Math.round(fight.avgProbRed * 100)}%</span>
                      </div>
                    </div>
                  </TooltipContent>
                </Tooltip>
              </TooltipProvider>
            )}
            <span className="text-xs text-muted-foreground/60 font-medium">vs</span>
            {wcAbbrev ? (
              wcLabel ? (
                <TooltipProvider>
                  <Tooltip>
                    <TooltipTrigger asChild>
                      <span
                        className="text-[9px] font-semibold uppercase px-1.5 py-0.5 rounded border cursor-help"
                        style={wcBadgeStyle(fight.weightClass!)}
                      >
                        {wcAbbrev}
                      </span>
                    </TooltipTrigger>
                    <TooltipContent side="bottom">
                      <span className="text-xs">{wcLabel}</span>
                    </TooltipContent>
                  </Tooltip>
                </TooltipProvider>
              ) : (
                <span className="text-[9px] font-semibold uppercase px-1.5 py-0.5 rounded border" style={wcBadgeStyle(fight.weightClass!)}>
                  {wcAbbrev}
                </span>
              )
            ) : null}
            {/* Triangle for bFighter pick: points RIGHT, blue */}
            {bPick && (
              <TooltipProvider>
                <Tooltip>
                  <TooltipTrigger asChild>
                    <span
                      className="absolute cursor-help"
                      style={{ left: '100%', top: '50%', transform: 'translateY(-50%)', paddingLeft: '6px' }}
                    >
                      <svg width={triW} height={triH} viewBox={`0 0 ${triW} ${triH}`}>
                        <polygon points={`${triW},${triH / 2} 0,0 0,${triH}`} fill="#60a5fa" opacity="0.85" />
                      </svg>
                    </span>
                  </TooltipTrigger>
                  <TooltipContent side="right" className="p-2 min-w-[150px]">
                    <div className="space-y-1 text-xs">
                      <div className="font-medium mb-1.5">
                        Model favors <span className="text-blue-400">{fight.bFighter.split(' ').pop()}</span>
                      </div>
                      <div className="flex justify-between gap-3">
                        <span className="text-muted-foreground">CatBoost</span>
                        <span className="tabular-nums">{Math.round((1 - fight.cbProbRed) * 100)}%</span>
                      </div>
                      <div className="flex justify-between gap-3">
                        <span className="text-muted-foreground">LightGBM</span>
                        <span className="tabular-nums">{Math.round((1 - fight.lgProbRed) * 100)}%</span>
                      </div>
                      <div className="flex justify-between gap-3 font-semibold border-t border-border/50 pt-1 mt-1">
                        <span className="text-muted-foreground">Average</span>
                        <span className="tabular-nums">{Math.round((1 - fight.avgProbRed) * 100)}%</span>
                      </div>
                    </div>
                  </TooltipContent>
                </Tooltip>
              </TooltipProvider>
            )}
          </div>

          {/* BLUE FIGHTER (right) */}
          <div className="flex flex-col items-center gap-1">
            <div className="w-12 h-12 rounded-full overflow-hidden flex-shrink-0">
              {bHeadshot ? (
                <img src={bHeadshot} alt={fight.bFighter} className="w-full h-full object-cover object-top" />
              ) : (
                <div className="w-full h-full bg-muted flex items-center justify-center rounded-full">
                  <User className="w-5 h-5 text-muted-foreground" />
                </div>
              )}
            </div>
            {bIsDebut ? (
              <span
                onClick={(e) => { e.stopPropagation(); onNameClick(fight.bFighter) }}
                className="font-semibold text-xl leading-tight cursor-pointer hover:underline text-center text-blue-400"
              >
                {fight.bFighter}
              </span>
            ) : (
              <FighterTooltip name={fight.bFighter} record={fight.bRecord} streak={fight.bStreak} weightClass={fight.weightClass} headshot={bHeadshot} tott={bTott}>
                <span
                  onClick={(e) => { e.stopPropagation(); onNameClick(fight.bFighter) }}
                  className="font-semibold text-xl leading-tight cursor-pointer hover:underline text-center text-blue-400"
                >
                  {fight.bFighter}
                </span>
              </FighterTooltip>
            )}
            {bIsDebut ? (
              <span className="text-xs text-muted-foreground">DEBUT</span>
            ) : (
              <span className="text-xs text-muted-foreground tabular-nums">
                {fight.bRecord}
                {fight.bStreak != null && fight.bStreak !== 0 && (
                  <span className={cn(
                    'px-1.5 py-0.5 rounded ml-1',
                    fight.bStreak > 0 
                      ? 'bg-emerald-500/10 text-emerald-500' 
                      : 'bg-red-500/10 text-red-500'
                  )}>
                    {fight.bStreak > 0 ? `W${fight.bStreak}` : `L${Math.abs(fight.bStreak)}`}
                  </span>
                )}
              </span>
            )}
            <OddsDisplay odds={fight.bOdds} isFav={!rIsFav} implied={fight.bImpliedNorm} showAsImplied={marketMode === 'exchange'} className="text-xl font-bold" />
          </div>
        </div>

        {/* Expand chevron */}
        <div className="flex justify-center mt-2">
          <ChevronDown className={cn('w-4 h-4 text-muted-foreground transition-transform', expanded && 'rotate-180')} />
        </div>
      </button>

      {/* Expanded details */}
      {expanded && (
        <div className="px-4 pb-4 space-y-4 border-t border-border/20">
          {isBet && <BetAdviceBanner fight={fight} mode={marketMode} />}

          {(Object.keys(fight.rPerBook).length > 0 || Object.keys(fight.bPerBook).length > 0) && (
            <PerBookOdds
              rPerBook={fight.rPerBook}
              bPerBook={fight.bPerBook}
              rFighter={fight.rFighter}
              bFighter={fight.bFighter}
              marketMode={marketMode}
              rEvolution={fight.rEvolution}
              bEvolution={fight.bEvolution}
            />
          )}
        </div>
      )}
    </div>
  )
}

// ─── Props Row ───

function PropRow({ prop }: { prop: PropBet }) {
  const [expanded, setExpanded] = useState(false)
  const books = [...new Set([...Object.keys(prop.rPerBook), ...Object.keys(prop.bPerBook)])].sort()
  return (
    <div className="border-b border-border/30 last:border-b-0">
      <button
        onClick={() => setExpanded(!expanded)}
        className="w-full px-4 py-2.5 grid gap-3 items-center text-left"
        style={{ gridTemplateColumns: '1fr auto auto auto' }}
      >
        <div className="min-w-0">
          <div className="text-sm text-foreground truncate">{prop.rName}</div>
          <div className="text-sm text-muted-foreground truncate mt-0.5">{prop.bName}</div>
        </div>
        <div className="font-mono text-xs text-right">
          <div className={cn(prop.rOdds != null && prop.rOdds > 0 ? 'text-blue-400' : 'text-amber-400')}>{formatOdds(prop.rOdds)}</div>
          <div className={cn('mt-0.5', prop.bOdds != null && prop.bOdds > 0 ? 'text-blue-400' : 'text-amber-400')}>{formatOdds(prop.bOdds)}</div>
        </div>
        <div className="text-[11px] text-muted-foreground/60 w-9 text-center">{books.length > 0 ? `${books.length}bk` : ''}</div>
        <ChevronDown className={cn('w-3.5 h-3.5 text-muted-foreground transition-transform', expanded && 'rotate-180')} />
      </button>
      {expanded && (
        <div className="px-4 pb-3">
          {books.length > 0 ? (
            <table className="w-full text-xs border-collapse">
              <thead>
                <tr className="text-muted-foreground">
                  <th className="text-left py-1 pr-2 font-medium">Book</th>
                  <th className="text-right py-1 px-2 font-medium text-foreground/80">Yes</th>
                  <th className="text-right py-1 pl-2 font-medium">No</th>
                </tr>
              </thead>
              <tbody>
                {books.map((book) => (
                  <tr key={book} className="border-t border-border/20">
                    <td className="py-1 pr-2 text-muted-foreground">{book}</td>
                    <td className="py-1 px-2 text-right tabular-nums font-mono text-foreground/80">{prop.rPerBook[book] != null ? formatOdds(prop.rPerBook[book]) : '—'}</td>
                    <td className="py-1 pl-2 text-right tabular-nums font-mono text-muted-foreground">{prop.bPerBook[book] != null ? formatOdds(prop.bPerBook[book]) : '—'}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          ) : (
            <p className="text-xs text-muted-foreground">No per-book breakdown available</p>
          )}
        </div>
      )}
    </div>
  )
}

// ─── Event Section ───

function EventSection({
  event,
  isFirst,
  marketMode,
  onNameClick,
  headshots,
  totts,
}: {
  event: EventPrediction
  isFirst: boolean
  marketMode: MarketMode
  onNameClick: (name: string) => void
  headshots: Map<string, string>
  totts: Map<string, TottEntry>
}) {
  const [collapsed, setCollapsed] = useState(!isFirst)
  const [tab, setTab] = useState<'fights' | 'props'>('fights')
  const n = event.fights.length
  const fightLabel = `${n} ${n === 1 ? 'fight' : 'fights'}`

  return (
    <div className={cn(
      'rounded-xl overflow-hidden mb-4 border border-border bg-card'
    )}>
      {/* Header */}
      <button
        onClick={() => setCollapsed(!collapsed)}
        className="w-full px-4 py-3.5 flex justify-between items-center border-b border-border/50 bg-card/50"
      >
        <div>
          <h2 className="text-base font-semibold text-foreground m-0 flex items-center gap-2">
            {event.url ? (
              <a
                href={event.url}
                target="_blank"
                rel="noopener noreferrer"
                onClick={(e) => e.stopPropagation()}
                className="text-foreground no-underline hover:underline"
              >
                {event.name}
                {event.date && (
                  <span className="text-muted-foreground text-sm font-normal mx-2">{event.date}</span>
                )}
                <span className="text-muted-foreground/50 text-xs">↗</span>
              </a>
            ) : (
              <>
                {event.name}
                {event.date && (
                  <span className="text-muted-foreground text-sm font-normal">{event.date}</span>
                )}
              </>
            )}
          </h2>
          <div className="flex gap-2 mt-1 items-center">
            {collapsed && (
              <>
                <span className="text-xs text-muted-foreground">{fightLabel}</span>
                {event.props.length > 0 && (
                  <span className="text-xs text-muted-foreground/60">· {event.props.length} props</span>
                )}
              </>
            )}
          </div>
        </div>
        <div className="flex gap-2 items-center">
          {event.activeBets > 0 && (
            <span className="bg-emerald-900/60 text-emerald-400 text-xs font-semibold px-2.5 py-1 rounded-full border border-emerald-800/50">
              {event.activeBets} BET{event.activeBets > 1 ? 'S' : ''}
            </span>
          )}
          <ChevronDown className={cn('w-4 h-4 text-muted-foreground transition-transform', collapsed && '-rotate-90')} />
        </div>
      </button>

      {!collapsed && (
        <>
          {/* Tabs: Fights / Props (perfectly aligned with columns below) */}
          <div className="border-b border-border/50 bg-muted/20 px-4">
            {event.props.length > 0 ? (
              <div className="grid grid-cols-[1fr_auto_1fr] gap-x-3 items-center">
                <button
                  onClick={() => setTab('fights')}
                  className={cn(
                    'py-2 text-xs font-medium transition-colors',
                    tab === 'fights' ? 'text-foreground bg-background/50 border-b-2 border-b-foreground' : 'text-muted-foreground hover:text-foreground'
                  )}
                >
                  Fights ({event.fights.length})
                </button>
                <div className="flex justify-center min-w-[52px]">
                  <div className="w-px h-4 bg-border/30" />
                </div>
                <button
                  onClick={() => setTab('props')}
                  className={cn(
                    'py-2 text-xs font-medium transition-colors',
                    tab === 'props' ? 'text-foreground bg-background/50 border-b-2 border-b-foreground' : 'text-muted-foreground hover:text-foreground'
                  )}
                >
                  Props ({event.props.length})
                </button>
              </div>
            ) : (
              <div className="grid grid-cols-[1fr_auto_1fr] gap-x-3 items-center">
                <div className="py-2 text-xs font-medium text-foreground bg-background/50 border-b-2 border-b-foreground text-center">
                  Fights ({event.fights.length})
                </div>
                <div className="min-w-[52px]" />
                <div />
              </div>
            )}
          </div>

          {tab === 'fights' && (
            <div>
              {event.fights.map((fight, i) => (
                <FightCard
                  key={`${fight.rFighter}-${fight.bFighter}-${i}`}
                  fight={fight}
                  marketMode={marketMode}
                  onNameClick={onNameClick}
                  headshots={headshots}
                  totts={totts}
                />
              ))}
            </div>
          )}

          {tab === 'props' && (
            <div>
              {event.props.map((prop, i) => (
                <PropRow key={`${prop.rName}-${prop.bName}-${i}`} prop={prop} />
              ))}
            </div>
          )}
        </>
      )}
    </div>
  )
}

// ─── Stat Card ───

function StatCard({
  label,
  value,
  accent,
  sub,
}: {
  label: string
  value: string | number
  accent?: string
  sub?: string
}) {
  return (
    <div className="bg-card border border-border rounded-lg px-4 py-3 flex flex-col items-center text-center justify-center min-h-[76px]">
      <div className="text-xs text-muted-foreground uppercase tracking-wide">{label}</div>
      <div className="text-2xl font-bold mt-1" style={{ color: accent || 'var(--foreground)' }}>{value}</div>
      {sub && <div className="text-[10px] text-muted-foreground/60 mt-0.5">{sub}</div>}
    </div>
  )
}

// ─── Loading Skeleton ───

function LoadingSkeleton() {
  return (
    <div className="space-y-6">
      {Array.from({ length: 3 }).map((_, i) => (
        <div key={i} className="space-y-2">
          <Skeleton className="h-10 rounded-lg" />
          {Array.from({ length: 4 }).map((_, j) => (
            <Skeleton key={j} className="h-16 rounded-lg" />
          ))}
        </div>
      ))}
    </div>
  )
}

// ─── Main Page ───

export default function UfcPredictionsPage() {
  const [data, setData] = useState<PredictionsData | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [marketMode, setMarketMode] = useState<MarketMode>('sportsbook')
  const [selectedFighter, setSelectedFighter] = useState<ChampionEntry | null>(null)
  const [headshots, setHeadshots] = useState<Map<string, string>>(new Map())
  const [totts, setTotts] = useState<Map<string, TottEntry>>(new Map())

  useEffect(() => {
    fetch('/api/ufc-predictions')
      .then((r) => {
        if (!r.ok) throw new Error(`HTTP ${r.status}`)
        return r.json()
      })
      .then((json) => {
        if (json.error) throw new Error(json.error)
        setData(json)
      })
      .catch((e) => setError(e.message))
      .finally(() => setLoading(false))
  }, [])

  // Fetch headshots + TOTT once data is loaded
  useEffect(() => {
    if (!data) return

    fetch('/api/ufc-tott')
      .then((r) => (r.ok ? r.json() : null))
      .then((json) => {
        if (!json?.data) return
        const map = new Map<string, TottEntry>()
        for (const r of json.data as { FIGHTER: string; HEIGHT?: string; WEIGHT?: string; REACH?: string; STANCE?: string; DOB?: string }[]) {
          map.set(normKey(r.FIGHTER), { height: r.HEIGHT, weight: r.WEIGHT, reach: r.REACH, stance: r.STANCE, dob: r.DOB })
        }
        setTotts(map)
      })
      .catch(() => { /* non-critical */ })

    // Per-fighter ilike queries (same approach as FighterDetailsModal — works through proxy)
    const allFighters = data.events.flatMap(e =>
      e.fights.flatMap(f => [f.rFighter, f.bFighter])
    )
    const uniqueFighters = [...new Set(allFighters)]

    Promise.all(uniqueFighters.map(name => {
      const parts = name.trim().split(/\s+/)
      const first = parts[0] ?? ''
      const last = parts[parts.length - 1] ?? ''
      return supabase
        .from('ufc_champion_images')
        .select('FIRST, LAST, HEADSHOT')
        .ilike('FIRST', `%${first}%`)
        .ilike('LAST', `%${last}%`)
        .limit(1)
        .then(({ data: r }: { data: { FIRST: string; LAST: string; HEADSHOT?: string }[] | null }) => ({ name, row: r?.[0] ?? null }))
    })).then(results => {
      const map = new Map<string, string>()
      for (const { name, row } of results) {
        if (!row?.HEADSHOT) continue
        const k = normKey(name)
        map.set(k, row.HEADSHOT)
        const lastOnly = normKey(name.split(' ').pop() ?? name)
        if (!map.has(lastOnly)) map.set(lastOnly, row.HEADSHOT)
      }
      if (map.size > 0) setHeadshots(map)
    }).catch(() => { /* non-critical */ })
  }, [data])

  const handleNameClick = useCallback((name: string) => {
    setSelectedFighter({ name } as ChampionEntry)
  }, [])

  if (loading) {
    return (
      <main className="max-w-[1200px] mx-auto px-4 pt-6">
        <LoadingSkeleton />
      </main>
    )
  }

  if (error) {
    return (
      <main className="max-w-[1200px] mx-auto px-4 pt-6">
        <div className="bg-destructive/10 border border-destructive/50 rounded-lg p-4 flex items-center gap-3">
          <AlertTriangle className="w-5 h-5 text-destructive flex-shrink-0" />
          <div>
            <p className="font-semibold text-destructive m-0">Failed to load</p>
            <p className="text-muted-foreground m-0 mt-1 text-sm">{error}</p>
          </div>
        </div>
      </main>
    )
  }

  if (!data) return null

  return (
    <main className="max-w-[1200px] mx-auto px-4 pt-4 pb-16">
      {/* Market toggle — sticky at top */}
      <div className="sticky top-4 z-40 w-full mb-6">
        <div className="w-full bg-background/80 backdrop-blur-md">
          <MarketToggle mode={marketMode} onChange={setMarketMode} />
        </div>
      </div>

      {/* Events */}
      <div>
        {data.events.map((event, i) => (
          <EventSection
            key={`${event.name}-${i}`}
            event={event}
            isFirst={i === 0}
            marketMode={marketMode}
            onNameClick={handleNameClick}
            headshots={headshots}
            totts={totts}
          />
        ))}
      </div>

      {/* Fighter Details Modal */}
      <AnimatePresence>
        {selectedFighter && (
          <FighterDetailsModal
            champion={selectedFighter}
            onClose={() => setSelectedFighter(null)}
            hideElo={false}
          />
        )}
      </AnimatePresence>
    </main>
  )
}
