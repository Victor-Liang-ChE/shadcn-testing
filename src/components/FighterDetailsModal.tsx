'use client';

import React, { useState, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { X, ArrowLeft, Loader2, Crown, User, ShieldCheck, ShieldOff } from 'lucide-react';
import { Skeleton } from '@/components/ui/skeleton';
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from '@/components/ui/tooltip';
import { supabase, FighterDetails, FightResult, UfcChampionImages } from '@/lib/supabaseClient';
import { cn } from '@/lib/utils';
import { wcBadgeStyle } from '@/lib/wc-colors';

/* ------------------------------------------------------------------ */
/*  TYPES                                                             */
/* ------------------------------------------------------------------ */
interface ProcessedFightResult extends FightResult {
  resultForFighter: 'win' | 'loss' | 'draw' | 'nc' | 'unknown';
  decisionDetails?: string;
  finishDetails?: string;
  defenseOutcome?: 'success' | 'fail' | null;
}

export interface ChampionEntry {
  name: string;
  reignStart?: string;
  reignEnd?: string | null;
  notes?: string;
  imageUrl?: string;    // full-body MAINSHOT (for champion cards)
  headshotUrl?: string; // square HEADSHOT (for modal circle)
  eloHistory?: Array<{ date: string; elo: number; event: string; eloBeforeFight?: number }>; // optional Elo data from UFC Elo page
}

/* ------------------------------------------------------------------ */
/*  HELPERS                                                           */
/* ------------------------------------------------------------------ */
function normalizeName(name: string): string {
  return name
    .normalize('NFD')
    .replace(/[\u0300-\u036f]/g, '')
    .replace(/[.'‑]/g, '')
    .trim();
}

function parseFighterName(full: string): { firstName: string; lastName: string } {
  const parts = full.trim().split(/\s+/);
  return { firstName: parts.slice(0, -1).join(' '), lastName: parts.at(-1) || '' };
}

function normalizeWeightClass(raw?: string | null): string {
  if (!raw) return '';
  let normalized = normalizeName(raw)
    .toLowerCase()
    .replace(/(title|championship|bout|world|undisputed|interim)/g, '')
    .replace(/^ufc\s*/, '')
    .replace(/\s+/g, ' ')
    .trim();
  // Map "womens X" → "women_X" so it matches the KNOWN_DIVISIONS/WC_ORDER keys
  normalized = normalized.replace(/^womens\s+/, 'women_');
  return normalized;
}

function extractOpponent(boutString: string, normalizedFighterName: string): string {
  if (!boutString || !normalizedFighterName) return 'Opponent Unknown';
  const parts = boutString.split(/\s+vs\.?\s+|\s+vs\s+/i);
  if (parts.length === 2) {
    const n1 = normalizeName(parts[0]).toLowerCase();
    const n2 = normalizeName(parts[1]).toLowerCase();
    const n = normalizedFighterName.toLowerCase();
    if (n1.includes(n) || n.includes(n1)) return parts[1].trim();
    if (n2.includes(n) || n.includes(n2)) return parts[0].trim();
  }
  return 'Opponent Unknown';
}

function extractDecisionDetails(details: string | null | undefined): string | undefined {
  if (!details) return undefined;
  const scores = details.match(/\d{1,2}\s*-\s*\d{1,2}/g);
  if (scores && scores.length >= 3) return scores.join(', ');
  return undefined;
}

function getResultForFighter(fight: FightResult, normalizedFighterName: string): ProcessedFightResult['resultForFighter'] {
  if (!fight.BOUT || !fight.OUTCOME) return 'unknown';
  if (fight.DETAILS?.toLowerCase().includes('overturned')) return 'nc';
  if (fight.METHOD?.toLowerCase().includes('no contest') || fight.METHOD?.toLowerCase() === 'nc') return 'nc';
  const outcomeRaw = fight.OUTCOME.toLowerCase().trim();
  if (outcomeRaw === 'nc' || outcomeRaw === 'nc/nc') return 'nc';
  if (outcomeRaw === 'd' || outcomeRaw === 'd/d' || outcomeRaw === 'draw') return 'draw';
  const outcomeParts = fight.OUTCOME.split('/');
  const fighters = fight.BOUT.split(/ vs\.? | vs /i);
  if (fighters.length !== 2) return 'unknown';
  const f1 = normalizeName(fighters[0].trim());
  const f2 = normalizeName(fighters[1].trim());
  const norm = normalizedFighterName.toLowerCase();
  if (f1.toLowerCase() === norm) {
    if (outcomeParts[0]?.toLowerCase() === 'w') return 'win';
    if (outcomeParts[0]?.toLowerCase() === 'l') return 'loss';
  } else if (f2.toLowerCase() === norm) {
    if (outcomeParts[1]?.toLowerCase() === 'w') return 'win';
    if (outcomeParts[1]?.toLowerCase() === 'l') return 'loss';
    if (outcomeParts.length === 1 && outcomeParts[0]?.toLowerCase() === 'w') return 'loss';
    if (outcomeParts.length === 1 && outcomeParts[0]?.toLowerCase() === 'l') return 'win';
  }
  return 'unknown';
}

function calculateAge(dob: string | null | undefined): number | null {
  if (!dob) return null;
  try {
    const birth = new Date(dob);
    const today = new Date();
    let age = today.getFullYear() - birth.getFullYear();
    const m = today.getMonth() - birth.getMonth();
    if (m < 0 || (m === 0 && today.getDate() < birth.getDate())) age--;
    return age;
  } catch { return null; }
}

function wcLbsLabel(key: string): string {
  const labels: Record<string, string> = {
    heavyweight: 'Heavyweight (265 lbs)',
    lightheavyweight: 'Light Heavyweight (205 lbs)',
    middleweight: 'Middleweight (185 lbs)',
    welterweight: 'Welterweight (170 lbs)',
    lightweight: 'Lightweight (155 lbs)',
    featherweight: 'Featherweight (145 lbs)',
    bantamweight: 'Bantamweight (135 lbs)',
    flyweight: 'Flyweight (125 lbs)',
    women_strawweight: "Women's Strawweight (115 lbs)",
    women_flyweight: "Women's Flyweight (125 lbs)",
    women_bantamweight: "Women's Bantamweight (135 lbs)",
    women_featherweight: "Women's Featherweight (145 lbs)",
    catchweight: 'Catchweight',
  };
  const k = key.toLowerCase().replace(/\s+/g, '');
  return labels[k] || key;
}

function abbreviateWeightClass(key: string): string {
  const abbr: Record<string, string> = {
    heavyweight: 'HW', lightheavyweight: 'LHW', middleweight: 'MW',
    welterweight: 'WW', lightweight: 'LW', featherweight: 'FeW',
    bantamweight: 'BW', flyweight: 'FLW', women_strawweight: 'W-SW',
    women_flyweight: 'W-FLW', women_bantamweight: 'W-BW', women_featherweight: 'W-FEW', catchweight: 'CW',
  };
  return abbr[key.toLowerCase().replace(/\s+/g, '')] || key.toUpperCase();
}

// Canonical weight classes in order from lightest to heaviest.
// Used to preserve consistent sorting in the UI.
const WEIGHT_CLASS_ORDER = [
  'women_strawweight',
  'women_flyweight',
  'women_bantamweight',
  'women_featherweight',
  'flyweight',
  'bantamweight',
  'featherweight',
  'lightweight',
  'welterweight',
  'middleweight',
  'lightheavyweight',
  'heavyweight',
] as const;

const KNOWN_DIVISIONS = new Set<string>(WEIGHT_CLASS_ORDER);
function isKnownDivision(key: string): boolean {
  return KNOWN_DIVISIONS.has(key.toLowerCase().replace(/\s+/g, ''));
}

function getWeightClassName(key: string): string {
  const names: Record<string, string> = {
    heavyweight: 'Heavyweight', lightheavyweight: 'Light Heavyweight',
    middleweight: 'Middleweight', welterweight: 'Welterweight',
    lightweight: 'Lightweight', featherweight: 'Featherweight',
    bantamweight: 'Bantamweight', flyweight: 'Flyweight',
    women_strawweight: "Women's Strawweight", women_flyweight: "Women's Flyweight",
    women_bantamweight: "Women's Bantamweight", women_featherweight: "Women's Featherweight", catchweight: 'Catchweight',
  };
  return names[key.toLowerCase().replace(/\s+/g, '')] || key;
}

function parseLocalDate(dateStr: string): Date {
  // Parse YYYY-MM-DD as local midnight (avoids UTC timezone shift for users west of UTC)
  const m = dateStr.match(/^(\d{4})-(\d{2})-(\d{2})/);
  if (m) return new Date(parseInt(m[1]), parseInt(m[2]) - 1, parseInt(m[3]));
  return new Date(dateStr);
}

function diffInDays(a: Date, b: Date): number {
  // Use floor so that "today at 4pm" vs "yesterday midnight" = 1 day, not 2
  return Math.floor((b.getTime() - a.getTime()) / 86_400_000);
}

function humanDuration(days: number): string {
  const y = Math.floor(days / 365);
  const mo = Math.floor((days % 365) / 30);
  const d = days % 30;
  return [y ? `${y} yr${y > 1 ? 's' : ''}` : '', mo ? `${mo} mo${mo > 1 ? 's' : ''}` : '', d ? `${d} day${d > 1 ? 's' : ''}` : '']
    .filter(Boolean).join(' ') || '0 days';
}

function createAthleteUrl(name: string): string {
  if (!name) return '#';
  let slug = name.normalize('NFD').replace(/[\u0300-\u036f]/g, '').replace(/\./g, '').trim().toLowerCase().replace(/\s+/g, '-');
  slug = slug.replace(/^b-?j-/, 'bj-');
  return `https://www.ufc.com/athlete/${slug}`;
}

/* ------------------------------------------------------------------ */
/*  COMPONENT                                                         */
/* ------------------------------------------------------------------ */
interface FighterDetailsModalProps {
  champion: ChampionEntry;
  onClose: () => void;
  onBack?: () => void; // present when navigated from another fighter modal
  eloLookup?: (name: string) => Array<{ date: string; elo: number; event: string }> | undefined;
  hideElo?: boolean;
}

export function FighterDetailsModal({ champion, onClose, onBack, eloLookup, hideElo }: FighterDetailsModalProps) {
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [fighterDetails, setFighterDetails] = useState<FighterDetails | null>(null);
  const [fighterImages, setFighterImages] = useState<UfcChampionImages | null>(null);
  const [fighterHistory, setFighterHistory] = useState<ProcessedFightResult[]>([]);
  const [eventDates, setEventDates] = useState<Record<string, string>>({});
  const [selectedOpponent, setSelectedOpponent] = useState<string | null>(null);
  const [fetchedEloHistory, setFetchedEloHistory] = useState<Array<{ date: string; elo: number; event: string; eloBeforeFight?: number }> | undefined>(undefined);
  // Start true so spinner shows immediately when Elo isn't already provided
  const [fetchingElo, setFetchingElo] = useState<boolean>(
    () => !champion.eloHistory || champion.eloHistory.length === 0
  );

  useEffect(() => {
    const original = champion.name.trim();
    const normName = normalizeName(original);
    // Strip only accent marks (not other punctuation) for DB ilike queries which are accent-sensitive in PostgreSQL
    const stripAccents = (s: string) => s.normalize('NFD').replace(/[\u0300-\u036f]/g, '').trim();
    const queryName = stripAccents(original);
    // Use first word and last word of the accent-stripped name so "Benoit Saint Denis" → FIRST=%Benoit%, LAST=%Denis%
    // This handles multi-word last/first names since DB FIRST/LAST ilike partial-matches
    const qParts = queryName.split(/\s+/);
    const firstName = qParts[0] || '';
    const lastName = qParts[qParts.length - 1] || '';

    setIsLoading(true);
    setError(null);

    Promise.all([
      supabase.from('ufc_fighter_tott').select('*').ilike('FIGHTER', `%${queryName}%`).limit(1).maybeSingle(),
      supabase.from('ufc_champion_images').select('*').ilike('FIRST', `%${firstName}%`).ilike('LAST', `%${lastName}%`).limit(1).maybeSingle(),
      supabase.from('ufc_fight_results').select('*').ilike('BOUT', `%${queryName}%`).order('id', { ascending: true }),
      supabase.from('ufc_event_details').select('EVENT,DATE').limit(2000),
    ])
      .then(([detailsRes, imagesRes, fightsRes, eventDatesRes]) => {
        setFighterDetails((detailsRes.data as FighterDetails | null) ?? null);

        // Prefer champion headshotUrl (circle crop), then HEADSHOT from DB
        const dbImages = imagesRes.data as UfcChampionImages | null;
        if (champion.headshotUrl && dbImages) {
          setFighterImages({ ...dbImages, HEADSHOT: champion.headshotUrl });
        } else if (champion.headshotUrl) {
          setFighterImages({ id: 0, HEADSHOT: champion.headshotUrl });
        } else {
          setFighterImages(dbImages);
        }

        const fights = ((fightsRes.data as FightResult[]) || []);
        const championStatus: Record<string, boolean> = {};
        const processed: ProcessedFightResult[] = [];

        for (let idx = fights.length - 1; idx >= 0; idx--) {
          const fight = fights[idx];
          const resultForFighter = getResultForFighter(fight, normName);
          const div = normalizeWeightClass(fight.WEIGHTCLASS);
          const isTitle = fight.WEIGHTCLASS?.toLowerCase().includes('title') ?? false;
          let defenseOutcome: 'success' | 'fail' | null = null;

          if (isTitle && div) {
            if (championStatus[div]) {
              if (resultForFighter === 'loss') { defenseOutcome = 'fail'; championStatus[div] = false; }
              else if (resultForFighter === 'win') { defenseOutcome = 'success'; }
            } else {
              if (resultForFighter === 'win') championStatus[div] = true;
            }
          }

          let decisionDetails: string | undefined;
          let finishDetails: string | undefined;
          if (fight.METHOD?.toLowerCase().includes('decision')) {
            decisionDetails = extractDecisionDetails(fight.DETAILS);
          } else if (fight.DETAILS?.trim()) {
            if (!/\d{1,2}\s*-\s*\d{1,2}/.test(fight.DETAILS)) finishDetails = fight.DETAILS.trim();
            else decisionDetails = extractDecisionDetails(fight.DETAILS);
          }

          processed.unshift({
            ...fight,
            OPPONENT: extractOpponent(fight.BOUT, normName),
            resultForFighter,
            decisionDetails,
            DETAILS: finishDetails,
            defenseOutcome,
          });
        }
        setFighterHistory(processed);

        const datesMap: Record<string, string> = {};
        if (!eventDatesRes.error && eventDatesRes.data) {
          for (const row of eventDatesRes.data as { EVENT: string; DATE: string | null }[]) {
            if (row.EVENT && row.DATE) {
              datesMap[row.EVENT] = row.DATE;
              datesMap[row.EVENT.toLowerCase().trim()] = row.DATE;
            }
          }
        }
        setEventDates(datesMap);
      })
      .catch((err) => setError(err?.message ?? 'Unknown error'))
      .finally(() => setIsLoading(false));
  }, [champion]);

  // Self-fetch Elo for fighters not in the pre-loaded batch (e.g. second-hop opponents)
  useEffect(() => {
    if (champion.eloHistory && champion.eloHistory.length > 0) {
      setFetchedEloHistory(undefined);
      setFetchingElo(false);
      return;
    }
    setFetchedEloHistory(undefined);
    setFetchingElo(true);
    const name = champion.name.trim();
    if (!name) { setFetchingElo(false); return; }
    let cancelled = false;
    fetch(`/api/ufc-elo?name=${encodeURIComponent(name)}`)
      .then(r => r.ok ? r.json() : null)
      .then(data => {
        if (!cancelled) {
          if (data?.history?.length > 0) setFetchedEloHistory(data.history);
          setFetchingElo(false);
        }
      })
      .catch(() => { if (!cancelled) setFetchingElo(false); });
    return () => { cancelled = true; };
  }, [champion.name, champion.eloHistory]);

  // When an opponent is selected, replace this modal — back arrow returns here, X closes all
  if (selectedOpponent) {
    return (
      <FighterDetailsModal
        key={selectedOpponent}
        champion={{ name: selectedOpponent, eloHistory: eloLookup?.(selectedOpponent) }}
        onClose={onClose}
        onBack={() => setSelectedOpponent(null)}
        eloLookup={eloLookup}
        hideElo={hideElo}
      />
    );
  }

  const activeEloHistory = (champion.eloHistory && champion.eloHistory.length > 0)
    ? champion.eloHistory
    : fetchedEloHistory;

  return (
    <AnimatePresence>
      <motion.div
        key="overlay-backdrop"
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        exit={{ opacity: 0 }}
        transition={{ duration: 0.25 }}
        className="fixed inset-0 z-40 flex items-center justify-center bg-black/60 backdrop-blur-sm p-4 md:p-8"
        onClick={onClose}
      >
        <motion.div
          initial={{ opacity: 0, scale: 0.96, y: 12 }}
          animate={{ opacity: 1, scale: 1, y: 0 }}
          exit={{ opacity: 0, scale: 0.96, y: 12 }}
          transition={{ duration: 0.2 }}
          className="bg-card text-card-foreground shadow-xl w-full max-w-6xl h-[90vh] max-h-[800px] flex flex-col relative overflow-hidden z-50 rounded-xl"
          onClick={(e) => e.stopPropagation()}
        >
          {onBack && (
            <button
              onClick={onBack}
              className="absolute top-3 left-3 text-muted-foreground hover:text-foreground z-50 p-1 rounded-full bg-card hover:bg-muted"
              aria-label="Go back"
            >
              <ArrowLeft size={20} />
            </button>
          )}
          <button
            onClick={onClose}
            className="absolute top-3 right-3 text-muted-foreground hover:text-foreground z-50 p-1 rounded-full bg-card hover:bg-muted"
            aria-label="Close details"
          >
            <X size={20} />
          </button>

          <div className="flex flex-1 h-full overflow-hidden pt-10">
            {isLoading && (
              <div className="absolute inset-0 bg-card z-20 flex overflow-hidden pt-10 rounded-xl">
                {/* Left panel skeleton */}
                <div className="w-1/3 border-r border-border p-4 md:p-6 flex flex-col items-center gap-3">
                  <Skeleton className="w-40 h-40 rounded-full" />
                  <Skeleton className="h-6 w-36 rounded" />
                  <Skeleton className="h-6 w-24 rounded-full" />
                  <div className="w-full space-y-2 mt-2">
                    <Skeleton className="h-4 w-3/4 mx-auto rounded" />
                    <Skeleton className="h-4 w-2/3 mx-auto rounded" />
                    <Skeleton className="h-4 w-3/4 mx-auto rounded" />
                    <Skeleton className="h-4 w-1/2 mx-auto rounded" />
                  </div>
                </div>
                {/* Right panel skeleton */}
                <div className="w-2/3 p-4 md:p-6 space-y-3 overflow-hidden">
                  {Array.from({ length: 7 }).map((_, i) => (
                    <Skeleton key={i} className="h-16 w-full rounded-lg" />
                  ))}
                </div>
              </div>
            )}
            {error && !isLoading && (
              <div className="absolute inset-0 flex flex-col items-center justify-center bg-card/80 z-20 p-4 text-center">
                <p className="text-destructive font-semibold mb-2">Error</p>
                <p className="text-sm text-destructive">{error}</p>
                <button onClick={onClose} className="mt-4 text-sm underline">Close</button>
              </div>
            )}
            {!isLoading && !error && (
              <>
                {/* Left Panel: Stats */}
                <div className="w-1/3 border-r border-border p-4 md:p-6 overflow-y-auto flex flex-col items-center">
                  <a
                    href={createAthleteUrl(fighterDetails?.FIGHTER || champion.name)}
                    target="_blank"
                    rel="noopener noreferrer"
                  >
                    <div
                      className="relative w-40 h-40 mb-4 rounded-full overflow-hidden bg-card flex items-center justify-center flex-shrink-0 hover:opacity-80 transition-opacity cursor-pointer"
                    >
                      {(fighterImages?.HEADSHOT || champion.imageUrl) ? (
                        // eslint-disable-next-line @next/next/no-img-element
                        <img
                          src={fighterImages?.HEADSHOT || champion.imageUrl || ''}
                          alt={champion.name}
                          className="absolute inset-0 h-full w-full object-cover object-top"
                          onError={(e) => { (e.currentTarget as HTMLImageElement).style.display = 'none'; }}
                        />
                      ) : (
                        <User className="h-16 w-16 text-gray-500 dark:text-gray-400" />
                      )}
                    </div>
                  </a>

                  <a
                    href={createAthleteUrl(fighterDetails?.FIGHTER || champion.name)}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="hover:underline hover:text-primary transition-colors cursor-pointer"
                  >
                    <h2 className="text-2xl font-bold mb-2 text-center">
                      {fighterDetails?.FIGHTER || champion.name}
                    </h2>
                  </a>

                  {/* Current Elo — gold ring; fetched from batch or self-fetch fallback */}
                  {!hideElo && (activeEloHistory && activeEloHistory.length > 0 ? (
                    <div className="inline-flex items-center gap-1 px-2.5 py-0.5 rounded-full border border-yellow-400/60 bg-yellow-400/10 mb-3">
                      <span className="text-xs text-foreground">Elo: <span className="font-bold text-sm">{activeEloHistory[activeEloHistory.length - 1].elo}</span></span>
                    </div>
                  ) : fetchingElo ? (
                    <div className="inline-flex items-center gap-1.5 px-2.5 py-0.5 rounded-full border border-border bg-muted/30 mb-3">
                      <span className="text-xs text-muted-foreground">Elo:</span>
                      <Loader2 className="w-3 h-3 animate-spin text-muted-foreground" />
                    </div>
                  ) : null)}

                  {fighterDetails ? (
                    <div className="space-y-2 text-sm w-full text-center">
                      {fighterDetails.HEIGHT && <p><strong>Height:</strong> {fighterDetails.HEIGHT}</p>}
                      {fighterDetails.WEIGHT && <p><strong>Weight:</strong> {fighterDetails.WEIGHT}</p>}
                      {fighterDetails.REACH && <p><strong>Reach:</strong> {fighterDetails.REACH}</p>}
                      {fighterDetails.STANCE && <p><strong>Stance:</strong> {fighterDetails.STANCE}</p>}
                      {fighterDetails.DOB && (
                        <p>
                          <strong>Born:</strong> {new Date(fighterDetails.DOB).toLocaleDateString()}
                          {calculateAge(fighterDetails.DOB) !== null && ` (Age: ${calculateAge(fighterDetails.DOB)})`}
                        </p>
                      )}
                      {/* Weight classes fought in (only known/named divisions) */}
                      {(() => {
                        const wcs = Array.from(new Set(
                          fighterHistory
                            .map(f => normalizeWeightClass(f.WEIGHTCLASS))
                            .filter(k => k && isKnownDivision(k))
                        ));

                        // Sort by canonical weight-class order (lightest → heaviest)
                        wcs.sort((a, b) => {
                          const ai = WEIGHT_CLASS_ORDER.indexOf(a as any);
                          const bi = WEIGHT_CLASS_ORDER.indexOf(b as any);
                          if (ai === -1) return 1;
                          if (bi === -1) return -1;
                          return ai - bi;
                        });
                        if (!wcs.length) return null;
                        return (
                          <div className="flex flex-wrap justify-center gap-1 pt-1">
                            {wcs.map(k => (
                              <TooltipProvider key={k} delayDuration={200}>
                                <Tooltip>
                                  <TooltipTrigger asChild>
                                    <span
                                      key={k}
                                      className="text-[10px] font-semibold uppercase px-1.5 py-0.5 rounded cursor-help border"
                                      style={wcBadgeStyle(k)}
                                    >
                                      {abbreviateWeightClass(k)}
                                    </span>
                                  </TooltipTrigger>
                                  <TooltipContent side="bottom">{wcLbsLabel(k)}</TooltipContent>
                                </Tooltip>
                              </TooltipProvider>
                            ))}
                          </div>
                        );
                      })()}
                    </div>
                  ) : (
                    <p className="text-sm text-muted-foreground mt-4">No detailed stats available.</p>
                  )}
                </div>

                {/* Right Panel: Fight History */}
                <div className="w-2/3 p-4 md:p-6 overflow-y-auto">
                  {fighterHistory.length > 0 ? (
                    <ul className="space-y-1">
                      {fighterHistory.map((fight) => {
                        const isTitleFight = fight.WEIGHTCLASS?.toLowerCase().includes('title') ?? false;
                        const normalizedKey = normalizeWeightClass(fight.WEIGHTCLASS);
                        const isNonStandard = !isKnownDivision(normalizedKey);
                        const abbreviation = normalizedKey ? abbreviateWeightClass(normalizedKey) : null;

                        // Elo delta for this fight
                        const eloEntryIdx = activeEloHistory?.findIndex(h => h.event === fight.EVENT) ?? -1;
                        const eloEntry = eloEntryIdx >= 0 ? activeEloHistory![eloEntryIdx] : null;
                        const prevElo = eloEntry?.eloBeforeFight !== undefined
                          ? eloEntry.eloBeforeFight
                          : (eloEntryIdx > 0 ? activeEloHistory![eloEntryIdx - 1].elo : 1500);
                        const eloDelta = eloEntry ? eloEntry.elo - prevElo : null;

                        return (
                          <li
                            key={fight.id}
                            className={cn(
                              'relative text-sm border-b border-border p-2 rounded-lg transition-colors duration-150 overflow-visible mb-2',
                              (fight.resultForFighter === 'draw' || fight.resultForFighter === 'nc') &&
                                'bg-gray-100 dark:bg-gray-800/30 border-gray-300 dark:border-gray-600',
                              fight.resultForFighter === 'win' &&
                                'bg-green-100 dark:bg-green-900/30 border-green-300 dark:border-green-700',
                              fight.resultForFighter === 'loss' &&
                                'bg-red-100 dark:bg-red-900/30 border-red-300 dark:border-red-700',
                              fight.resultForFighter === 'unknown' && 'bg-muted/50',
                              'hover:shadow-md cursor-pointer',
                            )}
                            onClick={() => { if (fight.URL) window.open(fight.URL, '_blank', 'noopener,noreferrer'); }}
                          >
                            <TooltipProvider>
                              {isTitleFight && !fight.defenseOutcome && (
                                <Tooltip delayDuration={300}>
                                  <TooltipTrigger asChild>
                                    <Crown size={16} className="absolute -top-1 -left-2 text-yellow-500 dark:text-yellow-400 opacity-70 -rotate-45 z-10 cursor-help" />
                                  </TooltipTrigger>
                                  <TooltipContent side="top">Title Fight</TooltipContent>
                                </Tooltip>
                              )}
                              {fight.defenseOutcome === 'success' && (
                                <Tooltip delayDuration={300}>
                                  <TooltipTrigger asChild>
                                    <ShieldCheck size={18} className="absolute -top-1 -left-2 text-blue-500 dark:text-blue-400 opacity-100 -rotate-45 z-10 cursor-help" />
                                  </TooltipTrigger>
                                  <TooltipContent side="top">Successful Title Defense</TooltipContent>
                                </Tooltip>
                              )}
                              {fight.defenseOutcome === 'fail' && (
                                <Tooltip delayDuration={300}>
                                  <TooltipTrigger asChild>
                                    <ShieldOff size={18} className="absolute -top-1 -left-2 text-red-500 dark:text-red-400 opacity-100 -rotate-45 z-10 cursor-help" />
                                  </TooltipTrigger>
                                  <TooltipContent side="top">Lost Title Defense</TooltipContent>
                                </Tooltip>
                              )}
                            </TooltipProvider>

                            <div className="flex justify-between items-center mb-1 pr-1">
                              <div className="flex items-center">
                                <span className={cn(
                                  'font-semibold mr-2 px-1.5 py-0.5 rounded text-xs uppercase',
                                  fight.resultForFighter === 'win' && 'bg-green-600 text-white',
                                  fight.resultForFighter === 'loss' && 'bg-red-600 text-white',
                                  (fight.resultForFighter === 'draw' || fight.resultForFighter === 'nc') && 'bg-gray-500 text-white',
                                  fight.resultForFighter === 'unknown' && 'bg-gray-400 text-white',
                                )}>
                                  {fight.resultForFighter === 'draw' ? 'DRAW' : fight.resultForFighter === 'nc' ? 'NC' : fight.resultForFighter}
                                </span>
                                <span>vs. </span>
                                <button
                                  className="font-bold hover:underline text-primary ml-1 cursor-pointer"
                                  onClick={(e) => { e.stopPropagation(); setSelectedOpponent(fight.OPPONENT || ''); }}
                                >
                                  {fight.OPPONENT || 'Opponent Unknown'}
                                </button>
                                {eloDelta !== null ? (
                                  <span className="inline-flex items-center ml-2 gap-1" style={{ minWidth: '2rem' }}>
                                    <span className={`text-xs font-semibold ${eloDelta >= 0 ? 'text-green-600 dark:text-green-400' : 'text-red-600 dark:text-red-400'}`}>
                                      {eloDelta >= 0 ? `+${eloDelta}` : String(eloDelta)}
                                    </span>
                                    <span className="text-[10px] text-muted-foreground">({eloEntry!.elo})</span>
                                  </span>
                                ) : !hideElo && fetchingElo ? (
                                  <span className="inline-flex items-center gap-0.5 ml-2 text-xs" style={{ minWidth: '2rem' }}>
                                    <span className={fight.resultForFighter === 'win' ? 'text-green-600 dark:text-green-400' : fight.resultForFighter === 'loss' ? 'text-red-600 dark:text-red-400' : 'text-muted-foreground'}>
                                      {fight.resultForFighter === 'win' ? '+' : fight.resultForFighter === 'loss' ? '−' : ''}
                                    </span>
                                    <Loader2 className="w-2.5 h-2.5 animate-spin text-muted-foreground" />
                                  </span>
                                ) : null}
                              </div>
                              <span className="text-xs text-muted-foreground text-right max-w-[45%] leading-tight">{fight.EVENT}</span>
                            </div>

                            <div className="flex justify-between items-baseline pr-1 mt-0.5">
                              <p className="text-xs text-muted-foreground flex-1">
                                {fight.METHOD}
                                {fight.ROUND && ` | Rd: ${fight.ROUND}`}
                                {fight.TIME && ` | Time: ${fight.TIME}`}
                              </p>
                              {fight.EVENT && (eventDates[fight.EVENT] || eventDates[fight.EVENT.toLowerCase().trim()]) && (
                                <TooltipProvider>
                                  <Tooltip delayDuration={300}>
                                    <TooltipTrigger asChild>
                                      <span className="text-xs text-muted-foreground/60 text-right leading-tight flex-shrink-0 ml-2 cursor-help">
                                        {eventDates[fight.EVENT] || eventDates[fight.EVENT.toLowerCase().trim()]}
                                      </span>
                                    </TooltipTrigger>
                                    <TooltipContent side="left">
                                      {humanDuration(diffInDays(
                                        parseLocalDate(eventDates[fight.EVENT] || eventDates[fight.EVENT.toLowerCase().trim()]),
                                        new Date()
                                      ))} ago
                                    </TooltipContent>
                                  </Tooltip>
                                </TooltipProvider>
                              )}
                            </div>

                            <div className="flex justify-between items-baseline pr-1">
                              <div className="flex-1">
                                {fight.decisionDetails && (
                                  <p className="text-xs text-muted-foreground italic mt-0.5">{fight.decisionDetails}</p>
                                )}
                                {fight.DETAILS && (
                                  <p className="text-xs text-muted-foreground mt-0.5">{fight.DETAILS}</p>
                                )}
                              </div>
                              {abbreviation && (
                                <TooltipProvider>
                                  <Tooltip delayDuration={300}>
                                    <TooltipTrigger asChild>
                                      <span
                                      className={cn(
                                        "text-[10px] font-semibold uppercase px-1 rounded-[3px] cursor-help mt-0.5 flex-shrink-0 ml-2",
                                        isNonStandard ? "bg-muted text-muted-foreground border border-border" : ""
                                      )}
                                      style={!isNonStandard ? wcBadgeStyle(normalizedKey) : {}}
                                    >
                                        {abbreviation}
                                      </span>
                                    </TooltipTrigger>
                                    <TooltipContent side="left">
                                      {wcLbsLabel(normalizedKey)}
                                    </TooltipContent>
                                  </Tooltip>
                                </TooltipProvider>
                              )}
                            </div>
                          </li>
                        );
                      })}
                    </ul>
                  ) : (
                    <p className="text-sm text-muted-foreground">No fight history found or data unavailable.</p>
                  )}
                </div>
              </>
            )}
          </div>
        </motion.div>
      </motion.div>
    </AnimatePresence>
  );
}
