"use client";

import React, {
  useState,
  useEffect,
  useMemo,
  useRef,
  useCallback,
} from "react";
import ReactECharts from "echarts-for-react";
import * as echarts from "echarts/core";
import type { EChartsOption } from "echarts";
import { LineChart, BarChart } from "echarts/charts";
import {
  TooltipComponent,
  GridComponent,
  MarkLineComponent,
} from "echarts/components";
import { CanvasRenderer } from "echarts/renderers";
import { Card, CardContent } from "@/components/ui/card";
import {
  ChevronUp,
  ChevronDown,
  ChevronsUpDown,
  BarChart2,
  Table2,
  LayoutGrid,
  Search,
  ArrowLeft,
  User,
  Loader2,
  Activity,
} from "lucide-react";
import { cn } from "@/lib/utils";
import { supabase } from "@/lib/supabaseClient";
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from "@/components/ui/tooltip";
import { Skeleton } from "@/components/ui/skeleton";
import { useTheme } from "next-themes";
import {
  FighterDetailsModal,
  type ChampionEntry,
} from "@/components/FighterDetailsModal";
import {
  BinFightersModal,
  type BinFighter,
} from "@/components/BinFightersModal";
import { AnimatePresence } from "framer-motion";
import { wcBadgeStyle } from "@/lib/wc-colors";

// Register ECharts components
echarts.use([
  TooltipComponent,
  GridComponent,
  MarkLineComponent,
  LineChart,
  BarChart,
  CanvasRenderer,
]);

/* ------------------------------------------------------------------ */
/*  TYPES                                                             */
/* ------------------------------------------------------------------ */
interface FighterHistoryPoint {
  date: string;
  elo: number;
  eloBeforeFight?: number; // post-decay, pre-fight — for accurate delta
  event: string;
  opponent: string;
  result: string;
  weightClass: string;
  url: string;
}

interface FighterSummary {
  name: string;
  currentElo: number;
  peakElo: number;
  peakDate: string;
  fights: number;
  primaryWeightClass: string;
  history: FighterHistoryPoint[];
}

interface EloData {
  candidatePool: FighterSummary[];
  byWeightClass: Record<string, FighterSummary[]>;
  wcDisplayNames: Record<string, string>;
  allFighters?: {
    name: string;
    currentElo: number;
    primaryWeightClass: string;
    lastFightDate: string;
  }[];
}

type SortCol =
  | "currentElo"
  | "peakElo"
  | "peakDate"
  | "fights"
  | "l5Wins"
  | "streak"
  | "wl";
type SortDir = "asc" | "desc";

/* ------------------------------------------------------------------ */
/*  COLOUR PALETTE                                                    */
/* ------------------------------------------------------------------ */
const PALETTE = [
  "#5470c6",
  "#91cc75",
  "#fac858",
  "#ee6666",
  "#73c0de",
  "#3ba272",
  "#fc8452",
  "#9a60b4",
  "#ea7ccc",
  "#ff9f7f",
  "#67e0e3",
  "#e690d1",
  "#e7bcf3",
  "#8378ea",
  "#96dee8",
];

const DIST_BIN = 20;

/* ------------------------------------------------------------------ */
/*  DIVISION ORDER & LABELS                                          */
/* ------------------------------------------------------------------ */
const DIVISION_ORDER = [
  "heavyweight",
  "lightheavyweight",
  "middleweight",
  "welterweight",
  "lightweight",
  "featherweight",
  "bantamweight",
  "flyweight",
  "women_featherweight",
  "women_bantamweight",
  "women_flyweight",
  "women_strawweight",
];

interface SearchResult {
  name: string;
  headshot?: string | null;
  currentElo?: number;
  peakElo?: number;
  record?: string;
  weightClass?: string;
  fights?: number;
  lastFightDate?: string;
  eloStatus: "loaded" | "loading" | "not_found";
}

const WC_ABBREV: Record<string, string> = {
  heavyweight: "HW",
  lightheavyweight: "LHW",
  middleweight: "MW",
  welterweight: "WW",
  lightweight: "LW",
  featherweight: "FEW",
  bantamweight: "BW",
  flyweight: "FLW",
  women_bantamweight: "W-BW",
  women_featherweight: "W-FEW",
  women_flyweight: "W-FLW",
  women_strawweight: "W-SW",
};

const TWO_YEARS_MS = 2 * 365.25 * 24 * 60 * 60 * 1000;

/* ------------------------------------------------------------------ */
/*  HELPERS                                                          */
/* ------------------------------------------------------------------ */
function getUFCEventUrl(eventName: string): string {
  const slug = eventName
    .toLowerCase()
    .normalize("NFD")
    .replace(/[\u0300-\u036f]/g, "")
    .replace(/[^a-z0-9\s]/g, "")
    .trim()
    .replace(/\s+/g, "-");
  return `https://www.ufc.com/event/${slug}`;
}

const WC_FULL_NAME: Record<string, string> = {
  heavyweight: "Heavyweight (265 lbs)",
  lightheavyweight: "Light Heavyweight (205 lbs)",
  middleweight: "Middleweight (185 lbs)",
  welterweight: "Welterweight (170 lbs)",
  lightweight: "Lightweight (155 lbs)",
  featherweight: "Featherweight (145 lbs)",
  bantamweight: "Bantamweight (135 lbs)",
  flyweight: "Flyweight (125 lbs)",
  women_bantamweight: "Women's Bantamweight (135 lbs)",
  women_flyweight: "Women's Flyweight (125 lbs)",
  women_strawweight: "Women's Strawweight (115 lbs)",
  women_featherweight: "Women's Featherweight (145 lbs)",
};

function computeFighterL5Streak(f: FighterSummary): {
  l5W: number;
  l5L: number;
  l5Str: string;
  streak: number;
} {
  const sorted = [...f.history].sort(
    (a, b) => new Date(a.date).getTime() - new Date(b.date).getTime(),
  );
  const last5 = sorted.slice(-5);
  const l5W = last5.filter((h) => h.result === "W").length;
  const l5L = last5.filter((h) => h.result === "L").length;
  const l5D = last5.filter((h) => h.result === "D").length;
  const l5Str = l5D > 0 ? `${l5W}-${l5L}-${l5D}` : `${l5W}-${l5L}`;
  let streak = 0;
  for (let i = sorted.length - 1; i >= 0; i--) {
    const r = sorted[i].result;
    if (streak === 0) {
      if (r === "W") streak = 1;
      else if (r === "L") streak = -1;
      else break;
    } else if (streak > 0 && r === "W") streak++;
    else if (streak < 0 && r === "L") streak--;
    else break;
  }
  return { l5W, l5L, l5Str, streak };
}

function getAthleteUrl(name: string): string {
  const slug = name
    .normalize("NFD")
    .replace(/[\u0300-\u036f]/g, "")
    .replace(/\./g, "")
    .trim()
    .toLowerCase()
    .replace(/\s+/g, "-");
  return `https://www.ufc.com/athlete/${slug}`;
}

/* Gaussian KDE helper */
function gaussianKDE(
  values: number[],
  bandwidth: number,
): (x: number) => number {
  return (x: number) => {
    let sum = 0;
    for (const v of values) {
      const u = (x - v) / bandwidth;
      sum += Math.exp(-0.5 * u * u);
    }
    return sum / (values.length * bandwidth * Math.sqrt(2 * Math.PI));
  };
}

function formatDob(dob: string): string {
  const d = new Date(dob);
  if (isNaN(d.getTime())) return dob;
  const month = d.getUTCMonth() + 1;
  const day = d.getUTCDate();
  const year = d.getUTCFullYear();
  const today = new Date();
  let age = today.getFullYear() - year;
  if (
    today.getMonth() < d.getUTCMonth() ||
    (today.getMonth() === d.getUTCMonth() && today.getDate() < day)
  )
    age--;
  return `${month}/${day}/${year} (Age: ${age})`;
}

/* ------------------------------------------------------------------ */
/*  COMPONENT                                                         */
/* ------------------------------------------------------------------ */
export default function UfcEloPage() {
  const [data, setData] = useState<EloData | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [selectedDivision, setSelectedDivision] = useState<string>("overall");
  const [sortCol, setSortCol] = useState<SortCol>("currentElo");
  const [sortDir, setSortDir] = useState<SortDir>("desc");
  const [lockedFighter, setLockedFighter] = useState<string | null>(null);
  const [selectedChampion, setSelectedChampion] =
    useState<ChampionEntry | null>(null);
  const [viewMode, setViewMode] = useState<
    "both" | "chart" | "table" | "distribution"
  >("both");
  const [activeOnly, setActiveOnly] = useState(false);
  const [searchQuery, setSearchQuery] = useState("");
  const [searchResults, setSearchResults] = useState<SearchResult[]>([]);
  const [isSearching, setIsSearching] = useState(false);
  const searchIdRef = useRef(0);
  const [distBinModal, setDistBinModal] = useState<{
    lo: number;
    hi: number;
    fighters: BinFighter[];
  } | null>(null);
  // Set of normalized non-interim champion names → name as it appears in the lineage
  const [championMap, setChampionMap] = useState<Map<string, ChampionEntry>>(
    new Map(),
  );
  const [candidateTott, setCandidateTott] = useState<
    Map<
      string,
      {
        height?: string | null;
        weight?: string | null;
        reach?: string | null;
        stance?: string | null;
        dob?: string | null;
      }
    >
  >(new Map());
  const [candidateHeadshots, setCandidateHeadshots] = useState<
    Map<string, string>
  >(new Map());
  const [tableRowHeight, setTableRowHeight] = useState<number | null>(null);
  const tableWrapperRef = useRef<HTMLDivElement>(null);
  const { resolvedTheme } = useTheme();

  // Refs for baseline-hover cursor tracking
  const chartRef = useRef<ReactECharts>(null);
  const chartContainerRef = useRef<HTMLDivElement>(null);
  const cursorDataYRef = useRef<number | null>(null);

  const handleChartMouseMove = useCallback(
    (e: React.MouseEvent<HTMLDivElement>) => {
      if (!chartRef.current || !chartContainerRef.current) return;
      try {
        const instance = (chartRef.current as any).getEchartsInstance();
        const rect = chartContainerRef.current.getBoundingClientRect();
        const pixelY = e.clientY - rect.top;
        const dataY = instance.convertFromPixel({ yAxisIndex: 0 }, pixelY);
        cursorDataYRef.current = typeof dataY === "number" ? dataY : null;
      } catch {
        cursorDataYRef.current = null;
      }
    },
    [],
  );

  const handleChartMouseLeave = useCallback(() => {
    cursorDataYRef.current = null;
  }, []);

  // Attach zrender handler once chart is ready so empty-area clicks clear isolation
  const handleChartReady = useCallback((instance: any) => {
    try {
      instance.getZr().on("click", (params: any) => {
        if (!params.target) {
          setLockedFighter(null);
        }
      });
    } catch {
      /* ignore */
    }
  }, []);

  // Fetch champion lineage once to determine who gets the popup
  useEffect(() => {
    fetch("/api/ufc-champions/lineage")
      .then((r) => r.json())
      .then((json) => {
        if (!json.lineage) return;
        const map = new Map<string, ChampionEntry>();
        for (const [, wc] of Object.entries(
          json.lineage as Record<string, { champions: ChampionEntry[] }>,
        )) {
          for (const champ of wc.champions) {
            const isInterim =
              champ.notes?.toLowerCase().includes("interim") ?? false;
            if (isInterim) continue;
            // Normalize: lowercase, strip accents/punctuation
            const key = champ.name
              .normalize("NFD")
              .replace(/[\u0300-\u036f]/g, "")
              .replace(/[.'\u2011]/g, "")
              .toLowerCase()
              .trim();
            if (!map.has(key)) map.set(key, champ);
          }
        }
        setChampionMap(map);
      })
      .catch(() => {
        /* non-critical */
      });
  }, []);

  // Pre-fetch TOTT physical stats + headshots for ALL displayed fighters (candidatePool + all byWeightClass)
  useEffect(() => {
    if (!data) return;
    const nk = (s: string) =>
      s
        .normalize("NFD")
        .replace(/[\u0300-\u036f]/g, "")
        .replace(/[.'\u2011]/g, "")
        .toLowerCase()
        .trim();
    const stripAccents = (s: string) =>
      s
        .normalize("NFD")
        .replace(/[\u0300-\u036f]/g, "")
        .trim();
    // Merge candidatePool + all byWeightClass fighters into a unique set
    const allPoolFighters = [
      ...(data.candidatePool || []),
      ...Object.values(data.byWeightClass || {}).flat(),
    ];
    const names = [...new Set(allPoolFighters.map((f) => f.name))];
    if (names.length === 0) return;
    // Fetch ALL TOTT rows via server-side API that paginates past PostgREST's 1000-row default cap
    fetch("/api/ufc-tott")
      .then((r) => (r.ok ? r.json() : null))
      .then((json) => {
        if (!json?.data) return;
        const map = new Map<
          string,
          {
            height?: string | null;
            weight?: string | null;
            reach?: string | null;
            stance?: string | null;
            dob?: string | null;
          }
        >();
        for (const r of json.data as any[]) {
          map.set(nk(r.FIGHTER), {
            height: r.HEIGHT,
            weight: r.WEIGHT,
            reach: r.REACH,
            stance: r.STANCE,
            dob: r.DOB,
          });
        }
        setCandidateTott(map);
      });
    // Fetch ALL champion images at once (table is small, ~200-400 rows) — avoids .in() truncation for large fighter pools
    supabase
      .from("ufc_champion_images")
      .select("FIRST, LAST, HEADSHOT")
      .limit(2000)
      .then(({ data: rows }) => {
        if (!rows) return;
        const map = new Map<string, string>();
        for (const r of rows as any[]) {
          if (r.HEADSHOT) map.set(nk(`${r.FIRST} ${r.LAST}`), r.HEADSHOT);
        }
        setCandidateHeadshots(map);
      });
  }, [data]);

  // Lock page scroll in single-view modes
  useEffect(() => {
    if (viewMode !== "both") {
      document.documentElement.style.overflow = "hidden";
      document.body.style.overflow = "hidden";
      // Also lock the layout's scrollable main container
      const layoutMain =
        document.querySelector<HTMLElement>("body > div > main");
      if (layoutMain) layoutMain.style.overflow = "hidden";
    } else {
      document.documentElement.style.overflow = "";
      document.body.style.overflow = "";
      const layoutMain =
        document.querySelector<HTMLElement>("body > div > main");
      if (layoutMain) layoutMain.style.overflow = "";
    }
    return () => {
      document.documentElement.style.overflow = "";
      document.body.style.overflow = "";
      const layoutMain =
        document.querySelector<HTMLElement>("body > div > main");
      if (layoutMain) layoutMain.style.overflow = "";
    };
  }, [viewMode]);

  // Reset lock when division changes
  useEffect(() => {
    setLockedFighter(null);
  }, [selectedDivision]);

  useEffect(() => {
    fetch("/api/ufc-elo")
      .then((r) => r.json())
      .then((json) => {
        if (json.error) throw new Error(json.error);
        setData(json);
      })
      .catch((e) => setError(e.message))
      .finally(() => setLoading(false));
  }, []);

  const isDark = resolvedTheme === "dark";
  const textColor = isDark ? "#e5e7eb" : "#374151";

  // Only real divisions, in canonical order
  const divisions = useMemo(() => {
    if (!data) return [];
    return Object.keys(data.byWeightClass)
      .filter((k) => DIVISION_ORDER.includes(k))
      .sort((a, b) => DIVISION_ORDER.indexOf(a) - DIVISION_ORDER.indexOf(b));
  }, [data]);

  // Full candidate pool for selected division
  const candidatePool = useMemo(() => {
    if (!data) return [];
    const base =
      selectedDivision === "overall"
        ? data.candidatePool
        : data.byWeightClass[selectedDivision] || [];
    if (!activeOnly) return base;
    const cutoff = Date.now() - TWO_YEARS_MS;
    return base.filter((f) => {
      const last = f.history[f.history.length - 1];
      return last && new Date(last.date).getTime() >= cutoff;
    });
  }, [data, selectedDivision, activeOnly]);

  const maxDisplay = selectedDivision === "overall" ? 20 : 10;

  // Flat name→history map across all divisions for opponent elo lookup
  const eloHistoryMap = useMemo(() => {
    if (!data) return new Map<string, FighterHistoryPoint[]>();
    const map = new Map<string, FighterHistoryPoint[]>();
    const add = (f: FighterSummary) => {
      const key = f.name
        .normalize("NFD")
        .replace(/[\u0300-\u036f]/g, "")
        .replace(/[.'\u2011]/g, "")
        .toLowerCase()
        .trim();
      if (!map.has(key)) map.set(key, f.history);
    };
    data.candidatePool.forEach(add);
    Object.values(data.byWeightClass).flat().forEach(add);
    return map;
  }, [data]);

  const eloLookup = useCallback(
    (name: string) => {
      const key = name
        .normalize("NFD")
        .replace(/[\u0300-\u036f]/g, "")
        .replace(/[.'\u2011]/g, "")
        .toLowerCase()
        .trim();
      return eloHistoryMap.get(key);
    },
    [eloHistoryMap],
  );

  // Fighter search
  useEffect(() => {
    const q = searchQuery.trim();
    if (!q) {
      setSearchResults([]);
      setIsSearching(false);
      return;
    }
    // Set searching immediately to show spinner during the debounce delay,
    // preventing "No fighters found" from flashing before results arrive.
    setIsSearching(true);
    const timer = setTimeout(async () => {
      const searchId = ++searchIdRef.current;
      try {
        const normKey = (s: string) =>
          s
            .normalize("NFD")
            .replace(/[\u0300-\u036f]/g, "")
            .replace(/[.'\u2011]/g, "")
            .toLowerCase()
            .trim();
        // Run prefix query (names starting with q) and contains query in parallel
        const normQ = normKey(q);
        const [{ data: prefixData }, { data: containsData }] =
          await Promise.all([
            supabase
              .from("ufc_fighter_tott")
              .select("FIGHTER")
              .ilike("FIGHTER", `${q}%`)
              .limit(100),
            supabase
              .from("ufc_fighter_tott")
              .select("FIGHTER")
              .ilike("FIGHTER", `%${q}%`)
              .limit(300),
          ]);
        if (searchIdRef.current !== searchId) return;

        // Merge and deduplicate
        const seen = new Set<string>();
        const allData: Array<{ FIGHTER: string }> = [];
        for (const row of [...(prefixData ?? []), ...(containsData ?? [])]) {
          const n = normKey((row as { FIGHTER: string }).FIGHTER);
          if (!seen.has(n)) {
            seen.add(n);
            allData.push(row as { FIGHTER: string });
          }
        }
        if (allData.length === 0) {
          setSearchResults([]);
          return;
        }

        // Score: full name starts with query (300) > first word starts with (200) > any word starts with (100) > contains anywhere (10)
        const scored = allData.map((f) => {
          const normName = normKey(f.FIGHTER);
          const words = normName.split(" ");
          let score = 0;
          if (normName.startsWith(normQ)) score = 300;
          else if (words[0]?.startsWith(normQ)) score = 200;
          else if (words.some((w: string) => w.startsWith(normQ))) score = 100;
          else score = 10;
          return { ...f, _score: score };
        });
        scored.sort(
          (a: any, b: any) =>
            b._score - a._score || a.FIGHTER.localeCompare(b.FIGHTER),
        );
        const fightersData = scored.slice(0, 14);

        // Fetch headshots for matched fighters — include accent-stripped variants of first words
        const stripA = (s: string) =>
          s
            .normalize("NFD")
            .replace(/[\u0300-\u036f]/g, "")
            .trim();
        const firstWords = [
          ...new Set(
            fightersData.slice(0, 50).flatMap((f: any) => {
              const word = f.FIGHTER.split(" ")[0];
              const stripped = stripA(word);
              return stripped !== word ? [word, stripped] : [word];
            }),
          ),
        ];
        const { data: imagesData } = await supabase
          .from("ufc_champion_images")
          .select("FIRST, LAST, HEADSHOT")
          .in("FIRST", firstWords)
          .limit(60);
        if (searchIdRef.current !== searchId) return;
        const imageMap = new Map<string, string>();
        if (imagesData) {
          for (const img of imagesData as any[]) {
            if (img.HEADSHOT)
              imageMap.set(normKey(`${img.FIRST} ${img.LAST}`), img.HEADSHOT);
          }
        }

        const buildResult = (
          name: string,
          hist: typeof eloHistoryMap extends Map<string, infer V> ? V : never,
        ): SearchResult => {
          const wins = hist.filter((h) => h.result === "W").length;
          const losses = hist.filter((h) => h.result === "L").length;
          const draws = hist.filter((h) => h.result === "D").length;
          const peakElo = Math.max(...hist.map((h) => h.elo));
          const wcCounts = new Map<string, number>();
          hist.forEach((h) => {
            if (h.weightClass)
              wcCounts.set(
                h.weightClass,
                (wcCounts.get(h.weightClass) ?? 0) + 1,
              );
          });
          const topWc = [...wcCounts.entries()].sort(
            (a, b) => b[1] - a[1],
          )[0]?.[0];
          return {
            name,
            headshot: imageMap.get(normKey(name)) ?? null,
            currentElo: hist[hist.length - 1]?.elo,
            peakElo,
            record:
              draws > 0 ? `${wins}-${losses}-${draws}` : `${wins}-${losses}`,
            weightClass: topWc,
            fights: hist.length,
            lastFightDate: hist[hist.length - 1]?.date,
            eloStatus: "loaded",
          };
        };

        // Show all fighters: those in batch get elo immediately, others get loading state
        const initial: SearchResult[] = fightersData.map((f: any) => {
          const key = normKey(f.FIGHTER);
          const hist = eloHistoryMap.get(key);
          if (hist && hist.length > 0) return buildResult(f.FIGHTER, hist);
          return {
            name: f.FIGHTER,
            headshot: imageMap.get(key) ?? null,
            eloStatus: "loading" as const,
          };
        });
        if (searchIdRef.current !== searchId) return;
        setSearchResults(initial);

        // Fetch elo for fighters not in the pre-loaded batch
        const toFetch = initial.filter((r) => r.eloStatus === "loading");
        if (toFetch.length > 0) {
          const promises = toFetch.map(async (r) => {
            try {
              const res = await fetch(
                `/api/ufc-elo?name=${encodeURIComponent(r.name)}`,
              );
              if (!res.ok)
                return { name: r.name, status: "not_found" as const };
              const d = await res.json();
              const hist: Array<{
                date: string;
                elo: number;
                event: string;
                result: string;
                weightClass: string;
              }> = d.history || [];
              if (hist.length === 0)
                return { name: r.name, status: "not_found" as const };
              const wins = hist.filter((h) => h.result === "W").length;
              const losses = hist.filter((h) => h.result === "L").length;
              const draws = hist.filter((h) => h.result === "D").length;
              const peakElo = Math.max(...hist.map((h) => h.elo));
              const wcCounts = new Map<string, number>();
              hist.forEach((h) => {
                if (h.weightClass)
                  wcCounts.set(
                    h.weightClass,
                    (wcCounts.get(h.weightClass) ?? 0) + 1,
                  );
              });
              const topWc = [...wcCounts.entries()].sort(
                (a, b) => b[1] - a[1],
              )[0]?.[0];
              return {
                name: r.name,
                status: "loaded" as const,
                currentElo: hist[hist.length - 1]?.elo,
                peakElo,
                record:
                  draws > 0
                    ? `${wins}-${losses}-${draws}`
                    : `${wins}-${losses}`,
                weightClass: topWc,
                fights: hist.length,
                lastFightDate: hist[hist.length - 1]?.date,
              };
            } catch {
              return { name: r.name, status: "not_found" as const };
            }
          });
          for (const p of promises) {
            p.then((result) => {
              if (searchIdRef.current !== searchId) return;
              setSearchResults((prev) =>
                prev.map((r) => {
                  if (normKey(r.name) !== normKey(result.name)) return r;
                  if (result.status === "not_found")
                    return { ...r, eloStatus: "not_found" as const };
                  return {
                    ...r,
                    eloStatus: "loaded" as const,
                    currentElo: result.currentElo,
                    peakElo: result.peakElo,
                    record: result.record,
                    weightClass: result.weightClass,
                    fights: result.fights,
                    lastFightDate: result.lastFightDate,
                  };
                }),
              );
            });
          }
        }
      } catch (e) {
        console.error("Fighter search error:", e);
      } finally {
        if (searchIdRef.current === searchId) setIsSearching(false);
      }
    }, 350);
    return () => clearTimeout(timer);
  }, [searchQuery, eloHistoryMap]);

  // Enrich the pool with precomputed L5/streak so we can sort by them
  const enrichedPool = useMemo(
    () => candidatePool.map((f) => ({ ...f, ...computeFighterL5Streak(f) })),
    [candidatePool],
  );

  // Sort pool and slice to display set — drives BOTH chart and table
  const displayFighters = useMemo(() => {
    const sorted = [...enrichedPool].sort((a, b) => {
      let vA: number, vB: number;
      if (sortCol === "currentElo") {
        vA = a.currentElo;
        vB = b.currentElo;
      } else if (sortCol === "peakElo") {
        vA = a.peakElo;
        vB = b.peakElo;
      } else if (sortCol === "peakDate") {
        vA = a.peakDate ? new Date(a.peakDate).getTime() : 0;
        vB = b.peakDate ? new Date(b.peakDate).getTime() : 0;
      } else if (sortCol === "l5Wins") {
        vA = a.l5W;
        vB = b.l5W;
      } else if (sortCol === "streak") {
        vA = a.streak;
        vB = b.streak;
      } else if (sortCol === "wl") {
        const getRatio = (f: { history: FighterHistoryPoint[] }) => {
          const w = f.history.filter((h) => h.result === "W").length;
          const l = f.history.filter((h) => h.result === "L").length;
          return l === 0 ? (w > 0 ? 1e9 : 0) : w / l;
        };
        vA = getRatio(a);
        vB = getRatio(b);
      } else {
        vA = a.fights;
        vB = b.fights;
      }
      return sortDir === "desc" ? vB - vA : vA - vB;
    });
    return sorted.slice(0, maxDisplay);
  }, [enrichedPool, sortCol, sortDir, maxDisplay]);

  // Measure table container to compute exact row heights that fill the viewport
  useEffect(() => {
    if (viewMode !== "table") return;
    const container = tableWrapperRef.current;
    if (!container) return;
    const update = () => {
      const thead = container.querySelector("thead") as HTMLElement | null;
      const theadH = thead ? thead.offsetHeight : 46;
      const borderTotal = Math.max(0, displayFighters.length - 1);
      const available = container.offsetHeight - theadH - borderTotal;
      setTableRowHeight(
        Math.max(28, Math.floor(available / displayFighters.length)),
      );
    };
    update();
    const ro = new ResizeObserver(update);
    ro.observe(container);
    return () => ro.disconnect();
  }, [viewMode, displayFighters.length]);

  const handleSort = (col: SortCol) => {
    if (sortCol === col) setSortDir((d) => (d === "desc" ? "asc" : "desc"));
    else {
      setSortCol(col);
      setSortDir(col === "peakDate" ? "asc" : "desc");
    }
  };

  // Per-series fight metadata for click handler [seriesIndex][dataIndex]
  const pointMeta = useMemo(
    () =>
      displayFighters.map((fighter) =>
        fighter.history
          .filter((h) => h.date && !isNaN(new Date(h.date).getTime()))
          .map((h) => ({
            event: h.event,
            opponent: h.opponent,
            result: h.result,
            date: h.date,
            url: h.url,
          })),
      ),
    [displayFighters],
  );

  const chartEvents = useMemo(
    () => ({
      click: (params: any) => {
        if (
          params.componentType === "series" &&
          params.componentSubType === "line"
        ) {
          // Toggle highlight/lock — setTimeout avoids getRawIndex errors during notMerge
          if (params.seriesName && params.seriesName !== "__baseline__") {
            const sn = params.seriesName;
            setTimeout(
              () => setLockedFighter((prev) => (prev === sn ? null : sn)),
              0,
            );
          }

          // Open UFC Stats link if available
          try {
            if (params.data && Array.isArray(params.data) && params.data[0]) {
              const clickTs = Number(params.data[0]);
              const seriesIdx = params.seriesIndex;
              const fights = pointMeta[seriesIdx];
              if (fights) {
                // Find fight exactly matching timestamp (ignore synthetic 1ms points)
                const fight = fights.find(
                  (f) => new Date(f.date).getTime() === clickTs,
                );
                if (fight?.url) {
                  window.open(fight.url, "_blank", "noopener,noreferrer");
                }
              }
            }
          } catch (e) {
            console.error("Chart click error:", e);
          }
        } else {
          // Click on empty chart area → clear isolation
          setLockedFighter(null);
        }
      },
    }),
    [pointMeta],
  );

  const SortArrow = ({ col }: { col: SortCol }) => {
    if (sortCol !== col)
      return <ChevronsUpDown className="inline w-3 h-3 ml-0.5 opacity-40" />;
    return sortDir === "desc" ? (
      <ChevronDown className="inline w-3 h-3 ml-0.5" />
    ) : (
      <ChevronUp className="inline w-3 h-3 ml-0.5" />
    );
  };

  // Build ECharts option — synced to displayFighters (same set as table)
  const chartOption: EChartsOption = useMemo(() => {
    if (!displayFighters.length) return {};

    const todayTs = Date.now();

    // Round y-axis min down to nearest 50 so axis starts clean (no 0)
    const allElos = displayFighters.flatMap((f) =>
      f.history.filter((h) => h.date).map((h) => h.elo),
    );
    const rawMin = allElos.length ? Math.min(...allElos) : 1000;
    const yMin = Math.floor((rawMin - 30) / 50) * 50;

    const series = displayFighters.map((fighter, i) => {
      const color = PALETTE[i % PALETTE.length];
      const isActive = lockedFighter === null || lockedFighter === fighter.name;
      const opacity = isActive ? 1 : 0.1;
      const points = fighter.history
        .filter((h) => h.date)
        .map((h) => {
          const ts = new Date(h.date).getTime();
          return isNaN(ts) ? null : [ts, h.elo];
        })
        .filter(Boolean) as [number, number][];

      // All timestamps we want to hide symbols for
      const originalLastTs =
        points.length > 0 ? points[points.length - 1][0] : 0;
      const ghostTs = points.length ? points[0][0] - 1 : null;
      const extendTs = points.length ? todayTs : null; // always extend to today

      // Build uniform [ts, elo] array — no mixed object/array format (prevents getRawIndex errors)
      const dataWithStart: [number, number][] = points.length
        ? [
            [points[0][0] - 1, 1500] as [number, number],
            ...points,
            [todayTs, points[points.length - 1][1]] as [number, number], // always extend to today
          ]
        : [];

      return {
        name: fighter.name,
        type: "line" as const,
        step: "end" as const,
        showSymbol: true,
        symbol: "circle",
        symbolSize: (value: number[]) => {
          const ts = value[0];
          if (ghostTs !== null && ts === ghostTs) return 0;
          if (extendTs !== null && ts === extendTs) return 0;
          return 5;
        },
        lineStyle: { width: isActive ? 1.5 : 0.8, color, opacity },
        itemStyle: { color, opacity },
        emphasis: { disabled: true },
        selectedMode: false,
        data: dataWithStart,
      };
    });

    // Add 1500 baseline as a proper series — avoids getRawIndex errors caused by markLine
    const allHistoryTs = displayFighters.flatMap((f) =>
      f.history
        .filter((h) => h.date && !isNaN(new Date(h.date).getTime()))
        .map((h) => new Date(h.date).getTime()),
    );
    const baselineMinTs = allHistoryTs.length
      ? Math.min(...allHistoryTs) - 1
      : todayTs - 10 * 365.25 * 24 * 60 * 60 * 1000;
    const baselineSeries = {
      name: "__baseline__",
      type: "line" as const,
      data: [
        [baselineMinTs, 1500],
        [todayTs, 1500],
      ] as [number, number][],
      lineStyle: {
        color: "rgba(239,68,68,0.4)",
        type: "dashed" as const,
        width: 1,
      },
      symbol: "none" as const,
      symbolSize: 0,
      silent: true,
      animation: false,
      z: -1,
      legendHoverLink: false,
      tooltip: { show: false },
      endLabel: {
        show: true,
        formatter: () => "1500",
        color: "rgba(239,68,68,0.65)",
        fontSize: 10,
        fontFamily: "Merriweather Sans",
      },
    };

    // Snapshot of history for the tooltip closure — pre-compute ts and prevElo for delta display
    const historySnapshot = displayFighters.map((f) => {
      const history = f.history
        .filter((h) => h.date && !isNaN(new Date(h.date).getTime()))
        .map((h, j, arr) => ({
          ts: new Date(h.date).getTime(),
          elo: h.elo,
          prevElo:
            h.eloBeforeFight !== undefined
              ? h.eloBeforeFight
              : j > 0
                ? arr[j - 1].elo
                : 1500,
          event: h.event,
        }));
      const originalLastTs =
        history.length > 0 ? history[history.length - 1].ts : 0;
      return { name: f.name, history, originalLastTs, isExtended: true }; // always extend to today
    });

    return {
      backgroundColor: "transparent",
      animation: true,
      animationDuration: 400,
      animationEasing: "cubicOut" as const,
      legend: { show: false },
      tooltip: {
        trigger: "axis",
        axisPointer: {
          type: "line",
          lineStyle: {
            color: isDark ? "rgba(255,255,255,0.3)" : "rgba(0,0,0,0.25)",
            type: "dashed" as const,
            width: 1,
          },
        },
        backgroundColor: isDark ? "#1f2937" : "#ffffff",
        borderColor: isDark ? "#374151" : "#d1d5db",
        textStyle: {
          color: textColor,
          fontSize: 11,
          fontFamily: "Merriweather Sans",
        },
        confine: true,
        formatter: (params: any) => {
          // Show baseline description when cursor is directly on the 1500 line (tight ±8 tolerance)
          const cursorY = cursorDataYRef.current;
          if (cursorY !== null && Math.abs(cursorY - 1500) <= 8) {
            return `<div style="padding:6px 8px;font-family:'Merriweather Sans',sans-serif;font-size:11px">Starting baseline — all fighters begin at 1500 Elo</div>`;
          }

          const cursorTs = Number(
            Array.isArray(params)
              ? (params[0]?.axisValue ?? params[0]?.value?.[0])
              : (params?.axisValue ?? params?.value?.[0]),
          );
          if (!cursorTs || isNaN(cursorTs)) return "";

          const isGhost = historySnapshot.some(({ history }) =>
            history.some((h) => h.ts - 1 === cursorTs),
          );
          if (isGhost) return null as any;

          const date = new Date(cursorTs).toLocaleDateString("en-US", {
            month: "short",
            day: "numeric",
            year: "numeric",
          });

          const rows: {
            name: string;
            elo: number;
            color: string;
            delta: number | null;
            fightEvent?: string;
          }[] = [];
          historySnapshot.forEach(
            ({ name, history, originalLastTs, isExtended }, i) => {
              if (lockedFighter !== null && name !== lockedFighter) return;
              let lastElo: number | null = null;
              let delta: number | null = null;
              let fightEvent: string | undefined;
              for (let j = 0; j < history.length; j++) {
                const h = history[j];
                if (h.ts <= cursorTs) {
                  lastElo = h.elo;
                  if (h.ts === cursorTs) {
                    delta = h.prevElo !== null ? h.elo - h.prevElo : null;
                    fightEvent = h.event;
                  } else {
                    delta = null;
                    fightEvent = undefined;
                  }
                } else break;
              }
              if (lastElo === null) return;
              rows.push({
                name,
                elo: lastElo,
                color: PALETTE[i % PALETTE.length],
                delta,
                fightEvent,
              });
            },
          );

          rows.sort((a, b) => b.elo - a.elo);
          if (!rows.length) return "";

          const eventName = rows.find((r) => r.fightEvent)?.fightEvent || "";

          const lines = rows.map((r) => {
            const deltaHtml =
              r.delta !== null
                ? r.delta >= 0
                  ? ` <span style="color:#22c55e;font-weight:600">(+${r.delta})</span>`
                  : ` <span style="color:#ef4444;font-weight:600">(${r.delta})</span>`
                : "";
            return (
              `<div style="display:flex;align-items:center;gap:5px;padding:1px 0">` +
              `<span style="width:8px;height:8px;border-radius:50%;background:${r.color};display:inline-block;flex-shrink:0"></span>` +
              `<span style="flex:1">${r.name}</span>` +
              `<span style="font-weight:600;padding-left:10px">${r.elo}</span>${deltaHtml}</div>`
            );
          });
          const headerHtml = eventName
            ? `<div style="font-size:10px;color:${textColor};margin-bottom:2px;text-transform:uppercase;letter-spacing:0.05em;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;max-width:260px">${eventName}</div>` +
              `<div style="font-weight:700;margin-bottom:6px;border-bottom:1px solid ${isDark ? "#374151" : "#e5e7eb"};padding-bottom:4px">${date}</div>`
            : `<div style="font-weight:700;margin-bottom:6px;border-bottom:1px solid ${isDark ? "#374151" : "#e5e7eb"};padding-bottom:4px">${date}</div>`;
          return (
            `<div style="padding:4px;font-size:11px;min-width:190px;font-family:'Merriweather Sans',sans-serif">` +
            headerHtml +
            lines.join("") +
            `</div>`
          );
        },
      },
      grid: {
        top: 8,
        bottom: viewMode === "chart" ? 50 : 28,
        left: "4%",
        right: "40px",
        containLabel: true,
      },
      xAxis: {
        type: "time",
        axisLine: { show: true, lineStyle: { color: textColor } },
        axisTick: {
          show: true,
          lineStyle: { color: textColor },
          length: 5,
          inside: false,
        },
        axisLabel: {
          color: textColor,
          fontSize: 11,
          fontFamily: "Merriweather Sans",
        },
        splitLine: { show: false },
      },
      yAxis: {
        type: "value",
        axisLine: { show: true, lineStyle: { color: textColor } },
        axisTick: {
          show: true,
          lineStyle: { color: textColor },
          length: 5,
          inside: false,
        },
        axisLabel: {
          color: textColor,
          fontSize: 11,
          fontFamily: "Merriweather Sans",
        },
        splitLine: { show: false },
        min: yMin,
      },
      color: PALETTE,
      series: [...series, baselineSeries] as any[],
    };
  }, [displayFighters, textColor, isDark, lockedFighter, viewMode]);

  /* ---------------------------------------------------------------- */
  /*  DISTRIBUTION CALCULATIONS                                       */
  /* ---------------------------------------------------------------- */
  const allFightersForDist = useMemo(() => {
    if (!data) return [];
    const all = data.allFighters || [];
    const pool =
      selectedDivision === "overall"
        ? all
        : all.filter((f) => f.primaryWeightClass === selectedDivision);
    if (!activeOnly) return pool;
    const cutoff = Date.now() - TWO_YEARS_MS;
    return pool.filter(
      (f) => f.lastFightDate && new Date(f.lastFightDate).getTime() >= cutoff,
    );
  }, [data, selectedDivision, activeOnly]);

  const distStats = useMemo(() => {
    if (allFightersForDist.length === 0) return null;
    const elos = [...allFightersForDist.map((f) => f.currentElo)].sort(
      (a, b) => a - b,
    );
    const n = elos.length;
    const mean = elos.reduce((a, b) => a + b, 0) / n;
    const median =
      n % 2 === 0
        ? (elos[n / 2 - 1] + elos[n / 2]) / 2
        : elos[Math.floor(n / 2)];
    const variance = elos.reduce((sum, x) => sum + (x - mean) ** 2, 0) / n;
    const std = Math.sqrt(variance);
    const q1 = elos[Math.floor(n * 0.25)];
    const q3 = elos[Math.floor(n * 0.75)];
    return {
      n,
      mean,
      median,
      std,
      q1,
      q3,
      min: elos[0],
      max: elos[n - 1],
      elos,
    };
  }, [allFightersForDist]);

  const distChartOption = useMemo((): EChartsOption => {
    if (!distStats) return {};
    const BIN = DIST_BIN;
    const { elos, mean, median, std } = distStats;
    const start = Math.floor(elos[0] / BIN) * BIN;
    const end = Math.ceil((elos[elos.length - 1] + 1) / BIN) * BIN;
    const barXs: number[] = [];
    const barYs: number[] = [];
    for (let lo = start; lo < end; lo += BIN) {
      barXs.push(lo + BIN / 2);
      barYs.push(elos.filter((e) => e >= lo && e < lo + BIN).length);
    }
    const maxY = Math.max(...barYs);

    // KDE evaluated at fine grid
    const kdeFn = gaussianKDE(elos, Math.max(std * 0.35, 15));
    const kdePoints: number[][] = [];
    const step = (end - start) / 120;
    for (let x = start; x <= end; x += step) {
      kdePoints.push([x, kdeFn(x) * elos.length * BIN]);
    }

    return {
      backgroundColor: "transparent",
      animation: false,
      grid: { top: 52, bottom: 60, left: 58, right: 24 },
      tooltip: {
        trigger: "axis",
        axisPointer: { type: "shadow" },
        backgroundColor: isDark ? "#1f2937" : "#ffffff",
        borderColor: isDark ? "#374151" : "#e5e7eb",
        textStyle: {
          color: textColor,
          fontSize: 12,
          fontFamily: "Merriweather Sans",
        },
        formatter: (params: any) => {
          const bar = Array.isArray(params)
            ? params.find((p: any) => p.seriesName === "Count")
            : null;
          if (!bar) return "";
          const center = bar.data[0] as number;
          const lo = center - BIN / 2;
          const hi = center + BIN / 2;
          const count = bar.data[1] as number;
          let html = `<div style="font-family:'Merriweather Sans',sans-serif"><span style="font-weight:600">${lo}–${hi}</span><br/>Fighters: <b>${count}</b>`;
          if (count <= 30) {
            const inBin = allFightersForDist
              .filter((f) => f.currentElo >= lo && f.currentElo < hi)
              .sort((a, b) => b.currentElo - a.currentElo);
            if (inBin.length > 0) {
              html += `<div style="margin-top:5px;border-top:1px solid ${isDark ? "#4b5563" : "#d1d5db"};padding-top:4px">`;
              html += inBin
                .map(
                  (f) =>
                    `<div style="font-size:11px;font-family:'Merriweather Sans',sans-serif">${f.name} <span style="opacity:0.65">(${f.currentElo})</span></div>`,
                )
                .join("");
              html += `</div>`;
            }
          }
          html += `</div>`;
          return html;
        },
      },
      xAxis: {
        type: "value",
        name: "Elo Rating",
        nameLocation: "middle",
        nameGap: 32,
        nameTextStyle: {
          color: textColor,
          fontSize: 12,
          fontFamily: "Merriweather Sans",
        },
        axisLabel: {
          color: textColor,
          fontSize: 11,
          fontFamily: "Merriweather Sans",
        },
        axisLine: {
          show: true,
          lineStyle: { color: isDark ? "#ffffff" : "#000000" },
        },
        splitLine: { show: false },
        min: start,
        max: end,
      },
      yAxis: {
        type: "value",
        name: "Fighters",
        nameTextStyle: {
          color: textColor,
          fontSize: 12,
          fontFamily: "Merriweather Sans",
        },
        axisLabel: {
          color: textColor,
          fontSize: 11,
          fontFamily: "Merriweather Sans",
        },
        axisLine: {
          show: true,
          lineStyle: { color: isDark ? "#ffffff" : "#000000" },
        },
        splitLine: {
          lineStyle: { color: isDark ? "#374151" : "#f3f4f6", type: "dashed" },
        },
        max: Math.ceil(maxY * 1.15),
      },
      series: [
        {
          name: "Count",
          type: "bar" as const,
          data: barXs.map((x, i) => [x, barYs[i]]),
          barWidth: "90%",
          cursor: "pointer",
          z: 10,
          itemStyle: {
            color: isDark ? "rgba(99,136,235,0.75)" : "rgba(84,112,198,0.75)",
            borderRadius: [3, 3, 0, 0],
          },
          markLine: {
            silent: true,
            symbol: ["none", "none"],
            lineStyle: { type: "dashed" as const, width: 2 },
            label: {
              show: true,
              position: "insideStartTop" as const,
              fontSize: 11,
              fontWeight: "bold" as const,
            },
            data: [
              {
                xAxis: mean,
                lineStyle: { color: "#fac858" },
                label: { show: false },
              },
              {
                xAxis: median,
                lineStyle: { color: "#ee6666" },
                label: { show: false },
              },
            ],
          } as any,
        },
        {
          name: "KDE",
          type: "line" as const,
          data: kdePoints,
          smooth: true,
          symbol: "none",
          silent: true,
          lineStyle: { color: "#73c0de", width: 2.5 },
          itemStyle: { color: "#73c0de" },
          areaStyle: {
            color: isDark ? "rgba(115,192,222,0.08)" : "rgba(115,192,222,0.15)",
          },
          tooltip: { show: false } as any,
        },
      ],
    };
  }, [distStats, textColor, isDark, allFightersForDist]);

  // Click handler for distribution chart bins
  const distChartEvents = useMemo(
    () => ({
      click: (params: any) => {
        if (params.componentType !== "series" || params.seriesName !== "Count")
          return;
        const center = params.data?.[0] as number | undefined;
        if (center == null) return;
        const lo = center - DIST_BIN / 2;
        const hi = center + DIST_BIN / 2;
        const fighters = [...allFightersForDist]
          .filter((f) => f.currentElo >= lo && f.currentElo < hi)
          .sort((a, b) => b.currentElo - a.currentElo);
        setDistBinModal({ lo, hi, fighters });
      },
    }),
    [allFightersForDist],
  );

  const FighterNameTooltip = ({
    fighter,
    onSelect,
    showWeightClass = true,
  }: {
    fighter: any;
    onSelect: (name: string) => void;
    showWeightClass?: boolean;
  }) => {
    const isInactive =
      !!fighter.lastFightDate &&
      Date.now() - new Date(fighter.lastFightDate).getTime() >
        2 * 365.25 * 24 * 60 * 60 * 1000;
    const abbrev = WC_ABBREV[fighter.primaryWeightClass];

    return (
      <TooltipProvider delayDuration={0}>
        <Tooltip>
          <TooltipTrigger asChild>
            <button
              onClick={() => onSelect(fighter.name)}
              className="flex items-center gap-2 text-left hover:text-primary transition-colors w-full group"
            >
              <span className="font-medium group-hover:underline truncate max-w-[170px]">
                {fighter.name}
              </span>
              {showWeightClass && abbrev && (
                <span
                  className="text-[10px] font-semibold px-1.5 py-0.5 rounded border flex-shrink-0"
                  style={wcBadgeStyle(fighter.primaryWeightClass)}
                >
                  {abbrev}
                </span>
              )}
              {isInactive && (
                <span className="text-[10px] font-semibold px-1.5 py-0.5 rounded border border-gray-400/40 bg-gray-400/10 text-gray-500 dark:text-gray-400 flex-shrink-0">
                  Inactive
                </span>
              )}
            </button>
          </TooltipTrigger>
          <TooltipContent
            side="right"
            className="p-0 overflow-hidden bg-card text-card-foreground border border-border shadow-lg"
          >
            {(() => {
              const normKey = fighter.name
                .normalize("NFD")
                .replace(/[\u0300-\u036f]/g, "")
                .replace(/[.'\u2011]/g, "")
                .toLowerCase()
                .trim();
              const champEntry = championMap.get(normKey);
              const headshotUrl =
                champEntry?.headshotUrl ?? candidateHeadshots.get(normKey);
              const tott = candidateTott.get(normKey);
              return (
                <div className="flex gap-3 p-3 min-w-[200px] max-w-[240px]">
                  <div className="w-14 h-14 rounded-full overflow-hidden bg-muted flex items-center justify-center flex-shrink-0 border border-border">
                    {headshotUrl ? (
                      <img
                        src={headshotUrl}
                        alt={fighter.name}
                        className="w-full h-full object-cover object-top"
                      />
                    ) : (
                      <User className="w-6 h-6 text-muted-foreground" />
                    )}
                  </div>
                  <div className="flex flex-col gap-0.5 text-xs min-w-0 justify-center">
                    <div className="font-bold text-sm truncate leading-tight">
                      {fighter.name}
                    </div>
                    {WC_ABBREV[fighter.primaryWeightClass] && (
                      <span
                        className="text-[10px] font-semibold uppercase px-1.5 py-0.5 rounded border w-fit"
                        style={wcBadgeStyle(fighter.primaryWeightClass)}
                      >
                        {WC_ABBREV[fighter.primaryWeightClass]}
                      </span>
                    )}
                    {fighter.currentElo !== undefined && (
                      <div className="flex items-center gap-1 my-0.5">
                        <span className="text-[10px] text-foreground/70">
                          Elo:
                        </span>
                        <span className="inline-flex items-center px-1.5 py-0.5 rounded-full border border-yellow-400/60 bg-yellow-400/10 font-bold text-xs text-yellow-500 dark:text-yellow-400">
                          {fighter.currentElo}
                        </span>
                        {fighter.peakElo !== undefined && (
                          <span className="text-[10px] text-foreground/50">
                            / {fighter.peakElo}
                          </span>
                        )}
                      </div>
                    )}
                    {tott?.height && (
                      <div>
                        <span className="font-medium text-[10px]">
                          Height:
                        </span>{" "}
                        {tott.height}
                      </div>
                    )}
                    {tott?.reach && (
                      <div>
                        <span className="font-medium text-[10px]">Reach:</span>{" "}
                        {tott.reach}
                      </div>
                    )}
                    {tott?.stance && (
                      <div>
                        <span className="font-medium text-[10px]">Stance:</span>{" "}
                        {tott.stance}
                      </div>
                    )}
                    {tott?.dob && (
                      <div>
                        <span className="font-medium text-[10px]">Born:</span>{" "}
                        {formatDob(tott.dob)}
                      </div>
                    )}
                  </div>
                </div>
              );
            })()}
          </TooltipContent>
        </Tooltip>
      </TooltipProvider>
    );
  };

  /* ---------------------------------------------------------------- */
  /*  RENDER                                                          */
  /* ---------------------------------------------------------------- */
  if (loading) {
    return (
      <main className="flex flex-col items-center pt-8 px-4 md:px-12 pb-12 w-full max-w-screen-2xl mx-auto gap-6">
        {/* Chart card skeleton */}
        <div className="w-full rounded-xl border border-border bg-card p-3 md:p-5">
          <div className="flex justify-center gap-2 mb-3">
            <Skeleton className="h-5 w-36" />
            <Skeleton className="h-5 w-4" />
            <Skeleton className="h-5 w-44" />
          </div>
          <Skeleton className="w-full h-[480px] rounded-lg" />
          <div className="flex flex-wrap justify-center gap-x-4 gap-y-1.5 pt-2">
            {Array.from({ length: 10 }).map((_, i) => (
              <div key={i} className="flex items-center gap-1.5">
                <Skeleton className="w-2.5 h-2.5 rounded-full" />
                <Skeleton className="h-3 w-20" />
              </div>
            ))}
          </div>
        </div>
        {/* Table card skeleton */}
        <div className="w-full rounded-xl border border-border bg-card p-3">
          <div className="space-y-2">
            <Skeleton className="h-8 w-full" />
            {Array.from({ length: 10 }).map((_, i) => (
              <Skeleton key={i} className="h-9 w-full" />
            ))}
          </div>
        </div>
      </main>
    );
  }

  if (error) {
    return (
      <main className="flex flex-col items-center justify-center pt-24 px-8 md:px-32 min-h-[60vh]">
        <p className="text-destructive font-semibold">Error: {error}</p>
      </main>
    );
  }

  const divisionLabel =
    selectedDivision === "overall"
      ? "All Weight Classes (P4P)"
      : data?.wcDisplayNames[selectedDivision] || selectedDivision;

  const thSort = (col: SortCol, label: string) => (
    <th
      className={`py-3 px-3 font-semibold text-center cursor-pointer select-none hover:text-primary transition-colors${sortCol === col ? " text-primary" : ""}`}
      onClick={() => handleSort(col)}
    >
      {label} <SortArrow col={col} />
    </th>
  );

  return (
    <main
      className={`flex flex-col px-4 md:px-12 w-full max-w-screen-2xl mx-auto ${
        viewMode === "both" && !searchQuery
          ? "pt-3 pb-12 items-center"
          : "h-full overflow-hidden pt-3 pb-0"
      }`}
    >
      {/* Top header row: view buttons LEFT (‘Both’ first), division pills (table), search RIGHT */}
      <div className="w-full flex items-center gap-2 mb-2 flex-shrink-0">
        <TooltipProvider delayDuration={0}>
          <div className="flex gap-1 flex-shrink-0">
            <Tooltip>
              <TooltipTrigger asChild>
                <button
                  onClick={() => setViewMode("both")}
                  className={`p-1.5 rounded border transition-colors ${viewMode === "both" ? "bg-primary text-primary-foreground border-primary" : "bg-card text-muted-foreground border-border hover:text-foreground hover:bg-muted"}`}
                >
                  <LayoutGrid className="w-3.5 h-3.5" />
                </button>
              </TooltipTrigger>
              <TooltipContent side="bottom">Chart and Table</TooltipContent>
            </Tooltip>
            <Tooltip>
              <TooltipTrigger asChild>
                <button
                  onClick={() => setViewMode("chart")}
                  className={`p-1.5 rounded border transition-colors ${viewMode === "chart" ? "bg-primary text-primary-foreground border-primary" : "bg-card text-muted-foreground border-border hover:text-foreground hover:bg-muted"}`}
                >
                  <BarChart2 className="w-3.5 h-3.5" />
                </button>
              </TooltipTrigger>
              <TooltipContent side="bottom">Chart only</TooltipContent>
            </Tooltip>
            <Tooltip>
              <TooltipTrigger asChild>
                <button
                  onClick={() => setViewMode("table")}
                  className={`p-1.5 rounded border transition-colors ${viewMode === "table" ? "bg-primary text-primary-foreground border-primary" : "bg-card text-muted-foreground border-border hover:text-foreground hover:bg-muted"}`}
                >
                  <Table2 className="w-3.5 h-3.5" />
                </button>
              </TooltipTrigger>
              <TooltipContent side="bottom">Table only</TooltipContent>
            </Tooltip>
            <Tooltip>
              <TooltipTrigger asChild>
                <button
                  onClick={() => setViewMode("distribution")}
                  className={`p-1.5 rounded border transition-colors ${viewMode === "distribution" ? "bg-primary text-primary-foreground border-primary" : "bg-card text-muted-foreground border-border hover:text-foreground hover:bg-muted"}`}
                >
                  <Activity className="w-3.5 h-3.5" />
                </button>
              </TooltipTrigger>
              <TooltipContent side="bottom">Elo Distribution</TooltipContent>
            </Tooltip>
          </div>
        </TooltipProvider>

        {/* Division pills — all modes, hidden when searching */}
        {!searchQuery && (
          <div
            className="flex-1 flex gap-1 justify-center overflow-x-auto min-w-0"
            style={{ scrollbarWidth: "none" }}
          >
            <button
              onClick={() => setSelectedDivision("overall")}
              className={cn(
                "px-2 py-0.5 rounded-full text-xs font-semibold whitespace-nowrap border-2 transition-all flex-shrink-0",
                selectedDivision === "overall"
                  ? "bg-primary text-primary-foreground border-primary"
                  : "bg-card border-border text-muted-foreground hover:text-foreground hover:bg-muted",
              )}
            >
              P4P
            </button>
            {divisions.map((key) => (
              <button
                key={key}
                onClick={() => setSelectedDivision(key)}
                className={cn(
                  "px-2 py-0.5 rounded-full text-xs font-semibold whitespace-nowrap border-2 transition-all flex-shrink-0",
                  selectedDivision === key
                    ? "shadow-sm"
                    : "opacity-60 hover:opacity-90",
                )}
                style={wcBadgeStyle(key)}
              >
                {WC_ABBREV[key] || key}
              </button>
            ))}
          </div>
        )}
        {searchQuery && <div className="flex-1" />}

        {/* Active toggle + Search bar — RIGHT */}
        <div className="flex items-center gap-1 flex-shrink-0">
          <button
            onClick={() => setActiveOnly((prev) => !prev)}
            className={cn(
              "h-7 px-2 text-xs rounded-lg border transition-colors whitespace-nowrap",
              activeOnly
                ? "bg-primary text-primary-foreground border-primary"
                : "bg-card text-muted-foreground border-border hover:text-foreground hover:bg-muted",
            )}
          >
            Active only
          </button>
          {searchQuery && (
            <button
              onClick={() => {
                setSearchQuery("");
                setSearchResults([]);
              }}
              className="p-1 text-muted-foreground hover:text-foreground rounded transition-colors"
              title="Back"
            >
              <ArrowLeft size={14} />
            </button>
          )}
          <div className="relative">
            <Search
              size={12}
              className="absolute left-2 top-1/2 -translate-y-1/2 text-muted-foreground pointer-events-none"
            />
            <input
              type="text"
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              placeholder="Search fighter…"
              className="pl-6 pr-2 h-7 text-xs rounded-lg border border-border bg-card focus:outline-none focus:ring-1 focus:ring-primary w-40"
            />
          </div>
        </div>
      </div>

      {/* Chart + Table — hidden while searching */}
      {!searchQuery && (
        <>
          {/* Chart Card */}
          {(viewMode === "chart" || viewMode === "both") && (
            <Card
              className={`w-full ${viewMode === "both" ? "mb-6" : ""} ${viewMode === "chart" ? "flex-1 min-h-0 flex flex-col overflow-hidden" : ""}`}
            >
              <CardContent
                className={`p-1.5 md:p-1 ${viewMode === "chart" ? "flex-1 min-h-0 flex flex-col overflow-hidden" : ""}`}
              >
                {/* Title */}
                <div className="flex items-center justify-center mb-1">
                  <TooltipProvider delayDuration={0}>
                    <Tooltip>
                      <TooltipTrigger asChild>
                        <span className="font-bold text-base whitespace-nowrap cursor-help underline decoration-dotted decoration-muted-foreground/50">
                          UFC Elo Ratings
                        </span>
                      </TooltipTrigger>
                      <TooltipContent
                        className="max-w-xs text-xs p-2"
                        side="bottom"
                      >
                        <div className="space-y-1">
                          <div>
                            <span className="font-semibold">Base Elo:</span>{" "}
                            1500
                          </div>
                          <div>
                            <span className="font-semibold">K-factor:</span> 40
                            (debut) → 24 (30+ fights), linear taper
                          </div>
                          <div>
                            <span className="font-semibold">Finish bonus:</span>{" "}
                            +10 pts
                          </div>
                        </div>
                      </TooltipContent>
                    </Tooltip>
                  </TooltipProvider>
                </div>

                {/* Chart */}
                <div
                  ref={chartContainerRef}
                  className={`flex justify-center w-full${viewMode === "chart" ? " flex-1 min-h-0" : ""}`}
                  onMouseMove={handleChartMouseMove}
                  onMouseLeave={handleChartMouseLeave}
                >
                  <ReactECharts
                    ref={chartRef}
                    option={chartOption}
                    style={{
                      width: "100%",
                      height: viewMode === "chart" ? "100%" : "480px",
                    }}
                    notMerge
                    lazyUpdate
                    onEvents={chartEvents}
                    onChartReady={handleChartReady}
                  />
                </div>

                {/* HTML legend — hidden in chart-only mode */}
                {viewMode !== "chart" && (
                  <div className="flex flex-wrap justify-center gap-x-4 gap-y-1.5 pt-2 pb-1">
                    {displayFighters.map((f, i) => {
                      const isActive =
                        lockedFighter === null || lockedFighter === f.name;
                      const lgKey = f.name
                        .normalize("NFD")
                        .replace(/[\u0300-\u036f]/g, "")
                        .replace(/[.'\u2011]/g, "")
                        .toLowerCase()
                        .trim();
                      const lgChamp = championMap.get(lgKey);
                      const lgHeadshot =
                        lgChamp?.headshotUrl ?? candidateHeadshots.get(lgKey);
                      const lgTott = candidateTott.get(lgKey);
                      const lgWins = f.history.filter(
                        (h) => h.result === "W",
                      ).length;
                      const lgLosses = f.history.filter(
                        (h) => h.result === "L",
                      ).length;
                      const lgDraws = f.history.filter(
                        (h) => h.result === "D",
                      ).length;
                      return (
                        <TooltipProvider key={f.name} delayDuration={400}>
                          <Tooltip>
                            <TooltipTrigger asChild>
                              <div
                                onClick={() =>
                                  setLockedFighter((prev) =>
                                    prev === f.name ? null : f.name,
                                  )
                                }
                                className={`flex items-center gap-1.5 text-xs transition-opacity select-none cursor-pointer hover:opacity-100 ${
                                  isActive
                                    ? "opacity-100 text-foreground"
                                    : "opacity-30 text-muted-foreground"
                                }`}
                              >
                                <div
                                  className="w-2.5 h-2.5 rounded-full flex-shrink-0"
                                  style={{
                                    backgroundColor:
                                      PALETTE[i % PALETTE.length],
                                  }}
                                />
                                <span>{f.name}</span>
                              </div>
                            </TooltipTrigger>
                            <TooltipContent
                              side="top"
                              className="p-0 overflow-hidden bg-card text-card-foreground border border-border shadow-lg"
                            >
                              <div className="flex gap-3 p-3 min-w-[200px] max-w-[240px]">
                                <div className="w-12 h-12 rounded-full overflow-hidden bg-muted flex items-center justify-center flex-shrink-0 border border-border">
                                  {lgHeadshot ? (
                                    <img
                                      src={lgHeadshot}
                                      alt={f.name}
                                      className="w-full h-full object-cover object-top"
                                    />
                                  ) : (
                                    <User className="w-6 h-6 text-muted-foreground" />
                                  )}
                                </div>
                                <div className="flex flex-col gap-0.5 text-xs min-w-0 justify-center">
                                  <div className="font-bold text-sm truncate leading-tight">
                                    {f.name}
                                  </div>
                                  {WC_ABBREV[f.primaryWeightClass] && (
                                    <span
                                      className="text-[10px] font-semibold uppercase px-1.5 py-0.5 rounded border w-fit"
                                      style={wcBadgeStyle(f.primaryWeightClass)}
                                    >
                                      {WC_ABBREV[f.primaryWeightClass]}
                                    </span>
                                  )}
                                  {/* Always show Elo with golden number only */}
                                  <div className="flex items-center gap-1 my-0.5">
                                    <span className="text-[10px] text-foreground/70">
                                      Elo:
                                    </span>
                                    <span className="inline-flex items-center px-1.5 py-0.5 rounded-full border border-yellow-400/60 bg-yellow-400/10 font-bold text-xs text-yellow-500 dark:text-yellow-400">
                                      {f.currentElo}
                                    </span>
                                    <span className="text-[10px] text-foreground/50">
                                      / {f.peakElo}
                                    </span>
                                  </div>
                                  {lgTott ? (
                                    <>
                                      {lgTott.height && (
                                        <div>
                                          <span className="font-medium">
                                            Height:
                                          </span>{" "}
                                          {lgTott.height}
                                        </div>
                                      )}
                                      {lgTott.weight && (
                                        <div>
                                          <span className="font-medium">
                                            Weight:
                                          </span>{" "}
                                          {lgTott.weight}
                                        </div>
                                      )}
                                      {lgTott.reach && (
                                        <div>
                                          <span className="font-medium">
                                            Reach:
                                          </span>{" "}
                                          {lgTott.reach}
                                        </div>
                                      )}
                                      {lgTott.stance && (
                                        <div>
                                          <span className="font-medium">
                                            Stance:
                                          </span>{" "}
                                          {lgTott.stance}
                                        </div>
                                      )}
                                      {lgTott.dob && (
                                        <div>
                                          <span className="font-medium">
                                            Born:
                                          </span>{" "}
                                          {formatDob(lgTott.dob)}
                                        </div>
                                      )}
                                    </>
                                  ) : (
                                    <div>
                                      {lgWins}-{lgLosses}
                                      {lgDraws > 0 ? `-${lgDraws}` : ""}
                                    </div>
                                  )}
                                </div>
                              </div>
                            </TooltipContent>
                          </Tooltip>
                        </TooltipProvider>
                      );
                    })}
                  </div>
                )}
              </CardContent>
            </Card>
          )}

          {/* Rankings table */}
          {(viewMode === "table" || viewMode === "both") && (
            <Card
              className={`w-full ${viewMode === "table" ? "flex-1 min-h-0 flex flex-col overflow-hidden" : ""}`}
            >
              <CardContent
                className={`p-0 md:p-0 ${viewMode === "table" ? "flex-1 min-h-0 flex flex-col overflow-hidden" : ""}`}
              >
                <div
                  ref={tableWrapperRef}
                  className={`overflow-x-auto ${viewMode === "table" ? "flex-1 min-h-0 overflow-y-hidden" : ""}`}
                >
                  <table
                    className={`w-full text-sm${viewMode === "table" ? " h-full" : ""}`}
                  >
                    <thead className="sticky top-0 bg-card z-10">
                      <tr className="border-b border-border text-left">
                        <th className="py-3 px-3 font-semibold w-8 text-center text-foreground">
                          #
                        </th>
                        <th className="py-3 px-3 font-semibold text-left">
                          Fighter
                        </th>
                        {thSort("currentElo", "Current Elo")}
                        {thSort("peakElo", "Peak Elo")}
                        {thSort("peakDate", "Peak Date")}
                        {thSort("fights", "UFC Fights")}
                        {thSort("l5Wins", "L5")}
                        {thSort("streak", "Strk")}
                        {thSort("wl", "W/L")}
                      </tr>
                    </thead>
                    <tbody>
                      {displayFighters.map((f, i) => {
                        const wins = f.history.filter(
                          (h) => h.result === "W",
                        ).length;
                        const losses = f.history.filter(
                          (h) => h.result === "L",
                        ).length;
                        const draws = f.history.filter(
                          (h) => h.result === "D",
                        ).length;
                        const record =
                          draws > 0
                            ? `${wins}-${losses}-${draws}`
                            : `${wins}-${losses}`;
                        const wlRatioRaw =
                          losses === 0
                            ? wins > 0
                              ? Infinity
                              : 0
                            : wins / losses;
                        const wlRatio =
                          losses === 0
                            ? wins > 0
                              ? "∞"
                              : "—"
                            : Number.isInteger(wlRatioRaw)
                              ? String(wlRatioRaw)
                              : wlRatioRaw.toFixed(2);
                        const { l5Str, streak } = f as typeof f & {
                          l5Str: string;
                          streak: number;
                        };
                        const strkDisplay =
                          streak === 0
                            ? "—"
                            : streak > 0
                              ? `W${streak}`
                              : `L${Math.abs(streak)}`;
                        const strkColor =
                          streak > 0
                            ? "text-green-500"
                            : streak < 0
                              ? "text-red-500"
                              : "text-foreground";
                        return (
                          <tr
                            key={f.name}
                            className="border-b border-border/50 hover:bg-muted/30 last:border-b-0"
                            style={
                              viewMode === "table"
                                ? { height: tableRowHeight ?? undefined }
                                : undefined
                            }
                          >
                            <td className="py-2 px-3 text-center text-foreground">
                              {i + 1}
                            </td>
                            <td className="py-2 px-3 font-medium text-left">
                              <TooltipProvider delayDuration={400}>
                                <Tooltip>
                                  <TooltipTrigger asChild>
                                    <button
                                      onClick={() => {
                                        const key = f.name
                                          .normalize("NFD")
                                          .replace(/[\u0300-\u036f]/g, "")
                                          .replace(/[.'\u2011]/g, "")
                                          .toLowerCase()
                                          .trim();
                                        const champ = championMap.get(key);
                                        setSelectedChampion(
                                          champ
                                            ? {
                                                ...champ,
                                                eloHistory: f.history,
                                              }
                                            : {
                                                name: f.name,
                                                eloHistory: f.history,
                                              },
                                        );
                                      }}
                                      className="hover:underline hover:text-primary cursor-pointer text-left font-medium"
                                    >
                                      {f.name}
                                    </button>
                                  </TooltipTrigger>
                                  <TooltipContent
                                    side="right"
                                    className="p-0 overflow-hidden bg-card text-card-foreground border border-border shadow-lg"
                                  >
                                    {(() => {
                                      const normKey = f.name
                                        .normalize("NFD")
                                        .replace(/[\u0300-\u036f]/g, "")
                                        .replace(/[.'\u2011]/g, "")
                                        .toLowerCase()
                                        .trim();
                                      const champEntry =
                                        championMap.get(normKey);
                                      const headshotUrl =
                                        champEntry?.headshotUrl ??
                                        candidateHeadshots.get(normKey);
                                      const tott = candidateTott.get(normKey);
                                      return (
                                        <div className="flex gap-3 p-3 min-w-[200px] max-w-[240px]">
                                          <div className="w-14 h-14 rounded-full overflow-hidden bg-muted flex items-center justify-center flex-shrink-0 border border-border">
                                            {headshotUrl ? (
                                              <img
                                                src={headshotUrl}
                                                alt={f.name}
                                                className="w-full h-full object-cover object-top"
                                              />
                                            ) : (
                                              <User className="w-6 h-6 text-muted-foreground" />
                                            )}
                                          </div>
                                          <div className="flex flex-col gap-0.5 text-xs min-w-0 justify-center">
                                            <div className="font-bold text-sm truncate leading-tight">
                                              {f.name}
                                            </div>
                                            {WC_ABBREV[
                                              f.primaryWeightClass
                                            ] && (
                                              <span
                                                className="text-[10px] font-semibold uppercase px-1.5 py-0.5 rounded border w-fit"
                                                style={wcBadgeStyle(
                                                  f.primaryWeightClass,
                                                )}
                                              >
                                                {
                                                  WC_ABBREV[
                                                    f.primaryWeightClass
                                                  ]
                                                }
                                              </span>
                                            )}
                                            {f.currentElo !== undefined && (
                                              <div className="flex items-center gap-1 my-0.5">
                                                <span className="text-[10px] text-foreground/70">
                                                  Elo:
                                                </span>
                                                <span className="inline-flex items-center px-1.5 py-0.5 rounded-full border border-yellow-400/60 bg-yellow-400/10 font-bold text-xs text-yellow-500 dark:text-yellow-400">
                                                  {f.currentElo}
                                                </span>
                                                {f.peakElo !== undefined && (
                                                  <span className="text-[10px] text-foreground/50">
                                                    / {f.peakElo}
                                                  </span>
                                                )}
                                              </div>
                                            )}
                                            {tott?.height && (
                                              <div>
                                                <span className="font-medium">
                                                  Height:
                                                </span>{" "}
                                                {tott.height}
                                              </div>
                                            )}
                                            {tott?.weight && (
                                              <div>
                                                <span className="font-medium">
                                                  Weight:
                                                </span>{" "}
                                                {tott.weight}
                                              </div>
                                            )}
                                            {tott?.reach && (
                                              <div>
                                                <span className="font-medium">
                                                  Reach:
                                                </span>{" "}
                                                {tott.reach}
                                              </div>
                                            )}
                                            {tott?.stance && (
                                              <div>
                                                <span className="font-medium">
                                                  Stance:
                                                </span>{" "}
                                                {tott.stance}
                                              </div>
                                            )}
                                            {tott?.dob && (
                                              <div>
                                                <span className="font-medium">
                                                  Born:
                                                </span>{" "}
                                                {formatDob(tott.dob)}
                                              </div>
                                            )}
                                          </div>
                                        </div>
                                      );
                                    })()}
                                  </TooltipContent>
                                </Tooltip>
                              </TooltipProvider>
                              {selectedDivision === "overall" &&
                                WC_ABBREV[f.primaryWeightClass] && (
                                  <TooltipProvider delayDuration={0}>
                                    <Tooltip>
                                      <TooltipTrigger asChild>
                                        <span
                                          className="ml-1.5 text-xs font-semibold rounded px-1 py-0.5 cursor-help border"
                                          style={wcBadgeStyle(
                                            f.primaryWeightClass,
                                          )}
                                        >
                                          {WC_ABBREV[f.primaryWeightClass]}
                                        </span>
                                      </TooltipTrigger>
                                      <TooltipContent>
                                        {WC_FULL_NAME[f.primaryWeightClass] ||
                                          f.primaryWeightClass}
                                      </TooltipContent>
                                    </Tooltip>
                                  </TooltipProvider>
                                )}
                              <span className="ml-1.5 text-xs text-foreground font-normal">
                                ({record})
                              </span>
                            </td>
                            <td className="py-2 px-3 text-center font-bold">
                              {f.currentElo}
                            </td>
                            <td className="py-1 px-3 text-center">
                              {f.peakElo}
                            </td>
                            <td className="py-1 px-3 text-center">
                              {f.peakDate
                                ? (() => {
                                    const peakFight = f.history.find(
                                      (h) => h.date === f.peakDate,
                                    );
                                    const dateStr = new Date(
                                      f.peakDate,
                                    ).toLocaleDateString();
                                    const eventName = peakFight?.event;
                                    const inner = peakFight?.url ? (
                                      <a
                                        href={peakFight.url}
                                        target="_blank"
                                        rel="noopener noreferrer"
                                        className="hover:underline hover:text-primary"
                                      >
                                        {dateStr}
                                      </a>
                                    ) : (
                                      <span>{dateStr}</span>
                                    );
                                    return eventName ? (
                                      <TooltipProvider delayDuration={0}>
                                        <Tooltip>
                                          <TooltipTrigger asChild>
                                            <span>{inner}</span>
                                          </TooltipTrigger>
                                          <TooltipContent
                                            side="top"
                                            className="text-xs"
                                          >
                                            {eventName}
                                          </TooltipContent>
                                        </Tooltip>
                                      </TooltipProvider>
                                    ) : (
                                      inner
                                    );
                                  })()
                                : "—"}
                            </td>
                            <td className="py-1 px-3 text-center">
                              {f.fights}
                            </td>
                            <td className="py-1 px-3 text-center text-foreground">
                              {l5Str}
                            </td>
                            <td
                              className={`py-1 px-3 text-center font-semibold ${strkColor}`}
                            >
                              {strkDisplay}
                            </td>
                            <td className="py-1 px-3 text-center text-foreground">
                              {wlRatio}
                            </td>
                          </tr>
                        );
                      })}
                    </tbody>
                  </table>
                </div>
              </CardContent>
            </Card>
          )}

          {/* Elo Distribution view */}
          {viewMode === "distribution" && (
            <Card className="w-full flex-1 min-h-0 flex flex-col overflow-hidden">
              <CardContent className="p-3 md:p-4 flex-1 min-h-0 flex flex-col gap-2 overflow-y-auto">
                {/* Title */}
                <div className="flex items-center justify-center">
                  <span className="font-bold text-base">Elo Distribution</span>
                </div>

                {/* Histogram chart */}
                <div className="flex-1 min-h-0" style={{ minHeight: 240 }}>
                  <ReactECharts
                    option={distChartOption}
                    style={{ width: "100%", height: "100%" }}
                    notMerge
                    onEvents={distChartEvents}
                  />
                </div>

                {/* Legend */}
                <div className="flex justify-center gap-6 text-xs text-muted-foreground">
                  <span className="flex items-center gap-1.5">
                    <span className="inline-block w-5 h-0.5 bg-[#fac858]" />{" "}
                    Mean
                  </span>
                  <span className="flex items-center gap-1.5">
                    <span className="inline-block w-5 h-0.5 bg-[#ee6666]" />{" "}
                    Median
                  </span>
                  <span className="flex items-center gap-1.5">
                    <span
                      className="inline-block w-5 h-0.5 rounded bg-[#73c0de]"
                      style={{ height: 2 }}
                    />{" "}
                    KDE curve
                  </span>
                </div>

                {/* Stats cards */}
                {distStats && (
                  <div className="grid grid-cols-4 sm:grid-cols-7 gap-2 w-full">
                    {(
                      [
                        {
                          label: "Fighters",
                          value: distStats.n,
                          tooltip: "Fighters with ≥5 UFC fights",
                        },
                        { label: "Mean", value: Math.round(distStats.mean) },
                        {
                          label: "Median",
                          value: Math.round(distStats.median),
                        },
                        { label: "Std Dev", value: Math.round(distStats.std) },
                        {
                          label: "Q1 / Q3",
                          value: `${distStats.q1} / ${distStats.q3}`,
                        },
                        { label: "Min", value: distStats.min },
                        { label: "Max", value: distStats.max },
                      ] as {
                        label: string;
                        value: number | string;
                        tooltip?: string;
                      }[]
                    ).map(({ label, value, tooltip }) => {
                      const cell = (
                        <div
                          className={cn(
                            "rounded-lg border border-border bg-muted/30 px-3 py-2 text-center",
                            tooltip && "cursor-help",
                          )}
                        >
                          <div className="text-[10px] text-muted-foreground uppercase tracking-wide font-semibold">
                            {label}
                          </div>
                          <div className="text-sm font-bold mt-0.5">
                            {value}
                          </div>
                        </div>
                      );
                      return tooltip ? (
                        <TooltipProvider key={label} delayDuration={0}>
                          <Tooltip>
                            <TooltipTrigger asChild>{cell}</TooltipTrigger>
                            <TooltipContent className="text-xs">
                              {tooltip}
                            </TooltipContent>
                          </Tooltip>
                        </TooltipProvider>
                      ) : (
                        <React.Fragment key={label}>{cell}</React.Fragment>
                      );
                    })}
                  </div>
                )}
              </CardContent>
            </Card>
          )}
        </>
      )}

      {/* Fighter search results */}
      {searchQuery && (
        <div className="w-full flex-1 overflow-y-auto">
          {isSearching ? (
            <div className="flex justify-center py-12">
              <Loader2 className="animate-spin h-6 w-6 text-muted-foreground" />
            </div>
          ) : searchResults.length === 0 ? (
            <p className="text-center text-muted-foreground py-12 text-sm">
              No fighters found matching &ldquo;{searchQuery}&rdquo;
            </p>
          ) : (
            <Card className="w-full rounded-xl overflow-hidden">
              <CardContent className="p-0">
                <ul className="divide-y divide-border">
                  {searchResults.map((result) => (
                    <li key={result.name}>
                      <button
                        className="w-full flex items-center gap-3 px-4 py-3 hover:bg-muted/50 transition-colors text-left"
                        onClick={() => {
                          const key = result.name
                            .normalize("NFD")
                            .replace(/[\u0300-\u036f]/g, "")
                            .replace(/[.'\u2011]/g, "")
                            .toLowerCase()
                            .trim();
                          const champEntry = championMap.get(key);
                          setSelectedChampion(
                            champEntry
                              ? {
                                  ...champEntry,
                                  eloHistory: eloLookup(result.name),
                                }
                              : {
                                  name: result.name,
                                  eloHistory: eloLookup(result.name),
                                },
                          );
                        }}
                      >
                        <div className="w-12 h-12 rounded-full overflow-hidden bg-muted flex-shrink-0 flex items-center justify-center">
                          {result.headshot ? (
                            <img
                              src={result.headshot}
                              alt={result.name}
                              className="w-full h-full object-cover object-top"
                            />
                          ) : (
                            <User className="w-6 h-6 text-muted-foreground" />
                          )}
                        </div>
                        <div className="flex-1 min-w-0">
                          <div className="flex items-center gap-2">
                            <p className="font-semibold text-sm truncate">
                              {result.name}
                            </p>
                            {result.weightClass &&
                              WC_ABBREV[result.weightClass] && (
                                <span
                                  className="text-[10px] font-semibold px-1.5 py-0.5 rounded border flex-shrink-0"
                                  style={wcBadgeStyle(result.weightClass)}
                                >
                                  {WC_ABBREV[result.weightClass]}
                                </span>
                              )}
                            {result.eloStatus === "loaded" &&
                              result.lastFightDate &&
                              Date.now() -
                                new Date(result.lastFightDate).getTime() >
                                TWO_YEARS_MS && (
                                <span className="text-[10px] font-semibold px-1.5 py-0.5 rounded border border-gray-400/40 bg-gray-400/10 text-gray-500 dark:text-gray-400 flex-shrink-0">
                                  Inactive
                                </span>
                              )}
                          </div>
                          <div className="flex items-center gap-2 mt-0.5">
                            {result.record && (
                              <span className="text-xs text-muted-foreground">
                                {result.record}
                              </span>
                            )}
                            {result.fights !== undefined && (
                              <span className="text-xs text-muted-foreground">
                                · {result.fights} UFC fights
                              </span>
                            )}
                          </div>
                        </div>
                        <div className="text-right flex-shrink-0 ml-2">
                          {result.eloStatus === "loading" ? (
                            <div className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full border border-border bg-muted/30">
                              <span className="text-xs text-muted-foreground">
                                Elo:
                              </span>
                              <Loader2 className="w-3 h-3 animate-spin text-muted-foreground" />
                            </div>
                          ) : result.eloStatus === "loaded" ? (
                            <div className="flex items-center gap-3">
                              <div>
                                <p className="text-[10px] text-muted-foreground">
                                  Current
                                </p>
                                <p className="font-bold text-sm">
                                  {result.currentElo}
                                </p>
                              </div>
                              {result.peakElo !== undefined && (
                                <div>
                                  <p className="text-[10px] text-muted-foreground">
                                    Peak
                                  </p>
                                  <p className="font-semibold text-sm">
                                    {result.peakElo}
                                  </p>
                                </div>
                              )}
                            </div>
                          ) : null}
                        </div>
                      </button>
                    </li>
                  ))}
                </ul>
              </CardContent>
            </Card>
          )}
        </div>
      )}

      {/* Fighter details popup — only for real champions */}
      <AnimatePresence>
        {selectedChampion && (
          <FighterDetailsModal
            champion={selectedChampion}
            onClose={() => setSelectedChampion(null)}
            eloLookup={eloLookup}
          />
        )}
      </AnimatePresence>

      {/* Elo distribution bin popup — shows all fighters in the clicked bin */}
      <AnimatePresence>
        {distBinModal && (
          <BinFightersModal
            lo={distBinModal.lo}
            hi={distBinModal.hi}
            fighters={distBinModal.fighters}
            onClose={() => setDistBinModal(null)}
            renderFighterName={(fighter) => (
              <FighterNameTooltip
                fighter={fighter}
                onSelect={(name) => {
                  setDistBinModal(null);
                  const key = name
                    .normalize("NFD")
                    .replace(/[\u0300-\u036f]/g, "")
                    .replace(/[.'\u2011]/g, "")
                    .toLowerCase()
                    .trim();
                  const champEntry = championMap.get(key);
                  setSelectedChampion(
                    champEntry
                      ? { ...champEntry, eloHistory: eloLookup(name) }
                      : { name, eloHistory: eloLookup(name) },
                  );
                }}
              />
            )}
            onSelectFighter={(name) => {
              setDistBinModal(null);
              const key = name
                .normalize("NFD")
                .replace(/[\u0300-\u036f]/g, "")
                .replace(/[.'\u2011]/g, "")
                .toLowerCase()
                .trim();
              const champEntry = championMap.get(key);
              setSelectedChampion(
                champEntry
                  ? { ...champEntry, eloHistory: eloLookup(name) }
                  : { name, eloHistory: eloLookup(name) },
              );
            }}
          />
        )}
      </AnimatePresence>
    </main>
  );
}
