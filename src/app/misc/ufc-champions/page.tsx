"use client";

import React, { useState, useEffect, useRef, useMemo, useCallback } from "react";
import { Card, CardContent } from "@/components/ui/card";
import { Checkbox } from "@/components/ui/checkbox";
import { Label } from "@/components/ui/label";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import {
  supabase,
  FighterDetails,
  FightResult,
  UfcChampionImages,
} from "@/lib/supabaseClient";
import {
  X,
  Loader2,
  Crown,
  User,
  ShieldCheck, // ✅ successful defence
  ShieldOff,  // ❌ broken defence
  Info,
} from "lucide-react"; // Updated imports
import { cn } from "@/lib/utils";
import { motion, AnimatePresence } from "framer-motion";
import { Skeleton } from "@/components/ui/skeleton"; // Import Skeleton
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from "@/components/ui/tooltip";
import { type Champion, type WeightClassData} from "@/data/ufcChampionsData";
import { FighterDetailsModal, type ChampionEntry } from "@/components/FighterDetailsModal";
import { wcBadgeStyle } from "@/lib/wc-colors";

/* ------------------------------------------------------------------ */
/*  TYPES                                                             */
/* ------------------------------------------------------------------ */
interface ProcessedFightResult extends FightResult {
  resultForFighter: "win" | "loss" | "draw" | "nc" | "unknown";
  decisionDetails?: string;
  finishDetails?: string;
  defenseOutcome?: "success" | "fail" | null; // Updated type
}

/* ------------------------------------------------------------------ */
/*  HELPERS                                                           */
/* ------------------------------------------------------------------ */

const parseDate = (iso: string | null): Date | null => {
  if (!iso) return null;
  if (iso === "Present") return new Date();
  const parts = iso.split("-");
  if (parts.length === 1) return new Date(`${iso}-01-01`);
  if (parts.length === 2) return new Date(`${iso}-01`);
  return new Date(iso);
};

const diffInDays = (a: Date, b: Date) =>
  Math.round((b.getTime() - a.getTime()) / 86_400_000);

const humanDuration = (days: number) => {
  const y = Math.floor(days / 365);
  const m = Math.floor((days % 365) / 30);
  const d = days % 30;
  return [
    y ? `${y} yr${y > 1 ? "s" : ""}` : "",
    m ? `${m} mo${m > 1 ? "s" : ""}` : "",
    d ? `${d} day${d > 1 ? "s" : ""}` : "",
  ]
    .filter(Boolean)
    .join(" ");
};

const normalizeName = (name: string) =>
  name
    .normalize("NFD")
    .replace(/[\u0300-\u036f]/g, "")
    .replace(/ł/g, "l")
    .replace(/Ł/g, "L")
    .replace(/[.'‑]/g, "")
    .trim();

/* NEW → collapse every weight‑class string to its plain division name */
const normalizeWeightClass = (raw?: string | null): string => {
  if (!raw) return "";
  return normalizeName(raw) // Apply normalizeName first (handles accents etc.)
    .toLowerCase()
    .replace(
      /(title|championship|bout|world|undisputed|interim)/g, // Remove keywords
      ""
    )
    .replace(/^ufc\s*/, "") // <<<--- ADD THIS LINE: Remove leading "ufc " prefix
    .replace(/\s+/g, " ") // Clean up multiple spaces
    .trim(); // Trim final whitespace
    // e.g., "UFC Welterweight title bout" → "ufc welterweight title bout" → "ufc welterweight " → "welterweight " → "welterweight " → "welterweight"
};

const parseFighterName = (full: string) => {
  const parts = full.trim().split(/\s+/);
  return { firstName: parts.slice(0, -1).join(" "), lastName: parts.at(-1) || "" };
};

const extractOpponent = (boutString: string, normalizedFighterName: string): string => {
  if (!boutString || !normalizedFighterName) return 'Opponent Unknown';
  const parts = boutString.split(/\s+vs\.?\s+|\s+vs\s+/i);
  if (parts.length === 2) {
    const normalized1 = normalizeName(parts[0]).toLowerCase();
    const normalized2 = normalizeName(parts[1]).toLowerCase();
    const normalized = normalizedFighterName.toLowerCase();
    if (normalized1.includes(normalized) || normalized.includes(normalized1)) {
      return parts[1].trim();
    } else if (normalized2.includes(normalized) || normalized.includes(normalized2)) {
      return parts[0].trim();
    }
  }
  return 'Opponent Unknown';
};

const extractDecisionDetails = (details: string | null | undefined): string | undefined => {
  if (!details) return undefined;
  const scorePattern = /\d{1,2}\s*-\s*\d{1,2}/g;
  const scores = details.match(scorePattern);
  if (scores && scores.length >= 3) {
    return scores.join(', ');
  }
  return undefined;
}

const createAthleteUrl = (name: string): string => {
  if (!name) return '#';
  let slug = name
    .normalize("NFD")
    .replace(/[\u0300-\u036f]/g, "")
    .replace(/\./g, "")
    .trim()
    .toLowerCase()
    .replace(/\s+/g, '-');
  slug = slug.replace(/^b-?j-/, "bj-");
  return `https://www.ufc.com/athlete/${slug}`;
};

const getResultForFighter = (fight: FightResult, normalizedFighterName: string): ProcessedFightResult['resultForFighter'] => {
    if (!fight.BOUT || !fight.OUTCOME) return 'unknown';
    if (fight.DETAILS?.toLowerCase().includes('overturned')) return 'nc';
    if (fight.METHOD?.toLowerCase().includes('no contest') || fight.METHOD?.toLowerCase() === 'nc') return 'nc';
    if (fight.OUTCOME.toLowerCase() === 'nc' || fight.OUTCOME.toLowerCase() === 'nc/nc') return 'nc';
    if (fight.OUTCOME.toLowerCase() === 'd' || fight.OUTCOME.toLowerCase() === 'd/d' || fight.OUTCOME.toLowerCase() === 'draw') return 'draw';
    const outcomeParts = fight.OUTCOME.split('/');
    const fighters = fight.BOUT.split(/ vs\.? | vs /i);
    if (fighters.length !== 2) return 'unknown';
    const fighter1 = normalizeName(fighters[0].trim());
    const fighter2 = normalizeName(fighters[1].trim());
    const fighterIsFighter1 = fighter1.toLowerCase() === normalizedFighterName.toLowerCase();
    const fighterIsFighter2 = fighter2.toLowerCase() === normalizedFighterName.toLowerCase();
    if (fighterIsFighter1) {
        if (outcomeParts[0]?.toLowerCase() === 'w') return 'win';
        if (outcomeParts[0]?.toLowerCase() === 'l') return 'loss';
    } else if (fighterIsFighter2) {
        if (outcomeParts[1]?.toLowerCase() === 'w') return 'win';
        if (outcomeParts[1]?.toLowerCase() === 'l') return 'loss';
        if (outcomeParts.length === 1 && outcomeParts[0]?.toLowerCase() === 'w') return 'loss';
        if (outcomeParts.length === 1 && outcomeParts[0]?.toLowerCase() === 'l') return 'win';
    }
    return 'unknown';
};

const calculateAge = (dobString: string | null | undefined): number | null => {
  if (!dobString) return null;
  try {
    const birthDate = new Date(dobString);
    const today = new Date();
    let age = today.getFullYear() - birthDate.getFullYear();
    const m = today.getMonth() - birthDate.getMonth();
    if (m < 0 || (m === 0 && today.getDate() < birthDate.getDate())) {
      age--;
    }
    return age;
  } catch (e) {
    console.error("Error parsing DOB:", e);
    return null;
  }
};

const getInitials = (name: string): string => {
  if (!name) return "?";
  const normalizedName = name.toLowerCase().replace(/[‑\s]+/g, ' ').trim();
  if (normalizedName === "b j penn") return "BJP";
  if (normalizedName === "georges st pierre") return "GSP"; // <<<--- ADDED THIS LINE
  if (normalizedName === "rafael dos anjos") return "RDA";
  if (normalizedName === "dricus du plessis") return "DDP";
  const parts = name
    .replace(/‑/g, "-")
    .split(/[\s-]+/)
    .filter(Boolean);
  if (parts.length === 0) return "?";
  if (parts.length === 1) return parts[0].substring(0, 2).toUpperCase();
  if (parts.length === 3) {
    const middle = parts[1].toLowerCase();
    if (["dos", "de", "du", "da", "von", "van", "del"].includes(middle)) {
      return (parts[0][0] + parts[2][0]).toUpperCase();
    }
  }
  return (parts[0][0] + (parts.length > 1 ? parts[parts.length - 1][0] : "")).toUpperCase();
};

const abbreviateWeightClass = (key: string): string => {
    const abbreviations: Record<string, string> = {
        heavyweight: "HW", lightheavyweight: "LHW", middleweight: "MW",
        welterweight: "WW", lightweight: "LW", featherweight: "FeW",
        bantamweight: "BW", flyweight: "FlW", women_strawweight: "W-SW",
        women_flyweight: "W-FLW", women_bantamweight: "W-BW",
        catchweight: "CW",
    };
    // Remove all spaces for matching (e.g., "light heavyweight" -> "lightheavyweight")
    const cleanKey = key.toLowerCase().replace(/\s+/g, "");
    return abbreviations[cleanKey] || key.toUpperCase();
};

const getWeightClassName = (key: string): string => {
    const names: Record<string, string> = {
        heavyweight: "Heavyweight",
        lightheavyweight: "Light Heavyweight",
        middleweight: "Middleweight",
        welterweight: "Welterweight",
        lightweight: "Lightweight",
        featherweight: "Featherweight",
        bantamweight: "Bantamweight",
        flyweight: "Flyweight",
        women_strawweight: "Women's Strawweight",
        women_flyweight: "Women's Flyweight",
        women_bantamweight: "Women's Bantamweight",
        catchweight: "Catchweight",
    };
    const cleanKey = key.toLowerCase().replace(/\s+/g, "");
    return names[cleanKey] || key;
};

 /**
 * Calculates background and text colors for a champion reign bar
 * based on the duration in days, using an HSL blue gradient.
 * @param days - Duration of the reign in days.
 * @returns Object with { bg, text, hoverBg } HSL color strings.
 */
const getReignColor = (days: number): { bg: string; text: string; hoverBg: string } => {
  const minDays = 0;
  const maxDays = 2550; // Approx 7 years, adjust upper limit as needed

  // Define HSL values for start (dark blue) and end (light blue)
  const startHue = 220; // Darker blue hue
  const startSaturation = 60;
  const startLightness = 35; // Darker

  const endHue = 200; // Lighter blue hue
  const endSaturation = 85;
  const endLightness = 88; // Lighter

  // Clamp days to prevent going outside the 0-1 range
  const clampedDays = Math.max(minDays, Math.min(maxDays, days));

  // Calculate interpolation factor (0 for minDays, 1 for maxDays)
  const t = (clampedDays - minDays) / (maxDays - minDays);

  // Interpolate HSL values
  const hue = Math.round(startHue + t * (endHue - startHue));
  const saturation = Math.round(startSaturation + t * (endSaturation - startSaturation));
  const lightness = Math.round(startLightness + t * (endLightness - startLightness));

  const bgColor = `hsl(${hue}, ${saturation}%, ${lightness}%)`;

  // Determine text color based on background lightness for contrast
  // Use a threshold slightly below the middle (50) as light text on mid-tones is often easier to read
  const textColor = lightness > 55 ? `hsl(${hue}, 70%, 20%)` : `hsl(${hue}, 30%, 95%)`; // Dark text on light bg, light text on dark bg

  // Calculate a slightly darker hover background color
  const hoverLightness = Math.max(15, lightness - 8); // Darken by 8%, ensure min 15% lightness
   const hoverBgColor = `hsl(${hue}, ${saturation}%, ${hoverLightness}%)`;


  return { bg: bgColor, text: textColor, hoverBg: hoverBgColor };
};


/* ------------------------------------------------------------------ */
/*  COMPONENT                                                         */
/* ------------------------------------------------------------------ */
// Canonical dropdown order: men heaviest→lightest, then women heaviest→lightest
const DIVISION_ORDER = [
  'heavyweight', 'lightheavyweight', 'middleweight', 'welterweight',
  'lightweight', 'featherweight', 'bantamweight', 'flyweight',
  'women_bantamweight', 'women_flyweight', 'women_strawweight',
];

const UfcChampionsDisplay: React.FC = () => {
  // --- STATE DECLARATIONS ---
  const [selected, setSelected] = useState("lightweight");
  const [selectedChampionForModal, setSelectedChampionForModal] = useState<ChampionEntry | null>(null);
  const [isTimelineView, setIsTimelineView] = useState(false);
  const [lineageData, setLineageData] = useState<Record<string, WeightClassData>>({});
  const [isLoadingLineage, setIsLoadingLineage] = useState(true);
  const [isTimelineCalculating, setIsTimelineCalculating] = useState(false);
  const [genderFilter, setGenderFilter] = useState<"all" | "men" | "women">("all");

  // Elo lookup — fetched in background so modals can show Elo ratings
  const eloMapRef = useRef<Map<string, Array<{ date: string; elo: number; event: string }>>>(new Map());
  useEffect(() => {
    fetch('/api/ufc-elo')
      .then(r => r.ok ? r.json() : null)
      .then((data: { candidatePool?: any[]; byWeightClass?: Record<string, any[]> } | null) => {
        if (!data) return;
        const map = eloMapRef.current;
        const norm = (s: string) => s.normalize('NFD').replace(/[\u0300-\u036f]/g, '').replace(/[.'\u2011]/g, '').toLowerCase().trim();
        const add = (f: { name: string; history: Array<{ date: string; elo: number; event: string }> }) => {
          const key = norm(f.name);
          if (!map.has(key)) map.set(key, f.history);
        };
        (data.candidatePool || []).forEach(add);
        Object.values(data.byWeightClass || {}).flat().forEach((f: any) => add(f));
      })
      .catch(() => {}); // silently fail — Elo is optional
  }, []);
  const eloLookup = useCallback((name: string) => {
    const key = name.normalize('NFD').replace(/[\u0300-\u036f]/g, '').replace(/[.'\u2011]/g, '').toLowerCase().trim();
    return eloMapRef.current.get(key);
  }, []);

  const scrollContainerRef = useRef<HTMLDivElement>(null);
  const timelineScrollContainerRef = useRef<HTMLDivElement>(null);

  // --- DRAG SCROLL REFS (Card View) ---
  const isDraggingRef = useRef(false);
  const startXRef = useRef(0);
  const scrollLeftStartRef = useRef(0);
  const didDragRef = useRef(false);

  // --- DRAG SCROLL REFS (Timeline View) --- <<<--- NEW REFS
  const isTimelineDraggingRef = useRef(false);
  const timelineStartXRef = useRef(0);
  const timelineScrollLeftStartRef = useRef(0);
  const timelineDidDragRef = useRef(false);
  // --- END STATE/REFS ---

  // Fetch full historical lineage from the API on mount
  useEffect(() => {
    fetch('/api/ufc-champions/lineage')
      .then(r => r.json())
      .then(json => {
        if (json.lineage) {
          // Sort by canonical division order before storing
          const raw = json.lineage as Record<string, WeightClassData>;
          const sorted: Record<string, WeightClassData> = {};
          for (const key of DIVISION_ORDER) {
            if (raw[key]) sorted[key] = raw[key];
          }
          // Append any divisions not in DIVISION_ORDER
          for (const key of Object.keys(raw)) {
            if (!DIVISION_ORDER.includes(key)) sorted[key] = raw[key];
          }
          setLineageData(sorted);
        }
      })
      .catch(e => console.warn('[ufc-champions] lineage fetch failed:', e))
      .finally(() => setIsLoadingLineage(false));
  }, []);

  const classes = lineageData;
  const data = classes[selected];

  // Lineage API is the single source of truth – use champions directly.
  const cards = useMemo(() => data?.champions ?? [], [data?.champions]);

  const filteredClasses = useMemo(() => {
    if (genderFilter === 'all') return classes;
    return Object.entries(classes)
      .filter(([key]) => {
        const isWomen = key.startsWith('women_');
        return genderFilter === 'women' ? isWomen : !isWomen;
      })
      .reduce((acc, [key, value]) => {
        acc[key] = value;
        return acc;
      }, {} as Record<string, WeightClassData>);
  }, [classes, genderFilter]);

  // Effect for card view scroll to end - runs when selection changes or switching from timeline view
  useEffect(() => {
    if (scrollContainerRef.current && data && !isTimelineView && !isDraggingRef.current) {
      const timer = setTimeout(() => {
        if (scrollContainerRef.current) {
          scrollContainerRef.current.scrollLeft = scrollContainerRef.current.scrollWidth;
        }
      }, 50);
      return () => clearTimeout(timer);
    }
  }, [selected, data, isTimelineView]); 

  // Effect for timeline view scroll to end
  useEffect(() => {
    if (isTimelineView && timelineScrollContainerRef.current) {
      timelineScrollContainerRef.current.scrollLeft = timelineScrollContainerRef.current.scrollWidth;
    }
  }, [isTimelineView]);

  /* ================================================================== */
  /*  TIMELINE CALCULATIONS                                             */
  /* ================================================================== */
  const timelineData = useMemo(() => {
    if (!isTimelineView || Object.keys(filteredClasses).length === 0) return null;
    setIsTimelineCalculating(true);
    console.log("Calculating timeline data...");
    let minDate = new Date();
    let maxDate = new Date(1990, 0, 1);
    const today = new Date();

    Object.values(filteredClasses).forEach((wc: WeightClassData) => {
      wc.champions.forEach((c: Champion) => {
        const start = parseDate(c.reignStart);
        const end = parseDate(c.reignEnd);
        if (start && start < minDate) minDate = start;
        if (end) {
          const cappedEnd = end > today ? today : end;
          if (cappedEnd > maxDate) maxDate = cappedEnd;
        } else {
            if (start && start > maxDate) maxDate = today;
            else if (today > maxDate) maxDate = today;
        }
      });
    });

    if (maxDate > today) maxDate = today;
    maxDate.setDate(maxDate.getDate() + 10);
    minDate.setDate(minDate.getDate() - 10);

    const totalDays = diffInDays(minDate, maxDate);
    const pixelsPerDay = 0.75;
    const todayX = diffInDays(minDate, today) * pixelsPerDay;

    const monthMarkers: { date: Date; x: number }[] = [];
    let currentMonth = new Date(minDate.getFullYear(), minDate.getMonth(), 1);
    while (currentMonth <= maxDate) {
      const x = diffInDays(minDate, currentMonth) * pixelsPerDay;
      monthMarkers.push({ date: new Date(currentMonth), x });
      currentMonth.setMonth(currentMonth.getMonth() + 1);
    }

    console.log("Timeline data calculated.");
    setIsTimelineCalculating(false);
    return { minDate, maxDate, totalDays, pixelsPerDay, totalWidth: totalDays * pixelsPerDay, monthMarkers, todayX };
  }, [isTimelineView, filteredClasses]);

  /* ================================================================== */
  /*  DRAG SCROLL HANDLERS (Card View - Updated)                        */
  /* ================================================================== */
  const handlePointerDown = (e: React.PointerEvent<HTMLDivElement>) => {
    // Only act if card view is active and ref exists
    if (isTimelineView || selectedChampionForModal || !scrollContainerRef.current) return;
    if (e.button !== 0) return; // Only main button
    e.preventDefault(); // Prevent text selection etc.
    (e.target as Element).setPointerCapture(e.pointerId); // Capture pointer

    isDraggingRef.current = true;
    didDragRef.current = false; // <<<--- Reset didDrag flag here
    startXRef.current = e.pageX - scrollContainerRef.current.offsetLeft;
    scrollLeftStartRef.current = scrollContainerRef.current.scrollLeft;
    scrollContainerRef.current.style.cursor = 'grabbing';
    scrollContainerRef.current.style.userSelect = 'none'; // Prevent text selection during drag
  };

  const handlePointerMove = (e: React.PointerEvent<HTMLDivElement>) => {
    if (!isDraggingRef.current || isTimelineView || selectedChampionForModal || !scrollContainerRef.current) return;
    e.preventDefault();

    const x = e.pageX - scrollContainerRef.current.offsetLeft;
    const walk = x - startXRef.current;

    // Set didDrag flag if movement exceeds a small threshold (e.g., 5px)
    if (Math.abs(walk) > 5) { // <<<--- Set didDrag flag if moved enough
        didDragRef.current = true;
    }

    scrollContainerRef.current.scrollLeft = scrollLeftStartRef.current - (walk * 1.5); // Adjust multiplier for scroll speed
  };

  const handlePointerUpOrLeave = (e: React.PointerEvent<HTMLDivElement>) => {
    if (!isDraggingRef.current || !scrollContainerRef.current) return;
    try { // Add try-catch as releasePointerCapture can sometimes throw if element is removed
        (e.target as Element).releasePointerCapture(e.pointerId);
    } catch (error) {
        console.warn("Could not release pointer capture:", error);
    }

    isDraggingRef.current = false;
    // No need to reset didDragRef here, onClick uses the value set during move/down
    scrollContainerRef.current.style.cursor = 'grab';
    scrollContainerRef.current.style.userSelect = '';
  };

  /* ================================================================== */
  /*  TIMELINE DRAG SCROLL HANDLERS                                     */
  /* ================================================================== */
  const handleTimelinePointerDown = (e: React.PointerEvent<HTMLDivElement>) => {
    // Only act if timeline view is active and ref exists
    if (!isTimelineView || selectedChampionForModal || !timelineScrollContainerRef.current) return;
    if (e.button !== 0) return; // Only main button
    e.preventDefault();
    (e.target as Element).setPointerCapture(e.pointerId); // Capture pointer

    isTimelineDraggingRef.current = true;
    timelineDidDragRef.current = false; // Reset didDrag flag for timeline
    timelineStartXRef.current = e.pageX - timelineScrollContainerRef.current.offsetLeft;
    timelineScrollLeftStartRef.current = timelineScrollContainerRef.current.scrollLeft;
    timelineScrollContainerRef.current.style.cursor = 'grabbing';
    timelineScrollContainerRef.current.style.userSelect = 'none'; // Prevent text selection
  };

  const handleTimelinePointerMove = (e: React.PointerEvent<HTMLDivElement>) => {
    if (!isTimelineDraggingRef.current || !isTimelineView || selectedChampionForModal || !timelineScrollContainerRef.current) return;
    e.preventDefault();

    const x = e.pageX - timelineScrollContainerRef.current.offsetLeft;
    const walk = x - timelineStartXRef.current;

    // Set timelineDidDrag flag if movement exceeds threshold
    if (Math.abs(walk) > 5) {
        timelineDidDragRef.current = true;
    }

    timelineScrollContainerRef.current.scrollLeft = timelineScrollLeftStartRef.current - (walk * 1.5); // Adjust multiplier as needed
  };

  const handleTimelinePointerUpOrLeave = (e: React.PointerEvent<HTMLDivElement>) => {
    if (!isTimelineDraggingRef.current || !timelineScrollContainerRef.current) return;
    try { (e.target as Element).releasePointerCapture(e.pointerId); } catch (error) { console.warn("Could not release pointer capture:", error); }

    isTimelineDraggingRef.current = false;
    // Don't reset timelineDidDragRef here
    timelineScrollContainerRef.current.style.cursor = 'grab';
    timelineScrollContainerRef.current.style.userSelect = ''; // Re-allow selection
  };

  /* ================================================================== */
  /*  WHEEL → HORIZONTAL SCROLL HANDLERS                               */
  /* ================================================================== */
  const handleCardWheel = (e: React.WheelEvent<HTMLDivElement>) => {
    if (!scrollContainerRef.current) return;
    e.preventDefault();
    scrollContainerRef.current.scrollLeft += e.deltaY + e.deltaX;
  };

  const handleTimelineWheel = (e: React.WheelEvent<HTMLDivElement>) => {
    if (!timelineScrollContainerRef.current) return;
    e.preventDefault();
    timelineScrollContainerRef.current.scrollLeft += e.deltaY + e.deltaX;
  };


  /* ------------------------------------------------------------------ */
  /*  UI                                                                */
  /* ------------------------------------------------------------------ */
  return (
    <main className="relative pt-8 md:pt-6 px-4 md:px-6 w-full flex flex-col h-full overflow-hidden">
      {/* --- View Toggle and Gender Filter --- */}
      <div className="w-full max-w-screen-2xl mx-auto flex justify-between items-start mb-4 px-0 md:px-0 min-h-[40px]">
        <div className="flex-1">
          {isTimelineView && (
            <div className="flex items-center space-x-2">
              <Select value={genderFilter} onValueChange={(value: "all" | "men" | "women") => setGenderFilter(value)} disabled={!!selectedChampionForModal}>
                <SelectTrigger className="w-180px]">
                  <SelectValue placeholder="Filter Divisions" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="men">Men's Divisions</SelectItem>
                  <SelectItem value="women">Women's Divisions</SelectItem>
                  <SelectItem value="all">All Divisions</SelectItem>
                </SelectContent>
              </Select>
            </div>
          )}
        </div>
        <div className="flex items-center space-x-2 pt-2">
          <Checkbox id="timeline-view-toggle" checked={isTimelineView} onCheckedChange={(checked) => { const isChecked = Boolean(checked); setIsTimelineView(isChecked); setGenderFilter(isChecked ? "men" : "all"); }} disabled={!!selectedChampionForModal} />
          <Label htmlFor="timeline-view-toggle" className="text-sm font-medium">Show Full Timeline View</Label>
        </div>
      </div>
      {/* --- Conditional Rendering: Card View or Timeline View --- */}
      {!isTimelineView ? (
        /* --- Main Content Card (Background) --- */
        (<Card className={cn( "w-full max-w-screen-2xl mx-auto transition-opacity duration-300 flex flex-col mb-auto", selectedChampionForModal ? "opacity-20 pointer-events-none" : "opacity-100" )}>
          <div>
            <CardContent className="space-y-6 pt-6">
               {/* selector */}
               <Select value={selected} onValueChange={setSelected} disabled={!!selectedChampionForModal || isLoadingLineage}>
                  <SelectTrigger className="w-full md:w-[300px] mx-auto">
                    <SelectValue placeholder={isLoadingLineage ? "Loading divisions…" : "Select a Weight Class"} />
                  </SelectTrigger>
                  <SelectContent> {Object.entries(classes).map(([key, val]) => ( <SelectItem key={key} value={key}> {val.displayName} </SelectItem> ))} </SelectContent>
               </Select>

                {/* champions timeline (Card View) */}
                {!data ? (
                   <div className="mt-6">
                     <div className="w-full overflow-x-auto pb-4">
                       <div className="flex space-x-3 items-stretch min-w-max">
                         {[...Array(8)].map((_, i) => (
                           <div key={i} className="flex-shrink-0 w-52 p-3 bg-muted rounded-md flex flex-col">
                             {/* image area skeleton — matches 3/4 aspect ratio */}
                             <Skeleton className="w-full rounded mb-2" style={{ aspectRatio: '3/4' }} />
                             {/* name line */}
                             <Skeleton className="h-4 w-3/4 rounded mb-1" />
                             {/* date line */}
                             <Skeleton className="h-3 w-full rounded mb-1" />
                             {/* duration line */}
                             <Skeleton className="h-3 w-1/2 rounded" />
                           </div>
                         ))}
                       </div>
                     </div>
                   </div>
                ) : (
                   <div className="mt-6">
                     {/* ---- Apply Drag Scroll Handlers Here ---- */}
                     <div
                       ref={scrollContainerRef}
                       className={cn(
                           "w-full overflow-x-auto overflow-y-visible pb-4",
                           !isTimelineView && !selectedChampionForModal && "cursor-grab"
                       )}
                       onPointerDown={handlePointerDown}
                       onPointerMove={handlePointerMove}
                       onPointerUp={handlePointerUpOrLeave}
                       onPointerLeave={handlePointerUpOrLeave}
                       onWheel={handleCardWheel}
                     >
                       <div className="flex space-x-3 items-stretch min-h-[200px] mx-auto min-w-max overflow-visible py-3 pr-6">
                         {cards.map((c, i) => {
                            const start = parseDate(c.reignStart)!;
                           const end = parseDate(c.reignEnd)!;
                           const reignDays = diffInDays(start, end);
                           const duration = humanDuration(reignDays);
                           const next = cards[i + 1];
                           let vacancyBlock: React.ReactNode = null;
                            if (next) {
                              const gap = diffInDays(end, parseDate(next.reignStart) || new Date()) - 1;
                              if (gap > 0) {
                                vacancyBlock = ( <div className="flex flex-col items-center justify-center px-2 text-center"> <div className="w-px h-16 bg-destructive/50 my-1"></div> <span className="text-xs text-destructive whitespace-normal w-16"> VACANT<br/>{humanDuration(gap)} </span> <div className="w-px h-16 bg-destructive/50 my-1"></div> </div> );
                              }
                            } else if (c.reignEnd && c.reignEnd !== "Present" && !c.notes?.toLowerCase().includes("interim")) {
                               vacancyBlock = ( <div className="flex flex-col items-center justify-center px-2 text-center"> <div className="w-px h-16 bg-destructive/50 my-1"></div> <span className="text-xs text-destructive whitespace-normal w-16"> VACANT </span> <div className="w-px h-16 bg-destructive/50 my-1"></div> </div> );
                            }

                           return (
                             <React.Fragment key={`${c.name}-${i}`}>
                               {/* ---- Add draggable=false and select-none to the card div ---- */}
                               <motion.div
                                 draggable={false}
                                 onClick={() => {
                                   if (!didDragRef.current) {
                                     setSelectedChampionForModal({ name: c.name, reignStart: c.reignStart, reignEnd: c.reignEnd, notes: c.notes, imageUrl: c.imageUrl, eloHistory: eloLookup(c.name) });
                                   }
                                 }}
                                 className="
                                   flex-shrink-0 w-52 p-3 bg-muted shadow-sm rounded-xl
                                   relative z-0 hover:z-10 flex flex-col
                                   select-none cursor-pointer group"
                                 whileHover={{ scale: 1.05, zIndex: 10 }}
                                 transition={{ duration: 0.2 }}
                               >
                                 <div
                                   draggable={false}
                                   className="relative w-full mb-2 bg-transparent rounded overflow-hidden flex items-center justify-center"
                                   style={{ aspectRatio: '3/4' }}
                                 >
                                   {/* Fallback icon only shown when no image */}
                                   {!c.imageUrl && <User className="h-12 w-12 text-gray-500 dark:text-gray-400" />}
                                   {c.imageUrl && (
                                     // eslint-disable-next-line @next/next/no-img-element
                                     <img
                                       src={c.imageUrl}
                                       alt={c.name}
                                       draggable={false}
                                       className="absolute inset-0 h-full w-full object-contain select-none"
                                       onError={(e) => { (e.currentTarget as HTMLImageElement).style.display = 'none'; (e.currentTarget.parentElement?.querySelector('.fallback-icon') as HTMLElement | null)?.style.setProperty('display', 'flex') }}
                                     />
                                   )}
                                 </div>
                                 {/* --- Text Block (Should be visible now) --- */}
                                 <div className="mt-1"> {/* <<<--- Added small top margin */}
                                     <p className="font-bold text-base leading-tight mb-1 truncate">{c.name}</p>
                                     <p className="text-xs"> {start.toLocaleDateString()} - {c.reignEnd === "Present" ? "Present" : end.toLocaleDateString()} </p>
                                     <p className="text-xs text-muted-foreground"> {duration} </p>
                                 </div>
                               </motion.div>
                               {vacancyBlock}
                             </React.Fragment>
                           );
                         })}
                       </div>
                     </div>
                     {/* ---- End Drag Scroll Container ---- */}
                   </div>
                 )}

             
            </CardContent>
          </div>
        </Card>)
      ) : (
        /* --- Timeline View --- */
        ((// Render actual timeline
      ((isTimelineCalculating || !timelineData) ? (<div className="w-full max-w-screen-2xl mx-auto border rounded-md bg-card p-4"> <div className="flex space-x-2"> <Skeleton className="h-[75vh] w-[2.5rem]" /> <div className="flex-1 space-y-2"> {[...Array(Object.keys(filteredClasses).length)].map((_, i) => ( <Skeleton key={i} className="h-[3.5rem] w-full" /> ))} <Skeleton className="h-[1.5rem] w-full mt-2" /> </div> </div> </div>) : (<div
          ref={timelineScrollContainerRef}
          onPointerDown={handleTimelinePointerDown}
          onPointerMove={handleTimelinePointerMove}
          onPointerUp={handleTimelinePointerUpOrLeave}
          onPointerLeave={handleTimelinePointerUpOrLeave}
          onWheel={handleTimelineWheel}
          className={cn(
            "w-full max-w-screen-2xl mx-auto overflow-x-auto border rounded-md bg-card text-card-foreground flex-1 overflow-y-hidden",
            isTimelineView && !selectedChampionForModal && "cursor-grab"
          )}
        >
        {/* Container for scrolling content (inner div, no handlers needed here) */}
        <div
          className="relative px-0 py-0 overflow-visible h-full select-none" // Add select-none here too
          style={{ width: `${timelineData.totalWidth + (2.5 * 16)}px`, paddingBottom: '1rem' }}
        >
          {/* Grid for rows, sidebar, and markers (no handlers needed here) */}
          <div
            className="relative z-10 overflow-visible grid h-full"
            style={{ gridTemplateRows: `1rem repeat(${Object.keys(filteredClasses).length}, 1fr) 1.5rem`, gridTemplateColumns: '2.5rem 1fr', width: '100%' }}
          >
             {/* Top Left Corner */}
             <div className="sticky left-0 bg-card z-30 border-r border-border/60" style={{ gridRow: 1, gridColumn: 1 }}></div>
             {/* Top Markers */}
             <div className="relative h-full border-b border-border/60" style={{ gridRow: 1, gridColumn: 2 }}>
              {timelineData?.monthMarkers.map(({ date, x }) => { const isYear = date.getMonth() === 0; return ( <div key={`top-month-${date.toISOString()}`} className="absolute bottom-0 flex flex-col items-center" style={{ left: `${x}px`, height: isYear ? '0.75rem' : '0.375rem' }}> <div className="w-px h-full bg-foreground/50"></div> </div> ); })}
              {timelineData && timelineData.todayX >= 0 && ( <div key="top-today-marker" className="absolute bottom-0 z-10" style={{ left: `${timelineData.todayX}px`, height: '1rem' }}> <div className="w-px h-full bg-primary/70"></div> </div> )}
            </div>
            {/* Weight Classes & Bars */}
            {Object.entries(filteredClasses).map(([key,wc],r)=>(
              <div key={key} className="contents">
                {/* abbrev name (add select-none) */}
                <div
                  className="sticky left-0 bg-card z-30 flex items-center justify-center px-1 border-r border-b border-border/60 select-none"
                  style={{ gridRow: r + 2, gridColumn: 1, width: '2.5rem' }}
                >
                  <TooltipProvider>
                    <Tooltip delayDuration={300}>
                      <TooltipTrigger asChild>
                        <span
                          className="text-[12px] font-bold uppercase text-center cursor-help px-0.5 rounded"
                          style={wcBadgeStyle(key, { bgOpacity: 0 })}
                        >
                          {abbreviateWeightClass(key)}
                        </span>
                      </TooltipTrigger>
                      <TooltipContent side="right">
                        {getWeightClassName(key)}
                      </TooltipContent>
                    </Tooltip>
                  </TooltipProvider>
                </div>
                {/* champion bars container (no handlers needed here) */}
                <div
                  className="relative h-full border-b border-border/60 overflow-visible"
                  style={{ gridRow: r + 2, gridColumn: 2, position: 'relative' }}
                >
                  {wc.champions.map((c,i)=>{
                    const s = parseDate(c.reignStart)!;
                    const e = parseDate(c.reignEnd)!;
                    const cs = s < timelineData!.minDate ? timelineData!.minDate : s;
                    const ce = e > timelineData!.maxDate ? timelineData!.maxDate : e;
                    const left = diffInDays(timelineData!.minDate,cs)*timelineData!.pixelsPerDay;
                    const width = diffInDays(cs,ce)*timelineData!.pixelsPerDay;
                    const reignYears = diffInDays(s, e) / 365;
                    const isLongReign = reignYears >= 1;

                    // <<<--- REMOVE old hue calculation (commented out for reference)
                    // let h=0; for(let k=0;k<c.name.length;k++) h=c.name.charCodeAt(k)+((h<<5)-h);
                    // const hue=h%360, bg=`hsl(${hue},70%,90%)`, hov=`hsl(${hue},70%,80%)`, col=`hsl(${hue},80%,25%)`;

                    // <<<--- ADD new color calculation:
                    const reignDaysTotal = diffInDays(s, e); // Calculate total reign days
                    const { bg, text, hoverBg } = getReignColor(reignDaysTotal);

                    return (
                      <motion.div
                        key={`${key}-${i}`}
                        draggable={false}
                        onClick={() => { if (!timelineDidDragRef.current) { setSelectedChampionForModal({ name: c.name, reignStart: c.reignStart, reignEnd: c.reignEnd, notes: c.notes, imageUrl: c.imageUrl, eloHistory: eloLookup(c.name) }); } }}
                        className="absolute top-[2px] bottom-[2px] flex items-center justify-center group cursor-pointer shadow-sm z-10 select-none"
                        style={{
                          left:`${left}px`,
                          width:`${Math.max(width,2)}px`,
                          backgroundColor: bg,
                          color: text,
                          borderRadius: '0.375rem'
                        }}
                        initial={{ scale:1 }}
                        whileHover={{
                          backgroundColor: hoverBg,
                          zIndex:20,
                          scale: 1.03,
                          borderRadius: '0.375rem'
                        }}
                        transition={{ duration:0.15 }}
                      >
                        {/* Inner spans will inherit text color */}
                        {isLongReign ? (
                          <span className="text-[20px] font-bold whitespace-nowrap px-1 overflow-hidden text-ellipsis max-w-full select-none">
                            {c.name}
                          </span>
                        ) : (
                          <span
                            className="text-[20px] font-bold whitespace-nowrap overflow-visible group-hover:opacity-0 transition-opacity duration-100 select-none"
                            style={{ position: 'relative', zIndex: 15 }}
                          >
                            {getInitials(c.name)}
                          </span>
                        )}
                        {!isLongReign && (
                          <span
                            className="absolute inset-0 flex items-center justify-center text-center text-[15px] font-semibold px-1 opacity-0 group-hover:opacity-100 transition-opacity duration-100 delay-50 pointer-events-none select-none"
                            style={{
                                background: hoverBg,
                                color: text,
                                borderRadius: '0.375rem',
                                zIndex: 21
                            }}
                          >
                            {c.name}
                          </span>
                        )}
                      </motion.div>
                    );
                  })}
                </div>
              </div>
            ))}

             {/* ... Bottom Left Corner ... */}
             {/* Bottom Markers (add select-none to text spans) */}
             <div className="relative h-full" style={{ gridRow: Object.keys(filteredClasses).length + 2, gridColumn: 2 }}>
               {timelineData?.monthMarkers.map(({ date, x }) => { const isYear = date.getMonth() === 0; return ( <div key={`bottom-month-${date.toISOString()}`} className="absolute top-0 flex flex-col items-center" style={{ left: `${x}px`, height: isYear ? '0.75rem' : '0.375rem' }}> <div className="w-px h-full bg-foreground/50"></div> {isYear && ( <span className="absolute top-full text-[12px] text-foreground transform -translate-x-1/2 select-none"> {date.getFullYear()} </span> )} </div> ); })}
               {timelineData && timelineData.todayX >= 0 && ( <div key="bottom-today-marker" className="absolute top-0 z-10 flex flex-col items-center" style={{ left: `${timelineData.todayX}px`, height: '1rem' }}> <div className="w-px h-full bg-primary/70"></div> <span className="absolute top-full text-[9px] text-primary font-semibold transform -translate-x-1/2 bg-card px-0.5 select-none"> Now </span> </div> )}
             </div>

          </div> {/* End of Grid */}
        </div> {/* End of scrollable content container */}
      </div>))) /* End of Timeline View */)
      )}
      {/* Fighter Details Modal */}
      {selectedChampionForModal && (
        <FighterDetailsModal
          champion={selectedChampionForModal}
          onClose={() => setSelectedChampionForModal(null)}
          eloLookup={eloLookup}
          hideElo={true}
        />
      )}
    </main>
  );
};

export default UfcChampionsDisplay;