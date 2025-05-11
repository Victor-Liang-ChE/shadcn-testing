"use client";

import React, { useState, useEffect, useRef, useMemo } from "react";
import Image from "next/image";
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
} from "lucide-react"; // Updated imports
import { cn } from "@/lib/utils";
import { motion, AnimatePresence } from "framer-motion";
import { Skeleton } from "@/components/ui/skeleton"; // Import Skeleton
import { ufcChampionsData, type Champion, type WeightClassData} from "@/data/ufcChampionsData";

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
        bantamweight: "BW", flyweight: "FlW", women_strawweight: "WSW",
        women_flyweight: "WFlW", women_bantamweight: "WBW",
    };
    return abbreviations[key.toLowerCase()] || key.toUpperCase();
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
const UfcChampionsDisplay: React.FC = () => {
  // --- STATE DECLARATIONS ---
  const [selected, setSelected] = useState("lightweight");
  const [showDetailsView, setShowDetailsView] = useState(false);
  const [selectedFighter, setSelectedFighter] = useState<Champion | null>(null);
  const [fighterDetails, setFighterDetails] = useState<FighterDetails | null>(null);
  const [fighterImages, setFighterImages] = useState<UfcChampionImages | null>(null);
  const [fighterHistory, setFighterHistory] = useState<ProcessedFightResult[]>([]);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [isTimelineView, setIsTimelineView] = useState(false);
  const [isTimelineCalculating, setIsTimelineCalculating] = useState(false);
  const [genderFilter, setGenderFilter] = useState<"all" | "men" | "women">("all");
  const dataCache = useRef<Record<string, { details: FighterDetails | null; images: UfcChampionImages | null; history: ProcessedFightResult[] }>>({});
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

  const classes = ufcChampionsData;
  const data = classes[selected];

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

  // Effect for card view scroll to end - ONLY runs when 'selected' changes
  useEffect(() => {
    // Check if the card view is actually active and not dragging
    if (scrollContainerRef.current && data && !isTimelineView && !showDetailsView && !isDraggingRef.current) {
      // Scroll to the end when a new weight class is selected
      console.log("Scrolling to end due to selection change:", selected);
      scrollContainerRef.current.scrollLeft = scrollContainerRef.current.scrollWidth;
    }
    // Depend ONLY on 'selected' and 'data' to trigger this scroll-to-end
  }, [selected, data]); // <<<--- Changed Dependency Array

  // Effect for timeline view scroll to end
  useEffect(() => {
    if (isTimelineView && timelineScrollContainerRef.current) {
      timelineScrollContainerRef.current.scrollLeft = timelineScrollContainerRef.current.scrollWidth;
    }
  }, [isTimelineView]);

  /* ================================================================== */
  /*  FETCH DATA                                                        */
  /* ================================================================== */
  const fetchFighterData = async (fighter: Champion) => {
    if (!fighter) return;
    const original = fighter.name.replace(/ /g, " ").trim();
    const normName = normalizeName(original);
    console.log(`Fetching data for: ${original} (Normalized: ${normName})`);

    if (dataCache.current[normName]) {
      console.log("Cache hit for:", normName);
      const cached = dataCache.current[normName];
      setFighterDetails(cached.details);
      setFighterImages(cached.images);
      setFighterHistory(cached.history);
      setSelectedFighter(fighter);
      setShowDetailsView(true);
      setIsLoading(false);
      setError(null);
      return;
    }
    console.log("Cache miss for:", normName);

    setIsLoading(true);
    setError(null);
    setFighterDetails(null);
    setFighterImages(null);
    setFighterHistory([]);
    setSelectedFighter(fighter);
    setShowDetailsView(true);

    try {
      const { firstName, lastName } = parseFighterName(original);
      console.log(`Querying images with FIRST: %${firstName}%, LAST: %${lastName}%`);

      const [detailsRes, imagesRes, fightsRes] = await Promise.all([
        supabase.from("ufc_fighter_tott").select("*").ilike("FIGHTER", `%${normName}%`).limit(1).maybeSingle(),
        supabase.from("ufc_champion_images").select("*").ilike("FIRST", `%${firstName}%`).ilike("LAST", `%${lastName}%`).limit(1).maybeSingle(),
        supabase.from("ufc_fight_results").select("*").ilike("BOUT", `%${normName}%`).order("id", { ascending: true }),
      ]);

      if (detailsRes.error) console.error("Supabase details fetch error:", detailsRes.error);
      const fighterDetailsFetched = detailsRes.data as FighterDetails | null;
      console.log("Fetched Details:", fighterDetailsFetched);

      if (imagesRes.error) console.error("Supabase images fetch error:", imagesRes.error);
      const fighterImagesFetched = imagesRes.data as UfcChampionImages | null;
      console.log("Fetched Images:", fighterImagesFetched);

      if (fightsRes.error) throw fightsRes.error;
      const fights = (fightsRes.data || []) as FightResult[];
      console.log("Fetched Raw Fight Results (Sorted):", fights);

      const processed: ProcessedFightResult[] = [];
      const championStatus: Record<string, boolean> = {};

      for (let idx = fights.length - 1; idx >= 0; idx--) {
        const fight = fights[idx];
        const resultForFighter = getResultForFighter(fight, normName);
        const div = normalizeWeightClass(fight.WEIGHTCLASS);
        const isTitle = fight.WEIGHTCLASS?.toLowerCase().includes("title") ?? false;
        let defenseOutcome: "success" | "fail" | null = null;

        if (isTitle && div) {
          if (championStatus[div]) {
            if (resultForFighter === "loss") {
              defenseOutcome = "fail"; championStatus[div] = false;
            } else if (resultForFighter === "win") {
              defenseOutcome = "success";
            }
          } else {
            if (resultForFighter === "win") championStatus[div] = true;
          }
        }

        let decisionDetails: string | undefined;
        let finishDetails: string | undefined;
        if (fight.METHOD?.toLowerCase().includes("decision")) {
          decisionDetails = extractDecisionDetails(fight.DETAILS);
        } else if (fight.DETAILS?.trim()) {
          if (!/\d{1,2}\s*-\s*\d{1,2}/.test(fight.DETAILS)) finishDetails = fight.DETAILS.trim();
          else decisionDetails = extractDecisionDetails(fight.DETAILS);
        }

        processed.unshift({
          ...fight,
          OPPONENT: extractOpponent(fight.BOUT, normName),
          resultForFighter, decisionDetails, DETAILS: finishDetails, defenseOutcome,
        });
      }
      console.log("Processed History:", processed);

      dataCache.current[normName] = { details: fighterDetailsFetched, images: fighterImagesFetched, history: processed };
      setFighterDetails(fighterDetailsFetched);
      setFighterImages(fighterImagesFetched);
      setFighterHistory(processed);
    } catch (err: any) {
      console.error("Error fetching fighter data:", err);
      setError(err.message ?? "Unknown error");
    } finally {
      setIsLoading(false);
    }
  };

  const handleCloseDetails = () => {
    setShowDetailsView(false);
    setError(null);
  };

  const getLayoutId = (champion: Champion | null) => {
      if (!champion) return undefined;
      return `champion-card-${normalizeName(champion.name)}-${champion.reignStart}`;
  }

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
    if (isTimelineView || showDetailsView || !scrollContainerRef.current) return;
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
    if (!isDraggingRef.current || isTimelineView || showDetailsView || !scrollContainerRef.current) return;
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
    if (!isTimelineView || showDetailsView || !timelineScrollContainerRef.current) return;
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
    if (!isTimelineDraggingRef.current || !isTimelineView || showDetailsView || !timelineScrollContainerRef.current) return;
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


  /* ------------------------------------------------------------------ */
  /*  UI                                                                */
  /* ------------------------------------------------------------------ */
  return (
    <main className="relative pt-8 md:pt-6 px-4 md:px-6 pb-20 w-full flex flex-col h-screen">
      {/* --- View Toggle and Gender Filter --- */}
      <div className="w-full max-w-screen-2xl mx-auto flex justify-between items-start mb-4 px-0 md:px-0 min-h-[40px]">
        <div className="flex-1">
          {isTimelineView && (
            <div className="flex items-center space-x-2">
              <Select value={genderFilter} onValueChange={(value: "all" | "men" | "women") => setGenderFilter(value)} disabled={showDetailsView}>
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
          <Checkbox id="timeline-view-toggle" checked={isTimelineView} onCheckedChange={(checked) => { const isChecked = Boolean(checked); setIsTimelineView(isChecked); setGenderFilter(isChecked ? "men" : "all"); }} disabled={showDetailsView} />
          <Label htmlFor="timeline-view-toggle" className="text-sm font-medium">Show Full Timeline View</Label>
        </div>
      </div>
      {/* --- Conditional Rendering: Card View or Timeline View --- */}
      {!isTimelineView ? (
        /* --- Main Content Card (Background) --- */
        (<Card className={cn( "w-full max-w-screen-2xl mx-auto transition-opacity duration-300 flex flex-col mb-auto", showDetailsView ? "opacity-20 pointer-events-none" : "opacity-100" )}>
          <div className="overflow-y-auto">
            <CardContent className="space-y-6 pt-6">
               {/* selector */}
               <Select value={selected} onValueChange={setSelected} disabled={showDetailsView}>
                  <SelectTrigger className="w-full md:w-[300px] mx-auto"> <SelectValue placeholder="Select a Weight Class" /> </SelectTrigger>
                  <SelectContent> {Object.entries(classes).map(([key, val]) => ( <SelectItem key={key} value={key}> {val.displayName} </SelectItem> ))} </SelectContent>
               </Select>

                {/* champions timeline (Card View) */}
                {!data ? (
                   <div className="mt-6">
                     <div className="w-full overflow-x-auto pb-4">
                       <div className="flex space-x-3 items-stretch min-h-[200px]"> {[...Array(5)].map((_, i) => ( <Skeleton key={i} className="flex-shrink-0 w-52 h-[180px] rounded-md" /> ))} </div>
                     </div>
                   </div>
                ) : (
                   <div className="mt-6">
                     {/* ---- Apply Drag Scroll Handlers Here ---- */}
                     <div
                       ref={scrollContainerRef}
                       className={cn(
                           "w-full overflow-x-auto pb-4",
                           !isTimelineView && !showDetailsView && "cursor-grab" // Add grab cursor
                       )}
                       onPointerDown={handlePointerDown}
                       onPointerMove={handlePointerMove}
                       onPointerUp={handlePointerUpOrLeave}
                       onPointerLeave={handlePointerUpOrLeave} // Handle leave while dragging
                     >
                       <div className="flex space-x-3 items-stretch min-h-[200px]">
                         {data.champions.map((c, i) => {
                           const start = parseDate(c.reignStart)!;
                           const end = parseDate(c.reignEnd)!;
                           const reignDays = diffInDays(start, end);
                           const duration = humanDuration(reignDays);
                           const next = data.champions[i + 1];
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
                                 layoutId={getLayoutId(c)}
                                 draggable={false} // Keep this on the motion.div
                                 onClick={() => {
                                   // <<<--- Check didDragRef before fetching
                                   if (!didDragRef.current) {
                                       fetchFighterData(c);
                                   }
                                   // didDragRef will be reset on the *next* pointer down
                                 }}
                                 className="
                                   flex-shrink-0 w-52 p-3 bg-muted rounded-md shadow-sm
                                   relative z-0 hover:z-10 flex flex-col
                                   overflow-hidden {/* Keep overflow-hidden */}
                                   select-none cursor-pointer group" // Added cursor-pointer back
                                 whileHover={{ scale: 1.05, zIndex: 10 }}
                                 transition={{ duration: 0.2 }}
                               >
                                 {/* ---- Add draggable=false ONLY to the container div ---- */}
                                 <div
                                   draggable={false} // Keep this on the container div
                                   className="relative h-24 w-full mb-2 bg-gradient-to-b from-gray-300 to-gray-400 dark:from-gray-600 dark:to-gray-700 rounded flex items-center justify-center"
                                 >
                                   {/* ---- REMOVE draggable={false} from User ---- */}
                                   <User className="h-12 w-12 text-gray-500 dark:text-gray-400" />
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

             <p className="text-xs text-center text-muted-foreground mt-8"> Missing UFC 294 and 308 Events. </p>
            </CardContent>
          </div>
        </Card>)
      ) : (
        /* --- Timeline View --- */
        ((isTimelineCalculating || !timelineData) ? // Render actual timeline
        (<div className="w-full max-w-screen-2xl mx-auto border rounded-md bg-card p-4"> <div className="flex space-x-2"> <Skeleton className="h-[75vh] w-[2.5rem]" /> <div className="flex-1 space-y-2"> {[...Array(Object.keys(filteredClasses).length)].map((_, i) => ( <Skeleton key={i} className="h-[3.5rem] w-full" /> ))} <Skeleton className="h-[1.5rem] w-full mt-2" /> </div> </div> </div>) : (<div
          ref={timelineScrollContainerRef}
          // Attach timeline-specific handlers
          onPointerDown={handleTimelinePointerDown}
          onPointerMove={handleTimelinePointerMove}
          onPointerUp={handleTimelinePointerUpOrLeave}
          onPointerLeave={handleTimelinePointerUpOrLeave}
          className={cn( // Add conditional cursor style
            "w-full max-w-screen-2xl mx-auto overflow-x-auto border rounded-md bg-card text-card-foreground flex-1 min-h-[70vh]",
            isTimelineView && !showDetailsView && "cursor-grab" // Grab cursor for timeline
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
                    <span className="text-[12px] font-bold uppercase text-center">
                      {abbreviateWeightClass(key)}
                    </span>
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
                          layoutId={getLayoutId(c)} // <<<--- ADD layoutId HERE
                          draggable={false}
                          onClick={() => { if (!timelineDidDragRef.current) { fetchFighterData(c) } }}
                          className="absolute top-[2px] bottom-[2px] rounded-sm flex items-center justify-center group cursor-pointer shadow-sm z-10 select-none"
                          // <<<--- UPDATE style prop
                          style={{
                            left:`${left}px`,
                            width:`${Math.max(width,2)}px`,
                            backgroundColor: bg, // Use calculated background
                            color: text,       // Use calculated text color
                            borderRadius: '0.125rem'
                          }}
                          initial={{ scale:1 }}
                          // <<<--- UPDATE whileHover prop
                          whileHover={{
                            backgroundColor: hoverBg, // Use calculated hover background
                            zIndex:20,
                            scale: 1.03,
                            borderRadius: '0.125rem'
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
                              className="absolute inset-0 flex items-center justify-center text-[15px] font-semibold px-1 opacity-0 group-hover:opacity-100 transition-opacity duration-100 delay-50 pointer-events-none select-none"
                              // <<<--- UPDATE hover span style
                              style={{
                                  background: hoverBg, // Use hover background
                                  color: text,        // Use base text color (should still contrast)
                                  borderRadius: '0.125rem',
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
        </div>) /* End of Timeline View */)
      )}
      {/* --- Fighter Details Overlay --- */}
      <AnimatePresence>
        {showDetailsView && selectedFighter && (
          <motion.div key="overlay-backdrop" initial={{ opacity: 0 }} animate={{ opacity: 1 }} exit={{ opacity: 0 }} transition={{ duration: 0.3 }} className="fixed inset-0 z-40 flex items-center justify-center bg-black/60 backdrop-blur-sm p-4 md:p-8" onClick={handleCloseDetails}>
            <motion.div layoutId={getLayoutId(selectedFighter)} className="bg-card text-card-foreground rounded-md shadow-xl w-full max-w-6xl h-[90vh] max-h-[800px] flex flex-col relative overflow-hidden z-50" onClick={(e) => e.stopPropagation()}>
              <motion.button initial={{ opacity: 0, scale: 0.5 }} animate={{ opacity: 1, scale: 1, transition: { delay: 0.2 } }} exit={{ opacity: 0 }} onClick={handleCloseDetails} className="absolute top-3 right-3 text-muted-foreground hover:text-foreground z-50 p-1 rounded-full bg-card hover:bg-muted" aria-label="Close details"> <X size={20} /> </motion.button>
              <div className="flex flex-1 h-full overflow-hidden pt-10">
                  {isLoading && ( <div className="absolute inset-0 flex items-center justify-center bg-card/80 z-20"> <Loader2 className="h-8 w-8 animate-spin text-primary" /> </div> )}
                  {error && !isLoading && ( <div className="absolute inset-0 flex flex-col items-center justify-center bg-card/80 z-20 p-4 text-center"> <p className="text-destructive font-semibold mb-2">Error</p> <p className="text-sm text-destructive">{error}</p> <button onClick={handleCloseDetails} className="mt-4 text-sm underline">Close</button> </div> )}
                  {!isLoading && !error && (
                    <>
                      {/* Left Panel: Stats */}
                      <div className="w-1/3 border-r border-border p-4 md:p-6 overflow-y-auto flex flex-col items-center">
                         <div className="relative h-40 w-40 mb-4 bg-gradient-to-b from-gray-300 to-gray-400 dark:from-gray-600 dark:to-gray-700 rounded-full flex items-center justify-center overflow-hidden">
                            {fighterImages?.HEADSHOT ? ( <Image src={fighterImages.HEADSHOT} alt={selectedFighter.name} layout="fill" objectFit="cover" /> ) : ( <User className="h-24 w-24 text-gray-500 dark:text-gray-400" /> )}
                         </div>
                         <h2 className="text-2xl font-bold mb-4 text-center">{fighterDetails?.FIGHTER || selectedFighter.name}</h2>
                         {fighterDetails ? (
                           <div className="space-y-2 text-sm w-full text-center">
                             {fighterDetails.HEIGHT && <p><strong>Height:</strong> {fighterDetails.HEIGHT}</p>} {fighterDetails.WEIGHT && <p><strong>Weight:</strong> {fighterDetails.WEIGHT}</p>} {fighterDetails.REACH && <p><strong>Reach:</strong> {fighterDetails.REACH}</p>} {fighterDetails.STANCE && <p><strong>Stance:</strong> {fighterDetails.STANCE}</p>}
                             {fighterDetails.DOB && ( <p> <strong>Born:</strong> {new Date(fighterDetails.DOB).toLocaleDateString()} {calculateAge(fighterDetails.DOB) !== null && ` (Age: ${calculateAge(fighterDetails.DOB)})`} </p> )}
                           </div>
                         ) : ( <p className="text-sm text-muted-foreground mt-4">No detailed stats available.</p> )}
                      </div>
                      {/* Right Panel: Fight History */}
                      <div className="w-2/3 p-4 md:p-6 overflow-y-auto">
                        {fighterHistory.length > 0 ? (
                          <ul className="space-y-1">
                            {fighterHistory.map((fight) => {
                              const isTitleFight = fight.WEIGHTCLASS?.toLowerCase().includes('title') ?? false;

                              // <<<--- Calculate Weight Class Abbreviation ---<<<
                              const normalizedKey = normalizeWeightClass(fight.WEIGHTCLASS);
                              const abbreviation = normalizedKey ? abbreviateWeightClass(normalizedKey) : null;
                              // <<<------------------------------------------<<<

                              return (
                                <li
                                  key={fight.id}
                                  // Make li relative to position the absolute span inside it
                                  className={cn(
                                    "relative text-sm border-b border-border p-2 rounded-lg transition-colors duration-150 overflow-visible mb-2",
                                    (fight.resultForFighter === 'draw' || fight.resultForFighter === 'nc') && 'bg-gray-100 dark:bg-gray-800/30 border-gray-300 dark:border-gray-600',
                                    fight.resultForFighter === 'win' && 'bg-green-100 dark:bg-green-900/30 border-green-300 dark:border-green-700',
                                    fight.resultForFighter === 'loss' && 'bg-red-100 dark:bg-red-900/30 border-red-300 dark:border-red-700',
                                    fight.resultForFighter === 'unknown' && 'bg-muted/50',
                                    "hover:shadow-md cursor-pointer"
                                  )}
                                  onClick={() => { if (fight.URL) window.open(fight.URL, '_blank', 'noopener,noreferrer'); }}
                                >
                                  {/* Icons (Crown, ShieldCheck, ShieldOff) - Unchanged */
                                  isTitleFight && !fight.defenseOutcome && ( <Crown size={16} className="absolute -top-1 -left-2 text-yellow-500 dark:text-yellow-400 opacity-70 -rotate-45 z-10" aria-label="Title Fight" /> )
                                  }
                                  {fight.defenseOutcome === 'success' && ( <ShieldCheck size={18} className="absolute -top-1 -left-2 text-blue-500 dark:text-blue-400 opacity-100 -rotate-45 z-10" aria-label="Successful Title Defence" /> )}
                                  {fight.defenseOutcome === 'fail' && ( <ShieldOff size={18} className="absolute -top-1 -left-2 text-red-500 dark:text-red-400 opacity-100 -rotate-45 z-10" aria-label="Lost Title Defence" /> )}

                                  {/* Result Badge and Opponent - Added mb-1 */}
                                  <div className="flex justify-between items-center mb-1"> {/* Added mb-1 for spacing */}
                                    <div>
                                      <span className={cn( "font-semibold mr-2 px-1.5 py-0.5 rounded text-xs uppercase", fight.resultForFighter === 'win' && 'bg-green-600 text-white', fight.resultForFighter === 'loss' && 'bg-red-600 text-white', (fight.resultForFighter === 'draw' || fight.resultForFighter === 'nc') && 'bg-gray-500 text-white', fight.resultForFighter === 'unknown' && 'bg-gray-400 text-white' )}> {fight.resultForFighter === 'draw' ? 'DRAW' : fight.resultForFighter === 'nc' ? 'NC' : fight.resultForFighter} </span>
                                      <span>vs. </span> <a href={createAthleteUrl(fight.OPPONENT || '')} target="_blank" rel="noopener noreferrer" className="font-bold hover:underline text-primary" onClick={(e) => e.stopPropagation()}> {fight.OPPONENT || 'Opponent Unknown'} </a>
                                    </div>
                                    <span className="text-xs text-muted-foreground">{fight.EVENT}</span>
                                  </div>

                                  {/* Method/Round/Time/Details - Added pr-8 */}
                                   <div className="pr-8"> {/* Added right padding to prevent overlap */}
                                    <p className="text-xs text-muted-foreground mt-0.5"> {fight.METHOD} {fight.ROUND && ` | Rd: ${fight.ROUND}`} {fight.TIME && ` | Time: ${fight.TIME}`} </p>
                                    {fight.decisionDetails && ( <p className="text-xs text-muted-foreground italic mt-0.5"> {fight.decisionDetails} </p> )}
                                    {fight.DETAILS && ( <p className="text-xs text-muted-foreground mt-0.5"> {fight.DETAILS} </p> )}
                                  </div>

                                  {/* <<<--- Add Weight Class Abbreviation Span ---<<< */}
                                  {abbreviation && (
                                    <span
                                      className="absolute bottom-[3px] right-[5px] text-[10px] font-semibold uppercase text-muted-foreground/80 bg-primary/5 px-1 rounded-[3px]"
                                      aria-label={`Weight Class: ${fight.WEIGHTCLASS || abbreviation}`} // More descriptive label
                                    >
                                      {abbreviation}
                                    </span>
                                  )}
                                  {/* <<<------------------------------------------<<< */}

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
        )}
      </AnimatePresence>
    </main>
  );
};

export default UfcChampionsDisplay;