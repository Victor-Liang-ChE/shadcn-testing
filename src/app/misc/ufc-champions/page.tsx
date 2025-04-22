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

/* ------------------------------------------------------------------ */
/*  TYPES                                                             */
/* ------------------------------------------------------------------ */
interface Champion {
  name: string;
  reignStart: string;
  reignEnd: string | null;
  notes?: string;
}
interface WeightClassData {
  displayName: string;
  champions: Champion[];
}
interface ProcessedFightResult extends FightResult {
  resultForFighter: "win" | "loss" | "draw" | "nc" | "unknown";
  decisionDetails?: string;
  finishDetails?: string;
  defenseOutcome?: "success" | "fail" | null; // Updated type
}

/* ------------------------------------------------------------------ */
/*  RAW DATA                                                          */
/* ------------------------------------------------------------------ */
const ufcChampionsData: Record<string, WeightClassData> = {
  /* --------------------------- MEN -------------------------------- */
  heavyweight: {
    displayName: "Heavyweight (206‑265 lbs)",
    champions: [
      { name: "Mark Coleman",   reignStart: "1997-02-07", reignEnd: "1997-07-27" },
      { name: "Maurice Smith",  reignStart: "1997-07-27", reignEnd: "1997-12-21" },
      { name: "Randy Couture",  reignStart: "1997-12-21", reignEnd: "1998-01-15", notes: "stripped (contract)" },
      { name: "Bas Rutten",     reignStart: "1999-05-07", reignEnd: "1999-06-09", notes: "vacated (drop to LHW)" },
      { name: "Kevin Randleman",reignStart: "1999-11-19", reignEnd: "2000-11-17" },
      { name: "Randy Couture",  reignStart: "2000-11-17", reignEnd: "2002-03-22" },
      { name: "Josh Barnett",   reignStart: "2002-03-22", reignEnd: "2002-07-26", notes: "stripped (drug test)" },
      { name: "Ricco Rodriguez",reignStart: "2002-09-27", reignEnd: "2003-02-28" },
      { name: "Tim Sylvia",     reignStart: "2003-02-28", reignEnd: "2003-10-15", notes: "stripped (drug test)" },
      { name: "Frank Mir",      reignStart: "2004-06-19", reignEnd: "2005-08-12", notes: "stripped (injury)" },
      { name: "Andrei Arlovski",reignStart: "2005-08-12", reignEnd: "2006-04-15" },
      { name: "Tim Sylvia",     reignStart: "2006-04-15", reignEnd: "2007-03-03" },
      { name: "Randy Couture",  reignStart: "2007-03-03", reignEnd: "2008-11-15" },
      { name: "Brock Lesnar",   reignStart: "2008-11-15", reignEnd: "2010-10-23" },
      { name: "Cain Velasquez", reignStart: "2010-10-23", reignEnd: "2011-11-12" },
      { name: "Junior dos Santos",reignStart: "2011-11-12",reignEnd: "2012-12-29" },
      { name: "Cain Velasquez", reignStart: "2012-12-29", reignEnd: "2015-06-13" },
      { name: "Fabrício Werdum",reignStart: "2015-06-13", reignEnd: "2016-05-14" },
      { name: "Stipe Miocic",   reignStart: "2016-05-14", reignEnd: "2018-07-07" },
      { name: "Daniel Cormier", reignStart: "2018-07-07", reignEnd: "2019-08-17" },
      { name: "Stipe Miocic",   reignStart: "2019-08-17", reignEnd: "2021-03-27" },
      { name: "Francis Ngannou",reignStart: "2021-03-27", reignEnd: "2023-01-14", notes: "vacated (free agent)" },
      { name: "Jon Jones",      reignStart: "2023-03-04", reignEnd: "Present" },
    ],
  },

  lightheavyweight: {
    displayName: "Light Heavyweight (186‑205 lbs)",
    champions: [
      { name: "Frank Shamrock",   reignStart: "1997-12-21", reignEnd: "1999-11-24" },
      { name: "Tito Ortiz",       reignStart: "2000-04-14", reignEnd: "2003-09-26" },
      { name: "Randy Couture",    reignStart: "2003-09-26", reignEnd: "2004-01-31" },
      { name: "Vitor Belfort",    reignStart: "2004-01-31", reignEnd: "2004-08-21" },
      { name: "Randy Couture",    reignStart: "2004-08-21", reignEnd: "2005-04-16" },
      { name: "Chuck Liddell",    reignStart: "2005-04-16", reignEnd: "2007-05-26" },
      { name: "Quinton Jackson",  reignStart: "2007-05-26", reignEnd: "2008-07-05" },
      { name: "Forrest Griffin",  reignStart: "2008-07-05", reignEnd: "2008-12-27" },
      { name: "Rashad Evans",     reignStart: "2008-12-27", reignEnd: "2009-05-23" },
      { name: "Lyoto Machida",    reignStart: "2009-05-23", reignEnd: "2010-05-08" },
      { name: "Maurício Rua",     reignStart: "2010-05-08", reignEnd: "2011-03-19" },
      { name: "Jon Jones",        reignStart: "2011-03-19", reignEnd: "2015-04-28" },
      { name: "Daniel Cormier",   reignStart: "2015-05-23", reignEnd: "2018-12-28" },
      { name: "Jon Jones",        reignStart: "2018-12-29", reignEnd: "2020-08-17" },
      { name: "Jan Błachowicz",   reignStart: "2020-09-27", reignEnd: "2021-10-30" },
      { name: "Glover Teixeira",  reignStart: "2021-10-30", reignEnd: "2022-06-12" },
      { name: "Jiří Procházka",   reignStart: "2022-06-12", reignEnd: "2022-11-23", notes: "vacated (injury)" },
      { name: "Jamahal Hill",     reignStart: "2023-01-21", reignEnd: "2023-07-14", notes: "vacated (injury)" },
      { name: "Alex Pereira",     reignStart: "2023-11-11", reignEnd: "2025-03-08" },
      { name: "Magomed Ankalaev", reignStart: "2025-03-08", reignEnd: "Present" },
    ],
  },

  middleweight: {
    displayName: "Middleweight (171‑185 lbs)",
    champions: [
      { name: "Dave Menne",        reignStart: "2001-09-28", reignEnd: "2002-01-11" },
      { name: "Murilo Bustamante", reignStart: "2002-01-11", reignEnd: "2002-10-05" },
      { name: "Evan Tanner",       reignStart: "2005-02-05", reignEnd: "2005-06-04" },
      { name: "Rich Franklin",     reignStart: "2005-06-04", reignEnd: "2006-10-14" },
      { name: "Anderson Silva",    reignStart: "2006-10-14", reignEnd: "2013-07-06" },
      { name: "Chris Weidman",     reignStart: "2013-07-06", reignEnd: "2015-12-12" },
      { name: "Luke Rockhold",     reignStart: "2015-12-12", reignEnd: "2016-06-04" },
      { name: "Michael Bisping",   reignStart: "2016-06-04", reignEnd: "2017-11-04" },
      { name: "Georges St‑Pierre", reignStart: "2017-11-04", reignEnd: "2017-12-07", notes: "vacated (colitis)" },
      { name: "Robert Whittaker",  reignStart: "2017-12-07", reignEnd: "2019-10-06" },
      { name: "Israel Adesanya",   reignStart: "2019-10-06", reignEnd: "2022-11-12" },
      { name: "Alex Pereira",      reignStart: "2022-11-12", reignEnd: "2023-04-08" },
      { name: "Israel Adesanya",   reignStart: "2023-04-08", reignEnd: "2023-09-10" },
      { name: "Sean Strickland",   reignStart: "2023-09-10", reignEnd: "2024-01-20" },
      { name: "Dricus du Plessis", reignStart: "2024-01-20", reignEnd: "Present" },
    ],
  },

  welterweight: {
    displayName: "Welterweight (156‑170 lbs)",
    champions: [
      { name: "Pat Miletich",  reignStart: "1998-10-16", reignEnd: "2001-05-04" },
      { name: "Carlos Newton", reignStart: "2001-05-04", reignEnd: "2001-11-02" },
      { name: "Matt Hughes",   reignStart: "2001-11-02", reignEnd: "2004-01-31" },
      { name: "B.J. Penn",     reignStart: "2004-01-31", reignEnd: "2004-05-17", notes: "stripped (left UFC)" },
      { name: "Matt Hughes",   reignStart: "2004-10-22", reignEnd: "2006-11-18" },
      { name: "Georges St-Pierre",reignStart: "2006-11-18",reignEnd: "2007-04-07" },
      { name: "Matt Serra",    reignStart: "2007-04-07", reignEnd: "2008-04-19" },
      { name: "Georges St-Pierre",reignStart: "2008-04-19",reignEnd: "2013-12-13", notes: "vacated (hiatus)" },
      { name: "Johny Hendricks",reignStart: "2014-03-15", reignEnd: "2014-12-06" },
      { name: "Robbie Lawler", reignStart: "2014-12-06", reignEnd: "2016-07-30" },
      { name: "Tyron Woodley", reignStart: "2016-07-30", reignEnd: "2019-03-02" },
      { name: "Kamaru Usman",  reignStart: "2019-03-02", reignEnd: "2022-08-20" },
      { name: "Leon Edwards",  reignStart: "2022-08-20", reignEnd: "2024-07-27" },
      { name: "Belal Muhammad",reignStart: "2024-07-27", reignEnd: "Present" },
    ],
  },

  lightweight: {
    displayName: "Lightweight (146‑155 lbs)",
    champions: [
      { name: "Jens Pulver",        reignStart: "2001-02-23", reignEnd: "2002-03-23", notes: "vacated (contract)" },
      { name: "Sean Sherk",         reignStart: "2006-10-14", reignEnd: "2007-12-08", notes: "stripped" },
      { name: "B.J. Penn",          reignStart: "2008-01-19", reignEnd: "2010-04-10" },
      { name: "Frankie Edgar",      reignStart: "2010-04-10", reignEnd: "2012-02-26" },
      { name: "Benson Henderson",   reignStart: "2012-02-26", reignEnd: "2013-08-31" },
      { name: "Anthony Pettis",     reignStart: "2013-08-31", reignEnd: "2015-03-14" },
      { name: "Rafael dos Anjos",   reignStart: "2015-03-14", reignEnd: "2016-07-07" },
      { name: "Eddie Alvarez",      reignStart: "2016-07-07", reignEnd: "2016-11-12" },
      { name: "Conor McGregor",     reignStart: "2016-11-12", reignEnd: "2018-04-07", notes: "stripped (inactivity)" },
      { name: "Khabib Nurmagomedov",reignStart: "2018-04-07", reignEnd: "2021-03-19", notes: "retired" },
      { name: "Charles Oliveira",   reignStart: "2021-05-15", reignEnd: "2022-05-07", notes: "stripped (weight)" },
      { name: "Islam Makhachev",    reignStart: "2022-10-22", reignEnd: "Present" },
    ],
  },

  featherweight: {
    displayName: "Featherweight (136‑145 lbs)",
    champions: [
      { name: "José Aldo",          reignStart: "2010-11-20", reignEnd: "2015-12-12" },
      { name: "Conor McGregor",     reignStart: "2015-12-12", reignEnd: "2016-11-26", notes: "stripped" },
      { name: "José Aldo",          reignStart: "2016-11-26", reignEnd: "2017-06-03" },
      { name: "Max Holloway",       reignStart: "2017-06-03", reignEnd: "2019-12-14" },
      { name: "Alexander Volkanovski",reignStart: "2019-12-14",reignEnd: "2024-02-17" },
      { name: "Ilia Topuria",       reignStart: "2024-02-17", reignEnd: "2025-02-19", notes: "vacated (move to LW)" },
      { name: "Alexander Volkanovski",reignStart: "2025-04-12",reignEnd: "Present" },
    ],
  },

  bantamweight: {
    displayName: "Bantamweight (126‑135 lbs)",
    champions: [
      { name: "Dominick Cruz",     reignStart: "2010-12-16", reignEnd: "2014-01-06" },
      { name: "Renan Barão",       reignStart: "2014-01-06", reignEnd: "2014-05-24" },
      { name: "T.J. Dillashaw",    reignStart: "2014-05-24", reignEnd: "2016-01-17" },
      { name: "Dominick Cruz",     reignStart: "2016-01-17", reignEnd: "2016-12-30" },
      { name: "Cody Garbrandt",    reignStart: "2016-12-30", reignEnd: "2017-11-04" },
      { name: "T.J. Dillashaw",    reignStart: "2017-11-04", reignEnd: "2019-03-20", notes: "vacated (drug test)" },
      { name: "Henry Cejudo",      reignStart: "2019-06-08", reignEnd: "2020-05-24", notes: "vacated (retired)" },
      { name: "Petr Yan",          reignStart: "2020-07-12", reignEnd: "2021-03-06" },
      { name: "Aljamain Sterling", reignStart: "2021-03-06", reignEnd: "2023-08-19" },
      { name: "Sean O'Malley",     reignStart: "2023-08-19", reignEnd: "2024-09-14" },
      { name: "Merab Dvalishvili", reignStart: "2024-09-14", reignEnd: "Present" },
    ],
  },

  flyweight: {
    displayName: "Flyweight (116‑125 lbs)",
    champions: [
      { name: "Demetrious Johnson",reignStart: "2012-09-22", reignEnd: "2018-08-04" },
      { name: "Henry Cejudo",      reignStart: "2018-08-04", reignEnd: "2020-02-29", notes: "vacated" },
      { name: "Deiveson Figueiredo",reignStart: "2020-07-18",reignEnd: "2021-06-12" },
      { name: "Brandon Moreno",    reignStart: "2021-06-12", reignEnd: "2022-01-22" },
      { name: "Deiveson Figueiredo",reignStart: "2022-01-22",reignEnd: "2023-01-21" },
      { name: "Brandon Moreno",    reignStart: "2023-01-21", reignEnd: "2023-07-08" },
      { name: "Alexandre Pantoja", reignStart: "2023-07-08", reignEnd: "Present" },
    ],
  },

  /* -------------------------- WOMEN ------------------------------- */
  women_strawweight: {
    displayName: "Women – Strawweight (106‑115 lbs)",
    champions: [
      { name: "Carla Esparza",   reignStart: "2014-12-12", reignEnd: "2015-03-14" },
      { name: "Joanna Jędrzejczyk",reignStart: "2015-03-14",reignEnd: "2017-11-04" },
      { name: "Rose Namajunas",  reignStart: "2017-11-04", reignEnd: "2019-05-11" },
      { name: "Jéssica Andrade", reignStart: "2019-05-11", reignEnd: "2019-08-31" },
      { name: "Zhang Weili",     reignStart: "2019-08-31", reignEnd: "2021-04-24" },
      { name: "Rose Namajunas",  reignStart: "2021-04-24", reignEnd: "2022-05-07" },
      { name: "Carla Esparza",   reignStart: "2022-05-07", reignEnd: "2022-11-12" },
      { name: "Zhang Weili",     reignStart: "2022-11-12", reignEnd: "Present" },
    ],
  },

  women_flyweight: {
    displayName: "Women – Flyweight (116‑125 lbs)",
    champions: [
      { name: "Nicco Montaño",      reignStart: "2017-12-01", reignEnd: "2018-09-08", notes: "stripped" },
      { name: "Valentina Shevchenko",reignStart: "2018-12-08", reignEnd: "2023-03-04" },
      { name: "Alexa Grasso",       reignStart: "2023-03-04", reignEnd: "2024-09-14" },
      { name: "Valentina Shevchenko",reignStart: "2024-09-14", reignEnd: "Present" },
    ],
  },

  women_bantamweight: {
    displayName: "Women – Bantamweight (126‑135 lbs)",
    champions: [
      { name: "Ronda Rousey",     reignStart: "2012-12-06", reignEnd: "2015-11-15" },
      { name: "Holly Holm",       reignStart: "2015-11-15", reignEnd: "2016-03-05" },
      { name: "Miesha Tate",      reignStart: "2016-03-05", reignEnd: "2016-07-09" },
      { name: "Amanda Nunes",     reignStart: "2016-07-09", reignEnd: "2021-12-11" },
      { name: "Julianna Peña",    reignStart: "2021-12-11", reignEnd: "2022-07-30" },
      { name: "Amanda Nunes",     reignStart: "2022-07-30", reignEnd: "2023-06-20", notes: "retired" },
      /* vacant 2023‑06‑20 → 2024‑01‑20 */
      { name: "Raquel Pennington",reignStart: "2024-01-20", reignEnd: "2024-10-05" },
      { name: "Julianna Peña",    reignStart: "2024-10-05", reignEnd: "Present" },
    ],
  },
};

/* ------------------------------------------------------------------ */
/*  HELPERS                                                           */
/* ------------------------------------------------------------------ */

// --- Re-added/Ensured Helper Functions ---

const parseDate = (iso: string | null): Date | null => {
  if (!iso) return null;
  if (iso === "Present") return new Date();
  // allow YYYY‑MM or YYYY
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
    .replace(/[.'‑]/g, "") // Removed hyphen from removal
    .trim();

/*  NEW → collapse every weight‑class string to its plain division name  */
const normalizeWeightClass = (raw?: string | null): string => {
  if (!raw) return "";
  return normalizeName(raw)
    .toLowerCase()
    .replace(
      /(title|championship|bout|world|undisputed|interim)/g,
      ""
    )
    .replace(/\s+/g, " ")
    .trim(); //  e.g.  "Lightweight title bout" → "lightweight"
};

const parseFighterName = (full: string) => {
  const parts = full.trim().split(/\s+/);
  return { firstName: parts.slice(0, -1).join(" "), lastName: parts.at(-1) || "" };
};

// Improved extractOpponent to fix "Alexander Volkano ki" issue
const extractOpponent = (boutString: string, normalizedFighterName: string): string => {
  if (!boutString || !normalizedFighterName) return 'Opponent Unknown';

  // Split on "vs" variations and trim spaces
  const parts = boutString.split(/\s+vs\.?\s+|\s+vs\s+/i);

  if (parts.length === 2) {
    // Get normalized versions of both names for comparison
    const normalized1 = normalizeName(parts[0]).toLowerCase();
    const normalized2 = normalizeName(parts[1]).toLowerCase();
    const normalized = normalizedFighterName.toLowerCase();

    // Return the non-matching part (the opponent)
    if (normalized1.includes(normalized) || normalized.includes(normalized1)) {
      return parts[1].trim();
    } else if (normalized2.includes(normalized) || normalized.includes(normalized2)) {
      return parts[0].trim();
    }
  }

  return 'Opponent Unknown';
};

// Updated extractDecisionDetails to remove parentheses around scores
const extractDecisionDetails = (details: string | null | undefined): string | undefined => {
  if (!details) return undefined;
  // Look for patterns like "XX - XX" scores
  const scorePattern = /\d{1,2}\s*-\s*\d{1,2}/g;
  const scores = details.match(scorePattern);
  // Return only scores if 3 sets are found, without parentheses
  if (scores && scores.length >= 3) {
    return scores.join(', '); // No parentheses
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

// Helper function to calculate age from DOB string (YYYY-MM-DD)
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

/* NEW Helper: Get Initials */
const getInitials = (name: string): string => {
  if (!name) return "?";
  
  // Special cases for specific fighters (case-insensitive, space/hyphen flexible)
  const normalizedName = name.toLowerCase().replace(/[‑\s]+/g, ' ').trim(); // Normalize spaces/hyphens
  if (normalizedName === "b j penn") return "BJ";
  if (normalizedName === "georges st pierre") return "GSP";
  if (normalizedName === "rafael dos anjos") return "RDA";
  if (normalizedName === "dricus du plessis") return "DDP";

  // Standard case - first letter of first part and first letter of last part
  const parts = name
    .replace(/‑/g, "-") // Normalize hyphens first
    .split(/[\s-]+/) // Split by space or hyphen
    .filter(Boolean); // Remove empty strings

  if (parts.length === 0) return "?";
  if (parts.length === 1) return parts[0].substring(0, 2).toUpperCase();

  // For three part names like "Rafael Dos Anjos", try to be smarter about it
  if (parts.length === 3) {
    const middle = parts[1].toLowerCase();
    // If middle part is a connector like "dos", "de", "du", etc., use first and last parts
    if (["dos", "de", "du", "da", "von", "van", "del"].includes(middle)) {
      return (parts[0][0] + parts[2][0]).toUpperCase();
    }
  }

  // Regular case: first letter of first name + first letter of last name
  return (parts[0][0] + (parts.length > 1 ? parts[parts.length - 1][0] : "")).toUpperCase();
};

/* NEW Helper: Abbreviate Weight Class */
const abbreviateWeightClass = (key: string): string => {
    const abbreviations: Record<string, string> = {
        heavyweight: "HW",
        lightheavyweight: "LHW",
        middleweight: "MW",
        welterweight: "WW",
        lightweight: "LW",
        featherweight: "FeW",
        bantamweight: "BW",
        flyweight: "FlW",
        women_strawweight: "WSW",
        women_flyweight: "WFlW",
        women_bantamweight: "WBW",
    };
    return abbreviations[key.toLowerCase()] || key.toUpperCase();
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
  // Add new state for gender filter
  const [genderFilter, setGenderFilter] = useState<"all" | "men" | "women">("all"); // Default to 'all' initially
  const dataCache = useRef<
    Record<
      string,
      { details: FighterDetails | null; images: UfcChampionImages | null; history: ProcessedFightResult[] }
    >
  >({});
  const scrollContainerRef = useRef<HTMLDivElement>(null);
  const timelineScrollContainerRef = useRef<HTMLDivElement>(null);
  // --- END REINSTATE STATE DECLARATIONS ---

  const classes = ufcChampionsData;
  const data = classes[selected];

  // NEW: Filter classes based on gender
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


  // Effect for card view scroll
  useEffect(() => {
    if (scrollContainerRef.current && data && !showDetailsView && !isTimelineView) {
      scrollContainerRef.current.scrollLeft = scrollContainerRef.current.scrollWidth;
    }
  }, [data, showDetailsView, isTimelineView]);

  // NEW Effect for timeline view scroll
  useEffect(() => {
    if (isTimelineView && timelineScrollContainerRef.current) {
      timelineScrollContainerRef.current.scrollLeft = timelineScrollContainerRef.current.scrollWidth;
    }
  }, [isTimelineView]);

  /* ================================================================== */
  /*  FETCH  +  **UPDATED DEFENCE LOGIC**                               */
  /* ================================================================== */
  const fetchFighterData = async (fighter: Champion) => {
    if (!fighter) return;
    const original = fighter.name.replace(/ /g, " ").trim();
    const normName = normalizeName(original);
    console.log(`Fetching data for: ${original} (Normalized: ${normName})`);

    // --- Check Cache ---
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
      // Log the names being used for the image query
      console.log(`Querying images with FIRST: %${firstName}%, LAST: %${lastName}%`);

      /* --- details + images fetch (Updated details query) ----------- */
      const [detailsRes, imagesRes] = await Promise.all([
        supabase
          .from("ufc_fighter_tott") // Use the new table name
          .select("*")
          .ilike("FIGHTER", `%${normName}%`) // Match using the FIGHTER column and normalized name
          .limit(1)
          .maybeSingle(),
        supabase
          .from("ufc_champion_images") // Keep image query as is, assuming it uses FIRST/LAST
          .select("*")
          .ilike("FIRST", `%${firstName}%`) // Use parsed name for images
          .ilike("LAST", `%${lastName}%`)   // Use parsed name for images
          .limit(1)
          .maybeSingle(),
      ]);

      // ... rest of the function remains the same ...
      if (detailsRes.error) console.error("Supabase details fetch error:", detailsRes.error);
      const fighterDetailsFetched = detailsRes.data as FighterDetails | null;
      console.log("Fetched Details:", fighterDetailsFetched);

      if (imagesRes.error) console.error("Supabase images fetch error:", imagesRes.error);
      const fighterImagesFetched = imagesRes.data as UfcChampionImages | null;
      console.log("Fetched Images:", fighterImagesFetched); // This log confirms if it's null


      /* --- fight‑history fetch  ------------------------------------ */
      const fightsRes = await supabase
        .from("ufc_fight_results")
        .select("*")
        .ilike("BOUT", `%${normName}%`)
        .order("id", { ascending: true }); // chronological

      if (fightsRes.error) throw fightsRes.error;
      const fights = (fightsRes.data || []) as FightResult[]; // Ensure fights is an array
      console.log("Fetched Raw Fight Results (Sorted):", fights);

      /* -------------------------------------------------------------- */
      /*  REAL TIME‑LINE WALK WITH “championStatus” MAP                 */
      /* -------------------------------------------------------------- */
      const processed: ProcessedFightResult[] = [];
      const championStatus: Record<string, boolean> = {}; // keyed by clean division

      // Reverse loop for newest-first display order
      for (let idx = fights.length - 1; idx >= 0; idx--) {
        const fight = fights[idx];
        const resultForFighter = getResultForFighter(fight, normName);
        const div = normalizeWeightClass(fight.WEIGHTCLASS);
        const isTitle = fight.WEIGHTCLASS?.toLowerCase().includes("title") ?? false;
        let defenseOutcome: "success" | "fail" | null = null;

        if (isTitle && div) {
          if (championStatus[div]) {
            // already champ → this is a defence attempt
            if (resultForFighter === "loss") {
              defenseOutcome = "fail";
              championStatus[div] = false; // belt lost
            } else if (resultForFighter === "win") { // *** CHANGED: Only count 'win' as success ***
              defenseOutcome = "success"; // win retains belt
              // retain champ flag
            }
            // NC or Draw means belt is retained, but not marked as 'success' or 'fail' defense
          } else {
            // challenger
            if (resultForFighter === "win") championStatus[div] = true; // new champ
            // otherwise championStatus[div] stays false
          }
        }

        /* --- optional decision / finish details (same as before) ---- */
        let decisionDetails: string | undefined;
        let finishDetails: string | undefined;
        if (fight.METHOD?.toLowerCase().includes("decision")) {
          decisionDetails = extractDecisionDetails(fight.DETAILS);
        } else if (fight.DETAILS?.trim()) {
          if (!/\d{1,2}\s*-\s*\d{1,2}/.test(fight.DETAILS))
            finishDetails = fight.DETAILS.trim();
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
      console.log("Processed History:", processed);

      /* --- cache + set state  -------------------------------------- */
      dataCache.current[normName] = { // Use dataCache here
        details: fighterDetailsFetched,
        images: fighterImagesFetched,
        history: processed,
      };
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

  // Generate a unique layoutId for a champion
  const getLayoutId = (champion: Champion | null) => {
      if (!champion) return undefined;
      // Use a combination of name and start date for uniqueness
      return `champion-card-${normalizeName(champion.name)}-${champion.reignStart}`;
  }

  /* ================================================================== */
  /*  TIMELINE VIEW CALCULATIONS                                        */
  /* ================================================================== */
  const timelineData = useMemo(() => {
    // Use filteredClasses here
    if (!isTimelineView || Object.keys(filteredClasses).length === 0) return null;
    setIsTimelineCalculating(true);
    console.log("Calculating timeline data...");
    let minDate = new Date();
    let maxDate = new Date(1990, 0, 1);
    const today = new Date();
    
    // Find min/max dates using filteredClasses
    Object.values(filteredClasses).forEach((wc: WeightClassData) => { // Explicit type for wc
      wc.champions.forEach((c: Champion) => { // Explicit type for c
        const start = parseDate(c.reignStart);
        const end = parseDate(c.reignEnd);
        if (start && start < minDate) minDate = start;
        if (end) {
          // If end date is in the future, cap it at today's date
          const cappedEnd = end > today ? today : end;
          if (cappedEnd > maxDate) maxDate = cappedEnd;
        } else { // Handle "Present" case for maxDate
            if (start && start > maxDate) maxDate = today; // If current champ started after current max, set max to today
            else if (today > maxDate) maxDate = today; // Ensure maxDate includes today if there's an active champ
        }
      });
    });

    // Don't extend into the future, use today as max date if needed
    if (maxDate > today) maxDate = today;
    
    // Add small buffer for display (don't set to first of month)
    maxDate.setDate(maxDate.getDate() + 10);
    minDate.setDate(minDate.getDate() - 10);

    const totalDays = diffInDays(minDate, maxDate);
    const pixelsPerDay = 0.75;
    const todayX = diffInDays(minDate, today) * pixelsPerDay;

    // Generate month markers
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
  }, [isTimelineView, filteredClasses]); // Add filteredClasses dependency


  /* ------------------------------------------------------------------ */
  /*  UI                                                                */
  /* ------------------------------------------------------------------ */
  return (
    <main className="relative pt-8 md:pt-6 px-4 md:px-6 pb-20 w-full flex flex-col h-screen"> {/* Added flex flex-col h-screen */}

      {/* --- View Toggle and Gender Filter --- */}
      {/* Added min-h-[...] to stabilize height */}
      <div className="w-full max-w-screen-2xl mx-auto flex justify-between items-start mb-4 px-0 md:px-0 min-h-[40px]">
        <div className="flex-1">
          {/* Conditionally render gender filter only when timeline is active */}
          {isTimelineView && (
            <div className="flex items-center space-x-2">
              <Select
                value={genderFilter}
                onValueChange={(value: "all" | "men" | "women") => setGenderFilter(value)}
                disabled={showDetailsView}
              >
                <SelectTrigger className="w-180px]">
                  <SelectValue placeholder="Filter Divisions" />
                </SelectTrigger>
                <SelectContent>
                  {/* Updated Order: Men, Women, All */}
                  <SelectItem value="men">Men's Divisions</SelectItem>
                  <SelectItem value="women">Women's Divisions</SelectItem>
                  <SelectItem value="all">All Divisions</SelectItem>
                </SelectContent>
              </Select>
            </div>
          )}
        </div>
        <div className="flex items-center space-x-2 pt-2"> {/* Added pt-2 for alignment */}
          <Checkbox
            id="timeline-view-toggle"
            checked={isTimelineView}
            onCheckedChange={(checked) => {
              const isChecked = Boolean(checked);
              setIsTimelineView(isChecked);
              // Set default filter to 'men' when enabling timeline, 'all' when disabling
              setGenderFilter(isChecked ? "men" : "all");
            }}
            disabled={showDetailsView}
          />
          <Label htmlFor="timeline-view-toggle" className="text-sm font-medium">
            Show Full Timeline View
          </Label>
        </div>
      </div>


      {/* --- Conditional Rendering: Card View or Timeline View --- */}
      {!isTimelineView ? (
        /* --- Main Content Card (Background) --- */
        /* Removed flex-1 h-full, added mb-auto to push footer down */
        <Card className={cn(
          "w-full max-w-screen-2xl mx-auto transition-opacity duration-300 flex flex-col mb-auto",
          showDetailsView ? "opacity-20 pointer-events-none" : "opacity-100"
        )}>
           {/* Removed flex-1 from wrapper div */}
           <div className="overflow-y-auto">
             <CardContent className="space-y-6 pt-6">
                {/* selector */}
                <Select value={selected} onValueChange={setSelected} disabled={showDetailsView}>
                   <SelectTrigger className="w-full md:w-[300px] mx-auto">
                    <SelectValue placeholder="Select a Weight Class" />
                  </SelectTrigger>
                  <SelectContent>
                    {Object.entries(classes).map(([key, val]) => (
                      <SelectItem key={key} value={key}>
                        {val.displayName}
                      </SelectItem>
                    ))}
                  </SelectContent>
              </Select>

                 {/* champions timeline (Card View) */}
                 {!data ? ( // Show skeleton if data for selected class isn't ready (unlikely but good practice)
                    <div className="mt-6">
                      <div className="w-full overflow-x-auto pb-4">
                        <div className="flex space-x-3 items-stretch min-h-[200px]">
                          {[...Array(5)].map((_, i) => (
                            <Skeleton key={i} className="flex-shrink-0 w-52 h-[180px] rounded-md" />
                          ))}
                        </div>
                      </div>
                    </div>
                 ) : (
                    <div className="mt-6">
                      <div ref={scrollContainerRef} className="w-full overflow-x-auto pb-4">
                        <div className="flex space-x-3 items-stretch min-h-[200px]">
                          {data.champions.map((c, i) => {
                            // ... champion card rendering logic ...
                            const start = parseDate(c.reignStart)!;
                            const end = parseDate(c.reignEnd)!;
                            const reignDays = diffInDays(start, end);
                            const duration = humanDuration(reignDays);
                            const next = data.champions[i + 1];
                            let vacancyBlock: React.ReactNode = null;
                             if (next) {
                              const gap =
                                diffInDays(
                                  end,
                                  parseDate(next.reignStart) || new Date()
                                ) - 1;
                              if (gap > 0) {
                                vacancyBlock = (
                                  <div className="flex flex-col items-center justify-center px-2 text-center">
                                     <div className="w-px h-16 bg-destructive/50 my-1"></div>
                                     <span className="text-xs text-destructive whitespace-normal w-16">
                                        VACANT<br/>{humanDuration(gap)}
                                     </span>
                                     <div className="w-px h-16 bg-destructive/50 my-1"></div>
                                  </div>
                                );
                              }
                            } else if (
                              c.reignEnd &&
                              c.reignEnd !== "Present" &&
                              !c.notes?.toLowerCase().includes("interim")
                            ) {
                               vacancyBlock = (
                                  <div className="flex flex-col items-center justify-center px-2 text-center">
                                     <div className="w-px h-16 bg-destructive/50 my-1"></div>
                                     <span className="text-xs text-destructive whitespace-normal w-16">
                                        VACANT
                                     </span>
                                     <div className="w-px h-16 bg-destructive/50 my-1"></div>
                                  </div>
                                );
                            }


                            return (
                              <React.Fragment key={`${c.name}-${i}`}>
                                <motion.div
                                  layoutId={getLayoutId(c)} // Keep layoutId for potential animation back
                                  onClick={() => fetchFighterData(c)}
                                  className="
                                    flex-shrink-0 w-52 p-3 bg-muted rounded-md shadow-sm
                                    cursor-pointer relative z-0 hover:z-10 flex flex-col justify-between group overflow-hidden"
                                  whileHover={{ scale: 1.05, zIndex: 10 }} // Example hover effect
                                  transition={{ duration: 0.2 }} // Adjust transition
                                >
                                  {/* ... Card content (Image Placeholder, Name, Dates, Duration, Notes) ... */}
                                  <div className="relative h-24 w-full mb-2 bg-gradient-to-b from-gray-300 to-gray-400 dark:from-gray-600 dark:to-gray-700 rounded flex items-center justify-center">
                                      <User className="h-12 w-12 text-gray-500 dark:text-gray-400" />
                                  </div>
                                  <div>
                                      <p className="font-bold text-base leading-tight mb-1 truncate">{c.name}</p>
                                      <p className="text-xs">
                                      {start.toLocaleDateString(undefined, { month: '2-digit', day: '2-digit', year: '2-digit' })}
                                      {' – '}
                                      {c.reignEnd === "Present"
                                          ? "Present"
                                          : end.toLocaleDateString(undefined, { month: '2-digit', day: '2-digit', year: '2-digit' })}
                                      </p>
                                      <p className="text-xs text-muted-foreground">
                                      {duration}
                                      </p>
                                  </div>
                                  <div className="absolute bottom-0 left-0 right-0 p-1">
                                      {c.notes && (
                                      <p className="text-xs mb-1 italic text-primary-foreground bg-primary/80 p-1 rounded opacity-0 max-h-0 overflow-hidden transition-all duration-300 ease-in-out group-hover:opacity-100 group-hover:max-h-40">
                                          {c.notes}
                                      </p>
                                      )}
                                  </div>
                                </motion.div>
                                {/* Vacancy Indicator */}
                                {vacancyBlock}
                              </React.Fragment>
                            );
                          })}
                        </div>
                      </div>
                    </div>
                  )}

              <p className="text-xs text-center text-muted-foreground mt-8">
                Missing UFC 294 and 308 Events.
              </p>
             </CardContent>
           </div> {/* End scrollable wrapper */}
         </Card>
      ) : (
        /* --- Timeline View --- */
        (isTimelineCalculating || !timelineData) ? ( // Show skeleton if calculating or data not ready
            <div className="w-full max-w-screen-2xl mx-auto border rounded-md bg-card p-4">
                <div className="flex space-x-2">
                    <Skeleton className="h-[75vh] w-[2.5rem]" /> {/* Sidebar skeleton */}
                    <div className="flex-1 space-y-2">
                        {/* Use filteredClasses for skeleton count */}
                        {[...Array(Object.keys(filteredClasses).length)].map((_, i) => (
                            <Skeleton key={i} className="h-[3.5rem] w-full" />
                        ))}
                        <Skeleton className="h-[1.5rem] w-full mt-2" /> {/* Marker row skeleton */}
                    </div>
                </div>
            </div>
        ) : ( // Render actual timeline
            <div
              ref={timelineScrollContainerRef}
              className="w-full max-w-screen-2xl mx-auto overflow-x-auto border rounded-md bg-card text-card-foreground flex-1 min-h-[70vh]" // Added flex-1
            >
              {/* Container for scrolling content */}
              <div
                className="relative px-0 py-0 overflow-visible h-full" // Added h-full
                style={{
                    // Removed minHeight
                    width: `${timelineData.totalWidth + (2.5 * 16)}px`, // 2.5rem sidebar width
                    paddingBottom: '1rem', // Keep padding to push content up from scrollbar
                }}
              >
                {/* Grid for rows, sidebar, and markers */}
                <div
                  className="relative z-10 overflow-visible grid h-full" // Ensure h-full here too
                  style={{
                    // Use filteredClasses for grid rows
                    // Use 1fr for champion rows to distribute height equally
                    gridTemplateRows: `1rem repeat(${Object.keys(filteredClasses).length}, 1fr) 1.5rem`, // 1fr rows will now expand
                    gridTemplateColumns: '2.5rem 1fr',
                    width: '100%'
                  }}
                >
                  {/* Empty cell for top-left corner */}
                  <div
                    className="sticky left-0 bg-card z-30 border-r border-border/60"
                    style={{ gridRow: 1, gridColumn: 1 }}
                  ></div>

                  {/* Top Markers Row */}
                  <div
                    className="relative h-full border-b border-border/60" // Border below top markers
                    style={{ gridRow: 1, gridColumn: 2 }}
                  >
                    {/* Top Month/Year Ticks */}
                    {timelineData?.monthMarkers.map(({ date, x }) => {
                       const isYear = date.getMonth() === 0;
                       return (
                         <div
                           key={`top-month-${date.toISOString()}`}
                           className="absolute bottom-0 flex flex-col items-center" // Align to bottom of this row
                           style={{ left: `${x}px`, height: isYear ? '0.75rem' : '0.375rem' }} // Control height via style
                         >
                           <div className="w-px h-full bg-foreground/50"></div>
                         </div>
                       );
                    })}
                    {/* Top Present Day Marker */}
                    {timelineData && timelineData.todayX >= 0 && (
                        <div
                            key="top-today-marker"
                            className="absolute bottom-0 z-10" // Align to bottom
                            style={{ left: `${timelineData.todayX}px`, height: '1rem' }} // Slightly taller
                        >
                            <div className="w-px h-full bg-primary/70"></div>
                        </div>
                    )}
                  </div>


                  {/* Weight Class Abbreviations and Champion Bars */}
                  {/* Use filteredClasses here */}
                  {Object.entries(filteredClasses).map(([key,wc],r)=>(
                    <div key={key} className="contents">
                      {/* abbrev name */}
                      <div
                        className="sticky left-0 bg-card z-30 flex items-center justify-center px-1 border-r border-b border-border/60"
                        // Row index needs to account for the top marker row (r + 2)
                        style={{ gridRow: r + 2, gridColumn: 1, width: '2.5rem' }}
                      >
                        <span className="text-[12px] font-bold uppercase text-center">
                          {abbreviateWeightClass(key)}
                        </span>
                      </div>
                      {/* champion bars */}
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
                          const reignYears = diffInDays(s, e) / 365; // Calculate reign length in years
                          const isLongReign = reignYears >= 1; // Check if reign is 1 year or longer
                          let h=0; for(let k=0;k<c.name.length;k++) h=c.name.charCodeAt(k)+((h<<5)-h);
                          const hue=h%360, bg=`hsl(${hue},70%,90%)`, hov=`hsl(${hue},70%,80%)`, col=`hsl(${hue},80%,25%)`;
                          
                          return (
                            <motion.div
                              key={`${key}-${i}`}
                              // Adjusted padding and margin for better fit within the row height
                              className="absolute top-[2px] bottom-[2px] rounded-sm flex items-center justify-center group cursor-pointer shadow-sm z-10" 
                              style={{ 
                                left:`${left}px`, 
                                width:`${Math.max(width,2)}px`, 
                                backgroundColor:bg, 
                                color:col,
                                borderRadius: '0.125rem' // Explicitly set borderRadius consistently
                              }}
                              onClick={() => fetchFighterData(c)}
                              initial={{ scale:1 }}
                              whileHover={{
                                backgroundColor:hov,
                                zIndex:20,
                                borderRadius: '0.125rem' // Keep border radius the same on hover
                              }}
                              transition={{ duration:0.15 }}
                            >
                              {/* Show full name for long reigns, otherwise show initials */}
                              {isLongReign ? (
                                <span className="text-[20px] font-bold whitespace-nowrap px-1 overflow-hidden text-ellipsis max-w-full">
                                  {c.name}
                                </span>
                              ) : (
                                <span
                                  className="text-[20px] font-bold whitespace-nowrap overflow-visible group-hover:opacity-0 transition-opacity duration-100"
                                  style={{ position: 'relative', zIndex: 15 }}
                                >
                                  {getInitials(c.name)}
                                </span>
                              )}
                              
                              {/* Full name on hover - only needed for short reigns */}
                              {!isLongReign && (
                                <span
                                  className="absolute inset-0 flex items-center justify-center text-[15px] font-semibold px-1 opacity-0 group-hover:opacity-100 transition-opacity duration-100 delay-50 pointer-events-none"
                                  style={{ background: hov, color: col, borderRadius: '0.125rem', zIndex: 21 }}
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

                  {/* Empty cell for bottom-left corner */}
                  <div
                    className="sticky left-0 bg-card z-30 border-r border-border/60"
                    // Row index is number of filtered classes + 2 (for top marker row)
                    style={{ gridRow: Object.keys(filteredClasses).length + 2, gridColumn: 1 }}
                  ></div>

                  {/* Bottom Markers Row - Adjusted height via gridTemplateRows */}
                  <div
                    className="relative h-full" // No border needed here if rows above have border-b
                    // Row index is number of filtered classes + 2
                    style={{ gridRow: Object.keys(filteredClasses).length + 2, gridColumn: 2 }}
                  >
                    {/* Bottom Month/Year Ticks */}
                    {timelineData?.monthMarkers.map(({ date, x }) => {
                       const isYear = date.getMonth() === 0;
                       return (
                         <div
                           key={`bottom-month-${date.toISOString()}`}
                           className="absolute top-0 flex flex-col items-center" // Align to top of this row
                           style={{ left: `${x}px`, height: isYear ? '0.75rem' : '0.375rem' }} // Control height
                         >
                           <div className="w-px h-full bg-foreground/50"></div>
                           {isYear && (
                             <span className="absolute top-full text-[12px] text-foreground transform -translate-x-1/2 mt-0.5"> {/* Centered label */}
                               {date.getFullYear()}
                             </span>
                           )}
                         </div>
                       );
                    })}
                    {/* Present Day Marker */}
                    {timelineData && timelineData.todayX >= 0 && (
                        <div
                            key="bottom-today-marker"
                            className="absolute top-0 z-10 flex flex-col items-center"
                            style={{ left: `${timelineData.todayX}px`, height: '1rem' }}
                        >
                            <div className="w-px h-full bg-primary/70"></div>
                             <span className="absolute top-full text-[9px] text-primary font-semibold transform -translate-x-1/2 mt-0.5 bg-card px-0.5"> {/* Centered label */}
                                Now
                            </span>
                        </div>
                    )}
                  </div>
                </div> {/* End of Grid */}
              </div> {/* End of scrollable content container */}
            </div> /* End of Timeline View */
        )
      )}


      {/* --- Fighter Details Overlay --- */}
      <AnimatePresence>
        {showDetailsView && selectedFighter && (
          <motion.div key="overlay-backdrop"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            transition={{ duration: 0.3 }}
            className="fixed inset-0 z-40 flex items-center justify-center bg-black/60 backdrop-blur-sm p-4 md:p-8"
            onClick={handleCloseDetails} // Close on backdrop click
          >
            {/*
              NOTE: The layoutId animation might look strange when transitioning
              from the Timeline View's small rectangle. Consider disabling layoutId
              or adjusting the animation if isTimelineView was true when opening.
              For simplicity, we keep it for now.
            */}
            <motion.div
              layoutId={getLayoutId(selectedFighter)}
              className="bg-card text-card-foreground rounded-md shadow-xl w-full max-w-6xl h-[90vh] max-h-[800px] flex flex-col relative overflow-hidden z-50"
              onClick={(e) => e.stopPropagation()} // Prevent closing when clicking inside
            >
              {/* Close Button */}
              <motion.button
                initial={{ opacity: 0, scale: 0.5 }}
                animate={{ opacity: 1, scale: 1, transition: { delay: 0.2 } }}
                exit={{ opacity: 0 }}
                onClick={handleCloseDetails}
                className="absolute top-3 right-3 text-muted-foreground hover:text-foreground z-50 p-1 rounded-full bg-card hover:bg-muted"
                aria-label="Close details"
              >
                <X size={20} />
              </motion.button>

              {/* Main content area for overlay */}
              <div className="flex flex-1 h-full overflow-hidden pt-10"> {/* Added pt-10 here */}
                  {isLoading && (
                    <div className="absolute inset-0 flex items-center justify-center bg-card/80 z-20">
                       <Loader2 className="h-8 w-8 animate-spin text-primary" />
                    </div>
                  )}
                  {error && !isLoading && (
                     <div className="absolute inset-0 flex flex-col items-center justify-center bg-card/80 z-20 p-4 text-center">
                       <p className="text-destructive font-semibold mb-2">Error</p>
                       <p className="text-sm text-destructive">{error}</p>
                       <button onClick={handleCloseDetails} className="mt-4 text-sm underline">Close</button>
                     </div>
                  )}
                  {!isLoading && !error && (
                    <>
                      {/* Left Panel: Stats */}
                      <div className="w-1/3 border-r border-border p-4 md:p-6 overflow-y-auto flex flex-col items-center">
                         {/* Image - Uses fighterImages.HEADSHOT */}
                         <div className="relative h-40 w-40 mb-4 bg-gradient-to-b from-gray-300 to-gray-400 dark:from-gray-600 dark:to-gray-700 rounded-full flex items-center justify-center overflow-hidden">
                            {fighterImages?.HEADSHOT ? (
                                <Image
                                    src={fighterImages.HEADSHOT}
                                    alt={selectedFighter.name}
                                    layout="fill"
                                    objectFit="cover"
                                    // Optionally add unoptimized={true} if Next.js optimization is suspected issue
                                    // unoptimized={true}
                                />
                            ) : (
                                <User className="h-24 w-24 text-gray-500 dark:text-gray-400" />
                            )}
                         </div>

                         {/* Name */}
                         <h2 className="text-2xl font-bold mb-4 text-center">{fighterDetails?.FIGHTER || selectedFighter.name}</h2>

                         {/* Details */}
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
                           </div>
                         ) : (
                           <p className="text-sm text-muted-foreground mt-4">No detailed stats available.</p>
                         )}
                      </div> {/* Closing tag for Left Panel */}


                      {/* Right Panel: Fight History */}
                      <div className="w-2/3 p-4 md:p-6 overflow-y-auto">
                        {fighterHistory.length > 0 ? (
                          <ul className="space-y-1"> {/* Keep space-y-1 for internal list item spacing */}
                            {fighterHistory.map((fight) => {
                              const isTitleFight = fight.WEIGHTCLASS?.toLowerCase().includes('title') ?? false;
                              return (
                                <li
                                  key={fight.id}
                                  className={cn(
                                    "relative text-sm border-b border-border p-2 rounded-lg transition-colors duration-150 overflow-visible",
                                    // Add margin-bottom for spacing between list items
                                    "mb-2", // Added margin-bottom here
                                    // ... background colors ...
                                    (fight.resultForFighter === 'draw' || fight.resultForFighter === 'nc') && 'bg-gray-100 dark:bg-gray-800/30 border-gray-300 dark:border-gray-600',
                                    fight.resultForFighter === 'win' && 'bg-green-100 dark:bg-green-900/30 border-green-300 dark:border-green-700',
                                    fight.resultForFighter === 'loss' && 'bg-red-100 dark:bg-red-900/30 border-red-300 dark:border-red-700',
                                    fight.resultForFighter === 'unknown' && 'bg-muted/50',
                                    "hover:shadow-md cursor-pointer"
                                  )}
                                  onClick={() => {
                                    if (fight.URL) window.open(fight.URL, '_blank', 'noopener,noreferrer');
                                  }}
                                >
                                  {/* Icons */}
                                  {isTitleFight && !fight.defenseOutcome && (
                                    <Crown
                                      size={16}
                                      className="absolute -top-1 -left-2 text-yellow-500 dark:text-yellow-400 opacity-70 -rotate-45 z-10"
                                      aria-label="Title Fight"
                                    />
                                  )}
                                  {fight.defenseOutcome === 'success' && (
                                    <ShieldCheck
                                      size={18}
                                      className="absolute -top-1 -left-2 text-blue-500 dark:text-blue-400 opacity-100 -rotate-45 z-10"
                                      aria-label="Successful Title Defence"
                                    />
                                  )}
                                  {fight.defenseOutcome === 'fail' && (
                                    <ShieldOff
                                      size={18}
                                      className="absolute -top-1 -left-2 text-red-500 dark:text-red-400 opacity-100 -rotate-45 z-10"
                                      aria-label="Lost Title Defence"
                                    />
                                  )}

                                  {/* Result Badge and Opponent */}
                                  <div className="flex justify-between items-center">
                                    <div>
                                      <span className={cn(
                                          "font-semibold mr-2 px-1.5 py-0.5 rounded text-xs uppercase",
                                          fight.resultForFighter === 'win' && 'bg-green-600 text-white',
                                          fight.resultForFighter === 'loss' && 'bg-red-600 text-white',
                                          (fight.resultForFighter === 'draw' || fight.resultForFighter === 'nc') && 'bg-gray-500 text-white',
                                          fight.resultForFighter === 'unknown' && 'bg-gray-400 text-white'
                                      )}>
                                          {fight.resultForFighter === 'draw' ? 'DRAW' : fight.resultForFighter === 'nc' ? 'NC' : fight.resultForFighter}
                                      </span>
                                      <span>vs. </span>
                                      <a
                                        href={createAthleteUrl(fight.OPPONENT || '')}
                                        target="_blank"
                                        rel="noopener noreferrer"
                                        className="font-bold hover:underline text-primary"
                                        onClick={(e) => e.stopPropagation()}
                                      >
                                        {fight.OPPONENT || 'Opponent Unknown'}
                                      </a>
                                    </div>
                                    <span className="text-xs text-muted-foreground">{fight.EVENT}</span>
                                  </div>

                                  {/* Method/Round/Time/Details */}
                                   <div>
                                    <p className="text-xs text-muted-foreground mt-0.5">
                                       {fight.METHOD}
                                       {fight.ROUND && ` | Rd: ${fight.ROUND}`}
                                       {fight.TIME && ` | Time: ${fight.TIME}`}
                                    </p>
                                    {fight.decisionDetails && (
                                      <p className="text-xs text-muted-foreground italic mt-0.5">
                                        {fight.decisionDetails}
                                      </p>
                                    )}
                                    {fight.DETAILS && (
                                      <p className="text-xs text-muted-foreground mt-0.5">
                                        {fight.DETAILS}
                                      </p>
                                    )}
                                  </div>
                                </li>
                              );
                            })}
                          </ul> // Closing tag for unordered list
                        ) : (
                          <p className="text-sm text-muted-foreground">No fight history found or data unavailable.</p>
                        )}
                      </div> {/* Closing tag for Right Panel */}
                    </>
                  )}
              </div> {/* Closing tag for main content area */}
            </motion.div> {/* Closing tag for motion.div overlay content */}
          </motion.div> // Closing tag for motion.div backdrop
        )}
      </AnimatePresence>
    </main> // Closing tag for main element
  );
};

export default UfcChampionsDisplay;




