import { createClient } from "@supabase/supabase-js";
import { NextResponse } from "next/server";

function createServerSupabaseClient() {
  const supabaseUrl = process.env.SUPABASE_URL;
  const supabaseAnonKey = process.env.SUPABASE_ANON_KEY;
  if (!supabaseUrl || !supabaseAnonKey)
    throw new Error("Supabase environment variables are missing.");
  return createClient(supabaseUrl, supabaseAnonKey);
}

/* ------------------------------------------------------------------ */
/*  ELO HELPERS                                                       */
/* ------------------------------------------------------------------ */
const DEFAULT_ELO = 1500;
const FINISH_BONUS = 5; // extra Elo for KO/TKO/Submission wins

// Smooth K decay: 40 at debut → 24 after 30 fights (continuous linear taper)
function getK(fightCount: number): number {
  return Math.max(24, 40 - (16 * Math.min(fightCount, 30)) / 30);
}

/**
 * Early finishes are more dominant than late finishes.
 * Infers 5-round fight if the ending round > 3.
 * R1 of 5-round fight → ×1.4 bonus; final round → ×1.0.
 */
function roundFinishMultiplier(round: number | null): number {
  if (!round || round < 1) return 1.0;
  const totalRounds = round > 3 ? 5 : 3;
  return 1 + (totalRounds - round) * 0.1;
}

/**
 * Split/majority decisions are more uncertain than unanimous ones.
 * Returns a discount multiplier applied to the full Elo change.
 */
function splitDecisionMultiplier(method: string): number {
  const lm = method.toLowerCase();
  if (lm.includes("split")) return 0.82;
  if (lm.includes("majority")) return 0.91;
  return 1.0;
}

/**
 * For decision wins, parse scorecard details (e.g. "29-28 48-47 30-27") and
 * return a multiplier that amplifies Elo changes for lopsided wins/losses.
 * A close 29-28 (avg margin 1) stays at ×1.0; a dominant 50-45 reaches ×1.4.
 */
function decisionMarginMultiplier(details: string | null | undefined): number {
  if (!details) return 1.0;
  const scoreRegex = /\b(\d{2})\s*-\s*(\d{2})\b/g;
  let totalMargin = 0;
  let count = 0;
  let match: RegExpExecArray | null;
  while ((match = scoreRegex.exec(details)) !== null) {
    const a = parseInt(match[1]);
    const b = parseInt(match[2]);
    if (a >= 27 && a <= 50 && b >= 27 && b <= 50) {
      totalMargin += Math.abs(a - b);
      count++;
    }
  }
  if (count === 0) return 1.0;
  const avgMargin = totalMargin / count;
  return 1 + 0.1 * Math.min(Math.max(avgMargin - 1, 0), 4);
}

function expectedScore(ratingA: number, ratingB: number): number {
  return 1 / (1 + Math.pow(10, (ratingB - ratingA) / 400));
}

function newRating(
  oldRating: number,
  k: number,
  expected: number,
  actual: number,
): number {
  return Math.round(oldRating + k * (actual - expected));
}

/**
 * Win/loss streak momentum: each consecutive win or loss adds +5% to K, capped at +25%.
 * A fight in a streak of N → multiplier = min(1 + N*0.05, 1.25).
 * Draws reset the counter to 0.
 */
function streakMultiplier(streak: number): number {
  return Math.min(1 + Math.abs(streak) * 0.05, 1.25);
}

/**
 * Ring-rust penalty: applied only when a fighter LOSES their first fight
 * back after a long layoff (>18 months). Amplifies the Elo loss.
 * 18 months → ×1.0 (no change); 36 months → ×1.75 (capped at ×2.0).
 */
const INACTIVITY_THRESHOLD_DAYS = 548; // ~18 months
function inactivityDecayMultiplier(
  lastDate: string,
  currentDate: string,
): number {
  if (!lastDate || !currentDate) return 1.0;
  const days =
    (new Date(currentDate).getTime() - new Date(lastDate).getTime()) /
    86_400_000;
  if (days < INACTIVITY_THRESHOLD_DAYS) return 1.0;
  return Math.min(1 + ((days - INACTIVITY_THRESHOLD_DAYS) / 548) * 0.75, 2.0);
}

/* ------------------------------------------------------------------ */
/*  WEIGHT-CLASS NORMALISER (mirrors the page helper)                 */
/* ------------------------------------------------------------------ */
function normalizeWeightClass(raw: string | null): string {
  if (!raw) return "unknown";
  const cleaned = raw
    .normalize("NFD")
    .replace(/[\u0300-\u036f]/g, "")
    .toLowerCase()
    .replace(/ultimate\s*fighter[\s\d]*/g, "")
    .replace(/\btournament\b/g, "")
    .replace(/\b(title|championship|bout|world|undisputed|interim)\b/g, "")
    .replace(/^ufc\s*/, "")
    .trim();

  // Check women's classes FIRST — must happen before any men's class check
  // to avoid "women's flyweight" matching the "flyweight" keyword below.
  if (/womens?|woman/.test(cleaned)) {
    if (/bantam/.test(cleaned)) return "women_bantamweight";
    if (/feather/.test(cleaned)) return "women_featherweight";
    if (/fly/.test(cleaned)) return "women_flyweight";
    if (/straw/.test(cleaned)) return "women_strawweight";
  }

  // Men's classes — use regex to handle spaced forms like "light heavyweight".
  // light heavyweight MUST come before heavyweight.
  if (/light\s*heavy/.test(cleaned)) return "lightheavyweight";
  if (/heavy/.test(cleaned)) return "heavyweight";
  if (/middle/.test(cleaned)) return "middleweight";
  if (/welter/.test(cleaned)) return "welterweight";
  // lightweight MUST come after light heavyweight (already returned above).
  if (/light/.test(cleaned)) return "lightweight";
  if (/feather/.test(cleaned)) return "featherweight";
  if (/bantam/.test(cleaned)) return "bantamweight";
  if (/fly/.test(cleaned)) return "flyweight";

  // Fallback: collapse all spaces
  return cleaned.replace(/\s+/g, "") || "unknown";
}

/* Map normalised key → display name */
function displayWeightClass(key: string): string {
  const map: Record<string, string> = {
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
    women_featherweight: "Women's Featherweight",
  };
  return map[key] || key;
}

/* ------------------------------------------------------------------ */
/*  MAIN HANDLER                                                      */
/* ------------------------------------------------------------------ */
export async function GET(request: Request) {
  const { searchParams } = new URL(request.url);
  const targetName = searchParams.get("name")?.trim() || null;
  try {
    const supabase = createServerSupabaseClient();

    // 1. Fetch ALL fight results ordered by event id (chronological)
    //    The id column is auto-incrementing and roughly chronological.
    let allFights: any[] = [];
    let offset = 0;
    const pageSize = 1000;
    while (true) {
      const { data, error } = await supabase
        .from("ufc_fight_results")
        .select(
          'id, "BOUT", "OUTCOME", "WEIGHTCLASS", "METHOD", "EVENT", "DETAILS", "ROUND", "TIME"',
        )
        .order("id", { ascending: true })
        .range(offset, offset + pageSize - 1);
      if (error) throw error;
      if (!data || data.length === 0) break;
      allFights = allFights.concat(data);
      if (data.length < pageSize) break;
      offset += pageSize;
    }

    // 2. Fetch event dates for timeline
    let allEvents: any[] = [];
    offset = 0;
    while (true) {
      const { data, error } = await supabase
        .from("ufc_event_details")
        .select('"EVENT", "DATE"')
        .range(offset, offset + pageSize - 1);
      if (error) throw error;
      if (!data || data.length === 0) break;
      allEvents = allEvents.concat(data);
      if (data.length < pageSize) break;
      offset += pageSize;
    }

    const eventDateMap: Record<string, string> = {};
    for (const e of allEvents) {
      if (e.EVENT && e.DATE) {
        eventDateMap[e.EVENT] = e.DATE;
        eventDateMap[e.EVENT.toLowerCase().trim()] = e.DATE;
      }
    }

    // Sort fights chronologically by event date (the id column is NOT chronological)
    const getEventDate = (fight: any): number => {
      const eventName = fight.EVENT || "";
      const dateStr =
        eventDateMap[eventName] ||
        eventDateMap[eventName.toLowerCase().trim()] ||
        "";
      if (!dateStr) return 0;
      return new Date(dateStr).getTime();
    };
    allFights.sort((a: any, b: any) => getEventDate(a) - getEventDate(b));

    // 2b. Fetch fight URLs from ufc_fight_details
    let allFightDetails: any[] = [];
    offset = 0;
    while (true) {
      const { data, error } = await supabase
        .from("ufc_fight_details")
        .select('"EVENT", "BOUT", "URL"')
        .range(offset, offset + pageSize - 1);
      if (error) throw error;
      if (!data || data.length === 0) break;
      allFightDetails = allFightDetails.concat(data);
      if (data.length < pageSize) break;
      offset += pageSize;
    }
    const fightUrlMap: Record<string, string> = {};
    for (const fd of allFightDetails) {
      if (fd.EVENT && fd.BOUT && fd.URL) {
        fightUrlMap[`${fd.EVENT}||${fd.BOUT}`] = fd.URL;
        fightUrlMap[
          `${fd.EVENT.toLowerCase().trim()}||${fd.BOUT.toLowerCase().trim()}`
        ] = fd.URL;
      }
    }

    // 3. Run Elo calculations chronologically
    const ratings: Record<string, number> = {}; // fighter → current elo
    const fightCounts: Record<string, number> = {}; // fighter → number of fights processed so far
    const lastFightDates: Record<string, string> = {}; // fighter → ISO date of most recent fight
    // positive = current win streak length, negative = current loss streak length, 0 = no streak
    const streaks: Record<string, number> = {};
    // Division tracker: current division and consecutive fight count for each fighter
    const divTracker: Record<string, { wc: string; count: number }> = {};
    // Per-fighter history: array of { date, elo, event, opponent, result, weightClass, url }
    const history: Record<
      string,
      Array<{
        date: string;
        elo: number;
        eloBeforeFight?: number;
        event: string;
        opponent: string;
        result: string;
        weightClass: string;
        url: string;
      }>
    > = {};

    for (const fight of allFights) {
      if (!fight.BOUT || !fight.OUTCOME) continue;

      const fighters = fight.BOUT.split(/ vs\.? | vs /i);
      if (fighters.length !== 2) continue;

      const fighterA = fighters[0].trim();
      const fighterB = fighters[1].trim();
      if (!fighterA || !fighterB) continue;

      // Determine outcomes
      const outcomeParts = fight.OUTCOME.split("/");
      let scoreA: number;
      let scoreB: number;

      const method = (fight.METHOD || "").toLowerCase();
      const outcomeRaw = fight.OUTCOME.toLowerCase().trim();

      if (
        method.includes("no contest") ||
        method === "nc" ||
        outcomeRaw === "nc" ||
        outcomeRaw === "nc/nc"
      ) {
        continue; // Skip no contests
      }

      if (outcomeRaw === "d" || outcomeRaw === "d/d" || outcomeRaw === "draw") {
        scoreA = 0.5;
        scoreB = 0.5;
      } else if (outcomeParts.length === 2) {
        const oA = outcomeParts[0].trim().toLowerCase();
        const oB = outcomeParts[1].trim().toLowerCase();
        if (oA === "w" && oB === "l") {
          scoreA = 1;
          scoreB = 0;
        } else if (oA === "l" && oB === "w") {
          scoreA = 0;
          scoreB = 1;
        } else {
          continue; // ambiguous
        }
      } else if (outcomeParts.length === 1) {
        const o = outcomeParts[0].trim().toLowerCase();
        if (o === "w") {
          scoreA = 1;
          scoreB = 0;
        } else if (o === "l") {
          scoreA = 0;
          scoreB = 1;
        } else continue;
      } else {
        continue;
      }

      // Initialise if new
      if (!(fighterA in ratings)) {
        ratings[fighterA] = DEFAULT_ELO;
        fightCounts[fighterA] = 0;
      }
      if (!(fighterB in ratings)) {
        ratings[fighterB] = DEFAULT_ELO;
        fightCounts[fighterB] = 0;
      }

      const fightDate =
        eventDateMap[fight.EVENT || ""] ||
        eventDateMap[(fight.EVENT || "").toLowerCase().trim()] ||
        "";
      const eloA = ratings[fighterA];
      const eloB = ratings[fighterB];
      const eloBeforeFightA = eloA; // pre-fight — used for accurate delta display
      const eloBeforeFightB = eloB;
      const prevDateA = lastFightDates[fighterA] || ""; // last fight date before this one (for inactivity check)
      const prevDateB = lastFightDates[fighterB] || "";
      const wc = normalizeWeightClass(fight.WEIGHTCLASS);
      // Division-change detection: fighter has ≥3 fights in their current division and is now in a different one
      const isDivChangeA = !!(
        divTracker[fighterA] &&
        divTracker[fighterA].wc !== wc &&
        divTracker[fighterA].count >= 3
      );
      const isDivChangeB = !!(
        divTracker[fighterB] &&
        divTracker[fighterB].wc !== wc &&
        divTracker[fighterB].count >= 3
      );
      const isTitleFight =
        (fight.WEIGHTCLASS || "").toLowerCase().startsWith("ufc ") &&
        (fight.WEIGHTCLASS || "").toLowerCase().includes("title") &&
        !(fight.EVENT || "").toLowerCase().includes("ultimate fighter");
      const titleMult = isTitleFight ? 1.5 : 1.0;
      const kA = Math.round(
        getK(fightCounts[fighterA]) *
          titleMult *
          streakMultiplier(streaks[fighterA] ?? 0),
      );
      const kB = Math.round(
        getK(fightCounts[fighterB]) *
          titleMult *
          streakMultiplier(streaks[fighterB] ?? 0),
      );

      const expectedA = expectedScore(eloA, eloB);
      const expectedB = expectedScore(eloB, eloA);

      let newEloA = newRating(eloA, kA, expectedA, scoreA);
      let newEloB = newRating(eloB, kB, expectedB, scoreB);

      // Finish bonus: scaled by round (early finish more impressive) and time (late stop discount)
      const isFinish = /\b(ko|tko|submission)\b/i.test(method);
      if (isFinish) {
        const finishBonus = Math.round(
          FINISH_BONUS * roundFinishMultiplier(fight.ROUND),
        );
        if (scoreA === 1) {
          newEloA += finishBonus;
          newEloB -= finishBonus;
        } else if (scoreB === 1) {
          newEloB += finishBonus;
          newEloA -= finishBonus;
        }
      }

      // Decision modifiers: margin-of-victory × split/majority discount
      const isDecision = method.includes("decision");
      if (isDecision && scoreA !== 0.5) {
        const combinedMult =
          decisionMarginMultiplier(fight.DETAILS) *
          splitDecisionMultiplier(method);
        if (combinedMult !== 1.0) {
          const changeA = newEloA - eloA;
          const changeB = newEloB - eloB;
          newEloA = Math.round(eloA + changeA * combinedMult);
          newEloB = Math.round(eloB + changeB * combinedMult);
        }
      }

      // Ring-rust penalty: amplify Elo loss when a fighter loses their first fight
      // back after an 18-month+ layoff. Only the loser is penalised; winner unchanged.
      if (scoreA === 0 && fightDate && prevDateA) {
        const inactMult = inactivityDecayMultiplier(prevDateA, fightDate);
        if (inactMult > 1.0)
          newEloA = Math.round(eloA + (newEloA - eloA) * inactMult);
      }
      if (scoreB === 0 && fightDate && prevDateB) {
        const inactMult = inactivityDecayMultiplier(prevDateB, fightDate);
        if (inactMult > 1.0)
          newEloB = Math.round(eloB + (newEloB - eloB) * inactMult);
      }

      // Division-change protection: fighters established in a division (≥3 fights) lose
      // 35% less Elo when moving to a different division and losing that fight.
      const DIV_CHANGE_PROTECT = 0.65;
      if (isDivChangeA && scoreA === 0)
        newEloA = Math.round(eloA + (newEloA - eloA) * DIV_CHANGE_PROTECT);
      if (isDivChangeB && scoreB === 0)
        newEloB = Math.round(eloB + (newEloB - eloB) * DIV_CHANGE_PROTECT);

      ratings[fighterA] = newEloA;
      ratings[fighterB] = newEloB;
      fightCounts[fighterA] = (fightCounts[fighterA] || 0) + 1;
      fightCounts[fighterB] = (fightCounts[fighterB] || 0) + 1;

      // Update win/loss streaks (positive = win streak length, negative = loss streak length)
      streaks[fighterA] =
        scoreA === 1
          ? Math.max(0, streaks[fighterA] ?? 0) + 1
          : scoreA === 0
            ? Math.min(0, streaks[fighterA] ?? 0) - 1
            : 0; // draw resets streak
      streaks[fighterB] =
        scoreB === 1
          ? Math.max(0, streaks[fighterB] ?? 0) + 1
          : scoreB === 0
            ? Math.min(0, streaks[fighterB] ?? 0) - 1
            : 0;

      // Update division tracker for each fighter
      if (!divTracker[fighterA]) divTracker[fighterA] = { wc, count: 1 };
      else if (divTracker[fighterA].wc === wc) divTracker[fighterA].count++;
      else divTracker[fighterA] = { wc, count: 1 };

      if (!divTracker[fighterB]) divTracker[fighterB] = { wc, count: 1 };
      else if (divTracker[fighterB].wc === wc) divTracker[fighterB].count++;
      else divTracker[fighterB] = { wc, count: 1 };

      const eventName = fight.EVENT || "";
      const date = fightDate || "";
      if (date) {
        lastFightDates[fighterA] = date;
        lastFightDates[fighterB] = date;
      }
      const fightUrl =
        fightUrlMap[`${eventName}||${fight.BOUT}`] ||
        fightUrlMap[
          `${eventName.toLowerCase().trim()}||${fight.BOUT.toLowerCase().trim()}`
        ] ||
        "";

      // Record history for both fighters
      if (!history[fighterA]) history[fighterA] = [];
      history[fighterA].push({
        date,
        elo: ratings[fighterA],
        eloBeforeFight: eloBeforeFightA,
        event: eventName,
        opponent: fighterB,
        result: scoreA === 1 ? "W" : scoreA === 0 ? "L" : "D",
        weightClass: wc,
        url: fightUrl,
      });

      if (!history[fighterB]) history[fighterB] = [];
      history[fighterB].push({
        date,
        elo: ratings[fighterB],
        eloBeforeFight: eloBeforeFightB,
        event: eventName,
        opponent: fighterA,
        result: scoreB === 1 ? "W" : scoreB === 0 ? "L" : "D",
        weightClass: wc,
        url: fightUrl,
      });
    }

    // Per-fighter lookup: return just this fighter's data without full batch processing
    if (targetName) {
      const norm = (s: string) =>
        s
          .normalize("NFD")
          .replace(/[\u0300-\u036f]/g, "")
          .replace(/[.'\u2011]/g, "")
          .toLowerCase()
          .trim();
      const tNorm = norm(targetName);
      const canonName =
        Object.keys(ratings).find((n) => norm(n) === tNorm) || targetName;
      return NextResponse.json({
        name: canonName,
        currentElo: ratings[canonName] ?? null,
        history: history[canonName] || [],
      });
    }

    // 4. Compute summary stats per fighter
    const fighterSummaries: Array<{
      name: string;
      currentElo: number;
      peakElo: number;
      peakDate: string;
      fights: number;
      primaryWeightClass: string;
      lastFightDate: string;
      history: (typeof history)[string];
    }> = [];

    for (const [name, hist] of Object.entries(history)) {
      if (hist.length < 5) continue; // Skip fighters with very few fights

      let peakElo = DEFAULT_ELO;
      let peakDate = "";
      const wcCounts: Record<string, number> = {};

      for (const h of hist) {
        if (h.elo > peakElo) {
          peakElo = h.elo;
          peakDate = h.date;
        }
        wcCounts[h.weightClass] = (wcCounts[h.weightClass] || 0) + 1;
      }

      // Primary weight class = most common
      const primaryWC =
        Object.entries(wcCounts).sort((a, b) => b[1] - a[1])[0]?.[0] ||
        "unknown";

      fighterSummaries.push({
        name,
        currentElo: ratings[name],
        peakElo,
        peakDate,
        fights: hist.length,
        primaryWeightClass: primaryWC,
        lastFightDate: hist[hist.length - 1]?.date || "",
        history: hist,
      });
    }

    // Sort by peak Elo descending
    fighterSummaries.sort((a, b) => b.peakElo - a.peakElo);

    // Valid weight classes only (exclude openweight, superfight, catch-weight, etc.)
    const VALID_WEIGHT_CLASSES = new Set([
      "heavyweight",
      "lightheavyweight",
      "middleweight",
      "welterweight",
      "lightweight",
      "featherweight",
      "bantamweight",
      "flyweight",
      "women_bantamweight",
      "women_flyweight",
      "women_strawweight",
      "women_featherweight",
    ]);

    // 5. Build weight-class grouped top fighters (top 20 per valid class)
    const byWeightClass: Record<string, typeof fighterSummaries> = {};
    for (const f of fighterSummaries) {
      const wc = f.primaryWeightClass;
      if (!VALID_WEIGHT_CLASSES.has(wc)) continue;
      if (!byWeightClass[wc]) byWeightClass[wc] = [];
      byWeightClass[wc].push(f);
    }

    // Weight class display name map
    const wcDisplayNames: Record<string, string> = {};
    for (const key of Object.keys(byWeightClass)) {
      wcDisplayNames[key] = displayWeightClass(key);
    }

    return NextResponse.json({
      // Top 100 overall by peak Elo (UI slices to 20 after sorting by chosen column)
      candidatePool: fighterSummaries.slice(0, 100),
      // Per weight class top 20
      byWeightClass,
      wcDisplayNames,
      // All fighters (≥5 UFC fights) for distribution analytics
      allFighters: fighterSummaries.map((f) => ({
        name: f.name,
        currentElo: f.currentElo,
        primaryWeightClass: f.primaryWeightClass,
        lastFightDate: f.lastFightDate,
      })),
    });
  } catch (err: any) {
    console.error("UFC Elo API Error:", err.message);
    return NextResponse.json({ error: err.message }, { status: 500 });
  }
}
