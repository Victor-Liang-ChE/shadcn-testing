import { NextResponse } from 'next/server'
import { createServerSupabaseClient } from '@/lib/supabaseServer'
import {
  predictFight,
  normalizeName,
  type FightInput,
  type FighterStats,
  type EventPrediction,
} from '@/lib/ufc-model'

// ─── Constants ───
const BUCKET = 'ufc-models'
const MODEL_CACHE_TTL_MS = 60 * 60 * 1000 // 1 hour

interface ModelMeta {
  features: string[]
  hp_label: string
  training_rows: number
}

// ─── Global Cache ───
let modelCache: {
  cb: any
  lg: any
  meta: ModelMeta
  source: 'storage'
  expiresAt: number
} | null = null

async function loadModels(
  supabase: ReturnType<typeof createServerSupabaseClient>
): Promise<NonNullable<typeof modelCache>> {
  if (modelCache && Date.now() < modelCache.expiresAt) return modelCache

  try {
    // Attempt download from Supabase Storage — no local JSON files anymore
    const [cbRes, lgRes, metaRes] = await Promise.all([
      supabase.storage.from(BUCKET).download('ufc-model-catboost.json'),
      supabase.storage.from(BUCKET).download('ufc-model-lightgbm.json'),
      supabase.storage.from(BUCKET).download('ufc-model-meta.json'),
    ])

    if (cbRes.error || !cbRes.data || lgRes.error || !lgRes.data || metaRes.error || !metaRes.data) {
      console.error('Model download failed:', { cb: cbRes.error, lg: lgRes.error, meta: metaRes.error })
      throw new Error('Supabase Storage models unavailable')
    }

    const [cb, lg, meta] = await Promise.all([
      cbRes.data.text().then(JSON.parse),
      lgRes.data.text().then(JSON.parse),
      metaRes.data.text().then(JSON.parse),
    ])

    modelCache = {
      cb,
      lg,
      meta,
      source: 'storage', // Always 'storage' now that we removed local bundles
      expiresAt: Date.now() + MODEL_CACHE_TTL_MS,
    }
    return modelCache
  } catch (err) {
    console.error('Error loading models:', err)
    throw err // Re-throw so the API returns a 500 until storage is fixed
  }
}

/** Paginated fetch from Supabase (1000 rows per page) */
async function fetchAll(
  supabase: ReturnType<typeof createServerSupabaseClient>,
  table: string,
  columns: string,
  orderBy = 'id'
) {
  const rows: Record<string, unknown>[] = []
  let offset = 0
  const pageSize = 1000
  while (true) {
    const { data, error } = await supabase
      .from(table)
      .select(columns)
      .order(orderBy, { ascending: true })
      .range(offset, offset + pageSize - 1)
    if (error) throw new Error(`${table}: ${error.message}`)
    if (!data || data.length === 0) break
    rows.push(...(data as unknown as Record<string, unknown>[]))
    if (data.length < pageSize) break
    offset += pageSize
  }
  return rows
}

/** Normalise a raw weight class string from ufc_fight_results to a canonical key. */
function extractWeightClass(raw: string | null | undefined): string | undefined {
  if (!raw) return undefined
  const s = raw.toLowerCase().normalize('NFD').replace(/[\u0300-\u036f]/g, '')
  // Women's classes must be checked before men's to avoid mismatches
  if (/womens?|woman/.test(s)) {
    if (/bantam/.test(s)) return 'women_bantamweight'
    if (/feather/.test(s)) return 'women_featherweight'
    if (/fly/.test(s)) return 'women_flyweight'
    if (/straw/.test(s)) return 'women_strawweight'
    return undefined
  }
  if (/light\s*heavy/.test(s)) return 'lightheavyweight'
  if (/heavy/.test(s)) return 'heavyweight'
  if (/middle/.test(s)) return 'middleweight'
  if (/welter/.test(s)) return 'welterweight'
  if (/light/.test(s)) return 'lightweight'
  if (/feather/.test(s)) return 'featherweight'
  if (/bantam/.test(s)) return 'bantamweight'
  if (/fly/.test(s)) return 'flyweight'
  return undefined
}

/**
 * Build a live fighter database from ufc_fight_results + ufc_event_details.
 * Computes current win streaks by walking through fights chronologically.
 * This replaces the static fighter_db.json.
 */
async function buildFighterDb(
  supabase: ReturnType<typeof createServerSupabaseClient>
): Promise<Record<string, FighterStats>> {
  // Fetch fight results and event dates in parallel
  const [fights, events] = await Promise.all([
    fetchAll(supabase, 'ufc_fight_results', 'id, "BOUT", "OUTCOME", "METHOD", "EVENT", "WEIGHTCLASS"'),
    fetchAll(supabase, 'ufc_event_details', '"EVENT", "DATE"'),
  ])

  // Build event → date map
  const eventDateMap: Record<string, string> = {}
  for (const e of events) {
    const ev = e.EVENT as string
    const dt = e.DATE as string
    if (ev && dt) {
      eventDateMap[ev] = dt
      eventDateMap[ev.toLowerCase().trim()] = dt
    }
  }

  // Sort fights chronologically
  fights.sort((a, b) => {
    const dateA = eventDateMap[a.EVENT as string] || eventDateMap[(a.EVENT as string || '').toLowerCase().trim()] || ''
    const dateB = eventDateMap[b.EVENT as string] || eventDateMap[(b.EVENT as string || '').toLowerCase().trim()] || ''
    return new Date(dateA).getTime() - new Date(dateB).getTime()
  })

  // Walk through fights chronologically to compute streaks
  const streaks: Record<string, number> = {} // fighter → current win streak (0 = reset)
  const wins: Record<string, number> = {}
  const losses: Record<string, number> = {}
  const weightClasses: Record<string, string> = {}

  for (const fight of fights) {
    const bout = fight.BOUT as string
    const outcome = fight.OUTCOME as string
    const method = ((fight.METHOD as string) || '').toLowerCase()
    if (!bout || !outcome) continue

    const parts = bout.split(/ vs\.? | vs /i)
    if (parts.length !== 2) continue

    const a = parts[0].trim()
    const b = parts[1].trim()
    if (!a || !b) continue

    // Skip no contests
    const oc = outcome.toLowerCase().trim()
    if (method.includes('no contest') || oc === 'nc' || oc === 'nc/nc') continue

    // Determine outcome
    const outParts = outcome.split('/')
    let aWon: boolean | null = null
    if (oc === 'd' || oc === 'd/d' || oc === 'draw') {
      aWon = null // draw resets both
    } else if (outParts.length === 2) {
      const oA = outParts[0].trim().toLowerCase()
      const oB = outParts[1].trim().toLowerCase()
      if (oA === 'w' && oB === 'l') aWon = true
      else if (oA === 'l' && oB === 'w') aWon = false
      else continue
    } else if (outParts.length === 1) {
      const o = outParts[0].trim().toLowerCase()
      if (o === 'w') aWon = true
      else if (o === 'l') aWon = false
      else continue
    } else continue

    const keyA = normalizeName(a)
    const keyB = normalizeName(b)

    // Track most-recent weight class for both fighters (chronological order means last seen wins)
    const wc = extractWeightClass(fight.WEIGHTCLASS as string | null)
    if (wc) {
      weightClasses[keyA] = wc
      weightClasses[keyB] = wc
    }

    if (aWon === null) {
      // Draw: reset both to 0
      streaks[keyA] = 0
      streaks[keyB] = 0
    } else if (aWon) {
      // keyA won, keyB lost
      streaks[keyA] = Math.max(0, streaks[keyA] ?? 0) + 1
      streaks[keyB] = Math.min(0, streaks[keyB] ?? 0) - 1
      wins[keyA] = (wins[keyA] ?? 0) + 1
      losses[keyB] = (losses[keyB] ?? 0) + 1
    } else {
      // keyB won, keyA lost
      streaks[keyB] = Math.max(0, streaks[keyB] ?? 0) + 1
      streaks[keyA] = Math.min(0, streaks[keyA] ?? 0) - 1
      wins[keyB] = (wins[keyB] ?? 0) + 1
      losses[keyA] = (losses[keyA] ?? 0) + 1
    }
  }

  // Build the fighter DB
  const db: Record<string, FighterStats> = {}
  const allKeys = new Set([...Object.keys(streaks), ...Object.keys(wins), ...Object.keys(losses)])
  for (const key of allKeys) {
    db[key] = {
      name: key,
      win_streak: streaks[key] ?? 0,
      wins: wins[key] ?? 0,
      losses: losses[key] ?? 0,
      weightClass: weightClasses[key],
    }
  }
  return db
}

export async function GET() {
  try {
    const supabase = createServerSupabaseClient()

    // Load models + fight odds + fighter DB + event dates all in parallel
    const [models, oddsResult, fighterDb, eventDatesRes] = await Promise.all([
      loadModels(supabase),
      supabase
        .from('ufc_betting_odds')
        .select('*')
        .order('id', { ascending: true }),
      buildFighterDb(supabase),
      fetchAll(supabase, 'ufc_event_details', '"EVENT", "DATE"'),
    ])

    // Build event name → date fallback map
    const eventDates: Record<string, string> = {}
    for (const e of eventDatesRes) {
      const ev = (e.EVENT as string || '').trim()
      const dt = e.DATE as string
      if (ev && dt) {
        eventDates[ev.toLowerCase()] = dt
      }
    }

    const { cbModel, lgModel, meta } = { cbModel: models.cb, lgModel: models.lg, meta: models.meta }

    const { data: oddsRows, error } = oddsResult
    if (error) throw new Error(`Supabase error: ${error.message}`)

    // Detect prop bets by name — fallback for rows that predate the is_prop column
    const PROP_RE = /\d|\bover\b|\bunder\b|\bdecision\b|ko\/tko|\btko\b|\bsubmission\b|\brds?\b|\bdistance\b|\binside\b|\bfinish\b|\bstarts?\b|\bwon'?t\b|\bdoesn'?t\b|\bgoes?\b|\bdraw\b|\bscorecards?\b|\baction\b|\bfight\b|\bno action\b/i
    const isRowProp = (row: Record<string, unknown>) =>
      row.is_prop === true ||
      PROP_RE.test(row.r_fighter as string) ||
      PROP_RE.test(row.b_fighter as string)

    // Group by event
    type PropRow = { rName: string; bName: string; rOdds: number | null; bOdds: number | null; rPerBook: Record<string, number>; bPerBook: Record<string, number> }
    const eventMap = new Map<
      string,
      { name: string; url?: string; date?: string; fights: FightInput[]; props: PropRow[] }
    >()

    for (const row of oddsRows ?? []) {
      const key = row.event_name as string
      if (!eventMap.has(key)) {
        const rawDate = row.event_date as string | undefined
        const lookupDate = eventDates[(row.event_name as string || '').toLowerCase().trim()]
        eventMap.set(key, {
          name: row.event_name as string,
          url: row.event_url as string | undefined,
          date: rawDate || lookupDate,
          fights: [],
          props: [],
        })
      }
      if (isRowProp(row)) {
        eventMap.get(key)!.props.push({
          rName: row.r_fighter as string,
          bName: row.b_fighter as string,
          rOdds: row.r_odds as number | null,
          bOdds: row.b_odds as number | null,
          rPerBook: (row.r_per_book as Record<string, number>) ?? {},
          bPerBook: (row.b_per_book as Record<string, number>) ?? {},
        })
      } else {
        eventMap.get(key)!.fights.push({
          rFighter: row.r_fighter as string,
          bFighter: row.b_fighter as string,
          rOdds: row.r_odds as number | null,
          bOdds: row.b_odds as number | null,
          rPerBook: (row.r_per_book as Record<string, number>) ?? {},
          bPerBook: (row.b_per_book as Record<string, number>) ?? {},
          rEvolution: row.r_evolution ?? undefined,
          bEvolution: row.b_evolution ?? undefined,
          matchupId: row.matchup_id as number | undefined,
        })
      }
    }

    // Run predictions for main fights only (props don't have fighter stats)
    const todayStr = new Date().toISOString().slice(0, 10) // e.g. "2026-03-24"
    const events: (EventPrediction & { props: PropRow[]; date?: string })[] = []
    for (const [, event] of eventMap) {
      // Skip events that have already passed
      if (event.date && event.date < todayStr) continue
      const predictions = event.fights.map((fight) =>
        predictFight(cbModel, lgModel, fight, fighterDb)
      )
      events.push({
        name: event.name,
        url: event.url,
        date: event.date,
        fights: predictions,
        activeBets: predictions.filter((p) => p.betSignal === 'bet').length,
        props: event.props,
      })
    }

    const totalBets = events.reduce((sum, e) => sum + e.activeBets, 0)
    const totalFights = events.reduce((sum, e) => sum + e.fights.length, 0)

    const scrapedAt =
      oddsRows && oddsRows.length > 0
        ? oddsRows[oddsRows.length - 1].scraped_at
        : null

    return NextResponse.json(
      {
        generatedAt: new Date().toISOString(),
        model: {
          name: 'Dual CatBoost + LightGBM',
          features: meta.features,
          hp: meta.hp_label,
          trainingFights: meta.training_rows,
          source: models.source,
        },
        events,
        totalBets,
        totalFights,
        oddsScrapedAt: scrapedAt,
      },
      {
        headers: {
          'Cache-Control': 'public, s-maxage=300, stale-while-revalidate=600',
        },
      }
    )
  } catch (err: unknown) {
    const message = err instanceof Error ? err.message : String(err)
    console.error('[ufc-predictions] Error:', message)
    return NextResponse.json({ error: message }, { status: 500 })
  }
}
