import { createClient } from '@supabase/supabase-js'
import { NextResponse } from 'next/server'

export const revalidate = 43200 // 12h cache

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------
interface FightRow {
  id: number
  BOUT: string | null
  OUTCOME: string | null
  WEIGHTCLASS: string | null
  EVENT: string | null
}

interface Classified extends FightRow {
  _date: string                            // ISO "YYYY-MM-DD"
  _key: string                             // e.g. "heavyweight"
  _isInterim: boolean
}

interface ChampEntry {
  name: string
  reignStart: string                       // ISO date
  reignEnd: string | null                  // ISO date, "Present", or null (= "Present")
  notes?: string
  imageUrl?: string
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
function norm(s: string | null | undefined): string {
  if (!s) return ''
  return s
    .normalize('NFD')
    .replace(/[\u0300-\u036f]/g, '')
    .replace(/[.'‑]/g, '')
    .toLowerCase()
    .trim()
}

function getWinner(bout: string | null, outcome: string | null): string | null {
  if (!bout || !outcome) return null
  const fighters = bout.split(/\s+vs\.?\s+/i)
  if (fighters.length !== 2) return null
  const parts = outcome.split('/')
  if (parts[0]?.toLowerCase() === 'w') return fighters[0].trim()
  if (parts.length > 1 && parts[1]?.toLowerCase() === 'w') return fighters[1].trim()
  return null // NC, Draw
}

function parseDateToISO(dateStr: string | null): string | null {
  if (!dateStr) return null
  try {
    const d = new Date(dateStr)
    if (isNaN(d.getTime())) return null
    return d.toISOString().split('T')[0]
  } catch { return null }
}

/**
 * Maps a WEIGHTCLASS string to a canonical division key.
 * Returns null for TUF tournament fights or unrecognised classes.
 */
function classifyWeightClass(
  raw: string | null,
): { key: string; isInterim: boolean } | null {
  if (!raw) return null
  const lc = raw.toLowerCase()

  // Exclude TUF / tournament title fights
  if (lc.includes('ultimate fighter') || lc.includes('tournament')) return null

  const isInterim = lc.includes('interim')
  const isWomen = lc.includes("women")

  // Women's divisions (check before men's to avoid "bantamweight" matching men's)
  if (isWomen && lc.includes('strawweight'))  return { key: 'women_strawweight',  isInterim }
  if (isWomen && lc.includes('flyweight'))    return { key: 'women_flyweight',    isInterim }
  if (isWomen && lc.includes('bantamweight')) return { key: 'women_bantamweight', isInterim }
  if (isWomen) return null

  // Men's – "light heavyweight" MUST come before "heavyweight" (substring overlap)
  if (lc.includes('light heavyweight')) return { key: 'lightheavyweight', isInterim }
  if (lc.includes('heavyweight'))       return { key: 'heavyweight',      isInterim }
  if (lc.includes('middleweight'))      return { key: 'middleweight',     isInterim }
  if (lc.includes('welterweight'))      return { key: 'welterweight',     isInterim }
  if (lc.includes('lightweight'))       return { key: 'lightweight',      isInterim }
  if (lc.includes('featherweight'))     return { key: 'featherweight',    isInterim }
  if (lc.includes('bantamweight'))      return { key: 'bantamweight',     isInterim }
  if (lc.includes('flyweight'))         return { key: 'flyweight',        isInterim }

  return null
}

const DIVISION_NAME: Record<string, string> = {
  heavyweight:        'Heavyweight (206-265 lbs)',
  lightheavyweight:   'Light Heavyweight (205 lbs)',
  middleweight:       'Middleweight (185 lbs)',
  welterweight:       'Welterweight (170 lbs)',
  lightweight:        'Lightweight (155 lbs)',
  featherweight:      'Featherweight (145 lbs)',
  bantamweight:       'Bantamweight (135 lbs)',
  flyweight:          'Flyweight (125 lbs)',
  women_strawweight:  "Women's Strawweight (115 lbs)",
  women_flyweight:    "Women's Flyweight (125 lbs)",
  women_bantamweight: "Women's Bantamweight (135 lbs)",
}

// ---------------------------------------------------------------------------
// Core lineage-building algorithm
// ---------------------------------------------------------------------------
/**
 * Walks through all title fights for ONE division (already sorted oldest→newest)
 * and builds a list of champion reigns.
 *
 * Key edge case handled: when the undisputed champion vacates/is stripped and
 * the interim champion is promoted WITHOUT a fight.  We detect this when the
 * current undisputed champion is absent from the next undisputed title fight
 * AND the current interim champion IS in that fight.
 */
function buildLineage(fights: Classified[]): ChampEntry[] {
  let undisputed: { name: string; reignStart: string } | null = null
  let interim: { name: string; reignStart: string } | null = null

  const lineage: ChampEntry[] = []

  for (const fight of fights) {
    const winner = getWinner(fight.BOUT, fight.OUTCOME)
    const fighters = (fight.BOUT ?? '').split(/\s+vs\.?\s+/i).map(s => s.trim())
    const date = fight._date

    // ------------------------------------------------------------------ INTERIM
    if (fight._isInterim) {
      if (!winner) continue
      if (interim && norm(interim.name) === norm(winner)) {
        // Same interim champ defending – keep original reignStart
      } else {
        interim = { name: winner, reignStart: date }
      }
      continue
    }

    // ---------------------------------------------------------------- UNDISPUTED – first ever
    if (!undisputed) {
      if (!winner) continue
      undisputed = { name: winner, reignStart: date }
      lineage.push({ name: winner, reignStart: date, reignEnd: null })
      continue
    }

    const champInFight = fighters.some(f => norm(f) === norm(undisputed!.name))

    if (champInFight) {
      if (!winner) continue  // NC / Draw – champ retains

      if (norm(winner) === norm(undisputed.name)) {
        // ---- Successful defense ----
        // If the interim champ was the opponent, they lost the unification bout → clear interim
        if (interim) {
          const opponent = fighters.find(f => norm(f) !== norm(undisputed!.name))
          if (opponent && norm(opponent) === norm(interim.name)) interim = null
        }
        continue
      }

      // ---- Champion LOST ----
      if (lineage.length > 0) lineage[lineage.length - 1].reignEnd = date
      undisputed = { name: winner, reignStart: date }
      lineage.push({ name: winner, reignStart: date, reignEnd: null })
      // If the new champ was the interim champ (unified), clear interim
      if (interim && norm(winner) === norm(interim.name)) interim = null

    } else {
      // ---- Current champion NOT in this fight → VACANCY ----
      if (lineage.length > 0) {
        lineage[lineage.length - 1].reignEnd = date
        lineage[lineage.length - 1].notes = 'vacated'
      }

      const interimInFight =
        interim !== null && fighters.some(f => norm(f) === norm(interim!.name))

      if (interimInFight) {
        // Interim champ was promoted and is fighting here as undisputed champion.
        // Use this fight's date as reignStart so there is NO overlap with the
        // previous champion (whose reignEnd is also set to `date` above).
        const promotedName = interim!.name

        if (!winner || norm(winner) === norm(promotedName)) {
          // NC/Draw OR win – promoted champ retains
          undisputed = { name: promotedName, reignStart: date }
          lineage.push({ name: promotedName, reignStart: date, reignEnd: null, notes: 'promoted from interim' })
        } else {
          // Promoted champ lost their first undisputed defence (rare)
          lineage.push({ name: promotedName, reignStart: date, reignEnd: date, notes: 'promoted from interim' })
          undisputed = { name: winner, reignStart: date }
          lineage.push({ name: winner, reignStart: date, reignEnd: null })
        }
        interim = null

      } else if (winner) {
        // Vacated belt – new champion crowned via a regular fight
        undisputed = { name: winner, reignStart: date }
        lineage.push({ name: winner, reignStart: date, reignEnd: null })
      }
      // NC on vacant belt → title stays vacant; nothing to push
    }
  }

  // Mark the last entry as current champion
  if (lineage.length > 0 && !lineage[lineage.length - 1].reignEnd) {
    lineage[lineage.length - 1].reignEnd = 'Present'
  }

  return lineage
}

// ---------------------------------------------------------------------------
// Route handler
// ---------------------------------------------------------------------------
export async function GET() {
  try {
    const supabaseUrl     = process.env.SUPABASE_URL
    const supabaseAnonKey = process.env.SUPABASE_ANON_KEY
    if (!supabaseUrl || !supabaseAnonKey) {
      return NextResponse.json({ error: 'Supabase env vars missing' }, { status: 500 })
    }
    const supabase = createClient(supabaseUrl, supabaseAnonKey)

    // Fetch every fight that has "title" in the weight-class column
    const { data: fights, error: fightsError } = await supabase
      .from('ufc_fight_results')
      .select('id, BOUT, OUTCOME, WEIGHTCLASS, EVENT')
      .ilike('WEIGHTCLASS', '%title%')
      .order('id', { ascending: true }) // low id = most recent in this DB

    if (fightsError) throw fightsError

    // Collect unique event names for the date lookup
    const eventNames = [
      ...new Set(
        (fights ?? []).map(f => f.EVENT?.trim()).filter((e): e is string => !!e),
      ),
    ]

    const { data: eventDetails } = await supabase
      .from('ufc_event_details')
      .select('EVENT, DATE')
      .in('EVENT', eventNames)

    const dateMap = new Map(
      (eventDetails ?? []).map(e => [e.EVENT?.trim() ?? '', e.DATE as string]),
    )

    // Classify each fight and attach its ISO date
    const classified: Classified[] = []
    for (const f of fights ?? []) {
      const cls = classifyWeightClass(f.WEIGHTCLASS)
      if (!cls) continue
      const rawDate = f.EVENT ? (dateMap.get(f.EVENT.trim()) ?? null) : null
      const isoDate = parseDateToISO(rawDate)
      if (!isoDate) continue
      classified.push({ ...f, _date: isoDate, _key: cls.key, _isInterim: cls.isInterim })
    }

    // Group by division, sort each group oldest → newest, then build lineage
    const byDivision = new Map<string, Classified[]>()
    for (const f of classified) {
      if (!byDivision.has(f._key)) byDivision.set(f._key, [])
      byDivision.get(f._key)!.push(f)
    }

    const lineage: Record<string, { displayName: string; champions: ChampEntry[] }> = {}

    for (const [key, divFights] of byDivision) {
      // Sort chronologically (ascending date string sorts correctly for ISO dates)
      divFights.sort((a, b) => a._date.localeCompare(b._date))
      const champions = buildLineage(divFights)
      if (champions.length > 0) {
        lineage[key] = { displayName: DIVISION_NAME[key] ?? key, champions }
      }
    }

    // Attach fighter image URLs from ufc_champion_images
    const { data: allImages } = await supabase
      .from('ufc_champion_images')
      .select('FIRST, LAST, MAINSHOT')

    const imageMap = new Map<string, string>()
    for (const img of allImages ?? []) {
      if (img.FIRST && img.LAST && img.MAINSHOT) {
        imageMap.set(norm(`${img.FIRST} ${img.LAST}`), img.MAINSHOT as string)
      }
    }

    for (const divData of Object.values(lineage)) {
      for (const champ of divData.champions) {
        const url = imageMap.get(norm(champ.name))
        if (url) champ.imageUrl = url
      }
    }

    return NextResponse.json({ lineage, fetchedAt: new Date().toISOString() })
  } catch (err: unknown) {
    const msg = err instanceof Error ? err.message : String(err)
    console.error('[ufc-champions/lineage]', msg)
    return NextResponse.json({ error: msg }, { status: 500 })
  }
}
