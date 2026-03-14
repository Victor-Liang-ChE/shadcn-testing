import { createClient } from '@supabase/supabase-js'
import { NextResponse } from 'next/server'

export const revalidate = 43200 // 12h cache

function parseWeightClass(raw: string | null): { key: string; display: string; isInterim: boolean } | null {
  if (!raw) return null
  const lower = raw.toLowerCase()
  const isInterim = lower.includes('interim')
  const clean = lower
    .replace(/ufc\s*/g, '')
    .replace(/interim\s*/g, '')
    .replace(/title\s*bout\s*/g, '')
    .replace(/championship\s*/g, '')
    .trim()

  if (lower.includes("women") && clean.includes("strawweight"))  return { key: 'women_strawweight',  display: "Women's Strawweight", isInterim }
  if (lower.includes("women") && clean.includes("flyweight"))    return { key: 'women_flyweight',    display: "Women's Flyweight",    isInterim }
  if (lower.includes("women") && clean.includes("bantamweight")) return { key: 'women_bantamweight', display: "Women's Bantamweight", isInterim }
  if (clean.includes('flyweight'))         return { key: 'flyweight',        display: 'Flyweight',        isInterim }
  if (clean.includes('bantamweight'))      return { key: 'bantamweight',     display: 'Bantamweight',     isInterim }
  if (clean.includes('featherweight'))     return { key: 'featherweight',    display: 'Featherweight',    isInterim }
  if (clean.includes('lightweight'))       return { key: 'lightweight',      display: 'Lightweight',      isInterim }
  if (clean.includes('welterweight'))      return { key: 'welterweight',     display: 'Welterweight',     isInterim }
  if (clean.includes('middleweight'))      return { key: 'middleweight',     display: 'Middleweight',     isInterim }
  if (clean.includes('light heavyweight')) return { key: 'lightheavyweight', display: 'Light Heavyweight', isInterim }
  if (clean.includes('heavyweight'))       return { key: 'heavyweight',      display: 'Heavyweight',      isInterim }
  return null
}

// Parse winner from "Fighter A vs. Fighter B" + "W/L" outcome
function getWinner(bout: string | null, outcome: string | null): string | null {
  if (!bout || !outcome) return null
  const fighters = bout.split(/\s+vs\.?\s+/i)
  if (fighters.length !== 2) return null
  const parts = outcome.split('/')
  if (parts[0]?.toLowerCase() === 'w') return fighters[0].trim()
  if (parts.length > 1 && parts[1]?.toLowerCase() === 'w') return fighters[1].trim()
  return null
}

// Convert "Month DD, YYYY" to "YYYY-MM-DD"
function parseDateString(dateStr: string | null): string | null {
  if (!dateStr) return null
  try {
    const d = new Date(dateStr)
    if (isNaN(d.getTime())) return null
    return d.toISOString().split('T')[0] // "YYYY-MM-DD"
  } catch { return null }
}

export async function GET() {
  try {
    const supabaseUrl = process.env.SUPABASE_URL
    const supabaseAnonKey = process.env.SUPABASE_ANON_KEY
    if (!supabaseUrl || !supabaseAnonKey) {
      return NextResponse.json({ error: 'Supabase env vars missing' }, { status: 500 })
    }
    const supabase = createClient(supabaseUrl, supabaseAnonKey)

    // Fetch all title fights. Lower id = more recent event in this DB.
    const { data: titleFights, error: fightsError } = await supabase
      .from('ufc_fight_results')
      .select('id, BOUT, OUTCOME, WEIGHTCLASS, EVENT')
      .ilike('WEIGHTCLASS', '%title%')
      .order('id', { ascending: true }) // lowest id = most recent

    if (fightsError) throw fightsError

    // Collect unique event names from those fights for date lookup
    const eventNames = [...new Set((titleFights ?? []).map(f => f.EVENT?.trim()).filter(Boolean))]

    // Fetch event dates (batch)
    const { data: eventDetails } = await supabase
      .from('ufc_event_details')
      .select('EVENT, DATE')
      .in('EVENT', eventNames)

    const dateMap = new Map<string, string>()
    for (const e of eventDetails ?? []) {
      if (e.EVENT && e.DATE) dateMap.set(e.EVENT.trim(), e.DATE)
    }

    // Scan newest-first (ASC by id) — first winner per class+type = current holder
    const champByClass = new Map<string, {
      name: string; event: string | null; date: string | null; isInterim: boolean; display: string
    }>()

    for (const fight of titleFights ?? []) {
      const wc = parseWeightClass(fight.WEIGHTCLASS)
      if (!wc) continue

      const mapKey = `${wc.key}:${wc.isInterim ? 'interim' : 'real'}`
      if (champByClass.has(mapKey)) continue

      const winner = getWinner(fight.BOUT, fight.OUTCOME)
      if (!winner) continue

      const eventName = fight.EVENT?.trim() ?? null
      const rawDate = eventName ? (dateMap.get(eventName) ?? null) : null
      const isoDate = parseDateString(rawDate) // "YYYY-MM-DD" or null

      champByClass.set(mapKey, {
        name: winner,
        event: eventName,
        date: isoDate,
        isInterim: wc.isInterim,
        display: wc.display,
      })
    }

    // Build output
    const champions = [...champByClass.entries()].map(([mapKey, c]) => ({
      name: c.name,
      weightClassKey: mapKey.split(':')[0],
      weightClassDisplay: c.display,
      isInterim: c.isInterim,
      reignStart: c.date,   // ISO date string or null
      event: c.event,
    }))

    return NextResponse.json({ champions, fetchedAt: new Date().toISOString() })
  } catch (err: any) {
    console.error('[ufc-champions API] Error:', err.message)
    return NextResponse.json({ error: err.message }, { status: 500 })
  }
}
