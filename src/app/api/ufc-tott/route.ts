import { createClient } from '@supabase/supabase-js'
import { NextResponse } from 'next/server'

function createServerSupabaseClient() {
  const supabaseUrl = process.env.SUPABASE_URL
  const supabaseAnonKey = process.env.SUPABASE_ANON_KEY
  if (!supabaseUrl || !supabaseAnonKey) throw new Error('Supabase environment variables are missing.')
  return createClient(supabaseUrl, supabaseAnonKey, {
    auth: { persistSession: false, autoRefreshToken: false, detectSessionInUrl: false },
  })
}

// Fetch all rows from ufc_fighter_tott in pages to bypass PostgREST's default 1000-row cap
export async function GET() {
  try {
    const supabase = createServerSupabaseClient()
    const PAGE = 1000
    const all: any[] = []
    let from = 0

    while (true) {
      const { data, error } = await supabase
        .from('ufc_fighter_tott')
        .select('FIGHTER, HEIGHT, WEIGHT, REACH, STANCE, DOB')
        .range(from, from + PAGE - 1)

      if (error) throw error
      if (!data || data.length === 0) break
      all.push(...data)
      if (data.length < PAGE) break
      from += PAGE
    }

    return NextResponse.json({ data: all })
  } catch (err: any) {
    console.error('UFC TOTT API Error:', err.message)
    return NextResponse.json({ error: err.message }, { status: 500 })
  }
}
