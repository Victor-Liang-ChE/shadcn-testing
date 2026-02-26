import { createClient } from '@supabase/supabase-js'
import { NextResponse } from 'next/server'

function createServerSupabaseClient() {
  const supabaseUrl = process.env.SUPABASE_URL
  const supabaseAnonKey = process.env.SUPABASE_ANON_KEY
  if (!supabaseUrl || !supabaseAnonKey) throw new Error('Supabase environment variables are missing.')
  return createClient(supabaseUrl, supabaseAnonKey)
}

export async function POST(request: Request) {
  try {
    const { table, select, filters, limit, order } = await request.json()
    const supabase = createServerSupabaseClient()

    let query: any = supabase.from(table).select(select || '*')

    if (filters && Array.isArray(filters)) {
      for (const f of filters) {
        if (f.type === 'eq') {
          query = query.eq(f.column, f.value)
        } else if (f.type === 'ilike') {
          query = query.ilike(f.column, f.value)
        } else if (f.type === 'in') {
          query = query.in(f.column, f.value)
        } else if (f.type === 'or') {
          query = query.or(f.value)
        }
      }
    }

    if (order) {
      query = query.order(order.column, { ascending: order.ascending })
    }

    if (limit) {
      query = query.limit(limit)
    }

    const { data, error } = await query

    if (error) throw error

    return NextResponse.json({ data })
  } catch (err: any) {
    console.error('Proxy Error:', err.message)
    return NextResponse.json({ error: err.message }, { status: 500 })
  }
}
