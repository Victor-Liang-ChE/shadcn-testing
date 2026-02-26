import { createClient } from '@supabase/supabase-js'

/**
 * Creates a server-side Supabase client using the private (non-public) env vars.
 * Use this only in server-only code ("use server" files, Route Handlers, etc.).
 * These env vars are never shipped to the browser.
 */
export function createServerSupabaseClient() {
  const supabaseUrl = process.env.SUPABASE_URL
  const supabaseAnonKey = process.env.SUPABASE_ANON_KEY
  if (!supabaseUrl || !supabaseAnonKey) {
    throw new Error('SUPABASE_URL or SUPABASE_ANON_KEY environment variable is missing.')
  }
  return createClient(supabaseUrl, supabaseAnonKey)
}
