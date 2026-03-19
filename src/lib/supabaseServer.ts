import { createClient } from '@supabase/supabase-js'

/**
 * Creates a server-side Supabase client using the private (non-public) env vars.
 * Use this only in server-only code ("use server" files, Route Handlers, etc.).
 * These env vars are never shipped to the browser.
 */
export function createServerSupabaseClient() {
  const supabaseUrl = process.env.SUPABASE_URL
  const supabaseAnonKey = process.env.SUPABASE_ANON_KEY
  const supabaseServiceKey = process.env.SUPABASE_SERVICE_KEY

  if (!supabaseUrl || (!supabaseAnonKey && !supabaseServiceKey)) {
    throw new Error('Supabase environment variables are missing.')
  }

  // Prefer Service Key for server-side operations if available (bypasses RLS/private bucket restrictions)
  return createClient(supabaseUrl, supabaseServiceKey || supabaseAnonKey!)
}
