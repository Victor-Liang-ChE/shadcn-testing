import { createClient } from '@supabase/supabase-js'

// The real credentials live server-side only (SUPABASE_URL / SUPABASE_ANON_KEY).
// Client-side code never calls Supabase directly — all queries go through
// /api/supabase-proxy.  We only need a dummy client here so TypeScript
// has the correct shape for the Proxy wrapper below.
const realSupabase = createClient(
  'https://placeholder.supabase.co',
  'placeholder-anon-key'
)

/**
 * PROXY WRAPPER: This mimics the Supabase client API but redirects 
 * calls through our server-side API proxy to hide the anon key 
 * from the Network tab.
 */
export const supabase = new Proxy(realSupabase, {
  get(target, prop) {
    if (prop === 'from') {
      return (table: string) => {
        const queryParams: any = { table, filters: [] };
        
        const builder = {
          select: (columns: string = '*') => {
            queryParams.select = columns;
            return builder;
          },
          ilike: (column: string, pattern: string) => {
            queryParams.filters.push({ type: 'ilike', column, value: pattern });
            return builder;
          },
          eq: (column: string, value: any) => {
            queryParams.filters.push({ type: 'eq', column, value });
            return builder;
          },
          in: (column: string, values: any[]) => {
            queryParams.filters.push({ type: 'in', column, value: values });
            return builder;
          },
          or: (condition: string) => {
            queryParams.filters.push({ type: 'or', value: condition });
            return builder;
          },
          order: (column: string, { ascending = true } = {}) => {
            queryParams.order = { column, ascending };
            return builder;
          },
          limit: (n: number) => {
            queryParams.limit = n;
            return builder;
          },
          // The terminal "then" allows this to be awaited like a real query
          then: async (onfulfilled: any, onrejected: any) => {
            try {
              const res = await fetch('/api/supabase-proxy', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(queryParams)
              });
              const json = await res.json();
              // Mirror the Supabase return object { data, error }
              const result = {
                  data: json.data || null,
                  error: json.error ? { message: json.error } : null,
                  status: res.status,
                  statusText: res.statusText
              };
              return onfulfilled ? onfulfilled(result) : result;
            } catch (err: any) {
              const errorResult = { data: null, error: { message: err.message }, status: 500 };
              return onrejected ? onrejected(errorResult) : Promise.resolve(errorResult);
            }
          }
        };
        return builder;
      };
    }
    return (target as any)[prop];
  }
}) as unknown as typeof realSupabase;

/**
 * Legacy proxy function (optional)
 */
export async function supabaseQueryProxy(params: {
  table: string,
  select?: string,
  filter?: any,
  limit?: number,
  ilike?: any,
  or?: string,
  order?: { column: string, ascending?: boolean }
}) {
  const response = await fetch('/api/supabase-proxy', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params)
  })
  if (!response.ok) {
    const errorJson = await response.json()
    throw new Error(errorJson.error || 'Proxy request failed')
  }
  return response.json()
}

// Define types based on your Supabase table schemas
export interface FighterDetails {
  id?: number;
  FIGHTER?: string | null; // Name of the fighter
  HEIGHT?: string | null;
  WEIGHT?: string | null;
  REACH?: string | null;
  STANCE?: string | null;
  DOB?: string | null; // Date string 'YYYY-MM-DD'
  URL?: string | null;
  created_at?: string;
  // Remove or comment out old fields if they don't exist in the new table
  // NICKNAME?: string | null;
  // FIRST?: string | null; // No longer used for fetching
  // LAST?: string | null;  // No longer used for fetching
}

export interface FightResult {
  id: number;
  BOUT: string;        // e.g., "Jon Jones vs. Ciryl Gane"
  OUTCOME: string;     // e.g., "win", "loss", "draw", "nc"
  WEIGHTCLASS: string | null;
  METHOD: string | null;
  ROUND: number | null;
  TIME: string | null;
  EVENT?: string | null;
  DETAILS?: string | null;
  URL?: string | null;
  // These might need to be derived or added to your table schema
  OPPONENT?: string;
  DATE?: string;
  FIGHTER_A?: string; // Consider adding columns like these to your DB
  FIGHTER_B?: string; // Consider adding columns like these to your DB
  WINNER?: string;    // Consider adding columns like these to your DB
}

// Add interface for the new images table
export interface UfcChampionImages {
    id: number;
    FIRST?: string | null; // Assuming name columns for matching
    LAST?: string | null;  // Assuming name columns for matching
    MAINSHOT?: string | null;
    HEADSHOT?: string | null;
}

