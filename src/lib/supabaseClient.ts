import { createClient } from '@supabase/supabase-js'

// Read environment variables
const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL
const supabaseAnonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY

// Check if variables are loaded (optional but good practice)
if (!supabaseUrl || !supabaseAnonKey) {
  throw new Error("Supabase URL or Anon Key is missing from environment variables.");
}

// Initialize and export the Supabase client
export const supabase = createClient(supabaseUrl, supabaseAnonKey)

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

