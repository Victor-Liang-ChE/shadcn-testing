export interface Champion {
  name: string;
  reignStart: string;
  reignEnd: string | null;
  notes?: string;
  imageUrl?: string;
}
export interface WeightClassData {
  displayName: string;
  champions: Champion[];
}

// Historical lineage is now built dynamically from Supabase via
// /api/ufc-champions/lineage – this constant is kept only as a fallback
// shape reference and is intentionally empty.
export const ufcChampionsData: Record<string, WeightClassData> = {};
