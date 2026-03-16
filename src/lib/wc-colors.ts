/**
 * Weight-class color palette.
 *
 * Gradient: Crystal Teal (#E0F7FA) → Vibrant Teal (#009688) → Soft Dark Teal (#00695C)
 * Ordered from lightest weight class to heaviest.
 */

// Ordered lightest → heaviest (men's divisions only; women map via WC_ALIAS)
const WC_ORDER = [
  'women_strawweight',   // 115 lbs (no men's equivalent; unique lightest colour)
  'flyweight',           // 125 lbs
  'bantamweight',        // 135 lbs
  'featherweight',       // 145 lbs
  'lightweight',         // 155 lbs
  'welterweight',        // 170 lbs
  'middleweight',        // 185 lbs
  'lightheavyweight',    // 205 lbs
  'heavyweight',         // 265 lbs
] as const;

// Women's divisions share the same colour as their men's weight equivalent.
// women_strawweight has no men's equivalent so it keeps its own entry in WC_ORDER above.
const WC_ALIAS: Record<string, string> = {
  women_flyweight:      'flyweight',
  women_bantamweight:   'bantamweight',
  women_featherweight:  'featherweight',
};

// Stop colors along the gradient: Crystal Teal → Vibrant Teal → Soft Dark Teal
const STOPS: [number, [number, number, number]][] = [
  [0,   [178, 235, 242]], // #B2EBF2 — Light Teal (Strawweights)
  [0.6, [0,   150, 136]], // #009688 — Vibrant Teal (Mid-weight classes)
  [1,   [0,   105, 92]],  // #00695C — Soft Dark Teal (Heavyweight)
];

function lerp(a: number, b: number, t: number) {
  return Math.round(a + (b - a) * t);
}

function gradientColor(t: number): string {
  // Find which two stops surround t
  let i = 0;
  while (i < STOPS.length - 2 && t > STOPS[i + 1][0]) i++;
  const [t0, c0] = STOPS[i];
  const [t1, c1] = STOPS[i + 1];
  const localT = t1 === t0 ? 0 : (t - t0) / (t1 - t0);
  const r = lerp(c0[0], c1[0], localT);
  const g = lerp(c0[1], c1[1], localT);
  const b = lerp(c0[2], c1[2], localT);
  return `rgb(${r},${g},${b})`;
}

/**
 * Returns `{ bg, text }` inline style values for a weight-class badge/chip.
 * `key` should be a normalized weight-class key (lowercase, no spaces, no "title").
 * `opacity` controls the background alpha (default 0.15 for light bg, text is the solid gradient colour).
 */
export function wcBadgeStyle(
  key: string,
  opts?: { bgOpacity?: number }
): { background: string; color: string; borderColor: string } {
  const cleanKey = key.toLowerCase().replace(/\s+/g, '');
  // Women's divisions share colour with their men's equivalent via alias
  const resolvedKey = WC_ALIAS[cleanKey] ?? cleanKey;
  const idx = (WC_ORDER as readonly string[]).indexOf(resolvedKey);
  const total = WC_ORDER.length - 1;
  const t = idx < 0 ? 0.5 : idx / total; // fallback: mid-point

  const solid = gradientColor(t);

  // Parse rgb values for alpha background
  const m = solid.match(/\d+/g)!;
  const [r, g, b] = m.map(Number);
  
  // Use a slightly higher default opacity for the text to ensure contrast 
  // with the very light sky blue on light backgrounds.
  const bgOpacity = opts?.bgOpacity ?? 0.12;
  const textOpacity = 0.95;

  return {
    background: `rgba(${r},${g},${b},${bgOpacity})`,
    color: `rgba(${r},${g},${b},${textOpacity})`,
    borderColor: `rgba(${r},${g},${b},0.35)`,
  };
}
