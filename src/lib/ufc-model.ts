/**
 * UFC Betting Model — Pure TypeScript tree evaluator
 * Evaluates CatBoost oblivious trees + LightGBM decision trees
 * No Python dependencies — runs anywhere JS/TS runs (Vercel, Edge, etc.)
 */

// ─── Types ───

interface CatBoostSplit {
  feature_idx: number
  border: number
}

interface CatBoostTree {
  splits: CatBoostSplit[]
  leaf_values: number[]
}

interface CatBoostModel {
  type: 'catboost'
  features: string[]
  scale_and_bias: [number | number[], number[]]
  trees: CatBoostTree[]
}

interface LightGBMLeaf {
  leaf_value: number
}

interface LightGBMSplit {
  feature: number
  threshold: number
  left: LightGBMNode
  right: LightGBMNode
  decision_type: string
}

type LightGBMNode = LightGBMLeaf | LightGBMSplit

interface LightGBMModel {
  type: 'lightgbm'
  features: string[]
  trees: LightGBMNode[]
}

// ─── Sigmoid ───

function sigmoid(x: number): number {
  if (x >= 0) {
    return 1 / (1 + Math.exp(-x))
  }
  const e = Math.exp(x)
  return e / (1 + e)
}

// ─── CatBoost Prediction ───

function predictCatBoostRaw(model: CatBoostModel, features: number[]): number {
  let total = 0
  for (const tree of model.trees) {
    let idx = 0
    for (let i = 0; i < tree.splits.length; i++) {
      const split = tree.splits[i]
      if (features[split.feature_idx] > split.border) {
        idx |= 1 << i
      }
    }
    total += tree.leaf_values[idx]
  }
  const scale =
    typeof model.scale_and_bias[0] === 'number'
      ? model.scale_and_bias[0]
      : model.scale_and_bias[0][0]
  const bias = model.scale_and_bias[1][0]
  return total * scale + bias
}

export function predictCatBoost(model: CatBoostModel, features: number[]): number {
  return sigmoid(predictCatBoostRaw(model, features))
}

// ─── LightGBM Prediction ───

function traverseLightGBM(node: LightGBMNode, features: number[]): number {
  if ('leaf_value' in node) {
    return node.leaf_value
  }
  const split = node as LightGBMSplit
  if (features[split.feature] <= split.threshold) {
    return traverseLightGBM(split.left, features)
  }
  return traverseLightGBM(split.right, features)
}

export function predictLightGBM(model: LightGBMModel, features: number[]): number {
  let total = 0
  for (const tree of model.trees) {
    total += traverseLightGBM(tree, features)
  }
  return sigmoid(total)
}

// ─── Odds Conversion ───

/** American odds → implied probability (with vig) */
export function americanToImplied(odds: number): number {
  if (odds === 0 || isNaN(odds)) return 0.5
  return odds < 0 ? Math.abs(odds) / (Math.abs(odds) + 100) : 100 / (odds + 100)
}

/** American odds → decimal odds */
export function americanToDecimal(odds: number): number {
  if (odds > 0) return 1 + odds / 100
  if (odds < 0) return 1 + 100 / Math.abs(odds)
  return 2.0
}

/** Decimal odds → American odds */
export function decimalToAmerican(dec: number): number {
  if (dec >= 2) return Math.round((dec - 1) * 100)
  if (dec > 1) return Math.round(-100 / (dec - 1))
  return 0
}

/** Remove vig: normalize implied probs to sum to 100% */
export function normalizeImplied(rImp: number, bImp: number): [number, number] {
  const total = rImp + bImp
  if (total === 0) return [0.5, 0.5]
  return [rImp / total, bImp / total]
}

// ─── Fighter Name Normalization ───

export function normalizeName(name: string): string {
  return name
    .toLowerCase()
    .replace(/\./g, '')
    .replace(/'/g, '')
    .replace(/-/g, ' ')
    .trim()
}

// ─── Fight Prediction ───

export interface FightInput {
  rFighter: string
  bFighter: string
  rOdds: number | null
  bOdds: number | null
  rPerBook: Record<string, number>
  bPerBook: Record<string, number>
  rEvolution?: OddsEvolution
  bEvolution?: OddsEvolution
  matchupId?: number
}

export interface OddsEvolution {
  mean: { timestamp: number; decimal: number; american: number }[]
  per_book: Record<string, { timestamp: number; decimal: number; american: number }[]>
}

export interface FighterStats {
  name: string
  win_streak: number
  wins: number
  losses: number
  weightClass?: string
}

/**
 * 'advised'     — historically +EV zone (dog odds ≤ +130, well-sampled)
 * 'caution'     — model finds edge but historically weak in this range (+130 to +175)
 * 'limited'     — insufficient historical sample to characterize (+175+)
 */
export type BetAdvice = 'advised' | 'caution' | 'limited'

export interface FightPrediction {
  rFighter: string
  bFighter: string
  rOdds: number | null
  bOdds: number | null
  rImplied: number
  bImplied: number
  rImpliedNorm: number
  bImpliedNorm: number
  rStreak: number
  bStreak: number
  cbProbRed: number
  lgProbRed: number
  avgProbRed: number
  modelPick: string
  confidence: number
  betSignal: 'bet' | 'none'
  betFighter: string
  betOdds: number
  /** Model's probability (0-1) that the bet fighter wins. 0 when no bet. */
  betProbModel: number
  /** Break-even probability at sportsbook odds (1/decimal). 0 when no bet. */
  breakEvenSB: number
  /** Break-even probability at fair/exchange odds (= de-vigged implied). 0 when no bet. */
  breakEvenEX: number
  /** Estimated $ EV per $100 flat bet at a sportsbook (with vig). 0 when no bet. */
  evSB: number
  /** Estimated $ EV per $100 flat bet at an exchange / no-vig market. 0 when no bet. */
  evEX: number
  /** Sportsbook-specific advice based on historical odds-range performance. */
  betAdviceSB: BetAdvice
  /** Exchange/prediction-market advice (same zones, slightly better EV). */
  betAdviceEX: BetAdvice
  rPerBook: Record<string, number>
  bPerBook: Record<string, number>
  rEvolution?: OddsEvolution
  bEvolution?: OddsEvolution
  rDebut: boolean
  bDebut: boolean
  rRecord: string
  bRecord: string
  weightClass?: string
  modelProb: number
  matchupId?: number
}

export interface EventPrediction {
  name: string
  url?: string
  fights: FightPrediction[]
  activeBets: number
}

export function predictFight(
  cbModel: CatBoostModel,
  lgModel: LightGBMModel,
  fight: FightInput,
  fighterDb: Record<string, FighterStats>
): FightPrediction {
  // Look up fighter stats
  const rInfo = lookupFighter(fighterDb, fight.rFighter)
  const bInfo = lookupFighter(fighterDb, fight.bFighter)

  const rStreak = rInfo?.win_streak ?? 0
  const bStreak = bInfo?.win_streak ?? 0

  // Compute implied probs
  const rImpRaw = fight.rOdds != null ? americanToImplied(fight.rOdds) : 0.5
  const bImpRaw = fight.bOdds != null ? americanToImplied(fight.bOdds) : 0.5
  const [rImpNorm, bImpNorm] = normalizeImplied(rImpRaw, bImpRaw)

  // Features: [ODDS_DIFF_NORM, WIN_STREAK_DIFF]
  // Use de-vigged (normalized) odds diff — the same feature the model was trained on
  const oddsDiff = rImpNorm - bImpNorm
  const streakDiff = rStreak - bStreak
  const features = [oddsDiff, streakDiff]

  // Model predictions
  const cbProbR = predictCatBoost(cbModel, features)
  const lgProbR = predictLightGBM(lgModel, features)
  const avgProbR = (cbProbR + lgProbR) / 2

  // Market favorite
  const favIsRed = rImpRaw >= bImpRaw

  // Both models disagree with market → bet signal
  const cbDisagrees = favIsRed !== (cbProbR >= 0.5)
  const lgDisagrees = favIsRed !== (lgProbR >= 0.5)
  const dualDisagree = cbDisagrees && lgDisagrees

  // NO weight class skipping — all fights get bet signals
  const betSignal: 'bet' | 'none' = dualDisagree ? 'bet' : 'none'

  const modelPickRed = avgProbR >= 0.5
  const confidence = Math.abs(avgProbR - 0.5) * 200

  const rRecord = rInfo ? `${rInfo.wins}-${rInfo.losses}` : '—'
  const bRecord = bInfo ? `${bInfo.wins}-${bInfo.losses}` : '—'
  const weightClass = rInfo?.weightClass ?? bInfo?.weightClass
  const modelProb = Math.round((modelPickRed ? avgProbR : 1 - avgProbR) * 100)

  // ─── Bet advice fields ───
  // betProbModel: model's probability of the underdog (bet fighter) winning
  const betProbModel = dualDisagree ? (favIsRed ? 1 - avgProbR : avgProbR) : 0
  // Fair implied probability for the bet fighter (used as break-even on exchange)
  const betFairImpNorm = dualDisagree ? (favIsRed ? bImpNorm : rImpNorm) : 0
  const betDecimalSB = dualDisagree ? americanToDecimal(favIsRed ? (fight.bOdds ?? 100) : (fight.rOdds ?? 100)) : 2
  const betDecimalEX = betFairImpNorm > 0 ? 1 / betFairImpNorm : 2

  const breakEvenSB = dualDisagree ? Math.round((1 / betDecimalSB) * 1000) / 10 : 0
  const breakEvenEX = dualDisagree ? Math.round(betFairImpNorm * 1000) / 10 : 0
  // EV per $100 flat bet: betProbModel * profit_if_win - (1-betProbModel) * 100
  const evSB = dualDisagree ? Math.round((betProbModel * (betDecimalSB - 1) - (1 - betProbModel)) * 100) : 0
  const evEX = dualDisagree ? Math.round((betProbModel * (betDecimalEX - 1) - (1 - betProbModel)) * 100) : 0

  // Advice tiers based on walk-forward performance (v160 analysis)
  // +EV range: -110 to +130. Caution: +130 to +175. Limited sample: +175+
  const betOddsVal = favIsRed ? (fight.bOdds ?? 0) : (fight.rOdds ?? 0)
  const betAdviceSB: BetAdvice = !dualDisagree ? 'limited'
    : betOddsVal <= 130 ? 'advised'
    : betOddsVal <= 175 ? 'caution'
    : 'limited'
  const betAdviceEX: BetAdvice = betAdviceSB  // same zones; EX is slightly better EV in each

  return {
    rFighter: fight.rFighter,
    bFighter: fight.bFighter,
    rOdds: fight.rOdds,
    bOdds: fight.bOdds,
    rImplied: Math.round(rImpRaw * 1000) / 10,
    bImplied: Math.round(bImpRaw * 1000) / 10,
    rImpliedNorm: Math.round(rImpNorm * 1000) / 10,
    bImpliedNorm: Math.round(bImpNorm * 1000) / 10,
    rStreak,
    bStreak,
    cbProbRed: Math.round(cbProbR * 10000) / 10000,
    lgProbRed: Math.round(lgProbR * 10000) / 10000,
    avgProbRed: Math.round(avgProbR * 10000) / 10000,
    modelPick: modelPickRed ? fight.rFighter : fight.bFighter,
    confidence: Math.round(confidence * 10) / 10,
    betSignal,
    betFighter: dualDisagree ? (favIsRed ? fight.bFighter : fight.rFighter) : '',
    betOdds: dualDisagree ? (favIsRed ? (fight.bOdds ?? 0) : (fight.rOdds ?? 0)) : 0,
    betProbModel: Math.round(betProbModel * 10000) / 10000,
    breakEvenSB,
    breakEvenEX,
    evSB,
    evEX,
    betAdviceSB,
    betAdviceEX,
    rPerBook: fight.rPerBook,
    bPerBook: fight.bPerBook,
    rEvolution: fight.rEvolution,
    bEvolution: fight.bEvolution,
    rDebut: rInfo ? (rInfo.wins === 0 && rInfo.losses === 0) : false,
    bDebut: bInfo ? (bInfo.wins === 0 && bInfo.losses === 0) : false,
    rRecord,
    bRecord,
    weightClass,
    modelProb,
    matchupId: fight.matchupId,
  }
}

/** Look up fighter by name: exact → last-name → partial */
function lookupFighter(
  db: Record<string, FighterStats>,
  name: string
): FighterStats | null {
  const key = normalizeName(name)
  if (db[key]) return db[key]

  // Try last name match
  const parts = key.split(/\s+/)
  if (parts.length > 0) {
    const last = parts[parts.length - 1]
    const matches = Object.entries(db).filter(([k]) => k.endsWith(last))
    if (matches.length === 1) return matches[0][1]
  }

  return null
}
