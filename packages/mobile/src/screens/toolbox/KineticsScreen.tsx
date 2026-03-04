import React, { useState, useMemo } from 'react';
import { ScrollView, View, Text, StyleSheet, TouchableOpacity, TextInput, useWindowDimensions } from 'react-native';
import { LineChart } from 'react-native-gifted-charts';
import SliderRow from '../../components/SliderRow';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

// ── RK4 ODE solver ──────────────────────────────────────────────────────────
function rk4(
  f: (t: number, y: number[]) => number[],
  y0: number[],
  tSpan: [number, number],
  nSteps: number,
): { t: number[]; y: number[][] } {
  const [t0, tf] = tSpan;
  const h = (tf - t0) / nSteps;
  const t: number[] = [t0];
  const y: number[][] = [y0.slice()];
  let yn = y0.slice();
  let tn = t0;
  for (let i = 0; i < nSteps; i++) {
    const k1 = f(tn, yn);
    const k2 = f(tn + h / 2, yn.map((v, j) => v + (h / 2) * k1[j]));
    const k3 = f(tn + h / 2, yn.map((v, j) => v + (h / 2) * k2[j]));
    const k4 = f(tn + h, yn.map((v, j) => v + h * k3[j]));
    yn = yn.map((v, j) => v + (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
    tn += h;
    t.push(tn);
    y.push(yn.slice());
  }
  return { t, y };
}

const COLORS = ['#60a5fa', '#f97316', '#22c55e', '#a855f7', '#ec4899', '#facc15'];

type ReactionType = 'A→B' | 'A→B→C' | 'A+B→C' | '2A→B';

interface ReactionConfig {
  label: ReactionType;
  species: string[];
  ode: (y: number[], k: number, k2: number) => number[];
  desc: string;
}

const REACTIONS: ReactionConfig[] = [
  {
    label: 'A→B',
    species: ['A', 'B'],
    ode: ([a, b], k) => [-k * a, k * a],
    desc: 'r = k·[A]  (1st order)',
  },
  {
    label: 'A→B→C',
    species: ['A', 'B', 'C'],
    ode: ([a, b, c], k, k2) => [-k * a, k * a - k2 * b, k2 * b],
    desc: 'r₁ = k₁·[A], r₂ = k₂·[B]',
  },
  {
    label: 'A+B→C',
    species: ['A', 'B', 'C'],
    ode: ([a, b, c], k) => [-k * a * b, -k * a * b, k * a * b],
    desc: 'r = k·[A]·[B]  (2nd order)',
  },
  {
    label: '2A→B',
    species: ['A', 'B'],
    ode: ([a, b], k) => [-2 * k * a * a, k * a * a],
    desc: 'r = k·[A]²  (2nd order)',
  },
];

export default function KineticsScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [rxnIdx, setRxnIdx] = useState(0);
  const [k1, setK1] = useState(0.3);
  const [k2, setK2] = useState(0.1);
  const [tMax, setTMax] = useState(20);
  const [C0, setC0] = useState<number[]>([2, 1, 0]);

  const rxn = REACTIONS[rxnIdx];
  const CHART_W = width - 64;
  const N_PTS = 150;

  const datasets = useMemo(() => {
    const y0 = rxn.species.map((_, i) => C0[i] ?? 0);
    const { t, y } = rk4((_, yn) => rxn.ode(yn, k1, k2), y0, [0, tMax], N_PTS);
    const step = Math.max(1, Math.floor(t.length / 5));
    return rxn.species.map((sp, si) =>
      y.map((row, i) => ({
        value: Math.max(0, row[si]),
        label: i % step === 0 ? t[i].toFixed(0) : '',
      })),
    );
  }, [rxn, k1, k2, tMax, C0]);

  function setC0i(i: number, val: number) {
    setC0(prev => {
      const next = [...prev];
      next[i] = val;
      return next;
    });
  }

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      {/* Reaction selector */}
      <View style={styles.chipRow}>
        {REACTIONS.map((r, i) => (
          <TouchableOpacity key={r.label} onPress={() => setRxnIdx(i)}
            style={[styles.chip, { backgroundColor: rxnIdx === i ? colors.blue : colors.card, borderColor: rxnIdx === i ? colors.blue : colors.border }]}>
            <Text style={[styles.chipText, { color: rxnIdx === i ? '#fff' : colors.foreground }]}>{r.label}</Text>
          </TouchableOpacity>
        ))}
      </View>
      <Text style={[styles.desc, { color: colors.mutedForeground }]}>{rxn.desc}</Text>

      <SectionCard title="Initial Concentrations (mol/L)">
        {rxn.species.map((sp, i) => (
          <SliderRow key={sp} label={`[${sp}]₀`} value={C0[i] ?? 0} min={0} max={5} step={0.1} decimals={1} onChange={v => setC0i(i, v)} />
        ))}
      </SectionCard>

      <SectionCard title="Rate Constants">
        <SliderRow label="k₁ (1/s or L/mol·s)" value={k1} min={0.01} max={2} step={0.01} decimals={2} onChange={setK1} />
        {rxnIdx === 1 && (
          <SliderRow label="k₂ (1/s)" value={k2} min={0.01} max={2} step={0.01} decimals={2} onChange={setK2} />
        )}
        <SliderRow label="Time Horizon (s)" value={tMax} min={1} max={100} step={1} decimals={0} onChange={setTMax} />
      </SectionCard>

      <SectionCard title="Concentration vs Time">
        <LineChart
          data={datasets[0]}
          width={CHART_W}
          height={220}
          color={COLORS[0]}
          thickness={2.5}
          dataSet={datasets.slice(1).map((d, i) => ({
            data: d,
            color: COLORS[i + 1],
            thickness: 2.5,
          }))}
          xAxisColor={colors.border}
          yAxisColor={colors.border}
          xAxisLabelTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          yAxisTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          dataPointsRadius={0}
          hideDataPoints
          initialSpacing={0}
          endSpacing={4}
          noOfSections={4}
          curved
        />
        <View style={styles.legend}>
          {rxn.species.map((sp, i) => (
            <View key={sp} style={styles.legendItem}>
              <View style={[styles.dot, { backgroundColor: COLORS[i] }]} />
              <Text style={[styles.legendText, { color: colors.mutedForeground }]}>[{sp}]</Text>
            </View>
          ))}
        </View>
      </SectionCard>
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  chipRow: { flexDirection: 'row', flexWrap: 'wrap', gap: 8, marginBottom: 8 },
  chip: { paddingHorizontal: 14, paddingVertical: 8, borderRadius: 20, borderWidth: 1 },
  chipText: { fontSize: 13, fontWeight: '600' },
  desc: { fontSize: 12, marginBottom: 14 },
  legend: { flexDirection: 'row', flexWrap: 'wrap', gap: 12, marginTop: 8, justifyContent: 'center' },
  legendItem: { flexDirection: 'row', alignItems: 'center', gap: 5 },
  dot: { width: 10, height: 10, borderRadius: 5 },
  legendText: { fontSize: 12 },
});
