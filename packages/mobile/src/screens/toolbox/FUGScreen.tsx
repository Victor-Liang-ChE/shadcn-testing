import React, { useState, useMemo } from 'react';
import { ScrollView, View, Text, StyleSheet, TextInput, TouchableOpacity, useWindowDimensions } from 'react-native';
import { BarChart } from 'react-native-gifted-charts';
import SectionCard from '../../components/SectionCard';
import SliderRow from '../../components/SliderRow';
import { useTheme } from '../../theme/ThemeContext';

// ─── FUG: Fenske-Underwood-Gilliland short-cut distillation design ───────────
// Uses relative volatility. Components: sorted heavy → light.

type Component = { name: string; z: number; alpha: number };

const DEFAULT_COMPS: Component[] = [
  { name: 'Propane', z: 0.25, alpha: 4.5 },
  { name: 'i-Butane (LK)', z: 0.35, alpha: 2.1 },
  { name: 'n-Butane (HK)', z: 0.25, alpha: 1.0 },
  { name: 'i-Pentane', z: 0.15, alpha: 0.45 },
];

// Fenske minimum stages
function fenske(LK_xD: number, LK_xB: number, HK_xD: number, HK_xB: number, alphaLH: number): number {
  return Math.log(LK_xD / LK_xB * HK_xB / HK_xD) / Math.log(alphaLH);
}

// Underwood theta (between alpha_HK=1 and alpha_LK)
function underwoodTheta(
  components: Component[], q: number
): number {
  // Σ [alpha_i * z_i / (alpha_i - theta)] = 1 - q
  // theta must lie between alpha_HK=1 and alpha_LK
  const target = 1 - q;
  let lo = 1.0001, hi = components.find(c => c.name.includes('LK') || c.alpha > 1.5)?.alpha ?? 2;
  hi = hi - 0.0001;
  for (let i = 0; i < 80; i++) {
    const theta = (lo + hi) / 2;
    const sum = components.reduce((s, c) => s + c.alpha * c.z / (c.alpha - theta), 0);
    if (sum < target) lo = theta; else hi = theta;
  }
  return (lo + hi) / 2;
}

// Min reflux ratio from Underwood
function underwood_Rmin(
  components: Component[], distillate_fracs: number[], theta: number
): number {
  const num = components.reduce((s, c, i) => s + c.alpha * distillate_fracs[i] / (c.alpha - theta), 0);
  return num - 1;
}

// Gilliland correlation (N vs R/Rmin)
function gilliland(Rmin: number, R: number, Nmin: number): number {
  const x = (R - Rmin) / (R + 1);
  const Y = 1 - Math.exp((1 + 54.4 * x) / (11 + 117.2 * x) * (x - 1) / Math.pow(x, 0.5));
  return (Nmin + Y) / (1 - Y);
}

const COLORS = ['#60a5fa', '#22c55e', '#f97316', '#a855f7', '#ec4899'];

export default function FUGScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [components, setComponents] = useState<Component[]>(DEFAULT_COMPS);
  const [LK_rec, setLK_rec] = useState(0.99); // LK recovery in distillate
  const [HK_rec, setHK_rec] = useState(0.01); // HK recovery in distillate
  const [q, setQ] = useState(1.0); // feed quality
  const [R_ratio, setR_ratio] = useState(1.5); // R / Rmin
  const CHART_W = width - 80;

  const lkIdx = components.findIndex(c => c.name.includes('LK') || (c.alpha > 1 && c.alpha < 3));
  const hkIdx = components.findIndex(c => c.name.includes('HK') || Math.abs(c.alpha - 1) < 0.1);
  const LK = components[lkIdx] ?? components[1];
  const HK = components[hkIdx] ?? components[2];

  const results = useMemo(() => {
    const totalZ = components.reduce((s, c) => s + c.z, 0) || 1;
    const norm = components.map(c => ({ ...c, z: c.z / totalZ }));
    const lk = norm.find(c => c.name === LK.name)!;
    const hk = norm.find(c => c.name === HK.name)!;

    // Assume total distillate D = sum of recovered fractions
    const D = lk.z * LK_rec + hk.z * HK_rec +
      norm.filter((_, i) => norm[i].alpha > lk.alpha).reduce((s, c) => s + c.z, 0);
    const B = 1 - D;
    const dFracs = norm.map(c => {
      if (c.name === lk.name) return lk.z * LK_rec / D;
      if (c.name === hk.name) return hk.z * HK_rec / D;
      return c.alpha > lk.alpha ? c.z / D : 0.001 * c.z / (B || 0.001);
    });
    const bFracs = norm.map((c, i) => (c.z - dFracs[i] * D) / (B || 0.001));

    const Nmin = fenske(dFracs[lkIdx], bFracs[lkIdx], dFracs[hkIdx], bFracs[hkIdx], lk.alpha / hk.alpha);
    const theta = underwoodTheta(norm, q);
    const Rmin = underwood_Rmin(norm, dFracs, theta);
    const R = Rmin * R_ratio;
    const N_actual = gilliland(Rmin, R, Nmin);
    const NF = 0.5 * N_actual; // approximate feed stage

    return { Nmin, Rmin, R, N_actual: Math.ceil(N_actual), NF: Math.ceil(NF), dFracs, bFracs, norm, D, B };
  }, [components, LK_rec, HK_rec, q, R_ratio]);

  const barData = results.dFracs.map((d, i) => ({
    value: parseFloat((d * 100).toFixed(1)),
    label: components[i]?.name.replace(' (LK)', '').replace(' (HK)', '').substring(0, 6),
    frontColor: COLORS[i % COLORS.length],
    topLabelComponent: () => (
      <Text style={{ color: colors.foreground, fontSize: 10, marginBottom: 2 }}>{(d * 100).toFixed(0)}%</Text>
    ),
  }));

  const botBarData = results.bFracs.map((b, i) => ({
    value: parseFloat((Math.max(0, b) * 100).toFixed(1)),
    label: components[i]?.name.replace(' (LK)', '').replace(' (HK)', '').substring(0, 6),
    frontColor: COLORS[i % COLORS.length],
    topLabelComponent: () => (
      <Text style={{ color: colors.foreground, fontSize: 10, marginBottom: 2 }}>{(Math.max(0, b) * 100).toFixed(0)}%</Text>
    ),
  }));

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      <SectionCard title="Feed Composition & Volatility">
        {components.map((c, i) => (
          <View key={c.name} style={[styles.compRow, { borderColor: colors.border }]}>
            <Text style={[styles.compName, { color: colors.foreground }]}>{c.name}</Text>
            <View style={styles.compSliders}>
              <SliderRow label={`z (feed)`} value={c.z} min={0.01} max={1} step={0.01} decimals={2}
                onChange={v => setComponents(p => p.map((x, j) => j === i ? { ...x, z: v } : x))} />
              <SliderRow label={`α (rel. vol.)`} value={c.alpha} min={0.1} max={20} step={0.05} decimals={2}
                onChange={v => setComponents(p => p.map((x, j) => j === i ? { ...x, alpha: v } : x))} />
            </View>
          </View>
        ))}
      </SectionCard>

      <SectionCard title="Design Specifications">
        <SliderRow label="LK Recovery in Distillate" value={LK_rec} min={0.8} max={0.9999} step={0.001} decimals={3} onChange={setLK_rec} />
        <SliderRow label="HK Recovery in Distillate" value={HK_rec} min={0.0001} max={0.2} step={0.001} decimals={3} onChange={setHK_rec} />
        <SliderRow label="q (feed quality)" value={q} min={0} max={2} step={0.05} decimals={2} onChange={setQ} />
        <SliderRow label="R / Rmin" value={R_ratio} min={1.05} max={5} step={0.05} decimals={2} onChange={setR_ratio} />
      </SectionCard>

      <SectionCard title="Results">
        <View style={styles.resultGrid}>
          {[
            { label: 'Nmin', value: results.Nmin.toFixed(1) },
            { label: 'Rmin', value: results.Rmin.toFixed(2) },
            { label: 'N actual', value: results.N_actual.toString() },
            { label: 'Feed Stage', value: results.NF.toString() },
            { label: 'R', value: results.R.toFixed(2) },
            { label: 'D/F', value: results.D.toFixed(3) },
          ].map(r => (
            <View key={r.label} style={[styles.resultItem, { backgroundColor: colors.card, borderColor: colors.border }]}>
              <Text style={[styles.resultLabel, { color: colors.mutedForeground }]}>{r.label}</Text>
              <Text style={[styles.resultValue, { color: colors.blue }]}>{r.value}</Text>
            </View>
          ))}
        </View>
      </SectionCard>

      <SectionCard title="Distillate Composition">
        <BarChart
          data={barData}
          width={CHART_W}
          height={160}
          barWidth={32}
          spacing={10}
          xAxisColor={colors.border}
          yAxisColor={colors.border}
          xAxisLabelTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          yAxisTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          noOfSections={4}
          maxValue={100}
          initialSpacing={8}
        />
        <Text style={[styles.chartNote, { color: colors.mutedForeground }]}>Mole % in distillate</Text>
      </SectionCard>

      <SectionCard title="Bottoms Composition">
        <BarChart
          data={botBarData}
          width={CHART_W}
          height={160}
          barWidth={32}
          spacing={10}
          xAxisColor={colors.border}
          yAxisColor={colors.border}
          xAxisLabelTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          yAxisTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          noOfSections={4}
          maxValue={100}
          initialSpacing={8}
        />
        <Text style={[styles.chartNote, { color: colors.mutedForeground }]}>Mole % in bottoms</Text>
      </SectionCard>
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  compRow: { borderBottomWidth: StyleSheet.hairlineWidth, paddingBottom: 8, marginBottom: 8 },
  compName: { fontSize: 14, fontWeight: '700', marginBottom: 4 },
  compSliders: {},
  resultGrid: { flexDirection: 'row', flexWrap: 'wrap', gap: 8 },
  resultItem: { borderWidth: 1, borderRadius: 8, padding: 10, minWidth: 80, alignItems: 'center' },
  resultLabel: { fontSize: 11, marginBottom: 2 },
  resultValue: { fontSize: 18, fontWeight: '700' },
  chartNote: { fontSize: 11, textAlign: 'center', marginTop: 4 },
});
