import React, { useState, useMemo } from 'react';
import { ScrollView, View, Text, StyleSheet, TouchableOpacity, useWindowDimensions } from 'react-native';
import { LineChart } from 'react-native-gifted-charts';
import SliderRow from '../../components/SliderRow';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

// ─── Reaction types ───────────────────────────────────────────────────────────
type RxnType = 'A→B' | '2A→B' | 'A→B→C' | 'A+B→C';
type RxnOrder = 1 | 2;

interface RxnConfig {
  label: RxnType;
  order: RxnOrder;
  desc: string;
  rA: (Fa: number, Fb: number, Ft: number, k: number, Cao: number) => number; // -rA in mol/(L·s)
}

const REACTIONS: RxnConfig[] = [
  { label: 'A→B', order: 1, desc: '1st order: -rA = k·CA', rA: (Fa, _Fb, Ft, k, Cao) => k * (Fa / Ft) * Cao },
  { label: '2A→B', order: 2, desc: '2nd order: -rA = k·CA²', rA: (Fa, _Fb, Ft, k, Cao) => k * Math.pow((Fa / Ft) * Cao, 2) },
  { label: 'A→B→C', order: 1, desc: 'Series: -rA = k₁·CA', rA: (Fa, _Fb, Ft, k, Cao) => k * (Fa / Ft) * Cao },
  { label: 'A+B→C', order: 2, desc: '2nd order: -rA = k·CA·CB', rA: (Fa, Fb, Ft, k, Cao) => k * (Fa / Ft) * Cao * (Fb / Ft) * Cao },
];

// RK4 for PFR
function rk4PFR(
  V0: number, Vf: number, Fa0: number, Fb0: number,
  rA: (Fa: number, Fb: number, Ft: number) => number,
  stoich_b: number, // stoichiometric coefficient of B consumed per A consumed (+1 for A+B→C, 0 for A→B)
): { V: number[]; X: number[]; Ca: number[] } {
  const N = 120;
  const dV = (Vf - V0) / N;
  const Ft0 = Fa0 + Fb0;
  let Fa = Fa0, Fb = Fb0;
  const Vs: number[] = [0];
  const Xs: number[] = [0];
  const Cas: number[] = [1]; // normalized

  for (let i = 0; i < N; i++) {
    const V = i * dV;
    const Ft = Fa + Fb + (Fa0 - Fa); // approximate total
    const step = (fa: number, fb: number) => {
      const ft = fa + Math.max(1e-9, fb) + Math.max(0, Fa0 - fa);
      return rA(fa, fb, ft);
    };

    const k1a = -step(Fa, Fb) * dV;
    const k1b = -stoich_b * step(Fa, Fb) * dV;
    const k2a = -step(Fa + k1a / 2, Fb + k1b / 2) * dV;
    const k2b = -stoich_b * step(Fa + k1a / 2, Fb + k1b / 2) * dV;
    const k3a = -step(Fa + k2a / 2, Fb + k2b / 2) * dV;
    const k3b = -stoich_b * step(Fa + k2a / 2, Fb + k2b / 2) * dV;
    const k4a = -step(Fa + k3a, Fb + k3b) * dV;
    const k4b = -stoich_b * step(Fa + k3a, Fb + k3b) * dV;

    Fa = Math.max(1e-9, Fa + (k1a + 2 * k2a + 2 * k3a + k4a) / 6);
    Fb = Math.max(1e-9, Fb + (k1b + 2 * k2b + 2 * k3b + k4b) / 6);

    const X = (Fa0 - Fa) / Fa0;
    Vs.push((i + 1) * dV);
    Xs.push(Math.min(1, Math.max(0, X)));
    Cas.push(Math.max(0, Fa / Fa0));
  }
  return { V: Vs, X: Xs, Ca: Cas };
}

// CSTR conversion: X = f(X) via Newton-Raphson (solve tau * (-rA) = Fa0*X/Fa0 = CAo*X)
function cstrConversion(
  tau: number,
  rA: (Fa: number, Fb: number, Ft: number) => number,
  Fa0: number, Fb0: number, Cao: number,
): number {
  let X = 0.5;
  for (let i = 0; i < 50; i++) {
    const Fa = Fa0 * (1 - X);
    const Fb = Math.max(1e-9, Fb0 - Fa0 * X);
    const Ft = Fa + Fb;
    const ra = rA(Fa, Fb, Ft);
    const F = tau * ra - Cao * X;
    const dFa = -Fa0;
    const dFb = -Fa0;
    const eps = 0.001;
    const dra = (rA(Fa + eps, Fb, Ft) - ra) / eps * dFa;
    const dF = tau * dra - Cao;
    const dX = -F / (dF || 1e-12);
    X = Math.min(0.999, Math.max(0.001, X + dX));
  }
  return X;
}

const COLORS = ['#60a5fa', '#f97316', '#22c55e'];

export default function ReactorDesignScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [rxnIdx, setRxnIdx] = useState(0);
  const [k, setK] = useState(0.5);
  const [Cao, setCao] = useState(2);
  const [Cbo, setCbo] = useState(1);
  const [Vmax, setVmax] = useState(10);
  const [reactorType, setReactorType] = useState<'PFR' | 'CSTR'>('PFR');
  const CHART_W = width - 64;

  const rxn = REACTIONS[rxnIdx];

  const { xData, leachingData, cstrX, pfr_X_at_Vmax } = useMemo(() => {
    const Fa0 = Cao; // flow rate = Cao (assume v0=1 L/s)
    const Fb0 = rxn.label === 'A+B→C' ? Cbo : 0;
    const stoich_b = rxn.label === 'A+B→C' ? 1 : 0;
    const ra = (Fa: number, Fb: number, Ft: number) =>
      rxn.rA(Fa, Fb, Ft, k, Cao);
    const { V, X } = rk4PFR(0, Vmax, Fa0, Fb0, ra, stoich_b);

    const stride = Math.max(1, Math.floor(V.length / 5));
    const xData = X.map((x, i) => ({
      value: parseFloat((x * 100).toFixed(2)),
      label: i % stride === 0 ? V[i].toFixed(1) : '',
    }));

    // CSTR: conversion at each volume (tau = V/v0, v0=1)
    const leachingData = V.map((v, i) => {
      const tau = v;
      const X_cstr = tau > 0 ? cstrConversion(tau, ra, Fa0, Fb0, Cao) : 0;
      return {
        value: parseFloat((X_cstr * 100).toFixed(2)),
        label: i % stride === 0 ? v.toFixed(1) : '',
      };
    });

    const cstrTauFinal = Vmax;
    const cstrX = cstrConversion(cstrTauFinal, ra, Fa0, Fb0, Cao);
    const pfr_X_at_Vmax = X[X.length - 1];
    return { xData, leachingData, cstrX, pfr_X_at_Vmax };
  }, [rxn, k, Cao, Cbo, Vmax, reactorType]);

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      <View style={styles.chipRow}>
        {REACTIONS.map((r, i) => (
          <TouchableOpacity key={r.label} onPress={() => setRxnIdx(i)}
            style={[styles.chip, { backgroundColor: rxnIdx === i ? colors.blue : colors.card, borderColor: rxnIdx === i ? colors.blue : colors.border }]}>
            <Text style={[styles.chipText, { color: rxnIdx === i ? '#fff' : colors.foreground }]}>{r.label}</Text>
          </TouchableOpacity>
        ))}
      </View>
      <Text style={[styles.desc, { color: colors.mutedForeground }]}>{rxn.desc}</Text>

      <SectionCard title="Parameters">
        <SliderRow label="Rate Constant k" value={k} min={0.01} max={5} step={0.01} decimals={2} onChange={setK} />
        <SliderRow label="[A]₀ (mol/L)" value={Cao} min={0.1} max={10} step={0.1} decimals={1} onChange={setCao} />
        {rxn.label === 'A+B→C' && (
          <SliderRow label="[B]₀ (mol/L)" value={Cbo} min={0.1} max={10} step={0.1} decimals={1} onChange={setCbo} />
        )}
        <SliderRow label="Max Volume (L)" value={Vmax} min={0.1} max={50} step={0.5} decimals={1} onChange={setVmax} />
      </SectionCard>

      <SectionCard title="Conversion vs Volume">
        <LineChart
          data={xData}
          width={CHART_W}
          height={220}
          color={colors.blue}
          thickness={2.5}
          secondaryData={leachingData}
          secondaryLineConfig={{ color: '#f97316', thickness: 2 }}
          xAxisColor={colors.border}
          yAxisColor={colors.border}
          xAxisLabelTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          yAxisTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          dataPointsRadius={0}
          hideDataPoints
          initialSpacing={0}
          endSpacing={4}
          noOfSections={4}
          maxValue={100}
        />
        <View style={styles.legend}>
          <View style={[styles.dot, { backgroundColor: colors.blue }]} />
          <Text style={[styles.legendText, { color: colors.mutedForeground }]}>PFR</Text>
          <View style={[styles.dot, { backgroundColor: '#f97316' }]} />
          <Text style={[styles.legendText, { color: colors.mutedForeground }]}>CSTR</Text>
        </View>
        <Text style={[styles.xLabel, { color: colors.mutedForeground }]}>Volume (L)</Text>
      </SectionCard>

      <SectionCard title={`Results at V = ${Vmax} L`}>
        <View style={styles.resultRow}>
          <View style={[styles.resultItem, { backgroundColor: colors.card, borderColor: colors.border }]}>
            <Text style={[styles.resultLabel, { color: colors.mutedForeground }]}>PFR Conversion</Text>
            <Text style={[styles.resultValue, { color: colors.blue }]}>{(pfr_X_at_Vmax * 100).toFixed(1)}%</Text>
          </View>
          <View style={[styles.resultItem, { backgroundColor: colors.card, borderColor: colors.border }]}>
            <Text style={[styles.resultLabel, { color: colors.mutedForeground }]}>CSTR Conversion</Text>
            <Text style={[styles.resultValue, { color: '#f97316' }]}>{(cstrX * 100).toFixed(1)}%</Text>
          </View>
        </View>
        <Text style={[styles.note, { color: colors.mutedForeground }]}>
          PFR outperforms CSTR for positive-order reactions (higher conversion at same volume).
        </Text>
      </SectionCard>
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  chipRow: { flexDirection: 'row', flexWrap: 'wrap', gap: 8, marginBottom: 8 },
  chip: { paddingHorizontal: 12, paddingVertical: 8, borderRadius: 20, borderWidth: 1 },
  chipText: { fontSize: 13, fontWeight: '600' },
  desc: { fontSize: 12, marginBottom: 14 },
  legend: { flexDirection: 'row', alignItems: 'center', gap: 6, marginTop: 8, justifyContent: 'center' },
  dot: { width: 10, height: 10, borderRadius: 5 },
  legendText: { fontSize: 12 },
  xLabel: { fontSize: 11, textAlign: 'center', marginTop: 4 },
  resultRow: { flexDirection: 'row', gap: 10, marginBottom: 8 },
  resultItem: { flex: 1, borderRadius: 10, borderWidth: 1, padding: 12, alignItems: 'center' },
  resultLabel: { fontSize: 12, marginBottom: 4 },
  resultValue: { fontSize: 24, fontWeight: '700' },
  note: { fontSize: 12, fontStyle: 'italic' },
});
