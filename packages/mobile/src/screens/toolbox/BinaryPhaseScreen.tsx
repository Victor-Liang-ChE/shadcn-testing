import React, { useState, useMemo } from 'react';
import { ScrollView, View, Text, StyleSheet, TouchableOpacity, useWindowDimensions } from 'react-native';
import { LineChart } from 'react-native-gifted-charts';
import SliderRow from '../../components/SliderRow';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

// ─── Same Antoine + Margules as AzeotropeFinder ──────────────────────────────
interface Antoine { A: number; B: number; C: number }
const SYSTEMS = [
  {
    label: 'Ethanol / Water', c1: 'Ethanol', c2: 'Water',
    a1: { A: 8.11220, B: 1592.864, C: 226.184 },
    a2: { A: 8.07131, B: 1730.63, C: 233.426 },
    A12: 1.6022, A21: 0.7947,
  },
  {
    label: 'Benzene / Toluene', c1: 'Benzene', c2: 'Toluene',
    a1: { A: 6.87987, B: 1196.760, C: 219.161 },
    a2: { A: 6.95087, B: 1342.310, C: 219.187 },
    A12: 0, A21: 0,
  },
  {
    label: 'Acetone / Water', c1: 'Acetone', c2: 'Water',
    a1: { A: 7.11714, B: 1210.595, C: 229.664 },
    a2: { A: 8.07131, B: 1730.63, C: 233.426 },
    A12: 2.0, A21: 1.7,
  },
  {
    label: 'n-Hexane / Toluene', c1: 'n-Hexane', c2: 'Toluene',
    a1: { A: 6.87601, B: 1171.17, C: 224.408 },
    a2: { A: 6.95087, B: 1342.310, C: 219.187 },
    A12: 0.3, A21: 0.3,
  },
];

function psat(a: Antoine, T_C: number) {
  return Math.pow(10, a.A - a.B / (T_C + a.C));
}
function margules(x1: number, A12: number, A21: number): [number, number] {
  const x2 = 1 - x1;
  return [
    Math.exp(x2 * x2 * (A12 + 2 * (A21 - A12) * x1)),
    Math.exp(x1 * x1 * (A21 + 2 * (A12 - A21) * x2)),
  ];
}

type DiagType = 'Txy' | 'Pxy';
const N = 60;

// Bubble T at given x1, P (mmHg)
function bubbleT(x1: number, P: number, a1: Antoine, a2: Antoine, A12: number, A21: number): number {
  let lo = 0, hi = 200;
  for (let i = 0; i < 60; i++) {
    const T = (lo + hi) / 2;
    const [g1, g2] = margules(x1, A12, A21);
    const sumK = x1 * g1 * psat(a1, T) + (1 - x1) * g2 * psat(a2, T);
    if (sumK > P) hi = T; else lo = T;
  }
  return (lo + hi) / 2;
}
// Dew T at given y1, P (mmHg)
function dewT(y1: number, P: number, a1: Antoine, a2: Antoine, A12: number, A21: number): number {
  let lo = 0, hi = 200;
  for (let i = 0; i < 60; i++) {
    const T = (lo + hi) / 2;
    // Initial x1 guess from ideal
    let x1 = y1 * P / psat(a1, T);
    x1 = Math.min(0.999, Math.max(0.001, x1));
    const [g1, g2] = margules(x1, A12, A21);
    const sumK = y1 / (g1 * psat(a1, T) / P) + (1 - y1) / (g2 * psat(a2, T) / P);
    if (sumK < 1) hi = T; else lo = T;
  }
  return (lo + hi) / 2;
}
// Bubble P at given x1, T (°C)
function bubbleP(x1: number, T: number, a1: Antoine, a2: Antoine, A12: number, A21: number): number {
  const [g1, g2] = margules(x1, A12, A21);
  return x1 * g1 * psat(a1, T) + (1 - x1) * g2 * psat(a2, T);
}
// Dew P at given y1, T (°C)
function dewP(y1: number, T: number, a1: Antoine, a2: Antoine, A12: number, A21: number): number {
  let lo = 0, hi = 5000;
  for (let i = 0; i < 60; i++) {
    const P = (lo + hi) / 2;
    let x1 = y1 * P / psat(a1, T);
    x1 = Math.min(0.999, Math.max(0.001, x1));
    const [g1, g2] = margules(x1, A12, A21);
    const sumInv = y1 / (g1 * psat(a1, T) / P) + (1 - y1) / (g2 * psat(a2, T) / P);
    if (sumInv > 1) lo = P; else hi = P;
  }
  return (lo + hi) / 2;
}

export default function BinaryPhaseScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [sysIdx, setSysIdx] = useState(0);
  const [diagType, setDiagType] = useState<DiagType>('Txy');
  const [P_mmHg, setP] = useState(760);
  const [T_C, setT_C] = useState(70);
  const CHART_W = width - 64;
  const sys = SYSTEMS[sysIdx];

  const { bubbleData, dewData, xlabel } = useMemo(() => {
    const { a1, a2, A12, A21 } = sys;
    const bubble: { value: number; label?: string }[] = [];
    const dew: { value: number; label?: string }[] = [];
    const stride = Math.floor(N / 5);

    for (let i = 0; i <= N; i++) {
      const z = i / N;
      const lbl = i % stride === 0 ? z.toFixed(1) : '';
      if (diagType === 'Txy') {
        const Tb = bubbleT(z, P_mmHg, a1, a2, A12, A21);
        const Td = dewT(z, P_mmHg, a1, a2, A12, A21);
        bubble.push({ value: parseFloat(Tb.toFixed(2)), label: lbl });
        dew.push({ value: parseFloat(Td.toFixed(2)), label: lbl });
      } else {
        const Pb = bubbleP(z, T_C, a1, a2, A12, A21);
        const Pd = dewP(z, T_C, a1, a2, A12, A21);
        bubble.push({ value: parseFloat(Pb.toFixed(2)), label: lbl });
        dew.push({ value: parseFloat(Pd.toFixed(2)), label: lbl });
      }
    }
    return { bubbleData: bubble, dewData: dew, xlabel: sys.c1 + ' mole fraction' };
  }, [sys, diagType, P_mmHg, T_C]);

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      <View style={styles.chipRow}>
        {SYSTEMS.map((s, i) => (
          <TouchableOpacity key={s.label} onPress={() => setSysIdx(i)}
            style={[styles.chip, { backgroundColor: sysIdx === i ? colors.blue : colors.card, borderColor: sysIdx === i ? colors.blue : colors.border }]}>
            <Text style={[styles.chipText, { color: sysIdx === i ? '#fff' : colors.foreground }]}>{s.label}</Text>
          </TouchableOpacity>
        ))}
      </View>

      <View style={styles.typeRow}>
        {(['Txy', 'Pxy'] as DiagType[]).map(t => (
          <TouchableOpacity key={t} onPress={() => setDiagType(t)}
            style={[styles.typeChip, { backgroundColor: diagType === t ? '#22c55e' : colors.card, borderColor: diagType === t ? '#22c55e' : colors.border }]}>
            <Text style={[styles.chipText, { color: diagType === t ? '#fff' : colors.foreground }]}>{t} Diagram</Text>
          </TouchableOpacity>
        ))}
      </View>

      <SectionCard title="Conditions">
        {diagType === 'Txy' ? (
          <SliderRow label="Pressure (mmHg)" value={P_mmHg} min={100} max={5000} step={50} decimals={0} onChange={setP} />
        ) : (
          <SliderRow label="Temperature (°C)" value={T_C} min={0} max={200} step={1} decimals={0} onChange={setT_C} />
        )}
      </SectionCard>

      <SectionCard title={diagType === 'Txy' ? `T-x-y at ${P_mmHg} mmHg` : `P-x-y at ${T_C} °C`}>
        <LineChart
          data={bubbleData}
          width={CHART_W}
          height={220}
          color={colors.blue}
          thickness={2.5}
          secondaryData={dewData}
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
        />
        <View style={styles.legend}>
          <View style={[styles.legendDot, { backgroundColor: colors.blue }]} />
          <Text style={[styles.legendText, { color: colors.mutedForeground }]}>Bubble curve</Text>
          <View style={[styles.legendDot, { backgroundColor: '#f97316' }]} />
          <Text style={[styles.legendText, { color: colors.mutedForeground }]}>Dew curve</Text>
        </View>
        <Text style={[styles.xLabel, { color: colors.mutedForeground }]}>{xlabel}</Text>
        <Text style={[styles.modelNote, { color: colors.mutedForeground }]}>
          Modified Raoult's Law + {sys.A12 === 0 ? 'Ideal (Raoult)' : 'Margules'}
        </Text>
      </SectionCard>
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  chipRow: { flexDirection: 'row', flexWrap: 'wrap', gap: 8, marginBottom: 8 },
  typeRow: { flexDirection: 'row', gap: 8, marginBottom: 14 },
  chip: { paddingHorizontal: 12, paddingVertical: 8, borderRadius: 20, borderWidth: 1 },
  typeChip: { paddingHorizontal: 14, paddingVertical: 8, borderRadius: 20, borderWidth: 1 },
  chipText: { fontSize: 12, fontWeight: '600' },
  legend: { flexDirection: 'row', alignItems: 'center', gap: 6, marginTop: 8, justifyContent: 'center' },
  legendDot: { width: 10, height: 10, borderRadius: 5 },
  legendText: { fontSize: 12 },
  xLabel: { fontSize: 11, textAlign: 'center', marginTop: 4 },
  modelNote: { fontSize: 10, textAlign: 'center', marginTop: 2, fontStyle: 'italic' },
});
