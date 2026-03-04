import React, { useState, useMemo } from 'react';
import { ScrollView, View, Text, StyleSheet, TouchableOpacity, useWindowDimensions } from 'react-native';
import { LineChart } from 'react-native-gifted-charts';
import SliderRow from '../../components/SliderRow';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

// ─── Modified Raoult's Law + Margules 1-suffix for azeotrope search ──────────
// y_i * P = x_i * gamma_i * Psat_i
// Psat via Antoine: log10(Psat_mmHg) = A - B/(T_C + C)
interface Antoine { A: number; B: number; C: number; Tmin: number; Tmax: number }
const SYSTEMS: { label: string; c1: string; c2: string; a1: Antoine; a2: Antoine; A12: number; A21: number }[] = [
  {
    label: 'Ethanol / Water',
    c1: 'Ethanol', c2: 'Water',
    a1: { A: 8.11220, B: 1592.864, C: 226.184, Tmin: 20, Tmax: 93 },
    a2: { A: 8.07131, B: 1730.63, C: 233.426, Tmin: 1, Tmax: 100 },
    A12: 1.6022, A21: 0.7947,
  },
  {
    label: 'Acetone / Water',
    c1: 'Acetone', c2: 'Water',
    a1: { A: 7.11714, B: 1210.595, C: 229.664, Tmin: -13, Tmax: 55 },
    a2: { A: 8.07131, B: 1730.63, C: 233.426, Tmin: 1, Tmax: 100 },
    A12: 2.0, A21: 1.7,
  },
  {
    label: 'Benzene / Toluene (ideal)',
    c1: 'Benzene', c2: 'Toluene',
    a1: { A: 6.87987, B: 1196.760, C: 219.161, Tmin: 8, Tmax: 80 },
    a2: { A: 6.95087, B: 1342.310, C: 219.187, Tmin: 6, Tmax: 137 },
    A12: 0, A21: 0,
  },
  {
    label: 'n-Hexane / Ethanol',
    c1: 'n-Hexane', c2: 'Ethanol',
    a1: { A: 6.87601, B: 1171.17, C: 224.408, Tmin: -25, Tmax: 92 },
    a2: { A: 8.11220, B: 1592.864, C: 226.184, Tmin: 20, Tmax: 93 },
    A12: 1.9, A21: 1.9,
  },
];

function psat(ant: Antoine, T_C: number): number {
  return Math.pow(10, ant.A - ant.B / (T_C + ant.C)); // mmHg
}

function margules(x1: number, A12: number, A21: number): [number, number] {
  const x2 = 1 - x1;
  const gamma1 = Math.exp(x2 * x2 * (A12 + 2 * (A21 - A12) * x1));
  const gamma2 = Math.exp(x1 * x1 * (A21 + 2 * (A12 - A21) * x2));
  return [gamma1, gamma2];
}

// Bubble temperature at given x1, P_mmHg via bisection
function bubbleT(x1: number, P: number, a1: Antoine, a2: Antoine, A12: number, A21: number): number {
  let Tlo = 20, Thi = 120;
  for (let i = 0; i < 50; i++) {
    const T = (Tlo + Thi) / 2;
    const [g1, g2] = margules(x1, A12, A21);
    const sumY = x1 * g1 * psat(a1, T) + (1 - x1) * g2 * psat(a2, T);
    if (sumY > P) Thi = T; else Tlo = T;
  }
  return (Tlo + Thi) / 2;
}

// Equilibrium y at given x1, T
function yEq(x1: number, T: number, P: number, a1: Antoine, a2: Antoine, A12: number, A21: number): number {
  const [g1] = margules(x1, A12, A21);
  return x1 * g1 * psat(a1, T) / P;
}

const N = 50;

export default function AzeotropeFinderScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [sysIdx, setSysIdx] = useState(0);
  const [P_mmHg, setP] = useState(760);
  const CHART_W = width - 64;

  const sys = SYSTEMS[sysIdx];

  const { xyData, diagData, azeo } = useMemo(() => {
    const { a1, a2, A12, A21 } = sys;
    const xyPts: { value: number; label?: string }[] = [];
    const diagPts: { value: number; label?: string }[] = [];
    let azeoX: number | null = null, azeoY: number | null = null;

    for (let i = 0; i <= N; i++) {
      const x1 = i / N;
      const T = bubbleT(x1, P_mmHg, a1, a2, A12, A21);
      const y1 = yEq(x1, T, P_mmHg, a1, a2, A12, A21);
      const lbl = i % 10 === 0 ? x1.toFixed(1) : '';
      xyPts.push({ value: parseFloat(Math.min(1, Math.max(0, y1)).toFixed(4)), label: lbl });
      diagPts.push({ value: x1, label: lbl });

      // Detect azeotrope: y crosses x
      if (i > 0) {
        const prevX = (i - 1) / N;
        const prevY = xyPts[i - 1].value;
        const curY = xyPts[i].value;
        // sign change in (y-x)
        if ((prevY - prevX) * (curY - x1) < 0 && prevX > 0.01 && x1 < 0.99) {
          azeoX = prevX + (prevX - prevX) / 2; // midpoint approx
          azeoX = (prevX + x1) / 2;
          azeoY = (prevY + curY) / 2;
        }
      }
    }

    return { xyData: xyPts, diagData: diagPts, azeo: azeoX !== null ? { x: azeoX, y: azeoY } : null };
  }, [sys, P_mmHg]);

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      <View style={styles.sysRow}>
        {SYSTEMS.map((s, i) => (
          <TouchableOpacity key={s.label} onPress={() => setSysIdx(i)}
            style={[styles.sysChip, { backgroundColor: sysIdx === i ? colors.blue : colors.card, borderColor: sysIdx === i ? colors.blue : colors.border }]}>
            <Text style={[styles.sysText, { color: sysIdx === i ? '#fff' : colors.foreground }]}>{s.label}</Text>
          </TouchableOpacity>
        ))}
      </View>

      <SectionCard title="Conditions">
        <SliderRow label="Pressure (mmHg)" value={P_mmHg} min={100} max={3000} step={50} decimals={0} onChange={setP} />
      </SectionCard>

      <SectionCard title={`x-y Diagram: ${sys.c1} / ${sys.c2}`}>
        <LineChart
          data={xyData}
          width={CHART_W}
          height={220}
          color={colors.blue}
          thickness={2.5}
          secondaryData={diagData}
          secondaryLineConfig={{ color: colors.border, thickness: 1 }}
          xAxisColor={colors.border}
          yAxisColor={colors.border}
          xAxisLabelTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          yAxisTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          dataPointsRadius={0}
          hideDataPoints
          initialSpacing={0}
          endSpacing={4}
          maxValue={1}
          noOfSections={4}
        />
        <Text style={[styles.axisLabel, { color: colors.mutedForeground }]}>x (liquid), y (vapor) — {sys.c1}</Text>
      </SectionCard>

      <SectionCard title="Azeotrope">
        {azeo ? (
          <View style={styles.azeoBox}>
            <Text style={[styles.azeoText, { color: colors.foreground }]}>
              Azeotrope detected at x ≈ {azeo.x?.toFixed(3)}, y ≈ {azeo.y?.toFixed(3)}
            </Text>
            {azeo.x !== null && (
              <Text style={[styles.azeoPressure, { color: colors.blue }]}>
                P = {P_mmHg} mmHg
              </Text>
            )}
          </View>
        ) : (
          <Text style={[styles.noAzeo, { color: colors.mutedForeground }]}>
            No azeotrope detected for this system at {P_mmHg} mmHg.
          </Text>
        )}
        <Text style={[styles.modelNote, { color: colors.mutedForeground }]}>
          Model: Modified Raoult's Law + {sys.A12 === 0 ? 'ideal (Raoult)' : 'Margules (1-suffix)'}
        </Text>
      </SectionCard>
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  sysRow: { flexDirection: 'row', flexWrap: 'wrap', gap: 8, marginBottom: 14 },
  sysChip: { paddingHorizontal: 12, paddingVertical: 8, borderRadius: 20, borderWidth: 1 },
  sysText: { fontSize: 12, fontWeight: '600' },
  axisLabel: { fontSize: 11, textAlign: 'center', marginTop: 4 },
  azeoBox: { padding: 4 },
  azeoText: { fontSize: 15, fontWeight: '700', marginBottom: 4 },
  azeoPressure: { fontSize: 13 },
  noAzeo: { fontSize: 13, fontStyle: 'italic' },
  modelNote: { fontSize: 11, marginTop: 8 },
});
