import React, { useState, useMemo } from 'react';
import { ScrollView, View, Text, StyleSheet, TouchableOpacity, useWindowDimensions } from 'react-native';
import { LineChart } from 'react-native-gifted-charts';
import SliderRow from '../../components/SliderRow';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

type Order = '1st' | '2nd';
type InputType = 'Step' | 'Ramp';

const N_PTS = 200;

function erf(x: number): number {
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
  const r = 1 - poly * Math.exp(-x * x);
  return x >= 0 ? r : -r;
}

export default function ProcessDynamicsScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [order, setOrder] = useState<Order>('1st');
  const [inputType, setInputType] = useState<InputType>('Step');
  const [K, setK] = useState(2);
  const [tau, setTau] = useState(5);
  const [M, setM] = useState(1);
  const [zeta, setZeta] = useState(0.5);
  const [tMax, setTMax] = useState(40);

  const CHART_W = width - 64;

  const { responseData, inputData, metrics } = useMemo(() => {
    const dt = tMax / N_PTS;
    const response: { value: number; label?: string }[] = [];
    const input: { value: number; label?: string }[] = [];
    const metrics: string[] = [];

    for (let i = 0; i <= N_PTS; i++) {
      const t = i * dt;
      const lbl = i % (N_PTS / 4) === 0 ? t.toFixed(0) : '';
      let y = 0, u = M;

      if (order === '1st') {
        if (inputType === 'Step') {
          y = K * M * (1 - Math.exp(-t / tau));
        } else {
          y = K * M * (t - tau * (1 - Math.exp(-t / tau)));
        }
      } else {
        if (zeta < 1 - 1e-6) {
          const wd = Math.sqrt(1 - zeta * zeta) / tau;
          const a = zeta / Math.sqrt(1 - zeta * zeta);
          if (inputType === 'Step') {
            y = K * M * (1 - Math.exp(-zeta * t / tau) * (Math.cos(wd * t) + a * Math.sin(wd * t)));
          } else {
            y = K * M * (t - 2 * zeta / (1 / tau) * (1 - Math.exp(-zeta * t / tau) * (Math.cos(wd * t) + a * Math.sin(wd * t))));
          }
        } else if (Math.abs(zeta - 1) < 1e-6) {
          y = K * M * (1 - (1 + t / tau) * Math.exp(-t / tau));
        } else {
          const r1 = (-zeta + Math.sqrt(zeta * zeta - 1)) / tau;
          const r2 = (-zeta - Math.sqrt(zeta * zeta - 1)) / tau;
          const A = K * M * r2 / (r2 - r1), B = -K * M * r1 / (r2 - r1);
          y = K * M + A * Math.exp(r1 * t) + B * Math.exp(r2 * t);
        }
      }

      response.push({ value: isFinite(y) ? y : 0, label: lbl });
      input.push({ value: u, label: lbl });
    }

    if (order === '2nd' && zeta < 1 && inputType === 'Step') {
      const wd = Math.sqrt(1 - zeta * zeta) / tau;
      const tp = Math.PI / wd;
      const OS = Math.exp(-Math.PI * zeta / Math.sqrt(1 - zeta * zeta));
      const ts = 4 / (zeta / tau);
      metrics.push(`Peak Time: ${tp.toFixed(2)} s`);
      metrics.push(`Overshoot: ${(OS * 100).toFixed(1)}%`);
      metrics.push(`Settling Time (4τ): ${ts.toFixed(2)} s`);
    }

    return { responseData: response, inputData: input, metrics };
  }, [order, inputType, K, tau, M, zeta, tMax]);

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      <View style={styles.rowGap}>
        {(['1st', '2nd'] as Order[]).map(o => (
          <TouchableOpacity key={o} onPress={() => setOrder(o)}
            style={[styles.chip, { backgroundColor: order === o ? colors.blue : colors.card, borderColor: order === o ? colors.blue : colors.border }]}>
            <Text style={[styles.chipText, { color: order === o ? '#fff' : colors.foreground }]}>{o} Order</Text>
          </TouchableOpacity>
        ))}
        {(['Step', 'Ramp'] as InputType[]).map(t => (
          <TouchableOpacity key={t} onPress={() => setInputType(t)}
            style={[styles.chip, { backgroundColor: inputType === t ? '#22c55e' : colors.card, borderColor: inputType === t ? '#22c55e' : colors.border }]}>
            <Text style={[styles.chipText, { color: inputType === t ? '#fff' : colors.foreground }]}>{t}</Text>
          </TouchableOpacity>
        ))}
      </View>

      <SectionCard title="Parameters">
        <SliderRow label="Gain K" value={K} min={0.1} max={10} step={0.1} decimals={1} onChange={setK} />
        <SliderRow label="Input Magnitude M" value={M} min={0.1} max={5} step={0.1} decimals={1} onChange={setM} />
        <SliderRow label="Time Constant τ (s)" value={tau} min={0.5} max={30} step={0.5} decimals={1} onChange={setTau} />
        {order === '2nd' && (
          <SliderRow label="Damping Ratio ζ" value={zeta} min={0.05} max={2} step={0.05} decimals={2} onChange={setZeta} />
        )}
        <SliderRow label="Time Horizon (s)" value={tMax} min={5} max={100} step={5} decimals={0} onChange={setTMax} />
      </SectionCard>

      <SectionCard title="System Response">
        <LineChart
          data={responseData}
          width={CHART_W}
          height={220}
          color={colors.blue}
          thickness={2.5}
          secondaryData={inputData}
          secondaryLineConfig={{ color: '#ef4444', thickness: 1.5 }}
          xAxisColor={colors.border}
          yAxisColor={colors.border}
          xAxisLabelTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          yAxisTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          dataPointsRadius={0}
          hideDataPoints
          curved
          initialSpacing={0}
          endSpacing={4}
          noOfSections={4}
        />
        <View style={styles.legendRow}>
          <View style={[styles.legendDot, { backgroundColor: colors.blue }]} />
          <Text style={[styles.legendText, { color: colors.mutedForeground }]}>Process Output</Text>
          <View style={[styles.legendDot, { backgroundColor: '#ef4444' }]} />
          <Text style={[styles.legendText, { color: colors.mutedForeground }]}>Input</Text>
        </View>
      </SectionCard>

      {metrics.length > 0 && (
        <SectionCard title="Performance Metrics">
          {metrics.map((m, i) => (
            <Text key={i} style={[styles.metric, { color: colors.foreground }]}>• {m}</Text>
          ))}
        </SectionCard>
      )}
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  rowGap: { flexDirection: 'row', flexWrap: 'wrap', gap: 8, marginBottom: 14 },
  chip: { paddingHorizontal: 14, paddingVertical: 8, borderRadius: 20, borderWidth: 1 },
  chipText: { fontSize: 13, fontWeight: '600' },
  legendRow: { flexDirection: 'row', alignItems: 'center', gap: 6, marginTop: 8, justifyContent: 'center' },
  legendDot: { width: 10, height: 10, borderRadius: 5 },
  legendText: { fontSize: 12 },
  metric: { fontSize: 14, marginBottom: 4 },
});
