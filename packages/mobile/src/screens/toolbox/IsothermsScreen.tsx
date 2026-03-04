import React, { useState, useMemo } from 'react';
import { ScrollView, View, Text, StyleSheet, TouchableOpacity, useWindowDimensions } from 'react-native';
import { LineChart } from 'react-native-gifted-charts';
import SliderRow from '../../components/SliderRow';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

type IsothermModel = 'Langmuir' | 'Freundlich' | 'Temkin';
type XVariable = 'Concentration' | 'Pressure';

export default function IsothermsScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [model, setModel] = useState<IsothermModel>('Langmuir');
  const [xVar, setXVar] = useState<XVariable>('Concentration');
  // Langmuir
  const [qmax, setQmax] = useState(50);
  const [KL, setKL] = useState(1);
  // Freundlich
  const [KF, setKF] = useState(5);
  const [n, setN] = useState(2);
  // Temkin
  const [A, setA] = useState(2);
  const [b, setB] = useState(10);

  const CHART_W = width - 64;
  const xLabel = xVar === 'Concentration' ? 'C (mmol/L)' : 'P (atm)';
  const xMax = xVar === 'Concentration' ? 10 : 100;
  const N_PTS = 80;
  const R = 8.314;
  const T = 298;

  const chartData = useMemo(() => {
    return Array.from({ length: N_PTS }, (_, i) => {
      const x = (i + 1) * xMax / N_PTS;
      let q = 0;
      if (model === 'Langmuir') {
        q = qmax * KL * x / (1 + KL * x);
      } else if (model === 'Freundlich') {
        q = KF * Math.pow(x, 1 / n);
      } else {
        q = (R * T / b) * Math.log(Math.max(A * x, 1e-10));
      }
      return { value: Math.max(q, 0), label: i % 20 === 0 ? x.toFixed(0) : '' };
    });
  }, [model, xVar, qmax, KL, KF, n, A, b]);

  const MODELS: IsothermModel[] = ['Langmuir', 'Freundlich', 'Temkin'];
  const XVARS: XVariable[] = ['Concentration', 'Pressure'];

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      <View style={styles.rowGap}>
        {MODELS.map(m => (
          <TouchableOpacity key={m} onPress={() => setModel(m)}
            style={[styles.chip, { backgroundColor: model === m ? colors.blue : colors.card, borderColor: model === m ? colors.blue : colors.border }]}>
            <Text style={[styles.chipText, { color: model === m ? '#fff' : colors.foreground }]}>{m}</Text>
          </TouchableOpacity>
        ))}
      </View>

      <View style={[styles.rowGap, { marginBottom: 14 }]}>
        {XVARS.map(v => (
          <TouchableOpacity key={v} onPress={() => setXVar(v)}
            style={[styles.chip, { backgroundColor: xVar === v ? '#22c55e' : colors.card, borderColor: xVar === v ? '#22c55e' : colors.border }]}>
            <Text style={[styles.chipText, { color: xVar === v ? '#fff' : colors.foreground }]}>{v}</Text>
          </TouchableOpacity>
        ))}
      </View>

      <SectionCard title="Parameters">
        {model === 'Langmuir' && <>
          <SliderRow label="q_max" value={qmax} min={1} max={200} step={1} unit=" mmol/g" decimals={0} onChange={setQmax} />
          <SliderRow label="K_L" value={KL} min={0.01} max={20} step={0.01} unit=" L/mmol" decimals={2} onChange={setKL} />
        </>}
        {model === 'Freundlich' && <>
          <SliderRow label="K_F" value={KF} min={0.1} max={50} step={0.1} unit="" decimals={2} onChange={setKF} />
          <SliderRow label="n" value={n} min={0.5} max={10} step={0.1} unit="" decimals={1} onChange={setN} />
        </>}
        {model === 'Temkin' && <>
          <SliderRow label="A" value={A} min={0.1} max={20} step={0.1} unit="" decimals={1} onChange={setA} />
          <SliderRow label="b (kJ/mol)" value={b} min={1} max={100} step={1} unit="" decimals={0} onChange={setB} />
        </>}
      </SectionCard>

      <SectionCard title={`q vs ${xLabel}`}>
        <LineChart
          data={chartData}
          width={CHART_W}
          height={220}
          color={colors.blue}
          thickness={2.5}
          xAxisColor={colors.border}
          yAxisColor={colors.border}
          xAxisLabelTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          yAxisTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          dataPointsColor={colors.blue}
          dataPointsRadius={0}
          hideDataPoints
          curved
          initialSpacing={4}
          endSpacing={4}
          noOfSections={4}
        />
        <Text style={[styles.axisLabel, { color: colors.mutedForeground }]}>q (mmol/g) vs {xLabel}</Text>
      </SectionCard>

      <SectionCard title="Model Equations">
        {model === 'Langmuir' && (
          <Text style={[styles.eq, { color: colors.foreground }]}>
            q = q_max · K_L · C / (1 + K_L · C){'\n'}
            Monolayer adsorption with fixed sites.
          </Text>
        )}
        {model === 'Freundlich' && (
          <Text style={[styles.eq, { color: colors.foreground }]}>
            q = K_F · C^(1/n){'\n'}
            Empirical model for heterogeneous surfaces.
          </Text>
        )}
        {model === 'Temkin' && (
          <Text style={[styles.eq, { color: colors.foreground }]}>
            q = (RT/b) · ln(A · C){'\n'}
            Assumes linear decrease in adsorption heat.
          </Text>
        )}
      </SectionCard>
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  rowGap: { flexDirection: 'row', gap: 8, marginBottom: 8 },
  chip: { paddingHorizontal: 14, paddingVertical: 8, borderRadius: 20, borderWidth: 1 },
  chipText: { fontSize: 13, fontWeight: '600' },
  axisLabel: { fontSize: 11, textAlign: 'center', marginTop: 6 },
  eq: { fontSize: 13, lineHeight: 20, fontFamily: 'monospace' },
});
