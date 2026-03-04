import React, { useState, useMemo, useCallback } from 'react';
import {
  ScrollView, View, Text, StyleSheet, TextInput, TouchableOpacity,
  KeyboardAvoidingView, Platform, useWindowDimensions,
} from 'react-native';
import { LineChart } from 'react-native-gifted-charts';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

type CalcMode = 'compound_interest' | 'EAR' | 'NPV' | 'annuity' | 'perpetuity' | 'inflation' | 'depreciation';

const MODES: { key: CalcMode; label: string }[] = [
  { key: 'compound_interest', label: 'Compound Interest' },
  { key: 'EAR', label: 'EAR' },
  { key: 'NPV', label: 'NPV' },
  { key: 'annuity', label: 'Annuity' },
  { key: 'perpetuity', label: 'Perpetuity' },
  { key: 'inflation', label: 'Inflation' },
  { key: 'depreciation', label: 'Depreciation' },
];

function fmt(n: number): string {
  if (Math.abs(n) >= 1e12) return (n / 1e12).toFixed(2) + 'T';
  if (Math.abs(n) >= 1e9) return (n / 1e9).toFixed(2) + 'B';
  if (Math.abs(n) >= 1e6) return (n / 1e6).toFixed(2) + 'M';
  if (Math.abs(n) >= 1e3) return (n / 1e3).toFixed(2) + 'K';
  return n.toFixed(2);
}

export default function ChemeEconScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [mode, setMode] = useState<CalcMode>('compound_interest');
  const [inputs, setInputs] = useState<Record<string, string>>({
    principal: '10000', rate: '5', periods: '10', n_compounds: '12',
    cashflow: '1000', discount: '8', payment: '500',
    growth: '3', initial_cost: '100000', salvage: '10000', life: '10',
  });

  const set = (k: string, v: string) => setInputs(prev => ({ ...prev, [k]: v }));

  const { result, chartData } = useMemo(() => {
    const p = (k: string) => parseFloat(inputs[k]) || 0;
    let result = '';
    let chartData: { value: number; label?: string }[] = [];

    if (mode === 'compound_interest') {
      const P = p('principal'), r = p('rate') / 100, n = p('n_compounds'), t = p('periods');
      for (let yr = 1; yr <= Math.min(t, 30); yr++) {
        chartData.push({ value: P * Math.pow(1 + r / n, n * yr), label: yr % 5 === 0 ? `${yr}y` : '' });
      }
      const FV = P * Math.pow(1 + r / n, n * t);
      result = `Future Value = $${fmt(FV)}\nInterest Earned = $${fmt(FV - P)}`;
    } else if (mode === 'EAR') {
      const nominal = p('rate') / 100, n = p('n_compounds');
      const ear = Math.pow(1 + nominal / n, n) - 1;
      result = `EAR = ${(ear * 100).toFixed(4)}%\n(Nominal: ${(nominal * 100).toFixed(2)}%, compounded ${n}×/yr)`;
    } else if (mode === 'NPV') {
      const cf = p('cashflow'), r = p('discount') / 100, t = p('periods'), I0 = p('principal');
      let npv = -I0;
      for (let i = 1; i <= Math.min(t, 30); i++) {
        const pv = cf / Math.pow(1 + r, i);
        npv += pv;
        chartData.push({ value: npv, label: i % 5 === 0 ? `${i}y` : '' });
      }
      result = `NPV = $${fmt(npv)}\n${npv > 0 ? '✓ Positive NPV – project profitable' : '✗ Negative NPV – project unprofitable'}`;
    } else if (mode === 'annuity') {
      const pmt = p('payment'), r = p('rate') / 100 / 12, n = p('periods') * 12;
      const PV = pmt * (1 - Math.pow(1 + r, -n)) / r;
      const FV = pmt * (Math.pow(1 + r, n) - 1) / r;
      for (let m = 1; m <= Math.min(n, 120); m += Math.ceil(n / 30)) {
        chartData.push({ value: pmt * (Math.pow(1 + r, m) - 1) / r, label: m % 12 === 0 ? `${m / 12}y` : '' });
      }
      result = `PV of Annuity = $${fmt(PV)}\nFV of Annuity = $${fmt(FV)}`;
    } else if (mode === 'perpetuity') {
      const pmt = p('payment'), r = p('rate') / 100, g = p('growth') / 100;
      const PV = r > g ? pmt / (r - g) : Infinity;
      result = `PV of Perpetuity = $${isFinite(PV) ? fmt(PV) : '∞'}\n(Growing perpetuity, g=${(g * 100).toFixed(1)}%)`;
    } else if (mode === 'inflation') {
      const P = p('principal'), inf = p('rate') / 100, t = p('periods');
      const realValue = P / Math.pow(1 + inf, t);
      for (let yr = 1; yr <= Math.min(t, 30); yr++) {
        chartData.push({ value: P / Math.pow(1 + inf, yr), label: yr % 5 === 0 ? `${yr}y` : '' });
      }
      result = `Real Value = $${fmt(realValue)}\nLoss of Purchasing Power = $${fmt(P - realValue)}`;
    } else if (mode === 'depreciation') {
      const cost = p('initial_cost'), salv = p('salvage'), life = p('life');
      const straightLine = (cost - salv) / life;
      let bookVal = cost;
      for (let yr = 1; yr <= Math.min(life, 30); yr++) {
        bookVal -= straightLine;
        chartData.push({ value: Math.max(bookVal, salv), label: yr % 2 === 0 ? `${yr}y` : '' });
      }
      result = `Straight-Line Depreciation = $${fmt(straightLine)}/yr\nFinal Book Value = $${fmt(Math.max(salv, cost - straightLine * life))}`;
    }

    return { result, chartData };
  }, [mode, inputs]);

  const InputRow = ({ label, field, hint }: { label: string; field: string; hint?: string }) => (
    <View style={styles.inputRow}>
      <Text style={[styles.inputLabel, { color: colors.foreground }]}>{label}</Text>
      <TextInput
        style={[styles.inputField, { color: colors.foreground, borderColor: colors.border, backgroundColor: colors.background }]}
        value={inputs[field]}
        onChangeText={v => set(field, v)}
        keyboardType="decimal-pad"
        placeholder={hint || '0'}
        placeholderTextColor={colors.mutedForeground}
      />
    </View>
  );

  const FIELD_MAP: Record<CalcMode, { label: string; field: string; hint?: string }[]> = {
    compound_interest: [
      { label: 'Principal ($)', field: 'principal' },
      { label: 'Annual Rate (%)', field: 'rate' },
      { label: 'Compounds/Year', field: 'n_compounds' },
      { label: 'Years', field: 'periods' },
    ],
    EAR: [
      { label: 'Nominal Rate (%)', field: 'rate' },
      { label: 'Compounds/Year', field: 'n_compounds' },
    ],
    NPV: [
      { label: 'Initial Investment ($)', field: 'principal' },
      { label: 'Annual Cash Flow ($)', field: 'cashflow' },
      { label: 'Discount Rate (%)', field: 'discount' },
      { label: 'Years', field: 'periods' },
    ],
    annuity: [
      { label: 'Payment/Month ($)', field: 'payment' },
      { label: 'Annual Rate (%)', field: 'rate' },
      { label: 'Years', field: 'periods' },
    ],
    perpetuity: [
      { label: 'Annual Payment ($)', field: 'payment' },
      { label: 'Discount Rate (%)', field: 'rate' },
      { label: 'Growth Rate (%)', field: 'growth' },
    ],
    inflation: [
      { label: 'Current Value ($)', field: 'principal' },
      { label: 'Inflation Rate (%)', field: 'rate' },
      { label: 'Years', field: 'periods' },
    ],
    depreciation: [
      { label: 'Initial Cost ($)', field: 'initial_cost' },
      { label: 'Salvage Value ($)', field: 'salvage' },
      { label: 'Useful Life (yr)', field: 'life' },
    ],
  };

  const CHART_W = width - 64;

  return (
    <KeyboardAvoidingView style={{ flex: 1, backgroundColor: colors.background }}
      behavior={Platform.OS === 'ios' ? 'padding' : undefined}>
      <ScrollView contentContainerStyle={styles.container} keyboardShouldPersistTaps="handled">
        {/* Mode picker */}
        <ScrollView horizontal showsHorizontalScrollIndicator={false} style={{ marginBottom: 14 }}>
          {MODES.map(m => (
            <TouchableOpacity key={m.key} onPress={() => setMode(m.key)}
              style={[styles.modeChip, {
                backgroundColor: mode === m.key ? colors.blue : colors.card,
                borderColor: mode === m.key ? colors.blue : colors.border,
              }]}>
              <Text style={[styles.modeText, { color: mode === m.key ? '#fff' : colors.foreground }]}>
                {m.label}
              </Text>
            </TouchableOpacity>
          ))}
        </ScrollView>

        <SectionCard title="Parameters">
          {FIELD_MAP[mode].map(f => <InputRow key={f.field} {...f} />)}
        </SectionCard>

        {result ? (
          <SectionCard title="Result">
            <Text style={[styles.resultText, { color: colors.blue }]}>{result}</Text>
          </SectionCard>
        ) : null}

        {chartData.length > 1 && (
          <SectionCard title="Chart">
            <LineChart
              data={chartData}
              width={CHART_W}
              height={200}
              color={colors.blue}
              thickness={2}
              dataPointsColor={colors.blue}
              dataPointsRadius={chartData.length > 15 ? 0 : 3}
              xAxisColor={colors.border}
              yAxisColor={colors.border}
              xAxisLabelTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
              yAxisTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
              hideDataPoints={chartData.length > 20}
              curved
              initialSpacing={4}
              endSpacing={4}
              noOfSections={4}
            />
          </SectionCard>
        )}
      </ScrollView>
    </KeyboardAvoidingView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  modeChip: { paddingHorizontal: 14, paddingVertical: 8, borderRadius: 20, borderWidth: 1, marginRight: 8 },
  modeText: { fontSize: 12, fontWeight: '600' },
  inputRow: { flexDirection: 'row', alignItems: 'center', marginBottom: 10, gap: 10 },
  inputLabel: { flex: 1, fontSize: 13 },
  inputField: { width: 120, borderWidth: 1, borderRadius: 8, padding: 8, fontSize: 14, textAlign: 'right' },
  resultText: { fontSize: 15, fontWeight: '600', lineHeight: 24 },
});
