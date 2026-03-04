import React, { useState, useMemo } from 'react';
import {
  ScrollView, View, Text, StyleSheet, TouchableOpacity, TextInput,
} from 'react-native';
import SliderRow from '../../components/SliderRow';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

type Mode = 'exact' | 'atleast' | 'range';

function logFactorial(n: number): number {
  if (n <= 1) return 0;
  let s = 0;
  for (let i = 2; i <= n; i++) s += Math.log(i);
  return s;
}

function logBinom(n: number, k: number): number {
  if (k > n || k < 0) return -Infinity;
  return logFactorial(n) - logFactorial(k) - logFactorial(n - k);
}

function binomProb(n: number, k: number, p: number): number {
  if (n > 1500) {
    // Normal approximation
    const mu = n * p;
    const sigma = Math.sqrt(n * p * (1 - p));
    if (sigma < 1e-10) return k === Math.round(mu) ? 1 : 0;
    const z = (k - mu) / sigma;
    return Math.exp(-0.5 * z * z) / (sigma * Math.sqrt(2 * Math.PI));
  }
  return Math.exp(logBinom(n, k) + k * Math.log(p + 1e-300) + (n - k) * Math.log(1 - p + 1e-300));
}

function probAtLeast(n: number, k: number, p: number): number {
  if (n > 1500) {
    const mu = n * p, sigma = Math.sqrt(n * p * (1 - p));
    if (sigma < 1e-10) return n * p >= k ? 1 : 0;
    const z = (k - 0.5 - mu) / sigma;
    return 0.5 * (1 + erf(-z / Math.SQRT2));
  }
  let prob = 0;
  for (let i = k; i <= n; i++) prob += binomProb(n, i, p);
  return prob;
}

function erf(x: number): number {
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
  const result = 1 - poly * Math.exp(-x * x);
  return x >= 0 ? result : -result;
}

export default function DropChanceScreen() {
  const { colors } = useTheme();
  const [mode, setMode] = useState<Mode>('atleast');
  const [dropChance, setDropChance] = useState(10);
  const [attempts, setAttempts] = useState(50);
  const [desired, setDesired] = useState(1);
  const [desired2, setDesired2] = useState(3);

  const result = useMemo(() => {
    const p = dropChance / 100;
    const n = Math.round(attempts);
    const k1 = Math.round(desired);
    const k2 = Math.round(desired2);
    if (mode === 'exact') {
      const prob = binomProb(n, k1, p);
      return `P(exactly ${k1} drop${k1 !== 1 ? 's' : ''} in ${n} attempts) = ${(prob * 100).toPrecision(4)}%`;
    } else if (mode === 'atleast') {
      const prob = probAtLeast(n, k1, p);
      return `P(at least ${k1} drop${k1 !== 1 ? 's' : ''} in ${n} attempts) = ${(prob * 100).toPrecision(4)}%`;
    } else {
      let prob = 0;
      const lo = Math.min(k1, k2), hi = Math.max(k1, k2);
      for (let i = lo; i <= hi; i++) prob += binomProb(n, i, p);
      return `P(between ${lo} and ${hi} drops in ${n} attempts) = ${(prob * 100).toPrecision(4)}%`;
    }
  }, [mode, dropChance, attempts, desired, desired2]);

  const TABS: { key: Mode; label: string }[] = [
    { key: 'exact', label: 'Exact' },
    { key: 'atleast', label: 'At Least' },
    { key: 'range', label: 'Range' },
  ];

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      {/* Mode tabs */}
      <View style={[styles.tabRow, { borderColor: colors.border }]}>
        {TABS.map(t => (
          <TouchableOpacity
            key={t.key}
            onPress={() => setMode(t.key)}
            style={[styles.tab, mode === t.key && { backgroundColor: colors.blue }]}
          >
            <Text style={[styles.tabText, { color: mode === t.key ? '#fff' : colors.foreground }]}>
              {t.label}
            </Text>
          </TouchableOpacity>
        ))}
      </View>

      <SectionCard title="Parameters">
        <SliderRow label="Drop Chance" value={dropChance} min={0.01} max={100} step={0.01}
          unit="%" decimals={2} onChange={setDropChance} />
        <SliderRow label="Attempts" value={attempts} min={1} max={1000} step={1}
          unit="" decimals={0} onChange={setAttempts} />
        <SliderRow label={mode === 'range' ? 'Min Desired' : 'Desired Drops'}
          value={desired} min={0} max={Math.min(Math.round(attempts), 100)} step={1}
          unit="" decimals={0} onChange={setDesired} />
        {mode === 'range' && (
          <SliderRow label="Max Desired" value={desired2} min={1}
            max={Math.min(Math.round(attempts), 100)} step={1}
            unit="" decimals={0} onChange={setDesired2} />
        )}
      </SectionCard>

      <SectionCard title="Result">
        <Text style={[styles.resultText, { color: colors.blue }]}>{result}</Text>
      </SectionCard>

      <SectionCard title="How it works">
        <Text style={[styles.infoText, { color: colors.mutedForeground }]}>
          Uses the binomial probability formula B(n,k,p). For very large n ({'>'} 1500), a normal
          approximation is applied for accuracy and performance.
        </Text>
      </SectionCard>
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  tabRow: { flexDirection: 'row', borderRadius: 10, borderWidth: 1, overflow: 'hidden', marginBottom: 14 },
  tab: { flex: 1, paddingVertical: 10, alignItems: 'center' },
  tabText: { fontSize: 14, fontWeight: '600' },
  resultText: { fontSize: 16, fontWeight: '700', textAlign: 'center', lineHeight: 24 },
  infoText: { fontSize: 12, lineHeight: 18 },
});
