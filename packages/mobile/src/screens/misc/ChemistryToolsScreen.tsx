import React, { useState, useCallback } from 'react';
import {
  ScrollView, View, Text, StyleSheet, TextInput, TouchableOpacity,
  KeyboardAvoidingView, Platform,
} from 'react-native';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

// ── Element data ──────────────────────────────────────────────────────────────
const ELEMENTS: Record<string, number> = {
  H:1,He:4,Li:7,Be:9,B:11,C:12,N:14,O:16,F:19,Ne:20,
  Na:23,Mg:24,Al:27,Si:28,P:31,S:32,Cl:35.5,Ar:40,K:39,Ca:40,
  Fe:55.8,Cu:63.5,Zn:65.4,Ag:108,Au:197,Pb:207,Br:80,I:127,Hg:200.6,
  Mn:54.9,Cr:52,Co:58.9,Ni:58.7,Ba:137.3,Sr:87.6,Sn:118.7,
  Ti:47.9,Mo:95.9,W:183.8,Pt:195,Pd:106.4,Al2:54,
};

function parseMolecule(formula: string): Record<string, number> | null {
  const result: Record<string, number> = {};
  const re = /([A-Z][a-z]?)(\d*)/g;
  let m;
  while ((m = re.exec(formula)) !== null) {
    const el = m[1], cnt = m[2] ? parseInt(m[2]) : 1;
    result[el] = (result[el] || 0) + cnt;
  }
  return Object.keys(result).length > 0 ? result : null;
}

function molWeight(formula: string): number {
  const atoms = parseMolecule(formula);
  if (!atoms) return 0;
  let mw = 0;
  for (const [el, cnt] of Object.entries(atoms)) {
    mw += (ELEMENTS[el] || 0) * cnt;
  }
  return mw;
}

// Very simplified equation balancer: just reports masses in/out for given stoich
interface ReactionResult {
  reactants: { formula: string; moles: number; mass: number }[];
  products: { formula: string; moles: number; mass: number }[];
  balanced: boolean;
}

function parseReaction(eq: string): ReactionResult | null {
  const sides = eq.split('->').map(s => s.trim());
  if (sides.length !== 2) return null;

  const parseSpecies = (s: string) =>
    s.split('+').map(p => {
      const m = p.trim().match(/^(\d*\.?\d*)\s*([A-Za-z0-9]+)$/);
      if (!m) return null;
      return { formula: m[2], moles: m[1] ? parseFloat(m[1]) : 1 };
    }).filter(Boolean) as { formula: string; moles: number }[];

  const reactants = parseSpecies(sides[0]);
  const products = parseSpecies(sides[1]);
  if (!reactants.length || !products.length) return null;

  const lhsMass = reactants.reduce((s, r) => s + r.moles * molWeight(r.formula), 0);
  const rhsMass = products.reduce((s, p) => s + p.moles * molWeight(p.formula), 0);
  const balanced = Math.abs(lhsMass - rhsMass) < 0.5;

  return {
    reactants: reactants.map(r => ({ ...r, mass: r.moles * molWeight(r.formula) })),
    products: products.map(p => ({ ...p, mass: p.moles * molWeight(p.formula) })),
    balanced,
  };
}

export default function ChemistryToolsScreen() {
  const { colors } = useTheme();
  const [reaction, setReaction] = useState('2 H2 + O2 -> 2 H2O');
  const [result, setResult] = useState<ReactionResult | null>(null);
  const [mwInput, setMwInput] = useState('C6H12O6');
  const [mwResult, setMwResult] = useState<number | null>(null);

  const calculate = useCallback(() => {
    setResult(parseReaction(reaction));
  }, [reaction]);

  return (
    <KeyboardAvoidingView
      style={{ flex: 1, backgroundColor: colors.background }}
      behavior={Platform.OS === 'ios' ? 'padding' : undefined}
    >
      <ScrollView contentContainerStyle={styles.container} keyboardShouldPersistTaps="handled">
        {/* Molecular Weight Calculator */}
        <SectionCard title="Molecular Weight Calculator">
          <View style={styles.row}>
            <TextInput
              style={[styles.input, { flex: 1, color: colors.foreground, borderColor: colors.border, backgroundColor: colors.background }]}
              value={mwInput}
              onChangeText={setMwInput}
              placeholder="e.g. C6H12O6"
              placeholderTextColor={colors.mutedForeground}
              autoCapitalize="none"
              autoCorrect={false}
            />
            <TouchableOpacity
              style={[styles.calcBtn, { backgroundColor: colors.blue }]}
              onPress={() => setMwResult(molWeight(mwInput))}
            >
              <Text style={styles.calcBtnText}>Calc</Text>
            </TouchableOpacity>
          </View>
          {mwResult !== null && (
            <Text style={[styles.resultBig, { color: colors.blue }]}>
              MW = {mwResult.toFixed(2)} g/mol
            </Text>
          )}
        </SectionCard>

        {/* Stoichiometry Calculator */}
        <SectionCard title="Stoichiometry Calculator">
          <Text style={[styles.hint, { color: colors.mutedForeground }]}>
            Format: {"  2 H2 + O2 -> 2 H2O  "}
          </Text>
          <TextInput
            style={[styles.input, { color: colors.foreground, borderColor: colors.border, backgroundColor: colors.background, marginBottom: 10 }]}
            value={reaction}
            onChangeText={setReaction}
            placeholder="e.g. 2 H2 + O2 -> 2 H2O"
            placeholderTextColor={colors.mutedForeground}
            autoCapitalize="none"
            autoCorrect={false}
          />
          <TouchableOpacity style={[styles.fullBtn, { backgroundColor: colors.blue }]} onPress={calculate}>
            <Text style={styles.calcBtnText}>Analyze Reaction</Text>
          </TouchableOpacity>

          {result && (
            <View style={{ marginTop: 14 }}>
              <View style={[styles.balancedBadge, { backgroundColor: result.balanced ? '#16a34a22' : '#ef444422' }]}>
                <Text style={{ color: result.balanced ? '#16a34a' : '#ef4444', fontWeight: '700' }}>
                  {result.balanced ? '✓ Mass balanced' : '⚠ Mass not balanced'}
                </Text>
              </View>

              <Text style={[styles.subheader, { color: colors.foreground }]}>Reactants</Text>
              {result.reactants.map((r, i) => (
                <View key={i} style={styles.speciesRow}>
                  <Text style={[styles.formula, { color: colors.foreground }]}>{r.moles} {r.formula}</Text>
                  <Text style={[styles.mass, { color: colors.mutedForeground }]}>
                    MW: {molWeight(r.formula).toFixed(2)} g/mol · Mass contribution: {r.mass.toFixed(2)} g
                  </Text>
                </View>
              ))}

              <Text style={[styles.subheader, { color: colors.foreground }]}>Products</Text>
              {result.products.map((p, i) => (
                <View key={i} style={styles.speciesRow}>
                  <Text style={[styles.formula, { color: colors.foreground }]}>{p.moles} {p.formula}</Text>
                  <Text style={[styles.mass, { color: colors.mutedForeground }]}>
                    MW: {molWeight(p.formula).toFixed(2)} g/mol · Mass contribution: {p.mass.toFixed(2)} g
                  </Text>
                </View>
              ))}
            </View>
          )}
        </SectionCard>
      </ScrollView>
    </KeyboardAvoidingView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  row: { flexDirection: 'row', gap: 8, marginBottom: 10 },
  input: { borderWidth: 1, borderRadius: 10, padding: 10, fontSize: 14 },
  calcBtn: { paddingHorizontal: 16, paddingVertical: 10, borderRadius: 10 },
  calcBtnText: { color: '#fff', fontWeight: '700', fontSize: 14 },
  fullBtn: { borderRadius: 10, paddingVertical: 12, alignItems: 'center' },
  resultBig: { fontSize: 18, fontWeight: '700', textAlign: 'center', marginTop: 10 },
  hint: { fontSize: 12, marginBottom: 6 },
  balancedBadge: { borderRadius: 8, padding: 10, marginBottom: 12, alignItems: 'center' },
  subheader: { fontSize: 14, fontWeight: '700', marginTop: 10, marginBottom: 6 },
  speciesRow: { marginBottom: 8 },
  formula: { fontSize: 15, fontWeight: '600' },
  mass: { fontSize: 12, marginTop: 2 },
});
