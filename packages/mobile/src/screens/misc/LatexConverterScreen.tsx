import React, { useState, useMemo } from 'react';
import {
  ScrollView, View, Text, StyleSheet, TextInput, TouchableOpacity,
  KeyboardAvoidingView, Platform,
} from 'react-native';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

// ── LaTeX → target language converters (ported from web app) ──────────────────

function replaceFrac(expr: string, open: string, div: string, close: string): string {
  let result = expr;
  let idx;
  while ((idx = result.indexOf('\\frac{')) !== -1) {
    // find matching brace for numerator
    let depth = 0, start = idx + 6, numEnd = -1;
    for (let i = start; i < result.length; i++) {
      if (result[i] === '{') depth++;
      else if (result[i] === '}') {
        if (depth === 0) { numEnd = i; break; }
        depth--;
      }
    }
    if (numEnd === -1) break;
    const num = result.slice(start, numEnd);
    if (result[numEnd + 1] !== '{') break;
    let depth2 = 0, denEnd = -1;
    for (let i = numEnd + 2; i < result.length; i++) {
      if (result[i] === '{') depth2++;
      else if (result[i] === '}') {
        if (depth2 === 0) { denEnd = i; break; }
        depth2--;
      }
    }
    if (denEnd === -1) break;
    const den = result.slice(numEnd + 2, denEnd);
    const replacement = open + replaceFrac(num, open, div, close) + div + replaceFrac(den, open, div, close) + close;
    result = result.slice(0, idx) + replacement + result.slice(denEnd + 1);
  }
  return result;
}

function latexToCommon(latex: string): string {
  return latex
    .replace(/\\sqrt\{([^}]+)\}/g, 'sqrt($1)')
    .replace(/\\sqrt\[(\d+)\]\{([^}]+)\}/g, '($2)^(1/$1)')
    .replace(/\\left\(/g, '(').replace(/\\right\)/g, ')')
    .replace(/\\left\[/g, '[').replace(/\\right\]/g, ']')
    .replace(/\\cdot/g, '*').replace(/\\times/g, '*')
    .replace(/\\div/g, '/').replace(/\\pm/g, '±')
    .replace(/\^{([^}]+)}/g, '^($1)').replace(/\^(\d)/g, '^$1')
    .replace(/_{([^}]+)}/g, '_$1')
    .replace(/\\pi/g, 'pi').replace(/\\e/g, 'e')
    .replace(/\\ln/g, 'ln').replace(/\\log/g, 'log')
    .replace(/\\exp/g, 'exp').replace(/\\sin/g, 'sin')
    .replace(/\\cos/g, 'cos').replace(/\\tan/g, 'tan')
    .replace(/\\alpha/g, 'alpha').replace(/\\beta/g, 'beta')
    .replace(/\\gamma/g, 'gamma').replace(/\\delta/g, 'delta')
    .replace(/\\theta/g, 'theta').replace(/\\lambda/g, 'lambda')
    .replace(/\\mu/g, 'mu').replace(/\\sigma/g, 'sigma')
    .replace(/\\rho/g, 'rho').replace(/\\omega/g, 'omega')
    .replace(/\\infty/g, 'Inf').replace(/{/g, '(').replace(/}/g, ')')
    .replace(/\s+/g, ' ').trim();
}

function toPython(latex: string): string {
  let expr = replaceFrac(latex, '(', ')/(', ')');
  expr = latexToCommon(expr)
    .replace(/\bpi\b/g, 'math.pi').replace(/\be\b/g, 'math.e')
    .replace(/\bsqrt\(/g, 'math.sqrt(').replace(/\bln\(/g, 'math.log(')
    .replace(/\blog\(/g, 'math.log10(').replace(/\bexp\(/g, 'math.exp(')
    .replace(/\bsin\(/g, 'math.sin(').replace(/\bcos\(/g, 'math.cos(')
    .replace(/\btan\(/g, 'math.tan(').replace(/\^/g, '**');
  return expr;
}

function toNumpy(latex: string): string {
  return toPython(latex)
    .replace(/math\./g, 'np.')
    .replace(/np\.log\(/g, 'np.log(').replace(/np\.log10\(/g, 'np.log10(');
}

function toJS(latex: string): string {
  let expr = replaceFrac(latex, '(', ')/(', ')');
  expr = latexToCommon(expr)
    .replace(/\bpi\b/g, 'Math.PI').replace(/\bInf\b/g, 'Infinity')
    .replace(/\bsqrt\(/g, 'Math.sqrt(').replace(/\bln\(/g, 'Math.log(')
    .replace(/\blog\(/g, 'Math.log10(').replace(/\bexp\(/g, 'Math.exp(')
    .replace(/\bsin\(/g, 'Math.sin(').replace(/\bcos\(/g, 'Math.cos(')
    .replace(/\btan\(/g, 'Math.tan(').replace(/\^(\([^)]+\))/g, '**$1')
    .replace(/\^(\w+)/g, '**$1');
  return expr;
}

function toMATLAB(latex: string): string {
  let expr = replaceFrac(latex, '(', ')./(', ')');
  return latexToCommon(expr).replace(/\^/g, '.^');
}

function toR(latex: string): string {
  let expr = replaceFrac(latex, '(', ')/(', ')');
  expr = latexToCommon(expr)
    .replace(/\bpi\b/g, 'pi').replace(/\bln\(/g, 'log(')
    .replace(/\blog\(/g, 'log10(').replace(/\^/g, '^');
  return expr;
}

function toExcel(latex: string): string {
  let expr = replaceFrac(latex, '(', ')/(', ')');
  expr = latexToCommon(expr)
    .replace(/\bsqrt\(/g, 'SQRT(').replace(/\bln\(/g, 'LN(')
    .replace(/\blog\(/g, 'LOG10(').replace(/\bexp\(/g, 'EXP(')
    .replace(/\bsin\(/g, 'SIN(').replace(/\bcos\(/g, 'COS(')
    .replace(/\btan\(/g, 'TAN(').replace(/\bpi\b/g, 'PI()')
    .replace(/\^/g, '^');
  return expr;
}

const TARGETS = [
  { key: 'python', label: 'Python', convert: toPython },
  { key: 'numpy', label: 'NumPy', convert: toNumpy },
  { key: 'js', label: 'JavaScript', convert: toJS },
  { key: 'matlab', label: 'MATLAB', convert: toMATLAB },
  { key: 'r', label: 'R', convert: toR },
  { key: 'excel', label: 'Excel', convert: toExcel },
];

export default function LatexConverterScreen() {
  const { colors } = useTheme();
  const [latex, setLatex] = useState('\\frac{-b \\pm \\sqrt{b^2 - 4ac}}{2a}');
  const [target, setTarget] = useState('python');

  const converted = useMemo(() => {
    const t = TARGETS.find(t => t.key === target);
    if (!t) return '';
    try { return t.convert(latex); } catch { return 'Parse error'; }
  }, [latex, target]);

  return (
    <KeyboardAvoidingView
      style={{ flex: 1, backgroundColor: colors.background }}
      behavior={Platform.OS === 'ios' ? 'padding' : undefined}
    >
      <ScrollView contentContainerStyle={styles.container} keyboardShouldPersistTaps="handled">
        <SectionCard title="LaTeX Expression">
          <TextInput
            style={[styles.textArea, {
              color: colors.foreground, borderColor: colors.border, backgroundColor: colors.background,
            }]}
            multiline
            numberOfLines={4}
            placeholder="Enter LaTeX expression…"
            placeholderTextColor={colors.mutedForeground}
            value={latex}
            onChangeText={setLatex}
            textAlignVertical="top"
            autoCapitalize="none"
            autoCorrect={false}
          />
        </SectionCard>

        <SectionCard title="Target Language">
          <View style={styles.targetGrid}>
            {TARGETS.map(t => (
              <TouchableOpacity key={t.key} onPress={() => setTarget(t.key)}
                style={[styles.targetChip, {
                  backgroundColor: target === t.key ? colors.blue : colors.card,
                  borderColor: target === t.key ? colors.blue : colors.border,
                }]}>
                <Text style={[styles.targetText, { color: target === t.key ? '#fff' : colors.foreground }]}>
                  {t.label}
                </Text>
              </TouchableOpacity>
            ))}
          </View>
        </SectionCard>

        <SectionCard title="Output">
          <Text style={[styles.output, { color: colors.blue, backgroundColor: colors.background }]}>
            {converted || '—'}
          </Text>
        </SectionCard>
      </ScrollView>
    </KeyboardAvoidingView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  textArea: {
    borderWidth: 1, borderRadius: 10, padding: 12, fontSize: 13,
    minHeight: 100, fontFamily: 'monospace',
  },
  targetGrid: { flexDirection: 'row', flexWrap: 'wrap', gap: 8 },
  targetChip: { paddingHorizontal: 14, paddingVertical: 8, borderRadius: 20, borderWidth: 1 },
  targetText: { fontSize: 13, fontWeight: '600' },
  output: {
    fontSize: 14, padding: 12, borderRadius: 10, fontFamily: 'monospace',
    lineHeight: 22, minHeight: 50,
  },
});
