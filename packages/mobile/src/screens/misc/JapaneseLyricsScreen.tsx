import React, { useState, useRef } from 'react';
import {
  ScrollView, View, Text, StyleSheet, TextInput, TouchableOpacity,
  ActivityIndicator, KeyboardAvoidingView, Platform,
} from 'react-native';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

interface WordEntry {
  word: string;
  reading?: string;
  pos?: string;
}

export default function JapaneseLyricsScreen() {
  const { colors } = useTheme();
  const [input, setInput] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [result, setResult] = useState<WordEntry[][]>([]);
  const [history, setHistory] = useState<{ input: string; result: WordEntry[][] }[]>([]);

  const analyze = async () => {
    if (!input.trim()) return;
    setLoading(true);
    setError('');
    try {
      const res = await fetch('https://victorliang.com/api/process-lyrics', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ lyrics: input }),
      });
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const data = await res.json();
      const parsed: WordEntry[][] = data.lines || data.result || [];
      setResult(parsed);
      setHistory(h => [{ input: input.trim(), result: parsed }, ...h.slice(0, 9)]);
    } catch (e: any) {
      setError(e.message || 'Analysis failed');
    } finally {
      setLoading(false);
    }
  };

  const renderLine = (line: WordEntry[], li: number) => (
    <View key={li} style={styles.line}>
      {line.map((word, wi) => (
        <View key={wi} style={styles.wordBlock}>
          {word.reading && word.reading !== word.word ? (
            <Text style={[styles.reading, { color: colors.mutedForeground }]}>{word.reading}</Text>
          ) : null}
          <Text style={[styles.word, { color: colors.foreground }]}>{word.word}</Text>
        </View>
      ))}
    </View>
  );

  return (
    <KeyboardAvoidingView
      style={{ flex: 1, backgroundColor: colors.background }}
      behavior={Platform.OS === 'ios' ? 'padding' : undefined}
    >
      <ScrollView contentContainerStyle={styles.container} keyboardShouldPersistTaps="handled">
        <SectionCard title="Japanese Lyrics Input">
          <TextInput
            style={[styles.textArea, {
              color: colors.foreground,
              borderColor: colors.border,
              backgroundColor: colors.background,
            }]}
            multiline
            numberOfLines={6}
            placeholder="Paste Japanese lyrics here…"
            placeholderTextColor={colors.mutedForeground}
            value={input}
            onChangeText={setInput}
            textAlignVertical="top"
          />
          <TouchableOpacity
            style={[styles.btn, { backgroundColor: colors.blue, opacity: loading ? 0.6 : 1 }]}
            onPress={analyze}
            disabled={loading}
          >
            {loading
              ? <ActivityIndicator color="#fff" size="small" />
              : <Text style={styles.btnText}>Analyze</Text>
            }
          </TouchableOpacity>
        </SectionCard>

        {error ? (
          <SectionCard>
            <Text style={{ color: '#ef4444', fontSize: 14 }}>{error}</Text>
          </SectionCard>
        ) : null}

        {result.length > 0 && (
          <SectionCard title="Furigana">
            {result.map((line, li) => renderLine(line, li))}
          </SectionCard>
        )}

        {history.length > 0 && (
          <SectionCard title="History">
            {history.map((h, i) => (
              <TouchableOpacity key={i} onPress={() => { setInput(h.input); setResult(h.result); }}
                style={[styles.historyItem, { borderBottomColor: colors.border }]}>
                <Text style={[styles.historyText, { color: colors.foreground }]} numberOfLines={2}>
                  {h.input}
                </Text>
              </TouchableOpacity>
            ))}
          </SectionCard>
        )}
      </ScrollView>
    </KeyboardAvoidingView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  textArea: {
    borderWidth: 1, borderRadius: 10, padding: 12, fontSize: 14,
    minHeight: 120, marginBottom: 12,
  },
  btn: { borderRadius: 10, paddingVertical: 12, alignItems: 'center' },
  btnText: { color: '#fff', fontSize: 15, fontWeight: '600' },
  line: { flexDirection: 'row', flexWrap: 'wrap', marginBottom: 8 },
  wordBlock: { alignItems: 'center', marginRight: 4, marginBottom: 6 },
  reading: { fontSize: 9, lineHeight: 12 },
  word: { fontSize: 16, lineHeight: 22 },
  historyItem: { paddingVertical: 10, borderBottomWidth: StyleSheet.hairlineWidth },
  historyText: { fontSize: 13 },
});
