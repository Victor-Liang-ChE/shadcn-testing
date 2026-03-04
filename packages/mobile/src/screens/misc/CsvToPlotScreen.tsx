import React, { useState, useCallback } from 'react';
import {
  ScrollView, View, Text, StyleSheet, TouchableOpacity, ActivityIndicator,
  useWindowDimensions, Platform,
} from 'react-native';
import * as DocumentPicker from 'expo-document-picker';
import { LineChart, BarChart } from 'react-native-gifted-charts';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

type ChartType = 'line' | 'bar';

function parseCSV(text: string): { headers: string[]; rows: number[][] } {
  const lines = text.trim().split('\n').filter(l => l.trim());
  if (lines.length < 2) return { headers: [], rows: [] };
  const headers = lines[0].split(',').map(h => h.trim().replace(/^"|"$/g, ''));
  const rows = lines.slice(1).map(l =>
    l.split(',').map(v => parseFloat(v.trim()) || 0)
  );
  return { headers, rows };
}

export default function CsvToPlotScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [headers, setHeaders] = useState<string[]>([]);
  const [rows, setRows] = useState<number[][]>([]);
  const [xCol, setXCol] = useState(0);
  const [yCol, setYCol] = useState(1);
  const [chartType, setChartType] = useState<ChartType>('line');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');

  const pickFile = useCallback(async () => {
    setError('');
    setLoading(true);
    try {
      const result = await DocumentPicker.getDocumentAsync({
        type: ['text/csv', 'text/comma-separated-values', 'application/csv', '*/*'],
        copyToCacheDirectory: true,
      });
      if (result.canceled || !result.assets?.[0]) { setLoading(false); return; }
      const asset = result.assets[0];
      const response = await fetch(asset.uri);
      const text = await response.text();
      const { headers: h, rows: r } = parseCSV(text);
      if (h.length === 0) throw new Error('Could not parse CSV');
      setHeaders(h);
      setRows(r);
      setXCol(0);
      setYCol(Math.min(1, h.length - 1));
    } catch (e: any) {
      setError(e.message || 'Failed to load file');
    } finally {
      setLoading(false);
    }
  }, []);

  const chartData = rows.map((row, i) => ({
    value: row[yCol] ?? 0,
    label: rows.length <= 20 ? String((row[xCol] ?? i).toFixed(1)) : '',
  }));

  const CHART_W = width - 64;

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      <SectionCard title="Load CSV File">
        <TouchableOpacity style={[styles.pickBtn, { backgroundColor: colors.blue }]} onPress={pickFile}>
          {loading
            ? <ActivityIndicator color="#fff" />
            : <Text style={styles.pickBtnText}>📂 Choose CSV File</Text>
          }
        </TouchableOpacity>
        {error ? <Text style={{ color: '#ef4444', marginTop: 8 }}>{error}</Text> : null}
      </SectionCard>

      {headers.length > 0 && (
        <>
          <SectionCard title="X Axis">
            <View style={styles.colGrid}>
              {headers.map((h, i) => (
                <TouchableOpacity key={i} onPress={() => setXCol(i)}
                  style={[styles.colChip, {
                    backgroundColor: i === xCol ? colors.blue : colors.card,
                    borderColor: i === xCol ? colors.blue : colors.border,
                  }]}>
                  <Text style={[styles.colText, { color: i === xCol ? '#fff' : colors.foreground }]}>{h}</Text>
                </TouchableOpacity>
              ))}
            </View>
          </SectionCard>

          <SectionCard title="Y Axis">
            <View style={styles.colGrid}>
              {headers.map((h, i) => (
                <TouchableOpacity key={i} onPress={() => setYCol(i)}
                  style={[styles.colChip, {
                    backgroundColor: i === yCol ? '#22c55e' : colors.card,
                    borderColor: i === yCol ? '#22c55e' : colors.border,
                  }]}>
                  <Text style={[styles.colText, { color: i === yCol ? '#fff' : colors.foreground }]}>{h}</Text>
                </TouchableOpacity>
              ))}
            </View>
          </SectionCard>

          <View style={[styles.tabRow, { borderColor: colors.border, marginBottom: 14 }]}>
            {(['line', 'bar'] as ChartType[]).map(t => (
              <TouchableOpacity key={t} onPress={() => setChartType(t)}
                style={[styles.tab, chartType === t && { backgroundColor: colors.blue }]}>
                <Text style={[styles.tabText, { color: chartType === t ? '#fff' : colors.foreground }]}>
                  {t.charAt(0).toUpperCase() + t.slice(1)} Chart
                </Text>
              </TouchableOpacity>
            ))}
          </View>

          {chartData.length > 0 && (
            <SectionCard title={`${headers[yCol]} vs ${headers[xCol]}`}>
              {chartType === 'line' ? (
                <LineChart
                  data={chartData}
                  width={CHART_W}
                  height={220}
                  color={colors.blue}
                  thickness={2}
                  dataPointsColor={colors.blue}
                  dataPointsRadius={chartData.length > 30 ? 0 : 3}
                  hideDataPoints={chartData.length > 50}
                  xAxisColor={colors.border}
                  yAxisColor={colors.border}
                  xAxisLabelTextStyle={{ color: colors.mutedForeground, fontSize: 9 }}
                  yAxisTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
                  curved
                  initialSpacing={4}
                  endSpacing={4}
                  noOfSections={4}
                />
              ) : (
                <BarChart
                  data={chartData.slice(0, 50)}
                  width={CHART_W}
                  height={220}
                  barWidth={Math.max(4, Math.floor(CHART_W / chartData.length) - 2)}
                  frontColor={colors.blue}
                  xAxisColor={colors.border}
                  yAxisColor={colors.border}
                  xAxisLabelTextStyle={{ color: colors.mutedForeground, fontSize: 9 }}
                  yAxisTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
                  noOfSections={4}
                  initialSpacing={4}
                />
              )}
              <Text style={[styles.rowCount, { color: colors.mutedForeground }]}>
                {rows.length} rows loaded
              </Text>
            </SectionCard>
          )}
        </>
      )}
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  pickBtn: { borderRadius: 10, paddingVertical: 14, alignItems: 'center' },
  pickBtnText: { color: '#fff', fontWeight: '700', fontSize: 15 },
  colGrid: { flexDirection: 'row', flexWrap: 'wrap', gap: 8 },
  colChip: { paddingHorizontal: 12, paddingVertical: 6, borderRadius: 16, borderWidth: 1 },
  colText: { fontSize: 12, fontWeight: '600' },
  tabRow: { flexDirection: 'row', borderRadius: 10, borderWidth: 1, overflow: 'hidden' },
  tab: { flex: 1, paddingVertical: 10, alignItems: 'center' },
  tabText: { fontSize: 14, fontWeight: '600' },
  rowCount: { fontSize: 11, textAlign: 'right', marginTop: 6 },
});
