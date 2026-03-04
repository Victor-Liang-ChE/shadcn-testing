import React, { useState } from 'react';
import {
  View,
  Text,
  ScrollView,
  StyleSheet,
  TouchableOpacity,
  useWindowDimensions,
} from 'react-native';
import { useNavigation } from '@react-navigation/native';
import { NativeStackNavigationProp } from '@react-navigation/native-stack';
import { useTheme } from '../theme/ThemeContext';
import { ALL_TOOLS, ToolItem } from '../constants/data';
import type { RootStackParamList } from '../navigation/types';

type NavProp = NativeStackNavigationProp<RootStackParamList>;

const CATEGORIES = [
  { key: 'all', label: 'All' },
  { key: 'cheme', label: 'Chemical Engineering' },
  { key: 'matsci', label: 'Materials Science' },
] as const;

type Category = (typeof CATEGORIES)[number]['key'];

export default function ToolboxScreen() {
  const navigation = useNavigation<NavProp>();
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [category, setCategory] = useState<Category>('all');

  const PADDING = 12;
  const GAP = 10;
  const cardWidth = (width - PADDING * 2 - GAP) / 2;

  const filtered = ALL_TOOLS.filter(
    t => category === 'all' || t.category === category,
  );

  const rows: ToolItem[][] = [];
  for (let i = 0; i < filtered.length; i += 2) {
    rows.push(filtered.slice(i, i + 2));
  }

  return (
    <ScrollView
      style={{ backgroundColor: colors.background }}
      contentContainerStyle={[styles.container, { paddingHorizontal: PADDING }]}
    >
      {/* Category filter chips */}
      <View style={styles.filterRow}>
        {CATEGORIES.map(cat => {
          const active = category === cat.key;
          return (
            <TouchableOpacity
              key={cat.key}
              onPress={() => setCategory(cat.key)}
              style={[
                styles.chip,
                {
                  backgroundColor: active ? colors.blue : colors.card,
                  borderColor: active ? colors.blue : colors.border,
                },
              ]}
            >
              <Text
                style={[
                  styles.chipText,
                  { color: active ? '#ffffff' : colors.foreground },
                ]}
              >
                {cat.label}
              </Text>
            </TouchableOpacity>
          );
        })}
      </View>

      {/* Tool grid */}
      {rows.map((row, idx) => (
        <View key={idx} style={[styles.row, { gap: GAP }]}>
          {row.map(tool => (
            <TouchableOpacity
              key={tool.path}
              style={[
                styles.card,
                {
                  backgroundColor: colors.card,
                  borderColor: colors.border,
                  width: cardWidth,
                },
              ]}
              activeOpacity={0.75}
              onPress={() =>
                tool.screenKey === 'WebView'
                  ? navigation.navigate('WebView', { url: `https://victorliang.com${tool.path}`, title: tool.name })
                  : navigation.navigate(tool.screenKey as any)
              }
            >
              <Text style={[styles.cardTitle, { color: colors.foreground }]}>
                {tool.name}
              </Text>
              <Text
                style={[styles.cardDesc, { color: colors.mutedForeground }]}
                numberOfLines={3}
              >
                {tool.description}
              </Text>
            </TouchableOpacity>
          ))}
          {row.length === 1 && <View style={{ width: cardWidth }} />}
        </View>
      ))}
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: {
    paddingTop: 16,
    paddingBottom: 32,
  },
  filterRow: {
    flexDirection: 'row',
    flexWrap: 'wrap',
    gap: 8,
    marginBottom: 16,
  },
  chip: {
    paddingHorizontal: 14,
    paddingVertical: 7,
    borderRadius: 20,
    borderWidth: 1,
  },
  chipText: {
    fontSize: 13,
    fontWeight: '600',
  },
  row: {
    flexDirection: 'row',
    marginBottom: 10,
  },
  card: {
    borderRadius: 14,
    borderWidth: 1,
    padding: 16,
    minHeight: 110,
    justifyContent: 'center',
  },
  cardTitle: {
    fontSize: 14,
    fontWeight: '700',
    marginBottom: 6,
  },
  cardDesc: {
    fontSize: 12,
    lineHeight: 17,
  },
});
