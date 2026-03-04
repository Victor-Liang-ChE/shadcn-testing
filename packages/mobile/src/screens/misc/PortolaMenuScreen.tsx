import React, { useState, useEffect } from 'react';
import {
  ScrollView, View, Text, StyleSheet, TouchableOpacity, ActivityIndicator,
} from 'react-native';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

type MealType = 'breakfast' | 'lunch' | 'dinner';
const MEALS: MealType[] = ['breakfast', 'lunch', 'dinner'];
const DAYS = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'];

interface MenuItem {
  name: string;
  station?: string;
  dietary?: string[];
}

interface DayMenu {
  day: string;
  items: MenuItem[];
}

export default function PortolaMenuScreen() {
  const { colors } = useTheme();
  const [meal, setMeal] = useState<MealType>('lunch');
  const [dayIdx, setDayIdx] = useState(0);
  const [menuData, setMenuData] = useState<Record<MealType, DayMenu[]> | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState('');

  useEffect(() => {
    const fetchMenu = async () => {
      setLoading(true);
      setError('');
      try {
        const res = await fetch('https://victorliang.com/api/portola-menu');
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        const data = await res.json();
        setMenuData(data);
      } catch (e: any) {
        setError(e.message || 'Failed to load menu');
      } finally {
        setLoading(false);
      }
    };
    fetchMenu();
  }, []);

  const currentItems: MenuItem[] = menuData
    ? (menuData[meal]?.[dayIdx]?.items ?? [])
    : [];

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      {/* Meal tabs */}
      <View style={[styles.tabRow, { borderColor: colors.border }]}>
        {MEALS.map(m => (
          <TouchableOpacity key={m} onPress={() => setMeal(m)}
            style={[styles.tab, meal === m && { backgroundColor: colors.blue }]}>
            <Text style={[styles.tabText, { color: meal === m ? '#fff' : colors.foreground }]}>
              {m.charAt(0).toUpperCase() + m.slice(1)}
            </Text>
          </TouchableOpacity>
        ))}
      </View>

      {/* Day selector */}
      <ScrollView horizontal showsHorizontalScrollIndicator={false} style={styles.dayScroll}>
        {DAYS.map((d, i) => (
          <TouchableOpacity key={d} onPress={() => setDayIdx(i)}
            style={[styles.dayChip, {
              backgroundColor: i === dayIdx ? colors.blue : colors.card,
              borderColor: i === dayIdx ? colors.blue : colors.border,
            }]}>
            <Text style={[styles.dayText, { color: i === dayIdx ? '#fff' : colors.foreground }]}>
              {d.slice(0, 3)}
            </Text>
          </TouchableOpacity>
        ))}
      </ScrollView>

      {loading && (
        <View style={styles.center}>
          <ActivityIndicator size="large" color={colors.blue} />
          <Text style={[styles.loadingText, { color: colors.mutedForeground }]}>Loading menu…</Text>
        </View>
      )}

      {error ? (
        <SectionCard>
          <Text style={[styles.errorText, { color: '#ef4444' }]}>{error}</Text>
        </SectionCard>
      ) : null}

      {!loading && !error && currentItems.length === 0 && (
        <SectionCard>
          <Text style={[styles.emptyText, { color: colors.mutedForeground }]}>
            No menu items available.
          </Text>
        </SectionCard>
      )}

      {!loading && currentItems.length > 0 && (
        <SectionCard title={`${DAYS[dayIdx]} ${meal.charAt(0).toUpperCase() + meal.slice(1)}`}>
          {currentItems.map((item, i) => (
            <View key={i} style={[styles.menuItem, { borderBottomColor: colors.border }]}>
              <Text style={[styles.itemName, { color: colors.foreground }]}>{item.name}</Text>
              {item.station ? (
                <Text style={[styles.itemStation, { color: colors.mutedForeground }]}>{item.station}</Text>
              ) : null}
              {item.dietary && item.dietary.length > 0 ? (
                <View style={styles.dietaryRow}>
                  {item.dietary.map((d, di) => (
                    <View key={di} style={[styles.dietaryBadge, { backgroundColor: colors.accent }]}>
                      <Text style={[styles.dietaryText, { color: colors.foreground }]}>{d}</Text>
                    </View>
                  ))}
                </View>
              ) : null}
            </View>
          ))}
        </SectionCard>
      )}
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  tabRow: { flexDirection: 'row', borderRadius: 10, borderWidth: 1, overflow: 'hidden', marginBottom: 12 },
  tab: { flex: 1, paddingVertical: 10, alignItems: 'center' },
  tabText: { fontSize: 14, fontWeight: '600' },
  dayScroll: { marginBottom: 14 },
  dayChip: { paddingHorizontal: 16, paddingVertical: 8, borderRadius: 20, borderWidth: 1, marginRight: 8 },
  dayText: { fontSize: 13, fontWeight: '600' },
  center: { alignItems: 'center', paddingVertical: 40 },
  loadingText: { marginTop: 10, fontSize: 14 },
  errorText: { fontSize: 14, textAlign: 'center' },
  emptyText: { fontSize: 14, textAlign: 'center' },
  menuItem: { paddingVertical: 10, borderBottomWidth: StyleSheet.hairlineWidth },
  itemName: { fontSize: 14, fontWeight: '500', marginBottom: 2 },
  itemStation: { fontSize: 12, marginBottom: 4 },
  dietaryRow: { flexDirection: 'row', flexWrap: 'wrap', gap: 4 },
  dietaryBadge: { paddingHorizontal: 6, paddingVertical: 2, borderRadius: 6 },
  dietaryText: { fontSize: 10, fontWeight: '600' },
});
