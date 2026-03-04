import React from 'react';
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
import { MISC_PROJECTS } from '../constants/data';
import type { RootStackParamList } from '../navigation/types';

type NavProp = NativeStackNavigationProp<RootStackParamList>;

export default function MiscScreen() {
  const navigation = useNavigation<NavProp>();
  const { colors } = useTheme();
  const { width } = useWindowDimensions();

  const PADDING = 12;
  const GAP = 10;
  const cardWidth = (width - PADDING * 2 - GAP) / 2;

  const rows = [];
  for (let i = 0; i < MISC_PROJECTS.length; i += 2) {
    rows.push(MISC_PROJECTS.slice(i, i + 2));
  }

  return (
    <ScrollView
      style={{ backgroundColor: colors.background }}
      contentContainerStyle={[styles.container, { paddingHorizontal: PADDING }]}
    >
      <Text style={[styles.pageDesc, { color: colors.mutedForeground }]}>
        A collection of various side projects, tools, and experiments.
      </Text>

      {rows.map((row, idx) => (
        <View key={idx} style={[styles.row, { gap: GAP }]}>
          {row.map(project => (
            <TouchableOpacity
              key={project.path}
              style={[
                styles.card,
                {
                  backgroundColor: colors.card,
                  borderColor: colors.border,
                  width: cardWidth,
                },
              ]}
              activeOpacity={0.75}
              onPress={() => navigation.navigate(project.screenKey as any)}
            >
              <Text style={[styles.cardTitle, { color: colors.foreground }]}>
                {project.name}
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
  pageDesc: {
    fontSize: 14,
    textAlign: 'center',
    lineHeight: 20,
    marginBottom: 18,
    paddingHorizontal: 8,
  },
  row: {
    flexDirection: 'row',
    marginBottom: 10,
  },
  card: {
    borderRadius: 14,
    borderWidth: 1,
    height: 100,
    justifyContent: 'center',
    alignItems: 'center',
    padding: 12,
  },
  cardTitle: {
    fontSize: 14,
    fontWeight: '600',
    textAlign: 'center',
  },
});
