import React from 'react';
import { View, Text, StyleSheet } from 'react-native';
import { Slider } from '@miblanchard/react-native-slider';
import { useTheme } from '../theme/ThemeContext';

interface Props {
  label: string;
  value: number;
  min: number;
  max: number;
  step?: number;
  unit?: string;
  decimals?: number;
  onChange: (v: number) => void;
}

export default function SliderRow({ label, value, min, max, step = 1, unit = '', decimals = 2, onChange }: Props) {
  const { colors } = useTheme();
  const display = Number.isInteger(value) && decimals === 0 ? value : value.toFixed(decimals);

  return (
    <View style={styles.row}>
      <View style={styles.header}>
        <Text style={[styles.label, { color: colors.foreground }]}>{label}</Text>
        <Text style={[styles.value, { color: colors.blue }]}>
          {display}{unit}
        </Text>
      </View>
      <Slider
        value={value}
        onValueChange={(v: number[]) => onChange(v[0])}
        minimumValue={min}
        maximumValue={max}
        step={step}
        minimumTrackTintColor={colors.blue}
        maximumTrackTintColor={colors.border}
        thumbTintColor={colors.blue}
      />
    </View>
  );
}

const styles = StyleSheet.create({
  row: { marginBottom: 10 },
  header: { flexDirection: 'row', justifyContent: 'space-between', marginBottom: 2 },
  label: { fontSize: 13 },
  value: { fontSize: 13, fontWeight: '700' },
});
