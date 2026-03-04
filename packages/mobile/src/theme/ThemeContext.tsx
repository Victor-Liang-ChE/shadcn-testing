import React, { createContext, useContext, useState } from 'react';
import { useColorScheme } from 'react-native';
import { Colors, AppColors } from './colors';

type ThemeContextType = {
  colorScheme: 'light' | 'dark';
  toggleTheme: () => void;
  colors: AppColors;
};

const ThemeContext = createContext<ThemeContextType>({
  colorScheme: 'dark',
  toggleTheme: () => {},
  colors: Colors.dark,
});

export function ThemeProvider({ children }: { children: React.ReactNode }) {
  const systemScheme = useColorScheme();
  const [override, setOverride] = useState<'light' | 'dark' | null>(null);

  const colorScheme: 'light' | 'dark' =
    override ?? (systemScheme === 'dark' ? 'dark' : 'light');

  const toggleTheme = () => {
    setOverride(prev => {
      if (prev === null) return systemScheme === 'dark' ? 'light' : 'dark';
      return prev === 'dark' ? 'light' : 'dark';
    });
  };

  return (
    <ThemeContext.Provider value={{ colorScheme, toggleTheme, colors: Colors[colorScheme] }}>
      {children}
    </ThemeContext.Provider>
  );
}

export const useTheme = () => useContext(ThemeContext);
