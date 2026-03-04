import React from 'react';
import { NavigationContainer, DefaultTheme, DarkTheme } from '@react-navigation/native';
import { SafeAreaProvider } from 'react-native-safe-area-context';
import { StatusBar } from 'expo-status-bar';
import { ThemeProvider, useTheme } from './src/theme/ThemeContext';
import RootNavigator from './src/navigation/RootNavigator';

function AppWithNavigation() {
  const { colorScheme, colors } = useTheme();

  const navTheme =
    colorScheme === 'dark'
      ? {
          ...DarkTheme,
          colors: {
            ...DarkTheme.colors,
            background: colors.background,
            card: colors.navbarBackground,
            text: colors.foreground,
            border: colors.border,
            primary: colors.blue,
            notification: colors.blue,
          },
        }
      : {
          ...DefaultTheme,
          colors: {
            ...DefaultTheme.colors,
            background: colors.background,
            card: colors.navbarBackground,
            text: colors.foreground,
            border: colors.border,
            primary: colors.blue,
            notification: colors.blue,
          },
        };

  return (
    <NavigationContainer theme={navTheme}>
      <StatusBar style={colorScheme === 'dark' ? 'light' : 'dark'} />
      <RootNavigator />
    </NavigationContainer>
  );
}

export default function App() {
  return (
    <SafeAreaProvider>
      <ThemeProvider>
        <AppWithNavigation />
      </ThemeProvider>
    </SafeAreaProvider>
  );
}
