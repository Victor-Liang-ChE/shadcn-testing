// Exact color tokens from globals.css, converted to hex/rgb for React Native

export type AppColors = {
  background: string;
  card: string;
  foreground: string;
  navbarBackground: string;
  mutedForeground: string;
  border: string;
  accent: string;
  primary: string;
  blue: string;
  tabBar: string;
  tabBarActive: string;
  tabBarInactive: string;
  searchBackground: string;
};

export const Colors: { light: AppColors; dark: AppColors } = {
  light: {
    // --background: hsl(0, 0%, 100%)
    background: '#ffffff',
    // --card: hsl(210, 40%, 96.1%)
    card: '#f0f4f8',
    // --foreground: hsl(222.2, 84%, 4.9%)
    foreground: '#0a0f1e',
    navbarBackground: '#f0f4f8',
    // --muted-foreground: hsl(215.4, 16.3%, 46.9%)
    mutedForeground: '#6b7a90',
    // --border: hsl(214.3, 31.8%, 91.4%)
    border: '#dde4ed',
    // --accent: hsl(210, 40%, 96.1%)
    accent: '#f0f4f8',
    // --primary: hsl(222.2, 47.4%, 11.2%)
    primary: '#1a2744',
    blue: '#60a5fa',
    tabBar: '#f0f4f8',
    tabBarActive: '#1a2744',
    tabBarInactive: '#6b7a90',
    searchBackground: '#e8edf3',
  },
  dark: {
    // --background: hsl(215, 35%, 18%) — Dark Navy Blue
    background: '#1e2d40',
    // --card: hsl(215, 40%, 30%) — Medium Blue
    card: '#2d4561',
    // --foreground: hsl(210, 40%, 98%)
    foreground: '#f9fafb',
    navbarBackground: '#2d4561',
    // --muted-foreground: hsl(215, 20%, 65%)
    mutedForeground: '#8fa3bc',
    // --border: hsl(215, 30%, 35%)
    border: '#3d5473',
    // --accent: hsl(215, 40%, 40%)
    accent: '#3b5b80',
    // --primary: hsl(210, 60%, 90%)
    primary: '#d0e4f5',
    blue: '#60a5fa',
    tabBar: '#2d4561',
    tabBarActive: '#60a5fa',
    tabBarInactive: '#8fa3bc',
    searchBackground: '#243650',
  },
};
