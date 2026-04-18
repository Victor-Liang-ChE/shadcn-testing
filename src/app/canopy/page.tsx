'use client';

import { useEffect } from 'react';
import { useCanopyStore } from '@/lib/canopy/store';
import CanopyNavbar from './components/CanopyNavbar';
import SetupScreen from './components/SetupScreen';
import FlowsheetScreen from './components/FlowsheetScreen';

export default function CanopyPage() {
  const { currentScreen, loadPreset, compounds } = useCanopyStore();

  // Auto-load the preset on first mount if no compounds selected
  useEffect(() => {
    if (compounds.length === 0) {
      loadPreset('cumene-production');
    }
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <CanopyNavbar />
      <div className="flex-1 flex flex-col overflow-hidden">
        {currentScreen === 'setup' ? <SetupScreen /> : <FlowsheetScreen />}
      </div>
    </>
  );
}
