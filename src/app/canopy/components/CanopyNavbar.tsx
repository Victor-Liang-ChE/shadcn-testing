'use client';

import { useState, useRef, useEffect } from 'react';
import Link from 'next/link';
import { Play, ArrowLeft, ArrowRight, Save, AlertTriangle, X } from 'lucide-react';
import { Button } from '@/components/ui/button';
import { ThemeToggle } from '@/components/ThemeToggle';
import { useCanopyStore } from '@/lib/canopy/store';
import CanopySaveSheet from './CanopySaveSheet';

export default function CanopyNavbar() {
  const { currentScreen, setScreen, solved, solveAll, solverErrors, dbWarnings } = useCanopyStore();
  const [showErrors, setShowErrors] = useState(false);
  const [saveSheetOpen, setSaveSheetOpen] = useState(false);
  const popoverRef = useRef<HTMLDivElement>(null);

  // Close popover when clicking outside
  useEffect(() => {
    if (!showErrors) return;
    const handler = (e: MouseEvent) => {
      if (popoverRef.current && !popoverRef.current.contains(e.target as Node)) {
        setShowErrors(false);
      }
    };
    document.addEventListener('mousedown', handler);
    return () => document.removeEventListener('mousedown', handler);
  }, [showErrors]);

  const allIssues = [
    ...solverErrors.map(e => ({ type: 'error' as const, msg: e })),
    ...dbWarnings.map(w => ({ type: 'warning' as const, msg: w })),
  ];

  return (
    <nav
      className="flex items-center justify-between px-4 py-2 border-b min-w-full text-navbar-foreground"
      style={{ backgroundColor: 'var(--navbar-background)' }}
    >
      {/* Left: Brand + Title */}
      <div className="flex items-center gap-3">
        <Link href="/" className="hover:opacity-80 transition-opacity">
          <img src="/favicon.ico" alt="Home" className="w-7 h-7" />
        </Link>
        <div className="h-5 w-px bg-border" />
        <span className="text-lg font-bold tracking-tight">
          🌿 Canopy
        </span>
        <span className="text-xs text-muted-foreground hidden sm:inline">
          Process Simulator
        </span>
      </div>

      {/* Center: Navigation + Solve */}
      <div className="flex items-center gap-2">
        <Button
          variant={currentScreen === 'setup' ? 'default' : 'ghost'}
          size="sm"
          onClick={() => setScreen('setup')}
          className="text-xs"
        >
          <ArrowLeft className="w-3 h-3 mr-1" />
          Setup
        </Button>

        <Button
          variant={currentScreen === 'flowsheet' ? 'default' : 'ghost'}
          size="sm"
          onClick={() => setScreen('flowsheet')}
          className="text-xs"
        >
          Flowsheet
          <ArrowRight className="w-3 h-3 ml-1" />
        </Button>

        <div className="h-5 w-px bg-border mx-1" />

        <Button
          variant="default"
          size="sm"
          onClick={solveAll}
          className="bg-green-600 hover:bg-green-700 text-white text-xs"
        >
          <Play className="w-3 h-3 mr-1" />
          Run
        </Button>

        {solved && (
          <span className="text-[10px] text-green-500 font-semibold">✓ Converged</span>
        )}

        {/* Clickable error/warning indicator with popup */}
        {allIssues.length > 0 && (
          <div className="relative" ref={popoverRef}>
            <button
              onClick={() => setShowErrors(!showErrors)}
              className="flex items-center gap-1 text-[10px] font-semibold px-1.5 py-0.5 rounded hover:bg-red-500/10 transition-colors cursor-pointer"
            >
              {solverErrors.length > 0 && (
                <span className="text-red-500">
                  ✗ {solverErrors.length} error{solverErrors.length > 1 ? 's' : ''}
                </span>
              )}
              {dbWarnings.length > 0 && (
                <span className="text-amber-500 flex items-center gap-0.5">
                  <AlertTriangle className="w-3 h-3" />
                  {dbWarnings.length}
                </span>
              )}
            </button>

            {showErrors && (
              <div className="absolute top-full right-0 mt-1 w-96 max-h-64 overflow-auto bg-popover border border-border rounded-lg shadow-xl z-50">
                <div className="flex items-center justify-between px-3 py-2 border-b border-border">
                  <span className="text-xs font-semibold">Diagnostics</span>
                  <button onClick={() => setShowErrors(false)} className="text-muted-foreground hover:text-foreground">
                    <X className="w-3.5 h-3.5" />
                  </button>
                </div>
                <div className="divide-y divide-border">
                  {allIssues.map((issue, i) => (
                    <div key={i} className="px-3 py-2 flex items-start gap-2">
                      {issue.type === 'error' ? (
                        <span className="text-red-500 text-xs mt-0.5 shrink-0">✗</span>
                      ) : (
                        <AlertTriangle className="w-3.5 h-3.5 text-amber-500 mt-0.5 shrink-0" />
                      )}
                      <p className={`text-xs leading-relaxed ${issue.type === 'error' ? 'text-red-600 dark:text-red-400' : 'text-amber-600 dark:text-amber-400'}`}>
                        {issue.msg}
                      </p>
                    </div>
                  ))}
                </div>
              </div>
            )}
          </div>
        )}
      </div>

      {/* Right: Theme + Future (save, account) */}
      <div className="flex items-center gap-2">
        <Button variant="ghost" size="sm" className="text-xs" onClick={() => setSaveSheetOpen(true)}>
          <Save className="w-3 h-3 mr-1" />
          Save
        </Button>
        <ThemeToggle />
      </div>
      <CanopySaveSheet open={saveSheetOpen} onOpenChange={setSaveSheetOpen} />
    </nav>
  );
}
