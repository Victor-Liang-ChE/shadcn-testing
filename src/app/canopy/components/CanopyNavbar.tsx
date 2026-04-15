'use client';

import Link from 'next/link';
import { Play, ArrowLeft, ArrowRight, ChevronDown, Save } from 'lucide-react';
import { Button } from '@/components/ui/button';
import { ThemeToggle } from '@/components/ThemeToggle';
import { useCanopyStore } from '@/lib/canopy/store';

export default function CanopyNavbar() {
  const { currentScreen, setScreen, solved, solveAll, solverErrors } = useCanopyStore();

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
        {solverErrors.length > 0 && (
          <span className="text-[10px] text-red-500 font-semibold">
            ✗ {solverErrors.length} error{solverErrors.length > 1 ? 's' : ''}
          </span>
        )}
      </div>

      {/* Right: Theme + Future (save, account) */}
      <div className="flex items-center gap-2">
        <Button variant="ghost" size="sm" className="text-xs" disabled>
          <Save className="w-3 h-3 mr-1" />
          Save
        </Button>
        <ThemeToggle />
      </div>
    </nav>
  );
}
