'use client';

import React from 'react';
import { motion } from 'framer-motion';
import { X } from 'lucide-react';
import { wcBadgeStyle } from '@/lib/wc-colors';

const TWO_YEARS_MS = 2 * 365.25 * 24 * 60 * 60 * 1000;

const WC_ABBREV: Record<string, string> = {
  heavyweight: 'HW',
  lightheavyweight: 'LHW',
  middleweight: 'MW',
  welterweight: 'WW',
  lightweight: 'LW',
  featherweight: 'FEW',
  bantamweight: 'BW',
  flyweight: 'FLW',
  women_bantamweight: 'W-BW',
  women_featherweight: 'W-FEW',
  women_flyweight: 'W-FLW',
  women_strawweight: 'W-SW',
};

export interface BinFighter {
  name: string;
  currentElo: number;
  primaryWeightClass: string;
  lastFightDate: string;
}

interface BinFightersModalProps {
  lo: number;
  hi: number;
  fighters: BinFighter[];
  onClose: () => void;
  onSelectFighter: (name: string) => void;
  renderFighterName?: (fighter: BinFighter) => React.ReactNode;
}

export function BinFightersModal({
  lo,
  hi,
  fighters,
  onClose,
  onSelectFighter,
  renderFighterName,
}: BinFightersModalProps) {
  return (
    <motion.div
      key="bin-overlay-backdrop"
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      exit={{ opacity: 0 }}
      transition={{ duration: 0.2 }}
      className="fixed inset-0 z-50 flex items-center justify-center bg-black/60 backdrop-blur-sm p-4 md:p-8"
      onClick={onClose}
    >
      <motion.div
        initial={{ opacity: 0, scale: 0.96, y: 12 }}
        animate={{ opacity: 1, scale: 1, y: 0 }}
        exit={{ opacity: 0, scale: 0.96, y: 12 }}
        transition={{ duration: 0.18 }}
        className="bg-card text-card-foreground shadow-xl w-full max-w-lg max-h-[80vh] flex flex-col rounded-xl overflow-hidden"
        onClick={(e) => e.stopPropagation()}
      >
        {/* Header */}
        <div className="flex items-center justify-between px-5 py-3.5 border-b border-border flex-shrink-0">
          <div>
            <h2 className="font-bold text-base" style={{ fontFamily: 'Merriweather Sans, sans-serif' }}>
              Elo {lo}–{hi}
            </h2>
            <p className="text-xs text-muted-foreground mt-0.5">
              {fighters.length} fighter{fighters.length !== 1 ? 's' : ''}
            </p>
          </div>
          <button
            onClick={onClose}
            className="text-muted-foreground hover:text-foreground p-1.5 rounded-full hover:bg-muted transition-colors"
            aria-label="Close"
          >
            <X size={18} />
          </button>
        </div>

        {/* Scrollable fighter list */}
        <div className="overflow-y-auto flex-1 min-h-0">
          {fighters.length === 0 ? (
            <div className="flex items-center justify-center py-12 text-muted-foreground text-sm">
              No fighters in this range.
            </div>
          ) : (
            <table className="w-full text-sm">
              <thead className="sticky top-0 bg-card z-10 border-b border-border">
                <tr>
                  <th className="py-2 px-3 text-center w-8 font-semibold text-muted-foreground text-xs">#</th>
                  <th className="py-2 px-3 text-left font-semibold text-xs text-muted-foreground">Fighter</th>
                  <th className="py-2 px-3 text-center font-semibold text-xs text-muted-foreground">Elo</th>
                </tr>
              </thead>
              <tbody>
                {fighters.map((f, i) => {
                  const isInactive =
                    !!f.lastFightDate &&
                    Date.now() - new Date(f.lastFightDate).getTime() > TWO_YEARS_MS;
                  const abbrev = WC_ABBREV[f.primaryWeightClass];
                  return (
                    <tr
                      key={f.name}
                      className="border-b border-border/40 last:border-0 hover:bg-muted/30 transition-colors"
                    >
                      <td className="py-2.5 px-3 text-center text-muted-foreground text-xs tabular-nums">
                        {i + 1}
                      </td>
                      <td className="py-2.5 px-3">
                        {renderFighterName ? (
                          renderFighterName(f)
                        ) : (
                          <button
                            className="flex items-center gap-2 text-left hover:text-primary transition-colors w-full group"
                            onClick={() => onSelectFighter(f.name)}
                          >
                            <span className="font-medium group-hover:underline truncate max-w-[160px]">
                              {f.name}
                            </span>
                            {abbrev && (
                              <span
                                className="text-[10px] font-semibold px-1.5 py-0.5 rounded border flex-shrink-0"
                                style={wcBadgeStyle(f.primaryWeightClass)}
                              >
                                {abbrev}
                              </span>
                            )}
                            {isInactive && (
                              <span className="text-[10px] font-semibold px-1.5 py-0.5 rounded border border-gray-400/40 bg-gray-400/10 text-gray-500 dark:text-gray-400 flex-shrink-0">
                                Inactive
                              </span>
                            )}
                          </button>
                        )}
                      </td>
                      <td className="py-2.5 px-3 text-center font-bold tabular-nums">
                        {f.currentElo}
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          )}
        </div>

        {/* Footer hint */}
        <div className="px-5 py-2.5 border-t border-border flex-shrink-0">
          <p className="text-xs text-muted-foreground text-center">
            Click a fighter to view their profile
          </p>
        </div>
      </motion.div>
    </motion.div>
  );
}
