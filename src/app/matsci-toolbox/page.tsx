'use client';

import Link from "next/link";
import { useTheme } from "next-themes";
import { useEffect, useState, useMemo } from "react";
import SearchBar from "@/components/SearchBar"
import FeedbackDialog from "@/components/FeedbackDialog"
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import MolecularDynamicsThumbnail from "@/components/MolecularDynamicsThumbnail";

export default function Page() {
  const { resolvedTheme } = useTheme();
  const [mounted, setMounted] = useState(false);
  const [searchQuery, setSearchQuery] = useState("");

  useEffect(() => {
    setMounted(true);
  }, []);

  const allSimulations = [
    {
      name: "Molecular Dynamics",
      path: "/toolbox/molecular-dynamics",
      description: "Interactive molecular dynamics simulator for studying particle interactions, energy evolution, and radial distribution functions.",
    },
  ];

  const filteredSimulations = useMemo(() => {
    if (!searchQuery) return allSimulations;
    const lowerQuery = searchQuery.toLowerCase();

    return allSimulations
      .filter((sim) =>
        sim.name.toLowerCase().includes(lowerQuery) ||
        sim.description.toLowerCase().includes(lowerQuery)
      )
      .sort((a, b) => {
        const nameA = a.name.toLowerCase();
        const nameB = b.name.toLowerCase();
        const descA = a.description.toLowerCase();
        const descB = b.description.toLowerCase();

        const getScore = (name: string, desc: string) => {
          if (name === lowerQuery) return 100;
          if (name.startsWith(lowerQuery)) return 80;
          if (name.includes(lowerQuery)) return 60;
          if (desc.startsWith(lowerQuery)) return 40;
          if (desc.includes(lowerQuery)) return 20;
          return 0;
        };

        const scoreA = getScore(nameA, descA);
        const scoreB = getScore(nameB, descB);

        return scoreB - scoreA;
      });
  }, [searchQuery, allSimulations]);

  return (
    <main className="flex flex-col items-center justify-between pt-6 px-6">
      {allSimulations.length > 0 && (
        <div className="w-full">
          {/* Search bar (left) and Submit Feedback (right) */}
          <div className="mb-6 flex justify-between items-start">
            <div>
              <SearchBar value={searchQuery} onChange={setSearchQuery} className="w-64" />
            </div>
            <div>
              <FeedbackDialog />
            </div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 w-full">
            {filteredSimulations.map((simulation, index) => (
              <Link
                key={index}
                href={simulation.path}
                className="h-full">
                <Card className="h-full hover:shadow-lg transition-shadow bg-card border border-border flex flex-col relative overflow-hidden">

                  <CardHeader>
                    <div className="relative mx-auto mb-4">
                      <div className="aspect-square w-64 h-64 rounded-3xl overflow-hidden relative">
                        <MolecularDynamicsThumbnail
                          className="w-full h-full"
                          isDark={mounted ? resolvedTheme !== 'light' : true}
                        />
                      </div>
                    </div>
                    <CardTitle>{simulation.name}</CardTitle>
                  </CardHeader>
                  <CardContent className="flex-grow">
                    <CardDescription>{simulation.description}</CardDescription>
                  </CardContent>
                </Card>
              </Link>
            ))}
          </div>
        </div>
      )}
    </main>
  );
}
