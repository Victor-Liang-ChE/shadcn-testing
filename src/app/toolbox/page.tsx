'use client';

import Link from "next/link";
import Image from "next/legacy/image";
import { useTheme } from "next-themes";
import { useEffect, useState, useMemo, Suspense } from "react";
import { useSearchParams } from "next/navigation";
import SearchBar from "@/components/SearchBar"
import FeedbackDialog from "@/components/FeedbackDialog"
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import CSTRVisualization from "@/components/CSTRVisualization";
import MolecularDynamicsThumbnail from "@/components/MolecularDynamicsThumbnail";

// Define the paths for the thumbnails
const mccabeThumbnailPath = "/thumbnails/mccabe-thiele-thumbnail.png";
const mccabeThumbnailLightPath = "/thumbnails/mccabe-thiele-thumbnail-light.png";
const azeotropeThumbnailPath = "/thumbnails/azeotrope-finder-thumbnail.png";
const azeotropeThumbnailLightPath = "/thumbnails/azeotrope-finder-thumbnail-light.png";
const kineticsThumbnailPath = "/thumbnails/kinetics-thumbnail.png";
const kineticsThumbnailLightPath = "/thumbnails/kinetics-thumbnail-light.png";
const dynamicsThumbnailPath = "/thumbnails/process-dynamics-thumbnail.png";
const dynamicsThumbnailLightPath = "/thumbnails/process-dynamics-thumbnail-light.png";
const pidThumbnailPath = "/thumbnails/pid-tuning-thumbnail.png";
const pidThumbnailLightPath = "/thumbnails/pid-tuning-thumbnail-light.png";
const compoundPropertiesThumbnailPath = "/thumbnails/compound-properties-thumbnail.png";
const compoundPropertiesThumbnailLightPath = "/thumbnails/compound-properties-thumbnail-light.png";
const residueCurveMapThumbnailPath = "/thumbnails/residue-curve-map-thumbnail.png";
const residueCurveMapThumbnailLightPath = "/thumbnails/residue-curve-map-thumbnail-light.png";
const reactorDesignThumbnailPath = "/thumbnails/reactor-design-thumbnail.png";
const reactorDesignThumbnailLightPath = "/thumbnails/reactor-design-thumbnail-light.png";
const isothermsThumbnailPath = "/thumbnails/isotherms-thumbnail.png";
const isothermsThumbnailLightPath = "/thumbnails/isotherms-thumbnail-light.png";
const binaryPhaseThumbnailPath = "/thumbnails/binary-phase-diagrams-thumbnail.png";
const binaryPhaseThumbnailLightPath = "/thumbnails/binary-phase-diagrams-thumbnail-light.png";
const heatTransferThumbnailPath = "/thumbnails/1d-heat-transfer-thumbnail.png";
const heatTransferThumbnailLightPath = "/thumbnails/1d-heat-transfer-thumbnail-light.png";
const unaryPhaseThumbnailPath = "/thumbnails/unary-phase-diagrams-thumbnail.png";
const unaryPhaseThumbnailLightPath = "/thumbnails/unary-phase-diagrams-thumbnail-light.png";
const molecularDynamicsThumbnailPath = "/thumbnails/molecular-dynamics-thumbnail.png";
const molecularDynamicsThumbnailLightPath = "/thumbnails/molecular-dynamics-thumbnail-light.png";
const fugThumbnailPath = "/thumbnails/FUG-thumbnail.png";
const fugThumbnailLightPath = "/thumbnails/FUG-thumbnail-light.png";

function ToolboxContent() {
  const { resolvedTheme } = useTheme();
  const searchParams = useSearchParams();
  const category = searchParams.get('category') || 'cheme';
  const [mounted, setMounted] = useState(false);
  const [searchQuery, setSearchQuery] = useState("");

  useEffect(() => {
    setMounted(true);
  }, []);

  const allSimulations = useMemo(() => [
    {
      name: "Reactor Design",
      path: "/toolbox/reactor-design",
      description: "Design and analyze chemical reactors. Calculate conversions and outlet flow rates for CSTR and PFR configurations.",
      thumbnailPath: reactorDesignThumbnailPath,
      thumbnailLightPath: reactorDesignThumbnailLightPath,
      category: "cheme"
    },
    {
      name: "Residue Curve Map",
      path: "/toolbox/residue-curve-map",
      description: "Visualize and analyze residue curve maps for ternary mixtures, aiding in distillation sequence design.",
      thumbnailPath: residueCurveMapThumbnailPath,
      thumbnailLightPath: residueCurveMapThumbnailLightPath,
      category: "cheme"
    },
    {
      name: "PID Tuning",
      path: "/toolbox/pid-tuning",
      description: "Interactive PID controller tuning simulation. Adjust parameters and observe system response in real-time.",
      thumbnailPath: pidThumbnailPath,
      thumbnailLightPath: pidThumbnailLightPath,
      category: "cheme"
    },
    {
      name: "McCabe-Thiele",
      path: "/toolbox/mccabe-thiele",
      description: "Select components and specify operating conditions to visualize distillation processes with accurate equilibrium diagrams.",
      thumbnailPath: mccabeThumbnailPath,
      thumbnailLightPath: mccabeThumbnailLightPath,
      category: "cheme"
    },
    {
      name: "Molecular Dynamics",
      path: "/toolbox/molecular-dynamics",
      description: "Interactive molecular dynamics simulator for studying particle interactions, energy evolution, and radial distribution functions.",
      thumbnailPath: molecularDynamicsThumbnailPath,
      thumbnailLightPath: molecularDynamicsThumbnailLightPath,
      category: "matsci"
    },
    {
      name: "Unary Phase Diagram",
      path: "/toolbox/unary-phase-diagrams",
      description: "Interactive phase diagrams for pure compounds showing vaporization, fusion, and sublimation curves with key thermodynamic points.",
      thumbnailPath: unaryPhaseThumbnailPath,
      thumbnailLightPath: unaryPhaseThumbnailLightPath,
      category: "cheme"
    },
    {
      name: "Binary Phase Diagram",
      path: "/toolbox/binary-phase-diagrams",
      description: "Generate and visualize binary phase diagrams using various thermodynamic models.",
      thumbnailPath: binaryPhaseThumbnailPath,
      thumbnailLightPath: binaryPhaseThumbnailLightPath,
      category: "cheme"
    },
    {
      name: "Compound Properties",
      path: "/toolbox/compound-properties",
      description: "Fetch, plot, and compare various physical and thermodynamic properties of chemical compounds.",
      thumbnailPath: compoundPropertiesThumbnailPath,
      thumbnailLightPath: compoundPropertiesThumbnailLightPath,
      category: "cheme"
    },
    {
      name: "Azeotrope Finder",
      path: "/toolbox/azeotrope-finder",
      description: "Predict and visualize azeotropic behavior of binary mixtures using various thermodynamic models.",
      thumbnailPath: azeotropeThumbnailPath,
      thumbnailLightPath: azeotropeThumbnailLightPath,
      category: "cheme"
    },
    {
      name: "Isotherms",
      path: "/toolbox/isotherms",
      description: "Explore and visualize adsorption isotherms including Langmuir, Freundlich, and Temkin models for surface chemistry analysis.",
      thumbnailPath: isothermsThumbnailPath,
      thumbnailLightPath: isothermsThumbnailLightPath,
      category: "cheme"
    },
    {
      name: "Reaction Kinetics",
      path: "/toolbox/kinetics",
      description: "Interactive simulator for chemical reaction kinetics. Model various reaction types and visualize concentration profiles over time.",
      thumbnailPath: kineticsThumbnailPath,
      thumbnailLightPath: kineticsThumbnailLightPath,
      category: "cheme"
    },
    {
      name: "Process Dynamics",
      path: "/toolbox/process-dynamics",
      description: "Simulate process dynamics with various inputs and understand system behavior in chemical processes.",
      thumbnailPath: dynamicsThumbnailPath,
      thumbnailLightPath: dynamicsThumbnailLightPath,
      category: "cheme"
    },
    {
      name: "1D Heat Transfer",
      path: "/toolbox/1d-heat-transfer",
      description: "Interactive heat transfer visualization through multiple layers with temperature controls.",
      thumbnailPath: heatTransferThumbnailPath,
      thumbnailLightPath: heatTransferThumbnailLightPath,
      category: "cheme"
    },
    {
      name: "FUG",
      path: "/toolbox/FUG",
      description: "Multicomponent flash and distillation design using Fenske, Underwood, and Gilliland relations.",
      thumbnailPath: fugThumbnailPath,
      thumbnailLightPath: fugThumbnailLightPath,
      category: "cheme"
    }
  ], []);

  const filteredSimulations = useMemo(() => {
    let simulations = allSimulations.filter(sim => sim.category === category);
    
    if (!searchQuery) return simulations;
    const lowerQuery = searchQuery.toLowerCase();

    return simulations
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
  }, [searchQuery, allSimulations, category]);

  return (
    <main className="flex flex-col items-center justify-between pt-6 px-6">
      {/* Search bar (left) and Submit Feedback (right) */}
      <div className="w-full mb-6 flex justify-between items-start">
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
                    {simulation.name === "Reactor Design" ? (
                      <CSTRVisualization className="w-full h-full" showLabel={true} />
                    ) : simulation.name === "Molecular Dynamics" ? (
                      <MolecularDynamicsThumbnail
                        className="w-full h-full"
                        isDark={mounted ? resolvedTheme !== 'light' : true}
                      />
                    ) : (
                      <Image
                        src={mounted && resolvedTheme === 'light' ? simulation.thumbnailLightPath : simulation.thumbnailPath}
                        alt={`${simulation.name} Thumbnail`}
                        layout="fill"
                        objectFit={simulation.name === "Residue Curve Map" ? "cover" : "contain"}
                      />
                    )}
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
    </main>
  );
}

export default function Page() {
  return (
    <Suspense fallback={<div className="flex justify-center items-center h-screen">Loading Toolbox...</div>}>
      <ToolboxContent />
    </Suspense>
  );
}

