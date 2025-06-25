'use client';

import Link from "next/link";
import Image from "next/image";
import { useTheme } from "next-themes";
import { useEffect, useState } from "react";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import CSTRVisualization from "@/components/CSTRVisualization";

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

export default function Page() {
  const { resolvedTheme } = useTheme();
  const [mounted, setMounted] = useState(false);

  useEffect(() => {
    setMounted(true);
  }, []);

  const allSimulations = [
    {
      name: "Reactor Design",
      path: "/toolbox/reactor-design",
      description: "Design and analyze chemical reactors. Calculate conversions and outlet flow rates for CSTR and PFR configurations.",
      thumbnailPath: reactorDesignThumbnailPath,
      thumbnailLightPath: reactorDesignThumbnailLightPath
    },
    {
      name: "Binary Phase Diagram",
      path: "/toolbox/binary-phase-diagrams",
      description: "Generate and visualize binary phase diagrams using various thermodynamic models.",
      thumbnailPath: binaryPhaseThumbnailPath,
      thumbnailLightPath: binaryPhaseThumbnailLightPath
    },
    {
      name: "PID Tuning",
      path: "/toolbox/pid-tuning",
      description: "Interactive PID controller tuning simulation. Adjust parameters and observe system response in real-time.",
      thumbnailPath: pidThumbnailPath,
      thumbnailLightPath: pidThumbnailLightPath
    },
    {
      name: "McCabe-Thiele",
      path: "/toolbox/mccabe-thiele",
      description: "Select components and specify operating conditions to visualize distillation processes with accurate equilibrium diagrams.",
      thumbnailPath: mccabeThumbnailPath,
      thumbnailLightPath: mccabeThumbnailLightPath
    },
    {
      name: "Reaction Kinetics",
      path: "/toolbox/kinetics",
      description: "Interactive simulator for chemical reaction kinetics. Model various reaction types and visualize concentration profiles over time.",
      thumbnailPath: kineticsThumbnailPath,
      thumbnailLightPath: kineticsThumbnailLightPath
    },
    {
      name: "Process Dynamics",
      path: "/toolbox/process-dynamics",
      description: "Simulate process dynamics with various inputs and understand system behavior in chemical processes.",
      thumbnailPath: dynamicsThumbnailPath,
      thumbnailLightPath: dynamicsThumbnailLightPath
    },
    {
      name: "Residue Curve Map",
      path: "/toolbox/residue-curve-map",
      description: "Visualize and analyze residue curve maps for ternary mixtures, aiding in distillation sequence design.",
      thumbnailPath: residueCurveMapThumbnailPath,
      thumbnailLightPath: residueCurveMapThumbnailLightPath
    },
    {
      name: "Azeotrope Finder",
      path: "/toolbox/azeotrope-finder",
      description: "Predict and visualize azeotropic behavior of binary mixtures using various thermodynamic models.",
      thumbnailPath: azeotropeThumbnailPath,
      thumbnailLightPath: azeotropeThumbnailLightPath
    },
    {
      name: "Compound Properties",
      path: "/toolbox/compound-properties",
      description: "Fetch, plot, and compare various physical and thermodynamic properties of chemical compounds.",
      thumbnailPath: compoundPropertiesThumbnailPath,
      thumbnailLightPath: compoundPropertiesThumbnailLightPath
    },
    {
      name: "Isotherms",
      path: "/toolbox/isotherms",
      description: "Explore and visualize adsorption isotherms including Langmuir, Freundlich, and Temkin models for surface chemistry analysis.",
      thumbnailPath: isothermsThumbnailPath,
      thumbnailLightPath: isothermsThumbnailLightPath
    }
  ];

  return (
    <main className="flex flex-col items-center justify-between pt-12 px-6">
      {/* Optionally, list all simulations below or link to a page with all of them if the list grows */}
      {allSimulations.length > 0 && (
        <div className="w-full">
          {/* <h2 className="text-xl font-bold mb-4 text-center">All Simulations</h2> */} {/* Title removed */}
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 w-full">
            {allSimulations.map((simulation, index) => (
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
        </div>
      )}
    </main>
  );
}
