'use client';

import Link from "next/link";
import Image from "next/image";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import CSTRVisualization from "@/components/CSTRVisualization";

// Define the paths for the thumbnails
const mccabeThumbnailPath = "/thumbnails/mccabe-thiele-thumbnail.png";
const azeotropeThumbnailPath = "/thumbnails/azeotrope-finder-thumbnail.png";
const kineticsThumbnailPath = "/thumbnails/kinetics-thumbnail.png";
const dynamicsThumbnailPath = "/thumbnails/process-dynamics-thumbnail.png";
const pidThumbnailPath = "/thumbnails/pid-tuning-thumbnail.png";
const compoundPropertiesThumbnailPath = "/thumbnails/compound-properties-thumbnail.png"; // Added new thumbnail path
const residueCurveMapThumbnailPath = "/thumbnails/residue-curve-map-thumbnail.png"; // Placeholder for new thumbnail
const reactorDesignThumbnailPath = "/thumbnails/reactor-design-thumbnail.png";

export default function Page() {
  const allSimulations = [
    {
      name: "McCabe-Thiele",
      path: "/toolbox/mccabe-thiele",
      description: "Select components and specify operating conditions to visualize distillation processes with accurate equilibrium diagrams.",
      thumbnailPath: mccabeThumbnailPath,
      isUpdated: true, // Added UPDATED badge
    },
    {
      name: "Reactor Design",
      path: "/toolbox/reactor-design",
      description: "Design and analyze chemical reactors. Calculate conversions and outlet flow rates for CSTR and PFR configurations.",
      thumbnailPath: reactorDesignThumbnailPath,
      isNew: true, // Added NEW badge
    },
    {
      name: "Azeotrope Finder",
      path: "/toolbox/azeotrope-finder",
      description: "Predict and visualize azeotropic behavior of binary mixtures using various thermodynamic models.",
      thumbnailPath: azeotropeThumbnailPath,
      isNew: true, // Existing NEW badge
    },
    {
      name: "Compound Properties",
      path: "/toolbox/compound-properties", // Assuming this is the correct path
      description: "Fetch, plot, and compare various physical and thermodynamic properties of chemical compounds.",
      thumbnailPath: compoundPropertiesThumbnailPath,
      isNew: true, // Added NEW badge
    },
    {
      name: "Residue Curve Map",
      path: "/toolbox/residue-curve-map",
      description: "Visualize and analyze residue curve maps for ternary mixtures, aiding in distillation sequence design.",
      thumbnailPath: residueCurveMapThumbnailPath, // Use the new thumbnail path
      isNew: true, 
    },
    {
      name: "Reaction Kinetics",
      path: "/toolbox/kinetics",
      description: "Interactive simulator for chemical reaction kinetics. Model various reaction types and visualize concentration profiles over time.",
      thumbnailPath: kineticsThumbnailPath
    },
    {
      name: "Process Dynamics",
      path: "/toolbox/process-dynamics",
      description: "Simulate process dynamics with various inputs and understand system behavior in chemical processes.",
      thumbnailPath: dynamicsThumbnailPath
    },
    {
      name: "PID Tuning",
      path: "/toolbox/pid-tuning",
      description: "Interactive PID controller tuning simulation. Adjust parameters and observe system response in real-time.",
      thumbnailPath: pidThumbnailPath
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
                  {simulation.isNew && (
                    <div className="absolute top-3 right-[-28px] transform rotate-45 bg-orange-500 text-white text-xs font-semibold py-2 px-10 shadow-lg z-10">
                      NEW
                    </div>
                  )}
                  {simulation.isUpdated && ( // Added UPDATED badge rendering
                    (<div className="absolute top-3 right-[-32px] transform rotate-45 bg-yellow-400 text-black text-xs font-semibold py-2 px-8 shadow-lg z-10">UPDATED
                                          </div>)
                  )}
                  <CardHeader>
                    <div className="relative mx-auto mb-4">
                      <div className="aspect-square w-64 h-64 rounded-3xl overflow-hidden relative">
                        {simulation.name === "Reactor Design" ? (
                          <CSTRVisualization className="w-full h-full" showLabel={true} />
                        ) : (
                          <Image
                            src={simulation.thumbnailPath}
                            alt={`${simulation.name} Thumbnail`}
                            layout="fill"
                            objectFit="contain"
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
