'use client';

import Link from "next/link";
import Image from "next/image";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";

// Define the paths for the thumbnails
const mccabeThumbnailPath = "/thumbnails/mccabe-thiele-thumbnail.png";
const azeotropeThumbnailPath = "/thumbnails/azeotrope-finder-thumbnail.png";
const kineticsThumbnailPath = "/thumbnails/kinetics-thumbnail.png";
const dynamicsThumbnailPath = "/thumbnails/process-dynamics-thumbnail.png";
const pidThumbnailPath = "/thumbnails/pid-tuning-thumbnail.png";

// Introduction text for the page
const pageDescription = "A portfolio of interactive chemical engineering simulations developed by me!\nThese tools demonstrate practical application of core concepts in McCabe-Thiele distillation, azeotrope prediction, reaction kinetics, process dynamics, and PID controller tuning.\n Under active development, with more simulations to come!";

export default function Page() {
  const allSimulations = [
    {
      name: "McCabe-Thiele",
      path: "/simulations/mccabe-thiele",
      description: "Select components and specify operating conditions to visualize distillation processes with accurate equilibrium diagrams.",
      thumbnailPath: mccabeThumbnailPath
    },
    {
      name: "Azeotrope Finder",
      path: "/simulations/azeotrope-finder",
      description: "Predict and visualize azeotropic behavior of binary mixtures using various thermodynamic models.",
      thumbnailPath: azeotropeThumbnailPath
    },
    {
      name: "Reaction Kinetics",
      path: "/simulations/kinetics",
      description: "Interactive simulator for chemical reaction kinetics. Model various reaction types and visualize concentration profiles over time.",
      thumbnailPath: kineticsThumbnailPath
    },
    {
      name: "Process Dynamics",
      path: "/simulations/process-dynamics",
      description: "Simulate process dynamics with various inputs and understand system behavior in chemical processes.",
      thumbnailPath: dynamicsThumbnailPath
    },
    {
      name: "PID Tuning",
      path: "/simulations/pid-tuning",
      description: "Interactive PID controller tuning simulation. Adjust parameters and observe system response in real-time.",
      thumbnailPath: pidThumbnailPath
    }
  ];

  // Display the first four simulations as "featured" according to the new order
  const featuredSimulations = allSimulations.slice(0, 4);

  return (
    <main className="flex flex-col items-center justify-between pt-12 px-6">
      <p className="mb-10 text-lg text-center text-muted-foreground max-w-3xl mx-auto whitespace-pre-line">
        {pageDescription}
      </p>
      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 w-full">
        {featuredSimulations.map((simulation, index) => (
          <Link
            key={index}
            href={simulation.path}
            className="block h-full"
            >
            <Card className="h-full hover:shadow-lg transition-shadow bg-card border border-border flex flex-col">
              <CardHeader>
                <div className="relative mx-auto mb-4">
                  <div className="aspect-square w-64 h-64 rounded-3xl overflow-hidden relative">
                    <Image
                      src={simulation.thumbnailPath}
                      alt={`${simulation.name} Thumbnail`}
                      layout="fill"
                      objectFit="contain"
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
      
      {/* Optionally, list all simulations below or link to a page with all of them if the list grows */}
      {/* For now, we are only showing the featured ones based on the current layout */}
      {allSimulations.length > featuredSimulations.length && (
        <div className="mt-12 w-full">
          <h2 className="text-xl font-bold mb-4 text-center">All Simulations</h2>
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 w-full">
            {allSimulations.map((simulation, index) => (
              <Link
                key={index}
                href={simulation.path}
                className="block h-full"
              >
                <Card className="h-full hover:shadow-lg transition-shadow bg-card border border-border flex flex-col">
                  <CardHeader>
                    <div className="relative mx-auto mb-4">
                      <div className="aspect-square w-64 h-64 rounded-3xl overflow-hidden relative">
                        <Image
                          src={simulation.thumbnailPath}
                          alt={`${simulation.name} Thumbnail`}
                          layout="fill"
                          objectFit="contain"
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
