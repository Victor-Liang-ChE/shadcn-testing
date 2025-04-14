'use client';

import Link from "next/link";
import Image from "next/image"; // Import Image component
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";

// Define the paths for the thumbnails
const kineticsThumbnailPath = "/thumbnails/kinetics-thumbnail.png";
const mccabeThumbnailPath = "/thumbnails/mccabe-thiele-thumbnail.png";
const dynamicsThumbnailPath = "/thumbnails/process-dynamics-thumbnail.png";
const pidThumbnailPath = "/thumbnails/pid-tuning-thumbnail.png";


export default function Page() {
  const simulations = [
    {
      name: "Reaction Kinetics",
      path: "/simulations/kinetics",
      description: "Interactive simulator for chemical reaction kinetics. Model various reaction types and visualize concentration profiles over time.",
      thumbnailPath: kineticsThumbnailPath // Added path
    },
    {
      name: "McCabe-Thiele",
      path: "/simulations/mccabe-thiele",
      description: "Select components and specify operating conditions to visualize distillation processes with accurate equilibrium diagrams.",
      thumbnailPath: mccabeThumbnailPath // Added path
    },
    {
      name: "Process Dynamics", // Corrected name if needed, ensure it matches image file
      path: "/simulations/process-dynamics", // Corrected path based on previous file structure
      description: "Simulate process dynamics with various inputs and understand system behavior in chemical processes.",
      thumbnailPath: dynamicsThumbnailPath // Added path
    },
    {
      name: "PID Tuning",
      path: "/simulations/pid-tuning",
      description: "Interactive PID controller tuning simulation. Adjust parameters and observe system response in real-time.",
      thumbnailPath: pidThumbnailPath // Added path
    }
  ];

  return (
    <main className="flex flex-col items-center justify-between pt-24 px-6">
      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 w-full">
        {simulations.map((simulation, index) => (
          <Link key={index} href={simulation.path} className="block h-full">
            <Card className="h-full hover:shadow-lg transition-shadow bg-card border border-border flex flex-col">
              <CardHeader>
                {/* Outer container for centering the block */}
                <div className="relative mx-auto mb-4">
                  {/* Inner container for aspect ratio, size, rounding, clipping, and relative positioning */}
                  <div className="aspect-square w-64 h-64 rounded-3xl overflow-hidden relative"> {/* Added bg-muted as fallback */}
                    {/* Use simulation.thumbnailPath for all images */}
                    <Image
                      src={simulation.thumbnailPath}
                      alt={`${simulation.name} Thumbnail`}
                      layout="fill"
                      objectFit="contain" // Use contain to ensure the whole image fits
                      // Optional: Add unoptimized if images are static and don't need Next.js optimization
                      // unoptimized={true}
                    />
                    {/* Removed the placeholder text div */}
                  </div>
                </div>
                <CardTitle>{simulation.name}</CardTitle>
              </CardHeader>
              {/* ...existing CardContent... */}
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
