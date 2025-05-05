'use client';

import Link from "next/link";
import { Card, CardContent } from "@/components/ui/card";

// Description text for the page - Added newline character \n
const pageDescription = "Just a collection of various side projects, tools, and experiments I've worked on.\nDon't take them too seriously!";

export default function Page() {
const projects = [
    { name: "Drop Chance", path: "/misc/drop-chance" },
    { name: "Japanese Lyrics Analyzer", path: "/misc/japanese-lyrics" },
    { name: "Portola Menu", path: "/misc/portola-menu" },
    { name: "Chemical Engineering Economics", path: "/misc/cheme-econ" },
    { name: "LaTeX Constructor & Converter", path: "/misc/latex" },
    { name: "Sandbox", path: "/misc/sandbox" },
    { name: "Chemistry Tools", path: "/misc/chemistry-tools" },
    { name: "Video/Audio Downloader", path: "/misc/downloader" },
    { name: "UFC Championship Lineage", path: "/misc/ufc-champions" },
    { name: "Fast .csv Plotter", path: "/misc/csv-to-plot" }, // Added this line
  ];

  return (
    // pt-24 still pushes everything down from the top nav bar
    <main className="flex flex-col items-center justify-between pt-12 px-8 md:px-32">

      {/* Page Description: Added whitespace-pre-line class */}
      <p className="mb-10 text-lg text-center text-muted-foreground max-w-3xl mx-auto whitespace-pre-line">
        {pageDescription}
      </p>

      <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-4 w-full">
        {projects.map((project, index) => (
          <Link key={index} href={project.path} className="block">
            <Card className="hover:shadow-lg transition-shadow bg-card text-card-foreground border border-border h-36 flex flex-col">
              <CardContent className="flex items-center justify-center h-full p-6 flex-grow">
                <span className="text-lg font-medium text-center">{project.name}</span>
              </CardContent>
            </Card>
          </Link>
        ))}
      </div>
    </main>
  );
}
