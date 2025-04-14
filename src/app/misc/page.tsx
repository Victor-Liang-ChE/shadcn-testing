'use client';

import Link from "next/link";
import { Card, CardContent } from "@/components/ui/card";

export default function Page() {
  const projects = [
    { name: "Drop Chance", path: "/misc/drop-chance" },
    { name: "Japanese Lyrics Analyzer", path: "/misc/japanese-lyrics" },
    { name: "Portola Menu", path: "/misc/portola-menu" },
    { name: "Chemical Engineering Economics", path: "/misc/cheme-econ" },
    { name: "LaTeX Constructor & Converter", path: "/misc/latex" },
    { name: "Sandbox", path: "/misc/sandbox" },
    { name: "Chemistry Tools", path: "/misc/chemistry-tools" },
    { name: "Video/Audio Downloader", path: "/misc/downloader" }
  ];

  return (
    <main className="flex flex-col items-center justify-between pt-24 px-8 md:px-32">
      <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-4">
        {projects.map((project, index) => (
          <Link key={index} href={project.path} className="block">
            <Card className="hover:shadow-lg transition-shadow bg-card text-card-foreground border border-border h-36">
              <CardContent className="flex items-center justify-center h-full p-6">
                <span className="text-lg font-medium text-center">{project.name}</span>
              </CardContent>
            </Card>
          </Link>
        ))}
      </div>
    </main>
  );
}
