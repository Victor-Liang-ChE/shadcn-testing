import Image from "next/image";
import Link from "next/link";
import {
  Card,
  CardContent,
  CardDescription,
  CardFooter,
  CardHeader,
  CardTitle,
} from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from "@/components/ui/tooltip";
import CSTRVisualization from "@/components/CSTRVisualization";

// Replace external URLs with local paths from public folder
const githubLogoPath = "/logos/github.svg";
const linkedinLogoPath = "/logos/linkedin.svg";
const ucsbLogoPath = "/logos/ucsb.svg";
const pythonLogoPath = "/logos/python.svg";
const reactLogoPath = "/logos/react.svg";
const nextjsLogoPath = "/logos/nextjs.svg"; // Added Next.js logo path
const javascriptLogoPath = "/logos/javascript.svg";
const html5LogoPath = "/logos/html5.svg";
const css3LogoPath = "/logos/css.svg";
const typescriptLogoPath = "/logos/typescript.svg";
const dockerLogoPath = "/logos/docker.svg";
const supabaseLogoPath = "/logos/supabase.svg"; // Added Supabase logo path
const tailwindLogoPath = "/logos/tailwind.svg";
const shadcnLogoPath = "/logos/shadcn.svg";
const kineticsThumbnailPath = "/thumbnails/kinetics-thumbnail.png";
const mccabeThumbnailPath = "/thumbnails/mccabe-thiele-thumbnail.png";
const dynamicsThumbnailPath = "/thumbnails/process-dynamics-thumbnail.png";
const azeotropeThumbnailPath = "/thumbnails/azeotrope-finder-thumbnail.png";
const labIllustrationPath = "/images/lab-illustration.png";
const compoundPropertiesThumbnailPath = "/thumbnails/compound-properties-thumbnail.png";
const residueCurveMapThumbnailPath = "/thumbnails/residue-curve-map-thumbnail.png"; // Added for Residue Curve Map
const reactorDesignThumbnailPath = "/thumbnails/reactor-design-thumbnail.png";

export default function Home() {
  const homeFeaturedSimulations = [
    {
      name: "McCabe-Thiele",
      path: "/toolbox/mccabe-thiele",
      description:
        "Select components and specify operating conditions to visualize distillation processes with accurate equilibrium diagrams.",
      thumbnailPath: mccabeThumbnailPath,
      isUpdated: true,
    },
    {
      name: "Reactor Design",
      path: "/toolbox/reactor-design",
      description:
        "Design and analyze chemical reactors. Calculate conversions and outlet flow rates for CSTR and PFR configurations.",
      thumbnailPath: reactorDesignThumbnailPath,
      isNew: true,
    },
    {
      name: "Azeotrope Finder",
      path: "/toolbox/azeotrope-finder",
      description:
        "Predict and visualize azeotropic behavior of binary mixtures using various thermodynamic models.",
      thumbnailPath: azeotropeThumbnailPath,
      isNew: true,
    },
    {
      name: "Compound Properties",
      path: "/toolbox/compound-properties",
      description:
        "Fetch, plot, and compare various physical and thermodynamic properties of chemical compounds.", // Updated description
      thumbnailPath: compoundPropertiesThumbnailPath,
      isNew: true,
    },
    {
      name: "Residue Curve Map",
      path: "/toolbox/residue-curve-map",
      description:
        "Visualize and analyze residue curve maps for ternary mixtures, aiding in distillation sequence design.", // Updated description
      thumbnailPath: residueCurveMapThumbnailPath,
      isNew: true,
    },
    {
      name: "Reaction Kinetics",
      path: "/toolbox/kinetics",
      description:
        "Interactive simulator for chemical reaction kinetics. Model various reaction types and visualize concentration profiles over time.",
      thumbnailPath: kineticsThumbnailPath,
    },
    {
      name: "Process Dynamics",
      path: "/toolbox/process-dynamics",
      description:
        "Simulate process dynamics with various inputs and understand system behavior in chemical processes.",
      thumbnailPath: dynamicsThumbnailPath,
    },
  ];

  return (
    <TooltipProvider>
      <div className="container mx-auto">
      <div className="bio-container pt-8">
        <Card className="w-full bg-card text-card-foreground px-8">
          <CardHeader>
            <CardTitle className="text-2xl font-bold">Hi, I'm Victor</CardTitle>
            <CardDescription className="text-2xl font-bold">
              Chemical Engineering Optimization & Predictive Modeling Enthusiast
            </CardDescription>
          </CardHeader>

          <CardContent className="grid grid-cols-1 md:grid-cols-5 gap-6 bio-section">
            <div className="md:col-span-3 space-y-4">
              <div className="bio-paragraph flex items-center flex-wrap gap-2">
                <span>Languages and Frameworks:</span>
                <span className="lang-icons flex flex-wrap gap-3 items-center">
                  <Tooltip>
                    <TooltipTrigger>
                      <Image
                        src={pythonLogoPath}
                        alt="Python Logo"
                        width={25}
                        height={25}
                        title="Python"
                      />
                    </TooltipTrigger>
                    <TooltipContent>
                      <p>Python</p>
                    </TooltipContent>
                  </Tooltip>
                  <Tooltip>
                    <TooltipTrigger>
                      <Image
                        src={reactLogoPath}
                        alt="React Logo"
                        width={25}
                        height={25}
                        title="React"
                      />
                    </TooltipTrigger>
                    <TooltipContent>
                      <p>React</p>
                    </TooltipContent>
                  </Tooltip>
                  <Tooltip>
                    <TooltipTrigger>
                      <Image
                        src={nextjsLogoPath}
                        alt="Next.js Logo"
                        width={25}
                        height={25}
                        title="Next.js"
                      />
                    </TooltipTrigger>
                    <TooltipContent>
                      <p>Next.js</p>
                    </TooltipContent>
                  </Tooltip>
                  <Tooltip>
                    <TooltipTrigger>
                      <Image
                        src={javascriptLogoPath}
                        alt="JavaScript Logo"
                        width={25}
                        height={25}
                        title="JavaScript"
                      />
                    </TooltipTrigger>
                    <TooltipContent>
                      <p>JavaScript</p>
                    </TooltipContent>
                  </Tooltip>
                  <Tooltip>
                    <TooltipTrigger>
                      <Image
                        src={typescriptLogoPath}
                        alt="TypeScript Logo"
                        width={25}
                        height={25}
                        title="TypeScript"
                      />
                    </TooltipTrigger>
                    <TooltipContent>
                      <p>TypeScript</p>
                    </TooltipContent>
                  </Tooltip>
                  <Tooltip>
                    <TooltipTrigger>
                      <Image
                        src={html5LogoPath}
                        alt="HTML5 Logo"
                        width={30}
                        height={30}
                        title="HTML5"
                      />
                    </TooltipTrigger>
                    <TooltipContent>
                      <p>HTML5</p>
                    </TooltipContent>
                  </Tooltip>
                  <Tooltip>
                    <TooltipTrigger>
                      <Image
                        src={css3LogoPath}
                        alt="CSS3 Logo"
                        width={35}
                        height={35}
                        title="CSS3"
                      />
                    </TooltipTrigger>
                    <TooltipContent>
                      <p>CSS3</p>
                    </TooltipContent>
                  </Tooltip>
                  <Tooltip>
                    <TooltipTrigger>
                      <Image
                        src={dockerLogoPath}
                        alt="Docker Logo"
                        width={25}
                        height={25}
                        title="Docker"
                      />
                    </TooltipTrigger>
                    <TooltipContent>
                      <p>Docker</p>
                    </TooltipContent>
                  </Tooltip>
                  <Tooltip>
                    <TooltipTrigger>
                      <Image
                        src={supabaseLogoPath}
                        alt="Supabase Logo"
                        width={25} 
                        height={25}
                        title="Supabase"
                      />
                    </TooltipTrigger>
                    <TooltipContent>
                      <p>Supabase</p>
                    </TooltipContent>
                  </Tooltip>
                  <Tooltip>
                    <TooltipTrigger>
                      <Image
                        src={tailwindLogoPath}
                        alt="Tailwind CSS Logo"
                        width={25}
                        height={25}
                        title="Tailwind CSS"
                      />
                    </TooltipTrigger>
                    <TooltipContent>
                      <p>Tailwind CSS</p>
                    </TooltipContent>
                  </Tooltip>
                  <Tooltip>
                    <TooltipTrigger>
                      <Image
                        src={shadcnLogoPath}
                        alt="shadcn/ui Logo"
                        width={25}
                        height={25}
                        title="shadcn/ui"
                      />
                    </TooltipTrigger>
                    <TooltipContent>
                      <p>shadcn/ui</p>
                    </TooltipContent>
                  </Tooltip>
                </span>
              </div>

              <div className="bio-paragraph flex items-center flex-wrap gap-2">
                <span>Python Packages:</span>
                <span className="">
                  NumPy, SciPy, Pandas, Matplotlib, Scikit-learn, Plotly, RegEx,
                  Control, BeautifulSoup4
                </span>
                <span className="ml-1">📦</span>
              </div>

              <div className="bio-paragraph flex items-center gap-2">
                <p>
                  Currently a senior at the University of California, Santa
                  Barbara.
                </p>
                <Image
                  src={ucsbLogoPath}
                  alt="UCSB Logo"
                  width={25}
                  height={25}
                  className="ucsb-logo dark:invert-[0.8]"
                />
              </div>

              <p className="bio-paragraph">
                Will pursue a masters degree in materials science next year. 🎓
              </p>

              <div className="bio-paragraph flex items-baseline gap-2">
                <span>Want to contact me?</span>
                <a
                  href="mailto:victorliang@ucsb.edu"
                  className="text-primary hover:underline"
                >
                  victorliang@ucsb.edu
                </a>
                <span>📧</span>
              </div>

              <div className="bio-paragraph flex items-center gap-2">
                <p>Under heavy construction, but take a look around! 😄</p>
              </div>
            </div>

            <div className="md:col-span-2 flex items-center justify-center h-full">
              <div className="relative w-full h-full rounded-lg overflow-hidden">
                <Image
                  src={labIllustrationPath}
                  alt="Lab Illustration"
                  layout="fill"
                  objectFit="contain"
                />
              </div>
            </div>
          </CardContent>

          <CardFooter className="flex flex-col items-start pt-0">
            <div className="social-links flex gap-2">
              <Button
                variant="outline"
                size="icon"
                asChild
                className="bg-card hover:bg-accent"
              >
                <a
                  href="https://github.com/Victor-Liang-ChE"
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  <Image
                    src={githubLogoPath}
                    alt="GitHub Logo"
                    width={30}
                    height={30}
                    className="dark:invert"
                  />
                  <span className="sr-only">GitHub</span>
                </a>
              </Button>
              <Button
                variant="outline"
                size="icon"
                asChild
                className="bg-card hover:bg-accent"
              >
                <a
                  href="https://www.linkedin.com/in/victor-liang-567238231/"
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  <Image
                    src={linkedinLogoPath}
                    alt="LinkedIn Logo"
                    width={30}
                    height={30}
                  />
                  <span className="sr-only">LinkedIn</span>
                </a>
              </Button>
            </div>
          </CardFooter>
        </Card>
      </div>
      {/* Featured Simulations Section */}
      <div className="simulations-showcase mt-8">
        <h2 className="text-2xl font-bold mb-4 text-center">
          Featured Simulations and Tools
        </h2>
        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 w-full">
          {homeFeaturedSimulations.map((simulation) => (
            <Link href={simulation.path} key={simulation.name}>
              <Card className="bg-card text-card-foreground hover:shadow-lg transition-shadow cursor-pointer h-full flex flex-col relative overflow-hidden">
                {simulation.isUpdated && (
                  <div className="absolute top-3 right-[-32px] transform rotate-45 bg-yellow-400 text-black text-xs font-semibold py-2 px-8 shadow-lg z-10">
                    UPDATED
                  </div>
                )}
                {simulation.isNew && (
                  <div className="absolute top-3 right-[-28px] transform rotate-45 bg-orange-500 text-white text-xs font-semibold py-2 px-10 shadow-lg z-10">
                    NEW
                  </div>
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
      </div>
    </TooltipProvider>
  );
}