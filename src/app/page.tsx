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
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";

// Replace external URLs with local paths from public folder
const githubLogoPath = "/logos/github.svg";
const linkedinLogoPath = "/logos/linkedin.svg";
const ucsbLogoPath = "/logos/ucsb.svg";
const pythonLogoPath = "/logos/python.svg";
const javascriptLogoPath = "/logos/javascript.svg";
const html5LogoPath = "/logos/html5.svg";
const css3LogoPath = "/logos/css.svg";
const typescriptLogoPath = "/logos/typescript.svg";
const dockerLogoPath = "/logos/docker.svg";
const reactLogoPath = "/logos/react.svg";
const dashLogoPath = "/logos/dash.svg"; // You'll need to add this logo
const kineticsThumbnailPath = "/thumbnails/kinetics-thumbnail.png"; // Added path
const mccabeThumbnailPath = "/thumbnails/mccabe-thiele-thumbnail.png";
const dynamicsThumbnailPath = "/thumbnails/process-dynamics-thumbnail.png"; // Added path
const pidThumbnailPath = "/thumbnails/pid-tuning-thumbnail.png"; // Added path
const azeotropeThumbnailPath = "/thumbnails/azeotrope-finder-thumbnail.png"; // Added for Azeotrope Finder
const labIllustrationPath = "/images/lab-illustration.png"; // Added path for the illustration

export default function Home() {
  const homeFeaturedSimulations = [
    {
      name: "McCabe-Thiele",
      path: "/simulations/mccabe-thiele",
      description:
        "Select components and specify operating conditions to visualize distillation processes with accurate equilibrium diagrams.",
      thumbnailPath: mccabeThumbnailPath,
      isUpdated: true, // To show the "UPDATED" badge
    },
    {
      name: "Azeotrope Finder",
      path: "/simulations/azeotrope-finder",
      description:
        "Predict and visualize azeotropic behavior of binary mixtures using various thermodynamic models.",
      thumbnailPath: azeotropeThumbnailPath,
    },
    {
      name: "Reaction Kinetics",
      path: "/simulations/kinetics",
      description:
        "Interactive simulator for chemical reaction kinetics. Model various reaction types and visualize concentration profiles over time.",
      thumbnailPath: kineticsThumbnailPath,
    },
    {
      name: "Process Dynamics",
      path: "/simulations/process-dynamics",
      description:
        "Simulate process dynamics with various inputs and understand system behavior in chemical processes.",
      thumbnailPath: dynamicsThumbnailPath,
    },
  ];

  return (
    <div className="container mx-auto">
      <div className="bio-container pt-8">
        <Card className="w-full bg-card text-card-foreground px-8">
          <CardHeader>
            <CardTitle className="text-2xl font-bold">Hi, I'm Victor</CardTitle>
            <CardDescription className="text-2xl font-bold">
              Chemical Engineering Optimization & Predictive Modeling Specialist
            </CardDescription>
          </CardHeader>

          {/* Modified CardContent to use 5-column grid */}
          <CardContent className="grid grid-cols-1 md:grid-cols-5 gap-6 bio-section">
            {/* Left Column for Bio Text - Takes 3/5ths width */}
            <div className="md:col-span-3 space-y-4">
              <div className="bio-paragraph flex items-center flex-wrap gap-2">
                <span>Languages and Frameworks:</span>
                <span className="lang-icons flex flex-wrap gap-3 items-center">
                  <Image
                    src={pythonLogoPath}
                    alt="Python Logo"
                    width={25}
                    height={25}
                    title="Python"
                  />
                  <Image
                    src={reactLogoPath}
                    alt="React Logo"
                    width={25}
                    height={25}
                    title="React"
                  />
                  <Image
                    src={javascriptLogoPath}
                    alt="JavaScript Logo"
                    width={25}
                    height={25}
                    title="JavaScript"
                  />
                  <Image
                    src={typescriptLogoPath}
                    alt="TypeScript Logo"
                    width={25}
                    height={25}
                    title="TypeScript"
                  />
                  <Image
                    src={html5LogoPath}
                    alt="HTML5 Logo"
                    width={30}
                    height={30}
                    title="HTML5"
                  />
                  <Image
                    src={css3LogoPath}
                    alt="CSS3 Logo"
                    width={35}
                    height={35}
                    title="CSS3"
                  />
                  <Image
                    src={dockerLogoPath}
                    alt="Docker Logo"
                    width={25}
                    height={25}
                    title="Docker"
                  />
                </span>
              </div>

              <div className="bio-paragraph flex items-center flex-wrap gap-2">
                <span>Python Packages:</span>
                <span className="text-muted-foreground">
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

            {/* Right Column for Image - Takes 2/5ths width */}
            <div className="md:col-span-2 flex items-center justify-center h-full">
              {/* Image container fills its parent column's height and width */}
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
            {/* Remove justify-center and w-full to align left */}
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
          Featured Simulations
        </h2>
        {/* Use consistent grid and card structure */}
        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 w-full">
          {homeFeaturedSimulations.map((simulation) => (
            <Link href={simulation.path} key={simulation.name}>
              <Card className="bg-card text-card-foreground hover:shadow-lg transition-shadow cursor-pointer h-full flex flex-col relative overflow-hidden">
                {simulation.isUpdated && (
                  <div className="absolute top-3 right-[-32px] transform rotate-45 bg-red-600 text-white text-xs font-semibold py-2 px-8 shadow-lg z-10">
                    UPDATED
                  </div>
                )}
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
    </div>
  );
}