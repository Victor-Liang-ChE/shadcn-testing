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

export default function Home() {
  return (
    <div className="container mx-auto px-2 py-8 flex flex-col items-center gap-8">
      <Card className="w-full max-w-5xl bg-card text-card-foreground">
        <CardHeader>
          <CardTitle className="text-2xl font-bold text-card-foreground/90 pt-1">Hi, I'm Victor Liang</CardTitle>
          <CardDescription className="text-2xl text-card-foreground/90 pt-1">
            Chemical Engineering Optimization & Predictive Modeling Specialist
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-0 equal-spacing">
          <div>
            <h3 className="flex items-center gap-4 flex-wrap">
              Languages, Frameworks, and Other Tools:
              <span className="flex flex-wrap gap-3 items-center ml-0">
                <Image src={pythonLogoPath} alt="Python Logo" width={24} height={24} title="Python" />
                <Image src={reactLogoPath} alt="React Logo" width={24} height={24} title="React" />
                <Image src={javascriptLogoPath} alt="JavaScript Logo" width={24} height={24} title="JavaScript" />
                <Image src={typescriptLogoPath} alt="TypeScript Logo" width={24} height={24} title="TypeScript" />
                <Image src={html5LogoPath} alt="HTML5 Logo" width={24} height={24} title="HTML5" />
                <Image src={css3LogoPath} alt="CSS3 Logo" width={24} height={24} title="CSS3" />
                <Image src={dockerLogoPath} alt="Docker Logo" width={24} height={24} title="Docker" />
              </span>
            </h3>
          </div>
          <div>
            <h3 className="flex items-center gap-4 flex-wrap">
              Python Packages:
              <span className="flex flex-wrap gap-2 ml-0">
                <Badge variant="secondary">NumPy</Badge>
                <Badge variant="secondary">SciPy</Badge>
                <Badge variant="secondary">Pandas</Badge>
                <Badge variant="secondary">Matplotlib</Badge>
                <Badge variant="secondary">Scikit-learn</Badge>
                <Badge variant="secondary">Plotly</Badge>
                <Badge variant="secondary">RegEx</Badge>
                <Badge variant="secondary">Control</Badge>
                <Badge variant="secondary">BeautifulSoup4</Badge>
              </span>
            </h3>
          </div>
          <div className="flex items-center gap-2">
            <p>Currently a senior at the University of California, Santa Barbara.</p>
            <Image src={ucsbLogoPath} alt="UCSB Logo" width={40} height={40} className="dark:invert-[0.8]" />
          </div>
          <p>Will pursue a masters degree in materials science next year.</p>
          <div className="flex items-baseline gap-2">
            <h3 className="font-semibold">Want to contact me?</h3>
            <a href="mailto:victorliang@ucsb.edu" className="text-primary hover:underline">
              victorliang@ucsb.edu
            </a>
          </div>
        </CardContent>
        <CardFooter className="flex flex-col items-center pt-6">
          <p className="text-sm text-muted-foreground mb-4">Under heavy construction, but take a look around!</p>
          <div className="flex gap-2 justify-center">
            <Button variant="outline" size="icon" asChild className="bg-card hover:bg-accent">
              <a href="https://github.com/Victor-Liang-ChE" target="_blank" rel="noopener noreferrer">
                <Image src={githubLogoPath} alt="GitHub Logo" width={20} height={20} className="dark:invert" />
                <span className="sr-only">GitHub</span>
              </a>
            </Button>
            <Button variant="outline" size="icon" asChild className="bg-card hover:bg-accent">
              <a href="https://www.linkedin.com/in/victor-liang-567238231/" target="_blank" rel="noopener noreferrer">
                <Image src={linkedinLogoPath} alt="LinkedIn Logo" width={20} height={20} />
                <span className="sr-only">LinkedIn</span>
              </a>
            </Button>
          </div>
        </CardFooter>
      </Card>

      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 w-full max-w-5xl">
        <Link href="/mccabe" passHref>
          <Card className="bg-card text-card-foreground hover:shadow-lg transition-shadow cursor-pointer h-full flex flex-col">
            <CardHeader>
              <div className="bg-muted h-32 w-full rounded-md mb-4 flex items-center justify-center text-muted-foreground">Thumbnail</div>
              <CardTitle>McCabe-Thiele</CardTitle>
            </CardHeader>
            <CardContent className="flex-grow">
              <CardDescription>Interactive McCabe-Thiele diagrams for binary distillation.</CardDescription>
            </CardContent>
          </Card>
        </Link>

        <Link href="/pid-tuning" passHref>
          <Card className="bg-card text-card-foreground hover:shadow-lg transition-shadow cursor-pointer h-full flex flex-col">
            <CardHeader>
              <div className="bg-muted h-32 w-full rounded-md mb-4 flex items-center justify-center text-muted-foreground">Thumbnail</div>
              <CardTitle>PID Tuning</CardTitle>
            </CardHeader>
            <CardContent className="flex-grow">
              <CardDescription>Tools and simulations for PID controller tuning methods.</CardDescription>
            </CardContent>
          </Card>
        </Link>

        <Link href="/kinetics" passHref>
          <Card className="bg-card text-card-foreground hover:shadow-lg transition-shadow cursor-pointer h-full flex flex-col">
            <CardHeader>
              <div className="bg-muted h-32 w-full rounded-md mb-4 flex items-center justify-center text-muted-foreground">Thumbnail</div>
              <CardTitle>Reaction Kinetics</CardTitle>
            </CardHeader>
            <CardContent className="flex-grow">
              <CardDescription>Simulate and analyze various chemical reaction kinetics.</CardDescription>
            </CardContent>
          </Card>
        </Link>

        <Link href="/process-control" passHref>
          <Card className="bg-card text-card-foreground hover:shadow-lg transition-shadow cursor-pointer h-full flex flex-col">
            <CardHeader>
              <div className="bg-muted h-32 w-full rounded-md mb-4 flex items-center justify-center text-muted-foreground">Thumbnail</div>
              <CardTitle>Process Control</CardTitle>
            </CardHeader>
            <CardContent className="flex-grow">
              <CardDescription>Explore concepts and simulations in process control systems.</CardDescription>
            </CardContent>
          </Card>
        </Link>
      </div>
    </div>
  );
}
