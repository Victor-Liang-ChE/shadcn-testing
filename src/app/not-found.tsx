import Link from "next/link";

export default function NotFound() {
  return (
    <div className="flex flex-col items-center justify-center min-h-screen bg-card text-card-foreground p-8">
      <h1 className="text-5xl font-bold mb-4">404</h1>
      <p className="text-xl mb-6 text-center max-w-xl">
        The heat is lost, the insulation is breached, your target destination cannot be reached.
      </p>
      <Link href="/">
        <span className="text-primary underline hover:text-blue-500 transition">Return Home</span>
      </Link>
    </div>
  );
}
