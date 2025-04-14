import type { Metadata } from "next";
// Remove Inter font import if it's not the intended default
// import { Inter } from "next/font/google";
import "./globals.css";
import { ThemeProvider } from "@/components/theme-provider";
import Navbar from "@/components/Navbar";
import Script from 'next/script';

// Remove Inter font initialization if not the default
// const inter = Inter({ subsets: ["latin"] });

export const metadata: Metadata = {
  title: "Victor Liang",
  description: "Personal Website",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en" suppressHydrationWarning>
      <head>
        {/* Add MathQuill CSS */}
        <link
          rel="stylesheet"
          href="https://cdnjs.cloudflare.com/ajax/libs/mathquill/0.10.1/mathquill.css"
        />
      </head>
      {/* Remove inter.className to use default font from globals.css */}
      <body>
        <ThemeProvider
          attribute="class"
          defaultTheme="system"
          enableSystem
          disableTransitionOnChange
        >
          <Navbar />
          <main>{children}</main>
        </ThemeProvider>
        {/* Change strategy to 'afterInteractive' */}
        <Script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js" strategy="afterInteractive" />
        <Script src="https://cdnjs.cloudflare.com/ajax/libs/mathquill/0.10.1/mathquill.js" strategy="afterInteractive" />
      </body>
    </html>
  );
}
