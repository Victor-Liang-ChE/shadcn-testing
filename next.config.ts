import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  /* config options here */
};

module.exports = {
  eslint: {
    // Warning: This allows production builds to successfully complete even if there are ESLint errors.
    ignoreDuringBuilds: true,
  },
  images: {
    domains: [
      "web.archive.org",
      'ufc.com',
    ],
  },
};

export default nextConfig;
