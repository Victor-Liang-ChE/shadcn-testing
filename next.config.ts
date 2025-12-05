import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  /* config options here */
};

module.exports = {
  images: {
    remotePatterns: [
      {
        protocol: 'https',
        hostname: 'web.archive.org',
      },
      {
        protocol: 'https',
        hostname: 'ufc.com',
      },
    ],
  },
};

export default nextConfig;
