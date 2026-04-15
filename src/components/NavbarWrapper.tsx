'use client';

import { usePathname } from 'next/navigation';
import Navbar from './Navbar';

/** Wraps the global Navbar and hides it on routes that provide their own (e.g. /canopy). */
export default function NavbarWrapper() {
  const pathname = usePathname();
  if (pathname.startsWith('/canopy')) return null;
  return <Navbar />;
}
