import type { Metadata } from 'next';

export const metadata: Metadata = {
  title: 'Canopy — Process Simulator',
  description: 'Chemical process simulator inspired by Aspen Plus. Build and solve flowsheets with rigorous thermodynamics.',
};

/** Canopy uses its own navbar — suppress the global layout navbar via a nested layout. */
export default function CanopyLayout({ children }: { children: React.ReactNode }) {
  return (
    <div className="flex flex-col h-screen overflow-hidden canopy-app">
      {children}
    </div>
  );
}
