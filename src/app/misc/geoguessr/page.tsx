// app/misc/geoguessr/page.tsx

'use client';

import { useState, useEffect, Suspense } from 'react';
import dynamic from 'next/dynamic';
import { LeafletMouseEvent, Layer } from 'leaflet';
import type { FeatureCollection, Feature } from 'geojson';

// Shadcn UI Components
import { Button } from '@/components/ui/button';
import {
  Sheet,
  SheetContent,
  SheetHeader,
  SheetTitle,
  SheetDescription,
} from '@/components/ui/sheet';
import {
  Drawer,
  DrawerContent,
  DrawerHeader,
  DrawerTitle,
  DrawerDescription,
  DrawerFooter,
  DrawerTrigger,
} from '@/components/ui/drawer';
import { Input } from '@/components/ui/input';
import { Card, CardContent, CardHeader } from '@/components/ui/card';
import { CardTitle } from '@/components/ui/card';

// Dynamically import react-leaflet components to prevent SSR errors
const MapContainer = dynamic(() => import('react-leaflet').then(mod => mod.MapContainer), { ssr: false });
const TileLayer = dynamic(() => import('react-leaflet').then(mod => mod.TileLayer), { ssr: false });
const GeoJSON = dynamic(() => import('react-leaflet').then(mod => mod.GeoJSON), { ssr: false });

// Mock database of hints. In a real app, this would come from a database.
const hintsDatabase = [
    {
        country: 'United States of America',
        id: 'usa-1',
        description: 'Yellow license plates are common in New York.',
        tags: ['license', 'plate', 'yellow', 'new', 'york', 'taxi'],
        imageUrl: 'https://via.placeholder.com/300x200.png?text=Yellow+Plate'
    },
    {
        country: 'Canada',
        id: 'can-1',
        description: 'Stop signs can be bilingual (English and French).',
        tags: ['stop', 'sign', 'red', 'street', 'french', 'english'],
        imageUrl: 'https://via.placeholder.com/300x200.png?text=Bilingual+Sign'
    }
];

export default function GeoGamePage() {
  // State variables declared once below

  const [isSheetOpen, setIsSheetOpen] = useState(false);

  const [geoJsonData, setGeoJsonData] = useState<FeatureCollection | null>(null);
  const [selectedCountry, setSelectedCountry] = useState<string | null>(null);
  const [countryHints, setCountryHints] = useState<any[]>([]);
  const [aiQuery, setAiQuery] = useState('');
  const [aiResults, setAiResults] = useState<any[]>([]);

  useEffect(() => {
    if (isSheetOpen && selectedCountry) {
      console.log('Sheet opened for country:', selectedCountry);
    } else if (!isSheetOpen) {
      console.log('Sheet closed, country cleared');
    }
  }, [isSheetOpen, selectedCountry]);

  useEffect(() => {
    if (selectedCountry) {
      console.log('Selected country updated to:', selectedCountry);
      // Optionally force a re-render or additional logic here if needed
    } else {
      console.log('Selected country cleared');
    }
  }, [selectedCountry]);

  useEffect(() => {
    fetch('https://raw.githubusercontent.com/datasets/geo-countries/master/data/countries.geojson')
      .then(response => response.json())
      .then(data => setGeoJsonData(data as FeatureCollection));
  }, []);

  useEffect(() => {
    if (selectedCountry) {
      console.log('Selected country updated to:', selectedCountry);
      // Optionally force a re-render or additional logic here if needed
    } else {
      console.log('Selected country cleared');
    }
  }, [selectedCountry]);
  
  const defaultStyle = {
    fillColor: "#3388ff",
    weight: 1,
    opacity: 1,
    color: '#333',
    fillOpacity: 0.2
  };

  const highlightStyle = {
    fillOpacity: 0.5,
    weight: 2
  };
  
  const onEachFeature = (feature: Feature, layer: Layer) => {
    layer.on({
      mouseover: (event: LeafletMouseEvent) => {
        event.target.setStyle(highlightStyle);
        event.target.bringToFront();
      },
      mouseout: (event: LeafletMouseEvent) => {
        event.target.setStyle(defaultStyle);
      },
      click: (event: LeafletMouseEvent) => {
        if (feature.properties) {
          const clickedCountry = feature.properties.ADMIN;
          if (isSheetOpen && selectedCountry === clickedCountry) {
            setIsSheetOpen(false);
          } else {
            setSelectedCountry(clickedCountry);
            setIsSheetOpen(true);
            const hints = hintsDatabase.filter(hint => hint.country === clickedCountry);
            setCountryHints(hints);
          }
        }
      },
    });
  };

  const handleAiSearch = () => {
    if (!aiQuery.trim()) {
        setAiResults([]);
        return;
    }
    const searchTerms = aiQuery.toLowerCase().split(' ');
    const results = hintsDatabase.filter(hint => 
      searchTerms.some(term => hint.tags.includes(term))
    );
    setAiResults(results);
  };

  if (!geoJsonData) {
    return <div className="h-screen w-screen flex items-center justify-center">Loading Country Data...</div>;
  }

  return (
    <Suspense fallback={<div>Loading UI...</div>}>
      <div className="relative h-screen w-screen">
        <MapContainer
          center={[30, 0]}
          zoom={3}
          style={{ height: '100%', width: '100%' }}
          scrollWheelZoom={true}
          worldCopyJump={true}
        >
          <TileLayer
            attribution='&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/attributions">CARTO</a>'
            url="https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png"
          />
          <GeoJSON 
              data={geoJsonData} 
              style={defaultStyle}
              onEachFeature={onEachFeature} 
          />
        </MapContainer>

        {/* Sidebar for Country Details */}
        <Sheet
          modal={false}
          open={isSheetOpen}
          onOpenChange={(open: boolean) => {
            setIsSheetOpen(open);
            if (!open) {
              setSelectedCountry(null);
            }
          }}
        >
          <SheetContent 
            className="w-[400px] sm:w-[540px] overflow-y-auto z-[1000] top-[85px] h-[calc(100%-64px)] [&>button]:hidden !border-none !shadow-none focus:outline-none focus:ring-0"
          >
            <SheetHeader>
              <SheetTitle>Hints for {selectedCountry}</SheetTitle>
              <SheetDescription>
                Here are some visual clues for identifying this country.
              </SheetDescription>
            </SheetHeader>
            {countryHints.length > 0 ? (
              <div className="grid gap-4 py-4">
                {countryHints.map((hint) => (
                  <Card key={hint.id}>
                    <CardHeader>
                      <CardTitle className="text-base font-medium">{hint.description}</CardTitle>
                    </CardHeader>
                    <CardContent>
                      <img src={hint.imageUrl} alt={hint.description} className="w-full rounded-md" />
                    </CardContent>
                  </Card>
                ))}
              </div>
            ) : (
              <div className="py-4 text-sm text-muted-foreground">
                No hints available for this country yet.
              </div>
            )}
          </SheetContent>
        </Sheet>

        {/* AI Assistant Drawer */}
        <div className="absolute bottom-4 right-4 z-[1000]">
          <Drawer>
              <DrawerTrigger asChild>
                  <Button size="lg" className="rounded-full h-16 w-16 shadow-lg">AI</Button>
              </DrawerTrigger>
              <DrawerContent>
                  <div className="mx-auto w-full max-w-sm">
                      <DrawerHeader>
                          <DrawerTitle>AI Geo Helper</DrawerTitle>
                          <DrawerDescription>Describe a clue and I'll find matching hints.</DrawerDescription>
                      </DrawerHeader>
                      <div className="p-4 pb-0">
                          <div className="flex items-center justify-center space-x-2">
                             <Input 
                                  id="ai-query" 
                                  value={aiQuery} 
                                  onChange={(e) => setAiQuery(e.target.value)}
                                  placeholder="e.g., yellow license plate"
                              />
                             <Button onClick={handleAiSearch}>Search</Button>
                          </div>
                          <div className="mt-4 max-h-48 overflow-y-auto">
                              {aiResults.map(hint => (
                                  <div key={hint.id} className="p-2 border rounded-md mb-2">
                                      <p className="font-bold">{hint.country}</p>
                                      <p>{hint.description}</p>
                                  </div>
                              ))}
                          </div>
                      </div>
                      <DrawerFooter>
                          <Button variant="outline" onClick={() => setAiResults([])}>Clear Results</Button>
                      </DrawerFooter>
                  </div>
              </DrawerContent>
          </Drawer>
        </div>
      </div>
    </Suspense>
  );
}