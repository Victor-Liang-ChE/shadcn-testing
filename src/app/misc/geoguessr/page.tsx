// app/misc/geoguessr/page.tsx

'use client'

import { useState, useEffect, Suspense, useRef } from 'react'
import type { LeafletMouseEvent, Layer } from 'leaflet'
import dynamic from 'next/dynamic'
import L from 'leaflet'
import type { FeatureCollection, Feature, Position } from 'geojson'

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

// Add this Set below your imports and mock database
const selectableCountries = new Set([
  // Europe - Western Europe
  'Andorra', 'Austria', 'Belgium', 'France', 'Germany', 'Greece', 'Ireland', 'Isle of Man', 'Italy', 'Luxembourg', 'Malta', 'Monaco', 'Netherlands', 'Portugal', 'Spain', 'Switzerland', 'United Kingdom',
  // Europe - Eastern Europe
  'Albania', 'Bulgaria', 'Croatia', 'Czechia', 'Hungary', 'Montenegro', 'North Macedonia', 'Poland', 'Romania', 'Russia', 'Republic of Serbia', 'Slovakia', 'Slovenia', 'Ukraine',
  // Europe - Nordics
  'Denmark', 'Faroe Islands', 'Finland', 'Greenland', 'Iceland', 'Norway', 'Sweden',
  // Europe - Baltics
  'Estonia', 'Latvia', 'Lithuania',
  // Americas - Latin America
  'Argentina', 'Bolivia', 'Brazil', 'Chile', 'Colombia', 'Costa Rica', 'Curaçao', 'Dominican Republic', 'Ecuador', 'Guatemala', 'Mexico', 'Panama', 'Peru', 'Puerto Rico', 'United States Virgin Islands', 'Uruguay',
  // Americas - North America
  'Bermuda', 'Canada', 'United States of America',
  // Asia - South & South-East Asia
  'Bangladesh', 'Bhutan', 'Cambodia', 'Indian Ocean Territories', 'India', 'Indonesia', 'Laos', 'Malaysia', 'Pakistan', 'Philippines', 'Singapore', 'Sri Lanka', 'Thailand', 'Vietnam',
  // Asia - Rest of Asia
  'China', 'Hong Kong S.A.R.', 'Japan', 'Kazakhstan', 'Kyrgyzstan', 'Mongolia', 'South Korea', 'Taiwan',
  // Asia - Middle East
  'Israel', 'Jordan', 'Palestine', 'Qatar', 'Tunisia', 'Turkey', 'United Arab Emirates',
  // Oceania
  'American Samoa', 'Australia', 'Guam', 'New Zealand', 'Northern Mariana Islands', 'United States Minor Outlying Islands',
  // Africa
  'Botswana', 'eSwatini', 'Ghana', 'Kenya', 'Lesotho', 'Madagascar', 'Nigeria', 'Rwanda', 'Senegal', 'South Africa', 'Uganda'
]);

// Define the bounding boxes for special territories
const specialRegions = [
  {
    name: 'Réunion',
    parent: 'France', // The name of the feature this territory belongs to
    bounds: {
      south: -21.44284571433252,
      north: -20.792980949405532,
      west: 55.16855111115359,
      east: 55.904345147111655,
    }
  },
  {
    name: 'Christmas Island',
    parent: 'Indian Ocean Territories', 
    bounds: {
      south: -10.663002117116974,
      north: -10.311947147775243,
      west: 105.43666677536476,
      east: 105.84598868861917,
    }
  }
];

// Add excludedRegions array for bounding boxes of regions to exclude from special handling.
const excludedRegions = [
  {
    // Cocos (Keeling) Islands
    parent: 'Indian Ocean Territories',
    bounds: {
      south: -12.281760044842969,
      north: -12.012616911899613,
      west: 96.71744895731854,
      east: 97.0129178271959,
    }
  },
  {
    // French Guiana
    parent: 'France',
    bounds: {
      south: 1.3461096780270339,
      north: 6.636884404557605,
      west: -55.25691176708574,
      east: -50.721209738287996,
    }
  }
];

// This helper function can go outside your component, after the imports
function findPolygonForRegion(feature: Feature, region: typeof specialRegions[0]) {
  if (feature.geometry.type !== 'MultiPolygon') return null;

  for (const polygonCoords of feature.geometry.coordinates) {
    let minLng = 180, maxLng = -180, minLat = 90, maxLat = -90;
    // The first element of polygonCoords is the outer ring of the polygon
    polygonCoords[0].forEach((point: Position) => {
      if (typeof point[0] === 'number' && typeof point[1] === 'number') {
        minLng = Math.min(minLng, point[0]);
        maxLng = Math.max(maxLng, point[0]);
        minLat = Math.min(minLat, point[1]);
        maxLat = Math.max(maxLat, point[1]);
      }
    });
    const centerLng = (minLng + maxLng) / 2;
    const centerLat = (minLat + maxLat) / 2;

    if (
      centerLat > region.bounds.south && centerLat < region.bounds.north &&
      centerLng > region.bounds.west && centerLng < region.bounds.east
    ) {
      return polygonCoords; // Return the coordinates of the matching polygon
    }
  }
  return null;
}

export default function GeoGamePage() {
  // AI search handler for Drawer
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
  // State variables declared once below

  const [isSheetOpen, setIsSheetOpen] = useState(false);

  const [geoJsonData, setGeoJsonData] = useState<FeatureCollection | null>(null);
  const [selectedCountry, setSelectedCountry] = useState<string | null>(null);
  const [countryHints, setCountryHints] = useState<any[]>([]);
  const [aiQuery, setAiQuery] = useState('');
  const [aiResults, setAiResults] = useState<any[]>([]);

  // Add this ref to store the temporary highlight layer
  const highlightLayerRef = useRef<L.Layer | null>(null);

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
      .then(data => {
        // --- Debugger Logic Starts Here ---

        // 1. Get all country names from the map data using the correct 'name' property.
        const namesFromData = new Set(
          data.features.map((feature: any) => feature.properties.name).filter((name: string) => name)
        );

        // 2. Compare your list against the map data to find mismatches.
        const mismatchedNames = [];
        for (const countryInYourList of selectableCountries) {
          if (!namesFromData.has(countryInYourList)) {
            mismatchedNames.push(countryInYourList);
          }
        }

        // 3. Log the definitive list of names you need to fix.
        if (mismatchedNames.length > 0) {
          console.log("MISMATCHED COUNTRIES:", mismatchedNames.sort());
          console.log("These names from your list were not found in the map data. You need to correct them.");
        } else {
          console.log("All country names in your list match the map data! ✅");
        }
        
        // --- Debugger Logic Ends ---

        setGeoJsonData(data as FeatureCollection);
      });
  }, []);

  useEffect(() => {
    if (selectedCountry) {
      console.log('Selected country updated to:', selectedCountry);
      // Optionally force a re-render or additional logic here if needed
    } else {
      console.log('Selected country cleared');
    }
  }, [selectedCountry]);
  
  const selectableStyle = {
    fillColor: "#3388ff",
    weight: 1,
    opacity: 1,
    color: '#333',
    fillOpacity: 0.2
  };

  const disabledStyle = {
    fillColor: "#666",
    weight: 1,
    opacity: 0.5,
    color: '#888',
    fillOpacity: 0.4,
    interactive: false // This disables all mouse events
  };

  const highlightStyle = {
    fillOpacity: 0.5,
    weight: 2
  };
  
  // This function decides the style for each country
  const styleFeature = (feature?: Feature) => {
    if (feature && selectableCountries.has(feature.properties!.name)) {
      return selectableStyle;
    }
    return disabledStyle;
  };

  const onEachFeature = (feature: Feature, layer: Layer) => {
    // This function now correctly applies the initial style to all parts of the feature
    const styleFeatureLayer = () => {
      const parentFeatureName = feature.properties!.name;
      if (typeof (layer as any).eachLayer === 'function') {
        (layer as any).eachLayer((polygonLayer: any) => {
          const layerCenter = polygonLayer.getBounds().getCenter();
          let isExcluded = false;
          for (const region of excludedRegions) {
            if (
              parentFeatureName === region.parent &&
              layerCenter.lat > region.bounds.south && layerCenter.lat < region.bounds.north &&
              layerCenter.lng > region.bounds.west && layerCenter.lng < region.bounds.east
            ) {
              isExcluded = true;
              break;
            }
          }
          polygonLayer.setStyle(isExcluded ? disabledStyle : selectableStyle);
        });
      } else {
        // Fallback for simple polygons
        (layer as L.Path).setStyle(styleFeature(feature));
      }
    };

    styleFeatureLayer(); // Apply the initial correct style

    layer.on({
      mouseover: (event: LeafletMouseEvent) => {
        const parentFeatureName = feature.properties!.name;
        const hoveredLat = event.latlng.lat;
        const hoveredLng = event.latlng.lng;

        // Check if the cursor is over an excluded region first.
        let isOverExcludedRegion = false;
        for (const region of excludedRegions) {
          if (
            parentFeatureName === region.parent &&
            hoveredLat > region.bounds.south && hoveredLat < region.bounds.north &&
            hoveredLng > region.bounds.west && hoveredLng < region.bounds.east
          ) {
            isOverExcludedRegion = true;
            break;
          }
        }

        // If it's an excluded region, simply do nothing and exit.
        if (isOverExcludedRegion) {
          return;
        }

        // --- The rest of your mouseover logic for valid areas ---
        if (highlightLayerRef.current) {
          highlightLayerRef.current.remove();
          highlightLayerRef.current = null;
        }

        let specialRegionHovered = null;
        for (const region of specialRegions) {
          if (
            parentFeatureName === region.parent &&
            hoveredLat > region.bounds.south && hoveredLat < region.bounds.north &&
            hoveredLng > region.bounds.west && hoveredLng < region.bounds.east
          ) {
            specialRegionHovered = region;
            break;
          }
        }

        if (specialRegionHovered) {
          const polygonToHighlight = findPolygonForRegion(feature, specialRegionHovered);
          if (polygonToHighlight) {
            import('leaflet').then(L => {
              const temporaryHighlightStyle = { ...highlightStyle, interactive: false };
              const highlight = L.geoJSON({
                type: 'Feature',
                properties: {},
                geometry: { type: 'Polygon', coordinates: polygonToHighlight },
              } as import('geojson').Feature, { style: temporaryHighlightStyle });
              highlight.addTo((event.target as any)._map);
              highlightLayerRef.current = highlight;
            });
          }
        } else {
          (event.target as L.Path).setStyle(highlightStyle);
        }
        (event.target as L.Path).bringToFront();
      },
      mouseout: (_event: LeafletMouseEvent) => {
        if (highlightLayerRef.current) {
          highlightLayerRef.current.remove();
          highlightLayerRef.current = null;
        }
        // Instead of re-applying a single style, re-run the logic that styles each part
        styleFeatureLayer();
      },
      click: (event: LeafletMouseEvent) => {
        const clickedLat = event.latlng.lat;
        const clickedLng = event.latlng.lng;
        const parentFeatureName = feature.properties!.name;

        // *** START FIX ***
        // Repeat the exclusion check for the click event.
        let isOverExcludedRegion = false;
        for (const region of excludedRegions) {
          if (
            parentFeatureName === region.parent &&
            clickedLat > region.bounds.south && clickedLat < region.bounds.north &&
            clickedLng > region.bounds.west && clickedLng < region.bounds.east
          ) {
            isOverExcludedRegion = true;
            break;
          }
        }

        if (isOverExcludedRegion) {
          return; // Do not proceed with the click action.
        }
        // *** END FIX ***
        
        let finalCountryName = parentFeatureName;
        for (const region of specialRegions) {
          if (
            parentFeatureName === region.parent &&
            clickedLat > region.bounds.south && clickedLat < region.bounds.north &&
            clickedLng > region.bounds.west && clickedLng < region.bounds.east
          ) {
            finalCountryName = region.name;
            break;
          }
        }

        if (isSheetOpen && selectedCountry === finalCountryName) {
          setIsSheetOpen(false);
        } else {
          setSelectedCountry(finalCountryName);
          setIsSheetOpen(true);
          const hints = hintsDatabase.filter(hint => hint.country === finalCountryName);
          setCountryHints(hints);
        }
      },
    });
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
              style={styleFeature}
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