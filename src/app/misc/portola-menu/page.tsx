'use client';

import { useState, useEffect } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Skeleton } from "@/components/ui/skeleton";

// We need to create the Tabs components, let's add this first

interface MenuItem {
  text: string;
  highlight?: boolean;
  class?: string;
}

interface MenuData {
  days: {
    title: string;
    items: MenuItem[];
  }[];
  dateRange: string;
}

// Custom useInterval hook with proper typing
function useInterval(callback: () => void, delay: number) {
  const savedCallback = useState(() => callback)[0];
  
  useEffect(() => {
    function tick() {
      savedCallback();
    }
    if (delay !== null) {
      const id = setInterval(tick, delay);
      return () => clearInterval(id);
    }
  }, [delay, savedCallback]);
}

export default function PortolaMenuPage() {
  const [menuData, setMenuData] = useState<MenuData | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [mealType, setMealType] = useState('dinner');

  const fetchMenu = async (meal: string) => {
    setIsLoading(true);
    try {
      const response = await fetch(`/api/portola-menu?meal=${meal}`);
      if (!response.ok) {
        throw new Error('Failed to fetch menu data');
      }
      const data = await response.json();
      setMenuData(data);
      setError(null);
    } catch (err) {
      setError('Error loading menu data. Please try again later.');
      console.error(err);
    } finally {
      setIsLoading(false);
    }
  };

  // Initial fetch
  useEffect(() => {
    fetchMenu(mealType);
  }, [mealType]);

  // Refresh every minute
  useInterval(() => {
    fetchMenu(mealType);
  }, 60000);

  if (error) {
    return (
      <div className="container mx-auto px-4 py-8">
        <h1 className="text-3xl font-bold mb-6">Portola Dining Menu</h1>
        <div className="bg-destructive/20 text-destructive border border-destructive/50 rounded-md p-4">
          {error}
        </div>
      </div>
    );
  }

  return (
    <div className="container mx-auto px-4 py-8">
      <div className="flex justify-between items-center mb-6 flex-col sm:flex-row gap-4">
        <h1 className="text-3xl font-bold">
          {menuData ? `Portola Dining ${mealType.charAt(0).toUpperCase() + mealType.slice(1)} Menu` : 'Portola Dining Menu'}
          {menuData && <span className="text-sm font-normal ml-2 block sm:inline">from {menuData.dateRange}</span>}
        </h1>
        
        <Tabs defaultValue={mealType} onValueChange={setMealType} className="w-full sm:w-auto">
          <TabsList>
            <TabsTrigger value="breakfast">Breakfast</TabsTrigger>
            <TabsTrigger value="lunch">Lunch/Brunch</TabsTrigger>
            <TabsTrigger value="dinner">Dinner</TabsTrigger>
          </TabsList>
        </Tabs>
      </div>

      {isLoading ? (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-5 gap-4">
          {[...Array(5)].map((_, i) => (
            <Card key={i}>
              <CardHeader>
                <Skeleton className="h-6 w-32" />
              </CardHeader>
              <CardContent>
                {[...Array(8)].map((_, j) => (
                  <Skeleton key={j} className="h-4 w-full my-2" />
                ))}
              </CardContent>
            </Card>
          ))}
        </div>
      ) : menuData ? (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-5 gap-4">
          {menuData.days.map((day, index) => (
            <Card key={index}>
              <CardHeader className="pb-2">
                <CardTitle className="text-lg">{day.title}</CardTitle>
              </CardHeader>
              <CardContent className="pt-2">
                <ul className="space-y-1">
                  {day.items.map((item, itemIndex) => {
                    if (item.class === 'course-row') {
                      return (
                        <li key={itemIndex} className="text-lg font-medium mt-2">
                          {item.text}
                        </li>
                      );
                    } else if (item.highlight) {
                      return (
                        <li key={itemIndex} className="font-bold text-red-600">
                          {item.text}
                        </li>
                      );
                    } else {
                      return <li key={itemIndex}>{item.text}</li>;
                    }
                  })}
                </ul>
              </CardContent>
            </Card>
          ))}
        </div>
      ) : (
        <div className="text-center p-8">
          <p className="text-xl text-destructive">No menu available for this week</p>
        </div>
      )}
    </div>
  );
}
