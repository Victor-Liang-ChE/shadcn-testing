'use client';

import React, { useState, useEffect } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Skeleton } from "@/components/ui/skeleton";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";

// Custom useInterval hook with proper typing
function useInterval(callback: () => void, delay: number | null) {
  // Provide an initial value (undefined) to useRef
  const savedCallback = React.useRef<(() => void) | undefined>(undefined);

  // Remember the latest callback.
  useEffect(() => {
    savedCallback.current = callback;
  }, [callback]);

  // Set up the interval.
  useEffect(() => {
    function tick() {
      // Check if savedCallback.current exists before calling
      if (savedCallback.current) {
        savedCallback.current(); // Correctly call the stored callback
      }
    }
    if (delay !== null) {
      const id = setInterval(tick, delay);
      return () => clearInterval(id);
    }
  }, [delay]);
}

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

// Dropdown options
const mealOptions = [
    { value: 'breakfast', label: 'Breakfast' },
    { value: 'lunch', label: 'Lunch/Brunch' },
    { value: 'dinner', label: 'Dinner' }
];

export default function PortolaMenuPage() {
  const [menuData, setMenuData] = useState<MenuData | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [mealType, setMealType] = useState('dinner');
  const [pageTitle, setPageTitle] = useState('Loading Menu...'); // State for title

  // Update page title based on state
  useEffect(() => {
    let titleText = 'Portola Dining Menu';
    if (isLoading) {
       titleText = 'Loading Menu...';
    } else if (error) {
       titleText = 'Error Loading Menu';
    } else if (menuData) {
        const mealCapitalized = mealType.charAt(0).toUpperCase() + mealType.slice(1);
        if (menuData.days.length > 0 && menuData.dateRange && menuData.dateRange !== "N/A") {
            titleText = `Portola Dining ${mealCapitalized} Menu (${menuData.dateRange})`;
        } else {
            titleText = `Portola Dining ${mealCapitalized} Menu (Unavailable)`;
        }
    }
    setPageTitle(titleText);
    // If you have a global title update mechanism, call it here
    // document.title = titleText; // Example for browser tab title
  }, [menuData, mealType, isLoading, error]);

  const fetchMenu = async (meal: string) => {
    console.log(`Fetching menu for: ${meal}`);
    setIsLoading(true);
    setError(null);
    // Keep previous data while loading for smoother transition? Optional.
    // setMenuData(null);

    try {
      // Use the local API endpoint
      const response = await fetch(`/api/portola-menu?meal=${meal}`);
      if (!response.ok) {
        let errorMsg = `HTTP error! status: ${response.status}`;
        try {
            const errorData = await response.json();
            errorMsg = errorData.error || errorMsg;
        } catch (e) { /* Ignore json parsing error if response not json */ }
        throw new Error(errorMsg);
      }

      const data: MenuData = await response.json();
      console.log("Menu data received:", data);

      if (!data || !data.days || typeof data.dateRange === 'undefined') {
        console.warn("Received unexpected data structure:", data);
        throw new Error('Invalid menu data structure received');
      }

      setMenuData(data);
    } catch (err: any) {
      console.error("Error fetching menu:", err);
      setError(err.message || 'Error loading menu data. Please try again later.');
      setMenuData(null); // Ensure data is null on error
    } finally {
      setIsLoading(false);
    }
  };

  // Initial fetch and fetch on mealType change
  useEffect(() => {
    fetchMenu(mealType);
  }, [mealType]);

  // Auto-refresh every 60 seconds
  useInterval(() => {
    fetchMenu(mealType);
  }, 60000); // Pass delay directly

  // Handler for Shadcn Select change
  const handleMealChange = (value: string) => {
    console.log("Changing meal type to:", value);
    setMealType(value);
  };

  if (error && !isLoading) { // Show error only if not loading
    return (
      <div className="container mx-auto px-4 py-8 pt-12">
        <h1 className="text-2xl font-bold mb-4 text-center">{pageTitle}</h1>
        <div className="bg-destructive/10 border border-destructive text-destructive px-4 py-3 rounded text-center">
          {error}
        </div>
      </div>
    );
  }

  return (
    <div className="container mx-auto px-4 py-8 pt-12">
      {/* Title and Select Dropdown */}
      <div className="flex flex-col sm:flex-row justify-between items-center mb-6 gap-4">
        <h1 className="text-xl md:text-2xl font-bold text-center sm:text-left flex-grow">
          {pageTitle}
        </h1>
        <Select value={mealType} onValueChange={handleMealChange}>
          <SelectTrigger className="w-full sm:w-[180px]">
            <SelectValue placeholder="Select Meal" />
          </SelectTrigger>
          <SelectContent>
            {mealOptions.map(option => (
              <SelectItem key={option.value} value={option.value}>
                {option.label}
              </SelectItem>
            ))}
          </SelectContent>
        </Select>
      </div>

      {/* Menu Grid or Loading Skeletons */}
      {isLoading ? (
        <div className="grid grid-cols-1 md:grid-cols-3 lg:grid-cols-5 gap-4">
          {[...Array(5)].map((_, i) => (
            <Card key={i}>
              <CardHeader>
                <Skeleton className="h-6 w-3/4" />
              </CardHeader>
              <CardContent className="space-y-3 pt-2">
                {[...Array(8)].map((_, j) => (
                  <Skeleton key={j} className="h-4 w-full" />
                ))}
              </CardContent>
            </Card>
          ))}
        </div>
      ) : menuData && menuData.days.length > 0 ? (
        <div className="grid grid-cols-1 md:grid-cols-3 lg:grid-cols-5 gap-4">
          {menuData.days.map((day, index) => (
            <Card key={index}>
              <CardHeader>
                <CardTitle className="text-lg">{day.title}</CardTitle>
              </CardHeader>
              <CardContent className="pt-0"> {/* Adjust padding */}
                <ul className="space-y-2.5 text-sm"> {/* Keep increased spacing */}
                  {day.items.map((item, itemIndex) => {
                    let itemClasses = "";
                    if (item.class === 'course-row') {
                      // Use Tailwind for styling course rows
                      itemClasses = "text-base font-medium mt-4 mb-2 list-none text-primary"; // Example styling
                    } else if (item.highlight) {
                      // Use Tailwind for highlighting
                      itemClasses = "font-semibold text-red-600 dark:text-red-400";
                    }

                    return (
                      <li key={itemIndex} className={itemClasses}>
                        {item.text}
                      </li>
                    );
                  })}
                </ul>
              </CardContent>
            </Card>
          ))}
        </div>
      ) : (
        <div className="text-center p-8 border rounded-md bg-muted">
          <p className="text-xl text-muted-foreground">No menu available for {mealType} this week.</p>
          {menuData?.dateRange && menuData.dateRange !== "N/A" && <p className="text-sm text-muted-foreground mt-1">Week: {menuData.dateRange}</p>}
        </div>
      )}
    </div>
  );
}
