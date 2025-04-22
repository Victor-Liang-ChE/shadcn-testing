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

// Interfaces for menu items and data structure
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

// Interface for the combined data fetched from the API
interface AllMealsData {
    breakfast: MenuData;
    lunch: MenuData;
    dinner: MenuData;
    dateRange: string;
}

// Dropdown options
const mealOptions = [
    { value: 'breakfast', label: 'Breakfast' },
    { value: 'lunch', label: 'Lunch/Brunch' },
    { value: 'dinner', label: 'Dinner' }
];

export default function PortolaMenuPage() {
  // State for the currently displayed meal's data
  const [menuData, setMenuData] = useState<MenuData | null>(null);
  // State to store data for all meals fetched initially
  const [allMenuData, setAllMenuData] = useState<AllMealsData | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [mealType, setMealType] = useState('dinner'); // Default meal type
  const [pageTitle, setPageTitle] = useState('Loading Menu...');

  // Update page title based on state
  useEffect(() => {
    let titleText = 'Portola Dining Menu';
    const topLevelDateRange = allMenuData?.dateRange;

    if (isLoading) {
       titleText = 'Loading Menu...';
    } else if (error) {
       titleText = 'Error Loading Menu';
    // Display the formatted date range if it's available and not "N/A"
    } else if (topLevelDateRange && topLevelDateRange !== "N/A") {
        try {
            const [startDateStr, endDateStr] = topLevelDateRange.split(' to ');
            const currentYear = new Date().getFullYear();
            // Append the current year for parsing
            const startDate = new Date(`${startDateStr}/${currentYear}`);
            const endDate = new Date(`${endDateStr}/${currentYear}`);

            // Adjust year for end date if it wraps around New Year
            if (endDate < startDate) {
                endDate.setFullYear(currentYear + 1);
            }

            const options: Intl.DateTimeFormatOptions = { weekday: 'short' };
            const startDay = startDate.toLocaleDateString('en-US', options);
            const endDay = endDate.toLocaleDateString('en-US', options);

            // Format: Day MM/DD - Day MM/DD
            titleText = `${startDay} ${startDateStr} - ${endDay} ${endDateStr}`;
        } catch (e) {
            console.error("Error formatting date range:", e);
            // Fallback to original string if formatting fails
            titleText = topLevelDateRange;
        }
    } else {
        // Fallback title if data is missing or date range is unavailable
        const mealCapitalized = mealType.charAt(0).toUpperCase() + mealType.slice(1).replace('lunch', 'Lunch/Brunch');
        titleText = `Portola Dining ${mealCapitalized} Menu (Unavailable)`;
    }
    setPageTitle(titleText);
  // Depend on allMenuData, mealType, isLoading, error
  }, [allMenuData, mealType, isLoading, error]);

  // Fetch all menu data once on component mount
  useEffect(() => {
    const fetchAllMenus = async () => {
      console.log(`Fetching all menus on page load`);
      setIsLoading(true);
      setError(null);
      setAllMenuData(null); // Clear previous all data
      setMenuData(null); // Clear current display data

      try {
        // Fetch without meal parameter to get all data
        const response = await fetch(`/api/portola-menu`);
        if (!response.ok) {
          let errorMsg = `HTTP error! status: ${response.status}`;
          try { const errorData = await response.json(); errorMsg = errorData.error || errorMsg; } catch (e) {}
          throw new Error(errorMsg);
        }

        const data: AllMealsData = await response.json();
        console.log("All menu data received:", data);

        if (!data || !data.breakfast || !data.lunch || !data.dinner || typeof data.dateRange === 'undefined') {
          console.warn("Received unexpected data structure for all meals:", data);
          throw new Error('Invalid menu data structure received');
        }

        setAllMenuData(data);
        // Set the initially displayed menu based on the default mealType
        setMenuData(data[mealType as keyof AllMealsData] as MenuData);

      } catch (err: any) {
        console.error("Error fetching all menus:", err);
        setError(err.message || 'Error loading menu data. Please try again later.');
        setAllMenuData(null);
        setMenuData(null);
      } finally {
        setIsLoading(false);
      }
    };

    fetchAllMenus();
  }, []); // Empty dependency array: runs only once on mount

  // Handler for Shadcn Select change - Update displayed data from allMenuData
  const handleMealChange = (value: string) => {
    console.log("Changing meal type to:", value);
    setMealType(value);
    if (allMenuData) {
      // Select the data for the chosen meal from the stored allMenuData
      setMenuData(allMenuData[value as keyof AllMealsData] as MenuData);
    } else {
      // Should ideally not happen if initial fetch worked, but handle defensively
      setMenuData(null);
    }
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
    <div className="container mx-auto px-4 py-8 pt-12"> {/* Adjusted top padding */}
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
                {/* Increased number of skeleton lines from 8 to 12 */}
                {[...Array(24)].map((_, j) => (
                  <Skeleton key={j} className="h-4 w-full" />
                ))}
              </CardContent>
            </Card>
          ))}
        </div>
      ) : menuData && menuData.days.length > 0 ? (
        <div className="grid grid-cols-1 md:grid-cols-3 lg:grid-cols-5 gap-4">
          {menuData.days.map((day, index) => (
            <Card key={index} className="gap-2"> {/* Changed gap-6 to gap-2 */}
              <CardHeader className="pb-2"> {/* Adjust padding */}
                {/* Display full day title */}
                <CardTitle className="text-lg">{day.title}</CardTitle>
                {/* Add horizontal line */}
                <hr className="my-2 border-border" />
              </CardHeader>
              <CardContent className="pt-0">
                <ul className="space-y-2.5 text-sm">
                  {/* Check for the special brunch notice */}
                  {day.items.length === 1 && day.items[0].class === 'brunch-notice' ? (
                     <li className="italic text-muted-foreground">{day.items[0].text}</li>
                  ) : (
                     day.items.map((item, itemIndex) => {
                        let itemClasses = "";
                        if (item.class === 'course-row') {
                          itemClasses = "text-base font-medium mt-4 mb-2 list-none text-primary";
                        } else if (item.highlight) {
                          itemClasses = "font-semibold text-red-600 dark:text-red-400";
                        }

                        return (
                          <li key={itemIndex} className={itemClasses}>
                            {item.text}
                          </li>
                        );
                     })
                  )}
                </ul>
              </CardContent>
            </Card>
          ))}
        </div>
      ) : (
        <div className="text-center p-8 border rounded-md bg-muted">
          <p className="text-xl text-muted-foreground">No menu available for {mealType} this week.</p>
          {allMenuData?.dateRange && allMenuData.dateRange !== "N/A" && <p className="text-sm text-muted-foreground mt-1">Week: {allMenuData.dateRange}</p>}
        </div>
      )}
    </div>
  );
}
