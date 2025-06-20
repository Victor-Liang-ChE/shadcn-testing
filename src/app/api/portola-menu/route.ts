import { NextResponse } from 'next/server';
import * as cheerio from 'cheerio';

// Interfaces remain the same as they match the frontend's needs
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
}

interface AllMealsData {
    breakfast: MenuData;
    lunch: MenuData;
    dinner: MenuData;
    dateRange: string; // We'll leave this for compatibility, but it's less important now
}

// The sections that trigger a highlight on their first item.
const SECTIONS_TO_HIGHLIGHT = new Set([
    "The Brick",
    "Chef's Choice",
    "International",
]);

// Headers to exclude from the menu display
const HEADERS_TO_EXCLUDE = new Set([
    "Greens & Grains",
    "The Brick",
    "Salad Bar Featured Items",
    "Bakery",
    "Chef's Choice",
    "International",
]);

/**
 * Cleans menu item text by removing dietary tags and extra whitespace
 * @param text - The raw menu item text
 * @returns Cleaned text with dietary tags removed
 */
function cleanMenuItemText(text: string): string {
    return text
        .replace(/\s*\(vgn\)\s*/gi, '') // Remove (vgn) tags
        .replace(/\s*\(v\)\s*/gi, '')   // Remove (v) tags
        .trim()                         // Remove leading/trailing whitespace
        .replace(/\s+/g, ' ');          // Normalize internal whitespace
}

/**
 * Parses a meal table from the new weekly menu layout.
 * @param $ - The Cheerio instance.
 * @param mealId - The ID for the meal's container (e.g., '#dinner-body').
 * @returns A MenuData object for the specified meal.
 */
function parseMealTable($: cheerio.CheerioAPI, mealId: string): MenuData {
  const mealTable = $(mealId);
  if (mealTable.length === 0) {
    return { days: [] };
  }

  // 1. Get all day titles from the table header
  const dayTitles: string[] = [];
  mealTable.find('thead th h4').each((_, th) => {
    // Clean up the title to remove extra whitespace
    const title = $(th).text().trim().replace(/\s+/g, ' ');
    dayTitles.push(title);
  });

  // 2. Initialize the data structure for the final output
  const days: MenuData['days'] = dayTitles.map(title => ({
    title: title,
    items: [],
  }));

  let currentStationTitle = '';

  // 3. Iterate over each row in the table body
  mealTable.find('tbody tr').each((_, rowElement) => {
    const row = $(rowElement);
    const courseHeader = row.find('td.course-row h4');

    if (courseHeader.length > 0) {
      // This is a course header row (e.g., "The Brick")
      currentStationTitle = courseHeader.text().trim();
      
      // Only add the station title if it's not in the exclusion list
      if (!HEADERS_TO_EXCLUDE.has(currentStationTitle)) {
        days.forEach(day => {
          day.items.push({ text: currentStationTitle, class: 'course-row' });
        });
      }

    } else {
      // This is a menu items row, corresponding to the previous header
      row.find('td').each((colIndex, cellElement) => {
        // Ensure we don't go out of bounds if the table is malformed
        if (colIndex >= days.length) return;

        $(cellElement).find('dd').each((itemIndex, itemElement) => {
          const rawItemText = $(itemElement).text().trim();
          const itemText = cleanMenuItemText(rawItemText);

          // Apply the new highlighting rule
          const shouldHighlight = itemIndex === 0 && SECTIONS_TO_HIGHLIGHT.has(currentStationTitle);

          // Add the menu item to the correct day
          days[colIndex].items.push({
            text: itemText,
            highlight: shouldHighlight,
          });
        });
      });
    }
  });

  return { days };
}

export async function GET() {
  try {
    // The new, correct URL for the weekly menu
    const menuUrl = 'https://apps.dining.ucsb.edu/menu/week?dc=portola&m=breakfast&m=brunch&m=lunch&m=dinner&m=late-night&food=';
    
    const response = await fetch(menuUrl, {
      next: { revalidate: 3600 } // Revalidate every hour
    });

    if (!response.ok) {
      throw new Error(`Failed to fetch menu. Status: ${response.status}`);
    }

    const html = await response.text();
    const $ = cheerio.load(html);

    // Parse each meal using the new table-based function
    const breakfastData = parseMealTable($, '#breakfast-body');
    const lunchData = parseMealTable($, '#lunch-body');
    const dinnerData = parseMealTable($, '#dinner-body');

    // The new page doesn't have a single date range, so we'll mark it N/A.
    // The frontend is already designed to construct the title from the day list.
    const allData: AllMealsData = {
      breakfast: breakfastData,
      lunch: lunchData,
      dinner: dinnerData,
      dateRange: "N/A",
    };

    return NextResponse.json(allData);

  } catch (error) {
    console.error('Error fetching or parsing menu:', error);
    const errorMessage = error instanceof Error ? error.message : 'An unknown error occurred';
    return NextResponse.json({ error: `Could not load menu data: ${errorMessage}` }, { status: 500 });
  }
}