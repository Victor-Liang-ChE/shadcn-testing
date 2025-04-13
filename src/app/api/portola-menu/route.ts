import { NextResponse } from 'next/server';
import { JSDOM } from 'jsdom';

// Define interfaces for menu items
interface MenuItem {
    text: string;
    highlight?: boolean;
    class?: string;
}

interface DayMenu {
    title: string;
    items: MenuItem[];
}

interface MenuData {
    days: DayMenu[];
    dateRange: string;
}

// Interface for the combined response
interface AllMealsResponse {
    breakfast: MenuData;
    lunch: MenuData;
    dinner: MenuData;
    dateRange: string; // Keep top-level date range
}

// Function to extract day titles - Return full name and date
function extractDayTitles(document: Document): string[] {
    const potentialHeaders = Array.from(document.querySelectorAll('h4, p.day-header, div.day-title'));
    const dayTitles: string[] = [];
    const dayNameRegex = /^(Mon|Tue|Wed|Thu|Fri|Sat|Sun)/i;

    potentialHeaders.forEach(el => {
        const fullText = el.textContent?.replace(/\s+/g, ' ').trim() || '';
        const match = fullText.match(dayNameRegex);
        if (match && fullText.match(/\d{1,2}\/\d{1,2}/)) {
            if (dayTitles.length < 5) {
                 // Return the full text like "Saturday, 04/12"
                 dayTitles.push(fullText);
            }
        }
    });

    if (dayTitles.length < 5) {
        console.warn(`DEBUG WARN: Found only ${dayTitles.length} potential day titles. Structure might have changed.`);
    }

    // Return full titles or empty array
    return dayTitles.length > 0 ? dayTitles : [];
}

// Function to process menu items for a section
function processTbody(tbody: Element | null, meal: string, highlightLogic: boolean): MenuItem[][] {
    console.log(`Processing ${meal} tbody with highlight: ${highlightLogic}`);
    const menuItems: MenuItem[][] = [[], [], [], [], []];
    if (!tbody) return menuItems;

    let trCount = 0; // Row counter (1-based index like in Python script)
    const rows = tbody.querySelectorAll('tr');
    rows.forEach(tr => {
        trCount++; // Increment row counter for each <tr>
        const tds = tr.querySelectorAll('td');
        const isCourseRow = tr.classList.contains('course-row') || tr.classList.contains('text-center');

        tds.forEach((td, dayIndex) => {
            if (dayIndex >= 5) return; // Only process first 5 columns (days)

            // Ensure the array for the day exists
            if (!menuItems[dayIndex]) {
                menuItems[dayIndex] = [];
            }

            if (isCourseRow) {
                const text = td.textContent?.replace(/\(v\)|\(vgn\)/g, '').trim() || '';
                if (text) {
                    menuItems[dayIndex].push({ text, class: 'course-row' });
                }
            } else {
                const dl = td.querySelector('dl');
                if (dl) {
                    const dds = dl.querySelectorAll('dd');
                    dds.forEach((dd, itemIndex) => { // itemIndex is 0-based index of <dd> within <dl>
                        const text = dd.textContent?.replace(/\(v\)|\(vgn\)/g, '').trim() || '';
                        if (text) {
                            let highlight = false;
                            if (highlightLogic) {
                                if (meal === 'lunch' && itemIndex === 0 && [6, 10, 12, 14].includes(trCount)) {
                                    highlight = true;
                                } else if (meal === 'dinner' && itemIndex === 0 && [6, 10, 12].includes(trCount)) { // Adjusted dinner rows based on previous logic
                                    highlight = true;
                                }
                            }
                            menuItems[dayIndex].push({ text, highlight });
                        }
                    });
                }
            }
        });
    });
    return menuItems;
}

// Helper to create MenuData structure
function createMenuDataObject(titles: string[], items: MenuItem[][], dateRange: string): MenuData {
    const daysData: DayMenu[] = titles.map((title, index) => ({
        title: title,
        items: items[index] || []
    }));
    return { days: daysData, dateRange };
}


export async function GET(request: Request) {
    console.log("DEBUG: Portola menu API route called (fetching all meals)");

    try {
        const url = 'https://apps.dining.ucsb.edu/menu/week?dc=portola&m=breakfast&m=brunch&m=lunch&m=dinner&m=late-night&food=';
        console.log(`DEBUG: Fetching menu from URL: ${url}`);
        const response = await fetch(url, { cache: 'no-store' });

        if (!response.ok) {
            console.error(`DEBUG ERROR: Failed to fetch menu HTML. Status: ${response.status}`);
            throw new Error(`Failed to fetch menu HTML: ${response.statusText}`);
        }
        const htmlContent = await response.text();
        if (!htmlContent) {
            throw new Error("Received empty HTML content.");
        }

        const dom = new JSDOM(htmlContent);
        const document = dom.window.document;

        // --- Extract Titles and Date Range Once ---
        const dayTitles = extractDayTitles(document);
        console.log('DEBUG: Extracted Day titles:', dayTitles);
        const finalDayTitles = dayTitles.length === 5 ? dayTitles : Array.from({ length: 5 }, (_, i) => `Day ${i + 1}`);

        let dateRange = "N/A";
        let dateElementFound = false;
        const allParagraphs = Array.from(document.querySelectorAll('p'));
        const dateRangeRegex = /\d{1,2}\/\d{1,2}\s+to\s+\d{1,2}\/\d{1,2}/;
        for (const p of allParagraphs) {
            if (p.textContent && dateRangeRegex.test(p.textContent)) {
                dateRange = p.textContent.trim();
                dateElementFound = true;
                break;
            }
        }
        if (!dateElementFound) {
             console.warn("DEBUG WARN: Date range element not found using pattern matching, generating fallback dates.");
             // Fallback date range generation if element not found
             const today = new Date(); // Consider timezone if needed
             const futureDate = new Date(today);
             futureDate.setDate(today.getDate() + 4);
             const formatDate = (date: Date) => `${(date.getMonth() + 1).toString().padStart(2, '0')}/${date.getDate().toString().padStart(2, '0')}`;
             dateRange = `${formatDate(today)} to ${formatDate(futureDate)}`;
        }
        console.log(`DEBUG: Date range determined as: ${dateRange}`);

        // --- Process Each Meal ---
        let breakfastItems: MenuItem[][] = [[], [], [], [], []];
        let lunchItems: MenuItem[][] = [[], [], [], [], []];
        let dinnerItems: MenuItem[][] = [[], [], [], [], []];

        // Breakfast
        const breakfastSection = document.querySelector('#breakfast-body');
        if (breakfastSection) {
            const breakfastTbody = breakfastSection.querySelector('tbody');
            breakfastItems = processTbody(breakfastTbody, 'breakfast', false); // No highlight for breakfast
            // Add "See Brunch" message for weekends
            finalDayTitles.forEach((title, index) => {
                if (title.startsWith('Sat') || title.startsWith('Sun')) {
                    if (!breakfastItems[index]) breakfastItems[index] = [];
                    // Clear existing items and add the message
                    breakfastItems[index] = [{ text: "See Lunch/Brunch section for Brunch menu.", class: 'brunch-notice' }];
                }
            });
        } else { console.log('DEBUG: No #breakfast-body section found'); }

        // Lunch/Brunch (Combined)
        const lunchSection = document.querySelector('#lunch-body');
        const brunchSection = document.querySelector('#brunch-body');
        if (lunchSection) {
            const lunchTbody = lunchSection.querySelector('tbody');
            const processedLunch = processTbody(lunchTbody, 'lunch', true); // Highlight lunch
            processedLunch.forEach((dayItems, dayIndex) => {
                 if (!lunchItems[dayIndex]) lunchItems[dayIndex] = [];
                 lunchItems[dayIndex].push(...dayItems);
            });
        } else { console.log('DEBUG: No #lunch-body section found'); }
        if (brunchSection) {
            const brunchTbody = brunchSection.querySelector('tbody');
            const processedBrunch = processTbody(brunchTbody, 'brunch', false); // No highlight brunch
            processedBrunch.forEach((dayItems, dayIndex) => {
                 if (!lunchItems[dayIndex]) lunchItems[dayIndex] = [];
                 lunchItems[dayIndex].push(...dayItems);
            });
        } else { console.log('DEBUG: No #brunch-body section found'); }


        // Dinner
        const dinnerSection = document.querySelector('#dinner-body');
        if (dinnerSection) {
            const dinnerTbody = dinnerSection.querySelector('tbody');
            dinnerItems = processTbody(dinnerTbody, 'dinner', true); // Highlight dinner
        } else { console.log('DEBUG: No #dinner-body section found'); }

        // --- Construct Response ---
        const allData: AllMealsResponse = {
            breakfast: createMenuDataObject(finalDayTitles, breakfastItems, dateRange),
            lunch: createMenuDataObject(finalDayTitles, lunchItems, dateRange),
            dinner: createMenuDataObject(finalDayTitles, dinnerItems, dateRange),
            dateRange: dateRange // Include top-level date range
        };

        console.log(`DEBUG: Successfully processed all meal data.`);
        return NextResponse.json(allData);

    } catch (error: any) {
        console.error(`DEBUG ERROR in API route (all meals):`, error);
        // Return a generic error structure if needed, or specific meal errors
        return NextResponse.json(
            { error: error.message || 'Failed to process menu data' },
            { status: 500 }
        );
    }
}
