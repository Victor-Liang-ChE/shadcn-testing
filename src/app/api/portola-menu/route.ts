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

// Function to extract day titles with dates and shorten them
function extractDayTitles(document: Document): string[] {
    const dayHeaders = Array.from(document.querySelectorAll('#week-menu-table thead th'));
    // Skip first header (empty cell), take the next 5
    return dayHeaders.slice(1, 6).map(th => {
        // Extract text, remove extra whitespace, should be like "Monday, 04/15"
        const fullText = th.textContent?.replace(/\s+/g, ' ').trim() || '';
        // Extract only the day name part and shorten it
        const dayNameMatch = fullText.match(/^(\w+)/);
        if (dayNameMatch && dayNameMatch[1]) {
            return dayNameMatch[1].substring(0, 3); // e.g., "Mon", "Tue"
        }
        return 'Day'; // Fallback
    });
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

                            // Apply highlighting logic based on Python script
                            // highlightLogic is true only for 'lunch' and 'dinner' main processing
                            if (highlightLogic) {
                                // Lunch highlighting: first item (itemIndex === 0) in specific rows (trCount)
                                if (meal === 'lunch' && itemIndex === 0 && [6, 10, 12, 14].includes(trCount)) {
                                    highlight = true;
                                    console.log(`Highlighted lunch item: ${text} (tr: ${trCount})`);
                                }
                                // Dinner highlighting: first item (itemIndex === 0) in specific rows (trCount)
                                else if (meal === 'dinner' && itemIndex === 0 && [4, 10, 12].includes(trCount)) {
                                    highlight = true;
                                    console.log(`Highlighted dinner item: ${text} (tr: ${trCount})`);
                                }
                            }
                            // Note: 'brunch' is processed with highlightLogic=false, so it won't highlight
                            // 'breakfast' is processed with highlightLogic=false, so it won't highlight

                            menuItems[dayIndex].push({ text, highlight });
                        }
                    });
                }
            }
        });
    });
    return menuItems;
}

export async function GET(request: Request) {
    console.log("DEBUG: Portola menu API route called");
    const { searchParams } = new URL(request.url);
    const meal = searchParams.get('meal') || 'dinner';
    console.log(`DEBUG: Requested meal type: ${meal}`);

    try {
        const url = 'https://apps.dining.ucsb.edu/menu/week?dc=portola&m=breakfast&m=brunch&m=lunch&m=dinner&m=late-night&food=';
        console.log(`DEBUG: Fetching menu from URL: ${url}`);
        const response = await fetch(url, { cache: 'no-store' }); // Disable cache

        if (!response.ok) {
            console.error(`DEBUG ERROR: Failed to fetch menu HTML. Status: ${response.status}`);
            throw new Error(`Failed to fetch menu HTML: ${response.statusText}`);
        }

        const htmlContent = await response.text();
        console.log(`DEBUG: Received HTML content (length: ${htmlContent.length})`);
        if (!htmlContent) {
            throw new Error("Received empty HTML content.");
        }

        const dom = new JSDOM(htmlContent);
        const document = dom.window.document;

        const dayTitles = extractDayTitles(document); // Use the updated function
        console.log('DEBUG: Extracted Day titles:', dayTitles);
        // Fallback if titles extraction fails
        const finalDayTitles = dayTitles.length === 5 ? dayTitles : Array.from({ length: 5 }, (_, i) => `Day ${i + 1}`);

        let processedMenuItems: MenuItem[][] = [[], [], [], [], []];
        let menuFound = false;

        if (meal === 'lunch') {
            console.log('DEBUG: Processing lunch/brunch menu');
            const lunchSection = document.querySelector('#lunch-body');
            const brunchSection = document.querySelector('#brunch-body');

            if (lunchSection) {
                const lunchTbody = lunchSection.querySelector('tbody');
                if (lunchTbody) {
                    console.log('DEBUG: Found lunch tbody');
                    // Process lunch with highlighting enabled (highlightLogic = true)
                    const lunchMenu = processTbody(lunchTbody, 'lunch', true);
                    lunchMenu.forEach((dayItems, dayIndex) => processedMenuItems[dayIndex].push(...dayItems));
                    menuFound = true;
                } else { console.log('DEBUG: No tbody found in #lunch-body'); }
            } else { console.log('DEBUG: No #lunch-body section found'); }

            if (brunchSection) {
                const brunchTbody = brunchSection.querySelector('tbody');
                if (brunchTbody) {
                    console.log('DEBUG: Found brunch tbody');
                    // Process brunch with highlighting disabled (highlightLogic = false)
                    const brunchMenu = processTbody(brunchTbody, 'brunch', false);
                    brunchMenu.forEach((dayItems, dayIndex) => {
                        if (!processedMenuItems[dayIndex]) {
                            processedMenuItems[dayIndex] = [];
                        }
                        processedMenuItems[dayIndex].push(...dayItems);
                    });
                    menuFound = true;
                } else { console.log('DEBUG: No tbody found in #brunch-body'); }
            } else { console.log('DEBUG: No #brunch-body section found'); }

        } else { // breakfast or dinner
            const sectionId = `${meal}-body`;
            console.log(`DEBUG: Processing ${meal} menu (selector: #${sectionId})`);
            const section = document.querySelector(`#${sectionId}`);

            if (section) {
                const tbody = section.querySelector('tbody');
                if (tbody) {
                    console.log(`DEBUG: Found ${meal} tbody`);
                    // Enable highlighting only for dinner (highlightLogic = meal === 'dinner')
                    processedMenuItems = processTbody(tbody, meal, meal === 'dinner');
                    menuFound = true;
                } else { console.log(`DEBUG: No tbody found in #${sectionId}`); }
            } else { console.log(`DEBUG: No #${sectionId} section found`); }
        }

        // Get date range
        let dateRange = "N/A";
        const dateElement = document.querySelector('.date-range'); // Assuming this class exists
        if (dateElement?.textContent) {
            dateRange = dateElement.textContent.trim();
        } else {
             // Fallback date range generation if element not found
             console.warn("DEBUG WARN: '.date-range' element not found, generating fallback dates.");
             const today = new Date();
             const futureDate = new Date(today);
             futureDate.setDate(today.getDate() + 4);
             const formatDate = (date: Date) => `${(date.getMonth() + 1).toString().padStart(2, '0')}/${date.getDate().toString().padStart(2, '0')}`;
             dateRange = `${formatDate(today)} to ${formatDate(futureDate)}`;
        }
        console.log(`DEBUG: Date range determined as: ${dateRange}`);

        // Check if any menu items were actually found across all days
        const totalItems = processedMenuItems.reduce((sum, day) => sum + (day?.length || 0), 0);
        if (!menuFound || totalItems === 0) {
            console.log(`DEBUG: No menu items found for ${meal}. Returning empty.`);
            // Return success but with empty data and the determined date range
            return NextResponse.json({ days: [], dateRange: dateRange });
        }

        const daysData: DayMenu[] = finalDayTitles.map((title, index) => ({
            title: title, // Use extracted short titles or fallback
            items: processedMenuItems[index] || [] // Ensure items array exists
        }));

        const menuData: MenuData = { days: daysData, dateRange };
        console.log(`DEBUG: Successfully processed ${meal} menu data. Highlighted items: ${daysData.flatMap(day => day.items).filter(item => item.highlight).length}`);
        return NextResponse.json(menuData);
    } catch (error: any) {
        console.error(`DEBUG ERROR in API route for meal "${meal}":`, error);
        return NextResponse.json(
            { error: error.message || 'Failed to process menu data' },
            { status: 500 }
        );
    }
}
