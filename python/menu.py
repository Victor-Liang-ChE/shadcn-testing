import sys
import json
import requests
from bs4 import BeautifulSoup
import datetime
import pytz

def process_tbody(tbody, meal, highlight):
    """
    Process the tbody for a particular meal.
    For lunch, highlight is True when requested;
    for brunch entries (meal=='brunch') highlight will be False.
    For breakfast and dinner, the appropriate highlighting logic is applied.
    """
    menu_items = []
    tr_count = 0
    for tr in tbody.find_all('tr'):
        tr_count += 1
        tds = tr.find_all('td')
        for i in range(min(5, len(tds))):
            if len(menu_items) <= i:
                menu_items.append([])
                
            if 'text-center course-row' in tr.get('class', []):
                text_content = tds[i].get_text(separator=' ', strip=True)
                text_content = text_content.replace('(v)', '').replace('(vgn)', '').strip()
                menu_items[i].append({'text': text_content, 'class': 'course-row'})
            else:
                first_dl = tds[i].find('dl')
                if first_dl:
                    for idx, dd in enumerate(first_dl.find_all('dd')):
                        text = dd.text.strip().replace('(v)', '').replace('(vgn)', '').strip()
                        # Logic for lunch with highlighting
                        if meal == 'lunch' and highlight:
                            if idx == 0 and tr_count in [6, 10, 12, 14]:
                                menu_items[i].append({'text': text, 'highlight': True})
                            else:
                                menu_items[i].append({'text': text})
                        # Breakfast: no highlighting
                        elif meal == 'breakfast':
                            menu_items[i].append({'text': text})
                        # Dinner highlighting logic
                        elif meal == 'dinner':
                            if idx == 0 and tr_count in [4, 10, 12]:
                                menu_items[i].append({'text': text, 'highlight': True})
                            else:
                                menu_items[i].append({'text': text})
                        else:
                            # For brunch entries, no highlighting even though processed under lunch
                            menu_items[i].append({'text': text})
    return menu_items

def scrape_menu(meal="dinner"):
    # Fetch the menu page (all meals are included)
    url = 'https://apps.dining.ucsb.edu/menu/week?dc=portola&m=breakfast&m=brunch&m=lunch&m=dinner&m=late-night&food='
    response = requests.get(url)
    soup = BeautifulSoup(response.content, 'html.parser')

    if meal == 'lunch':
        # For Lunch/Brunch mode, try to get both lunch and brunch sections
        merged_menu = []
        # Process lunch section (with highlighting)
        lunch_section = soup.find('div', id='lunch-body')
        if lunch_section:
            tbody_lunch = lunch_section.find('tbody')
            lunch_menu = process_tbody(tbody_lunch, meal='lunch', highlight=True)
            merged_menu = lunch_menu
        # Process brunch section (no highlighting)
        brunch_section = soup.find('div', id='brunch-body')
        if brunch_section:
            tbody_brunch = brunch_section.find('tbody')
            brunch_menu = process_tbody(tbody_brunch, meal='brunch', highlight=False)
            # Extend merged_menu with brunch_menu
            if not merged_menu:
                merged_menu = brunch_menu
            else:
                for i in range(min(len(merged_menu), len(brunch_menu))):
                    merged_menu[i].extend(brunch_menu[i])
        # If neither section is found, return None
        if not merged_menu:
            return None
        return merged_menu
    else:
        # For breakfast or dinner, use the respective section
        menu_section = soup.find('div', id=f'{meal}-body')
        if not menu_section:
            return None
        tbody = menu_section.find('tbody')
        # For breakfast force no highlight; for dinner use highlighting
        if meal == 'breakfast':
            return process_tbody(tbody, meal, highlight=False)
        else:
            return process_tbody(tbody, meal, highlight=True)

def format_menu_response(menu_items, meal_type):
    if menu_items is None:
        return {"error": "No menu available for this week"}
    
    # Get date range for the current week
    pst = pytz.timezone('America/Los_Angeles')
    today = datetime.datetime.now(pst)
    future_date = today + datetime.timedelta(days=4)
    date_range = f"{today.strftime('%m/%d')} to {future_date.strftime('%m/%d')}"
    
    # Format data for frontend
    days = []
    for i, items in enumerate(menu_items):
        day_title = (today + datetime.timedelta(days=i)).strftime('%A, %m/%d')
        days.append({
            "title": day_title,
            "items": items
        })
    
    return {
        "days": days,
        "dateRange": date_range
    }

if __name__ == "__main__":
    # Get meal type from command line argument
    meal_type = "dinner"  # default
    if len(sys.argv) > 1:
        meal_type = sys.argv[1]
    
    # Scrape menu for this meal type
    menu_items = scrape_menu(meal_type)
    
    # Format response and print as JSON
    response = format_menu_response(menu_items, meal_type)
    print(json.dumps(response))