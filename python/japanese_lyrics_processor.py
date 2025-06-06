import re
from sudachipy import dictionary, tokenizer
import pykakasi
import sys

# Initialize tokenizer and converter
try:
    tokenizer_obj = dictionary.Dictionary().create()
    mode = tokenizer.Tokenizer.SplitMode.C
    kks = pykakasi.kakasi()
    print("DEBUG: Sudachipy and PyKakasi initialized successfully.", file=sys.stderr)
except Exception as e:
    print(f"DEBUG ERROR: Failed to initialize Sudachipy/PyKakasi: {e}", file=sys.stderr)
    # Depending on the environment, you might want to exit or handle this differently
    # sys.exit(1) # Uncomment if initialization failure should stop the app

def katakana_to_hiragana(text):
    """Converts Katakana characters in a string to Hiragana."""
    result = []
    for ch in text:
        code = ord(ch)
        # Check if character is within the Katakana range
        if 0x30A1 <= code <= 0x30F6:
            # Convert to Hiragana by subtracting the offset
            result.append(chr(code - 0x60))
        else:
            result.append(ch)
    return "".join(result)

def add_furigana(text, format_type='hiragana'):
    """
    Adds furigana to Japanese text using <ruby> tags.
    format_type can be 'hiragana', 'katakana', or 'romanji'.
    """
    print(f"DEBUG: add_furigana called with format_type='{format_type}'", file=sys.stderr) # DEBUG
    if not text:
        print("DEBUG: add_furigana received empty text, returning []", file=sys.stderr) # DEBUG
        return []

    lines = text.splitlines()
    processed_lines = []
    kanji_pattern = re.compile(r'[\u4e00-\u9faf]') # Matches Kanji characters
    # Matches strings consisting entirely of Hiragana or Katakana (to avoid adding ruby to them)
    non_kanji_pattern = re.compile(r'^(?:[\u3040-\u309F\u30A0-\u30FF]+)$')

    for line_num, line in enumerate(lines): # DEBUG: Added line number
        print(f"DEBUG: Processing line {line_num + 1}: '{line}'", file=sys.stderr) # DEBUG
        if not line.strip(): # Keep empty lines
             processed_lines.append("")
             print(f"DEBUG: Line {line_num + 1} is empty, appending empty string.", file=sys.stderr) # DEBUG
             continue
        try:
            tokens = tokenizer_obj.tokenize(line, mode)
            result = []
            print(f"DEBUG: Tokens for line {line_num + 1}: {[t.surface() for t in tokens]}", file=sys.stderr) # DEBUG
            for token_num, token in enumerate(tokens): # DEBUG: Added token number
                surface = token.surface()
                reading_katakana = token.reading_form()
                print(f"DEBUG:   Token {token_num + 1}: Surface='{surface}', Reading(kata)='{reading_katakana}'", file=sys.stderr) # DEBUG

                # Determine the reading based on the requested format
                reading = ""
                if format_type == 'hiragana':
                    reading = katakana_to_hiragana(reading_katakana)
                elif format_type == 'katakana':
                    reading = reading_katakana
                elif format_type == 'romanji':
                    conversion = kks.convert(reading_katakana)
                    reading = "".join(item['hepburn'] for item in conversion)
                else: # Default to hiragana if format is unknown
                    reading = katakana_to_hiragana(reading_katakana)
                print(f"DEBUG:     Reading({format_type})='{reading}'", file=sys.stderr) # DEBUG

                # Add ruby tags only if:
                # 1. The surface contains Kanji.
                # 2. The surface is not purely Hiragana/Katakana.
                # 3. The surface form is different from the reading form.
                has_kanji = kanji_pattern.search(surface)
                is_not_kana_only = not non_kanji_pattern.fullmatch(surface)
                is_different = surface != reading
                print(f"DEBUG:     Conditions: has_kanji={bool(has_kanji)}, is_not_kana_only={is_not_kana_only}, is_different={is_different}", file=sys.stderr) # DEBUG

                if has_kanji and is_not_kana_only and is_different:
                    ruby_tag = f"<ruby>{surface}<rt>{reading}</rt></ruby>"
                    result.append(ruby_tag)
                    print(f"DEBUG:     Added ruby: {ruby_tag}", file=sys.stderr) # DEBUG
                else:
                    result.append(surface)
                    print(f"DEBUG:     Added surface: '{surface}'", file=sys.stderr) # DEBUG
            processed_line = "".join(result)
            processed_lines.append(processed_line)
            print(f"DEBUG: Finished line {line_num + 1}, result: '{processed_line}'", file=sys.stderr) # DEBUG
        except Exception as e:
             print(f"DEBUG ERROR: Error tokenizing/processing line: '{line}'. Error: {e}", file=sys.stderr)
             processed_lines.append(line) # Append original line on error

    print(f"DEBUG: add_furigana finished for format '{format_type}', returning {len(processed_lines)} lines.", file=sys.stderr) # DEBUG
    return processed_lines

def merge_adjacent_duplicates(lines):
    """Merges adjacent identical lines, appending 'xN' for N repetitions."""
    if not lines:
        return []
    merged = []
    count = 1
    # Use first line's processing status to handle potential initial empty lines
    prev_line = lines[0]

    for i in range(1, len(lines)):
        if lines[i] == prev_line:
            count += 1
        else:
            # Append the previous line (with count if > 1)
            merged.append(f"{prev_line} x{count}" if count > 1 else prev_line)
            # Reset for the new line
            prev_line = lines[i]
            count = 1
    # Append the last line sequence
    merged.append(f"{prev_line} x{count}" if count > 1 else prev_line)
    return merged

def process_lyrics(text):
    """
    Processes the input text to generate plain text and text with furigana
    in hiragana, katakana, and romanji formats. Merges adjacent duplicate lines.
    """
    print(f"DEBUG: process_lyrics called with text:\n---\n{text}\n---", file=sys.stderr) # DEBUG
    if not text:
        print("DEBUG: process_lyrics received empty text, returning empty dict.", file=sys.stderr) # DEBUG
        return {'plain': "", 'hiragana': "", 'katakana': "", 'romanji': ""}

    try:
        # Process for each format
        hiragana_lines = add_furigana(text, 'hiragana')
        katakana_lines = add_furigana(text, 'katakana')
        romanji_lines = add_furigana(text, 'romanji')
        plain_lines = text.splitlines() # Keep original lines for plain text

        # Merge duplicates for each format
        merged_plain = merge_adjacent_duplicates(plain_lines)
        merged_hiragana = merge_adjacent_duplicates(hiragana_lines)
        merged_katakana = merge_adjacent_duplicates(katakana_lines)
        merged_romanji = merge_adjacent_duplicates(romanji_lines)

        # Join lines with <br> for HTML display
        result_dict = {
            'plain': "<br>".join(merged_plain),
            'hiragana': "<br>".join(merged_hiragana),
            'katakana': "<br>".join(merged_katakana),
            'romanji': "<br>".join(merged_romanji)
        }
        print(f"DEBUG: process_lyrics finished successfully. Returning dict:", file=sys.stderr) # DEBUG
        # Print dictionary contents for easier debugging
        for key, value in result_dict.items():
             print(f"DEBUG:   {key}: {value[:100]}..." if len(value) > 100 else f"DEBUG:   {key}: {value}", file=sys.stderr) # DEBUG
        return result_dict
    except Exception as e:
        print(f"DEBUG ERROR: Failed during process_lyrics: {e}", file=sys.stderr)
        # Return error structure or raise exception
        return {"error": f"Processing failed: {e}"}

