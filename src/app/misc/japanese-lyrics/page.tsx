"use client";

import React, { useState, useEffect, useCallback } from 'react';
import { Button } from "@/components/ui/button";
import { Textarea } from "@/components/ui/textarea";
import {
    Select,
    SelectContent,
    SelectItem,
    SelectTrigger,
    SelectValue,
} from "@/components/ui/select";
import {
    Accordion,
    AccordionContent,
    AccordionItem,
    AccordionTrigger,
} from "@/components/ui/accordion";
import { Loader2 } from "lucide-react"; // Example loading icon

// Debounce function
function debounce<F extends (...args: any[]) => any>(func: F, waitFor: number) {
    let timeoutId: ReturnType<typeof setTimeout> | null = null;

    return (...args: Parameters<F>): Promise<ReturnType<F>> => {
        return new Promise((resolve) => {
            if (timeoutId) {
                clearTimeout(timeoutId);
            }

            timeoutId = setTimeout(() => {
                timeoutId = null;
                resolve(func(...args));
            }, waitFor);
        });
    };
}

interface ProcessedLyrics {
    plain: string;
    hiragana: string;
    katakana: string;
    romanji: string;
}

export default function JapaneseLyricsPage() {
    const [inputText, setInputText] = useState<string>('');
    const [processedLyrics, setProcessedLyrics] = useState<ProcessedLyrics | null>(null);
    const [showFurigana, setShowFurigana] = useState<boolean>(false);
    const [furiganaFormat, setFuriganaFormat] = useState<'hiragana' | 'katakana' | 'romanji'>('hiragana');
    const [showInput, setShowInput] = useState<boolean>(true);
    const [isLoading, setIsLoading] = useState<boolean>(false);
    const [error, setError] = useState<string | null>(null);

    // Debounced API call function
    const fetchProcessedLyrics = useCallback(
        debounce(async (text: string): Promise<ProcessedLyrics | null> => {
            if (!text.trim()) {
                return { plain: "", hiragana: "", katakana: "", romanji: "" };
            }
            setIsLoading(true);
            setError(null);
            try {
                // --- BACKEND API CALL ---
                // Replace with your actual API endpoint
                const response = await fetch('/api/process-lyrics', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify({ text }),
                });
                // --- END BACKEND API CALL ---

                if (!response.ok) {
                    throw new Error(`API error: ${response.statusText}`);
                }
                const data: ProcessedLyrics = await response.json();

                // Optional: Apply line merging on the client-side if backend doesn't do it
                // This assumes the backend returns arrays of strings per line
                /*
                const mergeLines = (content: string) => mergeAdjacentDuplicates(content.split('<br>')).join('<br>');
                return {
                    plain: mergeLines(data.plain),
                    hiragana: mergeLines(data.hiragana),
                    katakana: mergeLines(data.katakana),
                    romanji: mergeLines(data.romanji),
                };
                */
               // If backend already returns merged strings with <br>:
                return data;

            } catch (err) {
                console.error("Failed to fetch processed lyrics:", err);
                setError(err instanceof Error ? err.message : "An unknown error occurred");
                return null; // Indicate error or return previous state?
            } finally {
                setIsLoading(false);
            }
        }, 500), // Debounce API calls by 500ms
        []
    );

    // Effect to call API when input text changes
    useEffect(() => {
        // Define an async function inside useEffect
        const processText = async () => {
            // Await the result of the debounced function call
            const data = await fetchProcessedLyrics(inputText);
            if (data) {
                // Now 'data' is ProcessedLyrics | null, not a Promise
                setProcessedLyrics(data);
            }
            // Optional: Handle the case where data is null (e.g., API error)
            // else { setProcessedLyrics(null); }
        };

        // Call the async function
        processText();

        // No cleanup needed specifically for debounce in this pattern
    }, [inputText, fetchProcessedLyrics]); // Dependencies remain the same

    const handleInputChange = (event: React.ChangeEvent<HTMLTextAreaElement>) => {
        setInputText(event.target.value);
    };

    const toggleFurigana = () => {
        setShowFurigana(prev => !prev);
    };

    const toggleInputVisibility = () => {
        setShowInput(prev => !prev);
    };

    // Updated for Shadcn Select component
    const handleFormatChange = (value: string) => {
        setFuriganaFormat(value as 'hiragana' | 'katakana' | 'romanji');
    };

    const getDisplayText = () => {
        if (!processedLyrics) return "";
        if (showFurigana) {
            // Ensure the format exists, otherwise default to hiragana
            return processedLyrics[furiganaFormat] || processedLyrics.hiragana;
        }
        return processedLyrics.plain;
    };

    return (
        // Main container - relative positioning context
        <div className="relative h-screen w-screen overflow-hidden bg-background text-foreground">
            {/* Left Control Panel - Added overflow-y-auto */}
            <div
                className="absolute left-0 top-0 z-10 flex h-full w-[220px] flex-col gap-4 overflow-y-auto bg-background text-foreground p-4"
            >
                {/* Control Buttons and Select */}
                <Button onClick={toggleInputVisibility} variant="secondary">
                    {showInput ? "Hide Input" : "Show Input"}
                </Button>
                <Button onClick={toggleFurigana} variant="secondary">
                    Toggle Furigana
                </Button>
                <Select
                    value={furiganaFormat}
                    onValueChange={handleFormatChange}
                    disabled={!showFurigana}
                >
                    <SelectTrigger className="w-full">
                        <SelectValue placeholder="Select Format" />
                    </SelectTrigger>
                    <SelectContent>
                        <SelectItem value="hiragana">Hiragana</SelectItem>
                        <SelectItem value="katakana">Katakana</SelectItem>
                        <SelectItem value="romanji">Romanji</SelectItem>
                    </SelectContent>
                </Select>

                {/* Input Area Container - Conditional rendering/styling */}
                <div className={`flex flex-col transition-opacity duration-300 ${showInput ? 'opacity-100' : 'opacity-0 pointer-events-none h-0 overflow-hidden'}`}>
                    <Textarea
                        value={inputText}
                        onChange={handleInputChange}
                        placeholder="Paste the lyrics here..."
                        className="min-h-[200px] resize-none placeholder:text-foreground/60"
                    />
                    {isLoading && (
                        <div className="mt-4 flex items-center justify-center gap-2 text-sm text-primary-foreground">
                            <Loader2 className="h-4 w-4 animate-spin" />
                            Loading...
                        </div>
                    )}
                    {error && (
                        <div className="mt-4 text-center text-sm text-red-400">
                            Error: {error}
                        </div>
                    )}
                </div>

                {/* Journal Section - Removed mt-auto */}
                <div className="w-full"> {/* Removed mt-auto */}
                     <Accordion type="single" collapsible className="w-full">
                        <AccordionItem value="journal" className="border-b-0">
                            <AccordionTrigger className="px-0 py-2 text-sm font-medium hover:no-underline"> {/* Adjusted padding */}
                                Journal
                            </AccordionTrigger>
                            <AccordionContent className="max-h-[30vh] overflow-y-auto pb-1"> {/* Adjusted padding */}
                                {/* Using pre-wrap to preserve line breaks from the original markdown */}
                                {/* Changed text color for better contrast on primary background */}
                                <p className="whitespace-pre-wrap text-sm text-foreground/80">
{`- kinda been getting bored of music lately, i gotta listen to some old songs again
- ame to cappuccino is such a good song, she actually recommended the best songs ahhahhaha
- gonna try singing it by following the lyrics, but it sucks that i cant read most of the kanji lmao
- can i somehow enable the hirigana that sits on top of the kanji? oh i can cool
- oh its called furigana ok
- i wanna be able to toggle the furigana on and off so i can test my memory but spicetify doesnt let me toggle it on and off quickly
- ill just build my own furigana toggle thing then
- okok so i want a toggle button, probably a input box for copy and pasting the lyrics and then...
- detect the kanji? i know the kanji has multiple readings like the kanji for "one" can sound like "hito" or "ii"
- i want a dictionary that knows the context of the sentence and then give me the right kanji reading
- gonna look it up... something something tokenizer... like the tokens chatgpt uses??
- to split the sentence into words... ok
- gonna try the package MeCab i guess it looks like it has the best numbers on this tokenizer comparison chart
- nawww i know this doesnt sound right... 揺蕩 sounds like its sang: tayutau not youutouu...
- prob something wrong with the tokenizer, gonna try another one
- ngl maybe the slowest one on that chart, sudachi, might give me more accurate results
- like do i really care about how fast the tokenizer is? you not gonna troll me and load for a minute right?
- yea thats what i thought lmao, barely slower
- and it gave me the right reading too ok ill just stick with this then
- aight it looks like its coming together now! i should prob add more options like katakana and romanji for the furigana
- ahh sudachipy doesnt have a romanji conversion, can i hard code the romanji comversion in???
- jk lmao im gonna have to use another package for that
- pykakasi it is then...
- awesome its working!!! i kinda want to hide the input box when im done using it tho, its so big an distracting
- translation next??? i know the google translate is kinda ass... ill try anyway
- yea it is ass LMAO
- can i just take the translations from musixmatch? i know thats what spicetify does
- man they got a paid API that aint worth it for me
- ahh whatever ill just live without the translations ;-;`}
                                </p>
                            </AccordionContent>
                        </AccordionItem>
                    </Accordion>
                </div>
            </div>
            {/* Lyrics Display Area - Positioned next to left panel, scrolls */}
            {/* Changed bottom-[60px] to bottom-0 */}
            <div className="absolute bottom-0 left-[220px] right-0 top-0 overflow-y-auto p-5">
                <div className="mx-auto max-w-3xl text-center">
                    <div
                        className="whitespace-pre-wrap text-center text-2xl leading-[2.3] [&_ruby_rt]:text-[0.6em] [&_ruby_rt]:opacity-80"
                        dangerouslySetInnerHTML={{ __html: getDisplayText() }}
                    />
                </div>
            </div>
        </div>
    );
}
