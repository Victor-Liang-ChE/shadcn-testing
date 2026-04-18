'use client';

import { useEffect, useRef, useState } from 'react';
import {
  Cloud,
  Download,
  FolderOpen,
  Loader2,
  LogIn,
  LogOut,
  RefreshCw,
  Save,
  Trash2,
  Upload,
} from 'lucide-react';

import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import {
  Sheet,
  SheetContent,
  SheetDescription,
  SheetHeader,
  SheetTitle,
} from '@/components/ui/sheet';
import { useCanopyStore } from '@/lib/canopy/store';
import {
  CANOPY_BROWSER_DRAFT_KEY,
  CANOPY_BROWSER_DRAFT_META_KEY,
  buildFlowsheetDownloadFilename,
  parseCanopyFlowsheetJson,
  sanitizeCloudSaveName,
  type CanopyDraftMetadata,
} from '@/lib/canopy/save-utils';

interface CanopySaveSheetProps {
  open: boolean;
  onOpenChange: (open: boolean) => void;
}

interface CloudSaveSummary {
  id: string;
  name: string;
  created_at: string;
  updated_at: string;
}

interface SessionUser {
  id: string;
  email: string | null;
  fullName: string | null;
}

function formatTimestamp(value: string | null): string {
  if (!value) return 'Never';
  const date = new Date(value);
  return Number.isNaN(date.getTime()) ? 'Unknown' : date.toLocaleString();
}

export default function CanopySaveSheet({ open, onOpenChange }: CanopySaveSheetProps) {
  const exportFlowsheet = useCanopyStore(s => s.exportFlowsheet);
  const importFlowsheet = useCanopyStore(s => s.importFlowsheet);
  const compounds = useCanopyStore(s => s.compounds);
  const fluidPackage = useCanopyStore(s => s.fluidPackage);
  const streams = useCanopyStore(s => s.streams);
  const units = useCanopyStore(s => s.units);
  const nodePositions = useCanopyStore(s => s.nodePositions);
  const designSpecs = useCanopyStore(s => s.designSpecs);
  const sensitivities = useCanopyStore(s => s.sensitivities);
  const calculatorBlocks = useCanopyStore(s => s.calculatorBlocks);

  const [saveName, setSaveName] = useState('Canopy Flowsheet');
  const [statusMessage, setStatusMessage] = useState<string | null>(null);
  const [errorMessage, setErrorMessage] = useState<string | null>(null);
  const [dropActive, setDropActive] = useState(false);
  const [browserDraftMeta, setBrowserDraftMeta] = useState<CanopyDraftMetadata | null>(null);
  const [sessionUser, setSessionUser] = useState<SessionUser | null>(null);
  const [cloudSaves, setCloudSaves] = useState<CloudSaveSummary[]>([]);
  const [selectedCloudSaveId, setSelectedCloudSaveId] = useState<string | null>(null);
  const [busyAction, setBusyAction] = useState<string | null>(null);
  const fileInputRef = useRef<HTMLInputElement | null>(null);

  useEffect(() => {
    try {
      const raw = window.localStorage.getItem(CANOPY_BROWSER_DRAFT_META_KEY);
      if (raw) setBrowserDraftMeta(JSON.parse(raw) as CanopyDraftMetadata);
    } catch {
      setBrowserDraftMeta(null);
    }
  }, []);

  useEffect(() => {
    const timeout = window.setTimeout(() => {
      try {
        const metadata: CanopyDraftMetadata = {
          name: sanitizeCloudSaveName(saveName),
          savedAt: new Date().toISOString(),
        };
        window.localStorage.setItem(CANOPY_BROWSER_DRAFT_KEY, exportFlowsheet());
        window.localStorage.setItem(CANOPY_BROWSER_DRAFT_META_KEY, JSON.stringify(metadata));
        setBrowserDraftMeta(metadata);
      } catch {
        // Ignore local draft persistence failures.
      }
    }, 700);

    return () => window.clearTimeout(timeout);
  }, [
    exportFlowsheet,
    saveName,
    compounds,
    fluidPackage,
    streams,
    units,
    nodePositions,
    designSpecs,
    sensitivities,
    calculatorBlocks,
  ]);

  useEffect(() => {
    if (!open) return;
    void refreshSession();
  }, [open]);

  useEffect(() => {
    if (!sessionUser) {
      setCloudSaves([]);
      setSelectedCloudSaveId(null);
      return;
    }
    void refreshCloudSaves();
  }, [sessionUser]);

  const refreshSession = async () => {
    setErrorMessage(null);
    try {
      const response = await fetch('/api/canopy/auth/session', { cache: 'no-store' });
      const payload = await response.json();
      if (!response.ok) throw new Error(payload.error || 'Failed to read auth session.');
      setSessionUser(payload.user ?? null);
    } catch (error) {
      setErrorMessage(error instanceof Error ? error.message : 'Failed to read auth session.');
      setSessionUser(null);
    }
  };

  const refreshCloudSaves = async () => {
    setBusyAction('refresh');
    setErrorMessage(null);
    try {
      const response = await fetch('/api/canopy/cloud-saves', { cache: 'no-store' });
      const payload = await response.json();
      if (!response.ok) throw new Error(payload.error || 'Failed to load cloud saves.');
      setCloudSaves(payload.saves ?? []);
      if (selectedCloudSaveId && !(payload.saves ?? []).some((save: CloudSaveSummary) => save.id === selectedCloudSaveId)) {
        setSelectedCloudSaveId(null);
      }
    } catch (error) {
      setErrorMessage(error instanceof Error ? error.message : 'Failed to load cloud saves.');
    } finally {
      setBusyAction(null);
    }
  };

  const importFromText = (json: string, sourceLabel: string, fallbackName?: string) => {
    const parsed = parseCanopyFlowsheetJson(json);
    if (!parsed.ok) {
      setErrorMessage(parsed.error);
      return;
    }
    const ok = importFlowsheet(json);
    if (!ok) {
      setErrorMessage(`Failed to import ${sourceLabel}.`);
      return;
    }
    if (fallbackName) setSaveName(sanitizeCloudSaveName(fallbackName));
    setErrorMessage(null);
    setStatusMessage(`${sourceLabel} loaded.`);
  };

  const handleLocalDownload = () => {
    const json = exportFlowsheet();
    const blob = new Blob([json], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = buildFlowsheetDownloadFilename(saveName);
    a.click();
    URL.revokeObjectURL(url);
    setStatusMessage('Flowsheet downloaded as JSON.');
    setErrorMessage(null);
  };

  const handleFileList = async (files: FileList | null) => {
    const file = files?.[0];
    if (!file) return;
    try {
      const text = await file.text();
      importFromText(text, file.name, file.name.replace(/\.json$/i, ''));
    } catch (error) {
      setErrorMessage(error instanceof Error ? error.message : 'Failed to read the selected file.');
    }
  };

  const handleGoogleSignIn = () => {
    window.location.href = '/api/canopy/auth/google/start?next=/canopy';
  };

  const handleSignOut = async () => {
    setBusyAction('signout');
    setErrorMessage(null);
    try {
      const response = await fetch('/api/canopy/auth/session', { method: 'DELETE' });
      const payload = await response.json();
      if (!response.ok) throw new Error(payload.error || 'Failed to sign out.');
      setSessionUser(null);
      setCloudSaves([]);
      setSelectedCloudSaveId(null);
      setStatusMessage('Signed out.');
    } catch (error) {
      setErrorMessage(error instanceof Error ? error.message : 'Failed to sign out.');
    } finally {
      setBusyAction(null);
    }
  };

  const handleCloudSave = async (mode: 'new' | 'update') => {
    if (!sessionUser) {
      setErrorMessage('Sign in first to save to the cloud.');
      return;
    }

    const url = mode === 'new'
      ? '/api/canopy/cloud-saves'
      : `/api/canopy/cloud-saves/${selectedCloudSaveId}`;
    const method = mode === 'new' ? 'POST' : 'PUT';

    setBusyAction(mode === 'new' ? 'save-new' : 'save-update');
    setErrorMessage(null);
    try {
      const response = await fetch(url, {
        method,
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          name: saveName,
          flowsheetJson: exportFlowsheet(),
        }),
      });
      const payload = await response.json();
      if (!response.ok) throw new Error(payload.error || 'Cloud save failed.');
      setSaveName(payload.save.name);
      setSelectedCloudSaveId(payload.save.id);
      setStatusMessage(mode === 'new' ? 'Cloud save created.' : 'Cloud save updated.');
      await refreshCloudSaves();
    } catch (error) {
      setErrorMessage(error instanceof Error ? error.message : 'Cloud save failed.');
    } finally {
      setBusyAction(null);
    }
  };

  const handleCloudLoad = async (save: CloudSaveSummary) => {
    setBusyAction(`load-${save.id}`);
    setErrorMessage(null);
    try {
      const response = await fetch(`/api/canopy/cloud-saves/${save.id}`, { cache: 'no-store' });
      const payload = await response.json();
      if (!response.ok) throw new Error(payload.error || 'Failed to load cloud save.');
      importFromText(JSON.stringify(payload.save.flowsheet_json, null, 2), `Cloud save "${payload.save.name}"`, payload.save.name);
      setSelectedCloudSaveId(payload.save.id);
    } catch (error) {
      setErrorMessage(error instanceof Error ? error.message : 'Failed to load cloud save.');
    } finally {
      setBusyAction(null);
    }
  };

  const handleCloudDelete = async (save: CloudSaveSummary) => {
    setBusyAction(`delete-${save.id}`);
    setErrorMessage(null);
    try {
      const response = await fetch(`/api/canopy/cloud-saves/${save.id}`, { method: 'DELETE' });
      const payload = await response.json();
      if (!response.ok) throw new Error(payload.error || 'Failed to delete cloud save.');
      if (selectedCloudSaveId === save.id) setSelectedCloudSaveId(null);
      setStatusMessage(`Deleted "${save.name}".`);
      await refreshCloudSaves();
    } catch (error) {
      setErrorMessage(error instanceof Error ? error.message : 'Failed to delete cloud save.');
    } finally {
      setBusyAction(null);
    }
  };

  const handleRestoreBrowserDraft = () => {
    const draft = window.localStorage.getItem(CANOPY_BROWSER_DRAFT_KEY);
    if (!draft) {
      setErrorMessage('No browser draft is available.');
      return;
    }
    importFromText(draft, 'Browser draft', browserDraftMeta?.name);
  };

  const handleClearBrowserDraft = () => {
    window.localStorage.removeItem(CANOPY_BROWSER_DRAFT_KEY);
    window.localStorage.removeItem(CANOPY_BROWSER_DRAFT_META_KEY);
    setBrowserDraftMeta(null);
    setStatusMessage('Browser draft cleared.');
    setErrorMessage(null);
  };

  return (
    <Sheet open={open} onOpenChange={onOpenChange}>
      <SheetContent side="right" className="sm:max-w-xl overflow-y-auto">
        <SheetHeader>
          <SheetTitle>Save And Load</SheetTitle>
          <SheetDescription>
            Use local JSON files, browser draft restore, or a server-managed cloud library.
          </SheetDescription>
        </SheetHeader>

        <div className="mt-6 space-y-6">
          <section className="space-y-3 rounded-lg border border-border p-4">
            <div>
              <h3 className="text-sm font-semibold">Save Name</h3>
              <p className="text-xs text-muted-foreground">Used for downloaded filenames and cloud saves.</p>
            </div>
            <Input
              value={saveName}
              onChange={(event) => setSaveName(event.target.value)}
              placeholder="Canopy Flowsheet"
            />
          </section>

          <section className="space-y-3 rounded-lg border border-border p-4">
            <div>
              <h3 className="text-sm font-semibold">Local File</h3>
              <p className="text-xs text-muted-foreground">Download the current flowsheet or import a `.json` file by click or drag-and-drop.</p>
            </div>
            <div className="flex flex-wrap gap-2">
              <Button type="button" onClick={handleLocalDownload} className="text-xs">
                <Download className="w-3.5 h-3.5 mr-1" />
                Save Local
              </Button>
              <Button
                type="button"
                variant="outline"
                onClick={() => fileInputRef.current?.click()}
                className="text-xs"
              >
                <FolderOpen className="w-3.5 h-3.5 mr-1" />
                Load Local
              </Button>
              <input
                ref={fileInputRef}
                type="file"
                accept=".json,application/json"
                className="hidden"
                onChange={(event) => void handleFileList(event.target.files)}
              />
            </div>
            <div
              className={`rounded-lg border border-dashed p-4 text-sm transition-colors ${
                dropActive ? 'border-blue-500 bg-blue-500/5' : 'border-border'
              }`}
              onDragOver={(event) => {
                event.preventDefault();
                setDropActive(true);
              }}
              onDragLeave={(event) => {
                event.preventDefault();
                setDropActive(false);
              }}
              onDrop={(event) => {
                event.preventDefault();
                setDropActive(false);
                void handleFileList(event.dataTransfer.files);
              }}
            >
              <div className="flex items-center gap-2 text-muted-foreground">
                <Upload className="w-4 h-4" />
                <span>Drop a Canopy `.json` file here to load it.</span>
              </div>
            </div>
          </section>

          <section className="space-y-3 rounded-lg border border-border p-4">
            <div>
              <h3 className="text-sm font-semibold">Browser Draft</h3>
              <p className="text-xs text-muted-foreground">Canopy automatically keeps a local draft in this browser.</p>
            </div>
            <div className="text-xs text-muted-foreground">
              Last autosave: {formatTimestamp(browserDraftMeta?.savedAt ?? null)}
            </div>
            <div className="flex flex-wrap gap-2">
              <Button type="button" variant="outline" onClick={handleRestoreBrowserDraft} className="text-xs">
                Restore Draft
              </Button>
              <Button type="button" variant="outline" onClick={handleClearBrowserDraft} className="text-xs">
                Clear Draft
              </Button>
            </div>
          </section>

          <section className="space-y-3 rounded-lg border border-border p-4">
            <div className="flex items-start justify-between gap-3">
              <div>
                <h3 className="text-sm font-semibold">Cloud Library</h3>
                <p className="text-xs text-muted-foreground">
                  Google sign-in handled through server routes. The browser no longer needs a Supabase client key.
                </p>
              </div>
              <Cloud className="w-4 h-4 text-muted-foreground mt-0.5" />
            </div>

            {!sessionUser && (
              <Button type="button" onClick={handleGoogleSignIn} className="text-xs">
                <LogIn className="w-3.5 h-3.5 mr-1" />
                Sign In With Google
              </Button>
            )}

            {sessionUser && (
              <div className="space-y-3">
                <div className="flex items-center justify-between gap-3 rounded-md border border-border/70 p-3">
                  <div>
                    <div className="text-sm font-medium">{sessionUser.email ?? sessionUser.fullName ?? 'Signed in user'}</div>
                    <div className="text-xs text-muted-foreground">User ID: {sessionUser.id}</div>
                  </div>
                  <Button type="button" variant="outline" onClick={handleSignOut} disabled={busyAction === 'signout'} className="text-xs">
                    {busyAction === 'signout' ? <Loader2 className="w-3.5 h-3.5 mr-1 animate-spin" /> : <LogOut className="w-3.5 h-3.5 mr-1" />}
                    Sign Out
                  </Button>
                </div>

                <div className="flex flex-wrap gap-2">
                  <Button type="button" onClick={() => void handleCloudSave('new')} disabled={busyAction === 'save-new'} className="text-xs">
                    {busyAction === 'save-new' ? <Loader2 className="w-3.5 h-3.5 mr-1 animate-spin" /> : <Save className="w-3.5 h-3.5 mr-1" />}
                    Save New
                  </Button>
                  <Button
                    type="button"
                    variant="outline"
                    onClick={() => void handleCloudSave('update')}
                    disabled={!selectedCloudSaveId || busyAction === 'save-update'}
                    className="text-xs"
                  >
                    {busyAction === 'save-update' ? <Loader2 className="w-3.5 h-3.5 mr-1 animate-spin" /> : <Save className="w-3.5 h-3.5 mr-1" />}
                    Update Selected
                  </Button>
                  <Button type="button" variant="outline" onClick={() => void refreshCloudSaves()} disabled={busyAction === 'refresh'} className="text-xs">
                    {busyAction === 'refresh' ? <Loader2 className="w-3.5 h-3.5 mr-1 animate-spin" /> : <RefreshCw className="w-3.5 h-3.5 mr-1" />}
                    Refresh
                  </Button>
                </div>

                <div className="space-y-2">
                  {cloudSaves.length === 0 && (
                    <div className="rounded-md border border-border/70 p-3 text-xs text-muted-foreground">
                      No cloud saves yet.
                    </div>
                  )}
                  {cloudSaves.map((save) => (
                    <div
                      key={save.id}
                      className={`rounded-md border p-3 ${
                        save.id === selectedCloudSaveId ? 'border-blue-500 bg-blue-500/5' : 'border-border/70'
                      }`}
                    >
                      <div className="flex items-start justify-between gap-3">
                        <div className="min-w-0">
                          <button
                            type="button"
                            className="text-sm font-medium text-left hover:underline"
                            onClick={() => {
                              setSelectedCloudSaveId(save.id);
                              setSaveName(save.name);
                            }}
                          >
                            {save.name}
                          </button>
                          <div className="text-xs text-muted-foreground">
                            Updated {formatTimestamp(save.updated_at)}
                          </div>
                        </div>
                        <div className="flex gap-2">
                          <Button
                            type="button"
                            variant="outline"
                            size="sm"
                            className="h-8 px-2 text-xs"
                            onClick={() => void handleCloudLoad(save)}
                            disabled={busyAction === `load-${save.id}`}
                          >
                            {busyAction === `load-${save.id}` ? <Loader2 className="w-3.5 h-3.5 animate-spin" /> : 'Load'}
                          </Button>
                          <Button
                            type="button"
                            variant="outline"
                            size="sm"
                            className="h-8 px-2 text-xs text-red-600"
                            onClick={() => void handleCloudDelete(save)}
                            disabled={busyAction === `delete-${save.id}`}
                          >
                            {busyAction === `delete-${save.id}` ? <Loader2 className="w-3.5 h-3.5 animate-spin" /> : <Trash2 className="w-3.5 h-3.5" />}
                          </Button>
                        </div>
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}
          </section>

          {(statusMessage || errorMessage) && (
            <section className="space-y-2">
              {statusMessage && (
                <div className="rounded-md border border-green-500/40 bg-green-500/10 p-3 text-xs text-green-700 dark:text-green-300">
                  {statusMessage}
                </div>
              )}
              {errorMessage && (
                <div className="rounded-md border border-red-500/40 bg-red-500/10 p-3 text-xs text-red-700 dark:text-red-300">
                  {errorMessage}
                </div>
              )}
            </section>
          )}
        </div>
      </SheetContent>
    </Sheet>
  );
}
