'use client'

import React, { useState } from 'react'
import * as DialogPrimitive from '@radix-ui/react-dialog'
import { cn } from '@/lib/utils'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Label } from '@/components/ui/label'
import { Textarea } from '@/components/ui/textarea'

const Dialog = DialogPrimitive.Root;
const DialogTrigger = DialogPrimitive.Trigger;
const DialogPortal = DialogPrimitive.Portal;


const DialogOverlay = (({ className, ...props }: React.ComponentProps<typeof DialogPrimitive.Overlay>) => (
  <DialogPrimitive.Overlay
    className={cn(
      "data-[state=open]:animate-in data-[state=closed]:animate-out data-[state=closed]:fade-out-0 data-[state=open]:fade-in-0 fixed inset-0 z-50 bg-black/80 backdrop-blur-sm",
      className
    )}
    {...props}
  />
));

const DialogContent = (({ className, children, ...props }: React.ComponentProps<typeof DialogPrimitive.Content>) => (
  <DialogPortal>
    <DialogOverlay />
    <DialogPrimitive.Content
      className={cn(
        "bg-background data-[state=open]:animate-in data-[state=closed]:animate-out data-[state=closed]:fade-out-0 data-[state=open]:fade-in-0 fixed left-[50%] top-[50%] z-50 grid w-full max-w-lg translate-x-[-50%] translate-y-[-50%] gap-4 rounded-2xl border p-6 shadow-lg duration-200",
        className
      )}
      {...props}
    >
      {children}
    </DialogPrimitive.Content>
  </DialogPortal>
));

const DialogHeader = (({ className, ...props }: React.HTMLAttributes<HTMLDivElement>) => (
  <div className={cn("flex flex-col space-y-1.5 text-center sm:text-left", className)} {...props} />
));
const DialogFooter = (({ className, ...props }: React.HTMLAttributes<HTMLDivElement>) => (
  <div className={cn("flex flex-col-reverse sm:flex-row sm:justify-end sm:space-x-2", className)} {...props} />
));
const DialogTitle = (({ className, ...props }: React.ComponentProps<typeof DialogPrimitive.Title>) => (
  <DialogPrimitive.Title className={cn("text-lg font-semibold leading-none tracking-tight", className)} {...props} />
));
const DialogDescription = (({ className, ...props }: React.ComponentProps<typeof DialogPrimitive.Description>) => (
  <DialogPrimitive.Description className={cn("text-muted-foreground text-sm", className)} {...props} />
));

export default function FeedbackDialog({ triggerLabel = 'Submit Feedback' }: { triggerLabel?: string }) {
  const [open, setOpen] = useState(false);
  const [feedbackType, setFeedbackType] = useState('feedback');
  const [feedbackName, setFeedbackName] = useState('');
  const [feedbackEmail, setFeedbackEmail] = useState('');
  const [feedbackMessage, setFeedbackMessage] = useState('');
  const [feedbackError, setFeedbackError] = useState<string | null>(null);
  const [isSubmitting, setIsSubmitting] = useState(false);

  const handleSubmit = async () => {
    setFeedbackError(null);
    setIsSubmitting(true);
    try {
      const response = await fetch('/api/feedback', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ feedbackType, feedbackName, feedbackEmail, feedbackMessage }),
      });
      if (response.ok) {
        alert('Feedback submitted! Thank you.');
        setOpen(false);
        setFeedbackType('feedback');
        setFeedbackName('');
        setFeedbackEmail('');
        setFeedbackMessage('');
      } else {
        const errorData = await response.json();
        setFeedbackError(errorData.error || 'Failed to submit. Please try again.');
      }
    } catch (error) {
      console.error('Error submitting feedback:', error);
      setFeedbackError('Failed to submit. Please try again.');
    } finally {
      setIsSubmitting(false);
    }
  };

  return (
    <Dialog open={open} onOpenChange={setOpen}>
      <DialogTrigger asChild>
        <Button variant="outline">{triggerLabel}</Button>
      </DialogTrigger>
      <DialogContent>
        <DialogHeader>
          <DialogTitle>Feedback and Bug Reports</DialogTitle>
          <DialogDescription>I know I have bugs to fix and features to add and change!</DialogDescription>
        </DialogHeader>

        <div className="grid gap-4 py-4">
          <div className="grid gap-2">
            <Label htmlFor="fb-name">Name (Optional)</Label>
            <Input id="fb-name" value={feedbackName} onChange={(e) => setFeedbackName(e.target.value)} placeholder="Your name" />
          </div>
          <div className="grid gap-2">
            <Label htmlFor="fb-email">Email (Optional)</Label>
            <Input id="fb-email" type="email" value={feedbackEmail} onChange={(e) => setFeedbackEmail(e.target.value)} placeholder="your@email.com" />
          </div>
          <div className="grid gap-2">
            <Label htmlFor="fb-message">Message</Label>
            <Textarea id="fb-message" value={feedbackMessage} onChange={(e) => setFeedbackMessage(e.target.value)} placeholder="Describe your feedback or bug..." className="min-h-[120px] resize-none" />
          </div>
          {feedbackError && <div className="text-sm text-red-500">{feedbackError}</div>}
        </div>

        <DialogFooter>
          <Button onClick={handleSubmit} disabled={isSubmitting}>{isSubmitting ? 'Submitting...' : 'Submit'}</Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  );
}
