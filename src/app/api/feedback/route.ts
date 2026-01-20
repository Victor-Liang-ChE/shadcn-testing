import { NextRequest, NextResponse } from 'next/server';
import nodemailer from 'nodemailer';

// In-memory store for spam tracking (IP -> timestamp)
// In production, use Redis or a database
const spamTracker = new Map<string, number>();

// Clean up old entries every hour
setInterval(() => {
  const now = Date.now();
  const oneDay = 24 * 60 * 60 * 1000;
  for (const [ip, timestamp] of spamTracker.entries()) {
    if (now - timestamp > oneDay) {
      spamTracker.delete(ip);
    }
  }
}, 60 * 60 * 1000);

export async function POST(request: NextRequest) {
  try {
    const body = await request.json();
    const { feedbackType, feedbackName, feedbackEmail, feedbackMessage } = body;

    // Validate required fields
    if (!feedbackMessage || feedbackMessage.trim().length === 0) {
      return NextResponse.json(
        { error: 'Message is required.' },
        { status: 400 }
      );
    }

    // Get client IP for spam tracking
    const ip = request.headers.get('x-forwarded-for') || 
               request.headers.get('x-real-ip') || 
               'unknown';

    // Check if this IP has already submitted feedback in the last 24 hours
    if (spamTracker.has(ip)) {
      const lastSubmission = spamTracker.get(ip)!;
      const now = Date.now();
      const timeSinceLastSubmission = now - lastSubmission;
      const oneDay = 24 * 60 * 60 * 1000;

      if (timeSinceLastSubmission < oneDay) {
        const hoursRemaining = Math.ceil((oneDay - timeSinceLastSubmission) / (60 * 60 * 1000));
        return NextResponse.json(
          { 
            error: `You can submit feedback again in ${hoursRemaining} hour${hoursRemaining > 1 ? 's' : ''}. Please try again later.` 
          },
          { status: 429 }
        );
      }
    }

    // Update spam tracker with current timestamp
    spamTracker.set(ip, Date.now());

    // Configure nodemailer
    const transporter = nodemailer.createTransport({
      service: 'gmail',
      auth: {
        user: process.env.EMAIL_USER,
        pass: process.env.EMAIL_PASSWORD,
      },
    });

    // Prepare email content
    const subject = `${feedbackType === 'feedback' ? 'Feedback' : 'Bug Report'} Submission`;
    const emailBody = `
Type: ${feedbackType === 'feedback' ? 'Feedback' : 'Bug Report'}
Name: ${feedbackName || 'Not provided'}
Email: ${feedbackEmail || 'Not provided'}
Message: ${feedbackMessage}
    `.trim();

    // Send email
    await transporter.sendMail({
      from: process.env.EMAIL_USER,
      to: 'victorl1725@gmail.com',
      subject,
      text: emailBody,
      replyTo: feedbackEmail || process.env.EMAIL_USER,
    });

    return NextResponse.json(
      { message: 'Feedback submitted successfully.' },
      { status: 200 }
    );
  } catch (error) {
    console.error('Error processing feedback:', error);
    return NextResponse.json(
      { error: 'An error occurred while processing your feedback.' },
      { status: 500 }
    );
  }
}
