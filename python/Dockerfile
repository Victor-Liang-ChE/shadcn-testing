FROM python:3.12-slim

WORKDIR /app

# Expose the port the app runs on
ENV PORT=8080

# Copy requirements file
COPY requirements.txt ./

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy your application code from the "python" folder into /app
COPY . ./

# Run the web service
CMD ["sh", "-c", "exec gunicorn --bind :$PORT --workers 1 --threads 8 --timeout 0 app:app"]

# TO TEST LOCALLY:
# docker build -t mccabe-thiele-calculator .
# docker run -p 8080:8080 mccabe-thiele-calculator
# curl http://localhost:8080/health