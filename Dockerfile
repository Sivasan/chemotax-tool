# Use the official Python slim image as a base
FROM python:3.11-slim

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    FLASK_ENV=production

# Set up working directory
WORKDIR /app

# Copy requirements and install Python packages
COPY requirements.txt .
RUN pip install -r requirements.txt

# Copy all application files into the /app directory
COPY backend_api.py .
COPY index.html .
# We DO NOT copy gunicorn.conf.py anymore

# Expose the port Gunicorn will run on (this is just metadata)
EXPOSE 10000

# Start Gunicorn
# This command directly uses the $PORT variable provided by Render.
# This is more reliable.
CMD gunicorn --bind 0.0.0.0:$PORT --workers 3 --threads 2 --timeout 180 backend_api:app
