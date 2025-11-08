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
COPY gunicorn.conf.py .
COPY index.html .

# Expose the port Gunicorn will run on
EXPOSE 5000

# Start Gunicorn
# This is the command that a hosting service will use.
CMD ["gunicorn", "-c", "gunicorn.conf.py", "backend_api:app"]

