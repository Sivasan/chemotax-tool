# Use the official Python slim image as a base
FROM python:3.11-slim

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    FLASK_ENV=production

# Install system-level dependencies for document conversion
# This is the new part. It installs Pandoc and the LaTeX engine for PDFs.
# This step will take several minutes during deployment.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    pandoc \
    texlive-latex-base \
    texlive-fonts-recommended \
    lmodern \
    && rm -rf /var/lib/apt/lists/*

# Set up working directory
WORKDIR /app

# Copy requirements and install Python packages
COPY requirements.txt .
RUN pip install -r requirements.txt

# Copy all application files into the /app directory
COPY backend_api.py .
COPY index.html .

# Expose the port Gunicorn will run on
EXPOSE 10000

# Start Gunicorn
CMD gunicorn --bind 0.0.0.0:$PORT --workers 3 --threads 2 --timeout 180 backend_api:app
