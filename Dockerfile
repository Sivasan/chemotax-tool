# ---------- Corrected & Production-Ready Dockerfile ----------
FROM python:3.11-slim

# Avoid interactive prompts during apt installs
ENV DEBIAN_FRONTEND=noninteractive

# App/runtime defaults
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    FLASK_ENV=production \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PORT=10000

# System deps for Pandoc → PDF (Unicode/Tamil) and fonts
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      pandoc \
      texlive-latex-base \
      texlive-latex-recommended \
      texlive-latex-extra \
      texlive-plain-generic \
      texlive-xetex \
      texlive-fonts-recommended \
      texlive-fonts-extra \
      lmodern \
      fontconfig \
      fonts-noto-core \
      fonts-noto-cjk \
      fonts-noto-color-emoji \
      fonts-dejavu-core \
      fonts-lohit-tamil \
      curl \
    && fc-cache -f \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user
RUN useradd -m -u 10001 appuser

# Workdir
WORKDIR /app

# Install Python deps first for better layer caching
COPY requirements.txt .
RUN python -m pip install --upgrade pip wheel && \
    pip install --no-cache-dir -r requirements.txt

# Copy application (backend_api.py, index.html, etc.)
COPY . .

# Drop privileges
USER appuser

# Document the listening port
EXPOSE 10000

# Optional healthcheck (uncomment if you expose /healthz)
# HEALTHCHECK --interval=30s --timeout=5s --retries=3 \
#   CMD curl -fsS http://127.0.0.1:${PORT}/healthz || exit 1

# Start Gunicorn
CMD gunicorn --bind 0.0.0.0:${PORT} --workers 3 --threads 2 --timeout 180 backend_api:app
