# ===== Imagen base / Base image =====
FROM python:3.13-slim

# ===== Instalar gfortran 14.2 y utilidades necesarias / Install gfortran 14.2 and required tools =====
RUN apt-get update && apt-get install -y \
    gfortran-14 \
    build-essential \
    curl \
    && update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-14 100 \
    && rm -rf /var/lib/apt/lists/*

# ===== Instalar Poetry / Install Poetry =====
RUN pip install --upgrade pip \
    && pip install poetry

# ===== Crear directorio de trabajo / Set working directory =====
WORKDIR /app

# ===== Copiar pyproject.toml primero (para aprovechar cache de Docker) / Copy pyproject.toml first to leverage Docker cache =====
COPY pyproject.toml .

# ===== Instalar dependencias de Python sin virtualenv y sin instalar el proyecto root / Install Python dependencies without virtualenv and without installing the root project =====
RUN poetry config virtualenvs.create false \
    && poetry install --no-root --no-interaction --no-ansi

# ===== Copiar resto del proyecto / Copy the rest of the project =====
COPY data/ ./data
COPY fortran/ ./fortran
COPY python/ ./python
COPY web/ ./web

# ===== Exponer puerto para Flask / Expose port for Flask =====
EXPOSE 5000

# ===== Comando por defecto / Default command =====
CMD ["python", "python/main.py"]
