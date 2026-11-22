FROM python:3.11-slim
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
WORKDIR /app
COPY requirements.txt .
RUN python3 -m pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu \
 && python3 -m pip install --no-cache-dir -r requirements.txt
COPY data/ba1_exception_registry.json data/
COPY artifacts/equity_model_v10_clean.joblib artifacts/
COPY src/ src/
COPY CITATION.cff LICENSE ro-crate-metadata.json ./
EXPOSE 8080
# Run air-gapped in hospital VPC: docker run --network=none -p 8080:8080 genomic-equity:1.1.0
CMD ["python3", "src/api_server.py"]
