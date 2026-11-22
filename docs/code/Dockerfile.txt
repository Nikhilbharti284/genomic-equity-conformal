FROM python:3.11-slim
WORKDIR /app
COPY pyproject.toml README.md ./
COPY src/ ./src/
COPY data/ ./data/
COPY artifacts/equity_model_v10_clean.joblib ./artifacts/
COPY artifacts/v3_offline_clinical_ref.sqlite ./artifacts/
RUN pip install --no-cache-dir numpy scikit-learn joblib && pip install --no-cache-dir -e . --no-deps --no-build-isolation
EXPOSE 8000
CMD ["gec-api"]
