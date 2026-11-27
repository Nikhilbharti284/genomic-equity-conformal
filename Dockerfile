FROM python:3.11-slim
WORKDIR /app
COPY pyproject.toml README.md ./
COPY src/ ./src/
COPY data/ ./data/
COPY artifacts/ ./artifacts/
RUN pip install --no-cache-dir numpy scikit-learn pandas scipy psutil joblib && pip install --no-cache-dir -e . --no-deps --no-build-isolation
EXPOSE 8000
CMD ["gec-api"]
