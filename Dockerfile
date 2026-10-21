FROM python:3.11-slim
WORKDIR /workspace
COPY pyproject.toml README.md LICENSE ./
COPY src/ ./src/
RUN pip install --no-cache-dir -e '.[test]'
COPY scripts/ ./scripts/
COPY tests/ ./tests/
COPY config/ ./config/
COPY expected/ ./expected/
COPY figures/reference/ ./figures/reference/
ENTRYPOINT ["python", "scripts/run_all.py"]
