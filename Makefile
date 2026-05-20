.PHONY: install dev lint test ingest dbt-run dbt-test export clean

# ── Environment ───────────────────────────────────────────────────────────────
install:
	pip install -e ".[dev]"

dev: install

# ── Quality ───────────────────────────────────────────────────────────────────
lint:
	ruff check src/ ingest/ tests/
	mypy src/ ingest/

test:
	pytest tests/ --cov=src --cov=ingest --cov-report=term-missing

# ── Pipeline ──────────────────────────────────────────────────────────────────
ingest-transplants:
	python -m ingest.transplant_artifacts

ingest-pgx:
	python -m ingest.pgx_artifacts

ingest: ingest-transplants ingest-pgx

# ── dbt ───────────────────────────────────────────────────────────────────────
dbt-run:
	cd dbt_project && dbt run --profiles-dir .

dbt-test:
	cd dbt_project && dbt test --profiles-dir .

dbt-docs:
	cd dbt_project && dbt docs generate --profiles-dir . && dbt docs serve --profiles-dir .

# ── Export ────────────────────────────────────────────────────────────────────
export:
	python -m src.latam_equity.export

# ── Full pipeline ─────────────────────────────────────────────────────────────
pipeline: ingest dbt-run dbt-test export

# ── Cleanup ───────────────────────────────────────────────────────────────────
clean:
	rm -rf data/raw/*.json data/exports/*.json artifacts/*.json
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
