.PHONY: install test verify figures smoke data-help
install:
	python -m pip install -e '.[test]'
test:
	python -m pytest -q
verify:
	python scripts/validate_frozen.py
figures:
	python scripts/rebuild_figures_from_frozen.py
smoke:
	python scripts/run_rank_calibration.py --smoke

data-help:
	python scripts/fetch_sources.py --help
