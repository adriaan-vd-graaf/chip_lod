# chip_lod build. `make all` = python + r + compare + report.
# In the Lima VM set UV_PROJECT_ENVIRONMENT=$$HOME/.venvs/chip_lod so the host .venv is untouched.
UV      ?= uv
PY      ?= $(UV) run python
PARAMS  ?= config/params.json
RSCRIPT ?= Rscript

.PHONY: all python python-test r r-test compare report clean

all: python r compare report

# python must run before r: the R cross-language test compares against results/python

python-test:
	$(UV) run pytest -q

python: python-test
	$(PY) -m chip_lod run --params $(PARAMS) --out results/python
	$(PY) -m chip_lod figures --out results/python --figdir results/figures

r-test:
	cd r && $(RSCRIPT) -e 'testthat::test_dir("tests/testthat", stop_on_failure = TRUE)'

# results first, then the R tests (which include the cross-language byte comparison)
r:
	$(RSCRIPT) r/run_all.R $(PARAMS) results/r
	$(MAKE) r-test

compare:
	bash scripts/compare_outputs.sh

report:
	typst compile --root . report/chip_lod.typ report/chip_lod.pdf

clean:
	rm -rf results/python/* results/r/* results/figures/* results/hash_check.json report/chip_lod.pdf
