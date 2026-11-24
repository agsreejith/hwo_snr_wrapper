# Makefile for hwo_snr_wrapper (POLLUX SNR calculator + NiceGUI front end)
#
# Expected layout (run `make` from project_root):
#
#   project_root/
#   |-- Makefile
#   |-- pyproject.toml           defines the `pollux_etc` command
#   |-- data/                    input data
#   |-- output/                  generated SNR text files and plots
#   `-- src/
#       |-- get_pollux_spectra.py
#       |-- pollux_snr_util.py
#       |-- Telescope.py
#       |-- planet_sim.py
#       |-- hwo_snr_wrapper.py   (command-line script)
#       |-- app3.py               (NiceGUI front end)
#       `-- pollux_etc.py        (launcher behind the `pollux_etc` command)
#
# Quick start:
#   make install      # create the petc venv, install deps and the pollux_etc command
#   make link         # (optional) make `pollux_etc` available from any terminal
#   pollux_etc        # launch the front end on localhost
#
# Override settings on the command line, e.g.:  make app PORT=8080
 
SHELL    := /bin/bash
 
# ---- Settings ---------------------------------------------------------------
PYTHON   ?= python3
VENV     ?= petc
HOST     ?= 127.0.0.1
PORT     ?= 8000
BIN_DIR  ?= $(HOME)/.local/bin
SRC_DIR  := src
OUT_DIR  := output
DATA_DIR := data
 
VENV_PY  := $(abspath $(VENV))/bin/python
VENV_CMD := $(abspath $(VENV))/bin/pollux_etc
STAMP    := $(VENV)/.installed
 
REQUIRED_SRC := get_pollux_spectra.py pollux_snr_util.py Telescope.py planet_sim.py
 
.DEFAULT_GOAL := help
.PHONY: help venv install link unlink check launch app app-open script dirs \
        clean clean-output distclean
 
# ---- Help -------------------------------------------------------------------
help:
	@echo "HWO SNR wrapper - available targets:"
	@echo ""
	@echo "  make install       Create $(VENV)/, install dependencies and the 'pollux_etc' command"
	@echo "  make link          Symlink 'pollux_etc' into $(BIN_DIR) so it works from any terminal"
	@echo "  make unlink        Remove that symlink"
	@echo "  make launch        Run 'pollux_etc' (opens http://$(HOST):$(PORT))"
	@echo "  make check         Verify source files, directories and PYSYN_CDBS"
	@echo "  make app           Run the NiceGUI app in dev mode with auto-reload"
	@echo "  make script        Run the original command-line script (hwo_snr_wrapper.py)"
	@echo "  make clean-output  Delete generated files in $(OUT_DIR)/"
	@echo "  make clean         Remove __pycache__ and build folders"
	@echo "  make distclean     clean + remove the virtual environment"
	@echo ""
	@echo "Options: PORT=$(PORT)  HOST=$(HOST)  PYTHON=$(PYTHON)  VENV=$(VENV)  BIN_DIR=$(BIN_DIR)"
 
# ---- Environment ------------------------------------------------------------
venv: $(VENV)/bin/python
 
$(VENV)/bin/python:
	$(PYTHON) -m venv $(VENV)
	$(VENV_PY) -m pip install --upgrade pip
 
# Editable install: dependencies come from pyproject.toml, and the `pollux_etc`
# command runs the code in src/ directly (edits take effect without reinstalling).
$(STAMP): $(VENV)/bin/python pyproject.toml
	$(VENV_PY) -m pip install -e .
	@touch $(STAMP)
 
install: $(STAMP) dirs
	@echo ""
	@echo "Installed. Command: $(VENV_CMD)"
	@echo "  - Activate the environment:   source $(VENV)/bin/activate   then run: pollux_etc"
	@echo "  - Or run 'make link' to use 'pollux_etc' from any terminal."
	@echo "stsynphot also needs its reference data: follow the stsynphot docs,"
	@echo "set PYSYN_CDBS, then run 'make check'."
 
link: $(STAMP)
	@mkdir -p $(BIN_DIR)
	ln -sf $(VENV_CMD) $(BIN_DIR)/pollux_etc
	@echo "Linked $(BIN_DIR)/pollux_etc -> $(VENV_CMD)"
	@case ":$$PATH:" in *":$(BIN_DIR):"*) ;; \
	  *) echo "NOTE: $(BIN_DIR) is not on your PATH. Add this to ~/.bashrc or ~/.zshrc:"; \
	     echo "      export PATH=\"$(BIN_DIR):\$$PATH\"" ;; esac
 
unlink:
	rm -f $(BIN_DIR)/pollux_etc
 
dirs:
	@mkdir -p $(DATA_DIR) $(OUT_DIR)
 
# ---- Sanity checks ----------------------------------------------------------
check: dirs
	@ok=1; \
	for f in $(REQUIRED_SRC); do \
	  if [ ! -f $(SRC_DIR)/$$f ]; then echo "MISSING: $(SRC_DIR)/$$f"; ok=0; fi; \
	done; \
	for f in app3.py pollux_etc.py; do \
	  if [ ! -f $(SRC_DIR)/$$f ]; then echo "MISSING: $(SRC_DIR)/$$f"; ok=0; fi; \
	done; \
	if [ ! -f $(SRC_DIR)/hwo_snr_wrapper.py ]; then echo "NOTE: $(SRC_DIR)/hwo_snr_wrapper.py not found (only needed for 'make script')"; fi; \
	if [ -z "$$PYSYN_CDBS" ]; then echo "WARNING: PYSYN_CDBS is not set (stsynphot reference data)"; \
	elif [ ! -d "$$PYSYN_CDBS" ]; then echo "WARNING: PYSYN_CDBS=$$PYSYN_CDBS is not a directory"; \
	else echo "PYSYN_CDBS=$$PYSYN_CDBS"; fi; \
	if [ $$ok -eq 1 ]; then echo "Required source files found."; else exit 1; fi
 
# ---- Run --------------------------------------------------------------------
launch: $(STAMP) check
	$(VENV_CMD) --host $(HOST) --port $(PORT)
 
# Dev mode: run app3.py directly with auto-reload (from src/ so local imports work)
app: $(STAMP) check
	cd $(SRC_DIR) && $(VENV_PY) app3.py --reload --no-browser --host $(HOST) --port $(PORT)
 
app-open: $(STAMP) check
	cd $(SRC_DIR) && $(VENV_PY) app3.py --reload --host $(HOST) --port $(PORT)
 
script: $(STAMP) check
	cd $(SRC_DIR) && $(VENV_PY) hwo_snr_wrapper.py
 
# ---- Cleanup ----------------------------------------------------------------
clean-output:
	rm -f $(OUT_DIR)/hwo_*.txt $(OUT_DIR)/hwo_*.png
 
clean:
	find . -name "__pycache__" -type d -prune -exec rm -rf {} +
	rm -rf build *.egg-info $(SRC_DIR)/*.egg-info
 
distclean: clean unlink
	rm -rf $(VENV)
