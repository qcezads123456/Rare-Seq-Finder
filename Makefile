ifeq ($(firstword $(MAKECMDGOALS)),run)
RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
.PHONY: $(RUN_ARGS)
$(eval $(RUN_ARGS):;@:)
endif

ifeq ($(firstword $(MAKECMDGOALS)),generate_test_data)
GENERATE_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
.PHONY: $(GENERATE_ARGS)
$(eval $(GENERATE_ARGS):;@:)
endif


PYTHON := python3.12
VENV := .venv
VENV_PY := $(VENV)/bin/python
VENV_MS := $(VENV)/bin/mutation-simulator
MS_OUTPUT := data/MS_output
REF_DIR := data/ref
SELECTED_DIR := data/selected/selected.fasta
DEFAULT_MS_FASTA = $(firstword $(wildcard $(MS_OUTPUT)/*_ms.fasta))
DEFAULT_REF_FASTA = $(firstword $(wildcard $(REF_DIR)/ref.fasta))
.PHONY: run build ms_setup clean generate_test_data

$(VENV_PY):
	@$(PYTHON) -m venv $(VENV)

build : $(build/CMakeCache.txt)
	@cmake -S . -B build
$(build/CMakeCache.txt) : $(CmakeLists.txt)
	@cmake --build build

run : build
	@if [ -z "$(filter --reference -r,$(RUN_ARGS))" ] && [ -z "$(DEFAULT_REF_FASTA)" ]; then \
		echo "ERROR: ref.fasta not found"; \
		exit 1; \
	fi
	@if [ -z "$(filter --mutant -m,$(RUN_ARGS))" ] && [ -z "$(DEFAULT_MS_FASTA)" ]; then \
		echo "ERROR: *_ms.fasta not found"; \
		exit 1; \
	fi
	@RUN_VIA_MAKE="make run --" ./build/bin/rare-seq-finder \
	$(RUN_ARGS) \
	$(if $(filter -r,$(RUN_ARGS)),,-r $(DEFAULT_REF_FASTA)) \
	$(if $(filter -m,$(RUN_ARGS)),,-m $(DEFAULT_MS_FASTA))

ms_setup : $(VENV_MS)

$(VENV_MS) : $(VENV_PY)
	@$(VENV_PY) -m pip install Mutation-Simulator

generate_test_data : ms_setup
	@$(VENV_PY) scripts/select_seq.py $(GENERATE_ARGS)
	@$(VENV_MS) -o $(MS_OUTPUT)/mutant.fasta $(SELECTED_DIR) rmt config/mutation.rmt
	@echo "mutation simulation complete"
	
clean :
	@rm -rf $(VENV) build/*