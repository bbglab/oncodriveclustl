ROOT_DIR := $(shell echo $(dir $(lastword $(MAKEFILE_LIST))) | sed 's|/*$$||')

SHELL := /bin/bash

define version
$(shell uv run python -c "from oncodriveclustl import __version__; print(__version__)")
endef

define git_tag_or_sha
$(shell git describe --tags --exact-match 2>/dev/null || git rev-parse --short HEAD)
endef

define image
bbglab/oncodriveclustl:$(call version)
endef

BOLDRED := $(shell tput bold && tput setaf 1)
BOLDGREEN := $(shell tput bold && tput setaf 2)
BOLDYELLOW := $(shell tput bold && tput setaf 3)
BOLDBLUE := $(shell tput bold && tput setaf 4)
LIGHTBLUE := $(shell tput setaf 6)
WHITE := $(shell tput sgr0 && tput setaf 7)
RESET := $(shell tput sgr0)


.PHONY: help
help:
	@echo "$(BOLDYELLOW)Available targets:$(RESET)"
	@echo
	@echo "$(BOLDGREEN)  checks       $(WHITE)-> Run all the checks (format and lint)"
	@echo "$(BOLDGREEN)  check-format $(WHITE)-> Check for formatting errors"
	@echo "$(BOLDGREEN)  check-lint   $(WHITE)-> Check for lint errors"
	@echo "$(BOLDGREEN)  check-docker $(WHITE)-> Check the Dockerfile"
	@echo "$(BOLDGREEN)  format       $(WHITE)-> Format source code"
	@echo "$(BOLDGREEN)  build-dist   $(WHITE)-> Build source and wheel distribution files"
	@echo "$(BOLDGREEN)  build-image  $(WHITE)-> Build the Docker image"
	@echo "$(BOLDGREEN)  push-image   $(WHITE)-> Push the Docker image into DockerHub"
	@echo "$(BOLDGREEN)  run-example  $(WHITE)-> Run the included example using the Docker image"
	@echo "$(BOLDGREEN)  clean        $(WHITE)-> Clean the working directory (build files, virtual environments, caches)"
	@echo "$(RESET)"

.PHONY: uv-installed
uv-installed:
	@if ! which uv > /dev/null; then \
		echo "$(BOLDRED)This project build is managed by $(BOLDYELLOW)uv$(BOLDRED), which is not installed.$(RESET)"; \
		echo "$(LIGHTBLUE)Please follow these instructions to install it:$(RESET)"; \
		echo "$(LIGHTBLUE)--> $(BOLDBLUE)https://docs.astral.sh/uv/#getting-started$(RESET)"; \
		exit 1; \
	fi

.PHONY: ruff-installed
ruff-installed: uv-installed
	@if ! which ruff > /dev/null; then \
		echo "$(BOLDRED)This project requires $(BOLDYELLOW)ruff$(BOLDRED), which is not installed.$(RESET)"; \
		echo "$(LIGHTBLUE)Installing it with $(BOLDYELLOW)uv tool install ruff$(RESET)"; \
		uv tool install ruff; \
		ruff --version; \
	fi

.PHONY: checks
checks: check-format check-lint check-docker

.PHONY: check-format
check-format: ruff-installed
	@echo "$(BOLDGREEN)Checking code format ...$(RESET)"
	ruff format --check
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: check-lint
check-lint: ruff-installed
	@echo "$(BOLDGREEN)Checking lint ...$(RESET)"
	ruff check
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: check-docker
check-docker:
	@echo "$(BOLDGREEN)Checking Dockerfile ...$(RESET)"
	docker run --rm -i \
		-v $$(pwd):/project \
		hadolint/hadolint hadolint \
		--config /project/.hadolint.yaml \
		/project/Dockerfile
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: check-version
check-version: uv-installed
	@echo "$(BOLDGREEN)Checking that the version matches the tag ...$(RESET)"
	@if [ "$(call version)" != "$(call git_tag_or_sha)" ]; then \
	    echo "$(BOLDRED)==> Version $(BOLDYELLOW)$(call version)$(BOLDRED) doesn't match the git tag $(BOLDYELLOW)$(call git_tag_or_sha)$(BOLDRED) !!!$(RESET)"; \
		echo "$(BOLDRED)==> Please update the $(BOLDYELLOW)__version__$(BOLDRED) in $(BOLDYELLOW)oncodrivefml/__init__.py$(BOLDRED) and re-create the tag.$(RESET)"; \
	    exit 1; \
	fi
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: format
format: ruff-installed
	@echo "$(BOLDGREEN)Formatting code ...$(RESET)"
	ruff format

.PHONY: build-dist
build-dist: uv-installed
	@echo "$(BOLDGREEN)Building packages ...$(RESET)"
	uv build

.PHONY: publish-dist
publish-dist: uv-installed
	@echo "$(BOLDGREEN)Publishing OncodriveCLUSTL $(BOLDYELLOW)$(call version)$(BOLDGREEN) to PyPI ...$(RESET)"
	@if [ -z "$(PYPI_TOKEN)" ]; then \
		echo "$(BOLDRED)==> Missing PyPI token !!!$(RESET)"; \
		exit 1; \
	fi
	uv publish --token $(PYPI_TOKEN)

.PHONY: build-image
build-image: uv-installed
	@echo "$(BOLDGREEN)Building Docker image $(BOLDYELLOW)$(call image)$(BOLDGREEN) ...$(RESET)"
	docker build --progress=plain -t $(call image) .
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: save-image
save-image: build-image
	@echo "$(BOLDGREEN)Saving Docker image $(BOLDYELLOW)$(call image)$(BOLDGREEN) ...$(RESET)"
	docker save -o oncodriveclustl.tar $(call image)
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: load-image
load-image: uv-installed
	@echo "$(BOLDGREEN)Loading Docker image $(BOLDYELLOW)$(call image)$(BOLDGREEN) ...$(RESET)"
	docker load -i oncodriveclustl.tar
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: build-image
push-image: uv-installed
	@echo "$(BOLDGREEN)Pushing the Docker image into the DockerHub ...$(RESET)"
	docker push $(call image)
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: run-example
run-example: uv-installed
	@echo "$(BOLDGREEN)Running example ...$(RESET)"
	uv run oncodriveclustl \
		-i example/PAAD.tsv.gz -r example/cds.hg19.regions.gz -o example/output \
		-sw 15 -cw 15 -simw 35 -sim region_restricted --concatenate --clustplot -e KRAS
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: clean
clean:
	@echo "$(BOLDGREEN)Cleaning the repository ...$(RESET)"
	rm -rf ./oncodriveclustl.egg-info ./dist ./.ruff_cache ./.venv oncodriveclustl.tar
	find . -name "__pycache__" -type d -exec rm -r {} +
