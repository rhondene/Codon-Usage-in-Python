# Makefile for codon-usage-gui development

.PHONY: help install install-dev test test-coverage lint format clean build upload

help:  ## Show this help message
	@echo "Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

install:  ## Install the package
	pip install -e .

install-dev:  ## Install the package with development dependencies
	pip install -e ".[dev,test]"

test:  ## Run tests
	pytest

test-coverage:  ## Run tests with coverage report
	pytest --cov=codon_usage_gui --cov-report=html --cov-report=term

lint:  ## Run linting checks
	flake8 codon_usage_gui tests
	mypy codon_usage_gui

format:  ## Format code with black
	black codon_usage_gui tests

format-check:  ## Check if code is formatted correctly
	black --check codon_usage_gui tests

clean:  ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

build:  ## Build the package
	python -m build

check-build:  ## Check the built package
	python -m twine check dist/*

upload-test:  ## Upload to Test PyPI
	python -m twine upload --repository testpypi dist/*

upload:  ## Upload to PyPI
	python -m twine upload dist/*

run-gui:  ## Run the GUI application
	codon-usage-gui

run-dev:  ## Run the GUI application in development mode
	python -m codon_usage_gui.app

validate-install:  ## Validate package installation
	python scripts/validate_installation.py