#!/usr/bin/env python3
"""
Development setup script for codon-usage-gui

This script helps developers set up the development environment
and validates that everything is working correctly.
"""

import sys
import subprocess
import os
from pathlib import Path


def run_command(cmd, description, check=True):
    """Run a command and handle errors gracefully."""
    print(f"\n📦 {description}...")
    try:
        result = subprocess.run(cmd, shell=True, check=check, 
                              capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✅ {description} completed successfully")
            if result.stdout.strip():
                print(f"   Output: {result.stdout.strip()}")
        else:
            print(f"❌ {description} failed")
            if result.stderr.strip():
                print(f"   Error: {result.stderr.strip()}")
        return result.returncode == 0
    except Exception as e:
        print(f"❌ {description} failed with exception: {e}")
        return False


def check_git_repo():
    """Check if we're in a git repository."""
    return os.path.exists('.git')


def setup_development_environment():
    """Set up the development environment."""
    print("🧬 Codon Usage GUI - Development Setup")
    print("Setting up development environment...")
    
    # Check if we're in the right directory
    if not os.path.exists('pyproject.toml'):
        print("❌ pyproject.toml not found. Are you in the project root?")
        return False
    
    print("✅ Found pyproject.toml")
    
    # Install package in development mode
    success = True
    
    # Install development dependencies
    success &= run_command(
        "pip install -e .[dev,test]",
        "Installing package in development mode with dev dependencies"
    )
    
    # Install pre-commit hooks if in git repo
    if check_git_repo():
        success &= run_command(
            "pre-commit install",
            "Installing pre-commit hooks",
            check=False  # Don't fail if pre-commit not available
        )
    else:
        print("⚠️  Not in a git repository, skipping pre-commit setup")
    
    return success


def run_tests():
    """Run the test suite."""
    print("\n🧪 Running test suite...")
    
    # Run pytest with coverage
    success = run_command(
        "python -m pytest tests/ -v --cov=codon_usage_gui --cov-report=term-missing",
        "Running unit tests with coverage"
    )
    
    return success


def run_linting():
    """Run linting checks."""
    print("\n🔍 Running code quality checks...")
    
    success = True
    
    # Run flake8
    success &= run_command(
        "flake8 codon_usage_gui tests",
        "Running flake8 linting",
        check=False
    )
    
    # Run black check
    success &= run_command(
        "black --check codon_usage_gui tests",
        "Checking code formatting with black",
        check=False
    )
    
    # Run mypy
    success &= run_command(
        "mypy codon_usage_gui",
        "Running type checking with mypy",
        check=False
    )
    
    return success


def test_cli():
    """Test the CLI command."""
    print("\n💻 Testing CLI command...")
    
    success = run_command(
        "codon-usage-gui --help",
        "Testing CLI command availability",
        check=False
    )
    
    return success


def test_package_build():
    """Test package building."""
    print("\n📦 Testing package build...")
    
    # Clean previous builds
    run_command(
        "rm -rf dist/ build/ *.egg-info/",
        "Cleaning previous builds",
        check=False
    )
    
    # Build package
    success = run_command(
        "python -m build",
        "Building package distribution"
    )
    
    if success:
        # Check build artifacts
        dist_dir = Path("dist")
        if dist_dir.exists():
            files = list(dist_dir.glob("*"))
            print(f"✅ Built {len(files)} distribution files:")
            for file in files:
                print(f"   - {file.name}")
        
        # Validate build
        success &= run_command(
            "python -m twine check dist/*",
            "Validating built package"
        )
    
    return success


def main():
    """Main setup function."""
    print("🚀 Development Environment Setup for Codon Usage GUI")
    print("=" * 60)
    
    all_success = True
    
    # Setup development environment
    all_success &= setup_development_environment()
    
    # Run tests
    all_success &= run_tests()
    
    # Run linting
    linting_success = run_linting()
    if not linting_success:
        print("⚠️  Some linting checks failed (this is often normal during development)")
    
    # Test CLI
    all_success &= test_cli()
    
    # Test package building
    build_success = test_package_build()
    if not build_success:
        print("⚠️  Package build failed (check dependencies)")
    
    # Final summary
    print("\n" + "=" * 60)
    print("🏁 DEVELOPMENT SETUP SUMMARY")
    print("=" * 60)
    
    if all_success:
        print("🎉 Development environment setup completed successfully!")
        print("\nNext steps:")
        print("1. Start coding! The package is installed in development mode")
        print("2. Run tests frequently: make test")
        print("3. Check code quality: make lint")
        print("4. Format code: make format")
        print("5. Before committing: make test && make lint")
        
        print("\nUseful commands:")
        print("- make test          # Run tests")
        print("- make test-coverage # Run tests with coverage")
        print("- make lint          # Run linting")
        print("- make format        # Format code")
        print("- make clean         # Clean build artifacts")
        print("- make build         # Build package")
        
    else:
        print("❌ Some setup steps failed. Please check the errors above.")
        print("Common solutions:")
        print("- Make sure you have Python 3.8+ installed")
        print("- Update pip: python -m pip install --upgrade pip")
        print("- Install build tools: pip install build twine")
        
    return 0 if all_success else 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)