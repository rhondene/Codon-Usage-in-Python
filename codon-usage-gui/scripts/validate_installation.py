#!/usr/bin/env python3
"""
Installation validation script for codon-usage-gui package

This script validates that the package is correctly installed and all
dependencies are working properly. It's designed to be user-friendly
for biologists with limited technical experience.
"""

import sys
import subprocess
import importlib
import tempfile
import os
from pathlib import Path


def print_section(title):
    """Print a formatted section header."""
    print(f"\n{'='*60}")
    print(f" {title}")
    print('='*60)


def check_python_version():
    """Check if Python version is compatible."""
    print_section("Python Version Check")
    
    version = sys.version_info
    print(f"Python version: {version.major}.{version.minor}.{version.micro}")
    
    if version >= (3, 8):
        print("✅ Python version is compatible")
        return True
    else:
        print("❌ Python version is too old. Requires Python 3.8+")
        return False


def check_package_import():
    """Check if the package can be imported."""
    print_section("Package Import Check")
    
    try:
        import codon_usage_gui
        print(f"✅ codon_usage_gui imported successfully")
        print(f"   Version: {codon_usage_gui.__version__}")
        print(f"   Author: {codon_usage_gui.__author__}")
        return True
    except ImportError as e:
        print(f"❌ Failed to import codon_usage_gui: {e}")
        print("   Try installing with: pip install codon-usage-gui")
        return False


def check_dependencies():
    """Check if all required dependencies are available."""
    print_section("Dependency Check")
    
    required_packages = [
        'pandas',
        'numpy', 
        'streamlit',
        'plotly',
        'matplotlib',
        'seaborn'
    ]
    
    all_good = True
    
    for package in required_packages:
        try:
            module = importlib.import_module(package)
            version = getattr(module, '__version__', 'unknown')
            print(f"✅ {package} ({version})")
        except ImportError:
            print(f"❌ {package} - not installed")
            all_good = False
    
    return all_good


def check_core_functions():
    """Check if core analysis functions work correctly."""
    print_section("Core Functions Check")
    
    try:
        from codon_usage_gui.core.analysis import (
            parse_fasta_from_text,
            compute_codon_frequencies,
            compute_rscu_weights
        )
        
        # Test with simple data
        test_fasta = """>test_gene
ATGGCAAGCTAG"""
        
        headers, seqs = parse_fasta_from_text(test_fasta)
        print("✅ FASTA parsing works")
        
        df_codons, skipped = compute_codon_frequencies(headers, seqs)
        print("✅ Codon frequency calculation works")
        
        df_rscu = compute_rscu_weights(df_codons)
        print("✅ RSCU calculation works")
        
        return True
        
    except Exception as e:
        print(f"❌ Core functions test failed: {e}")
        return False


def check_cli_command():
    """Check if the CLI command is available."""
    print_section("Command Line Interface Check")
    
    try:
        result = subprocess.run(
            ['codon-usage-gui', '--help'],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode == 0:
            print("✅ CLI command 'codon-usage-gui' is available")
            return True
        else:
            print("❌ CLI command failed")
            print(f"   Error: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print("❌ CLI command timed out")
        return False
    except FileNotFoundError:
        print("❌ CLI command 'codon-usage-gui' not found")
        print("   Try reinstalling the package")
        return False
    except Exception as e:
        print(f"❌ CLI test failed: {e}")
        return False


def check_example_files():
    """Check if example files are accessible."""
    print_section("Example Files Check")
    
    try:
        # Try to find example files
        current_dir = Path(__file__).parent.parent
        examples_dir = current_dir / "examples"
        
        if examples_dir.exists():
            example_files = list(examples_dir.glob("*.fasta"))
            if example_files:
                print(f"✅ Found {len(example_files)} example FASTA files")
                for file in example_files:
                    print(f"   - {file.name}")
                return True
            else:
                print("⚠️  Example directory exists but no FASTA files found")
                return False
        else:
            print("⚠️  Example files not found (this is okay for pip installations)")
            return True
            
    except Exception as e:
        print(f"⚠️  Could not check example files: {e}")
        return True  # Not critical


def run_integration_test():
    """Run a complete integration test."""
    print_section("Integration Test")
    
    try:
        from codon_usage_gui import (
            parse_fasta_from_text,
            compute_codon_frequencies,
            compute_rscu_weights,
            compute_amino_acid_usage
        )
        
        # Create test data
        test_fasta = """>gene1
ATGGCAAGCGAATTTGCCGAAGCCCTGGACAAAGCGAAATTGCTGGAATAG
>gene2
ATGGCTCTGGAAATTGCAGAGGCTCTGGATAAAACTAAACTCCTGGAATAG
>gene3
ATGGCTCTTGAAATTGCAGAAGCACTTGATAAAACAAAACTTCTGGAATAG"""
        
        print("Running complete analysis pipeline...")
        
        # Parse sequences
        headers, seqs = parse_fasta_from_text(test_fasta)
        print(f"  ✅ Parsed {len(seqs)} sequences")
        
        # Compute frequencies
        df_codons, skipped = compute_codon_frequencies(headers, seqs)
        total_codons = df_codons['Obs_Freq'].sum()
        print(f"  ✅ Analyzed {total_codons} codons")
        
        # Compute RSCU
        df_rscu = compute_rscu_weights(df_codons)
        print(f"  ✅ Computed RSCU for {len(df_rscu)} codons")
        
        # Compute amino acid usage
        df_aa = compute_amino_acid_usage(df_rscu)
        print(f"  ✅ Analyzed {len(df_aa)} amino acids")
        
        print("✅ Integration test passed!")
        return True
        
    except Exception as e:
        print(f"❌ Integration test failed: {e}")
        return False


def print_installation_summary(results):
    """Print a summary of the installation validation."""
    print_section("Installation Summary")
    
    passed = sum(results.values())
    total = len(results)
    
    print(f"Tests passed: {passed}/{total}")
    print()
    
    for test_name, result in results.items():
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status} {test_name}")
    
    print()
    
    if passed == total:
        print("🎉 All tests passed! Your installation is working correctly.")
        print("   You can now use codon-usage-gui for your research.")
        print()
        print("To get started:")
        print("  1. Run 'codon-usage-gui' to launch the GUI")
        print("  2. Check the examples/ directory for sample data")
        print("  3. Read the documentation in README.md")
        
    elif passed >= total - 1:
        print("⚠️  Most tests passed. Your installation should work.")
        print("   Some minor features might not be available.")
        
    else:
        print("❌ Several tests failed. Your installation may have issues.")
        print("   Try reinstalling the package:")
        print("     pip uninstall codon-usage-gui")
        print("     pip install codon-usage-gui")


def main():
    """Run all validation tests."""
    print("🧬 Codon Usage GUI - Installation Validation")
    print("This script will test your installation to ensure everything works correctly.")
    
    # Run all tests
    results = {
        "Python Version": check_python_version(),
        "Package Import": check_package_import(),
        "Dependencies": check_dependencies(),
        "Core Functions": check_core_functions(),
        "CLI Command": check_cli_command(),
        "Example Files": check_example_files(),
        "Integration Test": run_integration_test()
    }
    
    # Print summary
    print_installation_summary(results)
    
    # Return appropriate exit code
    all_critical_passed = all([
        results["Python Version"],
        results["Package Import"], 
        results["Dependencies"],
        results["Core Functions"],
        results["Integration Test"]
    ])
    
    return 0 if all_critical_passed else 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)