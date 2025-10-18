"""
Command-line interface for launching the Streamlit GUI
"""

import subprocess
import sys
import os
from pathlib import Path


def main():
    """Main entry point for the CLI command."""
    print("🧬 Starting Codon Usage Analysis GUI...")

    # Get the path to the app.py file
    app_path = Path(__file__).parent / "app.py"

    # Launch Streamlit with the app
    try:
        cmd = [
            sys.executable, "-m", "streamlit", "run",
            str(app_path),
            "--server.headless", "true",
            "--browser.gatherUsageStats", "false"
        ]

        print(f"Launching Streamlit app at: {app_path}")
        print("The GUI will open in your default web browser...")

        subprocess.run(cmd)

    except KeyboardInterrupt:
        print("\n👋 Shutting down Codon Usage Analysis GUI...")
    except FileNotFoundError:
        print("❌ Error: Streamlit is not installed. Please install it with:")
        print("   pip install streamlit")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error launching GUI: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()