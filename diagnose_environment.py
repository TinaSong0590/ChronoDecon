#!/usr/bin/env python3
"""
Environment diagnostic script for ChronoDecon Dashboard
Run this script to check if all dependencies are properly installed.
"""

import sys
import os

def print_header(title):
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70)

def check_python():
    print_header("Python Environment")
    print(f"Python Version: {sys.version}")
    print(f"Python Executable: {sys.executable}")
    print(f"Working Directory: {os.getcwd()}")

def check_package(package_name, import_name=None):
    """Check if a package is installed and importable."""
    if import_name is None:
        import_name = package_name

    try:
        module = __import__(import_name)
        version = getattr(module, '__version__', 'N/A')
        location = getattr(module, '__file__', 'N/A')
        print(f"✅ {package_name:20} | Version: {version:10} | {location}")
        return True
    except ImportError:
        print(f"❌ {package_name:20} | NOT INSTALLED")
        return False

def check_matchms_details():
    print_header("matchms Detailed Check")
    try:
        import matchms
        print(f"✅ matchms version: {matchms.__version__}")
        print(f"   Location: {matchms.__file__}")

        submodules = [
            "matchms.importing.load_from_mgf",
            "matchms.filtering.normalize_intensities",
            "matchms.filtering.default_filters",
            "matchms.similarity.CosineGreedy",
            "matchms.Spectrum.Spectrum",
        ]

        print("\n   Checking matchms submodules:")
        for submodule in submodules:
            module_path, attr = submodule.rsplit(".", 1)
            try:
                module = __import__(module_path, fromlist=[attr])
                getattr(module, attr)
                print(f"   ✅ {attr}")
            except (ImportError, AttributeError) as e:
                print(f"   ❌ {attr} - {e}")
        return True
    except ImportError as e:
        print(f"❌ matchms not installed: {e}")
        return False

def check_chronodecon():
    print_header("ChronoDecon Modules Check")

    # Add current directory to path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    if current_dir not in sys.path:
        sys.path.insert(0, current_dir)

    modules = [
        ("decon", None),
        ("library_search", None),
        ("decon.deconvolute_mzml", None),
        ("library_search.run_local_search", None),
    ]

    all_ok = True
    for module_spec, _ in modules:
        try:
            __import__(module_spec)
            print(f"✅ {module_spec}")
        except ImportError as e:
            print(f"❌ {module_spec} - {e}")
            all_ok = False

    return all_ok

def main():
    print_header("ChronoDecon Environment Diagnostic")
    print("This script checks if all required dependencies are installed.")

    check_python()

    print_header("Required Packages")
    required_packages = [
        ("streamlit", "streamlit"),
        ("pandas", "pandas"),
        ("numpy", "numpy"),
        ("plotly", "plotly"),
        ("matchms", "matchms"),
    ]

    missing_packages = []
    for package, import_name in required_packages:
        if not check_package(package, import_name):
            missing_packages.append(package)

    check_matchms_details()
    check_chronodecon()

    print_header("Summary")
    if missing_packages:
        print("❌ Missing packages detected!")
        print("\nPlease install them with:")
        print(f"   pip install {' '.join(missing_packages)}")
        print("\nOr install all dependencies:")
        print("   pip install -r requirements.txt")
        return 1
    else:
        print("✅ All dependencies are installed!")
        print("\nYou can now run the dashboard with:")
        print("   bash run_dashboard.sh")
        return 0

if __name__ == "__main__":
    sys.exit(main())
