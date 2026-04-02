"""
RAW to mzML Conversion Module for ChronoDecon

Supports multiple conversion methods with automatic fallback:
1. PyOpenMS (pip-installable, primary method)
2. ThermoRawFileParser (external binary, auto-downloaded)
3. msconvert via Docker (optional)

Users can set proxy via environment variables:
  CHRONODECON_PROXY or http_proxy/https_proxy
"""

import subprocess
import sys
import os
import shutil
import tempfile
import zipfile
from pathlib import Path
from urllib.request import urlretrieve, urlopen, Request
from urllib.error import URLError


# ThermoRawFileParser download info
_THERMO_VERSION = "v2.0.0-dev"
_THERMO_LINUX_URL = (
    "https://github.com/CompOmics/ThermoRawFileParser/releases/download/"
    f"{_THERMO_VERSION}/ThermoRawFileParser-{_THERMO_VERSION}-linux.zip"
)
# Store the parser binary inside the package data directory
_THERMO_DIR = Path(__file__).parent / "bin" / "ThermoRawFileParser"


def _get_proxy_handler():
    """Get proxy settings from environment variables."""
    proxy_url = (
        os.environ.get("CHRONODECON_PROXY")
        or os.environ.get("https_proxy")
        or os.environ.get("HTTPS_PROXY")
        or os.environ.get("http_proxy")
        or os.environ.get("HTTP_PROXY")
        or None
    )
    return proxy_url


def _download_file(url, dest_path, desc=None):
    """Download a file with optional proxy support and progress."""
    proxy_url = _get_proxy_handler()
    if proxy_url:
        # Simple proxy-aware download using subprocess + curl (if available)
        if shutil.which("curl"):
            cmd = [
                "curl", "-L", "-x", proxy_url,
                "-o", str(dest_path), url
            ]
            label = desc or url.split("/")[-1]
            print(f"  Downloading {label} via proxy...")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if result.returncode != 0:
                raise RuntimeError(f"curl download failed: {result.stderr}")
            return
        # Fallback: set env and use urllib
        os.environ["http_proxy"] = proxy_url
        os.environ["https_proxy"] = proxy_url

    label = desc or url.split("/")[-1]
    print(f"  Downloading {label}...")

    def _report(block_num, block_size, total_size):
        downloaded = block_num * block_size
        pct = downloaded / total_size * 100 if total_size > 0 else 0
        sys.stdout.write(f"\r  Progress: {pct:.0f}%")
        sys.stdout.flush()

    urlretrieve(url, str(dest_path), reporthook=_report)
    print()  # newline after progress


def _ensure_thermo_parser():
    """Download and extract ThermoRawFileParser if not present.
    Returns path to the executable, or None if unavailable.
    """
    # Check if already installed
    exe_path = _THERMO_DIR / "ThermoRawFileParser"
    if exe_path.exists() and os.access(exe_path, os.X_OK):
        return str(exe_path)

    # Check system PATH
    if shutil.which("ThermoRawFileParser"):
        return "ThermoRawFileParser"

    print("  ThermoRawFileParser not found locally. Downloading...")

    try:
        _THERMO_DIR.mkdir(parents=True, exist_ok=True)

        # Download
        zip_path = Path(tempfile.gettempdir()) / "ThermoRawFileParser.zip"
        _download_file(_THERMO_LINUX_URL, zip_path, "ThermoRawFileParser")

        # Extract
        print("  Extracting...")
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(str(_THERMO_DIR))

        # Find the executable (it might be nested)
        for root, dirs, files in os.walk(str(_THERMO_DIR)):
            for f in files:
                if f == "ThermoRawFileParser":
                    found = Path(root) / f
                    found.chmod(0o755)
                    return str(found)

        print("  Warning: ThermoRawFileParser executable not found after extraction")
        return None

    except Exception as e:
        print(f"  Warning: Failed to download ThermoRawFileParser: {e}")
        return None


def convert_pyopenms(raw_path, output_dir):
    """Convert RAW to mzML using PyOpenMS (pip-installable).

    Requires: pip install pyopenms
    """
    try:
        from pyopenms import MzMLFile
    except ImportError:
        raise ImportError(
            "PyOpenMS is not installed. Install it with:\n"
            "  pip install pyopenms\n"
            "Or let ChronoDecon auto-download ThermoRawFileParser."
        )

    raw_file = Path(raw_path)
    output_path = Path(output_dir)
    mzml_file = output_path / f"{raw_file.stem}.mzML"

    try:
        fh = __import__("pyopenms", fromlist=["FileHandler"]).FileHandler()
        fh.loadExperiment(str(raw_file))
        exp = fh.getExperiment()
        MzMLFile().store(str(mzml_file), exp)

        if mzml_file.exists():
            return str(mzml_file)
        raise RuntimeError("mzML file was not created by PyOpenMS")
    except Exception as e:
        raise RuntimeError(f"PyOpenMS conversion failed: {e}")


def convert_thermo_parser(raw_path, output_dir):
    """Convert RAW to mzML using ThermoRawFileParser (auto-downloaded binary)."""
    exe = _ensure_thermo_parser()
    if not exe:
        raise RuntimeError("ThermoRawFileParser is not available")

    raw_file = Path(raw_path)
    output_path = Path(output_dir)
    mzml_file = output_path / f"{raw_file.stem}.mzML"

    cmd = [exe, "-i", str(raw_file), "-o", str(output_path), "-f", "2"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0 and not mzml_file.exists():
            raise RuntimeError(f"ThermoRawFileParser failed: {result.stderr}")

        if mzml_file.exists():
            return str(mzml_file)
        # Search for any mzML output
        for f in output_path.glob("*.mzML"):
            return str(f)
        raise RuntimeError("mzML file not found after conversion")
    except subprocess.TimeoutExpired:
        raise RuntimeError("Conversion timed out (10 min)")


def convert_docker(raw_path, output_dir):
    """Convert RAW to mzML using Docker + ProteoWizard msconvert."""
    if not shutil.which("docker"):
        raise RuntimeError("Docker is not available on this system")

    raw_file = Path(raw_path)
    mzml_file = Path(output_dir) / f"{raw_file.stem}.mzML"

    cmd = [
        "docker", "run", "--rm",
        "-v", f"{raw_file.parent}:/data",
        "proteme/proteowizard:latest",
        "msconvert", f"/data/{raw_file.name}", "--mzML", "-o", "/data/"
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            raise RuntimeError(f"Docker conversion failed: {result.stderr}")

        # Output may be in the original directory
        for candidate in [mzml_file, raw_file.parent / f"{raw_file.stem}.mzML"]:
            if candidate.exists():
                return str(candidate)
        raise RuntimeError("mzML file not found after Docker conversion")
    except subprocess.TimeoutExpired:
        raise RuntimeError("Conversion timed out (10 min)")


def convert_raw_to_mzml(raw_path, output_dir):
    """Convert RAW file to mzML with automatic method selection.

    Tries methods in order:
    1. ThermoRawFileParser (auto-downloaded, most reliable)
    2. PyOpenMS (if already installed)
    3. Docker ProteoWizard (if Docker available)

    Args:
        raw_path: Path to the .raw file
        output_dir: Directory for the output .mzML file

    Returns:
        str: Path to the converted .mzML file

    Raises:
        RuntimeError: If all conversion methods fail
    """
    errors = []

    # Method 1: ThermoRawFileParser (auto-download)
    try:
        print("Method 1: Trying ThermoRawFileParser (auto-download)...")
        return convert_thermo_parser(raw_path, output_dir)
    except Exception as e:
        errors.append(f"ThermoRawFileParser: {e}")
        print(f"  Failed: {e}")

    # Method 2: PyOpenMS
    try:
        print("Method 2: Trying PyOpenMS...")
        return convert_pyopenms(raw_path, output_dir)
    except Exception as e:
        errors.append(f"PyOpenMS: {e}")
        print(f"  Failed: {e}")

    # Method 3: Docker
    try:
        print("Method 3: Trying Docker ProteoWizard...")
        return convert_docker(raw_path, output_dir)
    except Exception as e:
        errors.append(f"Docker: {e}")
        print(f"  Failed: {e}")

    raise RuntimeError(
        "All RAW conversion methods failed:\n"
        + "\n".join(f"  - {e}" for e in errors)
        + "\n\nPlease install one of:\n"
        "  1. pip install pyopenms  (Python package)\n"
        "  2. Ensure internet access for auto-download of ThermoRawFileParser\n"
        "  3. Install Docker for ProteoWizard support\n"
        "  4. Or convert your .raw file to .mzML beforehand using Thermo Xcalibur / ProteoWizard"
    )
