#!/bin/bash
# ChronoDecon Environment Fix Script
# Diagnoses and fixes common dependency issues.

echo "========================================"
echo "  ChronoDecon Environment Fix"
echo "========================================"
echo ""

PYTHON_CMD="${PYTHON_CMD:-python3}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Step 1: Stop stale Streamlit processes
echo "[1/3] Stopping stale Streamlit processes..."
pkill -f "streamlit run" 2>/dev/null
echo "  Done."

# Step 2: Clear Streamlit cache
echo "[2/3] Clearing Streamlit cache..."
rm -rf ~/.streamlit/cache/ 2>/dev/null
echo "  Done."

# Step 3: Verify/install dependencies
echo "[3/3] Checking dependencies..."
echo ""

check_and_install() {
    local pkg=$1
    if $PYTHON_CMD -c "import $pkg" 2>/dev/null; then
        echo "  $pkg: OK"
    else
        echo "  $pkg: NOT FOUND - installing..."
        $PYTHON_CMD -m pip install "$pkg"
    fi
}

check_and_install "streamlit"
check_and_install "matchms"
check_and_install "pymzml"
check_and_install "plotly"
check_and_install "numpy"
check_and_install "scipy"
check_and_install "pandas"

echo ""
echo "========================================"
echo "  Fix Complete!"
echo "========================================"
echo ""
echo "Starting Dashboard..."
echo "Press Ctrl+C to stop."
echo ""

cd "$SCRIPT_DIR"
$PYTHON_CMD -m streamlit run app.py
