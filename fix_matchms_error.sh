#!/bin/bash
# Quick fix script for "No module named 'matchms'" error

echo "========================================"
echo "  ChronoDecon Matchms Error Fix"
echo "========================================"
echo ""

# Use miniconda3 environment (has matchms)
PYTHON_ENV="/home/knan/miniconda3/bin/python3"
PIP_CMD="/home/knan/miniconda3/bin/pip"
STREAMLIT_CMD="/home/knan/miniconda3/bin/streamlit"

# Step 1: Stop all Streamlit processes
echo "[1/4] Stopping Streamlit processes..."
pkill -f "streamlit run" 2>/dev/null
echo "✓ Streamlit processes stopped"
echo ""

# Step 2: Clear Streamlit cache
echo "[2/4] Clearing Streamlit cache..."
rm -rf ~/.streamlit/cache/ 2>/dev/null
echo "✓ Streamlit cache cleared"
echo ""

# Step 3: Verify matchms installation
echo "[3/4] Verifying matchms installation in miniconda3..."
if $PYTHON_ENV -c "import matchms; print(matchms.__version__)" 2>/dev/null; then
    echo "✓ matchms is installed"
else
    echo "✗ matchms is NOT installed"
    echo ""
    echo "Installing matchms..."
    $PIP_CMD install matchms
    echo ""
fi
echo ""

# Step 4: Run environment diagnostic
echo "[4/4] Running environment diagnostic..."
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
$PYTHON_ENV "$SCRIPT_DIR/diagnose_environment.py"

echo ""
echo "========================================"
echo "  Fix Complete!"
echo "========================================"
echo ""
echo "Starting Dashboard with miniconda3 environment..."
echo "Press Ctrl+C to stop"
echo ""

# Start the dashboard
cd "$SCRIPT_DIR"
$STREAMLIT_CMD run app.py
