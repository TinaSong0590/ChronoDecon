#!/bin/bash
# ChronoDecon Dashboard Launcher
# Works with any Python environment that has the required dependencies installed.

echo "================================"
echo "  ChronoDecon Dashboard v0.2.1"
echo "================================"
echo ""

# Switch to repo root (where this script lives)
REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_ROOT"

# Verify we're in the right place
if [ ! -f "app.py" ] || [ ! -d "chrono_decon" ]; then
    echo "Error: Please run this script from the repository root."
    echo "Expected: app.py and chrono_decon/ in the same directory."
    echo "Current directory: $(pwd)"
    exit 1
fi

# Use the current Python (whatever environment the user has activated)
PYTHON_CMD="${PYTHON_CMD:-python3}"
STREAMLIT_CMD="${STREAMLIT_CMD:-streamlit}"

echo "Python: $($PYTHON_CMD --version 2>&1)"
echo "Streamlit: $($STREAMLIT_CMD version 2>&1 | head -1)"

# Verify core dependencies
echo ""
echo "Checking dependencies..."

$PYTHON_CMD -c "import matchms" 2>/dev/null && echo "  matchms: OK" || {
    echo "  matchms: NOT FOUND - installing..."
    $PYTHON_CMD -m pip install matchms || {
        echo "  Error: Cannot install matchms. Please run: pip install -e ."
        exit 1
    }
}

$PYTHON_CMD -c "import streamlit" 2>/dev/null && echo "  streamlit: OK" || {
    echo "  streamlit: NOT FOUND - installing..."
    $PYTHON_CMD -m pip install streamlit || {
        echo "  Error: Cannot install streamlit."
        exit 1
    }
}

echo ""
echo "Starting ChronoDecon Dashboard..."
echo "Access: http://localhost:8501"
echo ""
echo "Press Ctrl+C to stop."
echo ""

cd "$REPO_ROOT"
$STREAMLIT_CMD run app.py
