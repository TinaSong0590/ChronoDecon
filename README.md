# ChronoDecon

Chronological Deconvolution for Mass Spectrometry - A Python package implementing the modified correlation algorithm from the 2021 Analytical Chemistry paper for single component reconstruction.

## Features

- **FFT-based Homologue Factor Calculation**: Uses cross-correlation to identify isotopic patterns
- **Adaptive Threshold Determination**: Implements D-response function for optimal HThr calculation
- **GNPS Library Integration**: Support for spectral matching against GNPS database
- **Elemental Analysis**: Formula validation with PPM error calculation
- **MCP Protocol Support**: Integration with AI agents and LangChain
- **Rich Visualization**: Interactive HTML reports using Plotly

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/chrono-decon.git
cd chrono-decon

# Install in development mode
pip install -e .
```

Or install with development dependencies:
```bash
pip install -e ".[dev]"
```

## Quick Start

### Command Line Interface

```bash
# Basic deconvolution
chrono-decon deconvolute input.mzML -o results.json

# With custom parameters
chrono-decon deconvolute input.mzML -h 0.9 -m 0.02 -o results.json

# With visualization
chrono-decon deconvolute input.mzML --visualize --viz-output plot.html
```

### Python API

```python
from chrono_decon import deconvolute_mzml

# Perform deconvolution
result = deconvolute_mzml(
    "path/to/file.mzML",
    h_threshold=0.87,
    min_intensity=0.01
)

# Access results
print(f"Found {result['num_components']} components")
```

## Core Algorithm

### Homologue Factor

The `homologue_factor()` function calculates the similarity between two isotopic distributions using FFT-based cross-correlation:

```python
from chrono_decon import homologue_factor
import numpy as np

psi_a = np.array([100, 50, 10, 2, 1])
psi_b = np.array([95, 48, 9, 3, 1])

h_factor = homologue_factor(psi_a, psi_b)
print(f"Homologue factor: {h_factor:.3f}")
```

### Deconvolution Workflow

1. Parse mzML file and extract spectrum data
2. Detect peaks in the spectrum
3. Calculate homologue factors between peaks
4. Group peaks into components using roulette wheel selection
5. Return deconvolution results

## Parameters

- `h_threshold`: Homologue factor threshold for grouping (default: 0.87)
- `min_intensity`: Minimum relative intensity threshold (0-1, default: 0.01)
- `max_isotope_diff`: Maximum isotope position difference to consider (default: 15)

## MCP Server Integration

Start the MCP server for AI agent integration:

```bash
# Start MCP server
python -m chrono_decon.mcp_tool
```

## Citation

If you use ChronoDecon in your research, please cite the original paper:

> "Single Component Reconstruction from Mass Spectrometry Data Using Modified Correlation Analysis", Analytical Chemistry, 2021.

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Acknowledgments

Based on the algorithm described in:
"Single Component Reconstruction from Mass Spectrometry Data Using Modified Correlation Analysis", Analytical Chemistry, 2021.
