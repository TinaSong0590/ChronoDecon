"""
MCP Tool for ChronoDecon - FastMCP integration for AI agents.

This module provides Model Context Protocol (MCP) tool interfaces
for seamless integration with Claude, LangChain, and other AI agents.
"""

import json
from typing import Optional
from pathlib import Path

from .decon import deconvolute_mzml, visualize_deconvolution


def deconvolute_mzml_tool(
    mzml_path: str,
    h_threshold: float = 0.87,
    min_intensity: float = 0.01,
    max_isotope_diff: int = 15,
    subtract_background: bool = True,
    bg_window_size: int = 20,
    bg_percentile: float = 10.0,
    symmetry_threshold: Optional[float] = None,
    auto_threshold: bool = True,
    simulate: bool = False,
    visualize: bool = False,
    viz_output: Optional[str] = None
) -> str:
    """
    Perform 2D deconvolution on an mzML file using the mXcorr algorithm (MCP protocol).

    This implements the full mXcorr workflow from the 2021 Analytical Chemistry paper:
    builds an all-ion chronogram (AIC) matrix, detects time shifts, and groups
    co-eluting peaks via the Dissymmetry Factor (DF).

    Parameters
    ----------
    mzml_path : str
        Path to the mzML file, or "simulate" to use simulated data
    h_threshold : float, default 0.87
        Legacy isotopic H-factor threshold (kept for backward compatibility)
    min_intensity : float, default 0.01
        Minimum relative intensity threshold (0-1)
    max_isotope_diff : int, default 15
        Maximum isotope position difference (legacy)
    subtract_background : bool, default True
        Whether to perform background subtraction on chronograms
    bg_window_size : int, default 20
        Window size for background estimation
    bg_percentile : float, default 10.0
        Percentile for background threshold
    symmetry_threshold : float, optional
        DF symmetry threshold for grouping. None = auto-determine.
    auto_threshold : bool, default True
        Whether to auto-determine the symmetry threshold
    simulate : bool, default False
        Use simulated data for testing (ignores file path)
    visualize : bool, default False
        Whether to generate visualization
    viz_output : str, optional
        Path to save visualization HTML

    Returns
    -------
    str
        JSON formatted deconvolution results
    """
    # Simulation mode for quick testing
    if simulate or mzml_path == "simulate":
        simulated_result = {
            'status': 'success',
            'mode': 'simulation',
            'message': 'Simulated deconvolution result for testing',
            'parameters': {
                'h_threshold': h_threshold,
                'min_intensity': min_intensity,
                'max_isotope_diff': max_isotope_diff
            },
            'components': [
                {
                    'component_id': 1,
                    'mz_values': [100.0, 101.0, 102.0, 103.0],
                    'intensities': [1.0, 0.65, 0.15, 0.02],
                    'average_mz': 100.0,
                    'total_intensity': 1.0
                },
                {
                    'component_id': 2,
                    'mz_values': [200.0, 201.0, 202.0],
                    'intensities': [1.0, 0.70, 0.12],
                    'average_mz': 200.0,
                    'total_intensity': 0.85
                },
                {
                    'component_id': 3,
                    'mz_values': [350.0, 351.0, 352.0, 353.0],
                    'intensities': [1.0, 0.68, 0.18, 0.03],
                    'average_mz': 350.0,
                    'total_intensity': 0.92
                }
            ],
            'num_components': 3,
            'num_peaks': 11
        }
        return json.dumps(simulated_result, indent=2, ensure_ascii=False)
    
    # Validate input file exists
    mzml_file = Path(mzml_path)
    if not mzml_file.exists():
        return json.dumps({
            'status': 'error',
            'message': f'File not found: {mzml_path}'
        }, indent=2, ensure_ascii=False)
    
    # Perform deconvolution (2D mXcorr approach)
    result = deconvolute_mzml(
        mzml_path=mzml_path,
        h_threshold=h_threshold,
        min_intensity=min_intensity,
        max_isotope_diff=max_isotope_diff,
        subtract_background=subtract_background,
        bg_window_size=bg_window_size,
        bg_percentile=bg_percentile,
        symmetry_threshold=symmetry_threshold,
        auto_threshold=auto_threshold
    )
    
    # Generate visualization if requested
    if visualize and result['status'] == 'success':
        try:
            if viz_output:
                visualize_deconvolution(result, output_path=viz_output)
                result['visualization'] = f'Generated: {viz_output}'
            else:
                visualize_deconvolution(result)
                result['visualization'] = 'Displayed in browser'
        except Exception as e:
            result['visualization'] = f'Error: {str(e)}'
    
    # Convert numpy types to Python types for JSON serialization
    def convert_to_serializable(obj):
        if isinstance(obj, dict):
            return {k: convert_to_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_to_serializable(v) for v in obj]
        elif hasattr(obj, 'dtype'):  # numpy types
            return obj.item() if obj.ndim == 0 else obj.tolist()
        else:
            return obj
    
    serializable_result = convert_to_serializable(result)
    return json.dumps(serializable_result, indent=2, ensure_ascii=False)


# FastMCP integration
try:
    from fastmcp import FastMCP
    
    # Create FastMCP instance
    mcp = FastMCP("chrono-decon")
    
    @mcp.tool()
    def deconvolute(
        mzml_path: str,
        symmetry_threshold: float = 0.0,
        min_intensity: float = 0.01,
        auto_threshold: bool = True
    ) -> str:
        """
        Deconvolute mass spectrometry mzML file using the mXcorr algorithm.

        Performs 2D deconvolution (time x m/z) by building all-ion chronograms
        and grouping co-eluting peaks via the Dissymmetry Factor.

        Args:
            mzml_path: Path to input mzML file
            symmetry_threshold: DF threshold for grouping (0=auto, default: 0)
            min_intensity: Minimum relative intensity threshold (default: 0.01)
            auto_threshold: Auto-determine symmetry threshold (default: True)

        Returns:
            JSON string with deconvolution results including:
            - Number of identified components
            - m/z values and intensities for each component
            - Dissymmetry Factor threshold used
            - Time window bounds
        """
        st = None if symmetry_threshold <= 0 else symmetry_threshold
        return deconvolute_mzml_tool(
            mzml_path=mzml_path,
            symmetry_threshold=st,
            min_intensity=min_intensity,
            auto_threshold=auto_threshold
        )
    
    HAS_FASTMCP = True

    @mcp.tool()
    def library_search(
        query_path: str,
        library_path: str = "",
        tolerance: float = 0.1,
        min_match: int = 3,
        top_n: int = 15,
        proxy: str = ""
    ) -> str:
        """
        Search GNPS spectral library for compound identification.

        Matches deconvoluted spectra against the GNPS public spectral library
        using CosineGreedy similarity scoring.

        Args:
            query_path: Path to query spectra (.mgf or .msp)
            library_path: Path to GNPS library (.mgf). Empty = auto-download.
            tolerance: m/z tolerance in Da (default: 0.1)
            min_match: Minimum matched peaks (default: 3)
            top_n: Top N matches per query (default: 15)
            proxy: HTTP proxy URL for downloading library

        Returns:
            JSON string with matching results including compound names,
            cosine scores, matched peaks, InChIKeys, and SMILES
        """
        from .library_search import msp_to_mgf, run_local_search

        query_file = Path(query_path)
        if not query_file.exists():
            return json.dumps({"status": "error", "message": f"File not found: {query_path}"}, indent=2)

        # Convert MSP to MGF if needed
        mgf_path = query_path
        if query_file.suffix.lower() == ".msp":
            mgf_path = str(query_file).replace(".msp", ".mgf")
            try:
                msp_to_mgf(query_path, mgf_path)
            except Exception as e:
                return json.dumps({"status": "error", "message": f"MSP conversion failed: {e}"}, indent=2)

        try:
            results = run_local_search(
                query_mgf=mgf_path,
                library_path=library_path if library_path else None,
                tolerance=tolerance,
                min_match=min_match,
                top_n=top_n,
                proxy=proxy if proxy else None,
            )
            results["status"] = "success"
            return json.dumps(results, indent=2, ensure_ascii=False)
        except Exception as e:
            return json.dumps({"status": "error", "message": str(e)}, indent=2, ensure_ascii=False)

except ImportError:
    HAS_FASTMCP = False


if __name__ == "__main__":
    print("🚀 Starting ChronoDecon MCP Server on http://0.0.0.0:8000")
    print("Tool registered: deconvolute_mzml_tool")
    print("Press Ctrl+C to stop the server...")
    mcp.run()
