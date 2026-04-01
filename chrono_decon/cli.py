"""
Command-line interface for ChronoDecon using Typer.

This module provides a CLI for performing deconvolution on mass spectrometry data.
"""

import typer
from pathlib import Path
from typing import Optional
import json
import sys
import logging

from .decon import deconvolute_mzml, visualize_deconvolution

app = typer.Typer(
    name="chrono-decon",
    help="Chronological Deconvolution for Mass Spectrometry",
    add_completion=False
)

logger = logging.getLogger("chrono_decon")


@app.command()
def deconvolute(
    mzml_path: str = typer.Argument(
        ...,
        help="Path to the input mzML file",
        exists=True
    ),
    output: Optional[str] = typer.Option(
        None,
        "-o",
        "--output",
        help="Output JSON file path for results"
    ),
    h_threshold: float = typer.Option(
        0.87,
        "-h",
        "--h-threshold",
        help="Homologue factor threshold for grouping (0-1, default: 0.87)"
    ),
    min_intensity: float = typer.Option(
        0.01,
        "-m",
        "--min-intensity",
        help="Minimum relative intensity threshold (0-1, default: 0.01)"
    ),
    max_isotope_diff: int = typer.Option(
        15,
        "-d",
        "--max-isotope-diff",
        help="Maximum isotope position difference to consider (default: 15)"
    ),
    visualize: bool = typer.Option(
        False,
        "-v",
        "--visualize",
        help="Generate interactive visualization"
    ),
    viz_output: Optional[str] = typer.Option(
        None,
        "--viz-output",
        help="Output path for visualization HTML file"
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-V",
        help="Show detailed output"
    )
) -> None:
    """
    Perform single component deconvolution on an mzML file.
    
    This command implements the modified correlation algorithm for identifying
    and deconvoluting single components from mass spectrometry data.
    
    Example usage:
        chrono-decon deconvolute input.mzML -o results.json --visualize
        chrono-decon deconvolute input.mzML -h 0.9 -m 0.02
    """
    if verbose:
        typer.echo(f"Processing file: {mzml_path}")
        typer.echo(f"H-threshold: {h_threshold}")
        typer.echo(f"Min intensity: {min_intensity}")
        typer.echo(f"Max isotope diff: {max_isotope_diff}")
    
    # Perform deconvolution
    result = deconvolute_mzml(
        mzml_path=mzml_path,
        h_threshold=h_threshold,
        min_intensity=min_intensity,
        max_isotope_diff=max_isotope_diff
    )
    
    # Check for errors
    if result['status'].startswith('error'):
        typer.echo(f"Error: {result['status']}", err=True)
        sys.exit(1)
    
    # Print summary
    typer.echo(f"✓ Deconvolution completed successfully")
    typer.echo(f"  File: {result['file_path']}")
    typer.echo(f"  Components identified: {result['num_components']}")
    
    if verbose:
        typer.echo(f"\nComponents:")
        for i, comp in enumerate(result['components']):
            typer.echo(f"  Component {i+1}:")
            typer.echo(f"    Peaks: {len(comp['peak_indices'])}")
            if comp['h_factors']:
                typer.echo(f"    H-factors: {[f'{h:.3f}' for h in comp['h_factors']]}")
    
    # Save results to JSON
    if output:
        output_path = Path(output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(result, f, indent=2)
        
        typer.echo(f"✓ Results saved to: {output}")
    
    # Generate visualization
    if visualize:
        try:
            if viz_output:
                visualize_deconvolution(result, output_path=viz_output)
                typer.echo(f"✓ Visualization saved to: {viz_output}")
            else:
                visualize_deconvolution(result)
        except Exception as e:
            typer.echo(f"Warning: Could not generate visualization: {e}", err=True)


@app.command("library-search")
def library_search(
    query: str = typer.Argument(
        ...,
        help="Path to query spectra file (.mgf or .msp). If .msp, will auto-convert to MGF."
    ),
    output_dir: Optional[str] = typer.Option(
        None,
        "-o",
        "--output-dir",
        help="Output directory for reports (default: same as query file)"
    ),
    library_path: Optional[str] = typer.Option(
        None,
        "-l",
        "--library",
        help="Path to GNPS library MGF file (default: auto-download to ~/.chrono_decon/)"
    ),
    tolerance: float = typer.Option(
        0.1,
        "-t",
        "--tolerance",
        help="m/z tolerance in Da for peak matching (default: 0.1)"
    ),
    min_match: int = typer.Option(
        3,
        "--min-match",
        help="Minimum number of matched peaks (default: 3)"
    ),
    top_n: int = typer.Option(
        15,
        "-n",
        "--top-n",
        help="Number of top matches to return per query (default: 15)"
    ),
    proxy: Optional[str] = typer.Option(
        None,
        "--proxy",
        help="HTTP proxy for downloading GNPS library, e.g. http://127.0.0.1:7890"
    ),
    cloud: bool = typer.Option(
        False,
        "--cloud",
        help="Upload to GNPS cloud for matching (requires internet)"
    ),
    force_download: bool = typer.Option(
        False,
        "--force-download",
        help="Force re-download of GNPS library even if cached"
    ),
) -> None:
    """
    Search GNPS spectral library for compound identification.

    Supports dual-mode: local matchms comparison (default) or GNPS cloud upload (--cloud).

    Examples:
        chrono-decon library-search query.mgf -o results/
        chrono-decon library-search query.msp -l /path/to/gnps_library.mgf --proxy http://127.0.0.1:7890
        chrono-decon library-search query.mgf --cloud
    """
    from .library_search import (
        msp_to_mgf, run_local_search, generate_html_report,
        generate_csv_report, print_summary,
    )

    query_path = Path(query)
    if not query_path.exists():
        typer.echo(f"Error: Query file not found: {query}", err=True)
        raise typer.Exit(1)

    # Determine output directory
    if output_dir:
        out_dir = Path(output_dir)
    else:
        out_dir = query_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = query_path.stem

    # Convert MSP to MGF if needed
    if query_path.suffix.lower() == ".msp":
        mgf_path = out_dir / f"{stem}.mgf"
        typer.echo(f"Converting MSP -> MGF: {query} -> {mgf_path}")
        try:
            msp_to_mgf(str(query_path), str(mgf_path))
        except Exception as e:
            typer.echo(f"Error converting MSP: {e}", err=True)
            raise typer.Exit(1)
    else:
        mgf_path = query_path

    # Cloud mode
    if cloud:
        typer.echo("Cloud mode: uploading to GNPS...")
        typer.echo(f"  Query file: {mgf_path}")
        typer.echo("  Please submit your file at: https://gnps.ucsd.edu/ProteoSAFe/index.jsp?params=#workflow=MASSBANK")
        typer.echo("  Or use the GNPS REST API for programmatic access.")
        return

    # Local mode
    typer.echo("=" * 60)
    typer.echo("  ChronoDecon - GNPS Library Search (Local Mode)")
    typer.echo("=" * 60)
    typer.echo(f"  Query:    {mgf_path}")
    typer.echo(f"  Library:  {library_path or 'auto-download'}")
    typer.echo(f"  Tolerance: {tolerance} Da")
    typer.echo(f"  Min peaks: {min_match}")
    typer.echo("=" * 60)

    try:
        results = run_local_search(
            query_mgf=str(mgf_path),
            library_path=library_path,
            tolerance=tolerance,
            min_match=min_match,
            top_n=top_n,
            proxy=proxy,
        )
    except Exception as e:
        typer.echo(f"Error running library search: {e}", err=True)
        if "download" in str(e).lower():
            typer.echo(f"Hint: Try downloading manually with proxy:")
            typer.echo(f"  wget -e use_proxy=yes -e http_proxy={proxy} -O ~/.chrono_decon/gnps_library.mgf \\")
            typer.echo(f"    https://gnps.ucsd.edu/ProteoSAFe/static/gnps_library.mgf", err=True)
        raise typer.Exit(1)

    # Save reports
    html_path = str(out_dir / f"{stem}_gnps_report.html")
    csv_path = str(out_dir / f"{stem}_gnps_results.csv")
    json_path = str(out_dir / f"{stem}_gnps_results.json")

    generate_html_report(results, html_path)
    generate_csv_report(results, csv_path)
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    # Print summary
    print_summary(results, top_n=10)

    typer.echo(f"\n  Reports saved:")
    typer.echo(f"    HTML: {html_path}")
    typer.echo(f"    CSV:  {csv_path}")
    typer.echo(f"    JSON: {json_path}")


@app.command()
def version() -> None:
    """Show version information."""
    from . import __version__
    typer.echo(f"ChronoDecon version: {__version__}")


@app.command()
def mcp():
    """Start MCP (Model Context Protocol) server for AI agent integration."""
    try:
        from .mcp_tool import mcp
        mcp.run()
    except ImportError as e:
        typer.echo(f"Error: {e}", err=True)
        typer.echo("Install dependencies with: pip install -e .[dev]", err=True)
        raise typer.Exit(1)


def main():
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()
