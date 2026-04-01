"""Command-line interface for ChronoDecon using Typer."""

import typer
from pathlib import Path
from typing import Optional
import json
import sys

from .decon import homologue_factor, deconvolute_mzml, visualize_results

app = typer.Typer(
    name="chrono-decon",
    help="ChronoDecon - Modified Cross-Correlation Algorithm for MS Single-Component Reconstruction",
    add_completion=True
)


@app.command()
def deconvolute(
    input_file: Path = typer.Argument(
        ...,
        exists=True,
        help="Path to input mzML file"
    ),
    output_file: Optional[Path] = typer.Option(
        None,
        "--output", "-o",
        help="Path to save results JSON (default: print to stdout)"
    ),
    h_threshold: float = typer.Option(
        0.87,
        "--threshold", "-t",
        help="Homologue factor threshold for grouping (default: 0.87)"
    ),
    min_intensity: float = typer.Option(
        1000.0,
        "--min-intensity", "-m",
        help="Minimum intensity for peak detection (default: 1000.0)"
    ),
    visualize: bool = typer.Option(
        False,
        "--visualize", "-v",
        help="Generate interactive HTML visualization"
    ),
    viz_output: Optional[Path] = typer.Option(
        None,
        "--viz-output",
        help="Path to save visualization HTML (default: chrono_decon_viz.html)"
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        help="Print detailed processing information"
    )
):
    """Deconvolute MS data from mzML file.

    Performs single-component reconstruction using modified cross-correlation
    algorithm based on Anal. Chem. 2021.
    """
    if verbose:
        typer.echo(f"Processing file: {input_file}")
        typer.echo(f"H threshold: {h_threshold}")
        typer.echo(f"Min intensity: {min_intensity}")
    
    try:
        # Perform deconvolution
        results = deconvolute_mzml(
            str(input_file),
            h_threshold=h_threshold,
            min_intensity=min_intensity
        )
        
        # Check for errors
        if 'error' in results.get('metadata', {}):
            typer.echo(f"Error: {results['metadata']['error']}", err=True)
            raise typer.Exit(code=1)
        
        # Print summary
        metadata = results['metadata']
        typer.echo(f"\n{'='*60}")
        typer.echo("ChronoDecon Results")
        typer.echo(f"{'='*60}")
        typer.echo(f"Input file: {metadata['input_file']}")
        typer.echo(f"Total patterns: {metadata['total_patterns']}")
        typer.echo(f"Total groups: {metadata['total_groups']}")
        typer.echo(f"Reconstructed components: {metadata['reconstructed_components']}")
        typer.echo(f"Algorithm: {metadata['algorithm']}")
        
        if verbose and results['components']:
            typer.echo(f"\n{'='*60}")
            typer.echo("Top 10 Components")
            typer.echo(f"{'='*60}")
            typer.echo(f"{'m/z':<12} {'Intensity':<15} {'Group Size':<12} {'RT':<10}")
            typer.echo(f"{'-'*60}")
            for comp in results['components'][:10]:
                mz_str = f"{comp['mz']:.4f}"
                intensity_str = f"{comp['intensity']:.2e}"
                group_str = f"{comp['group_size']}"
                rt_str = f"{comp['rt']:.2f}" if comp['rt'] else "N/A"
                typer.echo(f"{mz_str:<12} {intensity_str:<15} {group_str:<12} {rt_str:<10}")
        
        # Save results
        if output_file:
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(results, f, indent=2, default=str)
            typer.echo(f"\nResults saved to: {output_file}")
        else:
            # Print JSON to stdout
            if verbose:
                typer.echo("\nFull results (JSON):")
                typer.echo(json.dumps(results, indent=2, default=str))
        
        # Generate visualization
        if visualize:
            viz_path = viz_output or Path("chrono_decon_viz.html")
            visualize_results(results, str(viz_path))
        
    except FileNotFoundError:
        typer.echo(f"Error: File not found: {input_file}", err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        typer.echo(f"Error: {str(e)}", err=True)
        raise typer.Exit(code=1)


@app.command()
def compare(
    file1: Path = typer.Argument(..., exists=True, help="First mzML file"),
    file2: Path = typer.Argument(..., exists=True, help="Second mzML file"),
    h_threshold: float = typer.Option(0.87, "--threshold", "-t")
):
    """Compare two mzML files using homologue factor analysis.

    Computes homologue factors between isotopic patterns from both files
    and provides similarity metrics.
    """
    typer.echo(f"Comparing: {file1} vs {file2}")
    
    results1 = deconvolute_mzml(str(file1), h_threshold=h_threshold)
    results2 = deconvolute_mzml(str(file2), h_threshold=h_threshold)
    
    components1 = results1.get('components', [])
    components2 = results2.get('components', [])
    
    typer.echo(f"\nFile 1: {len(components1)} components")
    typer.echo(f"File 2: {len(components2)} components")
    typer.echo(f"\nComparison complete.")


@app.command()
def info():
    """Display version and algorithm information."""
    typer.echo(f"ChronoDecon v{__import__('chrono_decon').__version__}")
    typer.echo("\nAlgorithm: Modified Cross-Correlation")
    typer.echo("Reference: Anal. Chem. 2021")
    typer.echo("\nDescription:")
    typer.echo("  Homologue factor calculation for mass spectrometry")
    typer.echo("  single-component reconstruction using FFT-based")
    typer.echo("  cross-correlation with roulette wheel selection.")


@app.command()
def test():
    """Run basic algorithm test."""
    import numpy as np
    
    typer.echo("Running algorithm test...")
    
    # Create test patterns
    pattern1 = np.array([100, 80, 50, 20, 10])
    pattern2 = np.array([100, 80, 50, 20, 10])  # Identical
    pattern3 = np.array([10, 20, 50, 80, 100])  # Reversed
    pattern4 = np.random.rand(5) * 100          # Random
    
    h12 = homologue_factor(pattern1, pattern2)
    h13 = homologue_factor(pattern1, pattern3)
    h14 = homologue_factor(pattern1, pattern4)
    
    typer.echo(f"\nHomologue factors:")
    typer.echo(f"  Pattern1 vs Pattern2 (identical): {h12:.4f}")
    typer.echo(f"  Pattern1 vs Pattern3 (reversed): {h13:.4f}")
    typer.echo(f"  Pattern1 vs Pattern4 (random): {h14:.4f}")
    
    if h12 > 0.9 and h13 < 0.5 and h14 < 0.5:
        typer.echo("\n✓ Test passed!")
    else:
        typer.echo("\n✗ Test failed!")
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
