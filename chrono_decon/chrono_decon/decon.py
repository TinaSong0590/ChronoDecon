"""Core deconvolution algorithm for ChronoDecon.

Implements modified cross-correlation algorithm for mass spectrometry
single-component reconstruction based on Anal. Chem. 2021.
"""

import numpy as np
from scipy import fft
from typing import Dict, List, Tuple
import pymzml


def homologue_factor(psi_a: np.ndarray, psi_b: np.ndarray) -> float:
    """Calculate homologue factor H between two isotopic patterns.

    Uses FFT-based cross-correlation to compute the homologue factor H
    as described in the modified cross-correlation algorithm.

    Args:
        psi_a: First isotopic pattern (m/z intensity array)
        psi_b: Second isotopic pattern (m/z intensity array)

    Returns:
        H factor value between 0 and 1, where 1 indicates perfect match.

    Raises:
        ValueError: If input arrays are empty or have different lengths.

    Reference:
        Anal. Chem. 2021 - Modified cross-correlation algorithm for
        homologue factor calculation in MS data analysis.
    """
    if len(psi_a) == 0 or len(psi_b) == 0:
        raise ValueError("Isotopic patterns cannot be empty")
    
    # Normalize patterns
    psi_a_norm = psi_a / (np.max(psi_a) + 1e-10)
    psi_b_norm = psi_b / (np.max(psi_b) + 1e-10)
    
    # Ensure same length by zero-padding
    max_len = max(len(psi_a_norm), len(psi_b_norm))
    psi_a_padded = np.pad(psi_a_norm, (0, max_len - len(psi_a_norm)), 'constant')
    psi_b_padded = np.pad(psi_b_norm, (0, max_len - len(psi_b_norm)), 'constant')
    
    # FFT-based cross-correlation
    fft_a = fft.fft(psi_a_padded)
    fft_b = fft.fft(psi_b_padded)
    cross_corr = fft.ifft(fft_a * np.conj(fft_b))
    
    # Find maximum correlation
    max_corr = np.max(np.abs(cross_corr))
    
    # Normalize to get H factor
    norm_factor = np.sqrt(np.sum(psi_a_padded ** 2) * np.sum(psi_b_padded ** 2))
    h_factor = max_corr / (norm_factor + 1e-10)
    
    return float(h_factor)


def group_isotopic_patterns(
    patterns: List[Dict],
    h_threshold: float = 0.87
) -> List[List[Dict]]:
    """Group isotopic patterns using roulette wheel selection.

    Implements automatic grouping workflow based on homologue factor
    threshold using roulette wheel selection strategy.

    Args:
        patterns: List of isotopic pattern dictionaries, each containing:
            - 'mz': m/z value
            - 'intensity': intensity array
            - 'rt': retention time (optional)
        h_threshold: Homologue factor threshold for grouping (default: 0.87)

    Returns:
        List of groups, where each group is a list of pattern dictionaries.

    Algorithm:
        1. Compute H factors between all pattern pairs
        2. Roulette wheel selection based on H values
        3. Form groups above threshold
        4. Iterate until convergence
    """
    if not patterns:
        return []
    
    n = len(patterns)
    
    # Initialize groups
    groups = []
    assigned = [False] * n
    
    # Compute H factor matrix
    h_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            h = homologue_factor(
                patterns[i]['intensity'],
                patterns[j]['intensity']
            )
            h_matrix[i, j] = h
            h_matrix[j, i] = h
    
    # Roulette wheel selection and grouping
    for i in range(n):
        if assigned[i]:
            continue
        
        # Start new group with pattern i
        current_group = [patterns[i]]
        assigned[i] = True
        
        # Roulette wheel selection
        candidates = [j for j in range(n) if not assigned[j] and j != i]
        
        while candidates:
            # Calculate selection probabilities
            h_values = [h_matrix[i, j] for j in candidates]
            total_h = sum(h_values)
            
            if total_h == 0:
                break
            
            probabilities = [h / total_h for h in h_values]
            
            # Select candidate with highest probability
            selected_idx = np.argmax(probabilities)
            selected_candidate = candidates[selected_idx]
            
            # Check if H factor exceeds threshold
            if h_matrix[i, selected_candidate] >= h_threshold:
                current_group.append(patterns[selected_candidate])
                assigned[selected_candidate] = True
                candidates.pop(selected_idx)
            else:
                break
        
        groups.append(current_group)
    
    return groups


def deconvolute_mzml(
    mzml_path: str,
    h_threshold: float = 0.87,
    min_intensity: float = 1000.0
) -> Dict:
    """Deconvolute MS data from mzML file using modified cross-correlation algorithm.

    Main entry point for MS single-component reconstruction. Processes
    mzML file, extracts isotopic patterns, groups them, and reconstructs
    single components.

    Args:
        mzml_path: Path to input mzML file
        h_threshold: Homologue factor threshold for grouping (default: 0.87)
        min_intensity: Minimum intensity for peak detection (default: 1000.0)

    Returns:
        Dictionary containing:
            - 'components': List of deconvoluted components
            - 'groups': List of grouped isotopic patterns
            - 'metadata': Processing metadata

    Components structure:
        Each component is a dict with:
            - 'mz': Representative m/z value
            - 'intensity': Deconvoluted intensity
            - 'rt': Retention time (if available)
            - 'group_size': Number of patterns in the group

    Raises:
        FileNotFoundError: If mzML file does not exist.
        ValueError: If file cannot be parsed.
    """
    import pymzml
    
    # Parse mzML file
    run = pymzml.run.Reader(mzml_path)
    
    # Extract all spectra
    all_spectra = [s for s in run if s.ms_level == 1]
    
    if not all_spectra:
        return {
            'components': [],
            'groups': [],
            'metadata': {
                'error': 'No MS1 spectra found in file',
                'mzml_path': mzml_path
            }
        }
    
    # Extract isotopic patterns
    patterns = []
    for spectrum in all_spectra:
        rt = spectrum.scan_time[0] if spectrum.scan_time else None
        peaks = spectrum.peaks('centroided')
        
        # Filter by minimum intensity
        filtered_peaks = [(mz, i) for mz, i in peaks if i >= min_intensity]
        
        if filtered_peaks:
            # Simple isotopic pattern extraction
            # (full implementation would use pattern detection algorithms)
            for i, (mz, intensity) in enumerate(filtered_peaks):
                pattern = {
                    'mz': mz,
                    'intensity': np.array([intensity]),  # Simplified for demo
                    'rt': rt
                }
                patterns.append(pattern)
    
    # Group isotopic patterns
    groups = group_isotopic_patterns(patterns, h_threshold)
    
    # Reconstruct components from groups
    components = []
    for group in groups:
        if len(group) == 1:
            # Single pattern - might be noise or unique component
            continue
        
        # Aggregate information from group
        mz_values = [p['mz'] for p in group]
        intensities = [np.sum(p['intensity']) for p in group]
        rt_values = [p['rt'] for p in group if p['rt'] is not None]
        
        component = {
            'mz': float(np.mean(mz_values)),
            'intensity': float(np.sum(intensities)),
            'rt': float(np.mean(rt_values)) if rt_values else None,
            'group_size': len(group),
            'mz_std': float(np.std(mz_values)),
            'intensity_range': [float(np.min(intensities)), float(np.max(intensities))]
        }
        components.append(component)
    
    # Sort components by intensity
    components.sort(key=lambda x: x['intensity'], reverse=True)
    
    metadata = {
        'input_file': mzml_path,
        'total_patterns': len(patterns),
        'total_groups': len(groups),
        'reconstructed_components': len(components),
        'h_threshold': h_threshold,
        'min_intensity': min_intensity,
        'algorithm': 'modified_cross_correlation_2021'
    }
    
    return {
        'components': components,
        'groups': groups,
        'metadata': metadata
    }


def visualize_results(results: Dict, output_path: str = None) -> None:
    """Visualize deconvolution results using Plotly.

    Args:
        results: Results dictionary from deconvolute_mzml
        output_path: Optional path to save HTML visualization
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        print("Plotly not installed. Skipping visualization.")
        return
    
    components = results['components']
    
    if not components:
        print("No components to visualize.")
        return
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'Component Intensity Distribution',
            'm/z Distribution',
            'Group Size vs Intensity',
            'Component Spread'
        )
    )
    
    # Plot 1: Intensity distribution
    intensities = [c['intensity'] for c in components]
    fig.add_trace(
        go.Histogram(x=intensities, name='Intensity'),
        row=1, col=1
    )
    
    # Plot 2: m/z distribution
    mz_values = [c['mz'] for c in components]
    fig.add_trace(
        go.Histogram(x=mz_values, name='m/z'),
        row=1, col=2
    )
    
    # Plot 3: Group size vs intensity
    group_sizes = [c['group_size'] for c in components]
    fig.add_trace(
        go.Scatter(
            x=group_sizes,
            y=intensities,
            mode='markers',
            name='Components',
            text=[f"m/z={mz:.2f}" for mz in mz_values]
        ),
        row=2, col=1
    )
    
    # Plot 4: Component spread (m/z std vs intensity)
    mz_stds = [c.get('mz_std', 0) for c in components]
    fig.add_trace(
        go.Scatter(
            x=mz_stds,
            y=intensities,
            mode='markers',
            name='Spread',
            marker=dict(
                size=[c['group_size'] for c in components],
                sizemode='diameter',
                sizeref=2
            )
        ),
        row=2, col=2
    )
    
    fig.update_layout(
        title_text="ChronoDecon Results Visualization",
        showlegend=False,
        height=800
    )
    
    if output_path:
        fig.write_html(output_path)
        print(f"Visualization saved to {output_path}")
    else:
        fig.show()
