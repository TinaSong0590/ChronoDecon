"""
Core deconvolution algorithms for ChronoDecon (v2 - mXcorr).

Implements the modified cross-correlation (mXcorr) algorithm from the
2021 Analytical Chemistry paper for 2D (time x m/z) deconvolution.

Key concepts from author's MATLAB code (mXcorrForSeparation):
- getDF(): FFT-based Dissymmetry Factor
- determineHThr(): Automatic threshold determination
- Wavelet time shift detection
- Greedy grouping by signal intensity
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
import warnings


def _get_df(a, b):
    """Dissymmetry Factor. DF = -sum(|Im(FFT(a)*conj(FFT(b))|) / sum(|Re(...)|).
    Author's getDF() in b_2DProcessing.m. DF~0 = same component, DF~-1 = different."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if len(a) == 0 or len(b) == 0: return -1.0
    n = max(len(a), len(b))
    if len(a) < n: a = np.pad(a, (0, n - len(a)))
    if len(b) < n: b = np.pad(b, (0, n - len(b)))
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0: return -1.0
    a, b = a / na, b / nb
    cross = np.fft.fft(a) * np.conj(np.fft.fft(b))
    im_s, re_s = np.sum(np.abs(np.imag(cross))), np.sum(np.abs(np.real(cross)))
    return -im_s / re_s if re_s != 0 else -1.0


def _determine_threshold(aic, num_shifts=20, exclude_first_n=4):
    """Auto-determine symmetry threshold (determineHThr.m).

    Implements two approaches from the 2021 paper SI Section 8:
    1. Primary: D-response function R_D(mu) exponential decay fitting
    2. Fallback: Self-correlation DF exponential decay + percentile-based

    The D-response approach (Eq. S3-S4):
      R_D(mu)|psi_a = D(psi_a(t), psi_a(t+mu))
    At large mu, this fits mono-exponential: A*exp(-k*mu) + c
    The decay rate k relates to the discrimination factor D and thus to H*_Thr.
    """
    if aic.shape[1] < 5 or aic.shape[0] < 10:
        return 0.88

    # --- Approach 1: D-response function (paper Eq. S3/S4) ---
    try:
        return _determine_threshold_d_response(aic)
    except Exception:
        pass

    # --- Approach 2: Self-correlation DF (original fallback) ---
    nc = min(aic.shape[1], 300)
    idx = np.linspace(0, aic.shape[1]-1, nc, dtype=int)
    dcoef = np.zeros(num_shifts)
    cnt = 0
    for c0 in idx:
        t = aic[:, c0]
        if np.sum(np.abs(t)) < 1e-10: continue
        for c1 in range(1, num_shifts+1):
            dcoef[c1-1] += _get_df(t, np.roll(t, c1)); cnt += 1
    if cnt == 0: return 0.88
    dcoef /= cnt
    x, y = np.arange(1, num_shifts+1, dtype=np.float64), dcoef
    if len(x) > exclude_first_n: x, y = x[exclude_first_n:], y[exclude_first_n:]
    if len(x) < 2: return 0.88
    try:
        c_est = np.min(y)
        ys = np.maximum(y - c_est, 1e-10)
        A = np.vstack([np.ones(len(x)), -x]).T
        b1 = np.linalg.lstsq(A, np.log(ys), rcond=None)[0][1]
        thr = 1.0 - b1
        if thr < 0.75 or thr > 0.99:
            df_vals = []
            rng = np.random.default_rng(42)
            for _ in range(min(200, aic.shape[1]**2)):
                i, j = rng.choice(aic.shape[1], 2, replace=False)
                if np.sum(np.abs(aic[:, i])) < 1e-10 or np.sum(np.abs(aic[:, j])) < 1e-10: continue
                d = _get_df(aic[:, i], aic[:, j])
                if d > -1.0: df_vals.append(d)
            if len(df_vals) > 10:
                thr = -np.percentile(df_vals, 10)
            else:
                thr = 0.88
    except Exception:
        thr = 0.88
    return max(0.75, min(0.98, thr))


def _determine_threshold_d_response(aic, max_shift=50, min_shift=5, n_samples=100):
    """D-response based threshold from paper SI Section 8 (Eq. S3/S4).

    For each ion chronogram psi_a, compute:
      R_D(mu)|psi_a = D(psi_a(t), psi_a(t+mu))
    Average across ions, fit mono-exponential at large mu.
    The decay rate k relates to H*_Thr.

    Rationale: When the instrument transfer function is stable, the response
    of D at increasing temporal shifts decays exponentially. The decay rate
    reflects how quickly temporally shifted signals become "different",
    which is exactly what HThr controls.
    """
    nt, nmz = aic.shape
    shifts = np.arange(1, min(max_shift + 1, nt // 2))
    avg_response = np.zeros(len(shifts))

    # Sample columns to average over
    col_indices = np.linspace(0, nmz - 1, min(n_samples, nmz), dtype=int)
    n_valid = 0

    for ci in col_indices:
        chron = aic[:, ci]
        if np.sum(np.abs(chron)) < 1e-10:
            continue

        response = np.zeros(len(shifts))
        for si, mu in enumerate(shifts):
            shifted = np.roll(chron, mu)
            response[si] = _get_df(chron, shifted)

        # Only use chronograms that show decay behavior
        if response[0] > -0.99 and np.std(response) > 0.01:
            avg_response += response
            n_valid += 1

    if n_valid < 3:
        return 0.88

    avg_response /= n_valid

    # Fit mono-exponential at large shifts: y = A*exp(-k*x) + c
    large_mask = shifts >= min_shift
    x_fit = shifts[large_mask].astype(np.float64)
    y_fit = avg_response[large_mask]

    if len(x_fit) < 3:
        return 0.88

    try:
        # Estimate asymptote
        c_est = np.min(y_fit[-min(5, len(y_fit)):])
        ys = np.maximum(y_fit - c_est, 1e-10)

        # Linear fit to log: log(y-c) = log(A) - k*x
        valid = ys > 1e-8
        if np.sum(valid) < 2:
            return 0.88

        x_v, y_v = x_fit[valid], np.log(ys[valid])
        coeffs = np.polyfit(x_v, y_v, 1)
        k = -coeffs[0]  # decay rate

        # Convert decay rate to threshold
        # Paper states: H*_Thr is derived from D-H relation (Eq. S4)
        # Empirically: higher k (faster decay) = easier separation = lower HThr
        # k typically ranges 0.01-0.5 for real data
        thr = 1.0 / (1.0 + 2.0 * k)  # empirical mapping

        if 0.75 <= thr <= 0.98:
            return float(thr)

        # If out of range, use cross-correlation percentile approach
        df_between = []
        rng = np.random.default_rng(42)
        for _ in range(min(300, nmz * (nmz - 1) // 2)):
            i, j = rng.choice(nmz, 2, replace=False)
            if np.sum(np.abs(aic[:, i])) < 1e-10 or np.sum(np.abs(aic[:, j])) < 1e-10:
                continue
            d = _get_df(aic[:, i], aic[:, j])
            if d > -1.0:
                df_between.append(d)
        if len(df_between) > 20:
            return float(max(0.75, min(0.98, -np.percentile(df_between, 10))))
    except Exception:
        pass

    return 0.88


def _detect_time_shift(aic, wavelet_width=5.0):
    """Wavelet time shift detection (b_processingData.m)."""
    xv = np.arange(0, 21, dtype=np.float64)
    loz = wavelet_width / ((xv - 10)**2 + wavelet_width**2)
    wavelet = -np.gradient(np.gradient(loz))
    return np.array([np.argmax(np.convolve(aic[:, c], wavelet, mode='same'))
                     for c in range(aic.shape[1])])


def _select_time_window(aic, t_shift, min_fraction=0.5):
    """Select main elution window peaks (b_processingData.m)."""
    us, cts = np.unique(t_shift, return_counts=True)
    mask = cts >= np.mean(cts) * min_fraction
    us, cts = us[mask], cts[mask]
    if len(us) == 0:
        return np.ones(aic.shape[1], dtype=bool), (0, aic.shape[0]-1)
    bounds = (int(np.min(us)), int(np.max(us)))
    margin = max(1, (bounds[1]-bounds[0])//10)
    return (t_shift >= bounds[0]-margin) & (t_shift <= bounds[1]+margin), bounds


def _group_peaks_by_symmetry(aic, peak_list, symmetry_threshold=0.90, buffer_size=4096):
    """Greedy grouping by DF (b_processingData.m 82-163)."""
    np_count = aic.shape[1]
    if np_count == 0: return []
    buf = np.zeros((buffer_size, np_count))
    buf[:min(buffer_size, aic.shape[0])] = aic[:min(buffer_size, aic.shape[0])]
    pms = np.max(aic, axis=0)
    avail = np.ones(np_count, dtype=bool)
    groups, assigned, prev = [], 0, 0
    while assigned < np_count:
        ai = np.where(avail)[0]
        if len(ai) == 0: break
        ri = ai[np.argmax(pms[ai])]
        ref = buf[:, ri]
        grp = []
        for c1 in range(np_count):
            if not avail[c1]: continue
            cross = np.fft.fft(ref) * np.conj(np.fft.fft(buf[:, c1]))
            symm = np.sum(np.abs(np.imag(cross))) / np.sum(np.abs(np.real(cross)))
            if symm < symmetry_threshold:
                grp.append(peak_list[c1].tolist()); avail[c1] = False; assigned += 1
        if grp: groups.append(np.array(grp))
        if assigned == prev: break
        prev = assigned
    return groups


def homologue_factor(psi_a, psi_b):
    """Isotopic H-factor (backward compatible)."""
    a, b = np.asarray(psi_a), np.asarray(psi_b)
    if len(a) == 0 or len(b) == 0: return 0.0
    ac, bc = a - np.mean(a), b - np.mean(b)
    nc = len(a) + len(b) - 1
    cc = np.real(np.fft.ifft(np.fft.fft(ac, n=nc) * np.conj(np.fft.fft(bc, n=nc))))
    na, nb = np.linalg.norm(ac), np.linalg.norm(bc)
    if na == 0 or nb == 0: return 0.0
    return max(0.0, min(1.0, np.max(cc) / (na * nb)))


def _extract_isotopic_distribution(intensities, peak_idx, n_isotopes=5):
    return intensities[max(0, peak_idx-n_isotopes):min(len(intensities), peak_idx+n_isotopes+1)]


def _find_peaks_in_spectrum(mz_values, intensities, min_intensity=0.01, min_peak_distance=3):
    mx = np.max(intensities)
    if mx == 0: return np.array([])
    ni = intensities / mx
    peaks = [i for i in range(1, len(intensities)-1)
             if ni[i] > min_intensity and ni[i] > ni[i-1] and ni[i] > ni[i+1]]
    filt, last = [], -min_peak_distance
    for p in peaks:
        if p - last >= min_peak_distance: filt.append(p); last = p
    return np.array(filt)


def _subtract_background(intensities, window_size=20, percentile=10.0):
    if len(intensities) == 0: return intensities
    bg = np.array([np.percentile(intensities[max(0,i-window_size//2):min(len(intensities),i+window_size//2+1)], percentile)
                   for i in range(len(intensities))])
    c = intensities - bg; c[c<0] = 0; return c


def _group_homologous_peaks(peak_indices, intensities, h_threshold=0.87, max_isotope_diff=15):
    """Legacy m/z-space grouping (backward compatible)."""
    if len(peak_indices) < 2: return [list(peak_indices)]
    groups, assigned = [], set()
    for i, pi in enumerate(peak_indices):
        if pi in assigned: continue
        cg = [pi]; assigned.add(pi)
        pa = _extract_isotopic_distribution(intensities, pi)
        cands = []
        for j, oi in enumerate(peak_indices):
            if oi in assigned or abs(oi-pi) > max_isotope_diff: continue
            h = homologue_factor(pa, _extract_isotopic_distribution(intensities, oi))
            if h >= h_threshold: cands.append((oi, h))
        cands.sort(key=lambda x: x[1], reverse=True)
        for oi, _ in cands:
            if oi not in assigned: cg.append(oi); assigned.add(oi)
        cg.sort(); groups.append(cg)
    return groups


def deconvolute_mzml(mzml_path, h_threshold=0.87, min_intensity=0.01, max_isotope_diff=15,
                     subtract_background=True, bg_window_size=20, bg_percentile=10.0,
                     symmetry_threshold=None, auto_threshold=True, min_peaks_per_group=1,
                     sort_by='max_signal', mz_bin_size=0.005, min_spectra_ratio=0.5):
    """2D deconvolution on mzML using mXcorr algorithm.

    New parameters:
        symmetry_threshold: DF threshold (None=auto)
        auto_threshold: Auto-determine from data
        min_peaks_per_group: Min peaks per component
        sort_by: 'max_signal' or 'num_peaks'
        mz_bin_size: m/z binning width in Da (default 0.005 = 5 mDa)
        min_spectra_ratio: Min fraction of spectra a peak must appear in
    """
    try:
        import pymzml
    except ImportError:
        return {'file_path': mzml_path, 'num_components': 0, 'components': [],
                'parameters': {}, 'status': 'error: pymzml not installed'}
    try:
        reader = pymzml.run.Reader(mzml_path)
    except Exception as e:
        return {'file_path': mzml_path, 'num_components': 0, 'components': [],
                'parameters': {}, 'status': f'error: {e}'}

    spectra = [s for s in reader if s.ms_level == 1]
    if not spectra:
        return {'file_path': mzml_path, 'num_components': 0, 'components': [],
                'parameters': {}, 'status': 'error: no MS1 spectra'}

    nt = len(spectra)

    # --- Build AIC with m/z binning ---
    # Step 1: Collect all (m/z, intensity) pairs
    all_mz = []
    all_int = []
    for s in spectra:
        p = s.centroidedPeaks
        if len(p) > 0:
            all_mz.append(np.array([x[0] for x in p]))
            all_int.append(np.array([x[1] for x in p]))

    if not all_mz:
        return {'file_path': mzml_path, 'num_components': 0, 'components': [],
                'parameters': {}, 'status': 'error: no peaks'}

    # Step 2: Create uniform m/z grid (binning)
    global_min_mz = min(m.min() for m in all_mz)
    global_max_mz = max(m.max() for m in all_mz)
    mz_edges = np.arange(global_min_mz, global_max_mz + mz_bin_size, mz_bin_size)
    mz_centers = (mz_edges[:-1] + mz_edges[1:]) / 2
    nmz = len(mz_centers)

    # Step 3: Build AIC matrix by binning each spectrum
    aic = np.zeros((nt, nmz))
    for t_idx, (spec_mz, spec_int) in enumerate(zip(all_mz, all_int)):
        bin_idx = np.digitize(spec_mz, mz_edges) - 1
        valid = (bin_idx >= 0) & (bin_idx < nmz)
        np.add.at(aic[t_idx], bin_idx[valid], spec_int[valid])

    # Step 4: Filter by minimum presence (peak must appear in >= min_spectra_ratio of spectra)
    presence = np.sum(aic > 0, axis=0) / nt
    presence_mask = presence >= min_spectra_ratio
    aic = aic[:, presence_mask]
    mz_centers = mz_centers[presence_mask]
    nmz = aic.shape[1]

    if nmz < 2:
        return {'file_path': mzml_path, 'num_components': 0, 'components': [],
                'parameters': {}, 'status': f'error: only {nmz} m/z bins after filtering'}

    # Peak list: [mz, max_intensity]
    pl = np.column_stack([mz_centers, np.max(aic, axis=0)])

    # Filter by intensity
    mx = np.max(pl[:, 1])
    if mx == 0:
        return {'file_path': mzml_path, 'num_components': 0, 'components': [],
                'parameters': {}, 'status': 'error: zero intensities'}
    im = (pl[:, 1] / mx) >= min_intensity
    aic, pl, mz_centers = aic[:, im], pl[im, :], mz_centers[im]
    nmz = aic.shape[1]
    if nmz < 2:
        return {'file_path': mzml_path, 'num_components': 0, 'components': [],
                'parameters': {}, 'status': 'error: too few peaks'}

    # Time shift & window
    ts = _detect_time_shift(aic)
    sm, tb = _select_time_window(aic, ts)
    if np.sum(sm) < 2:
        sm = np.ones(aic.shape[1], dtype=bool); tb = (0, nt-1)
    aic_sel = aic[:, sm]
    pl_sel = pl[sm, :]
    mz_sel = mz_centers[sm]

    # Background subtraction
    if subtract_background:
        hw = max(1, bg_window_size // 2)
        for c1 in range(aic_sel.shape[1]):
            ch = aic_sel[:, c1]
            bl = np.array([np.percentile(ch[max(0,i-hw):min(len(ch),i+hw+1)], bg_percentile)
                           for i in range(len(ch))])
            aic_sel[:, c1] = np.maximum(ch - bl, 0)

    # Z-score normalization
    aic_zs = np.zeros_like(aic_sel)
    for c1 in range(aic_sel.shape[1]):
        c = aic_sel[:, c1]; s = np.std(c)
        if s > 0: aic_zs[:, c1] = (c - np.mean(c)) / s

    pl_proc = np.column_stack([pl_sel[:, 0], np.max(aic_sel, axis=0)])

    # Auto threshold
    ti = {}
    if auto_threshold and symmetry_threshold is None:
        try:
            symmetry_threshold = _determine_threshold(aic_zs)
            ti = {'method': 'auto', 'auto_threshold': round(symmetry_threshold, 4),
                  'note': 'Exponential decay fit on self-correlations'}
        except Exception as e:
            symmetry_threshold = 0.88
            ti = {'method': 'fallback', 'auto_threshold': 0.88, 'note': str(e)}
    elif symmetry_threshold is not None:
        ti = {'method': 'manual', 'auto_threshold': symmetry_threshold, 'note': 'User-specified'}
    else:
        symmetry_threshold = 0.88
        ti = {'method': 'default', 'auto_threshold': 0.88, 'note': 'Default'}

    # Group
    groups = _group_peaks_by_symmetry(aic_zs, pl_proc, symmetry_threshold)
    groups = [g for g in groups if len(g) >= min_peaks_per_group]
    if sort_by == 'num_peaks':
        groups.sort(key=lambda g: len(g), reverse=True)
    else:
        groups.sort(key=lambda g: np.max(g[:, 1]) if len(g) > 0 else 0, reverse=True)

    # Build result
    comps = []
    for idx, g in enumerate(groups):
        mz_v, int_v = g[:, 0].tolist(), g[:, 1].tolist()
        hfs = []
        if len(g) >= 2:
            sg = g[g[:, 0].argsort()]
            for i in range(len(sg)-1):
                if 0.9 < abs(sg[i+1, 0] - sg[i, 0]) < 1.1:
                    ia = np.argmin(np.abs(mz_sel - sg[i, 0]))
                    ib = np.argmin(np.abs(mz_sel - sg[i+1, 0]))
                    hfs.append(round(homologue_factor(aic_sel[:, ia], aic_sel[:, ib]), 4))
        comps.append({
            'component_id': idx+1,
            'mz_values': [round(m, 6) for m in mz_v],
            'intensities': [round(i, 4) for i in int_v],
            'num_peaks': len(g),
            'max_intensity': round(float(np.max(int_v)), 4) if int_v else 0,
            'total_intensity': round(float(np.sum(int_v)), 4) if int_v else 0,
            'h_factors': hfs,
            'avg_mz': round(float(np.mean(mz_v)), 4) if mz_v else 0})

    return {
        'file_path': mzml_path, 'num_components': len(comps), 'num_spectra': nt,
        'num_mz_bins': nmz, 'num_selected_mz': int(np.sum(sm)),
        'time_bounds': list(tb), 'components': comps,
        'parameters': {'symmetry_threshold': symmetry_threshold, 'min_intensity': min_intensity,
                       'subtract_background': subtract_background, 'bg_window_size': bg_window_size,
                       'bg_percentile': bg_percentile, 'auto_threshold': auto_threshold,
                       'mz_bin_size': mz_bin_size, 'min_spectra_ratio': min_spectra_ratio,
                       'sort_by': sort_by},
        'threshold_info': ti, 'status': 'success'}


def visualize_deconvolution(result, output_path=None):
    """Visualize with plotly."""
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        warnings.warn("plotly not installed"); return
    comps = result.get('components', [])[:20]
    if not comps: print("No components."); return
    fig = make_subplots(rows=len(comps), cols=1,
                        subplot_titles=[f"C{c['component_id']} ({c['num_peaks']}p, mz~{c['avg_mz']:.1f})" for c in comps],
                        vertical_spacing=0.04)
    colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf',
              '#999999','#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494',
              '#b3b3b3','#666666','#1b9e77','#d95f02']
    for i, c in enumerate(comps):
        fig.add_trace(go.Bar(x=c['mz_values'], y=c['intensities'], marker_color=colors[i%len(colors)],
                             showlegend=False, width=0.8), row=i+1, col=1)
    fig.update_layout(title=f"ChronoDecon mXcorr v2 - {result['num_components']} components",
                      height=250*len(comps), xaxis_title="m/z", yaxis_title="Intensity")
    if output_path:
        fig.write_html(output_path); print(f"Saved: {output_path}")
    else:
        fig.show()


def export_enriched_mgf(mzml_path, results, output_path,
                         min_peaks_per_spectrum=3, top_n_components=None,
                         intensity_percentile=95):
    """Export enriched MGF with full component spectra for library matching.

    Unlike the basic MGF which only lists grouped m/z bin centers,
    this function extracts complete peak data from the original mzML file
    at each component's elution time, producing richer spectra.

    Parameters
    ----------
    mzml_path : str
        Original mzML file path
    results : dict
        Deconvolution results from deconvolute_mzml()
    output_path : str
        Output MGF file path
    min_peaks_per_spectrum : int
        Minimum peaks to include a component (default: 3)
    top_n_components : int, optional
        Only export top N components by intensity
    intensity_percentile : float
        Use peaks above this percentile of intensity (0-100)
    """
    import pymzml

    reader = pymzml.run.Reader(mzml_path)
    spectra = [s for s in reader if s.ms_level == 1]
    components = results.get('components', [])

    if not spectra or not components:
        return 0

    # Sort by total intensity (descending), take top N
    components = sorted(components, key=lambda c: c.get('total_intensity', 0), reverse=True)
    if top_n_components:
        components = components[:top_n_components]

    nt = len(spectra)
    written = 0

    with open(output_path, 'w') as f:
        for comp in components:
            comp_id = comp['component_id']
            mz_values = comp.get('mz_values', [])
            intensities = comp.get('intensities', [])

            if len(mz_values) < min_peaks_per_spectrum:
                continue

            # Determine the component's elution time window:
            # Use the median m/z as the representative peak
            median_mz = float(np.median(mz_values))

            # Find which spectra have the strongest signal near median_mz
            best_spec_idx = 0
            best_signal = 0

            for si, spec in enumerate(spectra):
                peaks = spec.centroidedPeaks
                if len(peaks) == 0:
                    continue
                # Sum intensity of peaks within 1 Da of median_mz
                for mz_p, int_p in peaks:
                    if abs(mz_p - median_mz) < 1.0:
                        if int_p > best_signal:
                            best_signal = int_p
                            best_spec_idx = si

            # Extract full spectrum from best scan + nearby scans
            # Average across a window of +/- 5 scans around the best
            scan_window = 5
            start_idx = max(0, best_spec_idx - scan_window)
            end_idx = min(nt, best_spec_idx + scan_window + 1)

            peak_dict = {}
            for si in range(start_idx, end_idx):
                peaks = spectra[si].centroidedPeaks
                if len(peaks) == 0:
                    continue
                for mz_p, int_p in peaks:
                    # Round to 3 decimals for merging
                    mz_key = round(mz_p, 3)
                    peak_dict[mz_key] = peak_dict.get(mz_key, 0) + int_p

            if not peak_dict:
                continue

            # Convert to sorted lists
            peak_mz = sorted(peak_dict.keys())
            peak_int = [peak_dict[m] for m in peak_mz]

            # Filter by intensity percentile
            if len(peak_int) > min_peaks_per_spectrum:
                threshold = np.percentile(peak_int, intensity_percentile)
                filtered = [(m, i) for m, i in zip(peak_mz, peak_int) if i >= threshold]
            else:
                filtered = list(zip(peak_mz, peak_int))

            if len(filtered) < min_peaks_per_spectrum:
                continue

            # Normalize intensities to 0-999
            max_int = max(i for _, i in filtered)
            if max_int > 0:
                normalized = [(m, int(i / max_int * 999)) for m, i in filtered]
            else:
                normalized = filtered

            # Set PEPMASS: prefer highest-intensity peak (often molecular ion in ESI)
            pepmass_mz = max(normalized, key=lambda x: x[1])[0]

            # Write MGF entry
            f.write("BEGIN IONS\n")
            f.write(f"TITLE=Component_{comp_id}\n")
            f.write(f"PEPMASS={pepmass_mz:.6f}\n")
            f.write(f"SCANS={best_spec_idx+1}\n")
            f.write(f"MSLEVEL=2\n")
            f.write(f"COMMENT=Deconvoluted by ChronoDecon mXcorr, {len(normalized)} peaks\n")
            for mz, intensity in normalized:
                f.write(f"{mz:.6f} {intensity}\n")
            f.write("END IONS\n")
            written += 1

    return written


def export_simple_mgf(results, output_path, min_peaks=1):
    """Simple MGF export from deconvolution results (backward compatible)."""
    components = results.get('components', [])
    written = 0
    with open(output_path, 'w') as f:
        for comp in components:
            mz_values = comp.get('mz_values', [])
            intensities = comp.get('intensities', [])
            if len(mz_values) < min_peaks:
                continue
            max_int = max(intensities) if intensities else 1
            if max_int == 0:
                continue
            f.write("BEGIN IONS\n")
            f.write(f"TITLE=Component_{comp['component_id']}\n")
            pepmass = float(np.mean(mz_values)) if mz_values else 0
            pepmass_mz = pepmass
            f.write(f"SCANS=1\n")
            f.write(f"MSLEVEL=2\n")
            for mz, inten in zip(mz_values, intensities):
                norm_int = int(inten / max_int * 999)
                f.write(f"{mz:.6f} {norm_int}\n")
            f.write("END IONS\n")
            written += 1
    return written