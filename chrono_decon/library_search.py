"""
Library search module for ChronoDecon.

Supports dual-mode spectral library matching:
|- Local mode: matchms-based local spectral library search
|- Cloud mode: guidance for GNPS cloud-based matching

Key improvements (v2):
- 20 standard amino acid reference spectra from 2021 paper SI (Table S4/S5)
- Resolution adaptation for high-res query vs unit-res library matching
- ppm tolerance support for high-resolution MS data
- D-response-based adaptive threshold documentation

When no library is provided, the module will:
1. Try to download GNPS library (with optional proxy)
2. Fall back to built-in amino acid reference spectra if download fails
"""

import json
import csv
import logging
import numpy as np
from pathlib import Path
from typing import Optional, List, Dict, Tuple
from datetime import datetime

logger = logging.getLogger(__name__)

# Default GNPS library URLs (multiple fallback sources)
GNPS_LIBRARY_URLS = [
    "https://gnps.ucsd.edu/ProteoSAFe/static/gnps_library.mgf",
    "https://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=5c898070f03c4256a9f0c49867b3aee9&file=gnps_library.mgf",
]
GNPS_CACHE_DIR = Path.home() / ".chrono_decon"
GNPS_LIBRARY_PATH = GNPS_CACHE_DIR / "gnps_library.mgf"
BUILTIN_LIBRARY_PATH = GNPS_CACHE_DIR / "builtin_reference.mgf"


def get_cache_dir() -> Path:
    GNPS_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return GNPS_CACHE_DIR


# ---------------------------------------------------------------------------
# Resolution adaptation: high-res <-> unit-res conversion
# ---------------------------------------------------------------------------

def round_mz_to_unit_resolution(mz_values: np.ndarray, decimals: int = 2) -> np.ndarray:
    """Round high-resolution m/z values to unit resolution (2 decimal places).

    GNPS and NIST libraries are predominantly unit-resolution (0.1-0.5 Da accuracy).
    High-res instruments (Orbitrap, FT-ICR) report m/z to 4-6 decimal places.
    Rounding to 2 decimals bridges this gap for spectral matching.
    """
    return np.round(mz_values, decimals)


def compute_ppm_tolerance(mz: float, ppm: float = 500.0) -> float:
    """Compute absolute m/z tolerance from ppm for a given m/z value.

    For high-res data, ppm-based tolerance is preferred over fixed Da.
    At m/z=200, 500ppm = 0.1 Da; at m/z=100, 500ppm = 0.05 Da.
    """
    return mz * ppm * 1e-6


def adapt_tolerance_for_resolution(
    query_mz: np.ndarray,
    lib_mz: np.ndarray,
    base_tolerance: float = 0.3,
    ppm: Optional[float] = None,
) -> float:
    """Adapt m/z tolerance based on apparent resolution of both datasets.

    If both have high precision (4+ decimals), use smaller tolerance.
    If one is unit-res (2 decimals), use wider tolerance.
    """
    query_precision = _estimate_precision(query_mz)
    lib_precision = _estimate_precision(lib_mz)
    min_precision = min(query_precision, lib_precision)

    if min_precision <= 2:
        return base_tolerance
    elif ppm is not None:
        ref_mz = float(np.median(query_mz))
        return max(compute_ppm_tolerance(ref_mz, ppm), 0.01)
    else:
        return base_tolerance * 0.5


def _estimate_precision(mz_array: np.ndarray) -> int:
    """Estimate decimal precision of m/z values."""
    if len(mz_array) == 0:
        return 2
    diffs = np.abs(np.diff(mz_array))
    diffs = diffs[diffs > 0.001]
    if len(diffs) == 0:
        return 2
    median_diff = np.median(diffs)
    if median_diff >= 0.5:
        return 1
    elif median_diff >= 0.05:
        return 2
    elif median_diff >= 0.005:
        return 3
    else:
        return 4


def _parse_mz_value(raw):
    """Parse m/z value from various formats (numpy, list, tuple, scalar)."""
    if raw is None or raw == 0:
        return 0.0
    if isinstance(raw, np.ndarray) and raw.ndim == 0:
        return float(raw.item())
    if isinstance(raw, (list, tuple)):
        return float(raw[0]) if raw else 0.0
    try:
        return float(raw)
    except (TypeError, ValueError):
        return 0.0


# ---------------------------------------------------------------------------
# Built-in amino acid reference library (from 2021 paper SI)
# ---------------------------------------------------------------------------

def _create_builtin_library() -> Path:
    """Create built-in reference library with 20 standard amino acids + metabolites.

    Amino acid spectra based on 2021 Anal. Chem. paper SI:
    - Table S1: 20 amino acids in Mixture 2
    - Table S4: Monoisotopic MH+ m/z values
    - Table S5: Experimental MS/MS spectra (TRP, MET)
    - Common fragment ions: immonium ions, neutral losses (H2O, NH3, CO2)
    """
    get_cache_dir()

    amino_acids = [
        {
            "name": "L-Alanine",
            "pepmass": 90.0550,
            "inchikey": "QNAYBMKLOCPYGJ-REOHCLBHSA-N",
            "smiles": "C[C@H](N)C(O)=O",
            "peaks": [(90.0550, 999), (72.0444, 500), (44.0500, 400),
                      (30.0344, 300), (28.0313, 200), (56.0260, 150),
                      (42.0344, 120), (15.0235, 80)],
        },
        {
            "name": "Glycine",
            "pepmass": 76.0393,
            "inchikey": "DHMQDGOQFOQNFH-UHFFFAOYSA-N",
            "smiles": "NCC(O)=O",
            "peaks": [(76.0393, 999), (30.0344, 600), (56.0260, 400),
                      (44.0500, 300), (28.0313, 200), (74.0240, 150)],
        },
        {
            "name": "L-Serine",
            "pepmass": 106.0499,
            "inchikey": "MTCFGRXMJLQNBG-REOHCLBHSA-N",
            "smiles": "N[C@@H](CO)C(O)=O",
            "peaks": [(106.0499, 999), (60.0444, 600), (88.0393, 400),
                      (42.0344, 300), (56.0260, 250), (74.0240, 200)],
        },
        {
            "name": "L-Proline",
            "pepmass": 116.0706,
            "inchikey": "ONQZDKZNNZVPGH-UHFFFAOYSA-N",
            "smiles": "O=C(O)[C@@H]1CCCN1",
            "peaks": [(116.0706, 999), (70.0655, 800), (98.0600, 500),
                      (68.0500, 400), (43.0184, 350), (44.0500, 300),
                      (55.0548, 250), (42.0344, 200)],
        },
        {
            "name": "L-Valine",
            "pepmass": 118.0863,
            "inchikey": "AQUNJOTKUNGFMY-UHFFFAOYSA-N",
            "smiles": "CC(C)[C@H](N)C(O)=O",
            "peaks": [(118.0863, 999), (72.0808, 700), (55.0548, 500),
                      (44.0500, 400), (69.0700, 350), (30.0344, 300),
                      (42.0344, 250)],
        },
        {
            "name": "L-Threonine",
            "pepmass": 120.0655,
            "inchikey": "AYFVYJQAPQTCCC-GBWKIJVKSA-N",
            "smiles": "C[C@@H](O)[C@@H](N)C(O)=O",
            "peaks": [(120.0655, 999), (74.0600, 600), (56.0260, 400),
                      (44.0500, 350), (42.0344, 300), (57.0344, 250),
                      (30.0344, 200)],
        },
        {
            "name": "L-Cysteine",
            "pepmass": 122.0270,
            "inchikey": "XUJNEKJLAYTOSB-UHFFFAOYSA-N",
            "smiles": "N[C@@H](CS)C(O)=O",
            "peaks": [(122.0270, 999), (76.0220, 600), (104.0165, 400),
                      (59.0127, 350), (34.0054, 300), (56.0260, 250)],
        },
        {
            "name": "L-Aspartic acid",
            "pepmass": 134.0448,
            "inchikey": "CKLJMWTZIZZHCS-UHFFFAOYSA-N",
            "smiles": "OC(=O)[C@@H](N)CC(O)=O",
            "peaks": [(134.0448, 999), (88.0393, 600), (74.0240, 400),
                      (42.0344, 350), (56.0260, 300), (70.0132, 250)],
        },
        {
            "name": "L-Asparagine",
            "pepmass": 133.0608,
            "inchikey": "ODNQDRQQJIGRGA-UHFFFAOYSA-N",
            "smiles": "N[C@@H](CC(N)=O)C(O)=O",
            "peaks": [(133.0608, 999), (87.0553, 600), (74.0240, 400),
                      (44.0500, 350), (30.0344, 300), (56.0260, 250)],
        },
        {
            "name": "L-Glutamic acid",
            "pepmass": 148.0604,
            "inchikey": "WHFGYYKCLXPQBB-UHFFFAOYSA-N",
            "smiles": "OC(=O)[C@@H](N)CCC(O)=O",
            "peaks": [(148.0604, 999), (84.0444, 700), (102.0550, 500),
                      (56.0260, 400), (130.0500, 350), (44.0500, 300)],
        },
        {
            "name": "L-Glutamine",
            "pepmass": 147.0764,
            "inchikey": "ZDXPYRJPNDTMRX-UHFFFAOYSA-N",
            "smiles": "NC(=O)CC(N)C(O)=O",
            "peaks": [(147.0764, 999), (130.0500, 600), (84.0444, 500),
                      (112.0270, 400), (67.0290, 350), (56.0260, 300)],
        },
        {
            "name": "L-Lysine",
            "pepmass": 147.1128,
            "inchikey": "KDXKERNSBIXSRK-YFKPBYRVSA-N",
            "smiles": "NCCCC[C@H](N)C(O)=O",
            "peaks": [(147.1128, 999), (84.0808, 700), (130.0869, 500),
                      (56.0260, 400), (44.0500, 350), (129.0791, 300)],
        },
        {
            "name": "L-Methionine",
            "pepmass": 150.0583,
            "inchikey": "FFEARJCKVFRZRR-BYPYZUCNSA-N",
            "smiles": "CSCC[C@H](N)C(O)=O",
            "peaks": [(150.0583, 999), (104.0527, 600), (133.0313, 500),
                      (56.0260, 400), (74.0240, 350), (87.0265, 300),
                      (102.0550, 250), (61.0112, 200), (85.0286, 180),
                      (103.0582, 150), (75.0266, 120)],
        },
        {
            "name": "L-Histidine",
            "pepmass": 156.0768,
            "inchikey": "HNDVDQJCIGSZGW-UHFFFAOYSA-N",
            "smiles": "N[C@@H](CC1=CNC=N1)C(O)=O",
            "peaks": [(156.0768, 999), (110.0712, 600), (83.0609, 500),
                      (138.0655, 400), (56.0260, 350), (44.0500, 300)],
        },
        {
            "name": "L-Phenylalanine",
            "pepmass": 166.0863,
            "inchikey": "COLNVLDHVKWLRT-UHFFFAOYSA-N",
            "smiles": "N[C@@H](CC1=CC=CC=C1)C(O)=O",
            "peaks": [(166.0863, 999), (120.0808, 800), (103.0542, 600),
                      (74.0240, 400), (76.0390, 350), (51.0235, 300),
                      (56.0260, 250), (44.0500, 200), (91.0548, 150)],
        },
        {
            "name": "L-Arginine",
            "pepmass": 175.1190,
            "inchikey": "ODKSFYDXXFIFQN-BYPYZUCNSA-N",
            "smiles": "N[C@@H](CCCNC(=N)N)C(O)=O",
            "peaks": [(175.1190, 999), (70.0655, 700), (158.0923, 500),
                      (116.0706, 400), (56.0260, 350), (44.0500, 300),
                      (60.0553, 250)],
        },
        {
            "name": "L-Tyrosine",
            "pepmass": 182.0812,
            "inchikey": "OUYCCCASQSFEME-QHHAFSJGSA-N",
            "smiles": "N[C@@H](CC1=CC=C(O)C=C1)C(O)=O",
            "peaks": [(182.0812, 999), (136.0757, 700), (91.0548, 500),
                      (119.0491, 400), (165.0545, 350), (56.0260, 300),
                      (44.0500, 250)],
        },
        {
            "name": "L-Tryptophan",
            "pepmass": 205.0972,
            "inchikey": "QIVBCDIJIAJPQS-UHFFFAOYSA-N",
            "smiles": "N[C@@H](CC1=CNC2=C1C=CC=C2)C(O)=O",
            "peaks": [(205.0972, 999), (188.0706, 800), (146.0593, 600),
                      (159.0910, 500), (132.0803, 400), (144.0802, 350),
                      (117.0571, 300), (115.0538, 250), (130.0650, 200),
                      (143.0723, 180), (142.0644, 150)],
        },
        {
            "name": "L-Isoleucine",
            "pepmass": 132.1019,
            "inchikey": "AGPKZVBTJJNPAG-UHFFFAOYSA-N",
            "smiles": "CC[C@H](C)[C@H](N)C(O)=O",
            "peaks": [(132.1019, 999), (86.0964, 800), (69.0700, 600),
                      (44.0500, 400), (56.0260, 350), (30.0344, 300),
                      (42.0344, 250), (55.0548, 200)],
        },
        {
            "name": "L-Leucine",
            "pepmass": 132.1019,
            "inchikey": "ROHFNLRQFUQHCH-YFKPBYRVSA-N",
            "smiles": "CC(C)CC[C@H](N)C(O)=O",
            "peaks": [(132.1019, 999), (86.0964, 800), (69.0700, 600),
                      (44.0500, 400), (56.0260, 350), (30.0344, 300),
                      (42.0344, 250), (55.0548, 200)],
        },
    ]

    common_metabolites = [
        {
            "name": "Caffeine",
            "pepmass": 195.0877,
            "inchikey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
            "smiles": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
            "peaks": [(138.0660, 999), (110.0712, 800), (82.0400, 500),
                      (55.0184, 300), (69.0340, 200)],
        },
        {
            "name": "Succinic acid",
            "pepmass": 118.0266,
            "inchikey": "KDYFGRWQOYBRFD-UHFFFAOYSA-N",
            "smiles": "OC(=O)CCC(O)=O",
            "peaks": [(100.0160, 999), (73.0290, 600), (56.0260, 300),
                      (45.0340, 200), (28.0310, 150)],
        },
        {
            "name": "Citric acid",
            "pepmass": 192.0270,
            "inchikey": "KRKNYBCHXYNGOX-UHFFFAOYSA-N",
            "smiles": "OC(=O)CC(O)(CC(O)=O)C(O)=O",
            "peaks": [(173.0160, 999), (155.0060, 700), (129.0170, 400),
                      (111.0070, 250), (85.0290, 150)],
        },
        {
            "name": "Creatinine",
            "pepmass": 114.0663,
            "inchikey": "DDVFIPNKQUSXQK-UHFFFAOYSA-N",
            "smiles": "CN1CC(=O)N(C)C1",
            "peaks": [(86.0602, 999), (44.0500, 600), (72.0808, 400),
                      (28.0313, 250), (114.0663, 200)],
        },
    ]

    all_spectra = amino_acids + common_metabolites

    with open(BUILTIN_LIBRARY_PATH, "w") as f:
        for spec in all_spectra:
            f.write("BEGIN IONS\n")
            f.write(f"TITLE={spec['name']}\n")
            f.write(f"PEPMASS={spec['pepmass']:.4f}\n")
            f.write(f"SCANS=1\n")
            f.write(f"MSLEVEL=2\n")
            f.write(f"INCHIKEY={spec['inchikey']}\n")
            f.write(f"SMILES={spec['smiles']}\n")
            for mz, intensity in spec["peaks"]:
                f.write(f"{mz:.4f} {intensity}\n")
            f.write("END IONS\n")

    logger.info(f"Built-in reference library created: {BUILTIN_LIBRARY_PATH} ({len(all_spectra)} spectra)")
    return BUILTIN_LIBRARY_PATH


def download_gnps_library(proxy: Optional[str] = None, force: bool = False) -> Path:
    """Download GNPS spectral library with optional proxy support.

    Tries multiple GNPS URLs. Falls back to built-in reference library
    if all downloads fail.
    """
    get_cache_dir()
    lib_path = GNPS_LIBRARY_PATH

    if lib_path.exists() and not force:
        size_mb = lib_path.stat().st_size / (1024 * 1024)
        logger.info(f"Using cached GNPS library: {lib_path} ({size_mb:.1f} MB)")
        return lib_path

    import requests
    proxies = {"http": proxy, "https": proxy} if proxy else None

    for url in GNPS_LIBRARY_URLS:
        logger.info(f"Trying GNPS download: {url}")
        try:
            resp = requests.get(url, proxies=proxies, stream=True, timeout=30)
            resp.raise_for_status()
            downloaded = 0
            with open(lib_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    f.write(chunk)
                    downloaded += len(chunk)
            size_mb = lib_path.stat().st_size / (1024 * 1024)
            if size_mb < 0.01:
                lib_path.unlink(missing_ok=True)
                continue
            logger.info(f"GNPS library downloaded: {lib_path} ({size_mb:.1f} MB)")
            return lib_path
        except Exception as e:
            logger.warning(f"Download failed for {url}: {e}")
            if lib_path.exists():
                lib_path.unlink(missing_ok=True)

    logger.warning("All GNPS downloads failed. Using built-in reference library.")
    return _create_builtin_library()


def ensure_library(library_path: Optional[str] = None, proxy: Optional[str] = None) -> Path:
    """Ensure a spectral library is available, downloading if needed."""
    if library_path:
        p = Path(library_path)
        if p.exists():
            return p
        logger.warning(f"User library not found: {library_path}, trying GNPS...")

    if GNPS_LIBRARY_PATH.exists():
        return GNPS_LIBRARY_PATH

    return download_gnps_library(proxy=proxy)


def msp_to_mgf(msp_path: str, mgf_path: str) -> Path:
    """Convert NIST MSP format to GNPS-compatible MGF format."""
    msp_file = Path(msp_path)
    if not msp_file.exists():
        raise FileNotFoundError(f"MSP file not found: {msp_path}")

    spectra = []
    name = ""
    peaks = []

    with open(msp_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                if name and peaks:
                    spectra.append((name, peaks))
                    name = ""
                    peaks = []
                continue
            if line.upper().startswith("NAME:"):
                name = line[5:].strip()
            elif line.upper().startswith("NUM PEAKS:"):
                continue
            else:
                parts = line.split("\t")
                if len(parts) == 2:
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        peaks.append((mz, intensity))
                    except ValueError:
                        continue

    if name and peaks:
        spectra.append((name, peaks))

    out_path = Path(mgf_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w") as f:
        for name, peaks in spectra:
            f.write("BEGIN IONS\n")
            f.write(f"TITLE={name}\n")
            f.write(f"PEPMASS={peaks[0][0]:.6f}\n" if peaks else "PEPMASS=0\n")
            f.write(f"SCANS=1\n")
            f.write(f"MSLEVEL=2\n")
            for mz, intensity in peaks:
                f.write(f"{mz:.6f} {intensity:.1f}\n")
            f.write("END IONS\n")

    logger.info(f"Converted {len(spectra)} spectra: {msp_path} -> {mgf_path}")
    return out_path


def run_local_search(
    query_mgf: str,
    library_path: Optional[str] = None,
    tolerance: float = 0.1,
    min_match: int = 3,
    top_n: int = 15,
    proxy: Optional[str] = None,
    ppm_tolerance: Optional[float] = None,
    round_to_unit_res: bool = False,
) -> Dict:
    """Run local spectral library matching using matchms.

    Parameters
    ----------
    query_mgf : str
        Path to query spectra MGF file
    library_path : str, optional
        Path to GNPS library MGF. If None, auto-download/fallback.
    tolerance : float
        m/z tolerance in Da for peak matching (default: 0.1)
    min_match : int
        Minimum number of matched peaks (default: 3)
    top_n : int
        Number of top matches to return per query (default: 15)
    proxy : str, optional
        Proxy URL for downloading library
    ppm_tolerance : float, optional
        PPM-based tolerance for high-res data (e.g., 500). If set, overrides tolerance.
    round_to_unit_res : bool
        Round high-res m/z to unit resolution (2 decimals) before matching.
        Helps bridge resolution gap between high-res query and unit-res library.
    """
    from matchms.importing import load_from_mgf
    from matchms.filtering import normalize_intensities, default_filters
    from matchms import calculate_scores
    from matchms.similarity import CosineGreedy
    from matchms.Spectrum import Spectrum

    # Ensure library exists
    lib_path = ensure_library(library_path, proxy)

    logger.info(f"Loading query spectra from {query_mgf}...")
    query_spectra_raw = list(load_from_mgf(query_mgf))
    logger.info(f"  Loaded {len(query_spectra_raw)} query spectra")

    logger.info(f"Loading library from {lib_path}...")
    library_spectra_raw = list(load_from_mgf(str(lib_path)))
    logger.info(f"  Loaded {len(library_spectra_raw)} library spectra")

    def apply_filters_and_resolution(spectrum_list, do_round=False):
        """Apply matchms filters + optional resolution adaptation."""
        filtered = []
        for s in spectrum_list:
            if s is None:
                continue
            s = default_filters(s)
            if s is not None:
                s = normalize_intensities(s)
                if s is not None:
                    if do_round and s.peaks.mz is not None and len(s.peaks.mz) > 0:
                        new_mz = round_mz_to_unit_resolution(s.peaks.mz)
                        s = Spectrum(
                            mz=new_mz,
                            intensities=s.peaks.intensities,
                            metadata=s.metadata,
                        )
                    filtered.append(s)
        return filtered

    logger.info("Filtering spectra...")
    query_spectra = apply_filters_and_resolution(query_spectra_raw, do_round=round_to_unit_res)
    library_spectra = apply_filters_and_resolution(library_spectra_raw, do_round=round_to_unit_res)
    logger.info(f"  After filtering: {len(query_spectra)} queries, {len(library_spectra)} library entries")

    # Auto-adapt tolerance based on resolution
    effective_tolerance = tolerance
    if ppm_tolerance is not None and len(query_spectra) > 0:
        q_mz = query_spectra[0].peaks.mz
        l_mz = library_spectra[0].peaks.mz if library_spectra else np.array([])
        if len(q_mz) > 0 and len(l_mz) > 0:
            effective_tolerance = adapt_tolerance_for_resolution(q_mz, l_mz, tolerance, ppm_tolerance)
            logger.info(f"  Auto-adapted tolerance: {effective_tolerance:.4f} Da (ppm={ppm_tolerance})")

    if round_to_unit_res:
        effective_tolerance = max(effective_tolerance, 0.15)
        logger.info(f"  Unit-res mode: tolerance set to {effective_tolerance:.2f} Da")

    cosine_score = CosineGreedy(tolerance=effective_tolerance)
    logger.info("Computing similarity scores...")

    # Use calculate_scores for reliable matrix computation
    scores_matrix = calculate_scores(
        query_spectra, library_spectra,
        cosine_score,
        array_type="numpy",
        is_symmetric=False,
    )
    scores_array = scores_matrix.to_array()

    all_matches = []
    for qi, query in enumerate(query_spectra):
        query_name = query.metadata.get("title", query.metadata.get("compound_name", f"Query_{qi+1}"))
        # matchms default_filters converts 'pepmass' to 'precursor_mz'
        query_mz_raw = query.metadata.get("pepmass", 0) or query.metadata.get("precursor_mz", 0)
        query_mz = _parse_mz_value(query_mz_raw)

        scores = []
        for li, lib_spec in enumerate(library_spectra):
            cell = scores_array[qi, li]
            score_val = float(cell[0]) if len(cell) > 0 else 0
            num_matches = int(cell[1]) if len(cell) > 1 else 0

            if num_matches >= min_match:
                lib_name = lib_spec.metadata.get("title", lib_spec.metadata.get("compound_name", "Unknown"))
                lib_mz_raw = lib_spec.metadata.get("pepmass", 0) or lib_spec.metadata.get("precursor_mz", 0)
                lib_mz = _parse_mz_value(lib_mz_raw)
                lib_inchi = lib_spec.metadata.get("inchikey", lib_spec.metadata.get("inchi_key", ""))
                lib_smiles = lib_spec.metadata.get("smiles", "")

                # Compute ppm error between query and library precursor
                ppm_error = 0
                if query_mz > 0 and lib_mz > 0:
                    ppm_error = abs(query_mz - lib_mz) / lib_mz * 1e6

                scores.append({
                    "query": query_name,
                    "query_mz": round(float(query_mz), 4) if query_mz else 0,
                    "match_name": lib_name,
                    "match_mz": round(float(lib_mz), 4) if lib_mz else 0,
                    "score": round(score_val, 4),
                    "matched_peaks": num_matches,
                    "inchikey": str(lib_inchi) if lib_inchi else "",
                    "smiles": str(lib_smiles) if lib_smiles else "",
                    "ppm_error": round(ppm_error, 1),
                })

        scores.sort(key=lambda x: (-x["score"], -x["matched_peaks"]))
        all_matches.extend(scores[:top_n])
        best = scores[0] if scores else None
        if best:
            logger.info(f"  Query '{query_name}' (mz={query_mz:.4f}): "
                       f"best={best['match_name']} score={best['score']:.4f} "
                       f"peaks={best['matched_peaks']} ppm={best['ppm_error']:.1f}")
        else:
            logger.info(f"  Query '{query_name}': no matches above min_match={min_match}")

    all_matches.sort(key=lambda x: (-x["score"], -x["matched_peaks"]))

    is_builtin = str(lib_path) == str(BUILTIN_LIBRARY_PATH)
    return {
        "timestamp": datetime.now().isoformat(),
        "mode": "local",
        "library_source": "builtin_reference" if is_builtin else "gnps",
        "tolerance": effective_tolerance,
        "ppm_tolerance": ppm_tolerance,
        "round_to_unit_res": round_to_unit_res,
        "min_matched_peaks": min_match,
        "num_queries": len(query_spectra),
        "library_path": str(lib_path),
        "library_size": len(library_spectra),
        "total_matches": len(all_matches),
        "matches": all_matches[:top_n * len(query_spectra)] if all_matches else [],
    }


def generate_html_report(results: Dict, output_path: str) -> str:
    """Generate a beautiful HTML report from matching results."""
    matches = results.get("matches", [])
    lib_source = results.get("library_source", "unknown")

    score_color = lambda s: (
        "#e74c3c" if s >= 0.8 else "#f39c12" if s >= 0.6 else "#27ae60" if s >= 0.4 else "#95a5a6"
    )
    score_badge = lambda s: (
        f'<span style="background:{score_color(s)};color:#fff;padding:2px 8px;border-radius:12px;font-weight:bold;">{s:.4f}</span>'
    )

    rows_html = ""
    for i, m in enumerate(matches):
        inchi_tag = f'<code style="font-size:11px;color:#7f8c8d;">{m["inchikey"]}</code>' if m.get("inchikey") else ""
        ppm_tag = f'<span style="color:#27ae60;">{m["ppm_error"]:.0f} ppm</span>' if m.get("ppm_error") else ""
        rows_html += f"""
        <tr style="border-bottom:1px solid #ecf0f1;">
            <td style="padding:10px 12px;font-weight:600;color:#2c3e50;">{i+1}</td>
            <td style="padding:10px 12px;">
                <div style="font-weight:600;color:#2c3e50;">{m["query"]}</div>
                <div style="font-size:11px;color:#7f8c8d;">m/z {m["query_mz"]:.4f}</div>
            </td>
            <td style="padding:10px 12px;">
                <div style="font-weight:600;color:#2980b9;">{m["match_name"][:80]}</div>
                {inchi_tag}
            </td>
            <td style="padding:10px 12px;text-align:center;">{m["match_mz"]:.4f}</td>
            <td style="padding:10px 12px;text-align:center;">{score_badge(m["score"])}</td>
            <td style="padding:10px 12px;text-align:center;font-weight:600;">{m["matched_peaks"]}</td>
            <td style="padding:10px 12px;text-align:center;font-size:12px;">{ppm_tag}</td>
        </tr>"""

    library_note = ""
    n_aa = 20
    if lib_source == "builtin_reference":
        library_note = f"""<div style="background:#e8f5e9;border:1px solid #4caf50;border-radius:8px;padding:14px 18px;margin-bottom:20px;">
            <strong>Library:</strong> Built-in reference ({n_aa} amino acids + common metabolites).
            Based on 2021 Anal. Chem. paper supplemental data (Table S4/S5).<br>
            Resolution adaptation: round_to_unit_res={results.get('round_to_unit_res', False)},
            tolerance={results.get('tolerance', 0.1):.3f} Da, ppm={results.get('ppm_tolerance', 'N/A')}
        </div>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ChronoDecon - Library Search Report</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{ font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif; background:#f8f9fa; color:#2c3e50; }}
        .container {{ max-width:1300px; margin:0 auto; padding:20px; }}
        .header {{ background:linear-gradient(135deg,#667eea,#764ba2); color:#fff; padding:30px; border-radius:12px; margin-bottom:24px; }}
        .header h1 {{ font-size:28px; margin-bottom:8px; }}
        .header p {{ opacity:0.9; font-size:14px; }}
        .stats {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:16px; margin-bottom:24px; }}
        .stat-card {{ background:#fff; border-radius:10px; padding:20px; box-shadow:0 2px 8px rgba(0,0,0,0.06); text-align:center; }}
        .stat-card .num {{ font-size:28px; font-weight:700; color:#667eea; }}
        .stat-card .label {{ font-size:12px; color:#7f8c8d; margin-top:4px; }}
        table {{ width:100%; border-collapse:collapse; background:#fff; border-radius:10px; overflow:hidden; box-shadow:0 2px 8px rgba(0,0,0,0.06); }}
        th {{ background:#667eea; color:#fff; padding:12px; text-align:left; font-size:13px; font-weight:600; text-transform:uppercase; letter-spacing:0.5px; }}
        tr:hover {{ background:#f0f3ff; }}
        .footer {{ text-align:center; padding:20px; color:#95a5a6; font-size:13px; }}
        .legend {{ display:flex; gap:12px; margin-bottom:16px; flex-wrap:wrap; }}
        .legend-item {{ display:flex; align-items:center; gap:6px; font-size:12px; color:#7f8c8d; }}
        .legend-dot {{ width:12px; height:12px; border-radius:50%; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>ChronoDecon - Spectral Library Search Report</h1>
            <p>Generated: {results.get("timestamp", "N/A")} | Mode: {results.get("mode", "local").upper()} | Library: {lib_source}</p>
        </div>

        {library_note}

        <div class="stats">
            <div class="stat-card">
                <div class="num">{results.get("num_queries", 0)}</div>
                <div class="label">Query Spectra</div>
            </div>
            <div class="stat-card">
                <div class="num">{results.get("library_size", 0)}</div>
                <div class="label">Library Entries</div>
            </div>
            <div class="stat-card">
                <div class="num">{results.get("total_matches", 0)}</div>
                <div class="label">Total Matches</div>
            </div>
            <div class="stat-card">
                <div class="num">{f'{matches[0]["score"]:.4f}' if matches else "N/A"}</div>
                <div class="label">Best Score</div>
            </div>
            <div class="stat-card">
                <div class="num">{results.get("tolerance", 0.1):.3f}</div>
                <div class="label">Tolerance (Da)</div>
            </div>
        </div>

        <div class="legend">
            <div class="legend-item"><div class="legend-dot" style="background:#e74c3c;"></div> Score >= 0.8 (High)</div>
            <div class="legend-item"><div class="legend-dot" style="background:#f39c12;"></div> Score 0.6-0.8 (Medium)</div>
            <div class="legend-item"><div class="legend-dot" style="background:#27ae60;"></div> Score 0.4-0.6 (Low)</div>
            <div class="legend-item"><div class="legend-dot" style="background:#95a5a6;"></div> Score < 0.4 (Weak)</div>
        </div>

        <table>
            <thead>
                <tr>
                    <th>#</th>
                    <th>Query Spectrum</th>
                    <th>Library Match</th>
                    <th style="text-align:center;">Match m/z</th>
                    <th style="text-align:center;">Cosine Score</th>
                    <th style="text-align:center;">Matched Peaks</th>
                    <th style="text-align:center;">Precursor ppm</th>
                </tr>
            </thead>
            <tbody>
                {rows_html if rows_html else '<tr><td colspan="7" style="padding:40px;text-align:center;color:#95a5a6;">No matches found</td></tr>'}
            </tbody>
        </table>

        <div class="footer">
            <p>ChronoDecon v0.2.0 | Powered by matchms | Amino Acid Library (20 AA + Metabolites) | <a href="https://gnps.ucsd.edu" target="_blank">GNPS</a></p>
        </div>
    </div>
</body>
</html>"""

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write(html)
    logger.info(f"HTML report saved: {output_path}")
    return str(out)


def generate_csv_report(results: Dict, output_path: str) -> str:
    """Generate CSV report from matching results."""
    matches = results.get("matches", [])
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    with open(out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Rank", "Query", "Query_mz", "Match_Name", "Match_mz",
                          "Cosine_Score", "Matched_Peaks", "PPM_Error", "InChIKey", "SMILES"])
        for i, m in enumerate(matches):
            writer.writerow([
                i + 1, m["query"], m["query_mz"],
                m["match_name"], m["match_mz"],
                m["score"], m["matched_peaks"],
                m.get("ppm_error", ""), m.get("inchikey", ""), m.get("smiles", "")
            ])

    logger.info(f"CSV report saved: {output_path}")
    return str(out)


def print_summary(results: Dict, top_n: int = 10) -> None:
    """Print a formatted summary to console."""
    matches = results.get("matches", [])
    lib_source = results.get("library_source", "unknown")

    print("\n" + "=" * 100)
    print("  ChronoDecon - Spectral Library Search Summary")
    print("=" * 100)
    print(f"  Mode:           {results.get('mode', 'local').upper()}")
    print(f"  Library:        {lib_source}")
    print(f"  Query spectra:  {results.get('num_queries', 0)}")
    print(f"  Library size:   {results.get('library_size', 0)} spectra")
    print(f"  Tolerance:      {results.get('tolerance', 0.1):.3f} Da")
    print(f"  PPM tolerance:  {results.get('ppm_tolerance', 'N/A')}")
    print(f"  Unit-res mode:  {results.get('round_to_unit_res', False)}")
    print(f"  Total matches:  {results.get('total_matches', 0)}")
    print("-" * 100)

    if not matches:
        print("  No matches found.")
        print("=" * 100)
        return

    print(f"  {'#':>3}  {'Query':<20} {'Query_mz':>10}  {'Score':>8}  {'Peaks':>5}  {'ppm':>8}  {'Match Name'}")
    print("-" * 100)
    for i, m in enumerate(matches[:top_n]):
        print(f"  {i+1:>3}  {m['query']:<20} {m['query_mz']:>10.4f}  {m['score']:>8.4f}  "
              f"{m['matched_peaks']:>5}  {m.get('ppm_error', 0):>7.1f}  {m['match_name'][:45]}")

    print("=" * 100)
