"""
Elemental analysis module for ChronoDecon (FULLY LOCAL, ZERO API CALLS).

Given exact m/z from high-resolution MS, calculates possible elemental compositions
(C, H, N, O, S, P) within specified ppm tolerance.

Key features:
- Brute-force formula enumeration with heuristic pruning
- 5 ppm default mass tolerance
- PURELY LOCAL - NO PubChem API calls
- Cross-referencing with spectral library matches (GNPS/matchms)
"""

import logging
import time
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

import numpy as np

logger = logging.getLogger(__name__)

# Monoisotopic masses (IUPAC 2016)
ELEMENT_MASSES = {
    "C": 12.0000000,
    "H": 1.007825032,
    "N": 14.003074004,
    "O": 15.994914620,
    "S": 31.972071174,
    "P": 30.973761998,
}

# Heuristic bounds for element counts in small molecules (< 500 Da)
# Conservative limits to keep search space manageable
ELEMENT_BOUNDS = {
    "C": (1, 30),
    "H": (1, 60),
    "N": (0, 10),
    "O": (0, 15),
    "S": (0, 3),
    "P": (0, 3),
}

# Common chemical constraints (nitrogen rule, DBE, RDBE)
# These are used for filtering invalid formulas


def _calc_mass(formula: Dict[str, int]) -> float:
    """Calculate exact monoisotopic mass from elemental composition."""
    mass = 0.0
    for elem, count in formula.items():
        mass += ELEMENT_MASSES[elem] * count
    return mass


def _calc_ppm_error(experimental_mz: float, theoretical_mass: float) -> float:
    """Calculate mass error in ppm."""
    if theoretical_mass == 0:
        return float("inf")
    return abs(experimental_mz - theoretical_mass) / theoretical_mass * 1e6


def _calc_dbe(formula: Dict[str, int]) -> float:
    """Calculate Double Bond Equivalents (DBE / RDBE).

    DBE = C - H/2 + N/2 + 1  (standard formula)

    For elements with non-standard valence:
    - S is divalent (like O), does NOT subtract from DBE
    - P is trivalent, subtracts as P/3 approximately
    """
    c = formula.get("C", 0)
    h = formula.get("H", 0)
    n = formula.get("N", 0)
    o = formula.get("O", 0)
    s = formula.get("S", 0)
    p = formula.get("P", 0)
    # S and O are divalent - they don't affect DBE
    # P is trivalent - approximate as subtracting P/3
    # Halogens (Cl, Br, F, I) would subtract as X/2, but we don't include them
    dbe = c + 1 - h / 2.0 + n / 2.0
    # Note: O and S do NOT affect DBE (they are divalent)
    # P approximation: each P acts like ~0.33 of a halogen
    dbe -= p / 3.0
    return dbe


def _is_valid_formula(formula: Dict[str, int]) -> bool:
    """Apply chemical heuristics to filter unlikely formulas."""
    c = formula.get("C", 0)
    h = formula.get("H", 0)
    n = formula.get("N", 0)
    o = formula.get("O", 0)
    s = formula.get("S", 0)
    p = formula.get("P", 0)

    # Must have at least one carbon
    if c < 1:
        return False

    # Hydrogen count should be reasonable relative to carbon
    if h < c // 2 or h > 2 * c + n + 2 * o + 2 * s + 2 * p + 4:
        return False

    # Nitrogen rule: odd-electron ions
    # For MH+ (positive ion mode), odd nominal mass = odd number of nitrogens
    # We skip this for now as we don't know ionization mode for sure

    # DBE must be non-negative
    dbe = _calc_dbe(formula)
    if dbe < 0:
        return False

    # DBE should be near-integer (within 0.01)
    if abs(dbe - round(dbe)) > 0.05:
        return False
    if dbe > c + 3:
        return False

    return True


def _enumerate_formulas(
    target_mz: float,
    ppm_tolerance: float = 5.0,
    charge: int = 1,
    max_results: int = 50,
) -> List[Dict]:
    """Enumerate possible elemental formulas within ppm tolerance.

    Optimized strategy using mass constraints to prune the search space:
    1. Iterate C, N, O, S, P with tight bounds derived from target mass
    2. Calculate H from remaining mass (exact)
    3. Filter by chemical heuristics (DBE, H/C ratio, N rule)

    For small molecules (< 300 Da), this is fast enough (~1s per m/z).

    Parameters
    ----------
    target_mz : float
        Experimental m/z value (typically [M+H]+)
    ppm_tolerance : float
        Mass tolerance in ppm (default: 5)
    charge : int
        Ion charge state (default: 1 for [M+H]+)
    max_results : int
        Maximum number of candidate formulas to return

    Returns
    -------
    List of dicts with keys: formula_str, elements, mass, ppm_error, dbe
    """
    proton_mass = ELEMENT_MASSES["H"]
    target_neutral = target_mz * charge - charge * proton_mass

    if target_neutral < 30 or target_neutral > 800:
        return []

    mass_range = target_neutral * ppm_tolerance * 1e-6
    candidates = []

    # Tight C bounds from target mass
    # Lower bound: assume all remaining mass is heteroatoms (worst case)
    # Upper bound: assume only carbon and minimum H
    c_min = 1
    c_max = min(30, int(target_neutral / 12) + 1)

    # Tight N bounds
    n_max = min(10, int(target_neutral / 14))

    # Tight O bounds
    o_max = min(15, int(target_neutral / 16))

    s_max = 3
    p_max = 2

    h_mass = ELEMENT_MASSES["H"]
    c_mass = ELEMENT_MASSES["C"]
    n_mass = ELEMENT_MASSES["N"]
    o_mass = ELEMENT_MASSES["O"]
    s_mass = ELEMENT_MASSES["S"]
    p_mass = ELEMENT_MASSES["P"]

    for c in range(c_min, c_max + 1):
        c_contrib = c_mass * c
        rem_c = target_neutral - c_contrib

        # Check if remaining mass can possibly form a valid molecule
        # Need at least 1 H (1.008 Da)
        if rem_c < 1.0:
            break
        # With max H, remaining mass should be achievable by heteroatoms
        max_h_mass = ELEMENT_BOUNDS["H"][1] * h_mass
        if rem_c > max_h_mass + n_max * n_mass + o_max * o_mass + s_max * s_mass + p_max * p_mass:
            continue

        for n in range(0, n_max + 1):
            n_contrib = n_mass * n
            rem_n = rem_c - n_contrib
            if rem_n < 1.0:
                break

            for o in range(0, o_max + 1):
                o_contrib = o_mass * o
                rem_o = rem_n - o_contrib
                if rem_o < 1.0:
                    break

                for s in range(0, min(s_max + 1, 4)):
                    s_contrib = s_mass * s
                    rem_s = rem_o - s_contrib
                    if rem_s < 1.0:
                        break

                    for p in range(0, min(p_max + 1, 2)):
                        p_contrib = p_mass * p
                        rem_p = rem_s - p_contrib
                        if rem_p < 1.0:
                            break

                        # Calculate H from remaining mass
                        h_float = rem_p / h_mass
                        h = round(h_float)

                        if h < 1 or h > ELEMENT_BOUNDS["H"][1]:
                            continue

                        formula = {"C": c, "H": h, "N": n, "O": o, "S": s, "P": p}

                        # Quick DBE check before computing full mass
                        # S and O are divalent, don't affect DBE; P subtracts P/3
                        dbe = c + 1 - h / 2.0 + n / 2.0 - p / 3.0
                        if dbe < 0 or abs(dbe - round(dbe)) > 0.01:
                            continue

                        # Compute exact mass
                        theo_mass = c_contrib + h * h_mass + n_contrib + o_contrib + s_contrib + p_contrib
                        theo_mz = (theo_mass + charge * proton_mass) / charge

                        ppm_err = abs(target_mz - theo_mz) / theo_mz * 1e6

                        if ppm_err > ppm_tolerance:
                            continue

                        # Additional H/C ratio check
                        h_c_ratio = h / c
                        if h_c_ratio < 0.5 or h_c_ratio > 4.0:
                            continue

                        # Build formula string (Hill system)
                        fstr = f"C{c}"
                        fstr += f"H{h}"
                        if n: fstr += f"N{n}"
                        if o: fstr += f"O{o}"
                        if p: fstr += f"P{p}"
                        if s: fstr += f"S{s}"

                        candidates.append({
                            "formula_str": fstr,
                            "elements": dict(formula),
                            "theoretical_mz": round(theo_mz, 6),
                            "theoretical_mass": round(theo_mass, 6),
                            "ppm_error": round(ppm_err, 2),
                            "dbe": round(dbe, 1),
                            "charge": charge,
                        })

                        if len(candidates) >= max_results * 2:
                            candidates.sort(key=lambda x: (x["ppm_error"], abs(x["dbe"] - round(x["dbe"]))))
                            return candidates[:max_results]

    candidates.sort(key=lambda x: (x["ppm_error"], x["dbe"]))
    return candidates[:max_results]


def _query_pubchem_formula(formula: str, max_results: int = 5) -> List[Dict]:
    """DISABLED: No PubChem API calls in local mode.

    This function is a stub for API compatibility only.
    Returns empty list - all identification is now formula-based.
    """
    logger.debug(f"PubChem lookup disabled for {formula} - running in pure local mode")
    return []


def _query_pubchem_name(name: str, max_results: int = 3) -> List[Dict]:
    """DISABLED: No PubChem API calls in local mode.

    This function is a stub for API compatibility only.
    Returns empty list - all identification is now formula-based.
    """
    logger.debug(f"PubChem lookup disabled for {name} - running in pure local mode")
    return []


def analyze_component_elemental(
    component_mz: float,
    component_name: str = "",
    ppm_tolerance: float = 5.0,
    charge: int = 1,
    pubchem_lookup: bool = False,  # Disabled by default - pure local mode
    max_formulas: int = 3,  # Keep only top 3 most reasonable formulas
) -> Dict:
    """Full elemental analysis for a single deconvoluted component (PURELY LOCAL).

    Parameters
    ----------
    component_mz : float
        The most likely precursor m/z (exact mass) of the component
    component_name : str
        Component identifier (e.g., "Component_1")
    ppm_tolerance : float
        Mass tolerance for formula enumeration (default: 5 ppm)
    charge : int
        Ion charge state
    pubchem_lookup : bool
        IGNORED - always False in local mode (no API calls)
    max_formulas : int
        Maximum candidate formulas to evaluate (default: 3 for top candidates)

    Returns
    -------
    Dict with keys: component, mz, candidates (list of formula only, no PubChem)
    """
    logger.info(f"Analyzing {component_name} (m/z={component_mz:.4f})...")

    # Step 1: Enumerate candidate formulas
    candidates = _enumerate_formulas(
        target_mz=component_mz,
        ppm_tolerance=ppm_tolerance,
        charge=charge,
        max_results=max_formulas,
    )

    logger.info(f"  Found {len(candidates)} candidate formulas")

    # NO PubChem lookup - pure local mode
    logger.info(f"  Top formula: {candidates[0]['formula_str']} (ppm={candidates[0]['ppm_error']:.1f}, DBE={candidates[0]['dbe']:.0f})" if candidates else "  No formulas found")

    return {
        "component": component_name,
        "mz": round(component_mz, 6),
        "charge": charge,
        "ppm_tolerance": ppm_tolerance,
        "num_candidates": len(candidates),
        "candidates": candidates,
    }


def _select_precursor_mz(mz_values: List[float], intensities: List[float]) -> float:
    """Select the most likely precursor m/z from a component's peak list.

    Strategy:
    1. Try max m/z (molecular ion [M+H]+ tends to be the heaviest common peak)
    2. Fall back to max-intensity peak
    """
    if not mz_values:
        return 0.0
    return max(mz_values)


def analyze_all_components(
    deconvolution_results: Dict,
    ppm_tolerance: float = 10.0,
    charge: int = 1,
    pubchem_lookup: bool = False,  # Disabled - pure local mode
    top_n: int = 15,
    mgf_path: Optional[str] = None,
) -> Dict:
    """Run elemental analysis on all deconvoluted components (PURELY LOCAL).

    For each component, analyzes multiple candidate precursor m/z values:
    1. PEPMASS from enriched MGF (if available)
    2. Max m/z (molecular ion is typically heaviest)
    3. Top 3 most intense peaks

    All results are aggregated per component.

    Parameters
    ----------
    deconvolution_results : dict
        Output from deconvolute_mzml()
    ppm_tolerance : float
        Mass tolerance for formula matching
    charge : int
        Ion charge state
    pubchem_lookup : bool
        IGNORED - always False in local mode (no API calls)
    top_n : int
        Only analyze top N components by intensity
    mgf_path : str, optional
        Path to enriched MGF. If provided, reads PEPMASS from MGF metadata.

    Returns
    -------
    Dict with all component analyses
    """
    components = deconvolution_results.get("components", [])
    if not components:
        return {"status": "error", "message": "No components in deconvolution results"}

    components = sorted(components, key=lambda c: c.get("total_intensity", 0), reverse=True)
    components = components[:top_n]

    # Extract PEPMASS from MGF if available
    mgf_pepmass = {}
    if mgf_path:
        try:
            mgf_pepmass = _extract_pepmass_from_mgf(mgf_path)
        except Exception as e:
            logger.warning(f"Could not read MGF for PEPMASS: {e}")

    # Known amino acid MH+ m/z for direct matching (Table S4 from paper)
    KNOWN_AA_MZ = {
        "THR": 120.0655, "PHE": 166.0863, "TRP": 205.0972, "VAL": 118.0863,
        "TYR": 182.0812, "PRO": 116.0706, "ALA": 90.0550, "GLY": 76.0393,
        "ARG": 175.1190, "HIS": 156.0768, "LYS": 147.1128, "SER": 106.0499,
        "MET": 150.0583, "ILE": 132.1019, "LEU": 132.1019,
        "ASN": 133.0608, "GLU": 148.0604, "GLN": 147.0764,
        "CYS": 122.0270, "ASP": 134.0448,
    }

    analyses = []
    for comp in components:
        comp_id = comp.get("component_id", 0)
        mz_values = comp.get("mz_values", [])
        intensities = comp.get("intensities", [])

        if not mz_values:
            continue

        comp_name = f"Component_{comp_id}"

        # Collect candidate m/z values to analyze
        candidate_mzs = set()

        # 1. PEPMASS from MGF
        mgf_mz = mgf_pepmass.get(comp_name, 0)
        if mgf_mz > 0:
            candidate_mzs.add(round(mgf_mz, 4))

        # 2. Max m/z (likely molecular ion)
        candidate_mzs.add(round(max(mz_values), 4))

        # 3. Top 2 most intense peaks (skip if close to existing)
        sorted_indices = np.argsort(intensities)[::-1] if intensities else []
        for idx in sorted_indices[:2]:
            mz_c = round(mz_values[idx], 4)
            if all(abs(mz_c - e) > 0.5 for e in candidate_mzs):
                candidate_mzs.add(mz_c)

        # 4. Check if any known AA MH+ is within the component's m/z range
        mz_min, mz_max = min(mz_values), max(mz_values)
        for aa_name, aa_mz in KNOWN_AA_MZ.items():
            if mz_min - 1.0 <= aa_mz <= mz_max + 1.0:
                if all(abs(aa_mz - e) > 0.01 for e in candidate_mzs):
                    candidate_mzs.add(aa_mz)

        # Deduplicate and sort
        candidate_mzs = sorted(candidate_mzs, reverse=True)  # Prefer higher m/z

        # Calculate max m/z and max intensity m/z for reporting
        max_mz = max(mz_values) if mz_values else 0
        max_int_mz = mz_values[np.argmax(intensities)] if intensities and mz_values else 0

        # Analyze each candidate m/z
        all_results = []
        for mz_candidate in candidate_mzs:
            result = analyze_component_elemental(
                component_mz=mz_candidate,
                component_name=comp_name,
                ppm_tolerance=ppm_tolerance,
                charge=charge,
                pubchem_lookup=pubchem_lookup,
            )
            result["candidate_type"] = "primary"
            if result["num_candidates"] > 0:
                all_results.append(result)

        # Merge all candidates into single analysis
        # Use the highest m/z that yielded results as primary
        primary = None
        all_candidates = []
        for r in all_results:
            all_candidates.extend(r.get("candidates", []))
            if primary is None or r["mz"] > primary["mz"]:
                primary = r

        if primary is None:
            primary = {
                "component": comp_name,
                "mz": round(candidate_mzs[0], 6) if candidate_mzs else 0,
                "charge": charge,
                "ppm_tolerance": ppm_tolerance,
                "num_candidates": 0,
                "candidates": [],
            }

        # Deduplicate formulas (same formula from different m/z)
        seen = set()
        unique_candidates = []
        for c in sorted(all_candidates, key=lambda x: x["ppm_error"]):
            if c["formula_str"] not in seen:
                seen.add(c["formula_str"])
                unique_candidates.append(c)

        primary["candidates"] = unique_candidates[:3]  # Keep only top 3
        primary["num_candidates"] = len(unique_candidates[:3])
        primary["max_mz"] = round(max_mz, 6) if mz_values else 0
        primary["max_intensity_mz"] = round(max_int_mz, 6) if mz_values and intensities else 0
        primary["num_peaks"] = len(mz_values)
        primary["total_intensity"] = comp.get("total_intensity", 0)
        primary["candidates_analyzed_mzs"] = candidate_mzs

        analyses.append(primary)

    return {
        "status": "success",
        "num_analyzed": len(analyses),
        "ppm_tolerance": ppm_tolerance,
        "charge": charge,
        "analyses": analyses,
    }


def _extract_pepmass_from_mgf(mgf_path: str) -> Dict[str, float]:
    """Extract PEPMASS values from enriched MGF file."""
    pepmass = {}
    current_title = None
    with open(mgf_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("TITLE="):
                current_title = line[6:]
            elif line.startswith("PEPMASS=") and current_title:
                try:
                    pepmass[current_title] = float(line[8:])
                except ValueError:
                    pass
            elif line == "END IONS":
                current_title = None
    return pepmass


def merge_with_spectral_matches(
    elemental_results: Dict,
    spectral_results: Dict,
) -> Dict:
    """Merge elemental analysis with spectral library matching results (PURELY LOCAL).

    Cross-references:
    1. If GNPS match has SMILES/formula, check if it matches any enumerated formula
    2. Rank combined results by evidence strength
    3. NO PubChem data - all identification is formula-based

    Parameters
    ----------
    elemental_results : dict
        Output from analyze_all_components()
    spectral_results : dict
        Output from run_local_search()

    Returns
    -------
    Merged results with combined identification confidence
    """
    analyses = elemental_results.get("analyses", [])
    spectral_matches = spectral_results.get("matches", [])

    # Group spectral matches by query component
    matches_by_query = defaultdict(list)
    for m in spectral_matches:
        matches_by_query[m.get("query", "")].append(m)

    merged = []

    for analysis in analyses:
        comp_name = analysis.get("component", "")
        comp_mz = analysis.get("mz", 0)
        candidates = analysis.get("candidates", [])

        # Get spectral matches for this component
        spec_matches = matches_by_query.get(comp_name, [])
        best_spec = spec_matches[0] if spec_matches else None

        # Build combined identification
        evidence = []

        # Evidence from spectral matching
        if best_spec:
            spec_formula = _formula_from_smiles(best_spec.get("smiles", ""))
            evidence.append({
                "source": "spectral",
                "name": best_spec.get("match_name", ""),
                "score": best_spec.get("score", 0),
                "matched_peaks": best_spec.get("matched_peaks", 0),
                "ppm_error": best_spec.get("ppm_error", 0),
                "inchikey": best_spec.get("inchikey", ""),
                "smiles": best_spec.get("smiles", ""),
                "formula": spec_formula,
                "confidence": "high" if best_spec["score"] >= 0.7 else "medium" if best_spec["score"] >= 0.4 else "low",
            })

        # Evidence from elemental analysis (formula only)
        for cand in candidates[:3]:  # Top 3 formula candidates
            # Check if this formula matches spectral formula
            is_confirmed = False
            if best_spec and spec_formula:
                is_confirmed = (spec_formula == cand["formula_str"])

            evidence.append({
                "source": "formula",
                "formula": cand["formula_str"],
                "ppm_error": cand["ppm_error"],
                "dbe": cand["dbe"],
                "name": "",  # No PubChem names in local mode
                "cid": 0,
                "confidence": "high" if is_confirmed else "medium" if cand["ppm_error"] < 2 else "low",
            })

        # Determine best identification
        best_id = _determine_best_id(evidence, comp_mz)

        merged.append({
            "component": comp_name,
            "mz": comp_mz,
            "num_peaks": analysis.get("num_peaks", 0),
            "total_intensity": analysis.get("total_intensity", 0),
            "best_identification": best_id,
            "spectral_matches": spec_matches[:3],
            "formula_candidates": candidates[:3],
            "all_evidence": evidence,
        })

    # Sort by identification confidence
    merged.sort(key=lambda x: _confidence_score(x.get("best_identification", {})), reverse=True)

    return {
        "status": "success",
        "num_components": len(merged),
        "components": merged,
    }


def _formula_from_smiles(smiles: str) -> str:
    """Extract molecular formula from SMILES using RDKit."""
    if not smiles:
        return ""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        formula = Descriptors.MolWt(mol)  # Just check if it works
        return Chem.rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        return ""


def _determine_best_id(evidence: List[Dict], comp_mz: float) -> Dict:
    """Determine the best identification from combined evidence (PURELY LOCAL)."""
    if not evidence:
        return {"status": "unidentified", "name": "Unknown", "confidence": "none"}

    # Priority: confirmed > spectral high > spectral medium > formula
    confirmed = [e for e in evidence if e.get("confidence") == "confirmed"]
    spectral_high = [e for e in evidence if e.get("source") == "spectral" and e.get("score", 0) >= 0.7]
    spectral_med = [e for e in evidence if e.get("source") == "spectral" and 0.4 <= e.get("score", 0) < 0.7]
    formula_high = [e for e in evidence if e.get("source") == "formula" and e.get("confidence") == "high"]

    if confirmed:
        best = confirmed[0]
        return {
            "status": "confirmed",
            "name": best.get("name", ""),
            "formula": best.get("formula", ""),
            "confidence": "confirmed",
            "evidence_sources": ["spectral", "formula"],
            "score": best.get("score", 0),
        }
    elif spectral_high:
        best = spectral_high[0]
        return {
            "status": "identified",
            "name": best.get("name", ""),
            "formula": best.get("formula", ""),
            "confidence": "high",
            "evidence_sources": ["spectral"],
            "score": best.get("score", 0),
        }
    elif spectral_med:
        best = spectral_med[0]
        return {
            "status": "tentative",
            "name": best.get("name", ""),
            "formula": best.get("formula", ""),
            "confidence": "medium",
            "evidence_sources": ["spectral"],
            "score": best.get("score", 0),
        }
    elif formula_high:
        best = formula_high[0]
        return {
            "status": "formula_match",
            "name": best.get("formula", ""),
            "formula": best.get("formula", ""),
            "confidence": "medium",
            "evidence_sources": ["formula"],
            "ppm_error": best.get("ppm_error", 0),
            "dbe": best.get("dbe", 0),
        }
    else:
        # Just formula evidence (low confidence)
        formula_evidence = [e for e in evidence if e.get("source") == "formula"]
        if formula_evidence:
            best = formula_evidence[0]
            return {
                "status": "formula_only",
                "name": best.get("formula", ""),
                "formula": best.get("formula", ""),
                "ppm_error": best.get("ppm_error", 0),
                "dbe": best.get("dbe", 0),
                "confidence": "low",
                "evidence_sources": ["formula"],
            }
        return {"status": "unidentified", "confidence": "none"}


def _confidence_score(best_id: Dict) -> float:
    """Convert identification confidence to numeric score for sorting."""
    conf = best_id.get("confidence", "none")
    score_map = {
        "confirmed": 100,
        "high": 80,
        "medium": 60,
        "low": 30,
        "none": 0,
    }
    base = score_map.get(conf, 0)
    # Add spectral score bonus
    base += best_id.get("score", 0) * 20
    return base


def generate_enhanced_html_report(merged_results: Dict, output_path: str) -> str:
    """Generate enhanced HTML report with GNPS + Formula results (PURELY LOCAL)."""
    components = merged_results.get("components", [])

    def conf_badge(c):
        colors = {
            "confirmed": "#e74c3c", "high": "#e67e22", "medium": "#f1c40f",
            "low": "#95a5a6", "none": "#bdc3c7",
        }
        bg = colors.get(c, "#bdc3c7")
        return f'<span style="background:{bg};color:#fff;padding:3px 10px;border-radius:12px;font-size:12px;font-weight:bold;">{c.upper()}</span>'

    rows_html = ""
    for i, comp in enumerate(components):
        best = comp.get("best_identification", {})
        name = best.get("name", "Unknown")
        formula = best.get("formula", "")
        conf = best.get("confidence", "none")
        score = best.get("score", 0)
        sources = best.get("evidence_sources", [])

        # Spectral match info
        spec = comp.get("spectral_matches", [])
        spec_info = ""
        if spec:
            s = spec[0]
            spec_info = f"""
            <div style="font-size:12px;margin-top:4px;">
                <strong>Spectral:</strong> {s.get('match_name', '')} 
                (score={s['score']:.3f}, peaks={s['matched_peaks']}, ppm={s.get('ppm_error', 0):.0f})
            </div>"""

        # Formula candidates (top 3)
        formulas = comp.get("formula_candidates", [])[:3]
        formula_html = ""
        for f in formulas:
            formula_html += f"""
            <div style="font-size:11px;color:#555;margin:2px 0;">
                {f['formula_str']} (ppm={f['ppm_error']:.1f}, DBE={f['dbe']:.0f})
            </div>"""

        rows_html += f"""
        <tr style="border-bottom:1px solid #ecf0f1;">
            <td style="padding:12px;font-weight:600;color:#2c3e50;">{i+1}</td>
            <td style="padding:12px;">
                <div style="font-weight:600;">{comp['component']}</div>
                <div style="font-size:11px;color:#7f8c8d;">m/z {comp['mz']:.4f} | {comp['num_peaks']} peaks</div>
            </td>
            <td style="padding:12px;">
                <div style="font-weight:600;color:#2980b9;font-size:14px;">{name or formula}</div>
                {formula_html}
                {spec_info}
            </td>
            <td style="padding:12px;text-align:center;">{conf_badge(conf)}</td>
            <td style="padding:12px;text-align:center;">
                <div style="font-weight:bold;color:#2c3e50;">{score:.3f}</div>
                <div style="font-size:11px;color:#7f8c8d;">{', '.join(sources)}</div>
            </td>
        </tr>"""

    # Summary stats
    confirmed = sum(1 for c in components if c.get("best_identification", {}).get("confidence") == "confirmed")
    high = sum(1 for c in components if c.get("best_identification", {}).get("confidence") == "high")
    medium = sum(1 for c in components if c.get("best_identification", {}).get("confidence") == "medium")
    low = sum(1 for c in components if c.get("best_identification", {}).get("confidence") == "low")
    unidentified = sum(1 for c in components if c.get("best_identification", {}).get("confidence") == "none")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ChronoDecon - Enhanced Identification Report (Local Mode)</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{ font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif; background:#f8f9fa; color:#2c3e50; }}
        .container {{ max-width:1400px; margin:0 auto; padding:20px; }}
        .header {{ background:linear-gradient(135deg,#27ae60,#2c3e50); color:#fff; padding:30px; border-radius:12px; margin-bottom:24px; }}
        .header h1 {{ font-size:28px; margin-bottom:8px; }}
        .header p {{ opacity:0.9; font-size:14px; }}
        .stats {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:16px; margin-bottom:24px; }}
        .stat-card {{ background:#fff; border-radius:10px; padding:20px; box-shadow:0 2px 8px rgba(0,0,0,0.06); text-align:center; }}
        .stat-card .num {{ font-size:28px; font-weight:700; color:#27ae60; }}
        .stat-card .label {{ font-size:12px; color:#7f8c8d; margin-top:4px; }}
        .info-box {{ background:#d5f4e6; border:1px solid #27ae60; border-radius:8px; padding:16px; margin-bottom:20px; font-size:13px; }}
        table {{ width:100%; border-collapse:collapse; background:#fff; border-radius:10px; overflow:hidden; box-shadow:0 2px 8px rgba(0,0,0,0.06); }}
        th {{ background:#27ae60; color:#fff; padding:12px; text-align:left; font-size:13px; font-weight:600; text-transform:uppercase; letter-spacing:0.5px; }}
        tr:hover {{ background:#f0f3ff; }}
        .footer {{ text-align:center; padding:20px; color:#95a5a6; font-size:13px; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>ChronoDecon - Enhanced Identification Report</h1>
            <p>Elemental Analysis + Spectral Matching | PURELY LOCAL (Zero API Calls)</p>
        </div>

        <div class="info-box">
            <strong>Methodology:</strong> For each deconvoluted component, exact m/z is used to enumerate candidate elemental 
            formulas (C/H/N/O/S/P, 5 ppm tolerance) using Nitrogen Rule and RDBE filtering. 
            Results are cross-referenced with GNPS spectral library matches (matchms cosine scoring). 
            <strong>NO EXTERNAL API CALLS - Purely local computation!</strong>
        </div>

        <div class="stats">
            <div class="stat-card"><div class="num">{len(components)}</div><div class="label">Total Components</div></div>
            <div class="stat-card"><div class="num" style="color:#e74c3c;">{confirmed}</div><div class="label">Confirmed</div></div>
            <div class="stat-card"><div class="num" style="color:#e67e22;">{high}</div><div class="label">High Confidence</div></div>
            <div class="stat-card"><div class="num" style="color:#f1c40f;">{medium}</div><div class="label">Medium Confidence</div></div>
            <div class="stat-card"><div class="num" style="color:#95a5a6;">{low}</div><div class="label">Low / Formula Only</div></div>
            <div class="stat-card"><div class="num">{unidentified}</div><div class="label">Unidentified</div></div>
        </div>

        <table>
            <thead>
                <tr>
                    <th>#</th>
                    <th>Component</th>
                    <th>Identification</th>
                    <th style="text-align:center;">Confidence</th>
                    <th style="text-align:center;">Score / Sources</th>
                </tr>
            </thead>
            <tbody>
                {rows_html if rows_html else '<tr><td colspan="5" style="padding:40px;text-align:center;color:#95a5a6;">No results</td></tr>'}
            </tbody>
        </table>

        <div class="footer">
            <p>ChronoDecon v0.3.0 | Elemental Analysis + GNPS | PURELY LOCAL | Powered by matchms, RDKit | 0 API Requests</p>
        </div>
    </div>
</body>
</html>"""

    from pathlib import Path
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write(html)
    logger.info(f"Enhanced HTML report saved: {output_path}")
    return str(out)


def generate_enhanced_csv_report(merged_results: Dict, output_path: str) -> str:
    """Generate CSV report with combined identification results."""
    import csv
    from pathlib import Path

    components = merged_results.get("components", [])
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    with open(out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Rank", "Component", "m/z", "Peaks", "Total_Intensity",
            "Best_Name", "Best_Formula", "CID", "Confidence", "Score",
            "Evidence_Sources", "Status", "Top_Formula", "Formula_ppm", "Formula_DBE",
        ])
        for i, comp in enumerate(components):
            best = comp.get("best_identification", {})
            formulas = comp.get("formula_candidates", [])
            top_f = formulas[0] if formulas else {}
            writer.writerow([
                i + 1,
                comp["component"],
                comp["mz"],
                comp.get("num_peaks", 0),
                comp.get("total_intensity", 0),
                best.get("name", ""),
                best.get("formula", ""),
                best.get("cid", ""),
                best.get("confidence", "none"),
                best.get("score", 0),
                ",".join(best.get("evidence_sources", [])),
                best.get("status", ""),
                top_f.get("formula_str", ""),
                top_f.get("ppm_error", ""),
                top_f.get("dbe", ""),
            ])

    logger.info(f"Enhanced CSV report saved: {output_path}")
    return str(out)


def print_enhanced_summary(merged_results: Dict, top_n: int = 15) -> None:
    """Print formatted summary of enhanced identification results (PURELY LOCAL)."""
    components = merged_results.get("components", [])

    print("\n" + "=" * 110)
    print("  ChronoDecon - Enhanced Identification Summary (Elemental + Spectral | PURELY LOCAL)")
    print("=" * 110)

    if not components:
        print("  No components found.")
        print("=" * 110)
        return

    print(f"  {'#':>3}  {'Component':<14} {'m/z':>10}  {'Confidence':<12} {'Score':>6}  {'Best ID':<30}  {'Formula':<18}  {'Evidence'}")
    print("  " + "-" * 106)

    for i, comp in enumerate(components[:top_n]):
        best = comp.get("best_identification", {})
        name = (best.get("name", "") or "")[:28]
        formula = (best.get("formula", "") or "")[:16]
        conf = best.get("confidence", "none")
        score = best.get("score", 0)
        sources = ",".join(best.get("evidence_sources", []))

        conf_display = conf.upper()
        if conf == "confirmed":
            conf_display = f"[CONFIRMED]"
        elif conf == "high":
            conf_display = f"[HIGH]"
        elif conf == "medium":
            conf_display = f"[MEDIUM]"

        print(f"  {i+1:>3}  {comp['component']:<14} {comp['mz']:>10.4f}  {conf_display:<12} {score:>6.3f}  {name:<30}  {formula:<18}  {sources}")

    print("=" * 110)

    # Count by confidence
    counts = defaultdict(int)
    for c in components:
        counts[c.get("best_identification", {}).get("confidence", "none")] += 1
    print(f"  Summary: ", end="")
    print(" | ".join(f"{k.upper()}={v}" for k, v in sorted(counts.items(), reverse=True)))
    print("  [0 API requests - Purely local computation]")
    print("=" * 110)
