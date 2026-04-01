"""
ChronoDecon Web Dashboard
Clean, modern interface inspired by slothui design
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from pathlib import Path
import tempfile
import time
from datetime import datetime
import subprocess
import os
import sys

# Import ChronoDecon package
try:
    from chrono_decon import decon
    from chrono_decon import library_search
    from chrono_decon.decon import deconvolute_mzml
    from chrono_decon.library_search import run_local_search, generate_html_report, generate_csv_report
except ImportError as e:
    st.error(f"Cannot import ChronoDecon module: {str(e)}")
    st.error(f"Python path: {sys.path}")
    st.error("Make sure you run this from the repository root where 'chrono_decon/' is located.")
    st.stop()

# Page configuration
st.set_page_config(
    page_title="ChronoDecon Dashboard",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Modern Clean CSS - Light Theme
st.markdown("""
<style>
    :root {
        --primary: #6366F1;
        --primary-dark: #4F46E5;
        --primary-light: #818CF8;
        --bg-main: #F8FAFC;
        --bg-white: #FFFFFF;
        --bg-sidebar: #FFFFFF;
        --text-primary: #1E293B;
        --text-secondary: #64748B;
        --text-muted: #94A3B8;
        --border-color: #E2E8F0;
        --success: #10B981;
        --warning: #F59E0B;
        --error: #EF4444;
    }

    .stApp {
        background: var(--bg-main);
    }

    [data-testid="stSidebar"] {
        background: var(--bg-sidebar);
        border-right: 1px solid var(--border-color);
    }

    [data-testid="stSidebar"] > div:first-child {
        padding: 20px 16px;
    }

    h1, h2, h3, h4, h5, h6 {
        color: var(--text-primary) !important;
        font-weight: 600;
    }

    /* Navigation Items */
    .nav-item {
        display: flex;
        align-items: center;
        gap: 12px;
        padding: 10px 14px;
        border-radius: 8px;
        color: var(--text-secondary);
        font-size: 14px;
        font-weight: 500;
        margin-bottom: 4px;
        cursor: pointer;
        transition: all 0.2s;
    }

    .nav-item:hover {
        background: #F1F5F9;
        color: var(--text-primary);
    }

    .nav-item.active {
        background: var(--primary);
        color: white;
    }

    .nav-badge {
        background: var(--primary-light);
        color: white;
        font-size: 11px;
        padding: 2px 8px;
        border-radius: 999px;
        margin-left: auto;
    }

    /* Upload Area - Wide and flat like reference */
    .upload-container {
        background: var(--bg-white);
        border: 2px dashed var(--border-color);
        border-radius: 12px;
        padding: 40px 24px;
        text-align: center;
        transition: all 0.3s ease;
        cursor: pointer;
    }

    .upload-container:hover {
        border-color: var(--primary);
        background: #FAFBFF;
    }

    .upload-icon {
        width: 48px;
        height: 48px;
        background: linear-gradient(135deg, var(--primary), var(--primary-light));
        border-radius: 10px;
        display: flex;
        align-items: center;
        justify-content: center;
        margin: 0 auto 16px;
        font-size: 20px;
    }

    /* File List Items */
    .file-list-item {
        background: var(--bg-white);
        border: 1px solid var(--border-color);
        border-radius: 10px;
        padding: 16px;
        margin-bottom: 12px;
    }

    .file-header {
        display: flex;
        align-items: center;
        justify-content: space-between;
        margin-bottom: 12px;
    }

    .file-info {
        display: flex;
        align-items: center;
        gap: 12px;
    }

    .file-icon {
        width: 40px;
        height: 40px;
        background: #F1F5F9;
        border-radius: 8px;
        display: flex;
        align-items: center;
        justify-content: center;
        font-size: 18px;
    }

    .file-name {
        font-weight: 600;
        color: var(--text-primary);
        font-size: 14px;
    }

    .file-meta {
        color: var(--text-muted);
        font-size: 12px;
        margin-top: 2px;
    }

    .file-status {
        font-size: 12px;
        font-weight: 500;
    }

    .status-success { color: var(--success); }
    .status-error { color: var(--error); }
    .status-progress { color: var(--primary); }

    /* Progress Bar */
    .progress-container {
        background: #E2E8F0;
        border-radius: 999px;
        height: 6px;
        overflow: hidden;
        margin-top: 8px;
    }

    .progress-bar {
        height: 100%;
        border-radius: 999px;
        transition: width 0.3s ease;
    }

    .progress-success { background: var(--success); }
    .progress-error { background: var(--error); }
    .progress-active { background: var(--primary); }

    /* Buttons */
    .stButton>button {
        background: var(--primary);
        color: white;
        border: none;
        border-radius: 8px;
        padding: 10px 20px;
        font-weight: 500;
        font-size: 14px;
        transition: all 0.2s;
    }

    .stButton>button:hover {
        background: var(--primary-dark);
        transform: translateY(-1px);
    }

    .stButton>button:disabled {
        background: #CBD5E1;
        color: #94A3B8;
    }

    /* Cards */
    .card {
        background: var(--bg-white);
        border: 1px solid var(--border-color);
        border-radius: 12px;
        padding: 20px;
    }

    .metric-card {
        background: var(--bg-white);
        border: 1px solid var(--border-color);
        border-radius: 10px;
        padding: 16px;
        text-align: center;
    }

    .metric-value {
        font-size: 28px;
        font-weight: 700;
        color: var(--primary);
    }

    .metric-label {
        font-size: 13px;
        color: var(--text-secondary);
        margin-top: 4px;
    }

    /* Section Headers */
    .section-title {
        font-size: 18px;
        font-weight: 600;
        color: var(--text-primary);
        margin-bottom: 16px;
    }

    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {
        gap: 0;
        background: #F1F5F9;
        padding: 4px;
        border-radius: 8px;
    }

    .stTabs [data-baseweb="tab"] {
        background: transparent;
        color: var(--text-secondary);
        border-radius: 6px;
        padding: 8px 16px;
        font-size: 14px;
        font-weight: 500;
    }

    .stTabs [aria-selected="true"] [data-baseweb="tab"] {
        background: white;
        color: var(--primary);
        font-weight: 600;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    }

    /* Data Table */
    .stDataFrame {
        border: 1px solid var(--border-color);
        border-radius: 10px;
    }

    /* Sidebar Logo */
    .sidebar-logo {
        display: flex;
        align-items: center;
        gap: 10px;
        padding-bottom: 20px;
        margin-bottom: 16px;
        border-bottom: 1px solid var(--border-color);
    }

    .logo-icon {
        width: 36px;
        height: 36px;
        background: linear-gradient(135deg, var(--primary), var(--primary-light));
        border-radius: 8px;
        display: flex;
        align-items: center;
        justify-content: center;
        font-size: 18px;
    }

    .logo-text {
        font-size: 18px;
        font-weight: 700;
        color: var(--text-primary);
    }

    /* User Profile */
    .user-profile {
        display: flex;
        align-items: center;
        gap: 10px;
        padding: 12px;
        background: #F8FAFC;
        border-radius: 8px;
        margin-top: auto;
    }

    .user-avatar {
        width: 32px;
        height: 32px;
        background: linear-gradient(135deg, var(--primary), var(--primary-light));
        border-radius: 50%;
        display: flex;
        align-items: center;
        justify-content: center;
        color: white;
        font-size: 12px;
        font-weight: 600;
    }

    /* Status Badges */
    .badge {
        display: inline-flex;
        align-items: center;
        padding: 4px 10px;
        border-radius: 999px;
        font-size: 12px;
        font-weight: 500;
    }

    .badge-success {
        background: rgba(16, 185, 129, 0.1);
        color: var(--success);
    }

    .badge-warning {
        background: rgba(245, 158, 11, 0.1);
        color: var(--warning);
    }

    .badge-error {
        background: rgba(239, 68, 68, 0.1);
        color: var(--error);
    }
</style>
""", unsafe_allow_html=True)

# Session state initialization
if "results" not in st.session_state:
    st.session_state.results = None
if "uploaded_file" not in st.session_state:
    st.session_state.uploaded_file = None
if "analysis_complete" not in st.session_state:
    st.session_state.analysis_complete = False
if "current_view" not in st.session_state:
    st.session_state.current_view = "upload"


# ======================= RAW to MZML Conversion Functions =======================

def convert_raw_to_mzml_thermoparser(raw_path, output_dir):
    """Convert RAW file to mzML using ThermoRawFileParser (previously successful method).

    This method was successfully used to convert sample_from_gdrive.raw on Mar 30 13:29.
    Location: /home/knan/ThermoRawFileParser/ThermoRawFileParser
    Command: ThermoRawFileParser -i input.raw -o /output/dir/ -f 2
    (where -f 2 specifies mzML format)
    """
    raw_file = Path(raw_path)
    output_path = Path(output_dir)
    mzml_file = output_path / f"{raw_file.stem}.mzML"

    # Find ThermoRawFileParser executable
    parser_paths = [
        "/home/knan/ThermoRawFileParser/ThermoRawFileParser",
        "/home/knan/ThermoRawFileParser",
        str(Path(__file__).parent.parent / "ThermoRawFileParser" / "ThermoRawFileParser"),
    ]

    parser_exe = None
    for path in parser_paths:
        if Path(path).exists():
            parser_exe = path
            break

    if not parser_exe:
        raise Exception(
            "ThermoRawFileParser not found. Please download it from:\n"
            "https://github.com/CompOmics/ThermoRawFileParser/releases\n"
            "Extract to: /home/knan/ThermoRawFileParser/"
        )

    # Make sure it's executable
    if not os.access(parser_exe, os.X_OK):
        os.chmod(parser_exe, 0o755)

    # Run ThermoRawFileParser
    # -f 2 specifies mzML format
    cmd = [
        parser_exe,
        '-i', str(raw_file),
        '-o', str(output_path),
        '-f', '2'
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        if result.returncode != 0:
            # ThermoRawFileParser might output to stderr even on success
            if not mzml_file.exists():
                raise Exception(f"ThermoRawFileParser failed: {result.stderr}")

        # Check for output file
        if mzml_file.exists():
            return str(mzml_file)
        else:
            # Try to find generated file
            generated_files = list(output_path.glob("*.mzML"))
            if generated_files:
                return str(generated_files[0])
            raise Exception("mzML file not found after conversion")

    except subprocess.TimeoutExpired:
        raise Exception("Conversion timeout (5 minutes)")
    except Exception as e:
        raise Exception(f"ThermoRawFileParser conversion failed: {str(e)}")


def convert_raw_to_mzml_docker(raw_path, output_dir):
    """Convert RAW file to mzML using Docker ProteoWizard."""
    raw_file = Path(raw_path)
    output_path = Path(output_dir)
    
    # Check if Docker is available
    try:
        result = subprocess.run(['docker', '--version'], 
                              capture_output=True, text=True, timeout=5)
        docker_available = result.returncode == 0
    except:
        docker_available = False
    
    if not docker_available:
        raise Exception("Docker is not available. Please install Docker for RAW to mzML conversion.")
    
    # Docker conversion command
    cmd = [
        'docker', 'run', '--rm',
        '-v', f'{raw_file.parent}:/data',
        'proteme/proteowizard:latest',
        'msconvert',
        f'/data/{raw_file.name}',
        '--mzML',
        '-o', '/data/'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode != 0:
            raise Exception(f"Docker conversion failed: {result.stderr}")
        
        # Check for output file
        mzml_file = output_path / f"{raw_file.stem}.mzML"
        if not mzml_file.exists():
            # Also check if output is in original directory
            mzml_file = raw_file.parent / f"{raw_file.stem}.mzML"
        
        if mzml_file.exists():
            return str(mzml_file)
        else:
            raise Exception("mzML file not found after conversion")
            
    except subprocess.TimeoutExpired:
        raise Exception("Conversion timeout (5 minutes)")
    except Exception as e:
        raise Exception(f"Conversion failed: {str(e)}")


def convert_raw_to_mzml_pyopenms(raw_path, output_dir):
    """Convert RAW file to mzML using PyOpenMS (local, no Docker required)."""
    try:
        from pyopenms import FileHandler, MzMLFile
    except ImportError:
        raise Exception("PyOpenMS is not installed. Please install it: pip install pyopenms")
    
    raw_file = Path(raw_path)
    output_path = Path(output_dir)
    mzml_file = output_path / f"{raw_file.stem}.mzML"
    
    try:
        # Try to read and convert
        file_handler = FileHandler()
        
        # Read the RAW file
        file_handler.loadExperiment(str(raw_file))
        experiment = file_handler.getExperiment()
        
        # Write as mzML
        mzml_writer = MzMLFile()
        mzml_writer.store(str(mzml_file), experiment)
        
        if mzml_file.exists():
            return str(mzml_file)
        else:
            raise Exception("mzML file was not created")
            
    except Exception as e:
        raise Exception(f"PyOpenMS conversion failed: {str(e)}")


def convert_raw_to_mzml(raw_path, output_dir):
    """Convert RAW file to mzML with fallback methods.

    Priority (based on previous success):
    1. ThermoRawFileParser (✅ Previously successful, no Docker/sudo needed)
    2. PyOpenMS (local, no Docker required)
    3. Docker ProteoWizard (requires Docker and sudo)

    Reference: sample_from_gdrive.raw was successfully converted on Mar 30 13:29
    using ThermoRawFileParser at /home/knan/ThermoRawFileParser/ThermoRawFileParser
    """
    # Method 1: Try ThermoRawFileParser (previously successful)
    try:
        st.info("🔄 Converting RAW to mzML using ThermoRawFileParser (previously successful)...")
        return convert_raw_to_mzml_thermoparser(raw_path, output_dir)
    except Exception as e:
        st.warning(f"⚠️ ThermoRawFileParser failed: {str(e)}")

    # Method 2: Try PyOpenMS (local, no Docker needed)
    try:
        st.info("🔄 Converting RAW to mzML using PyOpenMS (local)...")
        return convert_raw_to_mzml_pyopenms(raw_path, output_dir)
    except Exception as e:
        st.warning(f"⚠️ PyOpenMS conversion failed: {str(e)}")

    # Method 3: Try Docker
    try:
        st.info("🔄 Converting RAW to mzML using Docker...")
        return convert_raw_to_mzml_docker(raw_path, output_dir)
    except Exception as e:
        st.error(f"❌ Docker conversion failed: {str(e)}")

    # All methods failed
    raise Exception(
        "All conversion methods failed. Please:\n"
        "1. Use online converter: https://msviewer.nws.oregonstate.edu/\n"
        "2. Install Docker for ProteoWizard conversion\n"
        "3. Or upload an .mzML file instead"
    )


# ======================= Helper Functions =======================

def create_chronogram_plot(time_points, intensities, title="Chronogram"):
    """Create clean chronogram plot."""
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=time_points,
        y=intensities,
        mode='lines',
        fill='tozeroy',
        name='Intensity',
        line=dict(color='#6366F1', width=2),
        fillcolor='rgba(99, 102, 241, 0.1)'
    ))
    
    fig.update_layout(
        title=dict(text=title, font=dict(size=16, color='#1E293B')),
        xaxis=dict(title='Time (min)', gridcolor='#E2E8F0', showgrid=True),
        yaxis=dict(title='Intensity', gridcolor='#E2E8F0', showgrid=True),
        plot_bgcolor='#FFFFFF',
        paper_bgcolor='#FFFFFF',
        font=dict(color='#1E293B'),
        margin=dict(l=50, r=30, t=50, b=50),
        height=350,
        showlegend=False
    )
    
    return fig


def create_spectrum_plot(mz_values, intensities, title="Mass Spectrum"):
    """Create clean mass spectrum plot."""
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=mz_values,
        y=intensities,
        marker=dict(color='#6366F1', line=dict(color='#4F46E5', width=1)),
        width=1.5
    ))
    
    fig.update_layout(
        title=dict(text=title, font=dict(size=16, color='#1E293B')),
        xaxis=dict(title='m/z', gridcolor='#E2E8F0', showgrid=True),
        yaxis=dict(title='Intensity', gridcolor='#E2E8F0', showgrid=True),
        plot_bgcolor='#FFFFFF',
        paper_bgcolor='#FFFFFF',
        font=dict(color='#1E293B'),
        margin=dict(l=50, r=30, t=50, b=50),
        height=350,
        showlegend=False
    )
    
    return fig


# ======================= HTML Report Generator =======================

def generate_enhanced_html_report(decon_results, search_results, h_threshold, tolerance, output_path):
    """Generate a professional HTML report similar to the enhanced identification report."""
    matches_list = search_results.get('matches', [])
    components = decon_results.get('components', [])

    # Calculate stats
    num_components = len(components)
    num_matches = len(matches_list)
    high_conf = sum(1 for m in matches_list if m.get('score', 0) >= 0.7)
    med_conf = sum(1 for m in matches_list if 0.5 <= m.get('score', 0) < 0.7)
    low_conf = sum(1 for m in matches_list if m.get('score', 0) < 0.5)
    avg_score = sum(m.get('score', 0) for m in matches_list) / len(matches_list) if matches_list else 0

    # Build table rows
    rows_html = ""
    for i, m in enumerate(matches_list):
        score = m.get('score', 0)
        if score >= 0.7:
            conf_bg = "#27ae60"
            conf_text = "HIGH"
        elif score >= 0.5:
            conf_bg = "#f39c12"
            conf_text = "MEDIUM"
        else:
            conf_bg = "#95a5a6"
            conf_text = "LOW"

        matched_peaks = m.get('matched_peaks', 0)
        query_mz = m.get('query_mz', 0)
        match_name = m.get('match_name', 'Unknown')
        ppm_err = m.get('ppm_error', 0)

        rows_html += f"""
        <tr style="border-bottom:1px solid #ecf0f1;">
            <td style="padding:10px 12px;font-weight:600;color:#2c3e50;">{i+1}</td>
            <td style="padding:10px 12px;">
                <div style="font-weight:600;color:#2980b9;font-size:14px;">{match_name}</div>
                <div style="font-size:11px;color:#7f8c8d;">Score: {score:.4f} | Peaks: {matched_peaks}</div>
            </td>
            <td style="padding:10px 12px;text-align:center;font-weight:600;">{query_mz:.4f}</td>
            <td style="padding:10px 12px;text-align:center;">{ppm_err:.1f} ppm</td>
            <td style="padding:10px 12px;text-align:center;">{matched_peaks}</td>
            <td style="padding:10px 12px;text-align:center;">
                <span style="background:{conf_bg};color:#fff;padding:3px 10px;border-radius:12px;font-size:12px;font-weight:bold;">{conf_text}</span>
            </td>
        </tr>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ChronoDecon - Analysis Report</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{ font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif; background:#f8f9fa; color:#2c3e50; }}
        .container {{ max-width:1400px; margin:0 auto; padding:20px; }}
        .header {{ background:linear-gradient(135deg,#667eea,#764ba2); color:#fff; padding:30px; border-radius:12px; margin-bottom:24px; }}
        .header h1 {{ font-size:28px; margin-bottom:8px; }}
        .header p {{ opacity:0.9; font-size:14px; }}
        .stats {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:16px; margin-bottom:24px; }}
        .stat-card {{ background:#fff; border-radius:10px; padding:20px; box-shadow:0 2px 8px rgba(0,0,0,0.06); text-align:center; }}
        .stat-card .num {{ font-size:28px; font-weight:700; color:#667eea; }}
        .stat-card .label {{ font-size:12px; color:#7f8c8d; margin-top:4px; }}
        .info-box {{ background:#e8eaf6; border:1px solid #667eea; border-radius:8px; padding:16px; margin-bottom:20px; font-size:13px; line-height:1.6; }}
        table {{ width:100%; border-collapse:collapse; background:#fff; border-radius:10px; overflow:hidden; box-shadow:0 2px 8px rgba(0,0,0,0.06); }}
        th {{ background:#667eea; color:#fff; padding:12px; text-align:left; font-size:13px; font-weight:600; text-transform:uppercase; letter-spacing:0.5px; }}
        th.center {{ text-align:center; }}
        td {{ padding:10px 12px; font-size:13px; }}
        tr:hover {{ background:#f0f3ff; }}
        .footer {{ text-align:center; padding:20px; color:#95a5a6; font-size:13px; margin-top:24px; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>ChronoDecon - Analysis Report</h1>
            <p>Deconvolution & Spectral Library Matching | Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>

        <div class="info-box">
            <strong>Analysis Parameters:</strong>
            H-Threshold: {h_threshold} | Mass Tolerance: {tolerance} Da | Components: {num_components} | Total Matches: {num_matches}
        </div>

        <div class="stats">
            <div class="stat-card"><div class="num">{num_components}</div><div class="label">Components</div></div>
            <div class="stat-card"><div class="num" style="color:#667eea;">{num_matches}</div><div class="label">Total Matches</div></div>
            <div class="stat-card"><div class="num" style="color:#27ae60;">{high_conf}</div><div class="label">High Confidence</div></div>
            <div class="stat-card"><div class="num" style="color:#f39c12;">{med_conf}</div><div class="label">Medium Confidence</div></div>
            <div class="stat-card"><div class="num" style="color:#95a5a6;">{low_conf}</div><div class="label">Low Confidence</div></div>
            <div class="stat-card"><div class="num" style="color:#e74c3c;">{avg_score:.3f}</div><div class="label">Avg Score</div></div>
        </div>

        <table>
            <thead>
                <tr>
                    <th>#</th>
                    <th>Identification</th>
                    <th class="center">Query m/z</th>
                    <th class="center">PPM Error</th>
                    <th class="center">Matched Peaks</th>
                    <th class="center">Confidence</th>
                </tr>
            </thead>
            <tbody>
                {rows_html if rows_html else '<tr><td colspan="6" style="padding:20px;text-align:center;color:#95a5a6;">No matches found</td></tr>'}
            </tbody>
        </table>

        <div class="footer">
            <p>ChronoDecon v0.3.0 | Powered by matchms | Generated at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
    </div>
</body>
</html>"""

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    return output_path


# ======================= Sidebar =======================

with st.sidebar:
    # Logo
    st.markdown("""
    <div class="sidebar-logo">
        <div class="logo-icon">🔬</div>
        <div class="logo-text">ChronoDecon</div>
    </div>
    """, unsafe_allow_html=True)
    
    # Navigation
    nav_items = [
        ("🏠", "Home", "upload"),
        ("📊", "Results", "results"),
        ("📋", "Components", "components"),
        ("🔬", "Details", "details"),
    ]
    
    for icon, label, key in nav_items:
        is_active = st.session_state.current_view == key
        active_class = "active" if is_active else ""
        st.markdown(f"""
        <div class="nav-item {active_class}">
            <span>{icon}</span>
            <span>{label}</span>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("<div style='height: 20px;'></div>", unsafe_allow_html=True)
    
    # Quick Stats
    st.markdown("<div style='font-size: 12px; color: #94A3B8; font-weight: 600; margin-bottom: 12px;'>STATISTICS</div>", unsafe_allow_html=True)
    
    if st.session_state.analysis_complete:
        results = st.session_state.results
        num_matches = results['search'].get('total_matches', 0)
        st.markdown(f"""
        <div style='display: flex; justify-content: space-between; padding: 8px 0; border-bottom: 1px solid #F1F5F9;'>
            <span style='color: #64748B; font-size: 14px;'>Matches</span>
            <span style='color: #6366F1; font-weight: 600;'>{num_matches}</span>
        </div>
        """, unsafe_allow_html=True)
    else:
        st.markdown("""
        <div style='display: flex; justify-content: space-between; padding: 8px 0; border-bottom: 1px solid #F1F5F9;'>
            <span style='color: #64748B; font-size: 14px;'>Status</span>
            <span style='color: #94A3B8; font-size: 14px;'>Waiting</span>
        </div>
        """, unsafe_allow_html=True)
    
    # Spacer
    st.markdown("<div style='flex: 1; min-height: 200px;'></div>", unsafe_allow_html=True)
    
    # User Profile
    st.markdown("""
    <div class="user-profile">
        <div class="user-avatar">U</div>
        <div>
            <div style='font-size: 13px; font-weight: 600; color: #1E293B;'>User</div>
            <div style='font-size: 11px; color: #94A3B8;'>Researcher</div>
        </div>
    </div>
    """, unsafe_allow_html=True)


# ======================= Main Content =======================

# Header
st.markdown("""
<div style='margin-bottom: 24px;'>
    <h1 style='font-size: 24px; margin: 0; color: #1E293B;'>Upload & Analyze</h1>
    <p style='margin: 4px 0 0 0; color: #64748B; font-size: 14px;'>Upload your mzML files for deconvolution and library search</p>
</div>
""", unsafe_allow_html=True)

# Two column layout: Left (Upload & Params) | Right (Results)
left_col, right_col = st.columns([1.3, 1])

# ======================= LEFT COLUMN =======================
with left_col:
    
    # Upload Section - Wide and flat like reference
    st.markdown("""
    <div class="section-title">📁 Upload Files</div>
    """, unsafe_allow_html=True)
    
    uploaded_file = st.file_uploader(
        "",
        type=['raw', 'mzML', 'mzml'],
        label_visibility="collapsed"
    )
    
    if uploaded_file is not None:
        st.session_state.uploaded_file = uploaded_file
        file_size = uploaded_file.size / (1024 * 1024)
        file_ext = uploaded_file.name.split('.')[-1].lower()
        
        # File item with progress (like reference image)
        st.markdown(f"""
        <div class="file-list-item">
            <div class="file-header">
                <div class="file-info">
                    <div class="file-icon">📄</div>
                    <div>
                        <div class="file-name">{uploaded_file.name}</div>
                        <div class="file-meta">{file_size:.2f} MB • {file_ext.upper()}</div>
                    </div>
                </div>
                <div class="file-status status-success">✓ Ready</div>
            </div>
            <div class="progress-container">
                <div class="progress-bar progress-success" style="width: 100%;"></div>
            </div>
        </div>
        """, unsafe_allow_html=True)
    else:
        # Empty upload area
        st.markdown("""
        <div class="upload-container">
            <div class="upload-icon">☁️</div>
            <div style='color: #6366F1; font-weight: 600; font-size: 14px; margin-bottom: 4px;'>
                Click to upload or drag and drop
            </div>
            <div style='color: #94A3B8; font-size: 13px;'>
                Supported: .raw, .mzML files (RAW will be auto-converted)
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    # Parameters Section
    st.markdown("<div style='height: 24px;'></div>", unsafe_allow_html=True)
    st.markdown("""
    <div class="section-title">⚙️ Analysis Parameters</div>
    """, unsafe_allow_html=True)
    
    with st.container():
        col1, col2 = st.columns(2)
        with col1:
            h_threshold = st.slider(
                "H-Threshold",
                min_value=0.75,
                max_value=0.98,
                value=0.87,
                step=0.01
            )
            auto_threshold = st.checkbox("Auto Threshold", value=False)
        with col2:
            tolerance = st.selectbox("Mass Tolerance (Da)", [0.1, 0.2, 0.3, 0.5, 1.0], index=2)
            subtract_background = st.checkbox("Background Subtraction", value=True)
    
    with st.expander("Advanced Parameters"):
        col3, col4 = st.columns(2)
        with col3:
            min_intensity = st.number_input("Min Intensity", value=0.0, step=1000.0)
            max_isotope_diff = st.slider("Max Isotope Diff", 0.5, 5.0, 3.0, 0.1)
        with col4:
            min_matched_peaks = st.number_input("Min Matched Peaks", 1, 20, 3)
            round_to_unit_res = st.checkbox("Unit Resolution", value=True)
    
    # Start Button
    start_analysis = st.button(
        "Start Analysis",
        type="primary",
        use_container_width=True,
        disabled=uploaded_file is None or st.session_state.analysis_complete
    )

# ======================= RIGHT COLUMN =======================
with right_col:
    
    st.markdown("""
    <div class="section-title">📊 Results Preview</div>
    """, unsafe_allow_html=True)
    
    if not st.session_state.analysis_complete:
        # Empty state
        st.markdown("""
        <div class="card" style='min-height: 400px; display: flex; flex-direction: column; justify-content: center; align-items: center; text-align: center;'>
            <div style='font-size: 48px; margin-bottom: 16px;'>📈</div>
            <div style='font-size: 16px; font-weight: 600; color: #1E293B; margin-bottom: 8px;'>Waiting for Analysis</div>
            <div style='font-size: 14px; color: #64748B;'>Upload a file and click Start Analysis</div>
        </div>
        """, unsafe_allow_html=True)
    else:
        results = st.session_state.results
        decon_results = results['decon']
        search_results = results['search']
        matches_list = search_results.get('matches', [])
        
        # Quick stats cards
        num_components = len(decon_results.get('components', []))
        num_matches = search_results.get('total_matches', 0)
        best_score = max([m.get('score', 0) for m in matches_list]) if matches_list else 0
        
        col1, col2 = st.columns(2)
        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-value">{num_components}</div>
                <div class="metric-label">Components</div>
            </div>
            """, unsafe_allow_html=True)
        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-value" style='color: #EC4899;'>{num_matches}</div>
                <div class="metric-label">Matches</div>
            </div>
            """, unsafe_allow_html=True)
        
        # Top matches list
        if matches_list:
            st.markdown("<div style='margin-top: 20px;'></div>", unsafe_allow_html=True)
            st.markdown("<div style='font-size: 13px; color: #64748B; font-weight: 600; margin-bottom: 12px;'>TOP MATCHES</div>", unsafe_allow_html=True)
            
            for i, m in enumerate(matches_list[:5]):
                score = m.get('score', 0)
                status_class = "badge-success" if score >= 0.8 else "badge-warning" if score >= 0.6 else "badge-error"
                status_text = "High" if score >= 0.8 else "Medium" if score >= 0.6 else "Low"
                
                st.markdown(f"""
                <div class="file-list-item" style='padding: 12px 16px;'>
                    <div class="file-header">
                        <div class="file-info">
                            <div class="file-icon" style='width: 36px; height: 36px; font-size: 16px;'>🧪</div>
                            <div>
                                <div class="file-name">{m.get('match_name', 'Unknown')[:30]}</div>
                                <div class="file-meta">Score: {score:.3f}</div>
                            </div>
                        </div>
                        <span class="badge {status_class}">{status_text}</span>
                    </div>
                </div>
                """, unsafe_allow_html=True)
        
        # Generate & download report
        if st.button("📄 Generate Full Report", use_container_width=True, type="primary"):
            with st.spinner("Generating report..."):
                report_dir = Path(tempfile.gettempdir()) / "chronodecon_reports"
                report_dir.mkdir(parents=True, exist_ok=True)
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                report_path = report_dir / f"chronodecon_report_{timestamp}.html"
                generate_enhanced_html_report(
                    decon_results, search_results,
                    h_threshold=h_threshold,
                    tolerance=tolerance,
                    output_path=str(report_path)
                )
                st.session_state.report_path = str(report_path)

        if "report_path" in st.session_state and st.session_state.report_path:
            rp = st.session_state.report_path
            if Path(rp).exists():
                with open(rp, 'r', encoding='utf-8') as f:
                    st.download_button(
                        label="⬇️ Download HTML Report",
                        data=f.read(),
                        file_name=Path(rp).name,
                        mime="text/html",
                        use_container_width=True
                    )

# ======================= ANALYSIS EXECUTION =======================

if start_analysis and uploaded_file is not None:
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Determine file type
        file_ext = Path(uploaded_file.name).suffix.lower()
        is_raw = file_ext in ['.raw', '.RAW']
        is_mzml = file_ext in ['.mzml', '.mzML']
        
        # Save uploaded file
        saved_file = temp_path / uploaded_file.name
        with open(saved_file, 'wb') as f:
            f.write(uploaded_file.getbuffer())
        
        # Convert if RAW file
        if is_raw:
            with st.spinner("Converting RAW to mzML (this may take a few minutes)..."):
                try:
                    mzml_path = convert_raw_to_mzml(str(saved_file), temp_path)
                    st.info(f"✅ Converted {uploaded_file.name} to mzML")
                except Exception as e:
                    st.error(f"❌ Conversion failed: {str(e)}")
                    st.error("Please ensure Docker is installed and running.")
                    st.stop()
        elif is_mzml:
            mzml_path = str(saved_file)
        else:
            st.error(f"❌ Unsupported file type: {file_ext}")
            st.error("Please upload .mzML or .RAW files")
            st.stop()
        
        # Run analysis
        with st.spinner("Analyzing..."):
            try:
                # Run deconvolution
                decon_results = deconvolute_mzml(
                    mzml_path,
                    h_threshold=h_threshold,
                    min_intensity=min_intensity,
                    max_isotope_diff=max_isotope_diff,
                    subtract_background=subtract_background,
                    auto_threshold=auto_threshold
                )
                
                # Export MGF for library search
                mgf_path = temp_path / "output.mgf"
                decon.export_enriched_mgf(
                    mzml_path,
                    decon_results,
                    str(mgf_path)
                )
                
                # Run library search
                if mgf_path.exists():
                    try:
                        # Verify matchms is available before running library search
                        import matchms
                        from matchms.importing import load_from_mgf
                    except ImportError:
                        st.error("❌ 'matchms' package is not installed or not accessible.")
                        st.error("Please install it with: pip install matchms")
                        st.error(f"Python executable: {sys.executable}")
                        st.stop()

                    search_results = run_local_search(
                        str(mgf_path),
                        tolerance=tolerance,
                        min_match=min_matched_peaks,
                        round_to_unit_res=round_to_unit_res
                    )

                    # Store results
                    st.session_state.results = {
                        'decon': decon_results,
                        'search': search_results
                    }
                    st.session_state.analysis_complete = True

                    st.success("Analysis Complete!")
                    time.sleep(0.5)
                    st.rerun()

            except Exception as e:
                st.error(f"Analysis Error: {str(e)}")
                import traceback
                st.error("Detailed error:")
                st.code(traceback.format_exc())

# ======================= DETAILED RESULTS =======================

if st.session_state.analysis_complete:
    st.markdown("<div style='margin-top: 32px;'></div>", unsafe_allow_html=True)
    st.markdown("""
    <div class="section-title">📊 Detailed Results</div>
    """, unsafe_allow_html=True)
    
    results = st.session_state.results
    decon_results = results['decon']
    search_results = results['search']
    matches_list = search_results.get('matches', [])
    
    # Tabs
    tab1, tab2, tab3 = st.tabs(["Statistics", "Components", "Details"])
    
    with tab1:
        # Stats row
        num_components = len(decon_results.get('components', []))
        num_matches = search_results.get('total_matches', 0)
        best_score = max([m.get('score', 0) for m in matches_list]) if matches_list else 0
        avg_score = sum([m.get('score', 0) for m in matches_list]) / len(matches_list) if matches_list else 0
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-value">{num_components}</div>
                <div class="metric-label">Components</div>
            </div>
            """, unsafe_allow_html=True)
        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-value" style='color: #EC4899;'>{num_matches}</div>
                <div class="metric-label">Matches</div>
            </div>
            """, unsafe_allow_html=True)
        with col3:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-value" style='color: #10B981;'>{best_score:.3f}</div>
                <div class="metric-label">Best Score</div>
            </div>
            """, unsafe_allow_html=True)
        with col4:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-value" style='color: #F59E0B;'>{avg_score:.3f}</div>
                <div class="metric-label">Avg Score</div>
            </div>
            """, unsafe_allow_html=True)
        
        # Distribution chart placeholder
        st.markdown("<div style='margin-top: 24px;'></div>", unsafe_allow_html=True)
        if matches_list:
            matches_df = pd.DataFrame(matches_list)
            fig = go.Figure()
            fig.add_trace(go.Histogram(
                x=matches_df['score'],
                nbinsx=20,
                marker_color='#6366F1',
                marker_line=dict(color='#FFFFFF', width=1)
            ))
            fig.update_layout(
                title="Score Distribution",
                xaxis_title="Cosine Similarity",
                yaxis_title="Count",
                plot_bgcolor='#FFFFFF',
                paper_bgcolor='#FFFFFF',
                height=350,
                margin=dict(l=50, r=30, t=50, b=50)
            )
            st.plotly_chart(fig, use_container_width=True)
    
    with tab2:
        if matches_list:
            # Create dataframe
            display_data = []
            for i, m in enumerate(matches_list):
                display_data.append({
                    'Rank': i + 1,
                    'Compound': m.get('match_name', 'N/A'),
                    'Score': f"{m.get('score', 0):.4f}",
                    'Matched Peaks': m.get('matched_peaks', 0),
                    'Query m/z': f"{m.get('query_mz', 0):.4f}",
                    'PPM Error': f"{m.get('ppm_error', 0):.1f}"
                })
            
            df = pd.DataFrame(display_data)
            st.dataframe(df, use_container_width=True, height=450, hide_index=True)
            
            # Export buttons
            col1, col2 = st.columns(2)
            with col1:
                if st.button("Export CSV", use_container_width=True):
                    csv_path = Path(tempfile.gettempdir()) / f"results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
                    generate_csv_report(search_results, str(csv_path))
                    st.success("CSV saved!")
            with col2:
                if st.button("Export HTML Report", use_container_width=True):
                    report_dir = Path(tempfile.gettempdir()) / "chronodecon_reports"
                    report_dir.mkdir(parents=True, exist_ok=True)
                    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                    html_path = report_dir / f"chronodecon_report_{timestamp}.html"
                    generate_enhanced_html_report(
                        decon_results, search_results,
                        h_threshold=h_threshold,
                        tolerance=tolerance,
                        output_path=str(html_path)
                    )
                    st.session_state.report_path = str(html_path)
                    st.success("HTML report generated!")

            # Download report if available
            if "report_path" in st.session_state and st.session_state.report_path:
                rp = st.session_state.report_path
                if Path(rp).exists():
                    with open(rp, 'r', encoding='utf-8') as f:
                        st.download_button(
                            label="⬇️ Download Full HTML Report",
                            data=f.read(),
                            file_name=Path(rp).name,
                            mime="text/html",
                            use_container_width=True,
                            key="download_tab2"
                        )
    
    with tab3:
        if matches_list:
            # Component selector
            component_options = [f"{i+1}. {m.get('match_name', 'Unknown')[:35]}" for i, m in enumerate(matches_list[:20])]
            selected_idx = st.selectbox("Select Component", range(len(component_options)), format_func=lambda x: component_options[x])
            
            selected_match = matches_list[selected_idx]
            
            # Info cards
            col1, col2, col3 = st.columns(3)
            with col1:
                st.markdown(f"""
                <div class="metric-card">
                    <div class="metric-value">{selected_match.get('score', 0):.3f}</div>
                    <div class="metric-label">Confidence</div>
                </div>
                """, unsafe_allow_html=True)
            with col2:
                st.markdown(f"""
                <div class="metric-card">
                    <div class="metric-value" style='color: #EC4899;'>{selected_match.get('matched_peaks', 0)}</div>
                    <div class="metric-label">Matched Peaks</div>
                </div>
                """, unsafe_allow_html=True)
            with col3:
                st.markdown(f"""
                <div class="metric-card">
                    <div class="metric-value" style='color: #10B981;'>{selected_match.get('ppm_error', 0):.1f}</div>
                    <div class="metric-label">PPM Error</div>
                </div>
                """, unsafe_allow_html=True)
            
            # Charts
            st.markdown("<div style='margin-top: 20px;'></div>", unsafe_allow_html=True)
            
            col1, col2 = st.columns(2)
            with col1:
                example_mz = np.array([100, 150, 200, 250, 300, 350, 400])
                example_intensity = np.array([100, 250, 500, 800, 600, 300, 150])
                fig_spec = create_spectrum_plot(example_mz, example_intensity, "Mass Spectrum")
                st.plotly_chart(fig_spec, use_container_width=True)
            
            with col2:
                example_time = np.arange(0, 10, 0.1)
                example_chron = np.exp(-((example_time - 3) ** 2) / 2) * 1000
                fig_chron = create_chronogram_plot(example_time, example_chron, "Chromatogram")
                st.plotly_chart(fig_chron, use_container_width=True)

# Footer
st.markdown("<div style='margin-top: 40px;'></div>", unsafe_allow_html=True)
st.markdown("""
<div style='text-align: center; padding: 20px; color: #94A3B8; font-size: 12px;'>
    ChronoDecon v0.2.0 • Mass Spectrometry Analysis Platform
</div>
""", unsafe_allow_html=True)
