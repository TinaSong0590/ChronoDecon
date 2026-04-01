![ChronoDecon Banner](https://github.com/TinaSong0590/ChronoDecon/raw/main/chrono_decon/banner.png)

# ChronoDecon

**AI-Powered Single-Component Reconstruction for High-Resolution Mass Spectrometry**  
**高分辨质谱单组分重建工具**

An open-source tool that implements the modified cross-correlation algorithm from the 2021 *Analytical Chemistry* paper (You et al.) to automatically extract clean single-component spectra from complex co-eluting LC-MS data.

---

## 🌐 Web Dashboard (Recommended)

No programming required — perfect for mass spectrometrists and lab researchers.

**How to run locally:**

```bash
git clone https://github.com/TinaSong0590/ChronoDecon.git
cd ChronoDecon
pip install -e .
streamlit run app.py
Features:

Drag-and-drop upload of .mzML files
One-click full analysis: deconvolution + GNPS library search + elemental analysis
Beautiful interactive visualization (summary, component list, detailed spectra)


✨ Key Features

Single-component deconvolution using time-domain morphology and homologue factor (H-factor)
Local GNPS spectral matching (zero external API calls)
Pure local elemental composition analysis (RDKit-based)
MCP Tool — can be directly called by AI agents (Claude, LangGraph, etc.)
Streamlit Web Dashboard — user-friendly graphical interface (v0.2)


🚀 Quick Start
1. Using the Web Dashboard (Recommended)
Bashstreamlit run app.py
2. Using Command Line
Bashchrono-decon library-search --input your_sample.mzML

📊 Example Results
Tested on real amino acid mixture data:

Successfully identified 15 major components
Main identifications include L-Methionine, L-Proline, L-Tryptophan, etc.
Fully local execution, no external API required


Installation
Bashgit clone https://github.com/TinaSong0590/ChronoDecon.git
cd ChronoDecon
pip install -e .

Roadmap

v0.3 — Direct support for Thermo .raw files with automatic conversion
v0.4 — Batch processing for multiple files
v0.5 — Docker one-click deployment + VS Code plugin


Citation
If you use ChronoDecon in your research, please cite the original paper:
You et al. (2021). Unsupervised Reconstruction of Analyte-Specific Mass Spectra Based on Time-Domain Morphology with a Modified Cross-Correlation Approach. Analytical Chemistry, 93(12), 5009–5014.

License: MIT License
GitHub: TinaSong0590/ChronoDecon
欢迎质谱研究者、AI + 科学计算爱好者一起参与贡献！
Any feedback or collaboration is highly appreciated!
text
