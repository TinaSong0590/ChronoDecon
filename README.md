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
