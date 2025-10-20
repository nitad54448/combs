# Combs – Powder Indexing Software

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)]()
[![Status](https://img.shields.io/badge/status-active-success.svg)]()
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](CONTRIBUTING.md)
[![Made with JavaScript](https://img.shields.io/badge/Made%20with-JavaScript-yellow.svg)]()

**Combs** is a modern, browser-based program for **ab initio powder diffraction indexing**.  
It determines crystal unit cell parameters directly from experimental powder X-ray diffraction (PXRD) data.  

Developed by **Nita Dragoe** at *Université Paris-Saclay (2024–2025)*, Combs builds on the author’s earlier software *Powder 4 (1999–2000)* and implements an efficient, exhaustive search algorithm tailored for scientists and researchers working in crystallography.

---

## ✨ Features
- 📂 Import standard diffraction data files (`.xy`, `.xrdml`, `.ras`, …)
- 🔍 Automated **peak detection** with adjustable background subtraction & smoothing
- ✏️ Manual peak list editing with interactive UI
- 📐 Multi-system **indexing algorithms** (Cubic, Tetragonal, Hexagonal, Orthorhombic, Monoclinic)
- 📊 **M(20) scoring** (de Wolff Figure of Merit) for ranking solutions
- 🧭 Automated **space group analysis** using systematic absences
- ⚙️ Advanced refinement: zero-point correction, least-squares fit
- 🖼️ Interactive visualization of observed vs calculated diffraction peaks
- 📝 Exportable reports (PDF with space group suggestions)

---

## 🚀 Quick Start

1. **Load Data File**  
   Import your PXRD dataset with `Select Data File`.

2. **Find Peaks**  
   Adjust peak-finding controls in the **Peaks tab**:
   - `Min peak (%)` – detection threshold  
   - `Radius (pts)` – background subtraction  
   - `Points` – smoothing window  

3. **Refine Peak List**  
   - Edit 2θ positions for precision  
   - Delete spurious peaks (noise, Kα₂, impurities)  
   - Add missing peaks via `Ctrl + Click`  

4. **Set Parameters** (in the **Parameters tab**):  
   - Confirm `Wavelength (Å)`  
   - Set a reasonable `Max Volume (Å³)`  
   - Adjust `2θ Error (°)` based on resolution  
   - Select crystal systems to search  

5. **Start Indexing**  
   Click `Start Indexing`.

6. **Analyze Solutions**  
   - Use the **Solutions tab**  
   - Sort by M(20)  
   - Compare calculated (blue) vs observed (red) peaks  
   - Review suggested space groups  

---

## 🖥️ User Interface

**Controls Panel (left)**  
- **Peaks Tab:** peak finding controls & editable peak table  
- **Parameters Tab:** wavelength, max volume, 2θ error, zero-point, symmetry checkboxes  
- **Solutions Tab:** sortable solution table  

**Results Area (right)**  
- Interactive diffraction chart (zoom, pan, reset, add peaks)  
- Overlay of observed vs calculated peak positions  

---

## 🔬 Algorithms & Methodology

- Converts observed peaks from 2θ to reciprocal space: `Q = 1/d²`  
- Solves indexing equations in quadratic form:  
  ```
  Q(hkl) = A·h² + B·k² + C·l² + D·kl + E·hl + F·hk
  ```
- Exhaustive search per crystal system (Cubic, Tetragonal, Hexagonal, Orthorhombic, Monoclinic)  
- Each trial cell refined & scored against peak list  

See [Technical Guide](https://nitad54448.github.io/combs/combs_help.html) for detailed methodology.

---

## 📊 Evaluating Solutions

Solutions ranked with **de Wolff Figure of Merit (M(20))**:

| M(20) | Interpretation |
|-------|----------------|
| > 20  | Very likely correct |
| > 10  | Likely correct if volume makes sense |
| 5–10  | Plausible, worth further study |
| < 5   | Probably spurious |

Refinements include **least-squares fit** and optional **zero-point error correction**.

---

## 🧭 Space Group Analysis

- Generates reflection list from refined cell  
- Tests observed peaks against extinction conditions  
- Determines probable lattice centering (P, I, F, A, B, C, R)  
- Filters possible space groups by systematic absences  
- Ranks by number of extinction violations  

⚠️ Note: Results are **suggestions**, final assignment requires chemical reasoning & refinement.

---

## ⚙️ Advanced Features
- Matrix-based cell transformations (primitive ↔ centered ↔ higher symmetry)  
- HKL divisor analysis (test halved axes)  
- Orthorhombic ↔ Hexagonal checks  
- Final sieving of similar-volume cells  

---

## 🛠️ Troubleshooting

**No solutions found?**
- Poor or noisy peak list  
- Incorrect wavelength / max volume  
- High zero-point error (enable correction)  
- Mixture data (single phase required)  

**Low M(20) but visually good fit?**
- Wrong cell multiple (sub/super-cell)  
- 2θ error tolerance too tight  
- Random/spurious solution  

---

## 📚 References

- de Wolff, P. M. (1968). *J. Appl. Cryst.*, **1**, 108–113  
- Klug, H. P. & Alexander, L. E. (1974). *X-Ray Diffraction Procedures*, 2nd ed. Wiley  
- Ito, T. (1949). *Nature*, 164, 755–756  
- Werner, P.-E. et al. (1985). *J. Appl. Cryst.*, 18, 367–370  
- Visser, J. W. (1969). *J. Appl. Cryst.*, 2, 89–95  
- Le Bail, A. (2004). *Powder Diffraction*, 19(3), 249–254  
- Boultif, A. & Louër, D. (2004). *J. Appl. Cryst.*, 37, 724–731  

---

## 📥 Installation & Usage

<a href="https://nitad54448.github.io/combs/combs.html" target="_blank">Open the main program in your browser</a>

> 💡 No compilation required — Combs is pure JavaScript/HTML and runs locally in any modern browser.

---

## 📜 License

This project is licensed under the [MIT License](LICENSE).

---

## 🧾 Citation

If you use Combs in academic work, please cite:

> Nita Dragoe (2024–2025). **Combs: A Browser-Based Powder Indexing Program**. Université Paris-Saclay.  

---

## 🧑‍💻 Author

- **Nita Dragoe** – Université Paris-Saclay (2024–2025)  
- Successor to *Powder 4 (1999–2000)*  

Help guide last updated by an AI assitant: **20 Oct 2025**
