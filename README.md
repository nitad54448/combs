# Combs Indexing - Technical Guide

This document provides a technical reference for the Combs Powder Indexing software, detailing its algorithms, parameters, and methodologies.

## How to Start

To start the program you can go [here](https://nitad54448.github.io/combs/combs.html). If it's the first time it is better to read a step-by-step guide on indexing your first pattern, see the [Quick Start Guide](#quick-start-guide).

## Table of Contents

* [Technical Overview & Methodology](#technical-overview--methodology)
* [Quick Start Guide](#quick-start-guide)
* [The User Interface](#the-user-interface)
* [Peak Finding in Detail](#peak-finding-in-detail)
* [Indexing Algorithm and Search Parameters](#indexing-algorithm-and-search-parameters)
* [Evaluating Solutions (M(20))](#evaluating-solutions-m20)
* [Space Group Analysis](#space-group-analysis)
* [Advanced Topics: Enhanced Search and Sieving](#advanced-topics-enhanced-search-and-sieving)
* [Troubleshooting & FAQ](#troubleshooting--faq)
* [References](#references)

---

<section id="introduction">

## Technical Overview & Methodology

This document provides a technical reference for the Combs Powder Indexing software. It details the program's algorithms, search parameters, and methodologies for scientists and researchers familiar with powder X-ray diffraction techniques.

### Core Methodology

The goal of *ab initio* powder indexing is to determine the unit cell parameters ($a, b, c, \alpha, \beta, \gamma$) from a list of observed diffraction peak positions ($2\theta$). This program implements a system-specific, exhaustive search algorithm to solve this problem.

The core principle is to assume that a small subset of the most intense, low-angle reflections must correspond to simple crystal planes with low-integer Miller indices $(hkl)$. The program solves the indexing equations by taking the exact number of unknown parameters for a given crystal system and solving a system of linear equations using that many observed peaks.

To linearize the indexing equations, all peak positions are first converted from $2\theta$ to Q-space, where $Q = 1/d^2$. The relationship between $Q$, the Miller indices, and the reciprocal cell parameters ($A, B, C, D, E, F$) is given by the general quadratic form:
$$Q_{hkl} = Ah^2 + Bk^2 + Cl^2 + Dkl + Ehl + Fhk$$
The software solves for these reciprocal parameters and then converts them back to the real-space cell parameters for the final solution. Each trial cell is immediately refined and scored against the full peak list.

</section>

---

<section id="quick-start">

## Quick Start Guide

Follow these steps for a standard indexing routine on a single-phase powder pattern.

1.  **Load Data File:** Use the `Select Data File` button. Supported formats include `.xy`, `.xrdml`, `.ras`, etc.
2.  **Find Peaks:** On the **Peaks** tab, adjust the `Min peak (%)`, `Radius (pts)`, and `Points` sliders to accurately capture your experimental peaks.
3.  **Review & Refine Peaks:** Critically examine the peak list. Edit 2θ positions for accuracy, delete spurious peaks (noise, Kα2), and add any missed reflections using `Ctrl + Click` on the chart. A clean list of 15-20 peaks is ideal.
4.  **Set Parameters:** On the **Parameters** tab:
    * Verify the X-ray `Wavelength (Å)`.
    * Set a chemically sensible `Max Volume (Å³)` to constrain the search.
    * Set the `2θ Error (°)` appropriate for your instrument's resolution (e.g., 0.02° for synchrotron, 0.05° for lab data).
    * Select the crystal systems to search. Start with higher symmetries first.
5.  **Start Indexing:** Click `Start Indexing`.
6.  **Analyze Solutions:** On the **Solutions** tab, review the results. Sort by M(20) and click on high-scoring solutions to visually compare the calculated (blue) and observed (red) peak markers on the chart. A good solution will also have plausible space group suggestions.

</section>

---

<section id="ui">

## The User Interface

The application is divided into the Controls Panel (left) and the Results Area (right).

### Controls Panel

This panel contains all inputs and controls, organized into three primary tabs.

#### 1. Peaks Tab
* **Peak Finding Sliders:** Control the background subtraction (`Radius`), data smoothing (`Points`), and peak detection threshold (`Min peak`, logarithmic scale).
* **2θ Range Sliders:** Define the angular range for peak detection, useful for excluding noisy regions.
* **Peak Table:** A list of all detected peaks. `2θ Obs (°)` values can be manually edited for precision.

#### 2. Parameters Tab
* **Wavelength (Å):** X-ray wavelength. Critical for d-spacing calculation.
* **Max Volume (Å³):** The upper limit for acceptable unit cell volume. A key parameter for constraining the search.
* **2θ Error (°):** Tolerance for matching calculated to observed peaks.
* **Impurity Peaks:** Number of allowed un-indexed peaks among the first 20 when calculating M(20).
* **Refine Zero-Point Error:** Enables simultaneous refinement of a zero-point correction with the lattice parameters. Highly recommended.
* **Crystal Systems to Search:** Checkboxes for Cubic, Tetragonal, Hexagonal, Orthorhombic, and Monoclinic systems.

#### 3. Solutions Tab
* **Solutions Table:** Lists all valid solutions found, with their system, parameters, volume, and M(20) score. Click a row to select it for visual analysis on the chart.

### Results Area
* **Chart:** Displays the diffraction pattern, observed peaks (red ticks), and calculated peaks for the selected solution (blue ticks). A good fit shows excellent alignment between red and blue ticks.
* **Chart Interaction:** Zoom (mouse wheel), Pan (click-drag), Reset Zoom (right-click), Add Peak (`Ctrl + Click`).

</section>

---

<section id="peak-finding">

## Peak Finding in Detail

Accurate peak positions are the most critical input for successful indexing. The program uses a multi-step process to identify peaks from raw data.

#### The Algorithm Steps
1.  **Background Subtraction:** A "rolling ball" algorithm estimates and subtracts the background signal. The `Radius` slider controls the size of the virtual ball.
2.  **Data Smoothing:** A Savitzky-Golay filter is applied to the background-subtracted data to reduce noise while preserving peak shape. The `Points` slider controls the smoothing window size.
3.  **Peak Detection:** The algorithm identifies local maxima in the smoothed data that are above the `Min peak (%)` threshold.
4.  **Position Refinement:** A parabolic fit to the top three points of each detected maximum is used to calculate a sub-pixel-accuracy peak position.

#### Practical Advice
* Start with default values and visually inspect the results.
* If weak but clear peaks are missed, lower the `Min peak (%)`. If picking up noise, increase it.
* For patterns with a broad amorphous background, increase the `Radius`.
* If data is very noisy, increase smoothing `Points`, but avoid over-smoothing, which can merge or shift peaks.
* **Always manually curate the final peak list.** Remove artifacts, Kα2 shoulders, and known impurity peaks.

</section>

---

<section id="indexing-method">

## Indexing Algorithm and Search Parameters

The program employs a dedicated search routine for each crystal system. This is an "exhaustive" or "brute-force" trial method that iterates through combinations of low-angle peaks and low-integer Miller indices. The number of peaks required to generate a trial cell depends on the number of unknown lattice parameters.

### System-by-System Search Logic

The search algorithm's goal is to solve a system of linear equations of the form $Q_{obs} = \sum P_i \cdot H_i$, where $Q_{obs}$ are the $1/d^2$ values from the observed peaks, $H_i$ are terms derived from the trial Miller indices (e.g., $h^2, k^2, l^2$), and $P_i$ are the reciprocal lattice parameters (e.g., $1/a^2, 1/b^2, 1/c^2$) we want to find.

* **Cubic (1 parameter, $A=1/a^2$):**
    Solves a 1x1 system:
    $$Q_{obs, 1} = (h_1^2 + k_1^2 + l_1^2) \cdot A$$
    The program iterates through the first 10 observed peaks, assigning each one a single trial $(hkl)$ vector (where $h,k,l$ are integers up to 8) to find a trial $a$-value.

* **Tetragonal & Hexagonal (2 parameters, $P_1, P_2$):**
    Solves a 2x2 system of equations using pairs of peaks from the first 10:
    $$Q_{obs, 1} = H_{1,a} \cdot P_1 + H_{1,c} \cdot P_2$$
    $$Q_{obs, 2} = H_{2,a} \cdot P_1 + H_{2,c} \cdot P_2$$
    It assigns pairs of trial $(hkl)$ vectors (where $h,k,l$ are integers up to 5) to solve for $P_1$ and $P_2$, which are then converted to $a$ and $c$.

* **Orthorhombic (3 parameters, $A=1/a^2, B=1/b^2, C=1/c^2$):**
    Solves a 3x3 system using triplets of peaks from the first 10:
    $$Q_{obs, 1} = h_1^2 A + k_1^2 B + l_1^2 C$$
    $$Q_{obs, 2} = h_2^2 A + k_2^2 B + l_2^2 C$$
    $$Q_{obs, 3} = h_3^2 A + k_3^2 B + l_3^2 C$$
    It assigns *triplets* of trial $(hkl)$ vectors from a pre-calculated basis set of 30 reflections to find $a, b, c$.

* **Monoclinic (4 parameters, $A, B, C, D$):**
    Solves a 4x4 system for the four reciprocal parameters ($A, B, C, D$) in the equation $Q_{hkl} = A h^2 + B k^2 + C l^2 + D hl$.
    It iterates through *quadruplets* of observed peaks (from the first 20), assigning *quadruplets* of trial $(hkl)$ vectors from a dynamically generated, curated basis set of 50 common, low-index reflections.

### Search Depth and HKL Limits

To ensure performance, the search depth and the range of Miller indices tested are internally limited. These hard-coded parameters represent a balance between exhaustive searching and computational time.

| Crystal System | Peak Subset Used | HKL Index Limits |
| :--- | :--- | :--- |
| Cubic | First 10 peaks | `h, k, l` up to 8 |
| Tetragonal | First 10 peaks | `h, k, l` up to 5 |
| Hexagonal | First 10 peaks | `h, k, l` up to 5 |
| Orthorhombic | First 10 peaks | Uses a basis set of the first 30 plausible low-index reflections. |
| Monoclinic | First 20 peaks | Uses a curated basis set of 50 common, low-index reflections. |

</section>

---

<section id="evaluating-solutions">

## Evaluating Solutions (M(20))

The indexing search often produces multiple candidate solutions. Distinguishing the correct one requires a reliable figure of merit and careful refinement.

### The de Wolff Figure of Merit: M(20)

Solutions are ranked by the **de Wolff Figure of Merit, M(20)**, which assesses the fit's accuracy and completeness. It is calculated using the first 20 observed reflections. A high M(20) value is a strong indicator of a correct solution.

| M(20) Value | Interpretation |
| :--- | :--- |
| > 20 | The solution is very likely to be correct. |
| > 10 | The solution is likely correct, especially if the volume is chemically sensible. |
| 5 - 10 | The solution is plausible and deserves further investigation. |
| < 5 | The solution is likely spurious and should be viewed with skepticism. |

### Least-Squares and Zero-Point Refinement

Every promising trial cell is refined using linear least-squares to minimize the difference between $Q_{obs}$ and $Q_{calc}$.

When the `Refine Zero-Point Error` option is enabled, the program uses a robust two-step approach:

1.  **Baseline Refinement:** It first performs a stable refinement of **only the lattice parameters**. This provides a reliable fallback solution.
2.  **Full Refinement:** It then attempts a more complex refinement solving for both the lattice parameters and the zero-point error simultaneously. This is done by adding a term to the least-squares fit that models the $2\theta$ shift.

If the full refinement is unstable or produces an invalid cell, the program gracefully falls back to the reliable baseline solution, ensuring a valid result is never discarded due to an unstable higher-order fit.

</section>

---

<section id="space-group-analysis">

## Space Group Analysis

After a high-scoring unit cell is found, the program provides an automated analysis to suggest the most probable space groups. This feature serves as a powerful guide for subsequent structure solution or Rietveld refinement.

The analysis relies on a comprehensive internal database of reflection conditions (`space_groups_complete.json`) derived from crystallographic libraries like **Pymatgen**.

### Methodology

The analysis is a systematic process of elimination based on observed systematic absences:

1.  **Index Observed Peaks:** The program takes the refined unit cell and generates a complete theoretical reflection list. It then indexes *all* high-quality observed peaks from the user's peak list, not just the first 20.
2.  **Build High-Confidence Set:** To avoid errors from peak overlap, the algorithm filters this list to find **unambiguous** reflections—observed peaks that correspond to a single, non-overlapping theoretical $(hkl)$ line.
3.  **Determine Centering:** This high-confidence set of *observed* $(hkl)$s is tested against the rules for all lattice centerings (P, I, F, A, B, C, R). The program determines the most likely centering by finding the one with the **minimum number of violations**. A violation is defined as an *observed reflection* that *should be absent* (e.g., observing a (100) reflection for an I-centered lattice).
4.  **Rank Space Groups:** The program filters its database for *all* space groups matching the solution's crystal system (e.g., 'monoclinic') and the most plausible centering(s).
5.  **Count Extinction Violations:** For each candidate space group, it checks the high-confidence observed reflections against that group's specific extinction rules for glide planes and screw axes (e.g., $0k0: k=2n$). Again, it only counts violations where an *observed peak* breaks a rule.

### How to Interpret the Results

The program ranks space groups by their number of violations. This ranking is displayed in the PDF report.

* **0 Violations:** This is the ideal result. It means that *no observed reflection* in the high-confidence set violated any of the space group's extinction rules. These are the strongest candidates.
* **1-2 Violations:** These are still plausible. A single violation could be caused by an experimental artifact, a mis-indexed ambiguous peak, or a weak, theoretically-forbidden reflection appearing due to experimental factors.
* **Many Violations:** These space groups are highly unlikely to be correct.

> **Important:** This analysis is a suggestion, not a definitive assignment. It is based *only* on the observed extinctions in a powder pattern, which can be incomplete. The final determination must be made with chemical knowledge and confirmed via structure refinement.

</section>

---

<section id="advanced-topics">

## Advanced Topics: Enhanced Search and Sieving

To improve the success rate, the program integrates several advanced search strategies after the initial indexing routine.

### 1. Matrix-based Cell Transformations
Promising solutions are transformed using crystallographic matrices to test for related primitive or higher-symmetry cells. For example, a body-centered cell ($I$) is tested for a valid primitive ($P$) equivalent.

### 2. HKL Divisor Analysis
The program examines the list of indexed Miller indices. If all indices in one direction (e.g., all $h$ values) share a common divisor, it tests a new cell with the corresponding axis halved.

### 3. Orthorhombic to Hexagonal Check
A hexagonal lattice can be indexed as a C-centered orthorhombic cell where $b/a \approx \sqrt{3}$. The program specifically checks all orthorhombic solutions for this condition and generates the equivalent hexagonal cell for evaluation.

### 4. Final Sieving
After all searches, a final sieving process is applied. If two solutions have very similar volumes, the one with higher symmetry and a better M(20) score is preferred, and the redundant solution is discarded.

</section>

---

<section id="troubleshooting">

## Troubleshooting & FAQ

#### Why were no solutions found?
* **Poor Peak List:** The most common cause. Ensure the peak list is accurate, with impurity peaks removed and positions refined. A minimum of 15-20 clean peaks is recommended.
* **Incorrect Parameters:** Verify the `Wavelength` and ensure the `Max Volume` is chemically reasonable and sufficiently large.
* **High Zero-Point Error:** A significant instrument misalignment may require enabling the `Refine Zero-Point Error` option.
* **Sample is a Mixture:** Indexing requires a peak list from a single crystalline phase.

#### Why is M(20) low for a visually good fit?
* **Incorrect Cell (Sub-multiple/Super-multiple):** The found cell may be a multiple or sub-multiple of the true cell. It indexes some lines but is penalized by M(20) for being too "empty" or "dense" with calculated reflections.
* **High Error / Low Resolution:** The `2θ Error` might be set too tight for broad peaks. Try slightly increasing this value.
* **Spurious Solution:** A random solution can sometimes fit a few lines by chance. Visual inspection is key; if major observed peaks are missed by the calculated pattern, the solution is incorrect regardless of the M(20) value.

</section>

---

<section id="references">

## References

This program was developed by Nita Dragoe at Université Paris-Saclay (2024-2025) as a successor to a prior software, Powder 4, originally created by the same author in 1999-2000. For a deeper understanding of the methodology, consulting the original scientific papers is highly recommended.

1.  **The M(20) Figure of Merit:**<br>
    de Wolff, P. M. (1968). "A Simplified Criterion for the Reliability of a Powder Pattern Indexing." *Journal of Applied Crystallography*, **1**, 108-113.

2.  **General Powder Diffraction Text:**<br>
    Klug, H. P. & Alexander, L. E. (1974). *X-Ray Diffraction Procedures for Polycrystalline and Amorphous Materials*, 2nd ed. New York: Wiley-Interscience.

3.  **Alternative Indexing Methods:**<br>
    Ito, T. (1949). "A General Powder X-ray Photography."*Nature*, 164, 755-756.<br>
    Werner, P.-E., Eriksson, L., & Westdahl, M. (1985). "TREOR, a semi-exhaustive trial-and-error powder indexing program for all symmetries." *Journal of Applied Crystallography*, **18**, 367-370.<br>
    Visser, J. W. (1969). "A fully automatic program for finding the unit cell from powder data." *Journal of Applied Crystallography*, **2**, 89-95.<br>
    Le Bail, A. (2004). "Monte Carlo Indexing with McMaille". *Powder Diffraction*, **19(3)**, 249-254.<br>
    Boultif, A. & Louër, D. (2004). "Powder pattern indexing with the dichotomy method." *Journal of Applied Crystallography*, **37**, 724-731.

---
*Help Guide updated on 20 Oct 2025 by an AI assistant.*
</section>
