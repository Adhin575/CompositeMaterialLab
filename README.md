# Composite Materials Virtual Laboratory (CMVL)
### v1.0 · University of Colorado Boulder

A production-grade desktop application for composite mechanics analysis, constitutive modeling, laminate design, failure prediction, and engineering visualization — built entirely in Python.

---

## Overview

CMVL is a single-file Python desktop tool that consolidates the core analytical workflows of composite materials engineering into one dark-themed, interactive GUI. It is designed for graduate-level coursework, research, and engineering practice in aerospace, mechanical, and materials science disciplines.

The application covers the full analytical pipeline: from raw fiber and matrix properties → micromechanics → laminate stacking → ABD matrices → failure indices → optimization — with interactive plots and export at every stage.

---

## Modules

| # | Module | Description |
|---|--------|-------------|
| 1 | **Material Database** | Browse and inspect built-in fiber and matrix property cards |
| 2 | **Constitutive Matrix** | Compute [C] and [S] for 7 symmetry classes; stability check; derived properties |
| 3 | **Micromechanics** | Four models (ROM, Halpin-Tsai, Mori-Tanaka, Self-Consistent); Vf sweep plots |
| 4 | **Laminate Analysis (CLT)** | Build stacking sequences; compute ABD matrix; ply-by-ply stresses; polar stiffness |
| 5 | **Failure Analysis** | Max Stress, Tsai-Hill, Tsai-Wu, Hashin; failure envelopes in σ₁-σ₂ space |
| 6 | **Stress-Strain Simulator** | Full 3D stress state → strains, principal stresses, von Mises, strain energy |
| 7 | **Optimization** | SciPy-based fiber volume fraction optimization; specific stiffness sweep plots |
| 8 | **Sensitivity Analysis** | Tornado charts and spider plots for parameter influence on E₁ |
| 9 | **Visualization Studio** | Directional modulus surface, Mohr's circle, failure envelope, stiffness heatmap |
| 10 | **Export Center** | Export stiffness matrices to Excel (.xlsx), material database to JSON and CSV |
| 11 | **Educational Mode** | Inline theory panels covering CLT, micromechanics, failure criteria, tensor notation |
| 12 | **Engineering Calculator** | Specific modulus/strength, safety factor, thin-wall pressure vessel, plate buckling, unit converter |

---

## Material Database

**Fibers**
- Carbon Fiber IM7, Carbon Fiber AS4
- E-Glass, S-Glass
- Kevlar-49, Basalt Fiber, Boron Fiber

**Matrices**
- Epoxy 3501-6, Epoxy 8552
- Polyester, Vinyl Ester
- PEEK, PEKK, Polyimide PMR-15

Each entry stores: E₁, E₂, ν₁₂, ν₂₃, G₁₂, G₂₃, density, Xₜ, Xc, Yₜ, Yc, S₁₂.

---

## Material Symmetry Classes (Constitutive Module)

| Class | Independent Constants | Parameters Required |
|-------|-----------------------|---------------------|
| Isotropic | 2 | E, ν |
| Transversely Isotropic | 5 | E₁, E₂, ν₁₂, ν₂₃, G₁₂, G₂₃ |
| Orthotropic | 9 | E₁–E₃, ν₁₂, ν₁₃, ν₂₃, G₁₂, G₁₃, G₂₃ |
| Monoclinic | 13 | Above + coupling terms C₁₆, C₂₆, C₃₆, C₄₅ |
| Cubic | 3 | C₁₁, C₁₂, C₄₄ |
| Hexagonal | 5 | C₁₁, C₁₂, C₁₃, C₃₃, C₄₄ |
| Anisotropic | 21 | Full upper-triangle of [C] |

The stability checker reports all six eigenvalues of [C] and flags whether the matrix is positive definite (Born stability criterion).

---

## Micromechanics Models

**Rule of Mixtures** — Voigt (longitudinal, exact upper bound) + Reuss (transverse, lower bound)

**Halpin-Tsai** — Semi-empirical model with reinforcement geometry parameter ξ (ξ = 2 for E₂, ξ = 1 for G₁₂)

**Mori-Tanaka** — Inclusion-based model using Eshelby tensor; more accurate at higher Vf

**Self-Consistent** — Iterative scheme embedding fiber in an effective medium; 50-iteration convergence

All four models are computed simultaneously at a single Vf, displayed in a comparative table, and swept across Vf = 0.01–0.80 in the plot panel.

---

## Classical Lamination Theory

The laminate module implements full CLT:

```
[N]   [A  B] [ε⁰]
[M] = [B  D] [κ ]
```

- **A** (extensional): Aᵢⱼ = Σ Q̄ᵢⱼ (zₖ − zₖ₋₁)
- **B** (coupling): Bᵢⱼ = ½ Σ Q̄ᵢⱼ (zₖ² − zₖ₋₁²)
- **D** (bending): Dᵢⱼ = ⅓ Σ Q̄ᵢⱼ (zₖ³ − zₖ₋₁³)

Engineering constants (Ex, Ey, Gxy, νxy) are back-computed from [A]⁻¹. Ply-by-ply stresses are resolved in both global and local (fiber) coordinates under user-defined in-plane loads and moments.

A default `[0/90/±45]s` 8-ply quasi-isotropic laminate is pre-loaded on startup.

---

## Failure Criteria

| Criterion | Type | Tension/Compression Distinction |
|-----------|------|---------------------------------|
| Maximum Stress | Non-interactive | Yes |
| Tsai-Hill | Interactive | No (single strength value per direction) |
| Tsai-Wu | Tensor polynomial | Yes |
| Hashin | Mode-based | Yes (separate fiber/matrix modes) |

All criteria output a **Failure Index (FI)** and **Reserve Factor (RF = 1/FI)**. The failure envelope plotter draws σ₁–σ₂ boundaries for a fixed τ₁₂.

---

## Requirements

```
Python 3.9+
numpy
scipy
matplotlib
openpyxl
tkinter   (bundled with standard Python distributions)
```

Optional (imported but gracefully skipped if absent):
```
sympy
pandas
scikit-learn
python-docx
reportlab
```

Install all at once:
```bash
pip install numpy scipy matplotlib openpyxl
```

---

## Running the Application

```bash
python GUI_Advanced.py
```

Minimum window size: **1200 × 750 px**. Recommended: **1500 × 950 px** or larger.

---

## File Structure

```
GUI_Advanced.py          ← entire application (single file, ~2500 lines)
README.md                ← this file
```

All logic, GUI, physics, and data are self-contained in the single `.py` file. No external data files or configuration are required.

---

## Architecture

The application follows a lazy-instantiation page pattern:

```
CMVLApp (root window)
├── Header bar
├── Sidebar navigation (12 buttons)
├── Content frame  ← pages are built on demand when selected
│   ├── BasePage (base class)
│   │   ├── MaterialDBPage
│   │   ├── ConstitutivePage
│   │   ├── MicromechanicsPage
│   │   ├── LaminatePage
│   │   ├── FailurePage
│   │   ├── StressStrainPage
│   │   ├── OptimizationPage
│   │   ├── SensitivityPage
│   │   ├── VisualizationPage
│   │   ├── ExportPage
│   │   ├── EducationalPage
│   │   └── CalculatorPage
└── Status bar
```

Core physics functions are pure (no GUI dependencies) and sit at module level, making them independently testable and reusable:

```
compliance_from_engineering()   → 6×6 compliance matrix [S]
stiffness_isotropic()           → isotropic [C] from E, ν
build_stiffness()               → dispatch for all 7 symmetry classes
stability_check()               → eigenvalue analysis of [C]
rule_of_mixtures()              → Voigt + Reuss
halpin_tsai()                   → Halpin-Tsai model
mori_tanaka()                   → Mori-Tanaka model
self_consistent()               → self-consistent iterative model
Q_matrix()                      → reduced stiffness [Q]
Qbar_matrix()                   → transformed [Q̄] at angle θ
ABD_matrix()                    → full CLT ABD matrices
laminate_engineering_constants()→ Ex, Ey, Gxy, νxy from ABD
failure_tsai_hill()             → Tsai-Hill FI
failure_tsai_wu()               → Tsai-Wu FI
failure_max_stress()            → Max Stress FI
failure_hashin()                → Hashin FI
stress_transform()              → global → local stress transformation
```

---

## Export Formats

| Format | Contents |
|--------|----------|
| `.xlsx` | Color-coded 6×6 stiffness matrix with diagonal/off-diagonal highlighting |
| `.json` | Full fiber and matrix databases |
| `.csv` | Tabular material database (fibers + matrices) |

---

## Keyboard / Navigation

- Click any sidebar button to switch modules instantly
- All input fields accept direct keyboard entry; press **▶ COMPUTE** (or equivalent) to update
- Matplotlib toolbar (zoom, pan, save) is embedded below each plot
- The paned divider between sidebar and content is draggable

---

## Known Limitations

- Ply-by-ply stress results use a linearized stress transformation; for thick laminates, a full 3D analysis would be more accurate
- The Mori-Tanaka implementation treats fibers as isotropic inclusions (E_avg of E₁ and E₂); anisotropic Eshelby tensors are not implemented
- No nonlinear material behavior (progressive damage, plasticity) is currently implemented
- The self-consistent model converges the shear modulus only; bulk modulus uses Halpin-Tsai as a starting point

---

## Potential Extensions

- Import custom material properties from CSV or JSON
- Progressive failure analysis (ply-by-ply degradation)
- Genetic algorithm / PSO optimization for stacking sequence
- PyVista 3D visualization of elastic surfaces and laminate stackups
- PDF/DOCX report generation with embedded figures

---

## Author

**Adhindra VS** — University of Colorado Boulder  
Composite mechanics research · Aerospace engineering

---
