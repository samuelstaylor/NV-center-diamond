Input files and some output files for the 2x2x2 supercell calculations. 

Analyzed results (from running `analyze.py` can be found in analyze.out)

nv_spectra plots can be found in `images`. `nv_spectra.py` contains scripts to generate all graphs along with tutorial-style comments to explain what each plot is and how it works
Here’s a version formatted for a **README.md**, ready to copy-paste. I used Markdown formatting for headings, code, and lists so it’s clear and readable:

---

# Plot Descriptions

This project generates multiple plots analyzing Huang–Rhys factors, spectral densities, and photoluminescence (PL) spectra for NV⁻ centers in diamond using PyPL. Below is a description of each plot, what it represents, and a brief explanation of how it is computed.

---

### **1. `plot_spectral_density.py`**

**What it plots:**

* Phonon spectral density ( S(\hbar \omega) ) for the system.
* Multiple curves correspond to different Gaussian broadening parameters.

**What it shows:**

* Peaks indicate phonon modes strongly coupled to the electronic transition.
* Height of each peak is proportional to the Huang–Rhys factor (HRF) for that mode.
* Different broadenings illustrate the effect of smoothing delta functions into Gaussians.

**How it works:**

* Computes ( S(\hbar \omega) = \sum_k S_k \delta(\hbar \omega - \hbar \omega_k) ).
* Uses **displacement-based HR factors**.
* Delta functions are replaced with Gaussians for visualization.

---

### **2. `plot_lineshape_function.py`**

**What it plots:**

* The **lineshape function** ( A(\hbar \omega) ) derived from the time-domain autocorrelation of the optical transition.

**What it shows:**

* Distribution of emission or absorption intensity as a function of energy.
* Shows zero-phonon line (ZPL) and phonon sidebands.
* Peaks are broadened by temperature and Gaussian parameters.

**How it works:**

* Uses HR factors (displacements or forces).
* Constructs thermal autocorrelation of the transition dipole.
* Fourier transform converts the time-domain correlation to energy-domain spectrum.

---

### **3. `plot_pl_spectrum.py`**

**What it plots:**

* Photoluminescence (PL) spectrum.

**What it shows:**

* Main ZPL peak and phonon sidebands in the emission spectrum.
* Produces an experimental-like PL spectrum using computed HR factors.

**How it works:**

* Uses the computed lineshape function ( A(\hbar \omega) ).
* Adds ZPL energy (e.g., 637 nm ≈ 1831 meV for NV⁻).
* Applies thermal occupation and additional broadening if needed.

---

### **4. `plot_pl_forces.py`**

**What it plots:**

* Comparison of PL spectra computed using **displacements vs. forces** for HR calculation.

**What it shows:**

* Solid curve: HR from atomic displacements (Δr).
* Dashed curve: HR from forces (ΔF).
* Demonstrates that both methods produce similar PL spectra.

**How it works:**

* Forces at the excited-state geometry are projected onto normal modes.
* Lineshape and PL spectrum computed from both HR dictionaries.

---

### **5. `plot_abs_spectrum.py`**

**What it plots:**

* Absorption spectral density ( S(\hbar \omega) ) from excited-state phonons.

**What it shows:**

* Peaks correspond to phonon modes contributing to **absorption**.
* Highlights differences between absorption and emission spectra.

**How it works:**

* Uses excited-state phonon frequencies and displacement-based HR factors.
* Gaussian broadening applied for visualization.

---

### **6. `plot_pl_absorption.py`**

**What it plots:**

* Comparison of **PL spectrum vs. absorption spectral density** on the same plot.

**What it shows:**

* Visual comparison of emission (PL) and absorption spectra.
* Illustrates differences in phonon contributions and energy distribution.

**How it works:**

* HR factors computed for emission (ground → excited) and absorption (excited → ground).
* Converts HR factors into spectral densities and plots both on the same energy axis.

---

### **7. `plot_td_pl.py`**

**What it plots:**

* Temperature-dependent PL spectra for multiple temperatures.

**What it shows:**

* ZPL broadening increases with temperature.
* Phonon sidebands change intensity with thermal population.
* Stacked curves illustrate evolution from cryogenic to room temperature.

**How it works:**

* Forces-based HR factors are used.
* Lineshape computed for each temperature.
* ZPL broadening adjusted per temperature.
* Curves are vertically offset for clarity in visualization.

---

Do you want me to also **add these descriptions as top-of-file docstrings in each Python script** so they appear when someone opens the file?
