# NV-center-diamond
NV center in diamond calculations using [ASE](https://ase-lib.org/), [QE](https://www.quantum-espresso.org/), and [WEST](https://west-code.org/).

### Overview

This project models the negatively charged nitrogen-vacancy (NV⁻) center in diamond — a point defect formed by substituting one carbon atom with nitrogen and removing a neighboring carbon atom. The NV⁻ center exhibits unique electronic and spin properties, making it a leading system for quantum sensing and quantum information applications (see:[ https://doi.org/10.1038/s41578-021-00306-y]( https://doi.org/10.1038/s41578-021-00306-y), [https://doi.org/10.1103/PhysRevB.74.104303](https://doi.org/10.1103/PhysRevB.74.104303), [https://doi.org/10.1088/1367-2630/13/2/025025](https://doi.org/10.1088/1367-2630/13/2/025025), [https://doi.org/10.1088/1367-2630/13/2/025019](https://doi.org/10.1088/1367-2630/13/2/025019)).



The electronic ground state of NV⁻ is a spin-triplet ³A₂ state, while the lowest optically accessible excited state is another triplet, ³E. The transition between these states gives rise to the zero-phonon line (ZPL) around [1.945 eV](https://doi.org/10.1098/rspa.1976.0039).
In this workflow:

- A DFT relaxation is performed to obtain the ground-state (³A₂) geometry.

- A BFGS TDDFT optimization is then run to find the excited-state (³E) geometry.

- The energy difference between the two relaxed states yields the adiabatic excitation energy, directly related to the ZPL.

### Step 1: build the diamond supercell and then introduce NV center

- use the script in the `build_structure` directory to generate, $2\times2\times2$ (64 atoms) and $3\times3\times3$ (216 atoms) supercells composed of diamond conventional cubic unit cell
- The script uses [ASE (Atomic Simulation Environment)](https://ase-lib.org/)

xz-plane 3x3x3 supercell (visualized using [xcrysden](http://www.xcrysden.org/))
![alt text](images/xz_3x3x3.png)

angled 3x3x3 supercell
![alt text](images/angled_3x3x3.png)

### Step 2: Quantum Espresso "relaxation" calculation to get GS geometry and GS energy
- Use the provided input file `2x2x2/relax/pw_nv_diamond_relax.in` to see how it is set up.
- A standard QE DFT planewave (PW) calculation is done by submiting `job_script.sh`
- Input files referenced the parameters used by Yu Jin for his calculations in the following paper: 
Yu Jin, Victor Wen-zhe Yu, Marco Govoni, Andrew C. Xu, and Giulia Galli
Journal of Chemical Theory and Computation 2023 19 (23), 8689-8705
DOI: [10.1021/acs.jctc.3c00986](doi.org/10.1021/acs.jctc.3c00986)
  - His data can be found here: 
    - https://notebook.rcc.uchicago.edu/files/arXiv.2309.03513/Datasets/PointDefects/NV-_diamond/PBE_216_Relaxation/

### Step 3: BFGS algorithm to find ES geometry and energies

- The BFGS (Broyden–Fletcher–Goldfarb–Shanno) algorithm is used to relax the excited-state (ES) geometry.

- It iteratively updates atomic positions using forces from the total ES energy (GS energy + excitation energy).
    - [WEST](https://west-code.org/) is used to calculated the excited state energies at each iteration using LR-TDDFT.

- [WESTpy](https://west-code.org/doc/westpy/latest/) calls pw.x and wbse.x at each step to compute these quantities.

- Convergence is reached when total forces and energy changes fall below thresholds.

- The optimized ES structure is used to compute the adiabatic excitation energy (AEE)-- which we use to approximate the zero-phonon line (ZPL) by assuming that the zero point energies are similar in the electronic ground and excited states.
- 
### Step 4: ΔSCF (Delta Self-Consistent Field) Calculation

The ΔSCF method provides an alternative way to compute excitation energies using ground-state DFT.

Instead of relying on the linear-response formalism of TDDFT, ΔSCF explicitly constructs an excited-state configuration by modifying the Kohn–Sham orbital occupations—typically promoting an electron from a filled state to an empty one while conserving spin and total charge.

A new self-consistent DFT calculation is then performed on this excited configuration to allow the charge density to relax in the presence of the excitation.

The excitation energy is obtained as the total energy difference between the excited-state SCF calculation and the ground-state SCF calculation.

This method often provides accurate results for localized excitations, such as the NV⁻ center’s ³A₂ → ³E transition, and serves as a useful benchmark against TDDFT predictions.

NOTE: DEPENDENTING ON WHAT YOU WANT THERE ARE TWO TYPES OF CALCULATIONS TO RUN:
- If you want a vertical ΔSCF (excited-state energy/forces at the ground-state geometry) — do NOT move ions; keep calculation = 'scf' and do not run ionic relaxation.

- If you want the adiabatic excited-state geometry (i.e. relax the structure on the excited-state PES), then allow ions to move and use calculation = 'relax' (or vc-relax if you want cell relaxation) with the ΔSCF occupations fixed at every ionic step.

### Step 5: Analyze
- run the `2x2x2/analyze.py` script to analyze the data and print out the important values
- the output can be found in `2x2x2/analyze.out`

### Step 6: Calculate Phonon Modes (see https://miccompy.github.io/pypl/tutorial.html)

**a. Ground-state phonon modes:**

- These modes are used to compute the **spectral density**, **lineshape function**, and **photoluminescence (PL) spectrum**.
- In the `phonon` folder, create a file `gs-dft-pw.in` containing the *relaxed ground-state geometry* obtained from Step 2 (QE relaxation).
  - Note: using an unrelaxed geometry is not recommended and will generally produce imaginary phonon frequencies since the structure is not at a minimum of the ground-state potential energy surface.
- Follow the instructions in `phonon/README.md` and `NOTES.md`.
- Run **Phonopy** to generate atomic displacements for each atom along ±x, ±y, ±z directions, taking crystal symmetry into account.
- Use `header_merge.py` to add control headers to each `supercell-***.in` file.
- Submit the jobs using `sbatch run_all_jobs.sh` to perform SCF calculations on all displaced structures. This yields **forces on atoms for each displacement**.
- Collect all forces and use **Phonopy** to compute **phonon frequencies and normal modes**, constructing the ground-state phonon mesh.

---

**b. Excited-state geometry phonon analysis (vertical approximation):**

- These calculations are used to compute **PL spectra via forces**, **absorption spectra**, and **temperature-dependent PL (TD-PL)**.
- In the `phonon` folder, create `es-dft-pw.in` containing the **relaxed excited-state geometry** obtained from TDDFT or ΔSCF (Step 3/4).
- Repeat the same Phonopy workflow: generate displacements, run SCF calculations, collect forces, and construct the phonon mesh.

- **Important clarification:**
  This procedure does *not* compute true excited-state phonon modes or the excited-state Hessian.

  Instead, it evaluates **ground-state DFT forces at the excited-state equilibrium geometry**, which are then projected onto a phonon basis.

---

**Equal Mode Approximation (EMA):**

- In practical vibronic calculations, it is commonly assumed that the vibrational normal modes of the ground and excited electronic states are approximately identical.

- This assumption is called the **Equal Mode Approximation (EMA)** (also known as the *parallel-mode approximation*).

- Mathematically, EMA assumes:
  - The phonon frequencies are approximately equal:
    \( \omega_k^{(g)} \approx \omega_k^{(e)} \)
  - The normal mode eigenvectors are approximately identical:
    \( \mathbf{e}_k^{(g)} \approx \mathbf{e}_k^{(e)} \)

- Physically, this means the excited-state potential energy surface is approximated as a **rigidly shifted version of the ground-state harmonic potential**, with unchanged curvature.

- Under EMA:
  - A single set of ground-state phonon modes is used as a common vibrational basis for both electronic states.
  - Electron–phonon coupling is described entirely by the displacement between equilibrium geometries projected onto these shared modes.

- This approximation enables efficient computation of vibrationally resolved spectra (PL and absorption) without requiring explicit excited-state phonon calculations.

- Fully computing excited-state phonons via TDDFT or linear-response methods is typically computationally prohibitive for large defect supercells.

- Note: WEST can compute excited-state forces, but full excited-state phonon spectra are not currently practical for routine use.

### STEP 7: Plot results.
- run each plot script for individual plots.
- or run `nv_spectra.py` for all plots.
- See details on each plot in `2x2x2/README.md`

### Step 8: Compute Photoluminescence (PL) and Absorption Spectra (1D-CCD Method) 
- calculations in `github_repos/NV-center-diamond/2x2x2/phonon/002_nv_diamond_1d_ccd`
- using `generate.py` generate all the images of the movement from the gs geometry to the es geometry.
- calculate es and gs at each image their using `./run_all.sh`
- This step uses the ground-state and excited-state geometries to compute ZPL, PL, and absorption spectra through the 1D configurational coordinate diagram (1D-CCD) method.
- The GS and ES phonon calculations (Step 6) generate a series of displaced geometries (“images”), each representing motion along the dominant vibrational coordinate.
- For each image, forces and energies are computed; from these, the mass-weighted displacement ΔQ, phonon frequencies, and vibrational modes are extracted.
- Using GS ↔ ES energy differences along ΔQ, a 1D potential energy curve is built for each electronic state.
- The overlap between GS and ES vibrational wavefunctions yields Franck–Condon factors, which determine:
- the PL spectrum (ES → GS transitions),
- the absorption spectrum (GS → ES transitions),
- and the ZPL position.
- With these FC factors and the chosen linewidths/broadening parameters, the code constructs the final PL and absorption spectra, which can then be plotted and compared to experiment or other theoretical methods.


### STEP 9: Plot results.
- run each plot script for individual plots.
- or run `nv_spectra_1d_ccd.py` for all plots.
- See details on each plot in `2x2x2/README.md`

NOTE:
- $\Delta$SCF would not converge for 2x2x2, 3x3x3 supercells. See the attempts to get convergence in their respective directories.
