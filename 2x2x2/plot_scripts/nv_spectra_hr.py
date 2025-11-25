# -- MASTER PLOTTING FILE -- BASED OFF OF PYPL TUTORIAL 1 Huang-Rhys Theory BUT FULLY COMMENTED
# cd to:        NV-center-diamond/2x2x2/
# use command:  python plot_scripts/nv_spectra_hr.py
# ---------------------------------------------------------------------------
# Tutorial-style, fully commented version of the PL (photoluminescence)
# processing script using PyPL's Huang-Rhys solver.
#
# WHAT THIS SCRIPT DOES (high-level):
#  1. Load phonon data (from phonopy HDF5) and atomic geometries (QE xml).
#  2. Compute Huang–Rhys factors (HRFs) that quantify electron-phonon
#     coupling per normal mode for NV⁻ transitions.
#  3. Build spectral densities S(ħω), phonon-broadened lineshapes A(ħω),
#     and PL / absorption spectra.
#  4. Compare displacement-based and force-based HR calculations.
#  5. Produce and save several diagnostic plots.
#
# IMPORTANT NOTES BEFORE RUNNING:
#  - This script assumes PyPL is installed and that helper parsers
#    (parse_phonopy_h5, parse_atoms_qexml, parse_forces_qexml) exist and
#    behave as used below.
#  - Carefully check units for phonon frequencies returned by parse_phonopy_h5.
#    Many libraries use THz (cycles/s) while some use angular frequency (rad/s).
#    Conversions to energy use hbar and may require a 2π factor if freqs are cycles/s.
#  - Make sure all referenced files exist relative to your working directory.
# ---------------------------------------------------------------------------

import numpy as np
from scipy import constants    # contains physical constants: hbar, eV, k, etc.
import matplotlib.pyplot as plt
import os

# Make plot font slightly bigger for readability in saved images
plt.rcParams.update({'font.size': 12})

# ---------------------------
# Optional: color palette
# ---------------------------
# Simple hex colors used for plotting. Not essential for functionality.
blue = '#4285F4'
red = '#DB4437'
deep_violet = '#9D80B8'
vibrant_purple = '#D391C2'
soft_coral = '#F3B1BA'
ocean_blue = '#81BFD6'
goldenrod_yellow = '#EBD68F'
emerald_green = '#02A650'

# ------------------------------------------------------------------------
# PyPL imports: hr_solver and utilities
# - hr_solver: high-level class providing HR calculations, spectral densities,
#   and lineshape computation.
# - pypl.utils: expected to provide parsers for phonopy and QE XML outputs.
# NOTE: If these imports raise an error, ensure the pypl package is installed
#       in your active Python environment and is on PYTHONPATH.
# ------------------------------------------------------------------------
from pypl.hr_solver import hr_solver
from pypl.utils import *   # convenience parsers (parse_phonopy_h5, parse_atoms_qexml, parse_forces_qexml)

# Ensure working directory is the location from which the script was called.
# This is defensive: relative file paths below are resolved from os.getcwd().
path = os.getcwd()
os.chdir(path)

# ------------------------------------------------------------------------
# FILE PATHS (user configurable)
# - gs_phonon_file: phonopy hdf5 containing ground-state frequencies & eigenvectors
# - gs_file: ground-state geometry XML (Quantum ESPRESSO's pwscf.xml)
# - es_file: excited-state geometry XML (computed separately, e.g. constrained-DFT)
#
# NOTE: The strings '³A₂' and '³E' in comments are **state labels**:
#       - superscript 3: triplet spin multiplicity
#       - A₂ / E: irreducible representations of C3v (symmetry labels)
#       These labels describe symmetry/spin, not file content or units.
# ------------------------------------------------------------------------
gs_phonon_file = 'phonon/gs-relaxed/gs_ph_mesh.hdf5'#NOTE: CHANGE THIS TO GS-RELAXED when finished
gs_file = 'relax/pwscf.xml'
es_file = 'tddft/final_geo/pwscf.xml'

# -------------------------
# PARSE INPUT FILES
# -------------------------
# parse_phonopy_h5(hdf5_file) should return:
#   - freqs : array-like of phonon frequencies (units must be checked; script
#             author noted "THz" but that could be cycles/s)
#   - modes : array-like of eigenvectors (mode shapes) corresponding to freqs
#
# parse_atoms_qexml(xml_file) should return:
#   - atomic_symbols : list of atomic labels (strings) in atom-ordering matching modes
#   - coords : Nx3 array of Cartesian coordinates in Angstrom
#   - cell_parameters : 3x3 matrix containing lattice vectors in Angstrom
#
# WARNING about units:
#   - If freqs are in THz (cycles/s), convert to angular freq ω = 2π f before
#     converting to energy using E = ħω. If freqs are already ω (rad/s),
#     do not multiply by 2π.
# -------------------------------------------------------------------------
gs_phonon_freqs, gs_phonon_modes = parse_phonopy_h5(gs_phonon_file)
atomic_symbols, gs_coord, cell_parameters = parse_atoms_qexml(gs_file)
atomic_symbols_2, es_coord, cell_parameters_2 = parse_atoms_qexml(es_file)

# -------------------------
# SANITY CHECKS
# -------------------------
# These assertions ensure that the ground and excited structures refer to the
# same atoms (same order in the arrays) and identical cell geometry. If these
# fail, the HR projection will be invalid because the mode basis is not aligned.
# -------------------------------------------------------------------------
assert atomic_symbols == atomic_symbols_2, (
    'Mismatch in atomic symbols between ground-state and excited-state structures.'
)
# cell_parameters are compared with a very tight tolerance. For defect
# supercells they should match exactly; if you performed geometry relaxation
# that changed the cell, you must handle that explicitly.
assert np.max(np.abs(cell_parameters - cell_parameters_2)) < 1e-12, (
    'Mismatch in cell parameters between ground-state and excited-state structures.'
)

# -------------------------
# MASSES (atomic mass dictionary)
# -------------------------
# Provide masses in atomic mass units (amu). Only include species present.
# PyPL will convert to the units it requires internally to compute
# mass-weighted displacements and HR factors.
# -------------------------------------------------------------------------
mass_list = {'C': 12.0107, 'N': 14.0067}

# ============================================================
# COMPUTE HUANG–RHYS FACTORS USING DISPLACEMENTS (Δr method)
# ============================================================
# Approach (displacement-based):
#   1. Take phonon normal modes from ground-state phonopy output.
#   2. Compute the geometric difference Δr = (r_excited - r_ground).
#   3. Project Δr onto each normal mode (mass-weighted) to obtain mode
#      displacements and the dimensionless Huang–Rhys factor S_k for each mode.
#
# hr_solver.compute_hrf_dis(...) should perform these steps and return a dict
# containing at minimum 'hr_factors' and 'freqs'. Check PyPL docs for more keys.
# ============================================================
pl_use_dis = hr_solver()
hrf_dict_pl_dis = pl_use_dis.compute_hrf_dis(
    gs_phonon_freqs,      # phonon frequencies from phonopy
    gs_phonon_modes,      # phonon eigenvectors from phonopy
    atomic_symbols,       # list of atomic labels (must match order used in modes)
    gs_coord,             # ground-state Cartesian coordinates (Angstrom)
    es_coord,             # excited-state Cartesian coordinates (Angstrom)
    cell_parameters,      # lattice vectors (Angstrom)
    mass_list=mass_list   # masses by atomic species (amu)
)

# The returned dictionary (hrf_dict_pl_dis) typically contains:
#   - 'hr_factors' : numpy array (S_k) for each normal mode k
#   - 'freqs'      : numpy array frequencies consistent with input (units need attention)
#   - possibly 'modes', 'dm' (mode-projected displacements), and reorg energies, etc.

# --------------------------------------------------------------------
# SHOW TOP CONTRIBUTING MODES (largest Huang–Rhys factors)
# --------------------------------------------------------------------
# We sort the HR factors, select the largest 10 contributors, and print them.
# For human readability we convert a frequency value into an energy (meV).
#
# IMPORTANT: The expression below multiplies frequency * ħ / eV to get eV (then *1e3 → meV).
# This is correct only if `freqs` are angular frequencies ω in units of 1/s.
# If `freqs` are cycles/s (Hz), multiply by 2π: ω = 2π f.
# If `freqs` are in THz (10^12 Hz), account for that scaling.
# --------------------------------------------------------------------
indices = np.argsort(hrf_dict_pl_dis['hr_factors'])[-10:]  # indices of top 10 HRF
for i in indices[::-1]:
    # i+1 to print 1-based index (human-friendly).
    freq = hrf_dict_pl_dis['freqs'][i]  # check units! (THz? Hz? rad/s?)
    # Convert frequency to energy (meV). If freqs are cycles/s, do 2*np.pi*freq here.
    energy_meV = freq * constants.hbar / constants.eV * 1e3
    print('index %4d    freq: % .4f meV    hrf: % .4e' %
          (i + 1, energy_meV, hrf_dict_pl_dis['hr_factors'][i]))

# --------------------------------------------------------------------
# BUILD SPECTRAL DENSITIES FOR DIFFERENT BROADENING PARAMETERS
# --------------------------------------------------------------------
# Spectral density S(ħω) is typically a sum over modes:
#   S(ħω) = sum_k S_k * delta(ħω - ħω_k)
# compute_spectral_density creates a smoothed representation by replacing
# delta functions with Gaussians or other broadeners controlled by 'sigma'.
#
# Here `sigmas` is a list of parameters; their exact meaning depends on PyPL:
#   - first element often sets energy-axis Gaussian width (meV)
#   - second could control mode-specific smoothing or HR smoothing
# Always consult PyPL docs for precise interpretation.
# --------------------------------------------------------------------
sigmas = [
    [3, 2], [3, 4], [3, 6],
    [6, 2], [6, 4], [6, 6],
]

spectral_density_pl_dis = []
for i, sigma in enumerate(sigmas):
    energy_axis = np.linspace(0, 200, 201)  # 0 → 200 meV in 1-meV steps
    tmp = pl_use_dis.compute_spectral_density(
        hrf_dict_pl_dis,
        energy_axis=energy_axis,
        sigma=sigma
    )
    # Many PyPL functions return (axis, values) tuple. The original script
    # treats tmp[0] as axis and tmp[1] as spectral density values.
    spectral_density_pl_dis.append(tmp)

# ---------------------------
# PLOT SPECTRAL DENSITIES
# ---------------------------
# Create a single plot showing multiple broadenings for comparison.
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 1, figsize=(7, 4))

colors = [
    deep_violet,
    vibrant_purple,
    soft_coral,
    ocean_blue,
    goldenrod_yellow,
    emerald_green,
]


max_val = 0.0

for i in range(len(spectral_density_pl_dis)):
    axis = spectral_density_pl_dis[i][0]
    values = spectral_density_pl_dis[i][1]
    ax.plot(axis, values, color=colors[i],
            label=r'$\sigma=[%.1f, %.1f]$' % (sigmas[i][0], sigmas[i][1]))
    max_val = max(max_val, np.max(values))  # track the highest y value

ax.set_xlabel('$\\hbar\\omega$ (meV)')
ax.set_ylabel('$S(\\hbar\\omega)$ (1/meV)')
ax.set_xlim(0,200)
ax.set_ylim(0, max_val * 1.05)  # set y-limit slightly above max
ax.legend()
os.makedirs("images/huang_rhys", exist_ok=True)
plt.savefig("images/huang_rhys/spectral_density.png", bbox_inches='tight', dpi=200)

# -------------------------------------------------------------------------
# COMPUTE LINESHAPE VIA FFT (time-domain autocorrelation → Fourier Transform)
# -------------------------------------------------------------------------
# compute_lineshape_fft commonly:
#  - constructs the thermal time autocorrelation of the optical transition
#    using HR factors and Boltzmann weights (temperature-dependent)
#  - Fourier transforms the correlation to obtain the frequency-domain
#    lineshape A(ħω).
#
# Parameters:
#  - temp: temperature in Kelvin (affects phonon occupation factors)
#  - sigma: broadening parameters (same caveat as before)
#  - zpl_broadening: additional broadening applied to zero-phonon line (meV)
# -------------------------------------------------------------------------
linshape_fft_pl_dis = pl_use_dis.compute_lineshape_fft(
    hrf_dict_pl_dis,
    temp=4,              # low temperature calculation (4 K)
    sigma=[6, 2],        # chosen broadening parameters (tweak & test)
    zpl_broadening=0.3   # broaden the ZPL by 0.3 meV
)

# Plot the computed lineshape function A(ħω)
fig, ax = plt.subplots(1, 1, figsize=(7, 4))
ax.plot(linshape_fft_pl_dis[0], linshape_fft_pl_dis[1], color=blue)
ax.set_xlabel('$\\hbar\\omega$ (meV)')
ax.set_ylabel('$A(\\hbar\\omega)$ (arb. units; check docs)')   # units are implementation-dependent
ax.set_ylim([0, 10])
ax.set_xlim([-100, 800])
plt.savefig("images/huang_rhys/lineshape_function.png", bbox_inches='tight', dpi=200)

# -------------------------------------------------------------------------
# COMPUTE PHOTOLUMINESCENCE (PL) SPECTRUM USING LINESHAPE
# -------------------------------------------------------------------------
# ezpl is the zero-phonon line energy (meV). For NV⁻ the ZPL ~ 637 nm ≈ 1945 meV.
# compute_spectrum uses the provided lineshape and ZPL location to construct
# the PL spectrum (possibly applying detailed balance factors for emission).
# -------------------------------------------------------------------------
ezpl = 1.831592 * 1000  # meV (≈637 nm experimentally) USING ANALYZE.OUT FOR ZPL
spectrum_pl_dis = pl_use_dis.compute_spectrum(
    ezpl,
    spectrum_type='PL',    # 'PL' for photoluminescence; some APIs accept 'Abs' for absorption
    lineshape=linshape_fft_pl_dis
)

fig, ax = plt.subplots(1, 1, figsize=(7, 4))
ax.plot(spectrum_pl_dis[0], spectrum_pl_dis[1], color=blue)
ax.set_xlabel('$\\hbar\\omega$ (meV)')
ax.set_ylabel('PL (arb. unit)')
ax.set_xlim([1300, 2000])
ax.set_ylim([0, 10])
plt.savefig("images/huang_rhys/pl_spectrum.png", bbox_inches='tight', dpi=200)

print("Saved plots: spectral_density.png, lineshape_function.png, pl_spectrum.png")

# -------------------------------------------------------------------------
# ALTERNATIVE HR CALCULATION: USING FORCES (ΔF method)
# -------------------------------------------------------------------------
# Two common ways to compute the mode coupling/reorganization:
#  1) Displacements: project geometric change Δr onto normal modes (used above)
#  2) Forces: compute Hellmann–Feynman forces at the excited-state geometry
#     using the ground-state electronic configuration (or vice versa), and
#     project these forces onto modes to obtain equivalent HR information.
#
# compute_hrf_forces(...) expects atomic forces at the excited-state geometry
# but computed with the ground-state potential (or the appropriate choice).
# The file forces_fname should contain forces for each atom in the same ordering.
# -------------------------------------------------------------------------
forces_fname = 'phonon/es/pwscf.xml'
atomic_symbols, gs_forces_es_coord = parse_forces_qexml(forces_fname)

# Create a new hr_solver instance (or reuse the previous one) and compute HR via forces
pl_using_forces = hr_solver()
hrf_dict_pl_forces = pl_using_forces.compute_hrf_forces(
    gs_phonon_freqs,
    gs_phonon_modes,
    atomic_symbols,
    gs_forces_es_coord,
    mass_list
)

# Build lineshape and spectrum from the forces-based HR dictionary
lineshape_fft_pl_forces = pl_using_forces.compute_lineshape_fft(
    hrf_dict_pl_forces,
    temp=4,
    sigma=[6, 2],
    zpl_broadening=0.3
)
spectrum_pl_forces = pl_using_forces.compute_spectrum(
    ezpl,
    lineshape=lineshape_fft_pl_forces,
    spectrum_type='PL'
)

# Plot comparison of displacement-based vs. forces-based PL
fig, ax = plt.subplots(1, 1, figsize=(7, 4))
ax.plot(spectrum_pl_dis[0], spectrum_pl_dis[1], label='use displacements', color=blue)
ax.plot(spectrum_pl_forces[0], spectrum_pl_forces[1], label='use forces', linestyle='--', color=red)
ax.legend()
ax.set_xlabel('$\\hbar\\omega$ (meV)')
ax.set_ylabel('PL (arb. unit)')
ax.set_xlim([1300, 2000])
ax.set_ylim([0, 10])
plt.savefig("images/huang_rhys/pl_forces.png", bbox_inches='tight', dpi=200)

# -------------------------------------------------------------------------
# ABSORPTION SPECTRUM (USING EXCITED-STATE PHONONS)
# -------------------------------------------------------------------------
# Optionally compute absorption spectrum using phonons/eigenvectors computed
# at the excited-state geometry. For this, we parse the excited-state phonon hdf5.
# This may be informative because excited-state geometry often has different
# normal modes and frequencies (mode softening/hardening).
# -------------------------------------------------------------------------
es_phonon_fname = 'phonon/es/es_ph_mesh.hdf5'
es_phonon_freqs, es_phonon_modes = parse_phonopy_h5(es_phonon_fname)

# Compute HR for absorption (projecting displacements appropriate for absorption)
abs_use_dis = hr_solver()
hrf_dict_abs_dis = abs_use_dis.compute_hrf_dis(
    es_phonon_freqs,
    es_phonon_modes,
    atomic_symbols,
    gs_coord,     # notice ordering: this may be ground → excited or vice versa based on API
    es_coord,
    cell_parameters,
    mass_list=mass_list
)

# Spectral density for absorption (single sigma choice)
spectral_density_abs_dis = abs_use_dis.compute_spectral_density(
    hrf_dict_abs_dis,
    energy_axis=np.linspace(0, 200, 201),
    sigma=[6.0, 2.0]
)

# Plot absorption spectral density vs PL spectral density (one of the PL sigmas)
fig, ax = plt.subplots(1, 1, figsize=(7, 4))
ax.plot(spectral_density_abs_dis[0], spectral_density_abs_dis[1], label='Abs, use displacements', color=blue)
# choose a representative PL spectral density (index 3 earlier was one sigma set)
ax.plot(spectral_density_pl_dis[3][0], spectral_density_pl_dis[3][1], label='PL, use displacements', color=red)
ax.legend()
ax.set_xlabel('$\\hbar\\omega$ (meV)')
ax.set_ylabel('$S(\\hbar\\omega)$ (meV)')
ax.set_xlim(0, 200)
ax.set_ylim(0.0, 0.1)
plt.savefig("images/huang_rhys/abs_spectrum.png", bbox_inches='tight', dpi=200)

# Compute lineshape and absorption spectrum (frequency-domain)
lineshape_fft_abs_dis = abs_use_dis.compute_lineshape_fft(
    hrf_dict_abs_dis,
    temp=4,
    sigma=[6.0, 2.0],
    zpl_broadening=0.3
)
spectrum_abs_dis = abs_use_dis.compute_spectrum(
    ezpl,
    spectrum_type='Abs',
    lineshape=lineshape_fft_abs_dis
)

# Plot absorption vs emission (PL) on one panel for comparison
fig, ax = plt.subplots(1, 1, figsize=(7, 4))
ax.plot(spectrum_abs_dis[0], spectrum_abs_dis[1], label='Absorption', color=blue)
ax.plot(spectrum_pl_dis[0], spectrum_pl_dis[1], label='PL', color=red)
ax.legend()
ax.set_xlabel('$\\hbar\\omega$ (meV)')
ax.set_ylabel('Absorption / PL (arb. unit)')
ax.set_xlim(1500, 2400)
ax.set_ylim([0, 10])
plt.savefig("images/huang_rhys/pl_absorption.png", bbox_inches='tight', dpi=200)

# -------------------------------------------------------------------------
# TEMPERATURE-DEPENDENT PL (USING FORCES-BASED HR)
# -------------------------------------------------------------------------
# Demonstrate how the lineshape & PL spectrum change with temperature by
# computing lineshapes at multiple temperatures. Often the ZPL broadens and
# phonon sidebands change intensity with temperature.
# -------------------------------------------------------------------------
temp_pl_forces = hr_solver()
hrf_dict_temp_pl_forces = temp_pl_forces.compute_hrf_forces(
    gs_phonon_freqs,
    gs_phonon_modes,
    atomic_symbols,
    gs_forces_es_coord,
    mass_list
)

# Example temperature set and corresponding ZPL broadenings (empirical choices)
# zpl_broadenings chosen here are illustrative; match to experiment if available.
temps = [8, 150, 200, 250, 300]
zpl_broadenings = [0.3, 0.8, 1.8, 3.2, 4.8]  # in meV, increasing with temperature
temp_pls = []

# Compute lineshape & spectrum for each temperature and store
for t, l in zip(temps, zpl_broadenings):
    tmp_lineshape = temp_pl_forces.compute_lineshape_fft(
        hrf_dict_temp_pl_forces,
        temp=t,
        sigma=[6.0, 2.0],
        zpl_broadening=l
    )
    tmp_spectrum = temp_pl_forces.compute_spectrum(
        ezpl,
        spectrum_type='PL',
        lineshape=tmp_lineshape
    )
    temp_pls.append(tmp_spectrum)

# Plot stacked/offset PL spectra to visualize temperature dependence
fig, ax = plt.subplots(1, 1, figsize=(7, 6))

colors = [deep_violet, vibrant_purple, soft_coral, ocean_blue, goldenrod_yellow]

# Plot from highest temperature first or last depending on visual preference.
# Here we loop from highest index to lowest to stack plots nicely.
for i in range(len(temps) - 1, -1, -1):
    axis = temp_pls[i][0]
    values = temp_pls[i][1]
    # vertically offset the curves by i units so they don't overlap; label includes T and ZPL broadening
    ax.fill_between(axis, 1 * i, values + 1 * i, color=colors[i],
                    label='$T = %d$ K, $\\lambda =$ %.1f meV' % (temps[i], zpl_broadenings[i]))
ax.legend(loc='upper left')
ax.set_xlabel('$\\hbar\\omega$ (meV)')
ax.set_ylabel('PL (arb. unit)')
ax.set_xlim([1400, 2100])
ax.set_ylim([0, 10 + len(temps)])  # allow vertical room for stacked plots
plt.savefig("images/huang_rhys/td_pl.png", bbox_inches='tight', dpi=200)

print("All done. Check the images/huang_rhys/ directory for generated plots.")

# ---------------------------------------------------------------------------
# QUICK TROUBLESHOOTING / DEBUG CHECKLIST (for you while learning)
# ---------------------------------------------------------------------------
# - If outputs are zero or NaN:
#     * confirm parse functions returned sensible arrays (print shapes, min/max)
#     * check that hr_factors are nonzero (if all zeros, the displacement projection failed)
# - If spectral energies are off by ~6.28x:
#     * frequency units are cycles/s → need factor 2π to convert to ω
# - If plots show no ZPL but only sidebands:
#     * check zpl_broadening and how compute_spectrum handles ZPL intensity
# - If parse_forces_qexml fails:
#     * ensure forces file corresponds to the same atom ordering and units (Ry/Bohr? eV/Å?)
# - If memory/time is large:
#     * reduce energy_axis resolution, reduce FFT grid sizes, or work with fewer modes for testing
# ---------------------------------------------------------------------------

# End of tutorial-style annotated script.
