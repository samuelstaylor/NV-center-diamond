#!/usr/bin/env python3
"""
Photoluminescence (PL) and absorption spectra of the nitrogen-vacancy (NV–) center
in diamond using the 1D configurational coordinate diagram (CCD) approach.

This script is a full Python version of PyPL tutorial 002, written as a
standalone, fully commented script instead of a Jupyter notebook.

Author: (your name)
"""

import numpy as np
from scipy import constants
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

# PyPL imports
from pypl.config_coord_1d_solver import config_coord_1d_solver
from pypl.utils import (
    parse_atoms_qexml,
    parse_total_energy_qexml,
    parse_phonopy_h5
)
from pypl.hr_solver import hr_solver


# ------------------------------------------------------------------------------
# Plot style
# ------------------------------------------------------------------------------
plt.rcParams.update({"font.size": 12})
blue = "#4285F4"
red = "#DB4437"


# ------------------------------------------------------------------------------
# 1. LOAD STRUCTURES FOR GROUND AND EXCITED STATES
# ------------------------------------------------------------------------------
print("\n=== Parsing atomic coordinates ===")

gs_xml = "relax/pwscf.xml"
es_xml = "tddft/final_geo/pwscf.xml"

atomic_symbols, gs_coord, cell_parameters = parse_atoms_qexml(gs_xml)
atomic_symbols_2, es_coord, cell_parameters_2 = parse_atoms_qexml(es_xml)

# Ensure consistent ordering between GS and ES structures
assert atomic_symbols == atomic_symbols_2, "GS/ES atomic symbols differ!"
assert np.max(np.abs(cell_parameters - cell_parameters_2)) < 1e-12, \
    "GS/ES cell parameters differ!"

# Atomic masses for mass-weighted coordinates
mass_list = {"C": 12.0107, "N": 14.0067}
all_masses = np.array([mass_list[s] for s in atomic_symbols])

# Compute mass-weighted displacement ΔQ
delta_q = np.linalg.norm((es_coord - gs_coord) * np.sqrt(all_masses[:, None]))
print(f"ΔQ = {delta_q:.10e} amu^(1/2) Å")


# ------------------------------------------------------------------------------
# 2. PARSE SINGLE-POINT ENERGIES ALONG INTERPOLATION PATH
# ------------------------------------------------------------------------------
print("\n=== Parsing energy profiles from GS and ES calculations ===")

relative_coordinates = np.linspace(-0.2, 1.2, 15)

gs_energies = np.zeros(len(relative_coordinates))
es_energies = np.zeros(len(relative_coordinates))

gs_rel_coord = []
es_rel_coord = []

for i in range(len(relative_coordinates)):
    # Ground state
    gs_file = f"002_nv_diamond_1d_ccd/gs_1d_ccd/Image-{i+1}/pwscf.xml"
    if os.path.exists(gs_file):
        gs_rel_coord.append(i)
        gs_energies[i] = parse_total_energy_qexml(gs_file)

    # Excited state
    es_file = f"002_nv_diamond_1d_ccd/es_1d_ccd/Image-{i+1}/pwscf.xml"
    if os.path.exists(es_file):
        es_rel_coord.append(i)
        es_energies[i] = parse_total_energy_qexml(es_file)

# Reset zero of energy to minimum GS energy
ref = np.min(gs_energies)
gs_energies = gs_energies[gs_rel_coord] - ref
es_energies = es_energies[es_rel_coord] - ref

# Convert to mass-weighted coordinate Q
gs_1d_coord = delta_q * relative_coordinates[gs_rel_coord]
es_1d_coord = delta_q * relative_coordinates[es_rel_coord]


# ------------------------------------------------------------------------------
# 3. FIT QUADRATIC POTENTIALS TO GS AND ES ENERGIES
# ------------------------------------------------------------------------------
print("\n=== Fitting parabolic potentials ===")

def gs_fit_fun(x, a, c):
    return a * x**2 + c

def es_fit_fun(x, a, c):
    return a * (x - delta_q)**2 + c

# Ground state fit (first ~5 points near minimum)
gs_params, _ = curve_fit(gs_fit_fun, gs_1d_coord[:5], gs_energies[:5])
print("GS fit parameters:", gs_params)

# Excited state fit (last ~10 points near minimum)
es_params, _ = curve_fit(es_fit_fun, es_1d_coord[4:], es_energies[4:])
print("ES fit parameters:", es_params)

# Evaluate fits
gs_fit_energies = gs_fit_fun(gs_1d_coord, *gs_params)
es_fit_energies = es_fit_fun(es_1d_coord, *es_params)


# ------------------------------------------------------------------------------
# 4. EXTRACT EFFECTIVE PHONON FREQUENCIES (IN meV)
# ------------------------------------------------------------------------------
def curvature_to_meV(a):
    """
    Convert quadratic coefficient 'a' [eV/(amu*Å^2)] to phonon energy in meV.
    """
    mu = constants.physical_constants["atomic mass constant"][0]  # kg
    factor = 2 * a * constants.eV / (1e-20 * mu)  # rad^2/s^2
    factor *= (constants.hbar**2 / constants.eV**2)
    return np.sqrt(factor) * 1e3  # meV

gs_phonon = curvature_to_meV(gs_params[0])
es_phonon = curvature_to_meV(es_params[0])

print(f"GS phonon = {gs_phonon:.5f} meV")
print(f"ES phonon = {es_phonon:.5f} meV")


# ------------------------------------------------------------------------------
# 5. COMPUTE HUANG–RHYS FACTORS
# ------------------------------------------------------------------------------
def compute_hr(delta_q, phonon_meV):
    mu = constants.physical_constants["atomic mass constant"][0]
    omega = phonon_meV * 1e-3 * constants.eV / constants.hbar
    return delta_q**2 * 1e-20 * mu * omega / (2 * constants.hbar)

gs_hrf = compute_hr(delta_q, gs_phonon)
es_hrf = compute_hr(delta_q, es_phonon)

print(f"GS Huang–Rhys factor: {gs_hrf:.5f}")
print(f"ES Huang–Rhys factor: {es_hrf:.5f}")


# ------------------------------------------------------------------------------
# 6. PLOT THE CONFIGURATIONAL COORDINATE DIAGRAM
# ------------------------------------------------------------------------------
plt.figure(figsize=(5, 7))
plt.plot(es_1d_coord, es_energies, "s", color=blue, label="ES energies")
plt.plot(gs_1d_coord, gs_energies, "o", color=red, label="GS energies")

plt.plot(es_1d_coord, es_fit_energies, "--", color=blue)
plt.plot(gs_1d_coord, gs_fit_energies, "--", color=red)

plt.axvline(0.0, color="gray", ls="--")
plt.axvline(delta_q, color="gray", ls="--")

plt.text(0.46, 0.4, f"GS phonon: {gs_phonon:.2f} meV\nGS HRF: {gs_hrf:.2f}",
         transform=plt.gca().transAxes, color=red)
plt.text(0.46, 0.6, f"ES phonon: {es_phonon:.2f} meV\nES HRF: {es_hrf:.2f}",
         transform=plt.gca().transAxes, color=blue)

plt.xlabel("Q (amu$^{1/2}$ Å)")
plt.ylabel("Energy (eV)")
plt.legend()
plt.tight_layout()
os.makedirs("images/1d_ccd", exist_ok=True)
plt.savefig("images/1d_ccd/plot_1.png", bbox_inches='tight', dpi=200)


# ------------------------------------------------------------------------------
# 7. COMPUTE PL AND ABSORPTION SPECTRA USING THE 1D CCD SOLVER
# ------------------------------------------------------------------------------
print("\n=== Computing PL and absorption spectra (1D CCD) ===")

order_es = 50
order_gs = 60
ene_range = [-200, 1200]  # meV
resol = 1401
temp = 5
gamma = 0.3
sigma = [8, 25]
ezpl = 1945  # meV

# PL spectrum: ES → GS
ccd_pl = config_coord_1d_solver(es_phonon, gs_phonon, delta_q)
ccd_pl.compute_franck_condon_integrals(ni=order_es, nf=order_gs)
ccd_pl.bulid_fc_lsp(
    eneaxis=np.linspace(ene_range[0], ene_range[1], resol),
    temp=temp, sigma=sigma, zpl_lorentzian=True, gamma=gamma
)
pl_spectrum = ccd_pl.compute_spectrum(tdm=1.0, zpl=ezpl, spectrum_type="PL")

# Absorption spectrum: GS → ES
ccd_abs = config_coord_1d_solver(gs_phonon, es_phonon, delta_q)
ccd_abs.compute_franck_condon_integrals(ni=order_gs, nf=order_es)
ccd_abs.bulid_fc_lsp(
    eneaxis=np.linspace(ene_range[0], ene_range[1], resol),
    temp=temp, sigma=sigma, zpl_lorentzian=True, gamma=gamma
)
abs_spectrum = ccd_abs.compute_spectrum(tdm=1.0, zpl=ezpl, spectrum_type="Abs")


# ------------------------------------------------------------------------------
# 8. PLOT PL AND ABSORPTION SPECTRA (1D CCD)
# ------------------------------------------------------------------------------
plt.figure(figsize=(12, 4))

# Normalize spectra
pl_norm = pl_spectrum[1] / (np.sum(pl_spectrum[1]) * abs(pl_spectrum[0][1] - pl_spectrum[0][0]))
abs_norm = abs_spectrum[1] / (np.sum(abs_spectrum[1]) * abs(abs_spectrum[0][1] - abs_spectrum[0][0]))

plt.plot(pl_spectrum[0] * 1e-3, pl_norm * 1e3, color=red, label="PL (1D CCD)")
plt.plot(abs_spectrum[0] * 1e-3, abs_norm * 1e3, color=blue, label="Abs (1D CCD)")

plt.xlim(1.5, 2.4)
plt.ylim(0, 6)

plt.xlabel("ℏω (eV)")
plt.ylabel("Intensity (arb. units)")
plt.legend()
plt.grid(ls="--", color="gray")
plt.tight_layout()
plt.savefig("images/1d_ccd/plot_2.png", bbox_inches='tight', dpi=200)




# ------------------------------------------------------------------------------
# 9. COMPUTE FULL PHONON H-R THEORY SPECTRA FOR COMPARISON
# ------------------------------------------------------------------------------
print("\n=== Computing full HR spectra (all phonon modes) ===")

gs_phonon_file = "001_nv_diamond_abs_pl/phonon/gs_ph_mesh.hdf5"
es_phonon_file = "001_nv_diamond_abs_pl/phonon/es_ph_mesh.hdf5"
gs_file = "001_nv_diamond_abs_pl/gs_dft/pwscf.xml"
es_file = "001_nv_diamond_abs_pl/es_cdft/pwscf.xml"

gs_freqs, gs_modes = parse_phonopy_h5(gs_phonon_file)
es_freqs, es_modes = parse_phonopy_h5(es_phonon_file)

atomic_symbols, gs_coord, cell_parameters = parse_atoms_qexml(gs_file)
_, es_coord, _ = parse_atoms_qexml(es_file)

hr_pl = hr_solver()
hr_abs = hr_solver()

hrf_pl = hr_pl.compute_hrf_dis(gs_freqs, gs_modes, atomic_symbols, gs_coord, es_coord, cell_parameters, mass_list=mass_list)
lin_pl = hr_pl.compute_lineshape_fft(hrf_pl, temp=4, sigma=[6, 2], zpl_broadening=0.3)
spectrum_pl_dis = hr_pl.compute_spectrum(ezpl, spectrum_type="PL", lineshape=lin_pl)

hrf_abs = hr_abs.compute_hrf_dis(es_freqs, es_modes, atomic_symbols, gs_coord, es_coord, cell_parameters, mass_list=mass_list)
lin_abs = hr_abs.compute_lineshape_fft(hrf_abs, temp=4, sigma=[6, 2], zpl_broadening=0.3)
spectrum_abs_dis = hr_abs.compute_spectrum(ezpl, spectrum_type="Abs", lineshape=lin_abs)


# ------------------------------------------------------------------------------
# 10. PLOT COMPARISON: 1D CCD vs FULL HR THEORY
# ------------------------------------------------------------------------------
plt.figure(figsize=(12, 4))

plt.plot(pl_spectrum[0]*1e-3, pl_norm*1e3, color=red, label="PL (1D CCD)")
plt.plot(abs_spectrum[0]*1e-3, abs_norm*1e3, color=blue, label="Abs (1D CCD)")

# Normalize HR spectra
pl_norm_hr = spectrum_pl_dis[1] / (np.sum(spectrum_pl_dis[1]) * abs(spectrum_pl_dis[0][1] - spectrum_pl_dis[0][0]))
abs_norm_hr = spectrum_abs_dis[1] / (np.sum(spectrum_abs_dis[1]) * abs(spectrum_abs_dis[0][1] - spectrum_abs_dis[0][0]))

plt.plot(spectrum_pl_dis[0]*1e-3, pl_norm_hr*1e3, "--", color=red, label="PL (HR full)")
plt.plot(spectrum_abs_dis[0]*1e-3, abs_norm_hr*1e3, "--", color=blue, label="Abs (HR full)")

plt.xlim(1.5, 2.4)
plt.ylim(0, 6)
plt.xlabel("ℏω (eV)")
plt.ylabel("Intensity (arb. units)")
plt.legend()
plt.grid(ls="--", color="gray")
plt.tight_layout()
plt.savefig("images/1d_ccd/plot_2.png", bbox_inches='tight', dpi=200)

print("\n=== DONE ===\n")
