import numpy as np
from scipy import constants
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})

blue = '#4285F4'
red = '#DB4437'

from pypl.config_coord_1d_solver import config_coord_1d_solver
from pypl.utils import *
import os

path = os.getcwd()
os.chdir(path)

# Unit of coordinates and cell_parameters is Angstrom
atomic_symbols, gs_coord, cell_parameters = parse_atoms_qexml('002_nv_diamond_1d_ccd/gs_dft/pwscf.xml')
atomic_symbols_2, es_coord, cell_parameters_2 = parse_atoms_qexml('002_nv_diamond_1d_ccd/es_cdft/pwscf.xml')

# Ensure that both sets of atomic coordinates have consistent atomic symbols and cell parameters.
assert atomic_symbols == atomic_symbols_2, 'Mismatch in atomic symbols between ground-state and excited-state structures.'
assert np.max(np.abs(cell_parameters - cell_parameters_2)) < 1e-12, 'Mismatch in cell parameters between ground-state and excited-state structures.'
mass_list = {'C': 12.0107, 'N': 14.0067}
all_masses = np.array([mass_list[sym] for sym in atomic_symbols])
delta_q = np.linalg.norm((es_coord - gs_coord) * all_masses[:, None]**0.5)
print('\Delta Q = % .10e amu^{0.5} \AA' % (delta_q))

relative_coordinates = np.linspace(-0.2, 1.2, 15)

gs_rel_coord = []
es_rel_coord = []
# energies are in eV
gs_energies = np.zeros(relative_coordinates.shape[0])
es_energies = np.zeros(relative_coordinates.shape[0])

for i in range(relative_coordinates.shape[0]):
    fileName = '002_nv_diamond_1d_ccd/gs_1d_ccd/Image-%d/pwscf.xml' % (i + 1)
    if os.path.exists(fileName):
        gs_rel_coord.append(i)
        gs_energies[i] = parse_total_energy_qexml(fileName)
    fileName = '002_nv_diamond_1d_ccd/es_1d_ccd/Image-%d/pwscf.xml' % (i + 1)
    if os.path.exists(fileName):
        es_rel_coord.append(i)
        es_energies[i] = parse_total_energy_qexml(fileName)

# reset the energy zero
ref_energy = np.min(gs_energies)
gs_energies = gs_energies[gs_rel_coord] - np.min(ref_energy)
es_energies = es_energies[es_rel_coord] - np.min(ref_energy)

gs_1d_coord = delta_q * relative_coordinates[gs_rel_coord]
es_1d_coord = delta_q * relative_coordinates[es_rel_coord]


from scipy.optimize import curve_fit

def gs_fit_fun(x, a, c):
    return a * x**2 + c

def es_fit_fun(x, a, c):
    return a * (x - delta_q) ** 2 + c

# Ground state phonon
gs_params = curve_fit(gs_fit_fun, gs_1d_coord[:5], gs_energies[:5])[0]
print('Parameters for GS curve: ', gs_params)
gs_fit_energies = gs_fit_fun(gs_1d_coord, gs_params[0], gs_params[1])

# (rad / s)^2
gs_phonon = 2 * gs_params[0] * constants.eV / (1e-10**2 * constants.physical_constants['atomic mass constant'][0])
# eV^2
gs_phonon *= (constants.hbar**2 / constants.eV**2)
# meV
gs_phonon = np.sqrt(gs_phonon) * 1000
print('GS phonon is %.5f meV' % gs_phonon)

# Excited state phonon
es_params = curve_fit(es_fit_fun, es_1d_coord[4:], es_energies[4:])[0]
print('Parameters for ES curve: ', es_params)
es_fit_energies = es_fit_fun(es_1d_coord, es_params[0], es_params[1])

# (rad / s)^2
es_phonon = 2 * es_params[0] * constants.eV / (1e-10**2 * constants.physical_constants['atomic mass constant'][0])
# eV^2
es_phonon *= (constants.hbar**2 / constants.eV**2)
# meV
es_phonon = np.sqrt(es_phonon) * 1000
print('ES phonon is %.5f meV' % es_phonon)


gs_hrf = (
    delta_q**2 * 1e-10**2 * constants.physical_constants['atomic mass constant'][0]
    * gs_phonon * 1e-3 * constants.eV / constants.hbar
    / (2 * constants.hbar)
)
print('HR for GS is %.5f' % gs_hrf)

es_hrf = (
    delta_q**2 * 1e-10**2 * constants.physical_constants['atomic mass constant'][0]
    * es_phonon * 1e-3 * constants.eV / constants.hbar
    / (2 * constants.hbar)
)
print('HR for ES is %.5f' % es_hrf)


ax = plt.subplots(1, 1, figsize=(5, 7))

ax.plot(es_1d_coord, es_energies, color=blue, linestyle='', marker='s', markersize=5, label='ES energy')
ax.plot(gs_1d_coord, gs_energies, color=red, linestyle='', marker='o', markersize=5, label='GS energy')

ax.plot(es_1d_coord, es_fit_energies, color=blue, linestyle='--', linewidth=1.0, marker='', label='ES fit')
ax.plot(gs_1d_coord, gs_fit_energies, color=red, linestyle='--', linewidth=1.0, marker='', label='GS fit')

ax.text(x=0.46, y=0.4, s='GS phonon: %.2f meV\nGS HRF: %.2f' % (gs_phonon, gs_hrf), color=red, transform=ax.transAxes)
ax.text(x=0.46, y=0.6, s='ES phonon: %.2f meV\nES HRF: %.2f' % (es_phonon, es_hrf), color=blue, transform=ax.transAxes)

ax.axvline(x=0.0, color='gray', linestyle='--', linewidth=1)
ax.axvline(x=delta_q, color='gray', linestyle='--', linewidth=1)

ax.legend(fontsize=12, loc='center left', edgecolor='black')
ax.set_xlim((-0.2, 0.9))
plt.xlabel('Q (amu$^{1/2}$ Å)')
plt.ylabel('Total Energy (eV)')
plt.tick_params(direction='in')

plt.show()


# order (must be large enough to converge the spectra and should be tested for each system)
order_es = 50
order_gs = 60

# energy range for plot (meV)
ene_range = [-200, 1200]

# resolution
resol = 1401

# temperature (K)
temp = 5

# broadening (empirically set to best match experiment)
gamma = 0.3
sigma = [8, 25]

# ZPL (meV)
ezpl = 1945



ccd_pl = config_coord_1d_solver(es_phonon, gs_phonon, delta_q)

ccd_pl.compute_franck_condon_integrals(ni=order_es, nf=order_gs)
ccd_pl.bulid_fc_lsp(eneaxis=np.linspace(ene_range[0], ene_range[1], resol),
                    temp=temp, sigma=sigma, zpl_lorentzian=True, gamma=gamma)
pl_spectrum = ccd_pl.compute_spectrum(tdm=1.0, zpl=ezpl, spectrum_type='PL')
ccd_abs = config_coord_1d_solver(gs_phonon, es_phonon, delta_q)

ccd_abs.compute_franck_condon_integrals(ni=order_gs, nf=order_es)
ccd_abs.bulid_fc_lsp(eneaxis=np.linspace(ene_range[0], ene_range[1], resol),
                        temp=temp, sigma=sigma, zpl_lorentzian=True, gamma=gamma)

abs_spectrum = ccd_abs.compute_spectrum(tdm=1.0, zpl=ezpl, spectrum_type='Abs')


fig, ax = plt.subplots(1, 1, figsize=(12, 4))

ax.plot(pl_spectrum[0] * 1e-3, pl_spectrum[1] / (np.sum(pl_spectrum[1]) * abs(pl_spectrum[0][1] - pl_spectrum[0][0])) * 1e3,
        color=red, linewidth=1, linestyle='-', label='PL')
ax.plot(abs_spectrum[0] * 1e-3, abs_spectrum[1] / (np.sum(abs_spectrum[1]) * abs(abs_spectrum[0][1] - abs_spectrum[0][0])) * 1e3,
        color=blue, linewidth=1, linestyle='-', label='Abs')

ax.set_xlim((1.5, 2.4))
ax.set_ylim((0.0, 6))

ax.legend(fontsize=12, loc='upper right', edgecolor='black')
ax.grid(color='gray', linestyle='--', linewidth=0.5)

ax.tick_params(direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlabel('$\hbar\omega$ (eV)')
ax.set_ylabel('PL (arb. unit.)')

plt.show()


from pypl.hr_solver import hr_solver

gs_phonon_file = '001_nv_diamond_abs_pl/phonon/gs_ph_mesh.hdf5'
gs_file = '001_nv_diamond_abs_pl/gs_dft/pwscf.xml'
es_file = '001_nv_diamond_abs_pl/es_cdft/pwscf.xml'

# Unit of freqs is THz
gs_phonon_freqs, gs_phonon_modes = parse_phonopy_h5(gs_phonon_file)
# Unit of coordinates and cell_parameters is Angstrom
atomic_symbols, gs_coord, cell_parameters = parse_atoms_qexml(gs_file)
atomic_symbols_2, es_coord, cell_parameters_2 = parse_atoms_qexml(es_file)

mass_list = {'C': 12.0107, 'N': 14.0067}

pl_use_dis = hr_solver()
hrf_dict_pl_dis = pl_use_dis.compute_hrf_dis(gs_phonon_freqs, gs_phonon_modes, atomic_symbols, gs_coord, es_coord, cell_parameters, mass_list=mass_list)
linshape_fft_pl_dis = pl_use_dis.compute_lineshape_fft(hrf_dict_pl_dis, temp=4, sigma=[6, 2], zpl_broadening=0.3)
spectrum_pl_dis = pl_use_dis.compute_spectrum(ezpl, spectrum_type='PL', lineshape=linshape_fft_pl_dis)


es_phonon_fname = '001_nv_diamond_abs_pl/phonon/es_ph_mesh.hdf5'
es_phonon_freqs, es_phonon_modes = parse_phonopy_h5(es_phonon_fname)

abs_use_dis = hr_solver()
hrf_dict_abs_dis = abs_use_dis.compute_hrf_dis(es_phonon_freqs, es_phonon_modes, atomic_symbols, gs_coord, es_coord, cell_parameters, mass_list=mass_list)
lineshape_fft_abs_dis = abs_use_dis.compute_lineshape_fft(hrf_dict_abs_dis, temp=4, sigma=[6.0, 2.0], zpl_broadening=0.3)
spectrum_abs_dis = abs_use_dis.compute_spectrum(ezpl, spectrum_type='Abs', lineshape=lineshape_fft_abs_dis)


fig, ax = plt.subplots(1, 1, figsize=(12, 4))

ax.plot(pl_spectrum[0] * 1e-3, pl_spectrum[1] / (np.sum(pl_spectrum[1]) * abs(pl_spectrum[0][1] - pl_spectrum[0][0])) * 1e3,
        color=red, linewidth=1, linestyle='-', label='PL, 1D CCD')
ax.plot(abs_spectrum[0] * 1e-3, abs_spectrum[1] / (np.sum(abs_spectrum[1]) * abs(abs_spectrum[0][1] - abs_spectrum[0][0])) * 1e3,
        color=blue, linewidth=1, linestyle='-', label='Abs, 1D CCD')

ax.plot(spectrum_pl_dis[0] * 1e-3, spectrum_pl_dis[1] / (np.sum(spectrum_pl_dis[1]) * abs(spectrum_pl_dis[0][1] - spectrum_pl_dis[0][0])) * 1e3,
        label='PL, HR', color=red, linestyle='--')
ax.plot(spectrum_abs_dis[0] * 1e-3, spectrum_abs_dis[1] / (np.sum(spectrum_abs_dis[1]) * abs(spectrum_abs_dis[0][1] - spectrum_abs_dis[0][0])) * 1e3,
        label='Abs, HR', color=blue, linestyle='--')

ax.set_xlim((1.5, 2.4))
ax.set_ylim((0.0, 6))

ax.legend(fontsize=12, loc='upper right', edgecolor='black')
ax.grid(color='gray', linestyle='--', linewidth=0.5)

ax.tick_params(direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlabel('$\hbar\omega$ (eV)')
ax.set_ylabel('PL (arb. unit.)')

plt.show()
