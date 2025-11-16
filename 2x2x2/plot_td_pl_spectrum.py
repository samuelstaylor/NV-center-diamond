import os
import matplotlib.pyplot as plt
from pypl.hr_solver import hr_solver
from pypl.utils import parse_phonopy_h5, parse_atoms_qexml

# Colors for plotting
deep_violet = '#9D80B8'
vibrant_purple = '#D391C2'
soft_coral = '#F3B1BA'
ocean_blue = '#81BFD6'
goldenrod_yellow = '#EBD68F'
colors = [deep_violet, vibrant_purple, soft_coral, ocean_blue, goldenrod_yellow]

# Load phonons and coordinates
gs_phonon_file = 'phonon/gs/gs-ph-mesh.hdf5'
gs_file = 'relax/pwscf.xml'
forces_fname = 'github_repos/NV-center-diamond/2x2x2/phonon/es/pwscf.xml'

gs_phonon_freqs, gs_phonon_modes = parse_phonopy_h5(gs_phonon_file)
atomic_symbols, gs_coord, cell_parameters = parse_atoms_qexml(gs_file)
atomic_symbols, gs_forces_es_coord = parse_forces_qexml(forces_fname)
mass_list = {'C': 12.0107, 'N': 14.0067}

# Forces-based HR
temp_pl_forces = hr_solver()
hrf_dict_temp_pl_forces = temp_pl_forces.compute_hrf_forces(gs_phonon_freqs, gs_phonon_modes,
                                                            atomic_symbols, gs_forces_es_coord,
                                                            mass_list)

# Temperature series
temps = [8, 150, 200, 250, 300]
zpl_broadenings = [0.3, 0.8, 1.8, 3.2, 4.8]
ezpl = 1831.592
temp_pls = []

# Compute spectra for each temperature
for t, l in zip(temps, zpl_broadenings):
    tmp_lineshape = temp_pl_forces.compute_lineshape_fft(hrf_dict_temp_pl_forces,
                                                         temp=t,
                                                         sigma=[6.0,2.0],
                                                         zpl_broadening=l)
    tmp_spectrum = temp_pl_forces.compute_spectrum(ezpl, spectrum_type='PL',
                                                   lineshape=tmp_lineshape)
    temp_pls.append(tmp_spectrum)

# Plot stacked spectra
fig, ax = plt.subplots(figsize=(7,6))
for i in range(len(temps)-1, -1, -1):
    axis, values = temp_pls[i]
    ax.fill_between(axis, 1*i, values+1*i, color=colors[i],
                    label='$T = %d$ K, $\\lambda = %.1f$ meV'%(temps[i], zpl_broadenings[i]))

ax.set_xlabel('$\\hbar\\omega$ (meV)')
ax.set_ylabel('PL (arb. unit)')
ax.set_xlim([1400,2100])
ax.set_ylim([0, 10 + len(temps)])
ax.legend(loc='upper left')

os.makedirs("images", exist_ok=True)
plt.savefig("images/td_pl.png", bbox_inches='tight', dpi=200)
print("Saved td_pl.png")
