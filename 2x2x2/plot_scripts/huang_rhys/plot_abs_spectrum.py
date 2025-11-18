import os
import matplotlib.pyplot as plt
from pypl.hr_solver import hr_solver
from pypl.utils import parse_phonopy_h5, parse_atoms_qexml
import numpy as np

# Colors
blue = '#4285F4'

# Load excited-state phonons
es_phonon_fname = 'phonon/es/es_ph_mesh.hdf5'
atomic_symbols, gs_coord, cell_parameters = parse_atoms_qexml('relax/pwscf.xml')
_, es_coord, _ = parse_atoms_qexml('tddft/final_geo/pwscf.xml')

es_phonon_freqs, es_phonon_modes = parse_phonopy_h5(es_phonon_fname)
mass_list = {'C': 12.0107, 'N': 14.0067}

# HR factors
abs_use_dis = hr_solver()
hrf_dict_abs_dis = abs_use_dis.compute_hrf_dis(es_phonon_freqs, es_phonon_modes,
                                               atomic_symbols, gs_coord, es_coord,
                                               cell_parameters, mass_list=mass_list)

# Spectral density
spectral_density_abs_dis = abs_use_dis.compute_spectral_density(hrf_dict_abs_dis,
                                                                energy_axis=np.linspace(0,200,201),
                                                                sigma=[6.0,2.0])

# Plot
fig, ax = plt.subplots(figsize=(7,4))
ax.plot(spectral_density_abs_dis[0], spectral_density_abs_dis[1], label='Abs, use displacements', color=blue)
ax.set_xlabel('$\\hbar\\omega$ (meV)')
ax.set_ylabel('$S(\\hbar\\omega)$ (meV)')
ax.set_xlim([0,200])
ax.set_ylim([0,0.1])
ax.legend()

os.makedirs("images", exist_ok=True)
plt.savefig("images/abs_spectrum.png", bbox_inches='tight', dpi=200)
print("Saved abs_spectrum.png")
