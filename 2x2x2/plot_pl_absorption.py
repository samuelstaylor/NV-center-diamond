import os
import matplotlib.pyplot as plt
from pypl.hr_solver import hr_solver
from pypl.utils import parse_phonopy_h5, parse_atoms_qexml

blue = '#4285F4'
red = '#DB4437'

# Load PL spectral density
gs_phonon_file = 'phonon/gs/gs-ph-mesh.hdf5'
gs_file = 'relax/pwscf.xml'
es_phonon_fname = '001_nv_diamond_abs_pl/phonon/es_ph_mesh.hdf5'
atomic_symbols, gs_coord, cell_parameters = parse_atoms_qexml(gs_file)
_, es_coord, _ = parse_atoms_qexml('tddft/final_geo/pwscf.xml')

gs_phonon_freqs, gs_phonon_modes = parse_phonopy_h5(gs_phonon_file)
es_phonon_freqs, es_phonon_modes = parse_phonopy_h5(es_phonon_fname)
mass_list = {'C': 12.0107, 'N': 14.0067}

# PL HR
pl_use_dis = hr_solver()
hrf_dict_pl_dis = pl_use_dis.compute_hrf_dis(gs_phonon_freqs, gs_phonon_modes,
                                             atomic_symbols, gs_coord, es_coord,
                                             cell_parameters, mass_list=mass_list)
spectral_density_pl_dis = pl_use_dis.compute_spectral_density(hrf_dict_pl_dis,
                                                              energy_axis=np.linspace(0,200,201),
                                                              sigma=[6.0,2.0])

# Absorption HR
abs_use_dis = hr_solver()
hrf_dict_abs_dis = abs_use_dis.compute_hrf_dis(es_phonon_freqs, es_phonon_modes,
                                               atomic_symbols, gs_coord, es_coord,
                                               cell_parameters, mass_list=mass_list)
spectral_density_abs_dis = abs_use_dis.compute_spectral_density(hrf_dict_abs_dis,
                                                                energy_axis=np.linspace(0,200,201),
                                                                sigma=[6.0,2.0])

# Plot
fig, ax = plt.subplots(figsize=(7,4))
ax.plot(spectral_density_abs_dis[0], spectral_density_abs_dis[1], label='Absorption', color=blue)
ax.plot(spectral_density_pl_dis[0], spectral_density_pl_dis[1], label='PL', color=red)
ax.set_xlabel('$\\hbar\\omega$ (meV)')
ax.set_ylabel('Absorption / PL (arb. unit)')
ax.set_xlim([1500,2400])
ax.set_ylim([0,10])
ax.legend()

os.makedirs("images", exist_ok=True)
plt.savefig("images/pl_absorption.png", bbox_inches='tight', dpi=200)
print("Saved pl_absorption.png")
