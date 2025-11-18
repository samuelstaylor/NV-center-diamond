import numpy as np
import matplotlib.pyplot as plt
import os
from pypl.hr_solver import hr_solver
from pypl.utils import parse_phonopy_h5, parse_atoms_qexml

plt.rcParams.update({'font.size': 12})

# Colors
deep_violet = '#9D80B8'
vibrant_purple = '#D391C2'
soft_coral = '#F3B1BA'
ocean_blue = '#81BFD6'
goldenrod_yellow = '#EBD68F'
emerald_green = '#02A650'

# File paths
gs_phonon_file = 'phonon/gs/gs-ph-mesh.hdf5'
gs_file = 'relax/pwscf.xml'
es_file = 'tddft/final_geo/pwscf.xml'

# Parse input
gs_phonon_freqs, gs_phonon_modes = parse_phonopy_h5(gs_phonon_file)
atomic_symbols, gs_coord, cell_parameters = parse_atoms_qexml(gs_file)
_, es_coord, _ = parse_atoms_qexml(es_file)

mass_list = {'C': 12.0107, 'N': 14.0067}

pl_use_dis = hr_solver()
hrf_dict_pl_dis = pl_use_dis.compute_hrf_dis(
    gs_phonon_freqs, gs_phonon_modes, atomic_symbols,
    gs_coord, es_coord, cell_parameters, mass_list=mass_list
)

# Spectral densities for multiple sigmas
sigmas = [[3,2],[3,4],[3,6],[6,2],[6,4],[6,6]]
spectral_density_pl_dis = []
for sigma in sigmas:
    energy_axis = np.linspace(0, 200, 201)
    spectral_density_pl_dis.append(
        pl_use_dis.compute_spectral_density(hrf_dict_pl_dis, energy_axis=energy_axis, sigma=sigma)
    )

colors = [deep_violet, vibrant_purple, soft_coral, ocean_blue, goldenrod_yellow, emerald_green]
fig, ax = plt.subplots(figsize=(7,4))

# Collect max value to adjust y-limit later
max_val = 0.0

for i in range(len(spectral_density_pl_dis)):
    axis, values = spectral_density_pl_dis[i]
    ax.plot(axis, values, color=colors[i], label=r'$\sigma=[%.1f, %.1f]$' % tuple(sigmas[i]))
    max_val = max(max_val, np.max(values))  # track the highest y value

ax.set_xlabel('$\\hbar\\omega$ (meV)')
ax.set_ylabel('$S(\\hbar\\omega)$ (1/meV)')
ax.set_xlim(0,200)
ax.set_ylim(0, max_val * 1.05)  # set y-limit slightly above max
ax.legend()
os.makedirs("images", exist_ok=True)
plt.savefig("images/spectral_density.png", bbox_inches='tight', dpi=200)

