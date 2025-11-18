import numpy as np
import matplotlib.pyplot as plt
import os
from pypl.hr_solver import hr_solver
from pypl.utils import parse_phonopy_h5, parse_atoms_qexml

plt.rcParams.update({'font.size': 12})
blue = '#4285F4'

gs_phonon_file = 'phonon/gs/gs-ph-mesh.hdf5'
gs_file = 'relax/pwscf.xml'
es_file = 'tddft/final_geo/pwscf.xml'


gs_phonon_freqs, gs_phonon_modes = parse_phonopy_h5(gs_phonon_file)
atomic_symbols, gs_coord, cell_parameters = parse_atoms_qexml(gs_file)
_, es_coord, _ = parse_atoms_qexml(es_file)
mass_list = {'C':12.0107,'N':14.0067}

pl_use_dis = hr_solver()
hrf_dict_pl_dis = pl_use_dis.compute_hrf_dis(gs_phonon_freqs, gs_phonon_modes,
                                             atomic_symbols, gs_coord, es_coord,
                                             cell_parameters, mass_list=mass_list)

lineshape_fft_pl_dis = pl_use_dis.compute_lineshape_fft(hrf_dict_pl_dis,
                                                       temp=4,
                                                       sigma=[6,2],
                                                       zpl_broadening=0.3)

fig, ax = plt.subplots(figsize=(7,4))
ax.plot(lineshape_fft_pl_dis[0], lineshape_fft_pl_dis[1], color=blue)
ax.set_xlabel('$\\hbar\\omega$ (meV)')
ax.set_ylabel('$A(\\hbar\\omega)$')
ax.set_xlim([-100,800])
ax.set_ylim([0,12])
os.makedirs("images", exist_ok=True)
plt.savefig("images/lineshape_function.png", bbox_inches='tight', dpi=200)
