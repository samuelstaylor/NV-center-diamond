import os
import matplotlib.pyplot as plt
from pypl.hr_solver import hr_solver
from pypl.utils import parse_phonopy_h5, parse_atoms_qexml, parse_forces_qexml

# Colors
blue = '#4285F4'
red = '#DB4437'

# Load phonons and coordinates
gs_phonon_file = 'phonon/gs/gs-ph-mesh.hdf5'
gs_file = 'relax/pwscf.xml'
es_file = 'tddft/final_geo/pwscf.xml'

gs_phonon_freqs, gs_phonon_modes = parse_phonopy_h5(gs_phonon_file)
atomic_symbols, gs_coord, cell_parameters = parse_atoms_qexml(gs_file)
_, es_coord, _ = parse_atoms_qexml(es_file)

# Load forces file
forces_fname = 'github_repos/NV-center-diamond/2x2x2/phonon/es/pwscf.xml'
atomic_symbols, gs_forces_es_coord = parse_forces_qexml(forces_fname)

mass_list = {'C': 12.0107, 'N': 14.0067}

# Displacement-based PL for comparison
pl_use_dis = hr_solver()
hrf_dict_pl_dis = pl_use_dis.compute_hrf_dis(gs_phonon_freqs, gs_phonon_modes,
                                             atomic_symbols, gs_coord, es_coord,
                                             cell_parameters, mass_list=mass_list)

# Forces-based HR
pl_using_forces = hr_solver()
hrf_dict_pl_forces = pl_using_forces.compute_hrf_forces(gs_phonon_freqs, gs_phonon_modes,
                                                        atomic_symbols, gs_forces_es_coord,
                                                        mass_list)

# Lineshape & spectrum
linshape_fft_pl_forces = pl_using_forces.compute_lineshape_fft(hrf_dict_pl_forces,
                                                               temp=4, sigma=[6,2],
                                                               zpl_broadening=0.3)
ezpl = 1831.592
spectrum_pl_forces = pl_using_forces.compute_spectrum(ezpl, spectrum_type='PL',
                                                      lineshape=linshape_fft_pl_forces)

spectrum_pl_dis = pl_use_dis.compute_spectrum(ezpl, spectrum_type='PL',
                                              lineshape=pl_use_dis.compute_lineshape_fft(hrf_dict_pl_dis,
                                                                                        temp=4, sigma=[6,2],
                                                                                        zpl_broadening=0.3))

# Plot
fig, ax = plt.subplots(figsize=(7,4))
ax.plot(spectrum_pl_dis[0], spectrum_pl_dis[1], label='use displacements', color=blue)
ax.plot(spectrum_pl_forces[0], spectrum_pl_forces[1], label='use forces', linestyle='--', color=red)
ax.set_xlabel('$\\hbar\\omega$ (meV)')
ax.set_ylabel('PL (arb. unit)')
ax.set_xlim([1300,2000])
ax.set_ylim([0,10])
ax.legend()

os.makedirs("images", exist_ok=True)
plt.savefig("images/pl_forces.png", bbox_inches='tight', dpi=200)
print("Saved pl_forces.png")
