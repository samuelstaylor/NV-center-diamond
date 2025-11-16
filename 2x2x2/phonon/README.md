Steps of phonon calculations using Phonopy and Quantum ESPRESSO

1. Generate displacements (from the ground state structure)
```
    phonopy --qe -d --dim="1 1 1" -c gs-dft-pw.in -v
```

2. With the `generated supercell-***.in` files, append the header for each scf calculation. Use script `header_merge.py` to do this
```
    python header_merge.py
```

3. Perform DFT calculations using the displaced structures. Do a "scf" calculation on every displaced structure (all `supercell-***.in` files). Computes the forces on each atom
```
    sbatch run_all_jobs.sh
```

4. Collect forces. creates the file: `FORCE_SETS`
```
    phonopy --qe -f supercell-{001..078}.out
```

5. Compute phonon frequencies and modes. Creates the file: `mesh.hdf5`. Rename it to `gs_ph_mesh.hdf5`
```
    phonopy --dim="1 1 1" --fc-symmetry --mesh="1 1 1" --eigenvectors --writefc --qe -c gs-dft-pw.in --mesh-format=hdf5
```

6. OPTIONAL: IF YOU WANT ABSORPTION SPECTRUM (vertical excitation):
- Repeat the same states but use the excited state geometry (calculated from tddft or dscf).
    - follow steps in `phonon/es/` directory.
