calculating the dft gs forces at the es coordinates.

This is an approximation for the ES phonon frequencies 
see note below:
```
True excited-state phonons would require:

For each small displacement of the atoms around the ES geometry, you’d need a full TDDFT (or CDFT) calculation to get the forces from the excited electronic structure.

This is extremely expensive because TDDFT is much more costly than regular DFT, especially if you have many phonon displacements.

Vertical (or “GS-forces-at-ES-geometry”) approximation:

You keep the atomic positions from the ES geometry, but compute forces using ground-state DFT.

This gives an approximate Hessian (second derivative matrix) for the phonons.

Result:

You get an estimate of the ES phonon frequencies that’s much cheaper.

It’s a good approximation if the excited-state forces are not drastically different from the ground-state forces at that geometry.

The accuracy depends on how much the electronic structure changes between GS and ES.

✅ So yes, the “gs_forces_at_es_coord” folder is doing exactly this: using the vertical approximation to get approximate ES phonons.
```

the informationon the final TDDFT positions is in:
NV-center-diamond/2x2x2/tddft/final_geo/pw.in

^ replaced the atom positions from gs-dft-pw.in with this ^

follow the steps outlined in steps