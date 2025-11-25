NOTE:

- notice that the calculation type is set to `calculation='relax'`
- this means that we want to do an adiabatic $\Delta$ SCF calculation. 
- this will find us the excited-state equilibrium geometry (minimize total energy for the fixed occupation/excited configuration)
- I noticed that previous calculations failed to converge. key notes
    - `estimated scf < 0.0342 Ry` did not get smaller across iterations
    - Total energy values mostly oscillate instead of monotonically decreasing
    - c_bands: eigenvalues not converged signs
Text book signs of SCF oscillation/instability caused by overly large mixing steps (mixing_beta too large)
- inappropriate mixing algorithm for a hard problem (charged + spin-polarized ΔSCF with fixed occupations)
- so for this calculation, mixing_beta is decreased from 0.7 to 0.3

Reducing β will slow each step's progress (smaller density step) but dramatically improves stability and usually reduces total iterations required to get to converged density for hard problems.