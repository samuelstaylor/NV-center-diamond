from westpy import *

# note, you don't want to do 96 since that is too many mpi processes for the size of the FFT grid
run_pw = 'mpirun -np 48 pw.x -in pw.in'
run_wbse = 'mpirun -np 48 wbse.x -in wbse.in'
pp = [f'{element}_ONCV_PBE-1.2.upf' for element in ['C', 'N']]

bfgs = bfgs_iter(run_pw=run_pw, run_wbse=run_wbse, pp=pp)
bfgs.solve()
