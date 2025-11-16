#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-gagalli
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1

module load intel/19.1.1
module load intelmpi/2019.up7+intel-19.1.1
module load mkl/2020.up1
module load python/anaconda-2020.11

export LD_LIBRARY_PATH=$PYTHON_DIR/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1

ulimit -s unlimited

# Loop over all supercell*.in files
for infile in supercell-*.in; do
    base=$(basename "$infile" .in)
    echo "Running $infile ..."
    mpirun -np 96 pw.x -in "$infile" > "${base}.out"
done
