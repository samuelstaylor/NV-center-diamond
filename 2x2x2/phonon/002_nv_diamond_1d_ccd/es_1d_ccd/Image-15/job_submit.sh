#!/bin/bash
#SBATCH --job-name=pl_es_Image-15
#SBATCH --time=24:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-gagalli
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --output=Image-15/slurm_%j.out
#SBATCH --error=Image-15/slurm_%j.err

module load intel/19.1.1
module load intelmpi/2019.up7+intel-19.1.1
module load mkl/2020.up1
module load python/anaconda-2020.11

export LD_LIBRARY_PATH=$PYTHON_DIR/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1

ulimit -s unlimited

cd /home/sstaylor/github_repos/NV-center-diamond/2x2x2/phonon/002_nv_diamond_1d_ccd/es_1d_ccd/Image-15
mpirun -np 96 pw.x -in pw.in > pw.out
