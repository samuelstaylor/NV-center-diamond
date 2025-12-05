#!/bin/bash

# Loop through all directory names that match Image-*
for img in Image-*; do
    if [ -d "$img" ]; then
        echo "Submitting job for $img"

        # Create a temporary job script inside the image directory
        cat <<EOF > $img/job_submit.sh
#!/bin/bash
#SBATCH --job-name=pl_es_${img}
#SBATCH --time=24:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-gagalli
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --output=$img/slurm_%j.out
#SBATCH --error=$img/slurm_%j.err

module load intel/19.1.1
module load intelmpi/2019.up7+intel-19.1.1
module load mkl/2020.up1
module load python/anaconda-2020.11

export LD_LIBRARY_PATH=\$PYTHON_DIR/lib:\$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1

ulimit -s unlimited

cd $PWD/$img
mpirun -np 96 pw.x -in pw.in > pw.out
EOF

        # Submit the job
        sbatch $img/job_submit.sh
    fi
done
