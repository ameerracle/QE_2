#!/bin/bash
#SBATCH --job-name=pdos_gen
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --account=def-peslherb

Set number of MPI processes based on SLURM allocation
export OMP_NUM_THREADS=1
NP=$SLURM_NTASKS
module load quantumespresso
Pass NP dynamically to the Python script

python3 prjwfc.py --prefix VN_Li2S4 --np $NP --run
