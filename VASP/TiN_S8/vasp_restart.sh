#!/bin/bash
#SBATCH --account=def-aragab
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3000M
#SBATCH --time=2-00:00:00
#SBATCH --job-name=TiN_S8_rst

module purge
module load StdEnv/2023 intel/2023.2.1 intelmpi/2021.9.0 imkl/2023.2.0 python/3.10

source ~/env/bin/activate

export OMP_NUM_THREADS=1

export ASE_VASP_COMMAND="srun /home/hoansar/vasp.6.4.2/bin/vasp_std"
export VASP_PP_PATH="/home/hoansar/potentials"

python TiN_S8_restart.py
