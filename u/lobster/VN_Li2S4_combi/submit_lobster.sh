#!/bin/bash
#SBATCH --job-name=LOB_VN_Li2S4_combi
#SBATCH --output=LOB_VN_Li2S4_combi_%j.out
#SBATCH --error=LOB_VN_Li2S4_combi_%j.err
#SBATCH --account=def-peslherb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=50G
#SBATCH --time=43:00:00

source ~/ase/bin/activate
module load quantumespresso
set -euo pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd /scratch/anizami/QE_2/u/lobster/VN_Li2S4_combi
rm -f lobster.out
/home/anizami/lobster/lobster-5.1.1 > lobster.out 2>&1
echo "$(date): VN_Li2S4_combi COMPLETED" >> /scratch/anizami/QE_2/u/lobster/automated_lobster_submissions.txt