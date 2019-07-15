#!/bin/bash
#SBATCH -J main_simtrial
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Number of CPUs; if 1, then requested cores will be on the same machine
#SBATCH -t 1000 # Runtime in minutes
#SBATCH -p murphy # Partition to submit to; alternative is murphy_secure
#SBATCH --mem=16000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o output/main_simtrial_%A_%a.out # Standard out goes to this file; %A is job name, %a is array_id for batch jobs
#SBATCH -e output/main_simtrial_%A_%a.err # Standard err goes to this filehostname
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=qiantianchen@fas.harvard.edu # Edit for your own email!
module load gcc/8.2.0-fasrc01 openmpi/3.1.1-fasrc01 R/3.6.1-fasrc01 # load R
export R_LIBS_USER=$HOME/apps/R_3.6.1:$R_LIBS_USER
R CMD BATCH /n/home00/tqian/gst_tmle/main_simtrial.R /n/home00/tqian/gst_tmle/Rout/main_simtrial.Rout.${SLURM_ARRAY_TASK_ID}
