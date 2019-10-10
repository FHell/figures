#!/bin/bash

#SBATCH --qos=standard
#SBATCH --job-name="ppplot"
#SBATCH --account=coen
#SBATCH --output=ppplot-%j-%N.out
#SBATCH --error=pplot-%j-%N.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=END
#SBATCH --mail-user=hellmann@pik-potsdam.de

echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "------------------------------------------------------------"

module load julia/1.1.0
module load hpc/2015
julia PhaseSpacePlots.jl
