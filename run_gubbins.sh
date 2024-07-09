#!/bin/bash

#SBATCH --job-name=gubbins
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=2
#SBATCH --mem=8GB
#SBATCH --error=gubbins.err
#SBATCH --output=gubbins.out

#conda activate ~/miniconda3/envs/gubbins

baps=$(awk -F '[\t]' '{print $0}' ./baps_list | awk "NR == $SLURM_ARRAY_TASK_ID")

cd ./${baps}

run_gubbins.py --prefix ${baps}_clincal ${baps}_clincal.aln
