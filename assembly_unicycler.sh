#!/bin/bash

#SBATCH --job-name=unicycler
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=1-222
#SBATCH --mem=64GB
#SBATCH --error=assembly_unicycler.err
#SBATCH --output=assembly_unicycler.out

source activate unicycler

run_acc=$(awk -F '[\t]' '{print $0}' ./tag_list_qc | awk "NR == $SLURM_ARRAY_TASK_ID")

mkdir ./${run_acc}/assembly_unicycler
unicycler -1 ./${run_acc}/${run_acc}_1_clean.fastq.gz -2 ./${run_acc}/${run_acc}_2_clean.fastq.gz -o ./${run_acc}/assembly_unicycler -t 32
