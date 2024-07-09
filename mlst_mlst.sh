#!/bin/bash

#SBATCH --job-name=mlst
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-239
#SBATCH --mem=8GB
#SBATCH --error=mlst.err
#SBATCH --output=mlst.out

source activate mlst

run_acc=$(awk -F '[\t]' '{print $0}' ./tag_list_qc | awk "NR == $SLURM_ARRAY_TASK_ID")

mlst ./${run_acc}/assembly_unicycler/assembly.fasta > ./${run_acc}/mlst_type.tsv
