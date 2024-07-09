#!/bin/bash

#SBATCH --job-name=beast
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --array=2-4
#SBATCH --mem=128GB
#SBATCH --error=beast_%j.err
#SBATCH --output=beast_%j.out

#module load java

baps=$(awk -F '[\t]' '{print $0}' ./baps_list | awk "NR == $SLURM_ARRAY_TASK_ID")

beast -java -working -threads 64 ./${baps}/beast/${baps}_city_mascot_1.xml
