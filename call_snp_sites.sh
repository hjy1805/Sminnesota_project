#!/bin/bash

#SBATCH --job-name=call
#SBATCH --time=60:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1
#SBATCH --mem=64GB
#SBATCH --error=call.err
#SBATCH --output=call.out

#conda activate snpsites

snp-sites -m -v -p -o sminnesota_final_all_clincal sminnesota_final_all_clincal.aln
