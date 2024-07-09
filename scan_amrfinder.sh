#!/bin/bash

#SBATCH --job-name=scan
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=2-260
#SBATCH --mem=8GB
#SBATCH --error=scan_amrfinder.err
#SBATCH --output=scan_amrfinder.out

#conda activate amrfinder

acc=$(awk -F '[\t]' '{print $0}' tag_list_qc | awk "NR == $SLURM_ARRAY_TASK_ID")

amrfinder --nucleotide ./${acc}/assembly_unicycler/assembly.fasta -O Salmonella --plus -o ./${acc}/${acc}_amrfinder.tsv

awk -v acc="$acc" 'BEGIN{FS=OFS="\t"} {if(NR==1) $0="Sample_name\t"$0; else $0=acc"\t"$0} 1' ./${acc}/${acc}_amrfinder.tsv > ./${acc}/${acc}_amrfinder_1.tsv

rm -rf ./${acc}/${acc}_amrfinder.tsv

mv ./${acc}/${acc}_amrfinder_1.tsv ./${acc}/${acc}_amrfinder.tsv
