#!/bin/bash

#SBATCH --job-name=mapping
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=260
#SBATCH --mem=64GB
#SBATCH --error=mapping.err
#SBATCH --output=mapping.out

#conda activate snippy

run_acc=$(awk -F '[\t]' '{print $0}' ./tag_list_qc | awk "NR == $SLURM_ARRAY_TASK_ID")

mkdir ./${run_acc}/mapping_K_099
snippy --cpus 16 --force --outdir ./${run_acc}/mapping_K_099 --ref ./K_099/pseudogenome.fasta --R1 ./${run_acc}/${run_acc}_1_clean.fastq.gz --R2 ./${run_acc}/${run_acc}_2_clean.fastq.gz
rm -rf ./${run_acc}/mapping_K_099/snps.bam
