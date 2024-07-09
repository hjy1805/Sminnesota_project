#!/bin/bash

#SBATCH --job-name=srst
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=1-260
#SBATCH --mem=64GB
#SBATCH --error=typing_srst.err
#SBATCH --output=typing_srst.out

#conda activate srst

acc=$(awk -F '[\t]' '{print $0}' tag_list_qc | awk "NR == $SLURM_ARRAY_TASK_ID")

cd ./${acc}

srst2 --input_pe ${acc}_1_clean.fastq.gz ${acc}_2_clean.fastq.gz --output ${acc} --forward _1_clean --reverse _2_clean --log --gene_db /ibex/project/c2205/databases/srst_db/VFDB_am_clean.fasta --threads 32 --min_coverage 60

#srst2 --input_pe ${acc}_1_clean.fastq.gz ${acc}_2_clean.fastq.gz --output ${acc} --forward _1_clean --reverse _2_clean --log --gene_db /ibex/project/c2205/databases/srst_db/CARD.fasta --threads 16

#srst2 --input_pe ${acc}_1_clean.fastq.gz ${acc}_2_clean.fastq.gz --output ${acc} --forward _1_clean --reverse _2_clean --log --gene_db /ibex/project/c2205/databases/srst_db/plasmidfinder.fasta --threads 16

#srst2 --input_pe ${acc}_1_clean.fastq.gz ${acc}_2_clean.fastq.gz --output ${acc} --forward _1_clean --reverse _2_clean --log --gene_db /ibex/project/c2205/databases/srst_db/ResFinder.fasta --threads 16
