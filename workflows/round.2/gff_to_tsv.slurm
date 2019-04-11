#!/bin/sh
#SBATCH -p priority
#SBATCH -J counts
#SBATCH -o gff_to_tsv_logs/run.o
#SBATCH -e gff_to_tsv_logs/run.e
#SBATCH -t 00-23:59:00
#SBATCH -c 1
#SBATCH --mem=16000

snakemake -p --stats gff_to_tsv_logs/run.stats -s gff_to_tsv.snakefile
