#!/bin/sh
#SBATCH -p priority
#SBATCH -J gff
#SBATCH -o run.o
#SBATCH -e run.e
#SBATCH -t 00-23:59:00
#SBATCH -c 1
#SBATCH --mem=16000

snakemake -p --stats run2.stats -s bam_to_gff.snakefile
