#!/bin/sh
#SBATCH -p priority
#SBATCH -J features
#SBATCH -o run.o
#SBATCH -e run.e
#SBATCH -t 00-73:59:00
#SBATCH -c 1
#SBATCH --mem=24000

snakemake -p --use-conda -s bam_to_feature.snakefile
