#!/bin/sh
#SBATCH -p medium
#SBATCH -J stats
#SBATCH -o tsv_to_stats_logs/run.o
#SBATCH -e tsv_to_stats_logs/run.e
#SBATCH -t 00-23:59:00
#SBATCH -c 1
#SBATCH --mem=94000
snakemake -s tsv_to_stats.snakefile --unlock 
snakemake -p --stats tsv_to_stats_logs/run.stats -s tsv_to_stats.snakefile --rerun-incomplete
