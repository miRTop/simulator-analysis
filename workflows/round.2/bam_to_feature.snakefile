from os.path import join, splitext, basename
import glob

FORMAT = dict()
FILES = dict()
OUTS = list()
TOOLS = ["star_default", "star_short", "bwa", "gsnap_m3", "gsnap_m5","shortstack"]
# TOOLS = ["star_default", "star_short"]
INPUTDIR = "../../output.round2/gffs"
for t in TOOLS:
    FORMAT[t] = format
    FILES[t] = join("../../output.round2/aligners_output", t,"simulated_miRNA_BF.bam")


print(FILES)
# Globals ---------------------------------------------------------------------


# Globals ---------------------------------------------------------------------

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = '../../raw/round2'
# A Snakemake regular expression matching the forward mate FASTQ files.
FASTQS = "simulated_miRNA_BF.fa" 

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
GFF = "../../reference22/hsa.hg38.mature.gff3"


rule all:
    input:
        counts = expand("../../output.round2/features/{sample}.txt", sample=TOOLS),
        counts_m = expand("../../output.round2/features_multiON/{sample}.txt", sample=TOOLS),
        counts_m_f = expand("../../output.round2/features_multiON_fractionON/{sample}.txt", sample=TOOLS)

rule counts:
    input:
        gff = GFF,
        bam = lambda wildcards: FILES[wildcards.sample]
    params: 
        options = lambda wildcards: FORMAT[wildcards.sample]
    conda: 
        "envs/bam_to_feature.yaml"
    output: ("../../output.round2/features/{sample}.txt")
    shell:
        'featureCounts -t miRNA -g Name -o {output} -a {input.gff}  -s 1 {input.bam}'

rule countsM:
    input:
        gff = GFF,
        bam = lambda wildcards: FILES[wildcards.sample]
    params: 
        options = lambda wildcards: FORMAT[wildcards.sample]
    conda: 
        "envs/bam_to_feature.yaml"
    output: ("../../output.round2/features_multiON/{sample}.txt")
    shell:
        'featureCounts -t miRNA -g Name -o {output} -a {input.gff} -M -s 1 {input.bam}'

rule countsMF:
    input:
        gff = GFF,
        bam = lambda wildcards: FILES[wildcards.sample]
    params: 
        options = lambda wildcards: FORMAT[wildcards.sample]
    conda: 
        "envs/bam_to_feature.yaml"
    output: ("../../output.round2/features_multiON_fractionON/{sample}.txt")
    shell:
        'featureCounts -t miRNA -g Name -o {output} -a {input.gff} -M -f  -s 1 {input.bam}'
