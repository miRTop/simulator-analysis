from os.path import join, splitext, basename
import glob

FORMAT = dict()
FILES = dict()
OUTS = list()
TOOLS = ["star_default", "star_short", "bwa", "gsnap_m3", "gsnap_m5","bowtie1", "shortstack", "manatee", "razers3", "seqbuster"]
# TOOLS = ["star_default", "star_short"]
INPUTDIR = "../../output.round2/gffs"
for t in TOOLS:
    ext  = "bam"
    format = "--genomic --low-memory"
    if t in ["manatee"]:
        ext = "sam"
        format = "--format manatee --genomic --low-memory"
    elif t == "seqbuster":
        ext = "mirna"
        format = "--format  seqbuster"
    elif t == "razers3":
        ext = "bam"
        format = "--low-memory"
    FORMAT[t] = format
    FILES[t] = join("../../output.round2/aligners_output", t,"simulated_miRNA_BF."+ext)

    output = splitext(basename(FILES[t]))[0] + ".gff"
    OUTS.append(join("../../output.round2/gffs", t, output))

print(FILES)
print(OUTS)
print(FORMAT)
# Globals ---------------------------------------------------------------------


# Globals ---------------------------------------------------------------------

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = '../../raw/round2'
# A Snakemake regular expression matching the forward mate FASTQ files.
FASTQS = "simulated_miRNA_BF.fa" 

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
FASTA = "../../reference22/hairpin.fa"
GFF = "../../reference22/hsa.hg38.gff3"


rule all:
    input:
        fixed = '../../reference22/hairpin_atcg_hsa.fa',
        gffs = expand("../../output.round2/gffs/{sample}/simulated_miRNA_BF.gff", sample=TOOLS)

rule reference:
    input:
        fasta = FASTA
    output:
        fixed = '../../reference22/hairpin_atcg_hsa.fa'
    shell:
        '''
        cat {input.fasta} | awk '{{if ($0~/>hsa/){{name=$0; print name}} else if ($0~/^>/){{name=0}};if (name!=0 && $0!~/^>/){{print $0;}}}}' | sed 's/U/T/g' > {output.fixed}
        '''

rule get_gffs:
    input:
        fasta = '../../reference22/hairpin_atcg_hsa.fa',
        gff = GFF,
        bam = lambda wildcards: FILES[wildcards.sample]
    params: 
        options = lambda wildcards: FORMAT[wildcards.sample],
        out = "../../output.round2/gffs/{sample}",
    output: ("../../output.round2/gffs/{sample}/simulated_miRNA_BF.gff")
    shell:
        'mirtop gff --keep-name {params.options} --hairpin {input.fasta} --gtf {input.gff} -o {params.out} {input.bam}'
