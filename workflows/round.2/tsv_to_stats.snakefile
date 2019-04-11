from os.path import join
import glob

FILES = dict()
OUTS = list()
TOOLS = ["star_default", "star_short", "bwa", "gsnap_m3", "gsnap_m5","bowtie1", "shortstack", "manatee", "razers3", "seqbuster"]
INPUTDIR = "../../output.round2/stats"
for t in TOOLS:
    FILES[t] = glob.glob(join("../../output.round2/counts", t, "*.tsv"))[0]
    OUTS.append(join("../../output.round2/stats", t + "_accurcy.tsv"))
(SAMPLES,FILENAMES) = glob_wildcards("../../output.round2/counts/{sample}/{filename}.tsv")
OUTS = expand("../../output.round2/stats/{sample}_accuracy.tsv", sample=TOOLS)
print(FILES)
print(OUTS)
print(SAMPLES)
print(FILENAMES)
# Globals ---------------------------------------------------------------------


# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
FASTA = "../../raw/round2/simulated_miRNA_BF.fa"

rule all:
    input: expand("../../output.round2/stats/{sample}_accuracy.tsv", sample=TOOLS)


rule getstats:
    input: 
        tsv = lambda wildcards: FILES[wildcards.sample],
        fasta = FASTA
    output: "../../output.round2/stats/{sample}_accuracy.tsv"
    params:
        prefix = "{sample}"
    shell:
        'python ../../scripts/parse_gff_count.py --counts {input.tsv} --input {input.fasta} --output ../../output.round2/stats/{params.prefix}'

