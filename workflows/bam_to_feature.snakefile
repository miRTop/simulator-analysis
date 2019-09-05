from os.path import join, splitext, basename
from os import mkdir
import glob

FORMAT = dict()
FILES = dict()
OUTS = list()
TOOLS = config['tools_features']
ADIR = snakemake.utils.format("../../output.{round}", round = config['round'])
IDIR = join('../../raw/', config['round'])
print(ADIR)

for folder in ['features', 'features_multiON', 'features_multiON_fractionON']:
    if not os.path.exists(join(ADIR, folder)):
        mkdir(join(ADIR, folder))

for t in TOOLS:
    FORMAT[t] = format
    FILES[t] = join(ADIR, "aligners_output", t, "simulated.bam")


print(FILES)
# Globals ---------------------------------------------------------------------

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
GFF = "../../reference22/hsa.hg38.mature.gff3"

rule all:
    input:
        counts = expand("{dir}/features/{sample}.txt", sample=TOOLS, dir = ADIR),
        counts_m = expand("{dir}/features_multiON/{sample}.txt", sample=TOOLS, dir = ADIR),
        counts_m_f = expand("{dir}/features_multiON_fractionON/{sample}.txt", sample=TOOLS, dir = ADIR)

rule counts:
    input:
        gff = GFF,
        bam = lambda wildcards: FILES[wildcards.sample]
    output: join(ADIR, "features/{sample}.txt")
    shell:
        'featureCounts -t miRNA -g Name -o {output} -a {input.gff}  -s 1 {input.bam}'

rule countsM:
    input:
        gff = GFF,
        bam = lambda wildcards: FILES[wildcards.sample]
    output: join(ADIR, "features_multiON/{sample}.txt")
    shell:
        'featureCounts -t miRNA -g Name -o {output} -a {input.gff} -M -s 1 {input.bam}'

rule countsMF:
    input:
        gff = GFF,
        bam = lambda wildcards: FILES[wildcards.sample]
    output: join(ADIR, "features_multiON_fractionON/{sample}.txt")
    shell:
        'featureCounts -t miRNA -g Name -o {output} -a {input.gff} -M -f  -s 1 {input.bam}'
