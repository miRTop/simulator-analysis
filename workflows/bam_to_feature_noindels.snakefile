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

for t in TOOLS:
    FORMAT[t] = format
    FILES[t] = join(ADIR, "aligners_output", t, "simulated.bam")

ADIR = join(ADIR, "noindels")
if not os.path.exists(ADIR):
    mkdir(ADIR)

for folder in ['features', 'features_multiON']:
    if not os.path.exists(join(ADIR, folder)):
        mkdir(join(ADIR, folder))

print(FILES)
# Globals ---------------------------------------------------------------------

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
GFF = "../../reference22/hsa.hg38.mature.gff3"

rule all:
    input:
        counts = expand("{dir}/features/{sample}.txt", sample=TOOLS, dir = ADIR),
        counts_m = expand("{dir}/features_multiON/{sample}.txt", sample=TOOLS, dir = ADIR),
        bam = expand("{dir}/{sample}.bam", sample=TOOLS, dir = ADIR)

rule clean:
    input:
        bam = lambda wildcards: FILES[wildcards.sample]
    output: join(ADIR, "{sample}.bam")
    shell:
        'samtools view -h {input.bam} | grep -v -i indel | samtools view -bS - >{output}'

rule counts:
    input:
        gff = GFF,
        bam = join(ADIR, "{sample}.bam")
    output: join(ADIR, "features/{sample}.txt")
    shell:
        'featureCounts -t miRNA -g Name -o {output} -a {input.gff}  -s 1 {input.bam}'

rule countsM:
    input:
        gff = GFF,
        bam = join(ADIR, "{sample}.bam")
    output: join(ADIR, "features_multiON/{sample}.txt")
    shell:
        'featureCounts -t miRNA -g Name -o {output} -a {input.gff} -M -s 1 {input.bam}'

