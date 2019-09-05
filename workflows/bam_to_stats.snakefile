from os.path import join, splitext, basename
from os import mkdir
import glob

FORMAT = dict()
FILES = dict()
OUTS = list()
TOOLS = config["tools"]
ADIR = snakemake.utils.format('../../output.{round}', round = config['round'])
print(ADIR)
for folder in ['gffs', 'counts', 'stats']:
    if not os.path.exists(join(ADIR, folder)):
        mkdir(join(ADIR, folder))

for t in TOOLS:
    ext = "bam"
    format = "--genomic --low-memory"
    if t in ["Manatee"]:
        ext = "sam"
        format = "--format manatee --genomic"
    elif t == "seqbuster":
        ext = "mirna"
        format = "--format  seqbuster"
    elif t == "seqbuster-quant":
        ext = "mirna"
        format = "--format  seqbuster"
    elif t == "razers3-pre":
        ext = "bam"
        format = "--low-memory"
    elif t == "bwa-pre":
        ext = "bam"
        format = "--low-memory"
    FORMAT[t] = format
    FILES[t] = join(ADIR, "aligners_output", t, "simulated." + ext)


print(FILES)
print(FORMAT)
# Globals ---------------------------------------------------------------------


# Globals ---------------------------------------------------------------------

# Full path to a folder that holds all of your FASTQ files.
IDIR = join('../../simulated_data/', config['round'])
SAMPLES = ["simulated"]
FASTQ = expand(join(IDIR, '{sample}.fa'), sample = SAMPLES)

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
FASTA = "../../reference22/hairpin.fa"
GFF = "../../reference22/hsa.hg38.gff3"

rule all:
    input:
        fixed = '../../reference22/hairpin_atcg_hsa.fa',
        gffs = expand("{dir}/gffs/{sample}/simulated.gff", sample=TOOLS, dir = ADIR),
        tsv = expand("{dir}/counts/{sample}/simulated.tsv", sample=TOOLS, dir = ADIR),
        stat = expand("{dir}/stats/{sample}_accuracy.tsv", sample=TOOLS, dir = ADIR)

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
        out = join(ADIR, "gffs/{sample}"),
    output: join(ADIR, "gffs/{sample}/simulated.gff")
    shell:
        'mirtop gff --keep-name {params.options} --hairpin {input.fasta} --gtf {input.gff} -o {params.out} {input.bam}'

rule get_tsvs:
    input:
        fasta = '../../reference22/hairpin_atcg_hsa.fa',
        mirbase = GFF,
        gff = join(ADIR, "gffs/{sample}/simulated.gff")
    params: 
        out = join(ADIR, "counts/{sample}")
    output: join(ADIR, "counts/{sample}/simulated.tsv")
    shell:
        'mirtop counts --hairpin {input.fasta} --gtf {input.mirbase} -o {params.out} --gff {input.gff}'

rule getstats:
    input: 
        tsv = join(ADIR, "counts/{sample}/simulated.tsv"),
        fasta = FASTQ
    output: join(ADIR, "stats/{sample}_accuracy.tsv")
    params:
        dir = join(ADIR, "stats"),
        prefix = join(ADIR, "stats/{sample}")
    run:
        if (config['noquant']):
            print("Skip parser")
            shell('mkdir -p {params.dir}')
            shell('touch {output}')
        else:
            shell('python ../../scripts/parse_gff_count.py --counts {input.tsv} --input {input.fasta} --output {params.prefix}')


