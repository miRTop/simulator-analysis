import glob
from os.path import join, splitext, basename

FILES = dict()
OUTS = list()
TOOLS = ["star_default", "star_short", "bwa", "gsnap_m3", "gsnap_m5","bowtie1", "shortstack", "manatee", "razers3", "seqbuster"]
# TOOLS = ["star_default", "star_short"]
INPUTDIR = "../../output.round2/gffs"
for t in TOOLS:
    FILES[t] = join("../../output.round2/gffs", t,"simulated_miRNA_BF.gff")

    output = splitext(basename(FILES[t]))[0] + ".tsv"
    OUTS.append(join("../../output.round2/counts", t, output))

print(FILES)
print(OUTS)
#from os.path import join

# Globals ---------------------------------------------------------------------

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
FASTA = "../../reference22/hairpin.fa"
GFF = "../../reference22/hsa.hg38.gff3"

rule all:
    input:
        fixed = '../../reference22/hairpin_atcg_hsa.fa',
        tsv = expand("../../output.round2/counts/{sample}/simulated_miRNA_BF.tsv", sample=TOOLS)

rule reference:
    input:
        fasta = FASTA
    output:
        fixed = '../../reference22/hairpin_atcg_hsa.fa'
    shell:
        '''
        cat {input.fasta} | awk '{{if ($0~/>hsa/){{name=$0; print name}} else if ($0~/^>/){{name=0}};if (name!=0 && $0!~/^>/){{print $0;}}}}' | sed 's/U/T/g' > {output.fixed}
        '''

rule get_tsvs:
    input:
        fasta = '../../reference22/hairpin_atcg_hsa.fa',
        mirbase = GFF,
        gff = lambda wildcards: FILES[wildcards.sample]
    params: 
        out = "../../output.round2/counts/{sample}"
    output: ("../../output.round2/counts/{sample}/simulated_miRNA_BF.tsv")
    shell:
        'mirtop counts --hairpin {input.fasta} --gtf {input.mirbase} -o {params.out} --gff {input.gff}'


