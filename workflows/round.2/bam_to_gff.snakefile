from os.path import join

# Globals ---------------------------------------------------------------------

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = '../../raw/round2'
# A Snakemake regular expression matching the forward mate FASTQ files.
SAMPLES = ["simulated_miRNA_BF"]
FASTQS = expand(join(FASTQ_DIR, '{sample}.fa'), sample = SAMPLES)

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
FASTA = "../../reference22/hairpin.fa"
GFF = "../../reference22/hsa.hg38.gff3"


rule all:
    input:
        fixed = '../../reference22/hairpin_atcg_hsa.fa',
        bowtie2 = '../../output.round2/gffs/bowtie1/bowtie.simulator.miRNA.combinations.gff',
        shortstack = '../../output.round2/gffs/shortstack/simulated_miRNA_BF.gff'

rule reference:
    input:
        fasta = FASTA
    output:
        fixed = '../../reference22/hairpin_atcg_hsa.fa'
    shell:
        '''
        cat {input.fasta} | awk '{{if ($0~/>hsa/){{name=$0; print name}} else if ($0~/^>/){{name=0}};if (name!=0 && $0!~/^>/){{print $0;}}}}' | sed 's/U/T/g' > {output.fixed}
        '''

rule bowtie2:
    input:
        fasta = '../../reference22/hairpin_atcg_hsa.fa',
        gff = GFF
    output:
        '../../output.round2/gffs/bowtie1/bowtie.simulator.miRNA.combinations.gff'
    shell:
        'mirtop gff --genomic --keep-name --low-memory --hairpin {input.fasta} --gtf {input.gff} -o ../../output.round2/gffs/bowtie1 ../../output.round2/aligners_output/Bowtie1/bowtie.simulator.miRNA.combinations.bam'

rule shortstack:
    input:
        fasta = '../../reference22/hairpin_atcg_hsa.fa',
        gff = GFF
    output:
        '../../output.round2/gffs/shortstack/simulated_miRNA_BF.gff'
    shell:
        'mirtop gff --genomic --keep-name --low-memory --hairpin {input.fasta} --gtf {input.gff} -o ../../output.round2/gffs/shortstack ../../output.round2/aligners_output/ShortStack/simulated_miRNA_BF.bam'

# 23h 64G RAM
rule manatee:
    input:
        fasta = '../../reference22/hairpin_atcg_hsa.fa',
        gff = GFF
    output:
        '../../output.round2/gffs/manatee/simulated_miRNA_BF_Manatee.gff'
    shell:
        'mirtop gff --format manatee --genomic --keep-name --low-memory --hairpin {input.fasta} --gtf {input.gff} -o ../../output.round2/gffs/manatee ../../output.round2/aligners_output/Manatee/simulated_miRNA_BF_Manatee.sam'
