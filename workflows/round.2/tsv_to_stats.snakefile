from os.path import join

# Globals ---------------------------------------------------------------------


# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
FASTA = "../../raw/round2/simulated_miRNA_BF.fa"


rule all:
    input:
        bowtie1 = '../../output.round2/stats/bowtie1.tsv',
        shortstack = '../../output.round2/stats/shortstack.tsv',
        manatee = "../../output.round2/stats/manatee.tsv",
        razers3 = '../../output.round2/stats/razers3.tsv'

rule bowtie1:
    input:
        fasta = FASTA,
        tsv = '../../output.round2/counts/bowtie1/bowtie.simulator.miRNA.combinations.tsv'
    output:
        '../../output.round2/stats/bowtie1.tsv'
    shell:
        'python ../../scripts/parse_gff_count.py --counts {input.tsv} --input {input.fasta} > {output}'

rule shortstack:
    input:
        fasta = FASTA,
        tsv = '../../output.round2/counts/shortstack/simulated_miRNA_BF.tsv'
    output:
        '../../output.round2/stats/shortstack.tsv'
    shell:
        'python ../../scripts/parse_gff_count.py --counts {input.tsv} --input {input.fasta} > {output}'

rule manatee:
    input:
        fasta = FASTA,
        tsv = '../../output.round2/counts/manatee/simulated_miRNA_BF_Manatee.tsv'
    output:
        '../../output.round2/stats/manatee.tsv'
    shell:
        'python ../../scripts/parse_gff_count.py --counts {input.tsv} --input {input.fasta} > {output}'

rule razers3:
    input:
        fasta = FASTA,
        tsv = '../../output.round2/counts/razers3/simulated_miRNA_BF.tsv'
    output:
        '../../output.round2/stats/razers3.tsv'
    shell:
        'python ../../scripts/parse_gff_count.py --counts {input.tsv} --input {input.fasta} > {output}'
