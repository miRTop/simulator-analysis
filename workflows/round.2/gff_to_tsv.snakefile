from os.path import join

# Globals ---------------------------------------------------------------------


# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
FASTA = "../../reference22/hairpin.fa"
GFF = "../../reference22/hsa.hg38.gff3"


rule all:
    input:
        fixed = '../../reference22/hairpin_atcg_hsa.fa',
        bowtie1 = '../../output.round2/counts/bowtie1/bowtie.simulator.miRNA.combinations.tsv',
        shortstack = '../../output.round2/counts/shortstack/simulated_miRNA_BF.tsv',
        manatee = '../../output.round2/counts/manatee/simulated_miRNA_BF_Manatee.tsv'

rule reference:
    input:
        fasta = FASTA
    output:
        fixed = '../../reference22/hairpin_atcg_hsa.fa'
    shell:
        '''
        cat {input.fasta} | awk '{{if ($0~/>hsa/){{name=$0; print name}} else if ($0~/^>/){{name=0}};if (name!=0 && $0!~/^>/){{print $0;}}}}' | sed 's/U/T/g' > {output.fixed}
        '''

rule bowtie1:
    input:
        fasta = '../../reference22/hairpin_atcg_hsa.fa',
        mirbase = GFF,
        gff = '../../output.round2/gffs/bowtie1/bowtie.simulator.miRNA.combinations.gff'
    output:
        '../../output.round2/counts/bowtie1/bowtie.simulator.miRNA.combinations.tsv'
    shell:
        'mirtop counts --hairpin {input.fasta} --gtf {input.mirbase} -o ../../output.round2/counts/bowtie1 --gff {input.gff}'

rule shortstack:
    input:
        fasta = '../../reference22/hairpin_atcg_hsa.fa',
        mirbase = GFF,
        gff = '../../output.round2/gffs/shortstack/simulated_miRNA_BF.gff'
    output:
        '../../output.round2/counts/shortstack/simulated_miRNA_BF.tsv'
    shell:
        'mirtop counts --hairpin {input.fasta} --gtf {input.mirbase} -o ../../output.round2/counts/shortstack --gff {input.gff}'

rule manatee:
    input:
        fasta = '../../reference22/hairpin_atcg_hsa.fa',
        mirbase = GFF,
        gff = '../../output.round2/gffs/manatee/simulated_miRNA_BF_Manatee.gff'
    output:
        '../../output.round2/counts/manatee/simulated_miRNA_BF_Manatee.tsv'
    shell:
        'mirtop counts --hairpin {input.fasta} --gtf {input.mirbase} -o ../../output.round2/counts/manatee --gff {input.gff}'
