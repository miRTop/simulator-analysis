from os.path import join

# Globals ---------------------------------------------------------------------

# Full path to a folder that holds all of your FASTQ files.
IDIR = join('../../simulated_data/', config['round'])
ADIR = snakemake.utils.format('../../output.{round}', round = config['round'])
# A Snakemake regular expression matching the forward mate FASTQ files.
SAMPLES = ["simulated"]
FASTQ = expand(join(IDIR, '{sample}.fa'), sample = SAMPLES)

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
FASTA = "../../reference22/hairpin.fa"
GFF = "../../reference22/hsa.gff3"

rule all:
    input:
        fixed = '../../reference22/hairpin_atcg_hsa.fa',
        collapsed = expand('{idir}/{sample}_trimmed.fastq', sample = SAMPLES, idir = IDIR),
        seqbuster_quant = expand('{adir}/aligners_output/seqbuster-quant/{sample}_trimmed.mirna', sample = SAMPLES, adir = ADIR),
        seqbuster = expand('{adir}/aligners_output/seqbuster/{sample}.mirna', sample = SAMPLES, adir = ADIR),
        srnabench = expand('{adir}/aligners_output/srnabench/{sample}.grouped', sample = SAMPLES, adir = ADIR),
        sam = expand(join(ADIR, 'aligners_output/razers3-pre', '{sample}.bam'), sample = SAMPLES),
        bwa = expand(join(ADIR, 'aligners_output/bwa-pre', '{sample}.bam'), sample = SAMPLES) 


rule reference:
    input:
        fasta = FASTA
    output:
        fixed = '../../reference22/hairpin_atcg_hsa.fa'
    shell:
        '''
        cat {input.fasta} | awk '{{if ($0~/>hsa/){{name=$0; print name}} else if ($0~/^>/){{name=0}};if (name!=0 && $0!~/^>/){{print $0;}}}} | sed 's/U/T/g' > {output.fixed}
        '''

rule razers3:
    threads: 4
    resources:
      mem_mb = 24000
    input:
        fasta = '../../reference22/hairpin_atcg_hsa.fa',
        reads = FASTQ
    output:
        join(ADIR, 'aligners_output/razers3-pre/{sample}.bam')
    params:
        dir = join(ADIR, 'aligners_output/razers3-pre')
    run:
        if (config['norazers3']):
            print("Skip razers3 with full file")
            shell('mkdir -p {params.dir}')
            shell('touch {output}')
        else:
            shell('mkdir -p {params.dir}')
            shell('razers3 -tc {threads} -f -o tmp.sam {input.fasta} {input.reads}')
            shell('samtools view -Sb tmp.sam > {output} && rm tmp.sam')

rule bwa:
    threads: 4
    resources:
        mem_mb = 24000
    input:
        fasta = '../../reference22/hairpin_atcg_hsa.fa',
        reads = FASTQ
    output:
        join(ADIR, 'aligners_output/bwa-pre/{sample}.bam')
    params:
        dir = join(ADIR, 'aligners_output/bwa-pre')
    shell:
        '''
        mkdir -p {params.dir}
        mkdir -p index
        bwa index -p index/pre -a is {input.fasta}
        bwa aln -t {threads} index/pre {input.reads} > tmp.sai
        bwa samse -f tmp.sam index/pre tmp.sai {input.reads}
        samtools view -Sb tmp.sam > {output}
        rm tmp.*
        '''


rule collapse:
    resources:
        mem_mb = 24000
    input:
        fasta = FASTQ
    output:
        join(IDIR, '{sample}_trimmed.fastq')
    params:
        dir = IDIR
    shell:
        '''
        mkdir -p {params.dir}
        seqcluster collapse -f {input.fasta} -m 0 --min_size 10 -o {params.dir}
        '''

rule seqbuster_quant:
    resources:
        mem_mb = 48000
    input:
        db = '../../reference22',
        reads = join(IDIR, '{sample}_trimmed.fastq')
    output:
        out = join(ADIR, "aligners_output/seqbuster-quant/{sample}_trimmed.mirna")
    params:
        dir = join(ADIR, "aligners_output/seqbuster-quant/{sample}_trimmed")
    shell:
        'mkdir -p {params.dir} && '
        'miraligner -Xms92m -Xmx{resources.mem_mb}m -sub 1 -trim 3 -add 3 -minl 16 -s hsa -i {input.reads} -db {input.db} -o {params.dir}'

rule seqbuster:
    resources:
        mem_mb = 124000
    input:
        db = '../../reference22',
        reads = FASTQ
    output:
        out = join(ADIR, "aligners_output/seqbuster/{sample}.mirna")
    params:
        dir = join(ADIR, "aligners_output/seqbuster/{sample}")
    run:
        if (config['noseqbusterfull']):
            print("Skip seqbuster with full file")
            shell('mkdir -p {params.dir}')
            shell('echo done > {output.out}')
        else:
            shell('mkdir -p {params.dir}')
            shell('miraligner -Xms92m -Xmx{resources.mem_mb}m -sub 1 -trim 3 -add 3 -minl 16 -s hsa -i {input.reads} -db {input.db} -o {params.dir}')

rule srnabench:
    resources:
        mem_mb = 48000
    input:
        db = "../../tools/sRNAtoolboxDB",
        reads = join(IDIR, '{sample}_trimmed.fastq')
    output:
        out = join(ADIR, "aligners_output/srnabench/{sample}.grouped")
    params:
        dir = join(ADIR, "aligners_output/srnabench")
    shell:
        '''
        mkdir -p {params.dir}
        sed  's/_x/#/g' {input.reads}  | sed  's/@/>/g' | sed  's/_//g' | grep -v IIII | grep -v '+' > {input.reads}.fasta
        java  -Xms92m -Xmx{resources.mem_mb}m -jar ../../tools/sRNAtoolboxDB/exec/sRNAbench.jar dbPath={input.db} microRNA=hsa input={input.reads}.fasta output={params.dir}
        mv {params.dir}/mature_sense.grouped {output.out}
        '''
