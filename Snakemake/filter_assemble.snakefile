rule all:
    input:
        "fastp/788707_20180129_S_R1.fastq.gz",
        "fastp/788707_20180129_S_R2.fastq.gz",
        "788707_20180129.bam.bai",
        "human/788707_20180129_S_R1.fastq.gz",
        "not_human/788707_20180129_S_R1.fastq.gz",

rule run_fastp:
    input:
        r1 = "788707_20180129_S_R1.fastq.gz",
        r2 = "788707_20180129_S_R2.fastq.gz",
        ad = "IlluminaAdapters.fa"
    output:
        r1 = "fastp/788707_20180129_S_R1.fastq.gz",
        r2 = "fastp/788707_20180129_S_R2.fastq.gz"
    shell:
        """
        fastp -n 1 -l 100 -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --adapter_fasta {input.ad}
        """

rule run_minimap:
    input:
        r1 = "fastp/788707_20180129_S_R1.fastq.gz",
        r2 = "fastp/788707_20180129_S_R2.fastq.gz",
        ref = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
    output:
        bam = "788707_20180129.bam"
    shell:
        """
        minimap2 -t 16 --split-prefix=tmp$$ -a -xsr {input.ref} {input.r1} {input.r2} | samtools view -bh | samtools sort -o {output.bam}
        """

rule index_bam:
    input:
        bam = "788707_20180129.bam"
    output:
        bai = "788707_20180129.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

rule human:
    input:
        bam = "788707_20180129.bam"
    output:
        r1 = "human/788707_20180129_S_R1.fastq.gz",
        r2 = "human/788707_20180129_S_R2.fastq.gz"
    shell:
        """
        samtools fastq -F 3588 -f 65 {input.bam} | gzip -c > {output.r1};
        samtools fastq -F 3588 -f 129 {input.bam} | gzip -c > {output.r2};
        """

rule not_human:
    input:
        bam = "788707_20180129.bam"
    output:
        r1 = "not_human/788707_20180129_S_R1.fastq.gz",
        r2 = "not_human/788707_20180129_S_R2.fastq.gz"
    shell:
        """
        samtools fastq -F 3584 -f 77 {input.bam} | gzip -c > {output.r1};
        samtools fastq -F 3584 -f 141 {input.bam} | gzip -c > {output.r2};
        """

rule assemble:
    input:
        r1 = "not_human/788707_20180129_S_R1.fastq.gz",
        r2 = "not_human/788707_20180129_S_R2.fastq.gz"
    output:
        directory("not_human_assembly")
    shell:
        """
        spades.py --meta -1 {input.r1} -2 {input.r2} -o not_human_assembly -t 8
        """
