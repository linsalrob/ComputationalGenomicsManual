"""

Process some metagenomes

Rob Edwards, September 2020

"""


import os
import sys


# set this to whatever the name of your directory
# with the reads is. If you are following along with the
# tutorial, you can leave this as fastq
READDIR = 'fastq/'

# Note that this example requires an 1 file AND an 2 file
# and that each file should match *_1.* and *_2.* 
SAMPLES,EXTENSIONS, = glob_wildcards(os.path.join(READDIR, '{sample}_1.{extensions}'))

# just get the first file extension as we don't need to iterate all of them
file_extension = EXTENSIONS[0]


# just check there is something to actually do!
if len(SAMPLES) == 0:
    sys.stderr.write("FATAL: We could not detect any samples at all.\n")
    sys.stderr.write(f"Do you have a directory called {READDIR} with some fastq files in it?\n")
    sys.stderr.write("Do those fastq files have _1. and _2.?\n")
    sys.exit()


rule all:
    input:
        expand(
            [
                os.path.join("results", "{sample}", "prinseq", "{sample}_R1.good.fastq"), 
                os.path.join("results", "{sample}", "prinseq", "{sample}_R2.good.fastq"),
                os.path.join("results", "{sample}", "focus", "output_All_levels.csv"),
                os.path.join("results", "{sample}", "superfocus", "output_all_levels_and_function.xls"),
                os.path.join("results", "{sample}", "kraken", "{sample}.report.tsv"),
                os.path.join("results", "{sample}", "kraken", "{sample}.output.tsv"),
                os.path.join("results", "{sample}", "assembly", "done"),
                os.path.join("results", "{sample}", "{sample}_contigs.bam"),
                os.path.join("results", "{sample}", "{sample}_contigs.tsv")
            ],
               sample=SAMPLES),
        os.path.join("results", "assemblies", "assembly", "done"),
        "results/coverage.tsv",
        "results/pearson.tsv"


rule run_prinseq:
    input:
        r1 = os.path.join(READDIR, "{sample}_1." + file_extension),
        r2 = os.path.join(READDIR, "{sample}_2." + file_extension)
    output:
        d = directory(os.path.join("results", "{sample}", "prinseq")),
        r1 = os.path.join("results", "{sample}", "prinseq", "{sample}_R1.good.fastq"),
        s1 = os.path.join("results", "{sample}", "prinseq", "{sample}_R1.single.fastq"),
        b1 = temporary(os.path.join("results", "{sample}", "prinseq", "{sample}_1.bad.fastq")),
        r2 = os.path.join("results", "{sample}", "prinseq", "{sample}_R2.good.fastq"),
        s2 = os.path.join("results", "{sample}", "prinseq", "{sample}_R2.single.fastq"),
        b2 = temporary(os.path.join("results", "{sample}", "prinseq", "{sample}_2.bad.fastq"))
    resources:
        cpus=2,
        mem_mb=4000
    shell:
        """
        mkdir -p {output.d} &&
        prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 1 -derep 1 \
            -out_format 0 -trim_tail_left 5 -trim_tail_right 5 \
            -trim_qual_type min -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10 \
            -out_good {output.r1} -out_single {output.s1} -out_bad {output.b1} \
            -out_good2 {output.r2} -out_single2 {output.s2} -out_bad2 {output.b2} \
            -fastq {input.r1} -fastq2 {input.r2} \
            -threads {resources.cpus}
        """

rule run_focus:
    input:
        os.path.join("results", "{sample}", "prinseq")
    output:
        d = directory(os.path.join("results", "{sample}", "focus")),
        a = os.path.join("results", "{sample}", "focus", "output_All_levels.csv"),
        f = os.path.join("results", "{sample}", "focus", "output_Family_tabular.csv"),
        k = os.path.join("results", "{sample}", "focus", "output_Kingdom_tabular.csv"),
        p = os.path.join("results", "{sample}", "focus", "output_Phylum_tabular.csv"),
        s = os.path.join("results", "{sample}", "focus", "output_Strain_tabular.csv"),
        c = os.path.join("results", "{sample}", "focus", "output_Class_tabular.csv"),
        g = os.path.join("results", "{sample}", "focus", "output_Genus_tabular.csv"),
        o = os.path.join("results", "{sample}", "focus", "output_Order_tabular.csv"),
        sp = os.path.join("results", "{sample}", "focus", "output_Species_tabular.csv")
    resources:
        cpus=2,
        mem_mb=4000
    shell:
        """
        focus -q {input} -o {output.d} -t {resources.cpus}
        """



rule run_superfocus:
    input:
        os.path.join("results", "{sample}", "prinseq")
    output:
        d = directory(os.path.join("results", "{sample}", "superfocus")),
        a = os.path.join("results", "{sample}", "superfocus", "output_all_levels_and_function.xls"),
        l1 = os.path.join("results", "{sample}", "superfocus", "output_subsystem_level_1.xls"),
        l2 = os.path.join("results", "{sample}", "superfocus", "output_subsystem_level_2.xls"),
        l3 = os.path.join("results", "{sample}", "superfocus", "output_subsystem_level_3.xls"),
    resources:
        cpus=4,
        mem_mb=16000
    shell:
        """
        superfocus -q {input} -dir {output.d} -a diamond -t {resources.cpus}
        """

rule run_kraken:
    input:
        r1 = os.path.join("results", "{sample}", "prinseq", "{sample}_R1.good.fastq"),
        r2 = os.path.join("results", "{sample}", "prinseq", "{sample}_R2.good.fastq")
    output:
        rt = os.path.join("results", "{sample}", "kraken", "{sample}.kraken_taxonomy.tsv"),
        ot = os.path.join("results", "{sample}", "kraken", "{sample}.kraken_output.tsv")
    resources:
        cpus=4,
        mem_mb=400000
    shell:
        """
        kraken2 --paired --report {output.rt} \
                --output {output.ot} \
                --threads {resources.cpus} \
                {input.r1} {input.r2}
        """


rule megahit_assemble:
    input:
        r1 = os.path.join("results", "{sample}", "prinseq", "{sample}_R1.good.fastq"),
        s1 = os.path.join("results", "{sample}", "prinseq", "{sample}_R1.single.fastq"),
        r2 = os.path.join("results", "{sample}", "prinseq", "{sample}_R2.good.fastq"),
        s2 = os.path.join("results", "{sample}", "prinseq", "{sample}_R2.single.fastq"),
    output:
        os.path.join("results", "{sample}", "assembly", "final.contigs.fa"),
        os.path.join("results", "{sample}", "assembly", "checkpoints.txt"),
        os.path.join("results", "{sample}", "assembly", "done"),
        directory(os.path.join("results", "{sample}", "assembly", "intermediate_contigs")),
        os.path.join("results", "{sample}", "assembly", "log"),
        os.path.join("results", "{sample}", "assembly", "options.json")
    params:
        odir = directory(os.path.join("results", "{sample}", "assembly"))
    resources:
        cpus=4,
        mem_mb=16000
    shell:
        """
        rmdir {params.odir} && \
        megahit -1 {input.r1} -2 {input.r2} -o {params.odir} -t {resources.cpus}
        """


rule combine_reads:
    input:
        r1 = expand(os.path.join("results", "{smpl}", "prinseq", "{smpl}_R1.good.fastq"), smpl=SAMPLES),
        s1 = expand(os.path.join("results", "{smpl}", "prinseq", "{smpl}_R1.single.fastq"), smpl=SAMPLES),
        r2 = expand(os.path.join("results", "{smpl}", "prinseq", "{smpl}_R2.good.fastq"), smpl=SAMPLES),
        s2 = expand(os.path.join("results", "{smpl}", "prinseq", "{smpl}_R2.single.fastq"), smpl=SAMPLES)
    output:
        r1 = os.path.join("results", "assemblies", "R1.fastq"),
        r2 = os.path.join("results", "assemblies", "R2.fastq"),
        s = os.path.join("results", "assemblies", "Singles.fastq"),
    shell:
        """
        cat {input.r1} > {output.r1} &&
        cat {input.r2} > {output.r2} &&
        cat {input.s1} {input.s2} > {output.s}
        """

rule assemble_all:
    input:
        r1 = os.path.join("results", "assemblies", "R1.fastq"),
        r2 = os.path.join("results", "assemblies", "R2.fastq"),
        s = os.path.join("results", "assemblies", "Singles.fastq"),
    output:
        os.path.join("results", "assemblies", "assembly", "final.contigs.fa"),
        os.path.join("results", "assemblies", "assembly", "checkpoints.txt"),
        os.path.join("results", "assemblies", "assembly", "done"),
        directory(os.path.join("results", "assemblies", "assembly", "intermediate_contigs")),
        os.path.join("results", "assemblies", "assembly", "log"),
        os.path.join("results", "assemblies", "assembly", "options.json")
    params:
        odir = directory(os.path.join("results", "assemblies", "assembly"))
    resources:
        cpus=16,
        mem_mb=40000
    shell:
        """
        rmdir {params.odir} && \
        megahit -1 {input.r1} -2 {input.r2} -o {params.odir} -t {resources.cpus} --min-contig-len 2000
        """

rule index_minimap:
    input:
        os.path.join("results", "assemblies", "assembly", "final.contigs.fa")
    output:
        os.path.join("results", "mapping", "contigs.mmi")
    params:
        o = os.path.join("results", "mapping")
    shell:
        """
        mkdir -p {params.o} && minimap2 -x sr -d {output} {input}
        """

rule run_minimap_all:
    input:
        mm = os.path.join("results", "mapping", "contigs.mmi"),
        r1 = os.path.join("results", "{sample}", "prinseq", "{sample}_R1.good.fastq"),
        r2 = os.path.join("results", "{sample}", "prinseq", "{sample}_R2.good.fastq"),
    output:
        os.path.join("results", "{sample}", "{sample}_contigs.bam")
    resources:
        cpus=16,
        mem_mb=40000
    shell:
        """

        minimap2 -ax sr {input.mm} \
                {input.r1} {input.r2} | \
                samtools view -@ {resources.cpus} -b -F 4  | \
                samtools sort -@ {resources.cpus} > \
                {output}
        """


rule indexed_bam:
    input:
        os.path.join("results", "{sample}", "{sample}_contigs.bam")
    output:
        os.path.join("results", "{sample}", "{sample}_contigs.bam.bai")
    shell:
        """
        samtools index {input}
        """

rule count_stats:
    input:
        bam = os.path.join("results", "{sample}", "{sample}_contigs.bam"),
        bai = os.path.join("results", "{sample}", "{sample}_contigs.bam.bai")
    output:
        os.path.join("results", "{sample}", "{sample}_contigs.tsv")
    shell:
        """
        samtools idxstats {input.bam} | cut -f 1,3 > {output}
        """

rule make_table:
    input:
        expand(os.path.join("results", "{smpl}", "{smpl}_contigs.tsv"), smpl=SAMPLES)
    output:
        "results/coverage.tsv"
    shell:
        """
        joinlists.pl -t {input} > {output}
        """

rule pearson:
    input:
        "results/coverage.tsv"
    output:
        "results/pearson.tsv"
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/bin/correlations.py  -l -r -f {input} > {output}
        """
