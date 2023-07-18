# Using snakemake

Snakemake gets its inspiration from the program _Make_.
_Make_ is an old program that was originally designed to compile and install software.
It was appropriated by bioinformaticians because the program is ideal for writing pipelines.

In other examples, we started with two fastq files, removed adapter sequences, mapped them to the human genome, and then separated out the human and not human sequences.

We can combine all of that into a single `snakemake` file, and it will do all of the steps for us.

[snakemake](https://snakemake.readthedocs.io/en/stable/) is a way to write reproducible code. There are lots of [tutorials about snakemake](https://www.google.com/search?q=snakemake+tutorial) and we have a lot of [Snakemake](https://edwards.flinders.edu.au/?s=snakemake) tutorials on the Edwards' lab website.

### The dataset

You can use any dataset you like for this step, but for the example, we'll use the CF data](../Datasets/CF)

## Writing the snakefile


You can use any text editor to write your `snakefile`. I recommend `nano` on your linux machine:

```bash
nano
```

The main commands are shown at the bottom of the screen. In these commands, the `^` means press the `ctrl` key.

```
^G Get Help      ^O Write Out     ^W Where Is      ^K Cut Text      ^J Justify       ^C Cur Pos       M-U Undo         M-A Mark Text    M-] To Bracket   M-Q Previous     ^B Back
^X Exit          ^R Read File     ^\ Replace       ^U Paste Text    ^T To Spell      ^_ Go To Line    M-E Redo         M-6 Copy Text    ^Q Where Was     M-W Next         ^F Forward
```


Snakemake uses standard [python](https://www.python.org/) python commands, and so if you have used python before you should find it familiar. There are some additions on top of python that you will need to learn.

Using the same concept as before, you start with your destination, _then_ you describe the journey.
__By default, Snakemake runs the first rule it encounters in the Snakefile.__
This is important because we use the first rule--`rule all`-- to declare the pipeline targets.
There are no outputs or shell command for rule all,
so the rule just tells Snakemake that we need the file "A.sam".

Create a file called `Snakefile` by using `nano`:

```
nano Snakefile
```

and add this text:


```python

rule all:
    input:
        "fastp/788707_20180129_S_R1.fastq.gz",
        "fastp/788707_20180129_S_R2.fastq.gz"

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
```

Exit `nano` by pressing `ctrl-x` and then pressing `y` to save the file.

This is your pipeline and Snakemake will recognise it when you run Snakemake.
The only thing you need to tell Snakemake at this stage is how many concurrent jobs to run.
For only 1 job at a time, you would run this:

```bash
snakemake -j 1
```

By default, `snakemake` looks for a file called `Snakefile` but that is not very imaginative, and you probably also want to call it something meaningful so you can find it again. You can use the `-s` flag to snakemake to get it to use a different file.

I normally name my snakefiles after the pipeline or process, and use a `.smk` file extension, so perhaps you might call this `preprocess.smk`. 


```bash 
snakemake -s preprocess.smk -j 1
```

# adding the minimap command

Edit your `Snakefile` it looks like this. 

*NOTE Make sure you add the additional line ("788707_20180129.bam.bai") to the `rule all` input!*

```
nano Snakefile
```

```python

rule all:
    input:
        "fastp/788707_20180129_S_R1.fastq.gz",
        "fastp/788707_20180129_S_R2.fastq.gz",
        "788707_20180129.bam.bai"

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

```


This modified snakefile adds three additional commands:

1. In the `rule all` we added the "788707_20180129.bam.bai" file that is the final output of the indexed bam file
2. The `rule run_minimap` performs the minimap alignment and creates the bam file
3. The `rule index_bam` indexes the bam file to make accessing it much faster


*But note!* In this command, we don't request the creation of the file `788707_20180129.bam`, but because it is required as an input to the command `index_bam`, snakemake knows that it has to make that as an output file. The _only_ way to make that file is to run the rule `run_minimap`, and that requires as input the two `fastq` files in the `fastp` folder, and so that needs the rule `run fastp` to complete! 

That's the beauty of snakemake, now you can always run the same commands.












