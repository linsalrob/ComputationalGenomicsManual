# Snakemake scripts

You should use these to run everything in one go, especially if you use a cluster!

- [process_metagenomesSAGC_conda.snakefile](process_metagenomesSAGC_conda.snakefile) A snakefile that uses conda (add the `--use-conda` flag) and installs the software you need. You will also need to download and extract [envs.zip](envs.zip), the conda environments.
- [process_metagenomesSAGC.snakefile](process_metagenomesSAGC.snakefile) A snakefile that assumes you have everything installed. This may not work!
