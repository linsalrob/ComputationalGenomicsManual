# Assignment 1

This is a test to make sure that you have the AWS image up and running and that you can access it and generate meaningful results.

It will also familiarize you with [NCBI EDirect](../Databases/NCBI_Edirect.html) that we will use during the course.

The file [genera.txt](https://raw.githubusercontent.com/linsalrob/ComputationalGenomicsManual/master/Assignments/Assignment1/genera.txt) is a plain text file that you can read on the command line using `less` or `more`. The file has 69 different organisms, listed as genus and species, with one entry per line. The assignment is to familiarize yourself with [NCBIs EDirect](../../Databases/NCBI_Edirect.md) tool and to use some advanced bash scripting.

In the first part of the assignment, you must identify which organism has the most number of genomes in the *assembly* database. You should use on of the `edirect` scripts provided on the AWS image to complete that task as shown in the manual.

In the second part of the assignment, you must calculate the **AVERAGE** (mean) genome size of the genomes associated with “**Prevotella buccalis**”. 

*Hint*: There is a command called countfasta.py that will likely help you with this step!