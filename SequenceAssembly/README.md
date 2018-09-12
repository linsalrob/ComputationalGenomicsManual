# Sequence Assembly

The essential problem that sequence assembly is trying to overcome is that the average microbial genome is 2,000,000 bp, while the typical sequence read length is 150-300 bp for Illumina sequences and upto 50,000 bp for PacBio or Nanopore sequences.

With such (relatively) short sequences, how can we assemble a whole, or nearly whole, genome?

The answer is by repetitively sequencing the same thing over and over again! If we start each sequence at a random location, and we have enough sequences, eventually we can join those sequences together to form what we call **contigs**.

