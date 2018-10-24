# Sample Metagenomes

<img src="images/abrolhos.png" alt="Abrolhos Islands" align="right" />

We have several sample metagenomes that we provide for use in the course. However, please feel free to bring your own metagenomes and data sets, and replace their filenames with the names we use here. Our sample data sets are from the Abrolhos region of Brazil, and are some low diversity metagenomes from an experiment where we tested the effects of Coral, Algae, CCA, or no treatment on the growth of microbes over time.

There are four groups of data:

*Algae treatment* &mdash; 4 replicates (Algae_11, Algae_12, Algae_13, Algae_14)
*CCA treatment* &mdash; 3 replicates (CCA_11, CCA_12, CCA_13)
*Control treatment* &mdash; 4 replicates (Control_11, Control_12, Control_13, Control_14)
*Coral treatment* &mdash; 4 replicates (Coral_11, Coral_12, Coral_13, Coral_14)

During the course, we are going to process these data sets as we normally process our metagenomes. We will start with quality control, perform some functional and taxonomic annotations, and then end up binning the samples to get individual metagenomes.

These samples were sequenced on an Ion Torrent, and so you will see quality differences and we’ll need to use the iontorrent flag when assembling them.

In the example text below, I use the Algae samples to demonstrate the commands. Be sure to switch the `Algae_12.fna` name to the file name that you are working on.
