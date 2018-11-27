Annotating the function in the bins with SUPER-FOCUS

Subsystems Profile by database reduction using FOCUS (SUPER-FOCUS) combines the taxonomic profiling using FOCUS, which is rapid, and the functional profiling using k-mers or blast. 

Download SUPER-FOCUS from “http://edwards.sdsu.edu/superfocus/” to run the bins through SUPERFOCUS to predict the subsystems in each bin. SUPER-FOCUS input can be the same as FOCUS, just the folder with all the bins. Example of the SUPER-FOCUS input file are saved under “/ec2-user/SUPER-FOCUS/SUPER-FOCUS_Input_Files/”

First download the database by running 

python superfocus__download.py rapsearch ( to download the database for rapserach algorithm)

Type the command , an example for Algae_MetaBatBins 

python superfocus.py -q Algae_MetaBatBins -m 1 -db DB_90 -dir Algae_MetaBatBins_Superfocus_results. 

The SUPER-FOCUS results for all the bins are saved in /ec2-user/SUPER-FOCUS/” 

