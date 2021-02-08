##################################################################
#                                                                #
# Minimal set up for a one day workshop. This only installs      #
# a few pieces of code for us, and is not a full install.        #
#                                                                #
# This assumes that we just become root to save typing sudo      #
#                                                                #
# Rob Edwards, Jan 2021                                          #
#                                                                #
##################################################################


# MANUAL SETUP
# You need to get access credentials from AWS before doing this. You need the Access key ID and	Secret access key)
 aws configure



##################################################################
#                                                                #
#                      Basic System Installation                 #
#                                                                #
##################################################################

# Much of this installation is just to make sure we have things for the class, should we need them.
# We are probably not going to use most of this stuff!

amazon-linux-extras install -y vim R4 python3.8 epel
yum-config-manager --enable epel
yum update -y

yum -y install gcc gcc-c++ libgcc ncurses ncurses-devel zlib-devel bzip2-devel xz-devel git  freetype.x86_64 freetype-devel.x86_64 libpng.x86_64 libpng-devel.x86_64 lapack.x86_64 lapack-devel.x86_64 gsl-devel.x86_64 gsl.x86_64 atlas.x86_64 atlas-devel.x86_64 boost  mesa-libGL-devel.x86_64 mesa-libGLU.x86_64 mesa-libGL.x86_64 mesa-libGLU-devel.x86_64 gmp.x86_64 gmp-devel.x86_64 libxslt-devel.x86_64 libxslt.x86_64 libxml2.x86_64 libxml2-devel.x86_64  libattr-devel perl-CPAN perl-IO-Socket-SSL perl-Net-SSLeay.x86_64   numpy tbb  tbb-devel cmake perl-Time-Piece perl-XML-Simple perl-Digest-MD5 git java perl-CPAN perl-Module-Build  gnuplot ghostscript   libcurl-devel boost-devel pigz  perl-LWP-Protocol-https  python38-devel tmux

# set python3 to work
update-alternatives --install /usr/bin/python3 python3 /bin/python3.8 1


# note that we do this in two steps to install dependencies first!
pip-3.8 install numpy matplotlib scipy cython virtualenv h5py checkm ete3 cherrypy bottle scikit-learn  lxml beautifulsoup4 scikit-learn pandas statsmodels MUSiCC biopython bottle 

pip-3.8 install checkm fishtaco numpy python_libsbml requests ipykernel jupyter gffutils metagenomics-focus superfocus pysam  snakemake


### This will require some manual input, so do this ONE LINE AT A TIME :)
## At the moment it seems like Prokka is the only thing that needs this.
# make sure we can set things up automatically
perl -MCPAN -e shell

perl -MCPAN -e 'my $c = "CPAN::HandleConfig"; $c->load(doit => 1, autoconfig => 1); $c->edit(prerequisites_policy => "follow"); $c->edit(build_requires_install_policy => "yes"); $c->commit'

#install bioperl 
cpan -i Bio::Perl 




##################################################################
#                                                                #
#    Individual Software Packagegs.                              #
#                                                                #
# These were chosen from the AWS_Setup.txt file in this          #
# directory because they may be useful for the workshop.         #
#                                                                #
##################################################################

# Set up bash fasta fail
set -euo pipefail

# fastq2fasta
mkdir -p /usr/local/genome/fastq2fasta
cd /usr/local/genome/fastq2fasta
wget https://raw.githubusercontent.com/linsalrob/EdwardsLab/master/bin/fastq2fasta.cpp
c++ -o fastq2fasta fastq2fasta.cpp 
rm -f fastq2fasta.cpp 
echo "export PATH=\$PATH:$PWD" > /etc/profile.d/fastqfasta.sh

# fastqc:
mkdir -p /usr/local/genome/fastqc
cd /usr/local/genome/fastqc
curl -Lo fastqc.zip https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc.zip 
cd FastQC/
chmod 755 fastqc
echo "export PATH=\$PATH:$PWD" > /etc/profile.d/fastqc.sh

# Jellyfish
mkdir  -p /usr/local/genome/jellyfish
cd /usr/local/genome/jellyfish 
curl -Lo jellyfish.tgz https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz
tar xf jellyfish.tgz
cd jellyfish-2.2.10
echo "include /usr/local/lib" >> /etc/ld.so.conf
ldconfig
./configure && make && make install


# SRA toolkit
mkdir -p /usr/local/genome/sra
cd /usr/local/genome/sra
curl -Lo sratoolkit.tgz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar xf sratoolkit.tgz 
rm -f sratoolkit.tgz 
ln -s sratoolkit.*centos_linux64 current
cd current/bin/
echo "export PATH=\$PATH:$PWD" > /etc/profile.d/sra.sh


# Prinseq
mkdir /usr/local/genome/prinseq
cd /usr/local/genome/prinseq
curl -Lo prinseq.tgz https://github.com/Adrian-Cantu/PRINSEQ-plus-plus/releases/download/v1.2/prinseq-plus-plus-1.2.tar.gz
tar xf prinseq.tgz 
rm -f prinseq.tgz 
ln -s prinseq-plus-plus-1.2/ current
cd current/
./configure  && make  && make install


# samtools 
# CHECK THE VERSION HERE: https://github.com/samtools/samtools/releases
mkdir /usr/local/genome/sam_bcf_hts
export VERSION=1.11
cd /usr/local/genome/sam_bcf_hts
curl -Lo samtools-$VERSION.tar.bz2  https://github.com/samtools/samtools/releases/download/$VERSION/samtools-$VERSION.tar.bz2
curl -Lo bcftools-$VERSION.tar.bz2 https://github.com/samtools/bcftools/releases/download/$VERSION/bcftools-$VERSION.tar.bz2
curl -Lo htslib-$VERSION.tar.bz2 https://github.com/samtools/htslib/releases/download/$VERSION/htslib-$VERSION.tar.bz2
tar xf htslib-$VERSION.tar.bz2 
tar xf bcftools-$VERSION.tar.bz2 
tar xf samtools-$VERSION.tar.bz2 
cd htslib-$VERSION
./configure --prefix=/usr/local/genome/sam_bcf_hts && make && make install
cd ../bcftools-$VERSION
./configure --prefix=/usr/local/genome/sam_bcf_hts && make && make install
cd ../samtools-$VERSION
./configure --prefix=/usr/local/genome/sam_bcf_hts && make && make install
cd ..
rm -rf bcftools-$VERSION  bcftools-$VERSION.tar.bz2 htslib-$VERSION  htslib-$VERSION.tar.bz2 samtools-$VERSION  samtools-$VERSION.tar.bz2
cd bin/
echo "export PATH=\$PATH:$PWD" > /etc/profile.d/samtools.sh
export PATH=$PATH:$PWD


# Diamond
# CHECK THE LATEST VERSION HERE
export DIA_VER=2.0.6
mkdir -p /usr/local/genome/diamond
cd /usr/local/genome/diamond
curl -LO https://github.com/bbuchfink/diamond/releases/download/v$DIA_VER/diamond-linux64.tar.gz
tar xf diamond-linux64.tar.gz
rm -f diamond-linux64.tar.gz
echo "export PATH=\$PATH:$PWD" > /etc/profile.d/diamond.sh

#minimap2
cd /usr/local/genome
git clone https://github.com/lh3/minimap2.git
cd minimap2
make -j
echo "export PATH=\$PATH:$PWD" > /etc/profile.d/minimap.sh



# Kraken2
mkdir -p /usr/local/genome/kraken/current
cd /usr/local/genome/kraken
K_VER=2.1.1
curl -LO https://github.com/DerrickWood/kraken2/archive/v${K_VER}.tar.gz
tar xf v${K_VER}.tar.gz
cd kraken2-${K_VER}/
./install_kraken2.sh /usr/local/genome/kraken/current/
cd ..
rm -rf kraken2-${K_VER}/ v${K_VER}.tar.gz
cd current
# This will need aws configure run before working
# this also saves disk by untarring the stream
aws s3 cp s3://genome-idx/kraken/k2_standard_8gb_20201202.tar.gz - | tar -zx
echo -e "export PATH=\$PATH:$PWD\nexport KRAKEN2_NUM_THREADS=4\nexport KRAKEN2_DB_PATH=$PWD\nexport KRAKEN2_DEFAULT_DB=$PWD" > /etc/profile.d/kraken.sh


# Megahit
MH_VER=1.2.9
mkdir -p /usr/local/genome/megahit
cd /usr/local/genome/megahit
curl -LO https://github.com/voutcn/megahit/releases/download/v${MH_VER}/MEGAHIT-${MH_VER}-Linux-x86_64-static.tar.gz
tar xf MEGAHIT-${MH_VER}-Linux-x86_64-static.tar.gz
ln -s MEGAHIT-${MH_VER}-Linux-x86_64-static current
cd current/bin
echo "export PATH=\$PATH:$PWD" > /etc/profile.d/megahit.sh


# Edwards lab
cd /usr/local/genome
git clone https://github.com/linsalrob/EdwardsLab.git
cd EdwardsLab/bin
echo "export PATH=\$PATH:$PWD" > /etc/profile.d/edwardslab.sh
cd ..
echo "export PYTHONPATH=\$PYTHONPATH:$PWD" >> /etc/profile.d/edwardslab.sh


# prokka
# Requires Bio::Perl from above
cd /usr/local/genome/
git clone https://github.com/tseemann/prokka.git
cd prokka
cd bin
./prokka --setupdb
echo "export PATH=\$PATH:$PWD" > /etc/profile.d/prokka.sh



exit

##################################################################
#                                                                #
#                          Manual Installation                   #
#                                                                #
##################################################################


## Finish the focus installation (you could probably do this automatically)
cd /usr/local/lib/python3.8/site-packages/focus_app/
unzip db.zip 
rm -f db.zip

## Finish the superfocus installation.
# Check the diamond version
/usr/local/genome/diamond/diamond --version
# if this is above 1 (most likely) you can continue with version 3 databases as shown here

export PATH=$PATH:/usr/local/bin
# check superfocus
which superfocus
# if this doesn't work, don't proceed!

mkdir ~/superfocus
cd ~/superfocus
FDIR=$(which superfocus | sed -e 's#bin/superfocus$#lib/python3.8/site-packages/superfocus_app/db/#')
curl -LO https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/diamond_v3/90_clusters.db.dmnd.zip
mkdir -p  $FDIR/static/diamond
unzip -d  $FDIR/static/diamond 90_clusters.db.dmnd.zip
if [ ! -e $FDIR/database_PKs.txt ]; then curl -Lo $FDIR/database_PKs.txt https://edwards.sdsu.edu/SUPERFOCUS/downloads/database_PKs.txt; fi
mkdir -p $FDIR/tmp 
chmod 777 $FDIR/tmp/




exit

##################################################################
#                                                                #
#   Install some data                                            #
#                                                                #
# NOTE: If you are using AWS in US N. Virginia (us-east-1)       #
#       you can download the data using AWS. Otherwise you       #
#       will have to download it from ncbi. Set up vdb-config    #
#       either way, since it needs to be done!                   #
#                                                                #
##################################################################

# Interactive
vdb-config -i

{keys are: c i a r s <enter> x}

# you can check whether you are using NCBI or Cloud. NCBI will have short path, cloud long path.
srapath SRR000001
mkdir -p /data/gut
cd /data/gut
for SRR_ID in SRR3466404 SRR3506419 SRR3506420 SRR3546776 SRR3546778 SRR3546779 SRR3546780 SRR3546781 SRR3546782; do fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip -N 5000 -X 255000 $SRR_ID; done

##################################################################
#                                                                #
#   Create some user accounts                                    #
#                                                                #
# NOTE: The new users text file here is not the one that we      #
#       really user because it has different passwords. Sorry!   #
#                                                                #
##################################################################

set up the newusers.txt file

~/GitHubs/EdwardsLab/bin/create_newusers.py 

then use newusers users.txt

This will add all the users.




edit /etc/ssh/sshd_config
change PasswordAuthentication yes
restart sshd: 
systemctl restart sshd


# set up accounts. This assumes you have your keys in /home/ec2-user/.ssh/authorized_keys and you've copied over your .vim stuff!

for USER in $(cut -f 1 -d : newusers.txt); do mkdir -p /home/$USER/.ssh; cp /home/ec2-user/.ssh/authorized_keys /home/$USER/.ssh; cp -r /home/ec2-user/.vim* /home/ec2-user/.bashrc /home/$USER/; chown -R $USER.$USER  /home/$USER/; chmod 600 /home/$USER/.ssh/authorized_keys; chmod 700 /home/$USER/.ssh; done

















