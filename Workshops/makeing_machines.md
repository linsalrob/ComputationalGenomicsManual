# Making machines for the workshops

We usually use [Google Cloud](https://console.google.com), [AWS](https://us-west-2.signin.aws.amazon.com/) or [Nectar](https://dashboard.rc.nectar.org.au/)

- When booting a new machine, make sure that you allow SSH access in the security groups. It needs a group that provides port 22 access to all IP address (0.0.0.0/0)
- Make sure you add additional space for user accounts, usually as an external volume:
    - mount the volume in the web server (e.g. as `/dev/vdb1`)
    - use `fdisk /dev/vdb` to make a new partition
    - use `mkfs.ext4 /dev/vdb1` to format that
    - mount it on `/storage`: `mkdir /storage && mount /dev/vdb1 /storage`
- make user and data directories: `mkdir /storage/users/ /storage/data`
- set the ubuntu password so you can enable passwords: `sudo passwd ubuntu`
- enable passwords:

```
sudo vi /etc/ssh/sshd_config
# Change PasswordAuthentication no to yes.
# Save and exit
sudo /etc/init.d/sshd reload

# for debian use
service sshd reload

# while you are at it:

apt update && apt install build-essential && apt -y dist-upgrade 

```





- create the users with [create_newusers.py](https://raw.githubusercontent.com/linsalrob/EdwardsLab/master/bin/create_newusers.py). Note that this version sets their home directory to `/storage/users/$USER`

```
python ~/GitHubs/EdwardsLab/bin/create_newusers.py -m storage/users -n 200 -s 5
```

- copy these files across to the ubuntu machine and use `newusers` to setup the accounts:

```
scp -i ~/.ssh/xxxxxxxxxxxxx_rsa users.tsv accounts.tsv ubuntu@aaa.bbb.ccc.ddd:
# login to the host
sudo newusers users.tsv
```

- check the you can log in with e.g. user number 100.


# Installing pony

We are trying [pony linux](https://github.com/NCGAS/PonyLinux/blob/master/installation_instructions.md) 

## install software

- install git/texinfo: `apt-get update && apt-get install -y git texinfo`
- install ponysay: 

```
cd /root
git clone https://github.com/erkin/ponysay.git
cd ponysay/
python3 setup.py --freedom=partial install
```

Now install PonyLinux somewhere everyone can use it:

```
cd /storage
mkdir git
cd git
git clone https://github.com/NCGAS/PonyLinux.git
```

Next, link this to each users home directory:

```
cd /storage/users/
for USER in user*; do ln -s /storage/git/PonyLinux $USER/PonyLinux; chown $USER:$USER $USER/PonyLinux; done 
```

Now when they log in, they can `cd PonyLinux; ./PonyLinux.sh`


# Install bioconda

Download miniforge:


```
curl -LO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash  Miniforge3-Linux-x86_64.sh
```


Make a ~/.bash_profile file

```
vi bash_profile

```

Paste this text:
```
if [ -f ~/.bashrc ]; then    
	. ~/.bashrc
fi
```

Then copy to all users:

```
for USER in user0*; do cp bash_profile $USER/.bash_profile; chown $USER:$USER $USER/.bash_profile; done
```


# Install binchicken and download databases


```
mamba create -n binchicken -y -c bioconda -c conda-forge "binchicken>=0.12.5"

mamba activate binchicken
binchicken build  --conda-prefix miniforge3  --singlem-metapackage metapackage --checkm2-db checkm2 --download-databases

```










