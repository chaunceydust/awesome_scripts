# update yum
sudo yum -y update

# install basic tools
sudo yum -y install make automake gcc gcc-c++ kernel-devel # https://unix.stackexchange.com/questions/1338/what-is-the-fedora-equivalent-of-the-debian-build-essential-package/259039
sudo yum -y install ncurses-devel # https://www.ostechnix.com/how-to-install-ncurses-library-in-linux/
sudo yum -y groupinstall "Development Tools" # https://www.ostechnix.com/install-development-tools-linux/
sudo yum -y install glibc libSM libICE libXpm libX11 # https://seisman.github.io/SAC_Docs_zh/introduction/linux-install/
sudo yum -y install zlib ncurses
sudo yum -y install zlib-devel.x86_64 libpng-devel.x86_64 libpng12-devel.x86_64 ncurses-devel.x86_64 gcc-c++ bzip2-devel xz-devel # https://octopus-toolkit2.readthedocs.io/en/latest/installation.html
sudo yum -y install python-devel.x86_64
sudo yum -y install git
sudo yum -y install cmake
sudo yum -y install java-1.8.0-openjdk* # https://stackoverflow.com/questions/5104817/how-to-install-java-sdk-on-centos
sudo yum -y install ant # https://docs.wso2.com/display/ESB460/Installing+Apache+Ant+on+Linux

# set up envrionment
curl http://data.biostarhandbook.com/install/bash_profile.txt >> ~/.bash_profile
curl http://data.biostarhandbook.com/install/bashrc.txt >> ~/.bashrc
source ~/.bashrc

# install miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
exit

conda config --add channels conda-forge
conda config --add channels bioconda
conda update conda -y
conda create -y --name bioinfo python=3.6

# activate conda image
conda activate bioinfo

# install softwares
curl http://data.biostarhandbook.com/install/conda.txt | xargs conda install -y

# check status
efetch -db nuccore -id 2 -format gb
mkdir -p ~/bin
curl http://data.biostarhandbook.com/install/doctor.py > ~/bin/doctor.py
chmod +x ~/bin/doctor.py
~/bin/doctor.py
doctor.py --fixme
mkdir -p ~/src
curl ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.zip > ~/src/edirect.zip
unzip -o ~/src/edirect.zip  -d ~/src
cecho 'export PATH=~/src/edirect:$PATH' >> ~/.bash_profile
source  ~/.bash_profile
conda install -y -f perl==5.26.2

# set users
chmod 755 ec2-user/
mkdir info
vim info/users.txt
for user in `cat /home/ec2-user/info/users.txt`; do sudo adduser $user  -g gbh0319; sudo cp /home/ec2-user/.bashrc /home/$user/.bashrc; sudo cp /home/ec2-user/.bash_profile /home/$user/.bash_profile; sudo su - $user -c "mkdir /home/$user/.ssh"; sudo su - $user -c "chmod 700 /home/$user/.ssh"; sudo su - $user -c "ssh-keygen -y -f /home/ec2-user/info/$user.pem > /home/$user/.ssh/authorized_keys"; sudo su - $user -c "chmod 600 /home/$user/.ssh/authorized_keys"; sudo su - $user -c "source /home/$user/.bashrc; source /home/$user/.bash_profile"; sudo su - $user -c "conda init"; done

wget -r http://data.biostarhandbook.com
mv data.biostarhandbook.com Data
