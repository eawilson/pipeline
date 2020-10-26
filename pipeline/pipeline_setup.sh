sudo apt-get -y update
sudo apt-get -y dist-upgrade

sudo apt-get -y install build-essential
sudo apt-get -y install zlib1g-dev # bwa, htslib
sudo apt-get -y install libbz2-dev # htslib
sudo apt-get -y install liblzma-dev # htslib
sudo apt-get -y install autoconf # htslib
sudo apt-get -y install libcurl4-gnutls-dev # htslib
sudo apt-get -y install libncurses5-dev # samtools
sudo apt-get -y install openjdk-8-jdk # varscan picard
sudo apt-get -y install libmysqlclient-dev # vep
sudo apt-get -y install python3-pip
sudo apt-get -y install python3-tk
sudo apt-get -y install unzip

sudo cpan App::cpanminus # vep
sudo cpanm Archive::Zip # vep
sudo cpanm DBD::mysql # vep
sudo cpanm Archive::Extract # vep
sudo cpanm Try::Tiny # vep
sudo cpanm JSON
sudo cpanm Module::Build

sudo pip3 install boto3
sudo pip3 install matplotlib


wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/\$latest/amd64-linux/bs?bt_package=latest" -O bs
chmod +x bs
sudo mv /usr/local/bin/bs

curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
rm awscliv2.zip
sudo rm -rf aws

git clone https://github.com/lh3/bwa.git
cd bwa
make
sudo mv bwa /usr/local/bin
cd ..

wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
tar xvfj htslib-1.10.2.tar.bz2
rm htslib-1.10.2.tar.bz2
cd htslib-1.10.2
autoreconf
./configue
make
sudo make install
cd ..

wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
tar xvfj samtools-1.10.tar.bz2
rm samtools-1.10.tar.bz2
cd samtools-1.10
autoreconf
./configue
make
sudo make install
cd ..

wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
tar xvfj bcftools-1.10.2.tar.bz2
rm bcftools-1.10.2.tar.bz2
cd bcftools-1.10.2
autoreconf
./configue
make
sudo make install
cd ..

# Not the latest version, latest version fails ?intermediate version
wget -O VarScan.v2.3.9.jar https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar/download
sudo mv VarScan.v2.3.9.jar /usr/local/bin
echo '#!/bin/bash' >varscan
echo 'java -jar /usr/local/bin/VarScan.v2.3.9.jar $@' >>varscan
chmod +x varscan
sudo mv varscan /usr/local/bin

sudo bash -c "$(curl -L https://basemount.basespace.illumina.com/install)"


sudo mkdir /usr/local/lib/site_perl
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
sudo perl INSTALL.pl --AUTO a -d /usr/local/lib/site_perl
sudo cp -r  modules/Bio/EnsEMBL /usr/local/lib/site_perl/Bio
sudo mv vep /usr/local/bin
cd ..



wget -O picard-2.23.4.jar https://github.com/broadinstitute/picard/releases/download/2.23.4/picard.jar
sudo mv picard-2.23.4.jar /usr/local/bin
echo '#!/bin/bash' >picard
echo 'java -jar /usr/local/bin/picard-2.23.4.jar $@' >>picard
chmod +x picard
sudo mv picard /usr/local/bin

git clone git@github.com:eawilson/pipeline.git
cd pipeline
sudo python3 setup.py develop
cd ..

git clone git@github.com:eawilson/dude.git
cd dude
make
sudo make install
cd ..

git clone https://github.com/eawilson/covermi.git
cd covermi
sudo python3 setup.py install
cd ..




