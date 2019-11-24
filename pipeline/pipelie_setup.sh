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
sudo apt-get install python3-tk

sudo cpan App::cpanminus # vep
sudo cpanm Archive::Zip # vep
sudo cpanm DBD::mysql # vep
sudo cpanm Archive::Extract # vep
sudo cpanm Try::Tiny # vep
sudo cpanm JSON

sudo pip3 install boto3
sudo pip3 install setuptools
sudo pip3 install matplotlib

mkdir workspace
sudo mkfs -t ext4 /dev/nvme0n1
sudo mount /dev/nvme0n1 workspace
sudo chmod go+rwx workspace



git clone https://github.com/lh3/bwa.git
cd bwa
make
sudo mv bwa /usr/local/bin
cd ..

wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar xvfj htslib-1.9.tar.bz2
rm htslib-1.9.tar.bz2
cd htslib-1.9
autoreconf
./configue
make
sudo make install
cd ..

wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar xvfj samtools-1.9.tar.bz2
rm samtools-1.9.tar.bz2
cd samtools-1.9
autoreconf
./configue
make
sudo make install
cd ..

wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar xvfj bcftools-1.9.tar.bz2
rm bcftools-1.9.tar.bz2
cd bcftools-1.9
autoreconf
./configue
make
sudo make install
cd ..

wget -O VarScan.v2.3.9.jar https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar/download
sudo mv VarScan.v2.3.9.jar /usr/local/bin
sudo ln -s /usr/local/bin/VarScan.v2.3.9.jar /usr/local/bin/varscan.jar

sudo bash -c "$(curl -L https://basemount.basespace.illumina.com/install)"

git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl --AUTO a
sudo cp vep /usr/local/bin/vep-98
sudo ln -s /usr/local/bin/vep-98 /usr/local/bin/vep
cd ..

wget -O picard-2.21.2.jar https://github.com/broadinstitute/picard/releases/download/2.21.2/picard.jar
sudo mv picard-2.21.2.jar /usr/local/bin
sudo ln -s /usr/local/bin/picard-2.21.2.jar /usr/local/bin/picard.jar

git clone git@github.com:eawilson/pipeline.git
cd pipeline
sudo python3 setup.py develop
cd ..

git clone https://github.com/eawilson/covermi.git
cd covermi
sudo python3 setup.py install
cd ..




