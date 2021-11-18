#!/usr/bin/env bash

set -e

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
sudo apt-get -y install python # gatk
sudo apt-get -y install unzip

yes | sudo cpan App::cpanminus # vep
sudo cpanm Archive::Zip # vep
sudo cpanm DBD::mysql # vep
sudo cpanm JSON # vep
sudo cpanm Set::IntervalTree # vep
sudo cpanm PerlIO::gzip # vep
sudo cpanm Try::Tiny # vep


sudo pip3 install boto3
sudo pip3 install matplotlib
sudo pip3 install scipy


wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O bs
chmod u+x bs
sudo mv bs /usr/local/bin


curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install --update
rm awscliv2.zip
sudo rm -rf aws


wget "https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2"
tar xvfj bwa-0.7.17.tar.bz2
rm bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
sudo mv bwa /usr/local/bin
cd ..
rm -rf bwa-0.7.17


wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
tar xvfj bwa-mem2-2.2.1_x64-linux.tar.bz2
rm bwa-mem2-2.2.1_x64-linux.tar.bz2
cd bwa-mem2-2.2.1_x64-linux
sudo chmod u+x bwa*
sudo mv bwa* /usr/local/bin
cd ..
rm -rf bwa-mem2-2.2.1_x64-linux


wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2
tar xvfj htslib-1.12.tar.bz2
rm htslib-1.12.tar.bz2
cd htslib-1.12
autoreconf
./configure
make
sudo make install
cd ..
rm -rf htslib-1.12


wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2
tar xvfj samtools-1.12.tar.bz2
rm samtools-1.12.tar.bz2
cd samtools-1.12
autoreconf
./configure
make
sudo make install
cd ..
rm -rf samtools-1.12


wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2
tar xvfj bcftools-1.12.tar.bz2
rm bcftools-1.12.tar.bz2
cd bcftools-1.12
autoreconf
./configure
make
sudo make install
cd ..
rm -rf bcftools-1.12


wget https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip
unzip gatk-4.2.0.0.zip
rm gatk-4.2.0.0.zip
cd gatk-4.2.0.0
sudo cp gatk /usr/local/bin
sudo cp gatk-package-4.2.0.0-local.jar /usr/local/bin
cd ..
rm -rf gatk-4.2.0.0


# Not the latest version, latest version fails ?intermediate version
wget https://github.com/dkoboldt/varscan/releases/download/2.4.2/VarScan.v2.4.2.jar
sudo mv VarScan.v2.4.2.jar /usr/local/bin
echo '#!/usr/bin/env bash
java -jar /usr/local/bin/VarScan.v2.4.2.jar "$@"' >varscan
chmod +x varscan
sudo mv varscan /usr/local/bin

sudo bash -c "$(curl -L https://basemount.basespace.illumina.com/install)"


wget https://github.com/eawilson/pipeline/archive/refs/tags/v1.1.0.tar.gz -O pipeline-1.1.0.tar.gz
tar xvzf pipeline-1.1.0.tar.gz
rm pipeline-1.1.0.tar.gz
cd pipeline-1.1.0
sudo python3 setup.py develop
cd ..


wget https://github.com/eawilson/elduderino/archive/refs/tags/v1.2.0.tar.gz -O elduderino-1.2.0.tar.gz
tar xvzf elduderino-1.2.0.tar.gz
rm elduderino-1.2.0.tar.gz
cd elduderino-1.2.0
make
sudo make install
cd ..


wget https://github.com/eawilson/udini/archive/refs/tags/v1.0.0.tar.gz -O udini-1.0.0.tar.gz
tar xvzf udini-1.0.0.tar.gz
rm udini-1.0.0.tar.gz
cd udini-1.0.0
make
sudo make install
cd ..


wget https://github.com/eawilson/trim/archive/refs/tags/v1.0.0.tar.gz -O trim-1.0.0.tar.gz
tar xvzf trim-1.0.0.tar.gz
rm trim-1.0.0.tar.gz
cd trim-1.0.0
make
sudo make install
cd ..


wget https://github.com/eawilson/CoverMi/archive/refs/tags/v6.0.0.tar.gz -O CoverMi-6.0.0.tar.gz
tar xvzf CoverMi-6.0.0.tar.gz
rm CoverMi-6.0.0.tar.gz
cd CoverMi-6.0.0
sudo python3 setup.py develop
cd ..


wget https://github.com/eawilson/downsample_sam/archive/refs/tags/v1.0.0.tar.gz -O downsample_sam-1.0.0.tar.gz
tar xvzf downsample_sam-1.0.0.tar.gz
rm downsample_sam-1.0.0.tar.gz
cd downsample_sam-1.0.0
make
sudo make install
cd ..


wget https://sourceforge.net/projects/fastuniq/files/FastUniq-1.1.tar.gz/download -O FastUniq-1.1.tar.gz
tar xvzf FastUniq-1.1.tar.gz
rm FastUniq-1.1.tar.gz
cd FastUniq/source
make
sudo mv fastuniq /usr/local/bin
cd ../..
rm -r FastUniq


VARDICT=VarDict-1.8.2
wget "https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/v1.8.2/$VARDICT.tar"
tar -xf $VARDICT.tar
rm $VARDICT.tar
sudo mv $VARDICT/bin/var2vcf_valid.pl /usr/local/bin
echo '#!/usr/bin/env bash
APP_HOME=/usr/local/bin/VarDict-1.8.2
CLASSPATH=$APP_HOME/lib/VarDict-1.8.2.jar:$APP_HOME/lib/commons-cli-1.2.jar:$APP_HOME/lib/commons-math3-3.6.1.jar:$APP_HOME/lib/jregex-1.2_01.jar:$APP_HOME/lib/htsjdk-2.21.1.jar
java -Xms768m -Xmx8g -classpath "$CLASSPATH" com.astrazeneca.vardict.Main -c 1 -S 2 -E 3 -g 4 "$@"' >vardictjava
chmod +x vardictjava
sudo mv vardictjava /usr/local/bin
sudo mv $VARDICT /usr/local/bin


sudo mkdir /usr/local/lib/site_perl
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
sudo perl INSTALL.pl --AUTO a -d /usr/local/lib/site_perl
sudo cp -r  modules/Bio/EnsEMBL /usr/local/lib/site_perl/Bio
sudo mv vep /usr/local/bin
cd ..
sudo rm -rf ensembl-vep


wget https://github.com/Ensembl/ensembl-xs/archive/refs/tags/2.3.2.tar.gz
tar xvfz 2.3.2.tar.gz
rm 2.3.2.tar.gz
cd ensembl-xs-2.3.2
perl Makefile.PL
make
make test
sudo make install
cd ..
rm -rf ensembl-xs-2.3.2


#############################################################################################################################
# Data files                                                                                                                #
#############################################################################################################################

# # https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
# wget "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
# gunzip hs37d5.fa.gz
# mkdir GRCh37_EBV_HPV_bwa_mem2
# chr_rename_fasta hs37d5.fa >GRCh37_EBV_HPV_bwa_mem2/GRCh37_EBV_HPV.fna
# rm hs37d5.fa
# cat EBV_HPV.fna >>GRCh37_EBV_HPV_bwa_mem2/GRCh37_EBV_HPV.fna
# bwa-mem2 index GRCh37_EBV_HPV_bwa_mem2/GRCh37_EBV_HPV.fna
# samtools faidx GRCh37_EBV_HPV_bwa_mem2/GRCh37_EBV_HPV.fna
# gatk CreateSequenceDictionary R=GRCh37_EBV_HPV_bwa_mem2/GRCh37_EBV_HPV.fna
# tar -cvzf GRCh37_EBV_HPV_bwa_mem2.tar.gz GRCh37_EBV_HPV_bwa_mem2
# 
# 
# wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
# gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
# mkdir GRCh38_EBV_HPV_bwa_mem2
# cp GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38_EBV_HPV_bwa_mem2/GRCh38_EBV_HPV.fna
# rm GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
# grep "HPV16" -A 10000 EBV_HPV.fna >>GRCh38_EBV_HPV_bwa_mem2/GRCh38_EBV_HPV.fna
# bwa-mem2 index GRCh38_EBV_HPV_bwa_mem2/GRCh38_EBV_HPV.fna
# samtools faidx GRCh38_EBV_HPV_bwa_mem2/GRCh38_EBV_HPV.fna
# gatk CreateSequenceDictionary R=GRCh38_EBV_HPV_bwa_mem2/GRCh38_EBV_HPV.fna
# tar -cvzf GRCh38_EBV_HPV_bwa_mem2.tar.gz GRCh38_EBV_HPV_bwa_mem2
# 
# 
# mkdir vep
# perl ~/ensembl-vep/INSTALL.pl --AUTO c --ASSEMBLY GRCh37 --SPECIES homo_sapiens_refseq --CACHEDIR ~/ephemoral/vep






















