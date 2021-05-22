FROM public.ecr.aws/lts/ubuntu:18.04_stable

RUN apt-get -y update && \
    apt-get -y dist-upgrade && \
    apt-get -y install --no-install-recommends python3-minimal python3-pip python3-tk && \
    apt-get -y install --no-install-recommends build-essential zlib1g-dev libbz2-dev liblzma-dev autoconf automake && \
    apt-get -y install --no-install-recommends openjdk-8-jdk libmysqlclient-dev wget curl && \
    cpan App::cpanminus && \
    cpanm DBD::mysql Archive::Zip Archive::Extract Try::Tiny JSON Module::Build && \
    mkdir /usr/local/lib/site_perl && \
    wget https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/101.0.tar.gz && \
    tar -xzf 101.0.tar.gz && \
    rm 101.0.tar.gz && \
    cd ensembl-vep-release-101.0 && \
    perl INSTALL.pl --NO_UPDATE --AUTO a -d /usr/local/lib/site_perl && \
    cp -r  modules/Bio/EnsEMBL /usr/local/lib/site_perl/Bio && \
    mv vep /usr/local/bin && \
    cd .. && \
    rm -rf ensembl-vep-release-101.0 && \
    pip3 install --no-cache-dir setuptools && \
    pip3 install --no-cache-dir boto3 matplotlib && \
    wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
    tar -xjf samtools-1.10.tar.bz2 && \
    rm samtools-1.10.tar.bz2 && \
    cd samtools-1.10 && \
    autoreconf && \
    ./configure --without-curses && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.10 && \
    wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar -xjf bwa-0.7.17.tar.bz2 && \
    rm bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    mv bwa /usr/local/bin && \
    cd .. && \
    rm -rf bwa-0.7.17 && \
    wget -O VarScan.v2.3.9.jar https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar/download && \
    mv VarScan.v2.3.9.jar /usr/local/bin && \
    echo '#!/usr/bin/env bash' >varscan && \
    echo 'java -jar /usr/local/bin/VarScan.v2.3.9.jar "$@"' >>varscan && \
    chmod +x varscan && \
    mv varscan /usr/local/bin && \
    wget 'https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/v1.8.2/VarDict-1.8.2.tar' && \
    tar -xf VarDict-1.8.2.tar && \
    rm VarDict-1.8.2.tar && \
    mv VarDict-1.8.2/bin/var2vcf_valid.pl /usr/local/bin && \
    echo '#!/usr/bin/env bash' >vardictjava && \
    echo 'APP_HOME=/usr/local/bin/VarDict-1.8.2' >>vardictjava && \
    echo 'CLASSPATH=$APP_HOME/lib/VarDict-1.8.2.jar:$APP_HOME/lib/commons-cli-1.2.jar:$APP_HOME/lib/commons-math3-3.6.1.jar:$APP_HOME/lib/jregex-1.2_01.jar:$APP_HOME/lib/htsjdk-2.21.1.jar' >>vardictjava && \
    echo 'java -Xms768m -Xmx8g -classpath "$CLASSPATH" com.astrazeneca.vardict.Main -c 1 -S 2 -E 3 -g 4 "$@"' >>vardictjava && \
    chmod +x vardictjava && \
    mv vardictjava /usr/local/bin && \
    mv VarDict-1.8.2 /usr/local/bin && \
    git clone git@github.com:eawilson/pipeline.git && \
    cd pipeline && \
    python3 setup.py develop && \
    cd .. && \
    git clone https://github.com/eawilson/covermi.git && \
    cd covermi && \
    python3 setup.py install && \
    cd ..

