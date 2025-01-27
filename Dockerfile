
ARG UBUNTU_VER=22.04
FROM ubuntu:${UBUNTU_VER}

# Install Python and other dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    wget \
    unzip \
    python3.10 \
    python3.10-venv \
    python3-pip \
    default-jre \
    libxml2 \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    pandoc \
    pandoc-citeproc \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev

# Install required Python packages
RUN python3.10 -m pip install --upgrade pip && \
    python3.10 -m pip install \
    numpy==1.23.1 \
    pandas==2.0.3 \
    natsort==8.3.1

# Install last Sqlite version
RUN wget https://www.sqlite.org/2024/sqlite-autoconf-3460000.tar.gz && \
    tar xvfz sqlite-autoconf-3460000.tar.gz && \
    cd sqlite-autoconf-3460000 && \
    export CFLAGS="-DSQLITE_ENABLE_FTS3 \
    -DSQLITE_ENABLE_FTS3_PARENTHESIS \
    -DSQLITE_ENABLE_FTS4 \
    -DSQLITE_ENABLE_FTS5 \
    -DSQLITE_ENABLE_JSON1 \
    -DSQLITE_ENABLE_LOAD_EXTENSION \
    -DSQLITE_ENABLE_RTREE \
    -DSQLITE_ENABLE_STAT4 \
    -DSQLITE_ENABLE_UPDATE_DELETE_LIMIT \
    -DSQLITE_SOUNDEX \
    -DSQLITE_TEMP_STORE=3 \
    -DSQLITE_USE_URI \
    -O2 \
    -fPIC" && \
    export PREFIX="/usr/local" && \
    LIBS="-lm" ./configure --disable-tcl --enable-shared --enable-tempstore=always --prefix="$PREFIX" && \
    make && \
    make install && \
    cd .libs && \
    mkdir backup && \
    cp /lib/x86_64-linux-gnu/libsqlite3.so.0 ./backup/libsqlite3.so.0 && \
    cp libsqlite3.so.0 /lib/x86_64-linux-gnu/libsqlite3.so.0 && \
    cd ../../


# Download FASTQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod 755 FastQC/fastqc && \
    mv FastQC /usr/local/bin/ && \
    ln -s /usr/local/bin/FastQC/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.11.9.zip 

# Download fastp
RUN wget -O /usr/local/bin/fastp http://opengene.org/fastp/fastp.0.23.2 && \
    chmod a+x /usr/local/bin/fastp

# Descargar e instalar Bowtie
RUN wget -c https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.1/bowtie-1.3.1-linux-x86_64.zip && \
    unzip bowtie-1.3.1-linux-x86_64.zip && \
    rm bowtie-1.3.1-linux-x86_64.zip
RUN cp -p bowtie-1.3.1-linux-x86_64/bowtie* /usr/local/bin

# Download multiqc
RUN pip3 install multiqc==1.22.2

# SRA-TOOLS (prefetch, fasterq-dump)
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz && \
    tar -vxzf sratoolkit.3.1.1-ubuntu64.tar.gz && \
    mv sratoolkit.3.1.1-ubuntu64 /usr/local/sra-tools

# Set up the PATH for SRA-TOOLS
ENV PATH /usr/local/sra-tools/bin:$PATH

# Install R
ARG r=4.2
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Madrid
RUN apt-get update -y && \
    apt-get install -y \
        dirmngr \
        apt-transport-https \
        software-properties-common && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
        tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    # NOTE: Only R >= 4.0 is available in this repo
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu '$(lsb_release -cs)'-cran40/' && \
    apt-get install -y \
        r-base=${r}* \
        r-recommended=${r}*

# Install R packages
RUN R -e "install.packages('curl', repos = 'http://cran.us.r-project.org', version='5.2.1')"
RUN R -e "install.packages('httr', repos = 'http://cran.us.r-project.org', version='1.4.7')"
RUN R -e "install.packages('usethis', repos = 'http://cran.us.r-project.org', version='2.2.3')"
RUN R -e "install.packages('renv', repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('devtools', repos = 'http://cran.us.r-project.org', version='2.4.5')"
RUN R -e "install.packages('BiocManager', repos='http://cran.us.r-project.org', version='1.30.23')"
RUN R -e "BiocManager::install('GO.db')"
RUN R -e "BiocManager::install('apeglm')"
RUN R -e "options(renv.settings.bioconductor.version = '3.19')"
RUN R -e "renv::install('bioc::DESeq2@1.44.0')"
RUN R -e "renv::install('bioc::genefilter@1.86.0')"
RUN R -e "renv::install('bioc::ComplexHeatmap@2.20.0')"
RUN Rscript -e 'devtools::install_github("PF2-pasteur-fr/SARTools")'
RUN R -e "install.packages('dendextend', repos = 'http://cran.us.r-project.org', version='1.17.1')"
RUN R -e "install.packages('ff', repos = 'http://cran.us.r-project.org', version='4.0.12')"
RUN R -e "install.packages('tidyverse', repos = 'http://cran.us.r-project.org', version='2.0.0')"
RUN R -e "install.packages('Ckmeans.1d.dp', repos = 'http://cran.us.r-project.org', version='4.3.5')"
RUN R -e "install.packages('ComplexUpset', repos = 'http://cran.us.r-project.org', version='1.3.5')"
RUN R -e "install.packages('ggthemes', repos = 'http://cran.us.r-project.org', version='5.1.0')"
RUN R -e "install.packages('argparse', repos = 'http://cran.us.r-project.org', version='2.2.3')"
RUN R -e "install.packages('plotly', repos = 'http://cran.us.r-project.org', version='4.10.3')"
