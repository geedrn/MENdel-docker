FROM continuumio/miniconda3

MAINTAINER Ryo NIWA

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev 

RUN conda create -n MENdel python=3.7

# install conda package
SHELL ["conda", "run", "-n", "MENdel", "/bin/bash", "-c"]
RUN conda install -c r r-essentials
RUN conda install -c r r
RUN pip install argparse 
RUN pip install pandas 
RUN pip install Bio 
RUN pip install numpy 
RUN pip install scipy 
RUN conda install -c bioconda r-rentrez

# Setting
RUN mkdir /MENdel_root
RUN cd /MENdel_root && \
    git clone https://github.com/shendurelab/Lindel.git && \
    git clone https://github.com/FriedbergLab/MENTHU-command-line.git && \
    apt-get update --allow-releaseinfo-change && \
    git clone https://github.com/FriedbergLab/MENdel-command-line.git

RUN Rscript -e "install.packages(c('xml2','rhandsontable', 'plyr', 'stringr', 'stringi', 'rlist', 'DT', 'devtools', 'curl', 'plyr', 'jsonlite', 'httr'), repos='https://cloud.r-project.org/')" 

RUN Rscript -e "install.packages(c('BiocManager'), repos='https://cloud.r-project.org/')" -e "BiocManager::install(c('Biostrings'))"

RUN echo "conda activate MENdel" >> ~/.bashrc
RUN mkdir /MENdel_root/MENdel_Output
WORKDIR /MENdel_root/MENdel-command-line/

