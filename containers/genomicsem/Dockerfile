FROM mambaorg/micromamba:1.5.10

# specific genSEM commit
ARG GENOMICSEM_COMMIT=da95d431a4709693d07eae55e172ade55526948b

LABEL org.opencontainers.image.title="sherlock-genomicsem" \
      org.opencontainers.image.description="GenomicSEM for nf-core/sherlock" \
      org.opencontainers.image.source="https://github.com/guillermocomesanacimadevila/sherlock" \
      org.opencontainers.image.licenses="GPL-3.0" \
      org.opencontainers.image.version="0.0.5-da95d43"

USER root

# install dependencies + all R deps
RUN micromamba install -y -n base -c conda-forge \
        r-base=4.4 \
        procps-ng \
        curl \
        r-plyr r-e1071 r-readr r-gdata r-lavaan \
        r-doparallel r-foreach r-iterators r-data.table \
        r-stringr r-splitstackshape r-matrix r-r.utils \
        r-dplyr r-rcpp r-mgsub r-sfsmisc r-optparse \
    && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH

RUN Rscript -e 'install.packages("simsalapar", repos = "https://cloud.r-project.org")'

RUN curl -sL https://github.com/GenomicSEM/GenomicSEM/archive/${GENOMICSEM_COMMIT}.tar.gz -o /tmp/gsem.tar.gz \
    && R CMD INSTALL /tmp/gsem.tar.gz \
    && rm /tmp/gsem.tar.gz

# verify!!!!!!
RUN Rscript -e 'library(GenomicSEM); library(optparse)'
