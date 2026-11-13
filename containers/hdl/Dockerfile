FROM mambaorg/micromamba:1.5.10

# specific HDL commit (repo has no tags)
ARG HDL_COMMIT=e94ba003349bffccc9b0ff23e0ff649bd84075a0

LABEL org.opencontainers.image.title="sherlock-hdl" \
      org.opencontainers.image.description="HDL for nf-core/sherlock" \
      org.opencontainers.image.source="https://github.com/guillermocomesanacimadevila/sherlock" \
      org.opencontainers.image.licenses="GPL-3.0" \
      org.opencontainers.image.version="1.4.3-e94ba00"

USER root

# install dependencies + all R deps
RUN micromamba install -y -n base -c conda-forge \
        r-base=4.4 \
        procps-ng \
        curl \
        r-dplyr r-data.table r-r.utils r-optparse \
    && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH

# HDL R package lives in the HDL/ subfolder of the repo
RUN curl -sL https://github.com/zhenin/HDL/archive/${HDL_COMMIT}.tar.gz | tar -xz -C /tmp \
    && R CMD INSTALL /tmp/HDL-${HDL_COMMIT}/HDL \
    && rm -rf /tmp/HDL-*

# verify!!!!!!
RUN Rscript -e 'library(HDL); library(optparse); library(R.utils)'
