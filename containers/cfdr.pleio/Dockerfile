FROM mambaorg/micromamba:1.5.10

# cfdr.pleio 0.0.0.9100, latest commit on master (repo has no tags)
ARG CFDRPLEIO_COMMIT=76d5085e6d3f3ca9576d5d7564d2acf11bcfd021

LABEL org.opencontainers.image.title="sherlock-cfdr.pleio" \
      org.opencontainers.image.description="cfdr.pleio (conjFDR) for nf-core/sherlock" \
      org.opencontainers.image.source="https://github.com/guillermocomesanacimadevila/sherlock" \
      org.opencontainers.image.licenses="GPL-3.0-only" \
      org.opencontainers.image.version="0.0.0.9100-76d5085"

USER root

# install dependencies + R deps from conda-forge (binom from CRAN below; cfdr.pleio is pure R -> no compiler needed)
RUN micromamba install -y -n base -c conda-forge \
        r-base=4.4 \
        r-r.utils \
        procps-ng \
        curl \
        r-data.table r-matrix r-matrixstats r-r6 \
        r-interp r-optparse \
    && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH

# r-binom has no linux-aarch64 build on conda-forge -> install from CRAN (pure R, no compiler)
RUN Rscript -e 'install.packages("binom", repos = "https://cloud.r-project.org")' \
    && Rscript -e 'library(binom)'

# cfdr.pleio at the pinned commit (package is the repo root)
RUN curl -sL https://github.com/alexploner/cfdr.pleio/archive/${CFDRPLEIO_COMMIT}.tar.gz | tar -xz -C /tmp \
    && R CMD INSTALL /tmp/cfdr.pleio-${CFDRPLEIO_COMMIT} \
    && rm -rf /tmp/cfdr.pleio-*

# verify!!!!!!
RUN Rscript -e 'library(cfdr.pleio); library(optparse); library(data.table); library(R.utils)'