FROM mambaorg/micromamba:1.5.10

# LAVA v0.1.5, latest commit on main (repo has no tags)
ARG LAVA_COMMIT=e729a245f7b6923967a96804fbf5246eadf2d6c6

LABEL org.opencontainers.image.title="sherlock-lava" \
      org.opencontainers.image.description="LAVA for nf-core/sherlock" \
      org.opencontainers.image.source="https://github.com/guillermocomesanacimadevila/sherlock" \
      org.opencontainers.image.licenses="NOASSERTION" \
      org.opencontainers.image.version="0.1.5-e729a24"

USER root

# install dependencies + all R deps (LAVA has C++ code -> compiler needed)
RUN micromamba install -y -n base -c conda-forge \
        r-base=4.4 \
        procps-ng \
        curl \
        cxx-compiler make \
        r-data.table r-r.utils r-optparse r-cpp11 \
    && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH

# CRAN-only deps (not on conda-forge), pinned via a dated CRAN snapshot
RUN Rscript -e 'install.packages(c("keep", "matrixsampling"), repos = "https://packagemanager.posit.co/cran/2026-11-01")'

# LAVA at the pinned commit (package is the repo root)
RUN curl -sL https://github.com/josefin-werme/LAVA/archive/${LAVA_COMMIT}.tar.gz | tar -xz -C /tmp \
    && R CMD INSTALL /tmp/LAVA-${LAVA_COMMIT} \
    && rm -rf /tmp/LAVA-*

# verify!!!!!!
RUN Rscript -e 'library(LAVA); library(matrixsampling); library(optparse); library(R.utils)'