FROM mambaorg/micromamba:1.5.10

# PLEIO, latest commit on master (repo has no tags, no licence)
ARG PLEIO_COMMIT=b18d507533ae39dcb6284573d3b36b91cdc13603
# LDSC, needed by ldsc_preprocess.py as a folder next to it (python 2.7)
ARG LDSC_COMMIT=2fdeeb3b44379408794154993dbd6101b8946b7e

LABEL org.opencontainers.image.title="sherlock-pleio" \
      org.opencontainers.image.description="PLEIO (+ LDSC for ldsc_preprocess) for nf-core/sherlock" \
      org.opencontainers.image.source="https://github.com/guillermocomesanacimadevila/sherlock" \
      org.opencontainers.image.licenses="NOASSERTION" \
      org.opencontainers.image.version="0.0.0-b18d507"

USER root

# tools for Nextflow (ps for metrics) + git to clone
RUN micromamba install -y -n base -c conda-forge \
        procps-ng \
        git \
        gzip \
    && micromamba clean -a -y

ENV PATH=/opt/conda/envs/pleio/bin:/opt/conda/bin:$PATH

# fix python version for LDSC preprocessing part
RUN micromamba create -y -n ldsc -c conda-forge -c bioconda \
        python=2.7 \
        ldsc=1.0.1 \
    && micromamba clean -a -y 

# fix python script to PY >=3.7 for pleio side of things
RUN micromamba create -y -n pleio -c conda-forge \
        python=3.10 \
        numpy=1.26 \
        pandas=1.5 \
        scipy=1.10 \
        polars=1.43.0 \
    && micromamba clean -a -y

# clone pleio repo 
RUN git clone https://github.com/cuelee/pleio.git /opt/pleio \
    && git -C /opt/pleio checkout ${PLEIO_COMMIT} \
    && rm -rf /opt/pleio/.git /opt/pleio/.ipynb_checkpoints /opt/pleio/__pycache__

# LDSC at the pinned commit, inside the pleio folder (where ldsc_preprocess.py looks for it)
RUN git clone https://github.com/bulik/ldsc.git /opt/pleio/ldsc \
    && git -C /opt/pleio/ldsc checkout ${LDSC_COMMIT} \
    && rm -rf /opt/pleio/ldsc/.git

# ldsc_preprocess.py writes ldsc/__init__.py on every run -> fails on read-only images (apptainer)
# so create it now and remove that line
RUN touch /opt/pleio/ldsc/__init__.py \
    && sed -i "s|^\(\s*\)open(ldsc_dir+'/__init__.py', 'a').close()|\1pass|" /opt/pleio/ldsc_preprocess.py

# wrappers so the module doesn't need to know about the two envs
RUN printf '#!/bin/bash\nexec /opt/conda/envs/ldsc/bin/python /opt/pleio/ldsc_preprocess.py "$@"\n' > /usr/local/bin/pleio_preprocess \
    && printf '#!/bin/bash\nexec /opt/conda/envs/pleio/bin/python /opt/pleio/pleio.py "$@"\n' > /usr/local/bin/pleio \
    && chmod +x /usr/local/bin/pleio_preprocess /usr/local/bin/pleio

# verify!!!!!!
RUN grep -q "^\s*pass$" /opt/pleio/ldsc_preprocess.py \
    && pleio_preprocess -h > /dev/null \
    && pleio -h > /dev/null \
    && python3 -c "import polars"
