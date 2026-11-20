# sherlock-pleio

PLEIO (commit `b18d507`) for nf-core/sherlock, with LDSC (commit `2fdeeb3`) for its `ldsc_preprocess.py` step.
Upstream: https://github.com/cuelee/pleio (old address hanlab-SNU/PLEIO now redirects here);
LDSC: https://github.com/bulik/ldsc

Used for PLEIO: a per-SNP test of whether a SNP affects more than one trait, in the pleiotropy step next to conjFDR.
`ldsc_preprocess.py` estimates the genetic covariance and environmental correlation between the traits with LDSC and writes
PLEIO's inputs (`metain`, `sg`, `ce`); `pleio.py` then runs the test (plus BLUP effect sizes with `--blup`).
The LD scores (`eur_w_ld_chr`) are not in the image: the pipeline passes the same `ldsc_ref` folder used by the rg step.

The image has two environments, because LDSC is Python 2.7 and PLEIO is Python 3:

| Command            | Runs                              | Environment                               |
| ------------------ | --------------------------------- | ----------------------------------------- |
| `pleio_preprocess` | `/opt/pleio/ldsc_preprocess.py`   | `ldsc`: python 2.7 + bioconda ldsc 1.0.1  |
| `pleio`            | `/opt/pleio/pleio.py`             | `pleio`: python 3.10, numpy 1.26, pandas 1.5, scipy 1.10 |

One change to upstream code: `ldsc_preprocess.py` creates `ldsc/__init__.py` next to itself on every run, which fails on
read-only images (Apptainer/Singularity). The file is created at build time and that line is replaced with `pass`.

Licence: none. The repo says "This project has no license currently" (files: Copyright (C) 2018 Cue Hyunkyu Lee), so the
image label is `NOASSERTION`. LDSC is GPL-3. All credit goes to the PLEIO and LDSC authors; please cite
Lee et al. 2021, AJHG 108(1):36-48 (https://doi.org/10.1016/j.ajhg.2020.11.017) and Bulik-Sullivan et al. 2015 when using it.

Build (from repo root; amd64 only; Python 2.7 is likely not available for linux-aarch64 on conda-forge):

```bash
docker build --platform linux/amd64 -t sherlock-pleio:0.0.0-b18d507 containers/pleio
```

Test:

```bash
docker run --rm sherlock-pleio:0.0.0-b18d507 pleio_preprocess -h
docker run --rm sherlock-pleio:0.0.0-b18d507 pleio -h
```
