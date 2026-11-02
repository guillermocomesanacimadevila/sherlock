# sherlock-genomicsem

GenomicSEM 0.0.5 (commit `da95d43`, GPL-3.0) for nf-core/sherlock.
Upstream: https://github.com/GenomicSEM/GenomicSEM

Build (from repo root):

```bash
docker build --platform linux/amd64 -t sherlock-genomicsem:0.0.5-da95d43 containers/genomicsem
```

Test:

```bash
docker run --rm sherlock-genomicsem:0.0.5-da95d43 Rscript -e 'library(GenomicSEM)'
```
