# sherlock-cfdr.pleio

cfdr.pleio 0.0.0.9100 (commit `76d5085`) for nf-core/sherlock.
Upstream: https://github.com/alexploner/cfdr.pleio (R re-implementation of pleiofdr, https://github.com/precimed/pleiofdr;
method: Andreassen et al. 2013, PLoS Genet)

Used for conjunctional FDR (conjFDR): SNPs associated with both traits of a pair, in the pleiotropy step next to PLEIO.
The LD reference is not in the image: download it once and pass it to the pipeline instead.

Licence: GPL-3 (upstream DESCRIPTION). All credit goes to the cfdr.pleio and pleiofdr authors; please cite
Andreassen et al. 2013 when using it.

Build (from repo root):

```bash
docker build --platform linux/amd64 -t sherlock-cfdr.pleio:0.0.0.9100-76d5085 containers/cfdr.pleio
```

Test:

```bash
docker run --rm sherlock-cfdr.pleio:0.0.0.9100-76d5085 Rscript -e 'library(cfdr.pleio); library(optparse); library(data.table)'
```
