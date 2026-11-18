# sherlock-lava

LAVA 0.1.5 (commit `e729a24`) for nf-core/sherlock.
Upstream: https://github.com/josefin-werme/LAVA (Werme et al. 2022, Nat Genet)

Used for local genetic correlation (rg) per locus, next to HDL-L. Reference data and locus files are not
in the image: pass the LAVA LD reference (UK Biobank recommended over 1000G for EUR) and the locus file
to the pipeline instead.

Licence: upstream DESCRIPTION says "All rights reserved" (no open licence), so the image is labelled
`NOASSERTION`. All credit goes to the LAVA authors; please cite their paper when using it.

Build (from repo root):

```bash
docker build --platform linux/amd64 -t sherlock-lava:0.1.5-e729a24 containers/lava
```

Test:

```bash
docker run --rm sherlock-lava:0.1.5-e729a24 Rscript -e 'library(LAVA); library(matrixsampling); library(optparse); library(R.utils)'
```
