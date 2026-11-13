# sherlock-hdl

HDL 1.4.3 (commit `e94ba00`, GPL-3.0) for nf-core/sherlock.
Upstream: https://github.com/zhenin/HDL (R package in the `HDL/` subfolder)

Used for global genetic correlation (rg) next to LDSC. Reference panels are not in the image: pass the
unpacked panel folder (e.g. `UKB_imputed_SVD_eigen99_extraction/`) to the pipeline instead.

Build (from repo root):

```bash
docker build --platform linux/amd64 -t sherlock-hdl:1.4.3-e94ba00 containers/hdl
```

Test:

```bash
docker run --rm sherlock-hdl:1.4.3-e94ba00 Rscript -e 'library(HDL); library(optparse); library(R.utils)'
```
