# sherlock-tools

Python library behind the [nf-core/sherlock](https://github.com/guillermocomesanacimadevila/sherlock) pipeline,
which looks for shared genetic signals between two traits and narrows them down to tissue, cell type, gene and SNP.

The pipeline's `bin/` scripts import it, but it also works on its own in Python.

## What it does

`QualityControl` cleans and harmonises one set of GWAS summary statistics (built on [polars](https://pola.rs)):

- fixes data types, strips `chr`, keeps autosomes (1-22) by default
- stops with an error if any P-value is outside 0-1 (e.g. a -log10 P column)
- removes missing values and duplicated SNP IDs
- removes indels, palindromic, multi-allelic and non-autosomal variants (each can be switched off:
  `drop_indels`, `drop_palindromes`, `drop_multiallelic`, `drop_non_autosomal`, all default `True`)
- MAF and (optional) INFO filters (INFO is only applied when both `info_col` and `info_thresh` are given)
- liftover between GRCh37 and GRCh38 (UCSC chain file)
- optional removal of the MHC (`drop_mhc`, default `True`) and APOE (`drop_apoe`, default `False`) regions
- aligns betas to the risk allele
- writes one sample size column `N` (effective N for case/control traits)

Output columns follow LDSC naming: `SNP CHR BP A1 A2 BETA SE P N [FRQ] [INFO]`.
`A1` is the risk allele: where BETA was negative, A1/A2 are swapped and BETA flipped (FRQ becomes 1 - FRQ),
so BETA is always positive.

Sample size is used as given and never rescaled. If your file has a half-size effective N column
(e.g. PGC `NEFFDIV2` / `NEFF_HALF`), double it first, or leave `n_col` empty and give `n_cases` + `n_controls`.
`n_col` must be a per-SNP sample size, not a count of studies.

## Install

```bash
pip install sherlock-tools
```

Requires Python >= 3.12.

## Example

```python
import polars as pl
from sherlock_tools.preprocessing import QualityControl

sumstats = pl.read_csv("AD.tsv", separator="\t")

qc = QualityControl(
    sumstats=sumstats,
    snp_col="SNP", chr_col="CHR", pos_col="BP",
    a1_col="A1", a2_col="A2",
    beta_col="BETA", se_col="SE", pval_col="P",
    n_col="N",
    maf_thresh=0.01,
    freq_col="FRQ",
    genome_build="GRCh37",
    target_genome_build="GRCh38",
    chain_file="hg19ToHg38.over.chain",
    trait_type="cc",
)

clean = qc.run()

# rows removed / left after each step
print(pl.DataFrame(qc.summary, schema=["step", "removed", "remaining"], orient="row"))
```

## Development install

```bash
pip install -e "sherlock-tools[dev]"
pytest sherlock-tools/tests
```

## Changes

- **0.1.3**: INFO filter is skipped when no threshold is set (previously removed every row);
  new P-value range check; per-step counts in `qc.summary`;
  indel / palindrome / multiallelic / non-autosomal filters can be switched off.
- **0.1.2**: package metadata and README.
- **0.1.1**: error when builds differ and no chain file is given.

## Licence

MIT
