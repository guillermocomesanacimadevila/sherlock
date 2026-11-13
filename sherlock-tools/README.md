# sherlock-tools

Python library behind the [nf-core/sherlock](https://github.com/guillermocomesanacimadevila/sherlock) pipeline,
which looks for shared genetic signals between two traits and narrows them down to tissue, cell type, gene and SNP.

The pipeline's `bin/` scripts import it, but it also works on its own in Python.

## What it does

`QualityControl` cleans and harmonises one set of GWAS summary statistics (built on [polars](https://pola.rs)):

- fixes data types, strips `chr`, keeps autosomes (1-22)
- removes missing values, indels, palindromic, duplicated and multi-allelic SNPs
- MAF and (optional) INFO filters
- liftover between GRCh37 and GRCh38 (UCSC chain file)
- optional removal of the MHC and APOE regions
- aligns betas to the risk allele
- writes one sample size column `N` (effective N for case/control traits)

Output columns follow LDSC naming: `SNP CHR BP A1 A2 BETA SE P N [FRQ] [INFO]`, with `A1` as the effect allele.

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
```

## Development install

```bash
pip install -e "sherlock-tools[dev]"
pytest sherlock-tools/tests
```

## Licence

MIT
