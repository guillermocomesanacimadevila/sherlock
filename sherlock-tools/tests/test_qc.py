import polars as pl 

from sherlock_tools.preprocessing import QualityControl

import pytest


# python3.12 -m pytest -q
# First bloody unit test
# generate toy df and dat
# one function to test each thing
# final funct to run

COLS = ["MarkerName", "Chrom", "Pos", "EA", "NEA", "Beta", "StdErr", "Pval", "Freq", "Info"]

ROWS = [
    # SNP     CHR     POS         A1    A2   BETA    SE      P       FRQ      INFO
    ("rs1",  "1",    "1000",     "a",  "g", "0.1",  "0.01", "1e-5", "0.30",  "0.95"),  # keep (lowercase alleles)
    ("rs2",  "chr2", "2000",     "A",  "C", "-0.2", "0.02", "0",    "0.40",  "0.90"),  # keep (chr prefix, P == 0, flip)
    ("rs3",  "3",    "3000",     "AT", "G", "0.1",  "0.01", "0.5",  "0.30",  "0.95"),  # indel
    ("rs4",  "4",    "4000",     "A",  "T", "0.1",  "0.01", "0.5",  "0.30",  "0.95"),  # palindrome
    ("rs5",  "5",    "5000",     "A",  "G", "0.1",  "0.01", "0.5",  "0.30",  "0.95"),  # duplicate ID
    ("rs5",  "5",    "5001",     "A",  "G", "0.1",  "0.01", "0.5",  "0.30",  "0.95"),  # duplicate ID
    ("rs6",  "7",    "7000",     "A",  "G", "0.1",  "0.01", "0.5",  "0.30",  "0.95"),  # multiallelic
    ("rs7",  "7",    "7000",     "A",  "C", "0.1",  "0.01", "0.5",  "0.30",  "0.95"),  # multiallelic
    ("rs8",  "X",    "8000",     "A",  "G", "0.1",  "0.01", "0.5",  "0.30",  "0.95"),  # non-autosomal
    ("rs9",  "9",    "9000",     "A",  "G", "0.1",  "0.01", "0.5",  "0.005", "0.95"),  # low MAF
    ("rs10", "10",   "10000",    "A",  "G", "0.1",  "0.01", "0.5",  "0.995", "0.95"),  # low MAF (FRQ > 0.5 side)
    ("rs11", "11",   "11000",    "A",  "G", "0.1",  "0.01", "0.5",  "0.30",  "0.50"),  # low INFO
    ("rs12", "12",   "12000",    "A",  "G", None,   "0.01", "0.5",  "0.30",  "0.95"),  # missing BETA
    ("rs13", "6",    "30000000", "A",  "G", "0.1",  "0.01", "0.5",  "0.30",  "0.95"),  # MHC (GRCh38)
    ("rs14", "19",   "44908684", "C",  "T", "0.1",  "0.01", "0.5",  "0.30",  "0.95"),  # APOE (GRCh38)
    ("rs15", "15",   "15000",    "C",  "T", "0.1",  "0.01", "0.5",  "0.30",  "0.95"),  # keep
]

@pytest.fixture
def toy_df():
    return pl.DataFrame(ROWS, schema=COLS, orient="row")

# print(toy_df())

@pytest.fixture
def make_qc(toy_df):
    # factory -> each test can override any arg
    def _make(**kwargs):
        args = dict(
            sumstats=toy_df,
            snp_col="MarkerName",
            chr_col="Chrom",
            pos_col="Pos",
            a1_col="EA",
            a2_col="NEA",
            beta_col="Beta",
            se_col="StdErr",
            pval_col="Pval",
            maf_thresh=0.01,
            genome_build="GRCh38",
            target_genome_build="GRCh38",
            trait_type="cc",
            n_cases=1000,
            n_controls=3000,
            freq_col="Freq",
            info_col="Info",
            info_thresh=0.8,
            drop_mhc=True,
            drop_apoe=True,
        )
        args.update(kwargs)
        return QualityControl(**args)
    return _make

def removed(before, after, col="MarkerName"):
    return set(before[col].to_list()) - set(after[col].to_list())

# one test per method!
def test_fix_data_types(make_qc):
    qc = make_qc()
    df = qc._fix_data_types()
    assert df["Pos"].dtype == pl.Int64
    assert df["Beta"].dtype == pl.Float64
    assert df["Freq"].dtype == pl.Float64
    assert df.filter(pl.col("MarkerName") == "rs2")["Chrom"].item() == "2"
    assert df.filter(pl.col("MarkerName") == "rs1")["EA"].item() == "A"

def test_fix_pval(make_qc):
    qc = make_qc()
    qc._fix_data_types()
    df = qc._fix_pval()
    assert df.filter(pl.col("MarkerName") == "rs2")["Pval"].item() == 1e-300
    assert (df["Pval"] > 0).all()

def test_correct_for_missingness(make_qc):
    qc = make_qc()
    before = qc._fix_data_types()
    after = qc._correct_for_missingness()
    assert removed(before, after) == {"rs12"}