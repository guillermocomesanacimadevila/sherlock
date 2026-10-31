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

def test_restrict_to_autosomal(make_qc):
    qc = make_qc()
    before = qc._fix_data_types()
    after = qc._restrict_to_autosomal()
    assert removed(before, after) == {"rs8"}

def test_remove_duplicates(make_qc):
    qc = make_qc()
    before = qc._fix_data_types()
    after = qc._remove_duplicates()
    assert removed(before, after) == {"rs5"}

def test_remove_multiallelic(make_qc):
    qc = make_qc()
    before = qc._fix_data_types()
    after = qc._remove_multiallelic()
    assert removed(before, after) == {"rs6", "rs7"}

def test_remove_indels(make_qc):
    qc = make_qc()
    before = qc._fix_data_types()
    after = qc._remove_indels()
    assert removed(before, after) == {"rs3"}

def test_remove_palindromes(make_qc):
    qc = make_qc()
    before = qc._fix_data_types()
    after = qc._remove_palindromes()
    assert removed(before, after) == {"rs4"}

def test_maf_filter(make_qc):
    qc = make_qc()
    before = qc._fix_data_types()
    after = qc._maf_filter()
    assert removed(before, after) == {"rs9", "rs10"}

def test_maf_filter_skips_without_freq(make_qc):
    qc = make_qc(freq_col=None)
    before = qc._fix_data_types()
    after = qc._maf_filter()
    assert after.height == before.height

def test_info_filter(make_qc):
    qc = make_qc()
    before = qc._fix_data_types()
    after = qc._info_filter()
    assert removed(before, after) == {"rs11"}

def test_remove_mhc(make_qc):
    qc = make_qc()
    before = qc._fix_data_types()
    after = qc._remove_mhc()
    assert removed(before, after) == {"rs13"}

def test_keep_mhc(make_qc):
    qc = make_qc(drop_mhc=False)
    before = qc._fix_data_types()
    after = qc._remove_mhc()
    assert after.height == before.height

def test_remove_apoe(make_qc):
    qc = make_qc()
    before = qc._fix_data_types()
    after = qc._remove_apoe()
    assert removed(before, after) == {"rs14"}

def test_keep_apoe(make_qc):
    qc = make_qc(drop_apoe=False)
    before = qc._fix_data_types()
    after = qc._remove_apoe()
    assert after.height == before.height

def test_align_betas_to_risk_allele(make_qc):
    qc = make_qc()
    qc._fix_data_types()
    df = qc._align_betas_to_risk_allele()
    rs2 = df.filter(pl.col("MarkerName") == "rs2").row(0, named=True)
    assert rs2["Beta"] == 0.2
    assert (rs2["EA"], rs2["NEA"]) == ("C", "A")
    assert rs2["Freq"] == pytest.approx(0.6)
    assert (df["Beta"].drop_nulls() >= 0).all()

# sample size
def test_effective_n_cc(make_qc):
    qc = make_qc()
    df = qc._impute_effective_n()
    assert (df["N"] == 3000).all()  # 4 * 1000 * 3000 / 4000

def test_effective_n_quant(make_qc):
    qc = make_qc(trait_type="quant", n=5000, n_cases=None, n_controls=None)
    df = qc._impute_effective_n()
    assert (df["N"] == 5000).all()

def test_effective_n_missing_raises(make_qc):
    qc = make_qc(trait_type="quant", n=None)
    with pytest.raises(ValueError):
        qc._impute_effective_n()

# liftover (fake chain file -> no real file needed)
class FakeChainFile:
    # chr1 -> +100 bp, chr2 -> no hit, chr15 -> lands on chr16
    def __init__(self, path, one_based=True):
        pass

    def __getitem__(self, chrom):
        def lookup(pos):
            if chrom == "chr1":
                return [("chr1", pos + 100, "+")]
            if chrom == "chr15":
                return [("chr16", pos, "+")]
            return []
        return _Lookup(lookup)

class _Lookup:
    def __init__(self, fn):
        self.fn = fn

    def __getitem__(self, pos):
        return self.fn(pos)

def test_liftover_skips_same_build(make_qc):
    qc = make_qc()
    before = qc._fix_data_types()
    after = qc._liftover()
    assert after.equals(before)

def test_liftover(make_qc, monkeypatch):
    monkeypatch.setattr("sherlock_tools.preprocessing.qc.ChainFile", FakeChainFile)
    qc = make_qc(genome_build="GRCh37", chain_file="fake.chain")
    qc._fix_data_types()
    qc.sumstats = qc.sumstats.filter(pl.col("MarkerName").is_in(["rs1", "rs2", "rs15"]))
    df = qc._liftover()
    assert df["MarkerName"].to_list() == ["rs1"]  # rs2 no hit, rs15 changed chr
    assert df["Pos"].item() == 1100

def test_liftover_missing_chain_raises(make_qc):
    qc = make_qc(genome_build="GRCh37", chain_file=None)
    with pytest.raises(ValueError):
        qc._liftover()

# full run
def test_run(make_qc):
    df = make_qc().run()
    assert df.columns == ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N", "FRQ", "INFO"]
    assert set(df["SNP"].to_list()) == {"rs1", "rs2", "rs15"}
    assert (df["BETA"] >= 0).all()
    assert (df["P"] > 0).all()
    assert (df["N"] == 3000).all()