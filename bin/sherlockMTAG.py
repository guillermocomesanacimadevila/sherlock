#!/usr/bin/env python3
import argparse
import subprocess
import polars as pl 


def run_mtag(
    trait1_sumstats: str,
    trait2_sumstats: str,
    ldsc_ref: str,
    n_cores: int,
    out: str
    ):

    """ Run MTAG on 2 traits (+ max FDR) -> mtag """

    cmd = [
        "mtag",
        "--sumstats", f"{trait1_sumstats},{trait2_sumstats}",
        "--snp_name", "SNP",
        "--chr_name", "CHR",
        "--bpos_name", "BP",
        "--a1_name", "A1",
        "--a2_name", "A2",
        "--eaf_name", "FRQ",
        "--n_name", "N",
        "--p_name", "P",
        "--z_name", "Z",
        "--ld_ref_panel", ldsc_ref.rstrip("/") + "/",
        "--fdr",
        "--n_approx",
        "--cores", str(n_cores),
        "--out", out,
        "--stream_stdout",
    ]
    subprocess.run(cmd, check=True)
    return f"{out}_trait_1.txt", f"{out}_trait_2.txt", f"{out}.log"
