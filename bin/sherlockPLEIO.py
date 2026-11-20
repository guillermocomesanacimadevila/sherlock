#!/usr/bin/env python3
import argparse
import subprocess
import polars as pl


def establish_type(
        sherlock_type: str
    ):
    if sherlock_type == "cc":
        return "binary"
    elif sherlock_type == "quant":
        return "quantitative"
    else:
        raise ValueError(f"trait_type must be 'cc' or 'quant', got '{sherlock_type}'")


def make_z_files(
        trait1,
        trait2,
        sumstats1: str, # QC file trait 1 (.qc.tsv.gz)
        sumstats2: str, # QC file trait 2
        out: str,
    ):

    """
    Make PLEIO Z files (SNP A1 A2 Z N P), both on trait 1's alleles.
    -----
    QC aligns each trait to its OWN risk allele, so A1 can differ between traits.
    ldsc_preprocess does not flip: it DROPS SNPs whose A1/A2 differ.
    -> flip trait 2 here (swap alleles + Z sign), drop real mismatches.
    """

    t1 = pl.read_csv(sumstats1, separator="\t", columns=["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"])
    t2 = pl.read_csv(sumstats2, separator="\t", columns=["SNP", "A1", "A2", "BETA", "SE", "P", "N"])

    t1 = t1.with_columns((pl.col("BETA") / pl.col("SE")).alias("Z"))
    t2 = t2.with_columns((pl.col("BETA") / pl.col("SE")).alias("Z"))

    # SNPs in both traits
    df = t1.join(t2, on="SNP", how="inner", suffix="_2")
    print(f"SNPs in both traits: {df.height}")

    # same alleles or swapped -> keep; anything else -> drop
    same    = (pl.col("A1_2") == pl.col("A1")) & (pl.col("A2_2") == pl.col("A2"))
    swapped = (pl.col("A1_2") == pl.col("A2")) & (pl.col("A2_2") == pl.col("A1"))

    n_before = df.height
    df = df.filter(same | swapped)
    print(f"Removed allele mismatch: {n_before - df.height}")

    # swapped -> flip trait 2's Z to trait 1's A1
    df = df.with_columns(
        pl.when(swapped).then(-pl.col("Z_2")).otherwise(pl.col("Z_2")).alias("Z_2")
    )

    # SE = 0 etc -> inf/NaN Z
    n_before = df.height
    df = df.filter(pl.col("Z").is_finite() & pl.col("Z_2").is_finite())
    print(f"Removed non-finite Z: {n_before - df.height}")
    print(f"Rows after Z files: {df.height}")

    zfile1 = f"{out}.{trait1}.z.txt"
    zfile2 = f"{out}.{trait2}.z.txt"

    df.select("SNP", "A1", "A2", "Z", "N", "P").write_csv(zfile1, separator="\t")
    df.select(
        "SNP", "A1", "A2",
        pl.col("Z_2").alias("Z"),
        pl.col("N_2").alias("N"),
        pl.col("P_2").alias("P"),
    ).write_csv(zfile2, separator="\t")

    # positions + both traits' Z and P, for the final table
    positions = df.select(
        "SNP", "CHR", "BP", "A1", "A2",
        pl.col("Z").alias("Z1"),
        pl.col("P").alias("P1"),
        pl.col("Z_2").alias("Z2"),
        pl.col("P_2").alias("P2"),
    )

    return zfile1, zfile2, positions


def gen_input_file(
        trait1,
        trait2,
        zfile1: str, # Z file (SNP A1 A2 Z N P), not the QC file
        zfile2: str,
        trait_type1,
        trait_type2,
        out_path: str, # e.g. <out>.pleio_input.txt
        pop_prev1 = None,
        pop_prev2 = None,
    ):

    """
    Generate input file for PLEIO!!!!
    -----
    # FILE	TYPE	SPREV	PPREV	NAME
    # als_A_sumstats_hg38.qc_S_clean.tsv	binary	0.25	8.50E-05	als_A
    # pkd_A_sumstats_hg38.qc_S_cleanNA.tsv	binary	0.061	0.017	pkd_A

    -----
    Types: (link to --trait-type)
    - binary
    - quantitative

    -----
    SPREV = cases / (cases + controls)
    sherlock always uses 0.5: QC writes N as Neff (= a balanced sample)
    """

    type1 = establish_type(str(trait_type1))
    type2 = establish_type(str(trait_type2))

    # binary needs both SPREV + PPREV, else PLEIO gives NaN effect sizes
    if type1 == "binary" and pop_prev1 is None:
        raise ValueError(f"{trait1} is cc: pop_prev is needed for PLEIO")
    if type2 == "binary" and pop_prev2 is None:
        raise ValueError(f"{trait2} is cc: pop_prev is needed for PLEIO")

    # quantitative -> nan for both
    sprev1, pprev1 = None, None
    sprev2, pprev2 = None, None
    if type1 == "binary":
        sprev1, pprev1 = 0.5, float(pop_prev1)
    if type2 == "binary":
        sprev2, pprev2 = 0.5, float(pop_prev2)

    df = pl.DataFrame({
        "FILE": [zfile1, zfile2],
        "TYPE": [type1, type2],
        "SPREV": [sprev1, sprev2],
        "PPREV": [pprev1, pprev2],
        "NAME": [str(trait1), str(trait2)],
    })

    # tab-separated .txt (ldsc_preprocess reads it with sep='\t')
    df.write_csv(out_path, separator="\t", null_value="nan")
    return df

def ldsc_preprocess(
        input_list: str, # path to <out>.pleio_input.txt
        ldsc_ref: str, # eur_w_ld_chr folder (params.ldsc_ref)
        out: str, # output folder -> metain/sg/ce .txt.gz
    ):

    # from Docker image (wrapper -> python2.7 ldsc_preprocess.py)
    # LDSC needs the trailing "/" (it appends 1.l2.ldscore.gz etc.)
    ref = ldsc_ref.rstrip("/") + "/"

    cmd = [
        "pleio_preprocess",
        "--input", input_list,
        "--ref-ld-chr", ref,
        "--w-ld-chr", ref,
        "--out", out,
    ]
    subprocess.run(cmd, check=True)
    return f"{out}/metain.txt.gz", f"{out}/sg.txt.gz", f"{out}/ce.txt.gz"

def run_pleio(
        metain_file: str, # output from ldsc_preprocess
        sg_file: str,
        ce_file: str,
        out: str,
        blup: bool,
        n_cores: int,
        nis: int = 100000,
    ):

    """ Run PLEIO + --blup arg! -> pleio """

    blup_cmd = [
        "pleio",
        "--metain", metain_file,
        "--sg", sg_file,
        "--ce", ce_file,
        "--nis", str(nis),
        "--parallel",
        "--ncores", str(n_cores),
        "--create",
        "--blup",
        "--out", out,
    ]

    no_blup_cmd = [
        "pleio",
        "--metain", metain_file,
        "--sg", sg_file,
        "--ce", ce_file,
        "--nis", str(nis),
        "--parallel",
        "--ncores", str(n_cores),
        "--create",
        "--out", out,
    ]

    if blup:
        subprocess.run(blup_cmd, check=True)
        return f"{out}.txt.gz", f"{out}.blup.gz"
    else:
        subprocess.run(no_blup_cmd, check=True)
        return f"{out}.txt.gz", None


def combine_results(
        pleio_file: str, # <out>.txt.gz from run_pleio
        blup_file, # <out>.blup.gz or None
        positions: pl.DataFrame, # from make_z_files
        trait1,
        trait2,
        out: str,
        thresh: float, # pleio_p cut-off (default 5e-8)
        trait_p: float, # per-trait P cut-off for "shared" (default 0.05)
    ):

    """
    Final table (same layout as conjFDR):
    SNP CHR BP A1 A2 Z1 P1 Z2 P2 tausq pleio_stat pleio_p LS_stat LS_p
    [BLUP1 BLUP1_SE BLUP2 BLUP2_SE] direction shared novel
    -----
    shared: P1 < trait_p AND P2 < trait_p (Guo et al. 2022/2023)
    novel:  P1 >= thresh AND P2 >= thresh -> not found by either GWAS alone (Lee et al. 2021)
    """

    res = pl.read_csv(pleio_file, separator="\t")
    res = positions.join(res, on="SNP", how="inner")

    if blup_file is not None:
        blup = pl.read_csv(blup_file, separator="\t")
        blup = blup.rename({
            str(trait1): "BLUP1",
            f"{trait1}_se": "BLUP1_SE",
            str(trait2): "BLUP2",
            f"{trait2}_se": "BLUP2_SE",
        })
        res = res.join(blup, on="SNP", how="left")

    # same sign in both traits (on trait 1's A1)?
    res = res.with_columns(
        pl.when((pl.col("Z1") * pl.col("Z2")) > 0)
        .then(pl.lit("concordant"))
        .otherwise(pl.lit("discordant"))
        .alias("direction")
    )

    # nominal in BOTH traits -> truly shared, not driven by one trait
    res = res.with_columns(
        ((pl.col("P1") < trait_p) & (pl.col("P2") < trait_p)).alias("shared")
    )

    # not genome-wide in either single GWAS -> new from combining
    res = res.with_columns(
        ((pl.col("P1") >= thresh) & (pl.col("P2") >= thresh)).alias("novel")
    )

    res = res.sort("CHR", "BP")

    passed = res.filter(pl.col("pleio_p") < thresh)
    print(f"SNPs tested: {res.height}")
    print(f"SNPs with pleio_p < {thresh}: {passed.height}")
    print(f"  of which shared (both P < {trait_p}): {passed.filter(pl.col('shared')).height}")
    print(f"  of which novel (neither P < {thresh}): {passed.filter(pl.col('novel')).height}")

    res.write_csv(f"{out}.pleio.all.tsv", separator="\t")
    passed.write_csv(f"{out}.pleio.pass.tsv", separator="\t")

    return res


def main():
    parser = argparse.ArgumentParser(description="PLEIO for one trait pair (sherlock pleio step)")

    # trait 1
    parser.add_argument("--trait1",      required=True)
    parser.add_argument("--sumstats1",   required=True, help="QC file (.qc.tsv.gz)")
    parser.add_argument("--trait-type1", required=True, choices=["cc", "quant"])
    parser.add_argument("--pop-prev1",   type=float, default=None)

    # trait 2
    parser.add_argument("--trait2",      required=True)
    parser.add_argument("--sumstats2",   required=True, help="QC file (.qc.tsv.gz)")
    parser.add_argument("--trait-type2", required=True, choices=["cc", "quant"])
    parser.add_argument("--pop-prev2",   type=float, default=None)

    # refs + settings
    parser.add_argument("--ldsc-ref", required=True, help="eur_w_ld_chr folder")
    parser.add_argument("--nis",      type=int, default=100000, help="importance sampling draws")
    parser.add_argument("--blup",     action="store_true", help="also estimate BLUP effects per trait")
    parser.add_argument("--thresh",   type=float, default=5e-8, help="pleio_p cut-off for the pass file")
    parser.add_argument("--trait-p",  type=float, default=0.05, help="per-trait P cut-off for the 'shared' flag")
    parser.add_argument("--cores",    type=int, default=1)
    parser.add_argument("--out",      required=True)

    args = parser.parse_args()

    # 1. Z files (aligned to trait 1)
    zfile1, zfile2, positions = make_z_files(
        args.trait1, args.trait2,
        args.sumstats1, args.sumstats2,
        args.out,
    )

    # 2. input list
    input_list = f"{args.out}.pleio_input.txt"
    gen_input_file(
        args.trait1, args.trait2,
        zfile1, zfile2,
        args.trait_type1, args.trait_type2,
        input_list,
        args.pop_prev1, args.pop_prev2,
    )

    # 3. ldsc_preprocess -> metain / sg / ce
    metain, sg, ce = ldsc_preprocess(input_list, args.ldsc_ref, f"{args.out}_pleio_in")

    # 4. PLEIO (+ BLUP)
    pleio_file, blup_file = run_pleio(metain, sg, ce, f"{args.out}.pleio", args.blup, args.cores, args.nis)

    # 5. positions + BLUP + direction + shared/novel -> all / pass
    combine_results(
        pleio_file, blup_file, positions,
        args.trait1, args.trait2,
        args.out, args.thresh, args.trait_p,
    )


if __name__ == "__main__":
    main()