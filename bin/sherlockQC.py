import argparse

import polars as pl
from sherlock_tools.preprocessing import QualityControl


def main():

    """
    Kinda self-explanatory script but..
    - Inherits QC class from sherlock-tools
    - Reads args from YAML file in params/ or inputs/
    - Reads X-separated sumstats (1 by 1) -> open as pl.DataFrame
    - Parse inputs as QC args into run()
    - Save as .parquet in out/
    """

    parser = argparse.ArgumentParser(description="sherlock GWAS/xQTL QC")
    parser.add_argument("--sumstats", required=True, help="input summary statistics (.tsv / .tsv.gz)")
    parser.add_argument("--out", required=True, help="output QC'd sumstats (.parquet)")
    parser.add_argument("--sep", default="\t", help="input separator")
    parser.add_argument("--snp-col", required=True)
    parser.add_argument("--chr-col", required=True)
    parser.add_argument("--pos-col", required=True)
    parser.add_argument("--a1-col", required=True, help="effect allele")
    parser.add_argument("--a2-col", required=True, help="other allele")
    parser.add_argument("--beta-col", required=True)
    parser.add_argument("--se-col", required=True)
    parser.add_argument("--pval-col", required=True)
    parser.add_argument("--freq-col", default=None)
    parser.add_argument("--info-col", default=None)
    parser.add_argument("--n-col", default=None)
    parser.add_argument("--trait-type", required=True, choices=["quant", "cc"])
    parser.add_argument("--n", type=int, default=None)
    parser.add_argument("--n-cases", type=int, default=None)
    parser.add_argument("--n-controls", type=int, default=None)
    parser.add_argument("--maf-thresh", type=float, default=0.01)
    parser.add_argument("--info-thresh", type=float, default=None)
    parser.add_argument("--genome-build", required=True, choices=["GRCh37", "GRCh38"])
    parser.add_argument("--target-genome-build", required=True, choices=["GRCh37", "GRCh38"])
    parser.add_argument("--chain-file", default=None)
    parser.add_argument("--drop-mhc", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--drop-apoe", action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()

    # read everything as strings -> _fix_data_types does the casting
    sumstats = pl.read_csv(
        args.sumstats,
        separator=args.sep,
        infer_schema=False,
        null_values=["", "NA", "NaN", "nan", "."],
    )

    qc = QualityControl(
        sumstats=sumstats,
        snp_col=args.snp_col,
        chr_col=args.chr_col,
        pos_col=args.pos_col,
        a1_col=args.a1_col,
        a2_col=args.a2_col,
        beta_col=args.beta_col,
        se_col=args.se_col,
        pval_col=args.pval_col,
        maf_thresh=args.maf_thresh,
        genome_build=args.genome_build,
        target_genome_build=args.target_genome_build,
        trait_type=args.trait_type,
        n_col=args.n_col,
        n=args.n,
        n_cases=args.n_cases,
        n_controls=args.n_controls,
        info_col=args.info_col,
        info_thresh=args.info_thresh,
        freq_col=args.freq_col,
        chain_file=args.chain_file,
        drop_mhc=args.drop_mhc,
        drop_apoe=args.drop_apoe,
    )

    df = qc.run()
    df.write_parquet(args.out)
    print(f"Written: {args.out}")

if __name__ == "__main__":
    main()