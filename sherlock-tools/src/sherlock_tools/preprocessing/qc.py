import polars as pl
from liftover import ChainFile


class QualityControl:

    """ Object compiling all functions pertaining to GWAS/xQTL QC """

    def __init__(
            self,
            sumstats: pl.DataFrame,
            snp_col: str,
            chr_col: str,
            pos_col: str,
            a1_col: str,
            a2_col: str,
            beta_col: str,
            se_col: str,
            pval_col: str,
            maf_thresh: float,
            genome_build: str,
            target_genome_build: str,
            trait_type: str,
            n_col: str | None = None,
            n: int | None = None,
            n_cases: int | None = None,
            n_controls: int | None = None,
            info_col: str | None = None,
            info_thresh: float | None = None,
            freq_col: str | None = None,
            chain_file: str | None = None, # although in params.yaml arg -> liftover_dir/
            drop_mhc: bool = True,
            drop_apoe: bool = False,
        ):

        # data
        self.sumstats = sumstats

        # column mapping
        self.snp_col = snp_col
        self.chr_col = chr_col
        self.pos_col = pos_col
        self.a1_col = a1_col
        self.a2_col = a2_col
        self.beta_col = beta_col
        self.se_col = se_col
        self.pval_col = pval_col
        self.info_col = info_col
        self.freq_col = freq_col

        # trait + sample size
        self.trait_type = trait_type
        self.n_col = n_col
        self.n = n
        self.n_cases = n_cases
        self.n_controls = n_controls

        # thresholds
        self.maf_thresh = maf_thresh
        self.info_thresh = info_thresh

        # builds
        self.genome_build = genome_build
        self.target_genome_build = target_genome_build
        self.chain_file = chain_file

        # region removal
        self.drop_mhc = drop_mhc
        self.drop_apoe = drop_apoe

    def _fix_data_types(self):
        # strict=False -> values that can't be cast become null (caught by _correct_for_missingness)
        df = self.sumstats.with_columns(
            pl.col(self.snp_col).cast(pl.String),
            pl.col(self.chr_col).cast(pl.String),
            pl.col(self.pos_col).cast(pl.Int64, strict=False),
            pl.col(self.a1_col.upper()).cast(pl.String),
            pl.col(self.a2_col.upper()).cast(pl.String),
            pl.col(self.beta_col).cast(pl.Float64, strict=False),
            pl.col(self.se_col).cast(pl.Float64, strict=False),
            pl.col(self.pval_col).cast(pl.Float64, strict=False),
        )

        # optional cols if given
        optional_cols = [c for c in [self.n_col, self.freq_col, self.info_col] if c is not None]
        df = df.with_columns(pl.col(c).cast(pl.Float64, strict=False) for c in optional_cols)
        self.sumstats = df
        return df

    def _correct_for_missingness(self):
        n_before = self.sumstats.height
        df = self.sumstats.filter(~pl.all_horizontal(pl.all().is_null()))
        n_after = df.height
        print(f"Removed empty rows: {n_before - n_after}")
        print(f"Rows after empty-row removal: {n_after}")

        # remove rows missing any key column (+ optional cols if given)
        key_cols = [
            self.snp_col, self.chr_col, self.pos_col,
            self.a1_col, self.a2_col,
            self.beta_col, self.se_col, self.pval_col,
        ]

        key_cols += [c for c in [self.n_col, self.freq_col, self.info_col] if c is not None]
        n_before = df.height
        df = df.drop_nulls(subset=key_cols)
        n_after = df.height
        print(f"Removed rows with missing key values: {n_before - n_after}")
        print(f"Rows after missing-value removal: {n_after}")
        self.sumstats = df
        return df

    def _remove_duplicates(self):
        df = self.sumstats.unique(subset=[self.snp_col])
        return df

    def _remove_multiallelic(self):
        # run after duplicate removal -> any repeated CHR:POS left = different alleles
        n_before = self.sumstats.height
        df = self.sumstats.filter(
            pl.len().over([self.chr_col, self.pos_col]) == 1
        )
        n_after = df.height
        print(f"Removed multiallelic variants: {n_before - n_after}")
        print(f"Rows after multiallelic removal: {n_after}")
        self.sumstats = df
        return df

    def _remove_indels(self):
        bases = ["A", "T", "C", "G"]
        a1 = pl.col(self.a1_col)
        a2 = pl.col(self.a2_col)
        ok_len = (a1.str.len_chars() == 1) & (a2.str.len_chars() == 1)
        ok_bases = a1.is_in(bases) & a2.is_in(bases)
        n_before = self.sumstats.height
        df = self.sumstats.filter(ok_len & ok_bases)
        n_after = df.height
        print(f"Removed INDELs / non-ACGT: {n_before - n_after}")
        print(f"Rows after INDEL removal: {n_after}")
        self.sumstats = df
        return df

    def _remove_palindromes(self):
        # A/T and C/G SNPs -> strand can't be resolved, always dropped
        a1 = pl.col(self.a1_col)
        a2 = pl.col(self.a2_col)
        is_palindrome = (
            ((a1 == "A") & (a2 == "T")) |
            ((a1 == "T") & (a2 == "A")) |
            ((a1 == "C") & (a2 == "G")) |
            ((a1 == "G") & (a2 == "C"))
        )
        n_before = self.sumstats.height
        df = self.sumstats.filter(~is_palindrome)
        n_after = df.height
        print(f"Removed palindromic SNPs: {n_before - n_after}")
        print(f"Rows after palindrome removal: {n_after}")
        self.sumstats = df
        return df
    
    def _maf_filter(self):
        if self.freq_col is None:
            print("Skipping MAF filter: no frequency column supplied.")
            return self.sumstats

        df = self.sumstats.filter(
            pl.col(self.freq_col) >= self.maf_thresh
        )
        self.sumstats = df
        return df

    def _info_filter(self):
        if self.info_col is None:
            print("Skipping INFO filter: no frequency column supplied.")
            return self.sumstats

        df = self.sumstats.filter(
            pl.col(self.info_col <= self.info_thresh)
        )
        self.sumstats = df
        return df

    def _impute_effective_n(self):
        # per-SNP N given -> keep it (_standardise_shape renames n_col -> N)
        if self.n_col is not None:
            print(f"Using per-SNP sample size from column: {self.n_col}")
            return self.sumstats

        # depends on whether type == cc / quant
        if self.trait_type == "quant":
            if self.n is None:
                raise ValueError("trait_type 'quant' needs n_col or n")
            n_eff = int(round(self.n))
        elif self.trait_type == "cc":
            if self.n_cases is None or self.n_controls is None:
                raise ValueError("trait_type 'cc' needs n_col or n_cases + n_controls")
            n_eff = int(round(4 * self.n_cases * self.n_controls / (self.n_cases + self.n_controls)))
        else:
            raise ValueError(f"Unsupported trait_type: {self.trait_type}")

        df = self.sumstats.with_columns(pl.lit(n_eff).alias("N"))
        print(f"N set to: {n_eff} ({self.trait_type})")
        self.sumstats = df
        return df
    

    def _liftover():
        return

    # --- --- --- --- --- --- --- --- --- --- --- 
    # --- These bottom ones -> after _liftover()
    # --- --- --- --- --- --- --- --- --- --- --- 

    def _align_betas(self):
        # align to risk allele: BETA < 0 -> flip sign + swap A1/A2 (+ FRQ -> 1 - FRQ)
        flip = pl.col(self.beta_col) < 0
        n_flip = self.sumstats.filter(flip).height

        cols = [
            pl.when(flip).then(-pl.col(self.beta_col)).otherwise(pl.col(self.beta_col)).alias(self.beta_col),
            pl.when(flip).then(pl.col(self.a2_col)).otherwise(pl.col(self.a1_col)).alias(self.a1_col),
            pl.when(flip).then(pl.col(self.a1_col)).otherwise(pl.col(self.a2_col)).alias(self.a2_col),
        ]
        if self.freq_col is not None:
            cols.append(
                pl.when(flip).then(1 - pl.col(self.freq_col)).otherwise(pl.col(self.freq_col)).alias(self.freq_col)
            )

        df = self.sumstats.with_columns(cols)
        print(f"Flipped to risk allele: {n_flip}")
        self.sumstats = df
        return df

    # removing specific regions
    def _remove_mhc(self):

        # to run after liftover
        if not self.drop_mhc:
            print("Keeping MHC region!")
            return self.sumstats

        if self.target_genome_build == "GRCh38":
            mhc_start = 24999772
            mhc_end = 34032223
        elif self.target_genome_build == "GRCh37":
            mhc_start = 25000000
            mhc_end = 34000000
        else:
            raise ValueError(f"Unsupported target_genome_build: {self.target_genome_build}")

        n_before = self.sumstats.height
        df = self.sumstats.filter(
            ~(
                (pl.col(self.chr_col) == "6") &
                (pl.col(self.pos_col) >= mhc_start) &
                (pl.col(self.pos_col) <= mhc_end)
            )
        )
        n_after = df.height
        print(f"Removed MHC variants: {n_before - n_after} ({self.target_genome_build})!")
        self.sumstats = df
        return df

    def _remove_apoe(self):

        if not self.drop_apoe:
            print("Keeping APOE region!")
            return self.sumstats

        if self.target_genome_build == "GRCh38":
            apoe_start = 43495848
            apoe_end = 45996742
        elif self.target_genome_build == "GRCh37":
            apoe_start = 44000000
            apoe_end = 46500000
        else:
            raise ValueError(f"Unsupported target_genome_build: {self.target_genome_build}")
        
        n_before = self.sumstats.height
        df = self.sumstats.filter(
            ~(
                (pl.col(self.chr_col) == "19") &
                (pl.col(self.pos_col) >= apoe_start) &
                (pl.col(self.pos_col) <= apoe_end)
            )
        )
        n_after = df.height
        print(f"Removed APOE variants: {n_before - n_after} ({self.target_genome_build})!")
        self.sumstats = df
        return df

    def _standardise_shape(self):
        # rename cols to the sherlock format (A1 = effect allele)
        cols = [
            pl.col(self.snp_col).alias("SNP"),
            pl.col(self.chr_col).alias("CHR"),
            pl.col(self.pos_col).alias("BP"),
            pl.col(self.a1_col).alias("A1"),
            pl.col(self.a2_col).alias("A2"),
            pl.col(self.beta_col).alias("BETA"),
            pl.col(self.se_col).alias("SE"),
            pl.col(self.pval_col).alias("P"),
        ]

        # sample size: always one N column
        if self.n_col is not None:
            cols.append(pl.col(self.n_col).alias("N"))
        elif "N" in self.sumstats.columns:
            cols.append(pl.col("N"))

        # optional cols if given
        if self.freq_col is not None:
            cols.append(pl.col(self.freq_col).alias("FRQ"))
        if self.info_col is not None:
            cols.append(pl.col(self.info_col).alias("INFO"))

        df = self.sumstats.select(cols)
        self.sumstats = df
        return df