#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(HDL)
    library(optparse)
    library(data.table)
})

# allocate HDL-L module to -> 'process_high'
get_args <- function() {
    
    # args for HDL-L    
    option_list <- list(

        # trait1
        make_option("--trait1",     type = "character", help = "Trait 1 name (meta.id)"),
        make_option("--sumstats1",  type = "character", help = "Trait 1 QC'd sumstats (<trait>.qc.tsv.gz)"),

        # trait 2
        make_option("--trait2",     type = "character", help = "Trait 2 name (meta.id)"),
        make_option("--sumstats2",  type = "character", help = "Trait 2 QC'd sumstats (<trait>.qc.tsv.gz)"),

        # reference stuff
        make_option("--ld",         type = "character", help = "HDL-L LD folder (with HDLL_LOC_snps.RData)"),
        make_option("--bim",        type = "character", help = "HDL-L bim folder"),
        make_option("--n-ref",      type = "double",    default = NULL, help = "Reference panel sample size (empty = HDL.L default)"),

        # settings
        # idk abt this one
        make_option("--n0",         type = "double",    default = 0,    help = "Number of overlapping samples between the two GWAS"),
        make_option("--eigen-cut",  type = "double",    default = 0.99, help = "Eigenvalue cut-off"),

        # Nextflow runs one chromosome per task
        make_option("--chr",        type = "integer",   default = NULL, help = "Chromosome to run (empty = whole genome)"),

        # output
        make_option("--out",        type = "character", help = "Output prefix (pair_id, e.g. AD_SCZ)")
    )

    args <- parse_args(OptionParser(option_list = option_list))
    for (a in c("trait1", "sumstats1", "trait2", "sumstats2", "ld", "bim", "out")) {
        if (is.null(args[[a]])) stop("--", a, " is required!")
    }
    return(args)
}

read_sumstats <- function(path) {
  df <- fread(path, select = c("SNP", "A1", "A2", "N", "BETA", "SE"))
  df[, Z := BETA / SE]
  return(as.data.frame(df[, .(SNP, A1, A2, N, Z)]))

# list of loci (chr + piece) from the HDL-L reference
get_loci <- function(args) {
    load(file.path(args$ld, "HDLL_LOC_snps.RData")) # NEWLOC
    loci <- data.frame(
        locus = as.integer(NEWLOC$LOC),
        chr = as.integer(NEWLOC$CHR),
        piece = as.integer(NEWLOC$piece),
        start_grch37 = NEWLOC$START,
        stop_grch37 = NEWLOC$STOP,
        n_snps_ref = unlist(nsnps.list.imputed)
    )

    if (!is.null(args$chr)) {
        loci <- loci[loci$chr == args$chr, ]
    }

    if (nrow(loci) == 0) stop("No loci found for chr ", args$chr)
    return(loci)
}

run_hdl_l <- function(args) {

    gwas1 <- read_sumstats(args$sumstats1)
    gwas2 <- read_sumstats(args$sumstats2)
    loci  <- get_loci(args)
    cat("Running HDL-L on", nrow(loci), "loci\n")

    # one locus -> one row (failed loci are skipped, not fatal)
    run_locus <- function(i) {
        tryCatch(
            HDL.L(
                gwas1.df   = gwas1,
                gwas2.df   = gwas2,
                Trait1name = args$trait1,
                Trait2name = args$trait2,
                LD.path    = args$ld,
                bim.path   = args$bim,
                Nref       = args$`n-ref`,
                N0         = args$n0,
                chr        = loci$chr[i],
                piece      = loci$piece[i],
                eigen.cut  = args$`eigen-cut`
            ),
            error = function(e) {
                cat("Locus chr", loci$chr[i], "piece", loci$piece[i], "failed:", conditionMessage(e), "\n")
                return(NULL)
            }
        )
    }
    results <- mclapply(seq_len(nrow(loci)), run_locus, mc.cores = args$cores)
    result  <- do.call(rbind, results)
    if (is.null(result)) stop("HDL-L failed for every locus!")
    result <- merge(loci, result, by = c("chr", "piece"))
    return(result[order(result$locus), ])
}

write_results <- function(result, args) {

    out <- data.frame(
        pair         = args$out,
        trait1       = args$trait1,
        trait2       = args$trait2,
        locus        = result$locus,
        chr          = result$chr,
        piece        = result$piece,
        start_grch37 = result$start_grch37,
        stop_grch37  = result$stop_grch37,
        n_snps_ref   = result$n_snps_ref,
        eigen_use    = result$eigen_use,
        h2_1         = result$Heritability_1,
        h2_1_p       = result$P_value_Heritability_1,
        h2_2         = result$Heritability_2,
        h2_2_p       = result$P_value_Heritability_2,
        gcov         = result$Genetic_Covariance,
        rg           = result$Genetic_Correlation,
        rg_lower     = result$Lower_bound_rg,
        rg_upper     = result$Upper_bound_rg,
        rg_p         = result$P
    )

    out$h2_1_p_fdr <- p.adjust(out$h2_1_p, method = "fdr")
    out$h2_2_p_fdr <- p.adjust(out$h2_2_p, method = "fdr")
    out$rg_p_fdr   <- p.adjust(out$rg_p, method = "fdr")
    write.table(out, paste0(args$out, ".hdll.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    return(out)
}

main <- function() {
    args <- get_args()
    result <- run_hdl(args = args)
    write_results(
        result = result,
        args = args
    )
}

main()