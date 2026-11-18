#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(HDL)
    library(optparse)
    library(data.table)
    library(parallel)
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
        make_option("--ld",         type = "character", help = "HDL-L LD folder (LD.path, *_LDSVD.rda files)"),
        make_option("--bim",        type = "character", help = "HDL-L bim folder (bim.path)"),
        make_option("--n-ref",      type = "double",    default = 335272, help = "Reference panel sample size (HDL-L UKB default: 335272)"),

        # settings
        make_option("--n0",         type = "double",    default = 0,    help = "Number of overlapping samples between the two GWAS (0 = no overlap)"),
        make_option("--eigen-cut",  type = "double",    default = 0.99, help = "Eigenvalue cut-off"),
        make_option("--univ-thresh", type = "double",   default = 0.05, help = "Local h2 p-value threshold: rg p is only corrected where both traits pass"),

        # whole genome by default; --chr for testing one chromosome
        make_option("--chr",        type = "integer",   default = NULL, help = "Chromosome to run (empty = whole genome)"),
        make_option("--cores",      type = "integer",   default = 1,    help = "Number of cores (task.cpus)"),

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
}

# loci = LD files in the panel (only these can run); coordinates from the matching bim
get_loci <- function(args) {
    ld_files <- grep("_LDSVD\\.rda$", list.files(args$ld), value = TRUE)
    loci <- data.frame(
        chr   = as.integer(sub("^ukb_(chr)?([0-9]+)\\.([0-9]+)_.*", "\\2", ld_files)),
        piece = as.integer(sub("^ukb_(chr)?([0-9]+)\\.([0-9]+)_.*", "\\3", ld_files))
    )
    if (!is.null(args$chr)) {
        loci <- loci[loci$chr == args$chr, ]
    }
    if (nrow(loci) == 0) stop("No loci found for chr ", args$chr)
    loci <- loci[order(loci$chr, loci$piece), ]
    snps <- vector("list", nrow(loci))

    # first/last SNP position (GRCh37) + number of reference SNPs + SNP ids
    for (i in seq_len(nrow(loci))) {
        bim_file <- list.files(args$bim, pattern = paste0("^ukb_(chr)?", loci$chr[i], "\\.", loci$piece[i], "_.*\\.bim$"), full.names = TRUE)
        bim <- fread(bim_file[1], select = c(2, 4))
        loci$start_grch37[i] <- min(bim[[2]])
        loci$stop_grch37[i]  <- max(bim[[2]])
        loci$n_snps_ref[i]   <- nrow(bim)
        snps[[i]]            <- bim[[1]]
    }
    return(list(loci = loci, snps = unique(unlist(snps))))
}

run_hdl_l <- function(args) {

    panel <- get_loci(args)
    loci  <- panel$loci

    # keep only panel SNPs once (HDL.L drops the rest anyway)
    gwas1 <- read_sumstats(args$sumstats1)
    gwas2 <- read_sumstats(args$sumstats2)
    gwas1 <- gwas1[gwas1$SNP %in% panel$snps, ]
    gwas2 <- gwas2[gwas2$SNP %in% panel$snps, ]
    cat("Panel SNPs in GWAS 1:", nrow(gwas1), " GWAS 2:", nrow(gwas2), "\n")
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
    # a killed worker (e.g. out of memory) silently returns NULL for all its loci;
    # R only warns "did not deliver results" -> turn that warning into an error
    results <- withCallingHandlers(
        mclapply(seq_len(nrow(loci)), run_locus, mc.cores = args$cores),
        warning = function(w) {
            if (grepl("did not deliver", conditionMessage(w))) {
                stop("a parallel worker died (likely out of memory): increase memory or lower --cores", call. = FALSE)
            }
        }
    )
    results <- Filter(is.data.frame, results)    # drop failed loci
    result  <- do.call(rbind, results)
    if (is.null(result)) stop("HDL-L failed for every locus!")
    result <- merge(loci, result, by = c("chr", "piece"))
    return(result[order(result$chr, result$piece), ])
}

write_results <- function(result, args) {

    out <- data.frame(
        pair         = args$out,
        trait1       = args$trait1,
        trait2       = args$trait2,
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

    # per-pair multiple testing (whole genome only): both local h2 > 0 AND p < univ-thresh, rg estimated
    # (h2 > 0 needed explicitly: HDL-L gives h2 p = 0 when h2 = 0)
    tested <- !is.na(out$rg_p) & out$h2_1 > 0 & out$h2_2 > 0 &
              out$h2_1_p < args$`univ-thresh` & out$h2_2_p < args$`univ-thresh`
    tested[is.na(tested)] <- FALSE
    out$tested    <- tested
    out$rg_p_fdr  <- NA_real_
    out$rg_p_bonf <- NA_real_
    if (is.null(args$chr) && any(tested)) {
        out$rg_p_fdr[tested]  <- p.adjust(out$rg_p[tested], method = "BH")
        out$rg_p_bonf[tested] <- p.adjust(out$rg_p[tested], method = "bonferroni")
    }

    write.table(out, paste0(args$out, ".hdll.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    return(out)
}

main <- function() {
    args <- get_args()
    result <- run_hdl_l(args = args)
    write_results(
        result = result,
        args = args
    )
}

main()