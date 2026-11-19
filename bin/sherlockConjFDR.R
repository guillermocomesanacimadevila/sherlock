#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(cfdr.pleio)
    library(optparse)
    library(data.table)
})

get_args <- function() {

    # compilation of args to run conjFDR (pairwise!!!)
    # vignette("Introduction", package = "cfdr.pleio")

    option_list <- list(

        # trait 1
        make_option("--trait1",       type = "character", help = "Trait 1 name (meta.id)"),
        make_option("--sumstats1",    type = "character", help = "Trait 1 QC'd sumstats (<trait>.qc.tsv.gz)"),

        # trait 2
        make_option("--trait2",       type = "character", help = "Trait 2 name (meta.id)"),
        make_option("--sumstats2",    type = "character", help = "Trait 2 QC'd sumstats (<trait>.qc.tsv.gz)"),

        # reference
        make_option("--ref",          type = "character", help = "cfdr.pleio reference folder (all_chr_perVariant.rds + chrXX_LDpairs.rds)"),

        # settings
        make_option("--n-iter",       dest = "n_iter",       type = "integer", default = 500,    help = "Number of random pruning iterations"),
        make_option("--seed",         type = "integer",      default = 154226, help = "Random seed for pruning"),
        make_option("--maf-min",      dest = "maf_min",      type = "double",  default = 0.005,  help = "Minimum MAF (reference panel)"),
        make_option("--exclude-mhc",  dest = "exclude_mhc",  action = "store_true", default = FALSE, help = "Exclude MHC (chr6 25.1-33.9 Mb, GRCh37)"),
        make_option("--exclude-apoe", dest = "exclude_apoe", action = "store_true", default = FALSE, help = "Exclude APOE +-1 Mb (chr19, GRCh37)"),
        make_option("--conjfdr-thresh", dest = "conjfdr_thresh", type = "double", default = 0.05, help = "conjFDR threshold for the pass file"),

        # output
        make_option("--out",          type = "character", help = "Output prefix (pair_id, e.g. AD_SCZ)")
    )

    # parse args
    args <- parse_args(OptionParser(option_list = option_list))

    # ensure essential args == there
    for (a in c("trait1", "sumstats1", "trait2", "sumstats2", "ref", "out")) {
        if (is.null(args[[a]])) stop("--", a, " is required!")
    }

    return(args)
}

mhc_range <- function() {
    # pleiofdr default
    return(data.frame(chr = 6, from = 25119106, to = 33854733))
}

apoe_range <- function() {
    # APOE (45,409,011-45,412,650) +- 1 Mb
    return(data.frame(chr = 19, from = 44409011, to = 46412650))
}

get_exclusion <- function(args) {
    # regions to drop (NULL if none)
    exclusion <- NULL
    if (args$exclude_mhc) {
        exclusion <- rbind(exclusion, mhc_range())
    }
    if (args$exclude_apoe) {
        exclusion <- rbind(exclusion, apoe_range())
    }
    return(exclusion)
}

read_trait <- function(path) {
    # cfdr.pleio only accepts SNP / BETA / PVAL as names
    dt <- fread(path, select = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "P"))
    setnames(dt, "P", "PVAL")
    return(dt)
}

run_conjfdr <- function(t1, t2, args) {

    cf <- cfdr_pleio$new()

    # cfdr.pleio == only 3 cols
    trait1    <- t1[, c("SNP", "BETA", "PVAL")]
    trait2    <- t2[, c("SNP", "BETA", "PVAL")]
    ref       <- refdata_location(args$ref)
    exclusion <- get_exclusion(args)

    if (is.null(exclusion)) {
        cf$init_data(
            trait1            = trait1,
            trait2            = trait2,
            trait_names       = c(args$trait1, args$trait2),
            refdat            = ref,
            local_refdat_path = "local_ref",
            filter_maf_min    = args$maf_min
        )
    } else {
        cf$init_data(
            trait1            = trait1,
            trait2            = trait2,
            trait_names       = c(args$trait1, args$trait2),
            refdat            = ref,
            local_refdat_path = "local_ref",
            filter_maf_min    = args$maf_min,
            exclusion_range   = exclusion
        )
    }

    # random pruning
    cf$initialize_pruning_index(n_iter = args$n_iter, seed = args$seed)

    # cond FDR both ways -> conj_fdr = max(cfdr12, cfdr21)
    cf$calculate_cond_fdr(fdr_trait = 1)
    cf$calculate_cond_fdr(fdr_trait = 2)
    cat("GC correction factors:", cf$gc_correct$corrfac, "\n")
    return(cf$get_trait_results())
}

write_results <- function(res, t1, t2, args) {

    # positions + alleles from trait 1 (target build)
    pos1 <- t1[, c("SNP", "CHR", "BP", "A1", "A2")]

    # alleles from trait 2 (renamed so they don't clash)
    alleles2 <- t2[, c("SNP", "A1", "A2")]
    setnames(alleles2, c("A1", "A2"), c("A1_t2", "A2_t2"))

    res <- merge(res, pos1,     by = "SNP")
    res <- merge(res, alleles2, by = "SNP")

    # same alleles (either order) -> keep; else drop (can't align)
    same    <- res$A1_t2 == res$A1 & res$A2_t2 == res$A2
    swapped <- res$A1_t2 == res$A2 & res$A2_t2 == res$A1
    keep    <- same | swapped
    cat("SNPs dropped (alleles don't match between traits):", sum(!keep), "\n")
    res <- res[keep, ]

    # flip BETA2 onto trait 1's A1 -> A1 = effect allele for both betas
    flip <- res$A1_t2 != res$A1
    res$BETA2[flip] <- -res$BETA2[flip]

    # same sign = concordant, else discordant
    res$direction <- ifelse(sign(res$BETA1) == sign(res$BETA2), "concordant", "discordant")

    # drop helper cols
    res$A1_t2 <- NULL
    res$A2_t2 <- NULL

    # order cols + rows
    setcolorder(res, c("SNP", "CHR", "BP", "A1", "A2"))
    res <- res[order(res$CHR, res$BP), ]

    # all SNPs
    fwrite(res, paste0(args$out, ".conjfdr.all.tsv"), sep = "\t")

    # SNPs passing conjFDR threshold -> loci step
    pass <- res[res$conj_fdr < args$conjfdr_thresh, ]
    cat("SNPs with conj_fdr <", args$conjfdr_thresh, ":", nrow(pass), "\n")
    fwrite(pass, paste0(args$out, ".conjfdr.pass.tsv"), sep = "\t")
}

main <- function() {
    args <- get_args()
    t1   <- read_trait(args$sumstats1)
    t2   <- read_trait(args$sumstats2)
    res  <- run_conjfdr(t1, t2, args)
    write_results(res, t1, t2, args)
}

main()