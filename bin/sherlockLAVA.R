#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(LAVA)
  library(optparse)
  library(parallel)
})

get_args <- function() {

    option_list <- list(

        # trait 1 stuff
        make_option("--trait1",         type = "character", help = "Trait 1 name (meta.id)"),
        make_option("--sumstats1",      type = "character", help = "Trait 1 QC'd sumstats (<trait>.qc.tsv.gz)"),
        make_option("--trait-type1",    dest = "trait_type1",    type = "character", help = "Trait 1 type (cc / quant)"),
        make_option("--pop-prev1",      dest = "pop_prev1",      type = "double",    default = NA, help = "Trait 1 population prevalence (cc only, optional)"),

        # trait 2 stuff
        make_option("--trait2",         type = "character", help = "Trait 2 name (meta.id)"),
        make_option("--sumstats2",      type = "character", help = "Trait 2 QC'd sumstats (<trait>.qc.tsv.gz)"),
        make_option("--trait-type2",    dest = "trait_type2",    type = "character", help = "Trait 2 type (cc / quant)"),
        make_option("--pop-prev2",      dest = "pop_prev2",      type = "double",    default = NA, help = "Trait 2 population prevalence (cc only, optional)"),

        # reference stuff
        make_option("--ref",            type = "character", help = "LAVA LD reference prefix (without _chr part)"),
        make_option("--loci",           type = "character", help = "Locus definition file (LOC, CHR, START, STOP [, SNPS])"),
        make_option("--sample-overlap", dest = "sample_overlap", type = "character", default = NULL, help = "Sample overlap file (empty = no overlap)"),

        # settings
        make_option("--univ-thresh",    dest = "univ_thresh",    type = "double",    default = 0.05, help = "Local h2 p-value threshold for running the bivariate test"),

        # whole genome by default; --chr for testing one chromosome
        make_option("--chr",            type = "integer",   default = NULL, help = "Chromosome to run (empty = whole genome)"),
        make_option("--cores",          type = "integer",   default = 1,    help = "Number of cores (task.cpus)"),

        # output
        make_option("--out",            type = "character", help = "Output prefix (pair_id, e.g. AD_SCZ)")
    )

    args <- parse_args(OptionParser(option_list = option_list))

    # ensure essential args == there
    for (a in c("trait1", "sumstats1", "trait_type1", "trait2", "sumstats2", "trait_type2", "ref", "loci", "out")) {
        if (is.null(args[[a]])) stop("--", a, " is required!")
    }

    # check for cc and quant status == actually only either of those 2 strings
    for (t in c(args$trait_type1, args$trait_type2)) {
        if (!t %in% c("quant", "cc")) stop("--trait-type must be quant or cc, got: ", t)
    }
    
    return(args)
}

make_info <- function(args) {
    cases <- c(NA, NA)
    if (args$trait_type1 == "cc") cases[1] <- 1
    if (args$trait_type2 == "cc") cases[2] <- 1

    # craft info df 
    info <- data.frame(
        phenotype  = c(args$trait1, args$trait2),
        cases      = cases,
        controls   = cases,
        prevalence = c(args$pop_prev1, args$pop_prev2),
        filename   = c(args$sumstats1, args$sumstats2)
    )

    info_file <- tempfile(fileext = ".txt")
    write.table(info, info_file, row.names = FALSE, quote = FALSE, sep = "\t")
    return(info_file)
}

# != sample_overlap_file <- function()
# We assume sample.overlap == NULL (maybe in future versions we account for it!!)

run_lava <- function(info_file, args) {

    # lava input
    input <- process.input(
        input.info.file     = info_file,
        sample.overlap.file = args$sample_overlap,
        ref.prefix          = args$ref,
        phenos              = c(args$trait1, args$trait2)
    )

    # now loci (one chromosome if --chr is used)
    loci <- read.loci(loc.file = args$loci)
    if (!is.null(args$chr)) loci <- loci[loci$CHR == args$chr, ]
    cat("Running LAVA on", nrow(loci), "loci with", args$cores, "cores\n")

    # one locus: univ (local h2) + bivar (local rg, only if both traits pass univ.thresh)
    run_locus <- function(i) {
        locus <- process.locus(
            locus = loci[i, ],
            input = input
        )

        if (is.null(locus)) return(NULL)

        loc_info <- data.frame(
            locus  = locus$id,
            chr    = locus$chr,
            start  = locus$start,
            stop   = locus$stop,
            n_snps = locus$n.snps,
            n_pcs  = locus$K
        )

        out <- run.univ.bivar(
            locus       = locus,
            univ.thresh = args$univ_thresh
        )

        univ  <- cbind(loc_info, out$univ)
        bivar <- NULL
        if (!is.null(out$bivar)) bivar <- cbind(loc_info, out$bivar)
        return(list(univ = univ, bivar = bivar))
    }

    # failed loci == skipped
    res <- mclapply(seq_len(nrow(loci)), function(i) {
        tryCatch(run_locus(i), error = function(e) {
            message("Locus ", loci$LOC[i], " failed: ", conditionMessage(e))
            return(NULL)
        })
    }, mc.cores = args$cores)
    cat("Loci skipped (too few SNPs or failed):", sum(sapply(res, is.null)), "\n")

    univ  <- do.call(rbind, lapply(res, function(r) r$univ))
    bivar <- do.call(rbind, lapply(res, function(r) r$bivar))
    return(list(univ = univ, bivar = bivar))
}

main <- function() {
    args <- get_args()

    # LAVA p-values/CIs are simulated -> fixed seed (one stream per core) so -resume gives same results
    RNGkind("L'Ecuyer-CMRG")
    set.seed(1)

    info_file <- make_info(args = args)
    out       <- run_lava(args = args, info_file = info_file)

    # always write both files (bivar is NULL if no locus passes univ.thresh in both traits)
    univ  <- out$univ
    bivar <- out$bivar
    if (is.null(univ))  univ  <- data.frame()
    if (is.null(bivar)) bivar <- data.frame()

    # per-pair multiple testing (whole genome only): loci with a bivar p
    # (both traits passed univ.thresh, rho estimated); h2 = 0 loci never reach bivar
    if (nrow(bivar) > 0) {
        bivar$rg_p_fdr  <- NA_real_
        bivar$rg_p_bonf <- NA_real_
        tested <- !is.na(bivar$p)
        if (is.null(args$chr) && any(tested)) {
            bivar$rg_p_fdr[tested]  <- p.adjust(bivar$p[tested], method = "BH")
            bivar$rg_p_bonf[tested] <- p.adjust(bivar$p[tested], method = "bonferroni")
        }
    }

    write.table(univ,  paste0(args$out, ".lava_univ.tsv"),  row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(bivar, paste0(args$out, ".lava_bivar.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
    cat("Done:", nrow(univ), "univ rows,", nrow(bivar), "bivar rows\n")
}

main()