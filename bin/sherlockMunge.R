#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicSEM)
})

# argument parser function!!!!
get_args <- function() {
  option_list <- list(
    make_option("--sumstats", type = "character",  help = "QC'd sumstats (<trait>.qc.tsv.gz)"),
    make_option("--trait", type = "character", help = "Trait name (meta.id)"),
    make_option("--hm3", type = "character", help = "HapMap3 SNP list (w_hm3.snplist)")
  )
  args <- parse_args(OptionParser(option_list = option_list))

  for (a in c("sumstats", "trait", "hm3")) {
    if (is.null(args[[a]])) stop("--", a, " is required")
  }
  return(args)
}

munge_sumstats <- function(trait, gwas_sumstats, hm3) {
  munge(
    files       = gwas_sumstats,
    trait.names = trait,
    hm3         = hm3,
    N           = NULL, # use N column from QC
    maf.filter  = 0,
    info.filter = 0,
    log.name    = trait
  )
}

main <- function() {
  args <- get_args()
  munge_sumstats(
    trait         = args$trait,
    gwas_sumstats = args$sumstats,
    hm3           = args$hm3
  )
}

main()