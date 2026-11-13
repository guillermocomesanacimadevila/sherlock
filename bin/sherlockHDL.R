#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(HDL)
  library(optparse)
  library(data.table)
})

get_args <- function() {
  option_list <- list(
    
    # trait 1
    make_option("--trait1",     type = "character", help = "Trait 1 name (meta.id)"),
    make_option("--sumstats1",  type = "character", help = "Trait 1 QC'd sumstats (<trait>.qc.tsv.gz)"),
    
    # trait 2
    make_option("--trait2",     type = "character", help = "Trait 2 name (meta.id)"),
    make_option("--sumstats2",  type = "character", help = "Trait 2 QC'd sumstats (<trait>.qc.tsv.gz)"),
    
    # reference + output
    make_option("--ref",        type = "character", help = "HDL reference panel folder (e.g. UKB_imputed_SVD_eigen99_extraction/)"),
    make_option("--n-ref",      type = "double",    default = 335265, help = "Sample size of the reference panel -> the one provided (UKB default: 335265)"),
    make_option("--out",        type = "character", help = "Output prefix (pair_id, e.g. AD_SCZ)")
  )
  
  args <- parse_args(OptionParser(option_list = option_list))
  
  # check if essential args == there!!
  for (a in c("trait1", "sumstats1", "trait2", "sumstats2", "ref", "out")) {
    if (is.null(args[[a]])) stop("--", a, " is required!")
  }
  return(args)
}

# from qc output
read_sumstats <- function(path) {
  df <- fread(path, select = c("SNP", "A1", "A2", "N", "BETA", "SE"))
  df[, Z := BETA / SE]
  return(as.data.frame(df[, .(SNP, A1, A2, N, Z)]))
}

run_hdl <- function(args) {
  
  gwas1 <- read_sumstats(args$sumstats1)
  gwas2 <- read_sumstats(args$sumstats2)

  result <- HDL.rg(
    gwas1.df    = gwas1,
    gwas2.df    = gwas2,
    LD.path     = args$ref,
    Nref        = args$`n-ref`,
    output.file = paste0(args$out, "_hdl.log")
  )
  return(result)
}

write_results <- function(result, args) {
  
  est <- result$estimates.df
  # rows: Heritability_1, Heritability_2, Genetic_Covariance, Genetic_Correlation
  
  rg    <- result$rg
  rg_se <- result$rg.se
  out <- data.frame(
    pair      = args$out,
    trait1    = args$trait1,
    trait2    = args$trait2,
    h2_1      = est["Heritability_1", "Estimate"],
    h2_1_se   = est["Heritability_1", "se"],
    h2_2      = est["Heritability_2", "Estimate"],
    h2_2_se   = est["Heritability_2", "se"],
    gcov      = est["Genetic_Covariance", "Estimate"],
    gcov_se   = est["Genetic_Covariance", "se"],
    rg        = rg,
    rg_se     = rg_se,
    rg_z      = rg / rg_se,
    rg_p      = result$P,
    eigen_use = result$eigen.use
  )
  
  write.table(out, paste0(args$out, ".hdl.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  return(out)
  
}

  
main <- function() {
  args   <- get_args()
  result <- run_hdl(args = args)
  write_results(result = result, args = args)
}

main()