#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicSEM)
})

get_args <- function() {
  # takes munged sumstats as input!!!!
  option_list <- list(
    
    # trait 1
    make_option("--trait1",       type = "character", help = "Trait 1 name (meta.id)"),
    make_option("--sumstats1",    type = "character", help = "Trait 1 munged sumstats (<trait>.sumstats.gz)"),
    make_option("--trait-type1",  type = "character", help = "Trait 1 type: quant/cc"),
    make_option("--pop-prev1",    type = "double",    default = NA, help = "Trait 1 population prevalence (cc only)"),
    
    # trait 2
    make_option("--trait2",       type = "character", help = "Trait 2 name (meta.id)"),
    make_option("--sumstats2",    type = "character", help = "Trait 2 munged sumstats (<trait>.sumstats.gz)"),
    make_option("--trait-type2",  type = "character", help = "Trait 2 type: quant/cc"),
    make_option("--pop-prev2",    type = "double",    default = NA, help = "Trait 2 population prevalence (cc only)"),
    
    # reference + output
    make_option("--ld",           type = "character", help = "LD score folder (eur_w_ld_chr/)"),
    make_option("--wld",          type = "character", help = "LD weights folder (usually same as --ld)"),
    make_option("--out",          type = "character", help = "Output prefix (pair_id, e.g. AD_SCZ)")
  )
  
  args <- parse_args(OptionParser(option_list = option_list))
  
  # check if essential args == there!!
  for (a in c("trait1", "sumstats1", "trait-type1", "trait2", "sumstats2", "trait-type2", "ld", "wld", "out")) {
    if (is.null(args[[a]])) stop("--", a, " is required!")
  }
  return(args)
}

run_ldsc <- function(args) {
  
  get_sample_prev <- function(trait_type) {
    if (trait_type == "cc") {
      return(0.5)
    }
    return(NA)
  }
  
  # population-wide float val for quant traits
  get_pop_prev <- function(trait_type, pop_prev) {
    if (trait_type == "cc") {
      return(pop_prev)
    }
    return(NA)
  }
  
  sample_prev <- c(
    get_sample_prev(args$`trait-type1`),
    get_sample_prev(args$`trait-type2`)
  )
  
  pop_prev <- c(
    get_pop_prev(args$`trait-type1`, args$`pop-prev1`),
    get_pop_prev(args$`trait-type2`, args$`pop-prev2`)
  )
  
  result <- ldsc(
    traits          = c(args$sumstats1, args$sumstats2), 
    trait.names     = c(args$trait1, args$trait2),
    sample.prev     = sample_prev,
    population.prev = pop_prev,
    ld              = args$ld,
    wld             = args$wld,
    ldsc.log        = args$out,
    stand           = TRUE
  )
  
  return(result)
}

# assemble LDSC outputs from log files into an actual dataframe
write_results <- function(result, args) {
  
  # matrix outputs from genSEM LDSC
  S       <- result$S
  V       <- result$V
  S_stand <- result$S_Stand
  V_stand <- result$V_Stand
  I       <- result$I
  rg      <- S_stand[1, 2]
  rg_se   <- sqrt(V_stand[2, 2])
  rg_z    <- rg / rg_se
  rg_p    <- 2 * pnorm(-abs(rg_z))
  
  out <- data.frame(
    pair         = args$out,
    trait1       = args$trait1,
    trait2       = args$trait2,
    h2_1         = S[1, 1],
    h2_1_se      = sqrt(V[1, 1]),
    h2_2         = S[2, 2],
    h2_2_se      = sqrt(V[3, 3]),
    gcov         = S[1, 2],
    gcov_se      = sqrt(V[2, 2]),
    rg           = rg,
    rg_se        = rg_se,
    rg_z         = rg_z,
    rg_p         = rg_p,
    intercept_1  = I[1, 1],
    intercept_2  = I[2, 2],
    intercept_12 = I[1, 2],
    n_snps       = result$m
  )
  
  write.table(out, paste0(args$out, ".ldsc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  return(out)

}

main <- function() {
  args   <- get_args()
  result <- run_ldsc(args = args)
  write_results(result = result, args = args)
}

main()