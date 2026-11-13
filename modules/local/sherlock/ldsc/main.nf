process SHERLOCK_LDSC {
    tag "${meta1.id}_${meta2.id}"
    label 'process_low'

    // GenomicSEM is not on conda -> own image only
    container "ghcr.io/guillermocomesanacimadevila/sherlock-genomicsem:0.0.5-da95d43"

    input:
    tuple val(meta1), path(sumstats1), val(meta2), path(sumstats2) // one pair: both traits' munged sumstats (sorted by trait name)
    path ld_ref                                                     // LD score folder (eur_w_ld_chr/)

    output:
    // one-row results table (h2, gcov, rg, intercepts) + GenomicSEM log
    tuple val(meta1), val(meta2), path("*.ldsc.tsv"),  emit: results
    tuple val(meta1), val(meta2), path("*_ldsc.log"),  emit: log

    // GenomicSEM version -> pipeline_info/nf_core_sherlock_software_versions.yml
    tuple val("${task.process}"),
          val('GenomicSEM'),
          eval("Rscript -e 'cat(as.character(packageVersion(\"GenomicSEM\")))'"),
          topic: versions,
          emit: versions_genomicsem

    // run unless modules.config sets ext.when = false
    when:
    task.ext.when == null || task.ext.when

    script:

    // ---- no conda build for GenomicSEM ----
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SHERLOCK_LDSC: GenomicSEM is not on conda. Use -profile docker or singularity.")
    }

    // ---- settings from conf/modules.config ----
    def args = ''
    if (task.ext.args) {
        args = task.ext.args
    }

    // ---- output name (default = pair id, e.g. AD_SCZ) ----
    def prefix = "${meta1.id}_${meta2.id}"
    if (task.ext.prefix) {
        prefix = task.ext.prefix
    }

    // ---- optional flags: prevalence only if the samplesheet cell is filled ----
    def optional = []
    if (meta1.pop_prev) { optional << "--pop-prev1 ${meta1.pop_prev}" }
    if (meta2.pop_prev) { optional << "--pop-prev2 ${meta2.pop_prev}" }
    optional = optional.join(' ')

    """
    sherlockLDSC.R \\
        --trait1        ${meta1.id} \\
        --sumstats1     ${sumstats1} \\
        --trait-type1   ${meta1.trait_type} \\
        --trait2        ${meta2.id} \\
        --sumstats2     ${sumstats2} \\
        --trait-type2   ${meta2.trait_type} \\
        --ld            ${ld_ref} \\
        --wld           ${ld_ref} \\
        --out           ${prefix} \\
        ${optional} \\
        ${args}
    """

    // ---- stub: fake output for quick dry runs (-stub) ----
    stub:
    def prefix = "${meta1.id}_${meta2.id}"
    if (task.ext.prefix) {
        prefix = task.ext.prefix
    }

    """
    touch ${prefix}.ldsc.tsv
    touch ${prefix}_ldsc.log
    """
}
