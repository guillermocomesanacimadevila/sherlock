process SHERLOCK_HDL {

    tag "${meta1.id}_${meta2.id}"
    label 'process_low'

    // avoiding conda here (** for now **)
    container "ghcr.io/guillermocomesanacimadevila/sherlock-hdl:1.4.3-e94ba00"

    input:
    tuple val(meta1), path(sumstats1), val(meta2), path(sumstats2)
    path hdl_ref // // HDL reference panel folder (e.g. UKB_imputed_SVD_eigen99_extraction/)

    output:
    // one-row results table (h2, gcov, rg, eigen_use) + HDL log
    tuple val(meta1), val(meta2), path("*.hdl.tsv"),  emit: results
    tuple val(meta1), val(meta2), path("*_hdl.log"),  emit: log

    // HDL version -> pipeline_info/nf_core_sherlock_software_versions.yml
    tuple val("${task.process}"),
          val('HDL'),
          eval("Rscript -e 'cat(as.character(packageVersion(\"HDL\")))'"),
          topic: versions,
          emit: versions_hdl

    // run unless modules.config sets ext.when = false
    when:
    task.ext.when == null || task.ext.when

    script:

    // ---- no conda build for HDL ----
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SHERLOCK_HDL: HDL is not on conda. Use -profile docker or singularity.")
    }

    // ---- settings from conf/modules.config (e.g. --n-ref) ----
    def args = ''
    if (task.ext.args) {
        args = task.ext.args
    }

    // ---- output name (default = pair id, e.g. AD_SCZ) ----
    def prefix = "${meta1.id}_${meta2.id}"
    if (task.ext.prefix) {
        prefix = task.ext.prefix
    }

    """
    sherlockHDL.R \\
        --trait1        ${meta1.id} \\
        --sumstats1     ${sumstats1} \\
        --trait2        ${meta2.id} \\
        --sumstats2     ${sumstats2} \\
        --ref           ${hdl_ref} \\
        --out           ${prefix} \\
        ${args}
    """

    // ---- stub: fake output for quick dry runs (-stub) ----
    stub:
    def prefix = "${meta1.id}_${meta2.id}"
    if (task.ext.prefix) {
        prefix = task.ext.prefix
    }

    """
    touch ${prefix}.hdl.tsv
    touch ${prefix}_hdl.log
    """
}