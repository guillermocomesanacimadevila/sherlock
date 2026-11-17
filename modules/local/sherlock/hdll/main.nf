process SHERLOCK_HDLL {

    tag "${meta1.id}_${meta2.id}"
    label 'process_medium'

    // avoiding conda here (** for now **)
    container "ghcr.io/guillermocomesanacimadevila/sherlock-hdl:1.4.3-e94ba00"

    input:
    tuple val(meta1), path(sumstats1), val(meta2), path(sumstats2)
    path hdll_ref // HDL-L reference folder (contains LD.path/ and bim.path/)

    output:
    // one row per locus (local h2, gcov, rg, 95% interval, p)
    tuple val(meta1), val(meta2), path("*.hdll.tsv"),  emit: results
    tuple val(meta1), val(meta2), path("*_hdll.log"),  emit: log

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
        error("SHERLOCK_HDLL: HDL is not on conda. Use -profile docker or singularity.")
    }

    // ---- settings from conf/modules.config (e.g. --n0) ----
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
    sherlockHDL-L.R \\
        --trait1        ${meta1.id} \\
        --sumstats1     ${sumstats1} \\
        --trait2        ${meta2.id} \\
        --sumstats2     ${sumstats2} \\
        --ld            ${hdll_ref}/LD.path \\
        --bim           ${hdll_ref}/bim.path \\
        --cores         ${task.cpus} \\
        --out           ${prefix} \\
        ${args} \\
        > ${prefix}_hdll.log 2>&1
    """

    // ---- stub: fake output for quick dry runs (-stub) ----
    stub:
    def prefix = "${meta1.id}_${meta2.id}"
    if (task.ext.prefix) {
        prefix = task.ext.prefix
    }

    """
    touch ${prefix}.hdll.tsv
    touch ${prefix}_hdll.log
    """
}