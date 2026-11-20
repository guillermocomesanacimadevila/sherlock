process SHERLOCK_PLEIO {

    tag "${meta1.id}_${meta2.id}"
    label 'process_medium'

    container "ghcr.io/guillermocomesanacimadevila/sherlock-cfdr.pleio:0.0.0.9100-76d5085"

    input:
    tuple val(meta1), path(sumstats1), val(meta2), path(sumstats2)
    path ldsc_ref // eur_w_ld_chr folder (same as the rg step's params.ldsc_ref)

    output:
    // AD_SCZ.pleio.all.tsv  (all SNPs)
    // AD_SCZ.pleio.pass.tsv (pleio_p < --thresh -> loci step)
    tuple val(meta1), val(meta2), path("*.pleio.all.tsv"),  emit: pleio_all
    tuple val(meta1), val(meta2), path("*.pleio.pass.tsv"), emit: pleio_pass
    tuple val(meta1), val(meta2), path("*_pleio.log"),      emit: log

    // PLEIO version -> pipeline_info (no tags upstream -> pinned commit)
    tuple val("${task.process}"),
          val('pleio'),
          val('b18d507'),
          topic: versions,
          emit: versions_pleio

    // run unless modules.config sets ext.when = false
    when:
    task.ext.when == null || task.ext.when

    script:

    // ---- no conda build for PLEIO ----
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SHERLOCK_PLEIO: PLEIO is not on conda. Use -profile docker or singularity.")
    }

    // ---- settings from conf/modules.config (--nis, --blup, --thresh, --trait-p) ----
    def args = ''
    if (task.ext.args) {
        args = task.ext.args
    }

    // ---- output name (default = pair id, e.g. AD_SCZ) ----
    def prefix = "${meta1.id}_${meta2.id}"
    if (task.ext.prefix) {
        prefix = task.ext.prefix
    }

    // ---- pop_prev only if set (needed for cc traits) ----
    def pop_prev1 = ''
    if (meta1.pop_prev) { pop_prev1 = "--pop-prev1 ${meta1.pop_prev}" }
    def pop_prev2 = ''
    if (meta2.pop_prev) { pop_prev2 = "--pop-prev2 ${meta2.pop_prev}" }

    """
    # keep python prints in order with the tools' output in the log
    export PYTHONUNBUFFERED=1

    sherlockPLEIO.py \\
        --trait1      ${meta1.id} \\
        --sumstats1   ${sumstats1} \\
        --trait-type1 ${meta1.trait_type} \\
        ${pop_prev1} \\
        --trait2      ${meta2.id} \\
        --sumstats2   ${sumstats2} \\
        --trait-type2 ${meta2.trait_type} \\
        ${pop_prev2} \\
        --ldsc-ref    ${ldsc_ref} \\
        --cores       ${task.cpus} \\
        --out         ${prefix} \\
        ${args} \\
        > ${prefix}_pleio.log 2>&1
    """

    // ---- stub: fake output for quick dry runs (-stub) ----
    stub:
    def prefix = "${meta1.id}_${meta2.id}"
    if (task.ext.prefix) {
        prefix = task.ext.prefix
    }

    """
    touch ${prefix}.pleio.all.tsv
    touch ${prefix}.pleio.pass.tsv
    touch ${prefix}_pleio.log
    """
}