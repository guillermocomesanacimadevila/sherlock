process SHERLOCK_CONJFDR {

    tag "${meta1.id}_${meta2.id}"
    label 'process_single'

    container "ghcr.io/guillermocomesanacimadevila/sherlock-cfdr.pleio:0.0.0.9100-76d5085"

    input:
    tuple val(meta1), path(sumstats1), val(meta2), path(sumstats2)
    path conjfdr_ref // cfdr.pleio reference folder (all_chr_perVariant.rds + chrXX_LDpairs.rds)

    output:
    // AD_SCZ.conjfdr.all.tsv  (all SNPs)
    // AD_SCZ.conjfdr.pass.tsv (conj_fdr < --conjfdr-thresh -> loci step)
    tuple val(meta1), val(meta2), path("*.conjfdr.all.tsv"),  emit: cfdr_all
    tuple val(meta1), val(meta2), path("*.conjfdr.pass.tsv"), emit: cfdr_pass
    tuple val(meta1), val(meta2), path("*_conjfdr.log"),      emit: log

    // conjFDR version -> pipeline_info
    tuple val("${task.process}"),
          val('cfdr.pleio'),
          eval("Rscript -e 'cat(as.character(packageVersion(\"cfdr.pleio\")))'"),
          topic: versions,
          emit: versions_cfdr

    // run unless modules.config sets ext.when = false
    when:
    task.ext.when == null || task.ext.when

    script:

    // ---- no conda build for cfdr.pleio ----
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SHERLOCK_CONJFDR: cfdr.pleio is not on conda. Use -profile docker or singularity.")
    }

    // ---- settings from conf/modules.config (--n-iter, --seed, --maf-min, --exclude-mhc, --exclude-apoe, --conjfdr-thresh) ----
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
    sherlockConjFDR.R \\
        --trait1    ${meta1.id} \\
        --sumstats1 ${sumstats1} \\
        --trait2    ${meta2.id} \\
        --sumstats2 ${sumstats2} \\
        --ref       ${conjfdr_ref} \\
        --out       ${prefix} \\
        ${args} \\
        > ${prefix}_conjfdr.log 2>&1
    """

    // ---- stub: fake output for quick dry runs (-stub) ----
    stub:
    def prefix = "${meta1.id}_${meta2.id}"
    if (task.ext.prefix) {
        prefix = task.ext.prefix
    }

    """
    touch ${prefix}.conjfdr.all.tsv
    touch ${prefix}.conjfdr.pass.tsv
    touch ${prefix}_conjfdr.log
    """
}
