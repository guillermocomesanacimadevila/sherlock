process SHERLOCK_LAVA {

    tag "${meta1.id}_${meta2.id}"
    label 'process_medium'

    container "ghcr.io/guillermocomesanacimadevila/sherlock-lava:0.1.5-e729a24"

    input:
    tuple val(meta1), path(sumstats1), val(meta2), path(sumstats2)
    path lava_ref  // LAVA reference folder (lava-ukb-v1.1_chr*.bcor/.info)
    path lava_loci // locus definition file (LOC CHR START STOP), same build as lava_ref

    output:
    // local h2 per locus per trait + local rg where both traits pass --univ-thresh
    tuple val(meta1), val(meta2), path("*.lava_univ.tsv"),  emit: univ
    tuple val(meta1), val(meta2), path("*.lava_bivar.tsv"), emit: bivar
    tuple val(meta1), val(meta2), path("*_lava.log"),       emit: log

    // LAVA version -> pipeline_info
    tuple val("${task.process}"),
          val('LAVA'),
          eval("Rscript -e 'cat(as.character(packageVersion(\"LAVA\")))'"),
          topic: versions,
          emit: versions_lava

    // run unless modules.config sets ext.when = false
    when:
    task.ext.when == null || task.ext.when

    script:

    // ---- no conda build for LAVA ----
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SHERLOCK_LAVA: LAVA is not on conda. Use -profile docker or singularity.")
    }

    // ---- settings from conf/modules.config (e.g. --univ-thresh) ----
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

    // ref prefix from the .info files (e.g. lava-ukb-v1.1_chr1.info -> lava-ukb-v1.1)
    """
    ref=\$(ls ${lava_ref}/*_chr*.info | head -n 1 | sed 's/_chr[0-9]*\\.info\$//')

    sherlockLAVA.R \\
        --trait1        ${meta1.id} \\
        --sumstats1     ${sumstats1} \\
        --trait-type1   ${meta1.trait_type} \\
        --trait2        ${meta2.id} \\
        --sumstats2     ${sumstats2} \\
        --trait-type2   ${meta2.trait_type} \\
        --ref           \$ref \\
        --loci          ${lava_loci} \\
        --cores         ${task.cpus} \\
        --out           ${prefix} \\
        ${optional} \\
        ${args} \\
        > ${prefix}_lava.log 2>&1
    """

    // ---- stub: fake output for quick dry runs (-stub) ----
    stub:
    def prefix = "${meta1.id}_${meta2.id}"
    if (task.ext.prefix) {
        prefix = task.ext.prefix
    }

    """
    touch ${prefix}.lava_univ.tsv
    touch ${prefix}.lava_bivar.tsv
    touch ${prefix}_lava.log
    """
}
