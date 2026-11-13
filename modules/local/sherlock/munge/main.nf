process SHERLOCK_MUNGE {

    tag "${meta.id}"
    label 'process_low'

    container "ghcr.io/guillermocomesanacimadevila/sherlock-genomicsem:0.0.5-da95d43"

    input:
    tuple val(meta), path(sumstats) // one trait: [ samplesheet row, <trait>.qc.tsv.gz ]
    path hm3

    // GenomicSEM version -> pipeline_info/nf_core_sherlock_software_versions.yml
    output:

    tuple val(meta), path("*.sumstats.gz"), emit: sumstats
    tuple val(meta), path("*_munge.log"),   emit: log

    tuple val("${task.process}"),
          val('GenomicSEM'),
          eval("Rscript -e 'cat(as.character(packageVersion(\"GenomicSEM\")))'"),
          topic: versions,
          emit: versions_genomicsem


    // when
    when:
    task.ext.when == null || task.ext.when

    script:

    // if user selects -> conda / mamba == cooked!
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SHERLOCK_MUNGE: GenomicSEM is not on conda. Use -profile docker or singularity.")
    }

    def args = ''
    if (task.ext.args) {
        args = task.ext.args
    }

    def prefix = meta.id
    if (task.ext.prefix) {
        prefix = task.ext.prefix
    }

    """
    sherlockMunge.R \\
        --sumstats  ${sumstats} \\
        --trait     ${prefix} \\
        --hm3       ${hm3} \\
        ${args}
    """

    // stub
    stub:
    def prefix = meta.id
    if (task.ext.prefix) {
        prefix = task.ext.prefix
    }

    """
    echo "" | gzip > ${prefix}.sumstats.gz
    touch ${prefix}_munge.log
    """
}
