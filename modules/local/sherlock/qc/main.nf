process SHERLOCK_QC {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // container: add once the Seqera container exists (step 8)

    input:
    tuple val(meta), path(sumstats) // one trait: [ samplesheet row, sumstats file ]
    path chain_file // liftover chain file, or [] if none

    output:
    // QC'd sumstats
    tuple val(meta), path("*.qc.parquet"), emit: parquet

    // sherlock-tools version -> pipeline_info/nf_core_sherlock_software_versions.yml
    tuple val("${task.process}"),
          val('sherlock-tools'),
          eval("pip show sherlock-tools | sed -n 's/^Version: //p'"),
          topic: versions,
          emit: versions_sherlock_tools

    // run unless modules.config sets ext.when = false
    when:
    task.ext.when == null || task.ext.when

    script:

    // ---- settings from conf/modules.config ----
    def args = ''
    if (task.ext.args) {
        args = task.ext.args
    }

    // ---- output name (default = trait name, e.g. AD) ----
    def prefix = meta.id
    if (task.ext.prefix) {
        prefix = task.ext.prefix
    }

    // ---- optional flags: only added if the samplesheet cell is filled ----
    def optional = []
    if (meta.freq_col)   { optional << "--freq-col ${meta.freq_col}" }
    if (meta.info_col)   { optional << "--info-col ${meta.info_col}" }
    if (meta.n_col)      { optional << "--n-col ${meta.n_col}" }
    if (meta.n)          { optional << "--n ${meta.n}" }
    if (meta.n_cases)    { optional << "--n-cases ${meta.n_cases}" }
    if (meta.n_controls) { optional << "--n-controls ${meta.n_controls}" }
    if (chain_file)      { optional << "--chain-file ${chain_file}" }
    optional = optional.join(' ')

    """
    sherlockQC.py \\
        --sumstats      ${sumstats} \\
        --out           ${prefix}.qc.parquet \\
        --trait-type    ${meta.trait_type} \\
        --genome-build  ${meta.build} \\
        --snp-col       ${meta.snp_col} \\
        --chr-col       ${meta.chr_col} \\
        --pos-col       ${meta.pos_col} \\
        --a1-col        ${meta.a1_col} \\
        --a2-col        ${meta.a2_col} \\
        --beta-col      ${meta.beta_col} \\
        --se-col        ${meta.se_col} \\
        --pval-col      ${meta.pval_col} \\
        ${optional} \\
        ${args}
    """

    // ---- stub: fake output for quick dry runs (-stub) ----
    stub:
    def prefix = meta.id
    if (task.ext.prefix) {
        prefix = task.ext.prefix
    }
    """
    touch ${prefix}.qc.parquet
    """
}
