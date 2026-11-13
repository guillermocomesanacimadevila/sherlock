/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_sherlock_pipeline'
include { SHERLOCK_QC            } from '../modules/local/sherlock/qc'
include { SHERLOCK_MUNGE         } from '../modules/local/sherlock/munge'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SHERLOCK {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    outdir

    main:

    def ch_versions = channel.empty()

    //
    // QC
    //

    // check liftover setup before running
    if (params.chain_file) {

        def chain_name = file(params.chain_file).name
        def expected   = 'ToHg38'
        if (params.target_build == 'GRCh37') {
            expected = 'ToHg19'
        }

        if (!chain_name.contains(expected)) {
            log.warn "chain_file '${chain_name}' doesn't look like a liftover to ${params.target_build} (expected '*${expected}*'). Check it's the right direction."
        }
    }

    ch_samplesheet
        .map { meta, sumstats  ->
            if (meta.build != params.target_build && !params.chain_file) {
                error("${meta.id} is ${meta.build} but target_build is ${params.target_build}: set chain_file in params.yaml for liftover.")
            }
            [ meta, sumstats ]
        }
        .set { ch_gwas }

    // chain file for the module: [] = none
    def ch_chain = channel.value([])
    if (params.chain_file) {
        ch_chain = channel.value(file(params.chain_file, checkIfExists: true))
    }

    SHERLOCK_QC(ch_gwas, ch_chain)

    //
    // GENETIC CORRELATION (rg)
    //

    def steps = params.steps.tokenize(',')

    // munged sumstats, stays empty if rg isn't run
    def ch_munged = channel.empty()

    if (steps.contains('rg')) {

        // hm3 is required for munge
        if (!params.hm3) {
            error("--steps includes 'rg' but hm3 is not set: add the HapMap3 SNP list (w_hm3.snplist) to params.yaml.")
        }
        def ch_hm3 = channel.value(file(params.hm3, checkIfExists: true))

        SHERLOCK_MUNGE(SHERLOCK_QC.out.sumstats, ch_hm3)
        ch_munged = SHERLOCK_MUNGE.out.sumstats
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'sherlock_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
    emit:
    qc             = SHERLOCK_QC.out.sumstats    // channel: [ meta, path(*.qc.tsv.gz) ]
    munged         = ch_munged                   // channel: [ meta, path(*.sumstats.gz) ]
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
