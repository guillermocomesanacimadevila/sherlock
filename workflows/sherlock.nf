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
include { SHERLOCK_LDSC          } from '../modules/local/sherlock/ldsc'
include { SHERLOCK_HDL           } from '../modules/local/sherlock/hdl'
include { SHERLOCK_HDLL          } from '../modules/local/sherlock/hdll'
include { SHERLOCK_LAVA          } from '../modules/local/sherlock/lava'
include { SHERLOCK_CONJFDR       } from '../modules/local/sherlock/conjfdr'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// [ id1, id2 ] per pair + [ meta, sumstats ] per trait  ->  [ meta1, sumstats1, meta2, sumstats2 ] per pair
def pairFiles(ch_pairs, ch_files) {
    def ch_by_id = ch_files.map { meta, sumstats -> [ meta.id, meta, sumstats ] }
    return ch_pairs
        .combine(ch_by_id, by: 0)                                                    // match trait 1
        .map { _id1, id2, meta1, sumstats1 -> [ id2, meta1, sumstats1 ] }
        .combine(ch_by_id, by: 0)                                                    // match trait 2
        .map { _id2, meta1, sumstats1, meta2, sumstats2 -> [ meta1, sumstats1, meta2, sumstats2 ] }
}

workflow SHERLOCK {

    take:
    ch_traits // channel: [ meta, sumstats ], one per trait
    ch_pairs  // channel: [ id1, id2 ], one per trait pair
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

    ch_traits
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

    def steps      = params.steps.tokenize(',')
    def rg_methods = params.rg_methods.tokenize(',')

    // munged sumstats + rg results, stay empty if rg (or that method) isn't run
    def ch_munged  = channel.empty()
    def ch_ldsc    = channel.empty()
    def ch_hdl     = channel.empty()
    def ch_hdll    = channel.empty()
    def ch_lava    = channel.empty()
    def ch_conjfdr = channel.empty()

    if (steps.contains('rg')) {

        // ---- LDSC: munge -> rg, needs hm3 + ldsc_ref ----
        if (rg_methods.contains('ldsc')) {

            if (!params.hm3) {
                error("rg_methods includes 'ldsc' but hm3 is not set: add the HapMap3 SNP list (w_hm3.snplist) to params.yaml.")
            }
            if (!params.ldsc_ref) {
                error("rg_methods includes 'ldsc' but ldsc_ref is not set: add the LD score folder (eur_w_ld_chr/) to params.yaml.")
            }
            def ch_hm3      = channel.value(file(params.hm3, checkIfExists: true))
            def ch_ldsc_ref = channel.value(file(params.ldsc_ref, checkIfExists: true))

            SHERLOCK_MUNGE(SHERLOCK_QC.out.sumstats, ch_hm3)
            ch_munged = SHERLOCK_MUNGE.out.sumstats

            // one tuple per pair: [ meta1, sumstats1, meta2, sumstats2 ] (trait names sorted, e.g. AD, SCZ)
            pairFiles(ch_pairs, ch_munged).set { ch_pair_munged }

            SHERLOCK_LDSC(ch_pair_munged, ch_ldsc_ref)
            ch_ldsc = SHERLOCK_LDSC.out.results
        }

        // ---- HDL: reads QC output directly (no munge), needs hdl_ref ----
        if (rg_methods.contains('hdl')) {

            if (!params.hdl_ref) {
                error("rg_methods includes 'hdl' but hdl_ref is not set: add the HDL panel folder to params.yaml, or set rg_methods: ldsc to skip HDL.")
            }
            def ch_hdl_ref = channel.value(file(params.hdl_ref, checkIfExists: true))

            // same pairing as LDSC, but on the QC'd sumstats
            pairFiles(ch_pairs, SHERLOCK_QC.out.sumstats).set { ch_pair_qc }

            SHERLOCK_HDL(ch_pair_qc, ch_hdl_ref)
            ch_hdl = SHERLOCK_HDL.out.results
        }
    }

    //
    // LOCAL GENETIC CORRELATION (lrg)
    //

    def lrg_methods = params.lrg_methods.tokenize(',')

    if (steps.contains('lrg')) {

        // one tuple per pair on the QC'd sumstats (shared by HDL-L and LAVA)
        def ch_pair_qc_lrg = pairFiles(ch_pairs, SHERLOCK_QC.out.sumstats)

        // ---- HDL-L: reads QC output directly (no munge), needs hdll_ref ----
        if (lrg_methods.contains('hdll')) {

            if (!params.hdll_ref) {
                error("lrg_methods includes 'hdll' but hdll_ref is not set: add the HDL-L reference folder (with LD.path/ and bim.path/) to params.yaml.")
            }
            def ch_hdll_ref = channel.value(file(params.hdll_ref, checkIfExists: true))

            SHERLOCK_HDLL(ch_pair_qc_lrg, ch_hdll_ref)
            ch_hdll = SHERLOCK_HDLL.out.results
        }

        // ---- LAVA: reads QC output directly (no munge), needs lava_loci, lava_ref ----
        if (lrg_methods.contains('lava')) {

            if (!params.lava_ref) {
                error("lrg_methods includes 'lava' but lava_ref is not set: add the LAVA reference folder (with <prefix>_chr*.bcor/.info files, e.g. lava-ukb-v1.1) to params.yaml.")
            }

            if (!params.lava_loci) {
                error("lrg_methods includes 'lava' but lava_loci is not set: add the LAVA locus file (LOC CHR START STOP, same build as lava_ref, e.g. blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile) to params.yaml.")
            }

            def ch_lava_ref  = channel.value(file(params.lava_ref,  checkIfExists: true))
            def ch_lava_loci = channel.value(file(params.lava_loci, checkIfExists: true))

            // sherlock lava
            SHERLOCK_LAVA(
                ch_pair_qc_lrg,
                ch_lava_ref,
                ch_lava_loci
            )
            ch_lava = SHERLOCK_LAVA.out.bivar
        }
    }

    //
    // PLEIOTROPY (pleio)
    //

    def pleio_methods = params.pleio_methods.tokenize(',')

    if (steps.contains('pleio')) {

        // one tuple per pair on the QC'd sumstats (shared by conjFDR and PLEIO later)
        def ch_pair_qc_pleio = pairFiles(ch_pairs, SHERLOCK_QC.out.sumstats)

        // ----- conjFDR
        if (pleio_methods.contains('conjfdr')) {

            if (!params.conjfdr_ref) {
                error("pleio_methods includes 'conjfdr' but conjfdr_ref is not set: add the cfdr.pleio reference folder (with all_chr_perVariant.rds and chr01..22_LDpairs.rds, unzipped genref.zip from Zenodo 5750318) to params.yaml.")
            }

            def ch_conjfdr_ref = channel.value(file(params.conjfdr_ref, checkIfExists: true))

            // execute conjFDR module
            SHERLOCK_CONJFDR(
                ch_pair_qc_pleio,
                ch_conjfdr_ref
            )
            ch_conjfdr = SHERLOCK_CONJFDR.out.cfdr_pass
        }
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
    ldsc           = ch_ldsc                     // channel: [ meta1, meta2, path(*.ldsc.tsv) ]
    hdl            = ch_hdl                      // channel: [ meta1, meta2, path(*.hdl.tsv) ]
    hdll           = ch_hdll                     // channel: [ meta1, meta2, path(*.hdll.tsv) ]
    lava           = ch_lava                     // channel: [ meta1, meta2, path(*.lava_bivar.tsv) ]
    conjfdr        = ch_conjfdr                  // channel: [ meta1, meta2, path(*.conjfdr.pass.tsv) ]
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
