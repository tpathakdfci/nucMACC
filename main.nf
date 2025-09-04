#!/usr/bin/env nextflow

 /*
 ===============================================================================
                      nextflow based nucMACC pipeline
 ===============================================================================
Authors:
Uwe Schwartz <uwe.schwartz@ur.de>
 -------------------------------------------------------------------------------
 */

nextflow.enable.dsl = 2

 //                           show settings
 if (!params.help) {
         include{settings} from './modules/setting'
         settings()
 }

 //                       help message
 // Show help message
 if (params.help) {
     include{helpMessage} from './modules/help'
     helpMessage()
     exit 0
 }


// Check mandatory parameters

if (params.genome && params.genomes[params.genome]) {
    def genome_config = params.genomes[params.genome]
    params.genome = file(genome_config.genome, checkIfExists: false)
    params.genomeIdx = file(genome_config.genomeIdx, checkIfExists: false)
} else {
    exit 1, 'Error: You must provide a valid --genome key (e.g., GRCh38) that is configured in igenomes.config.'
}


if (params.csvInput) {
    ch_csv = file(params.csvInput)
} else {
    exit 1, 'Input samplesheet not found!'
}

if (params.TSS) {
    def ch_TSS = file(params.TSS)
    if (ch_TSS.isEmpty()) {
        exit 1, 'TSS file not found!'
    }
}

if (params.blacklist) {
    def ch_blacklist = file(params.blacklist)
    if (ch_blacklist.isEmpty()) {
        exit 1, 'Blacklist file not found!'
    }
}

if (params.analysis == 'MNaseQC' || params.analysis == 'nucMACC') {
    if (params.bamEntry == false) {
        if (!params.genome) {
            exit 1, 'Genome FASTA not found!'
        }
        if (!params.genomeIdx) {
            exit 1, 'Bowtie2 index folder not found!'
        }
    }
}

// Read csv file
if (params.test) {
    if (params.bamEntry) {
        Channel.fromPath(params.csvInput)
            .splitCsv(header: true)
            .map { row -> tuple(row.Sample_Name, file(row.path_mono)) }
            .set { bamEntry_mono }
        Channel.fromPath(params.csvInput)
            .splitCsv(header: true)
            .map { row -> tuple(row.Sample_Name, file(row.path_sub)) }
            .set { bamEntry_sub }
    } else {
        Channel.fromPath(params.csvInput)
            .splitCsv(header: true)
            .map { row -> tuple(row.Sample_Name, file(params.project.concat(row.path_fwdReads))) }
            .set { samples_fwd_ch }
        Channel.fromPath(params.csvInput)
            .splitCsv(header: true)
            .map { row -> tuple(row.Sample_Name, file(params.project.concat(row.path_revReads))) }
            .set { samples_rev_ch }
    }
} else if (params.bamEntry == true) {
    Channel.fromPath(params.csvInput)
        .splitCsv(header: true)
        .map { row -> tuple(row.Sample_Name, file(row.path_mono)) }
        .set { bamEntry_mono }
    Channel.fromPath(params.csvInput)
        .splitCsv(header: true)
        .map { row -> tuple(row.Sample_Name, file(row.path_sub)) }
        .set { bamEntry_sub }
    println "BamEntry csv part"
} else {
    Channel.fromPath(params.csvInput)
        .splitCsv(header: true)
        .map { row -> tuple(row.Sample_Name, file(row.path_fwdReads)) }
        .set { samples_fwd_ch }
    Channel.fromPath(params.csvInput)
        .splitCsv(header: true)
        .map { row -> tuple(row.Sample_Name, file(row.path_revReads)) }
        .set { samples_rev_ch }
}

if (params.bamEntry == false) {
    samples_fwd_ch.mix(samples_rev_ch).set { sampleSingle_ch }
    samples_fwd_ch.join(samples_rev_ch).set { samplePair_ch }
}

Channel.fromPath(params.csvInput)
    .splitCsv(header: true)
    .map { row -> tuple(row.MNase_U.toDouble(), row.Sample_Name) }
    .set { samples_conc }

include{MNaseQC} from './workflows/MNaseQC'
include{sub_bamEntry; sub_FASTQ_entry; common_nucMACC} from './workflows/nucMACC'

workflow {
    if (params.analysis == 'MNaseQC') {
        MNaseQC(sampleSingle_ch, samplePair_ch, samples_conc)
    }
    if (params.analysis == 'nucMACC') {
        if (params.bamEntry == true) {
            sub_bamEntry(bamEntry_mono, bamEntry_sub)
            common_nucMACC(sub_bamEntry.out[0], sub_bamEntry.out[1], samples_conc)
        } else {
            sub_FASTQ_entry(sampleSingle_ch, samplePair_ch, samples_conc)
            common_nucMACC(sub_FASTQ_entry.out[0], sub_FASTQ_entry.out[1], samples_conc)
        }
    }
}
