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

 //                      workflow

/// Read the CSV file once and create all required channels
Channel
    .fromPath(params.csvInput)
    .splitCsv(header: true)
    .set { samples_ch }

// Create channels based on analysis type 
if (params.bamEntry) {
    samples_ch
        .map { row -> tuple(row.Sample_Name, file(row.path_mono)) }
        .set { bamEntry_mono }

    samples_ch
        .map { row -> tuple(row.Sample_Name, file(row.path_sub)) }
        .set { bamEntry_sub }
} else {
    // creates the paired-end channel for alignment
    samples_ch
        .map { row -> tuple(row.Sample_Name, file(row.path_fwdReads), file(row.path_revReads)) }
        .set { samplePair_ch }

    samplePair_ch.view() { it -> "Item in samplePair_ch: ${it}" }

    // creates the single-end channel for fastqc
    samples_ch
        .map { row -> tuple(row.Sample_Name, file(row.path_fwdReads)) }
        .mix(samples_ch.map { row -> tuple(row.Sample_Name, file(row.path_revReads)) })
        .set { sampleSingle_ch }
}

//  Create the MNase concentration channel 
samples_ch
    .map { row -> tuple(row.MNase_U.toDouble(), row.Sample_Name) }
    .set { samples_conc }



// load workflows
// generate profiles
include{MNaseQC} from './workflows/MNaseQC'
include{sub_bamEntry; sub_FASTQ_entry; common_nucMACC} from './workflows/nucMACC'


// Check mandatory parameters
if (params.csvInput) { ch_csv = file(params.csvInput) } else { exit 1, 'Input samplesheet not found!' }

if (params.TSS) {ch_TSS = file (params.TSS)}
      if(params.TSS){
        if (ch_TSS.isEmpty()) { exit 1, 'TSS file not found!'}
        }

if (params.blacklist) {ch_blacklist = file (params.blacklist)}
      if(params.blacklist){
        if (ch_blacklist.isEmpty()) { exit 1, 'Blacklist file not found!'}
        }


  Channel.fromPath(params.genome).first().set { genome_ch }
  Channel.fromPath(params.genomeIdx).set{genome_tarball_ch}


workflow{
        if(params.analysis=='MNaseQC'){
                MNaseQC(sampleSingle_ch,samplePair_ch,samples_conc)
        }
        if(params.analysis=='nucMACC'){
                if (params.bamEntry == true) {
                  sub_bamEntry(bamEntry_mono,bamEntry_sub)
                  common_nucMACC(sub_bamEntry.out[0], sub_bamEntry.out[1], samples_conc)
                  }
                else {
                  sub_FASTQ_entry(sampleSingle_ch,samplePair_ch,samples_conc, genome_ch, genome_tarball_ch)
                  common_nucMACC(sub_FASTQ_entry.out[0], sub_FASTQ_entry.out[1], samples_conc, genome_ch)
                }
        }
}