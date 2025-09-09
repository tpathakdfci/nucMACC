process qualimap {

  container 'uschwartz/qualimap'
  label 'mid'
  publishDir "${params.outDir}/QC/03_QUALIMAP/${sampleID}", mode: 'copy'

  input:
  tuple val(sampleID), file(bam), file(idx)

  output:
  file "${sampleID}"
  tuple val(sampleID), file("${sampleID}_qualimap_genome_results.txt")


  script:
  """
  qualimap bamqc --java-mem-size=14G -bam $bam -c -outdir ${sampleID}
  mv ${sampleID}/genome_results.txt ${sampleID}_qualimap_genome_results.txt
  """
}
