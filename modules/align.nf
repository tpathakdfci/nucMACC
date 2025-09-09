process alignment{

 
  label 'bigCPU'
  memory { params.genomeSize > 200000000 ? params.high_memory : params.low_memory}
  publishDir "${params.outDir}/QC/02_ALIGNMENT/${sampleID}", mode: 'copy', pattern: "*_alignment_stats.txt"
  publishDir "${params.outDir}/RUN/00_ALIGNMENT/${sampleID}", mode: 'copy', pattern: "*_aligned.bam"


  input:
  tuple val(sampleID), file(read1), file(read2)
  path genome
  path genome_tarball

  output:
  file "*_alignment_stats.txt"
  tuple val(sampleID), file("*_aligned.bam"), file("*_aligned.bam.bai")
  tuple val(sampleID), file("*_alignment_stats.txt")

   script:
  """

  mkdir genome_idx
  tar -xzvf ${genome_tarball} -C genome_idx


  bowtie2 -t \\
  --threads $task.cpus \\
  --very-sensitive-local \\
  --no-discordant \\
  --no-mix \\
  --dovetail \\
  -x genome_idx/${genome.baseName} \\
  -1 $read1 \\
  -2 $read2 \\
  2> ${sampleID}_alignment_stats.txt \\
  | samtools view -bS -q 30 -f 2 -@ $task.cpus - | samtools sort -@ $task.cpus - > ${sampleID}_aligned.bam

  samtools index -b ${sampleID}_aligned.bam
  """

}
