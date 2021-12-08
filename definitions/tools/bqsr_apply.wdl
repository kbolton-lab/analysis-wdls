version 1.0

task bqsrApply {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bam_bai
    String output_name = "final"
    Array[File] known_sites
    Array[File] known_sites_tbi  # secondaryFiles...
    Array[String] intervals = ["chr1", "chr2", "chr3", "chr4", "chr5","chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
  }

  Int cores = 1
  Int space_needed_gb = 10 + round(size(known_sites, "GB") + size(known_sites_tbi, "GB") + size([reference, reference_fai, reference_dict], "GB") + size([bam, bam_bai], "GB") * 2)
  Int preemptible = 1
  Int maxRetries = 0
  runtime {
    cpu: cores
    docker: "broadinstitute/gatk:4.1.8.1"
    memory: "18GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    preemptible: preemptible
    maxRetries: maxRetries
  }

  command <<<
    /gatk/gatk --java-options -Xmx16g BaseRecalibrator -O bqsr.table ~{sep=" " prefix("-L ", intervals)} -R ~{reference} -I ~{bam} ~{sep=" " prefix("--known-sites ", known_sites)}
    /gatk/gatk --java-options -Xmx16g ApplyBQSR -O ~{output_name}.bam ~{sep=" " prefix("--static-quantized-quals ", [10, 20, 30])} -R ~{reference} -I ~{bam} -bqsr bqsr.table
  >>>

  output {
    File bqsr_bam = "~{output_name}.bam"
    File bqsr_bam_bai = "~{output_name}.bai"
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bam_bai
    String output_name = "final"

    Array[File] known_sites
    Array[File] known_sites_tbi
    Array[String]? intervals
  }
  call bqsrApply {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    bam_bai=bam_bai,
    output_name=output_name,
    known_sites=known_sites,
    known_sites_tbi=known_sites_tbi,
    intervals=intervals
  }
}
