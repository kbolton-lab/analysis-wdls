version 1.0

task bqsr {
  input {
    File reference
    File reference_fai   # secondaryFiles...
    File reference_dict  # secondaryFiles...

    File bam
    File bam_bai  # secondaryFiles...

    Array[File] known_sites
    Array[File] known_sites_tbi  # secondaryFiles...

    Array[String] intervals = ["chr1", "chr2", "chr3", "chr4", "chr5","chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
  }

  Int cores = 1
  Int space_needed_gb = 10 + round(size(flatten([known_sites,known_sites_tbi]), "GB") + size([reference, reference_fai, reference_dict, bam, bam_bai], "GB"))
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

  String outfile = "bqsr.table"
  command <<<
    /gatk/gatk --java-options -Xmx16g BaseRecalibrator -O ~{outfile} ~{sep=" " prefix("-L ", intervals)} -R ~{reference} -I ~{bam} ~{sep=" " prefix("--known-sites ", known_sites)}
  >>>

  # TODO: how much space to allocate for bqsr_table? what file does it scale with?
  output {
    File bqsr_table = outfile
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict

    File bam
    File bam_bai

    Array[File] known_sites
    Array[File] known_sites_tbi

    Array[String]? intervals
  }
  call bqsr {
    input:
    reference = reference,
    reference_fai = reference_fai,
    reference_dict = reference_dict,
    bam = bam,
    bam_bai = bam_bai,
    known_sites = known_sites,
    known_sites_tbi = known_sites_tbi,
    intervals = intervals
  }
}
