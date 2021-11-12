version 1.0

task samtoolsFlagstat {
  input {
    File bam
    File bam_bai
  }

  Int space_needed_gb = 5 + round(size([bam, bam_bai], "GB")*2)
  Int preemptible = 1
  Int maxRetries = 0
  Int cores = 1
  runtime {
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  String outfile = basename(bam) + ".flagstat"
  command <<<
    /usr/local/bin/samtools flagstat ~{bam} > ~{outfile}
  >>>

  output {
    File flagstats = outfile
  }
}
