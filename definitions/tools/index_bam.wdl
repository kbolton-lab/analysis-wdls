version 1.0

task indexBam {
  input { File bam }

  Int space_needed_gb = 10 + round(size(bam, "GB")*2)
  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    cpu: cores
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    memory: "6GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    preemptible: preemptible
    maxRetries: maxRetries
  }
  command <<<
    mv ~{bam} ~{basename(bam)}
    /usr/local/bin/samtools index ~{basename(bam)} ~{basename(bam)}.bai
  >>>
  output {
    File indexed_bam = basename(bam)
    File indexed_bam_bai = "~{basename(bam)}.bai"
  }
}

workflow wf {
  input { File bam }
  call indexBam { input: bam=bam }
}
