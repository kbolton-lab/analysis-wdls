version 1.0

task indexVcf {
  input {
    File vcf
  }

  Int space_needed_gb = 10 + round(3*size(vcf, "GB"))
  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    memory: "6GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  command <<<
    cp ~{vcf} ~{basename(vcf)}
    /usr/local/bin/tabix -p vcf ~{basename(vcf)}
  >>>
  output {
    File indexed_vcf = basename(vcf)
    File indexed_vcf_tbi = basename(vcf) + ".tbi"
  }
}

workflow wf {
  input { File vcf }
  call indexVcf {
    input: vcf=vcf
  }
}
