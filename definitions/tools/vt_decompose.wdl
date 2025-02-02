version 1.0

task vtDecompose {
  input {
    File vcf
    File vcf_tbi
  }

  Int space_needed_gb = 5 + round(size([vcf, vcf_tbi], "GB")*2)
  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    memory: "6GB"
    docker: "kboltonlab/vt"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  command <<<
    /opt/vt/vt decompose -s -o decomposed.vcf.gz ~{vcf}

    /usr/bin/tabix decomposed.vcf.gz
  >>>

  output {
    File decomposed_vcf = "decomposed.vcf.gz"
    File decomposed_vcf_tbi = "decomposed.vcf.gz.tbi"
  }
}

workflow wf {
  input {
    File vcf
    File vcf_tbi
  }

  call vtDecompose {
    input:
    vcf=vcf,
    vcf_tbi=vcf_tbi
  }
}
