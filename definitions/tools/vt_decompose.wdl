version 1.0

task vtDecompose {
  input {
    File vcf
    File vcf_tbi
    String variant_caller = "caller"
  }

  Int space_needed_gb = 10 + round(size([vcf, vcf_tbi], "GB")*2)
  runtime {
    memory: "4GB"
    docker: "kboltonlab/vt"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    /opt/vt/vt decompose -s -o ~{variant_caller}.decomposed.vcf.gz ~{vcf}

    /usr/bin/local/tabix ~{variant_caller}.decomposed.vcf.gz
  >>>

  output {
    File decomposed_vcf = "~{variant_caller}.decomposed.vcf.gz"
    File decomposed_vcf_tbi = "~{variant_caller}.decomposed.vcf.gz.tbi"
  }
}

workflow wf {
  input {
    File vcf
    File vcf_tbi
    String variant_caller = "caller"
  }

  call vtDecompose {
    input:
    vcf=vcf,
    vcf_tbi=vcf_tbi,
    variant_caller=variant_caller
  }
}
