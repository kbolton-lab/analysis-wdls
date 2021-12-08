version 1.0

task fgbioFastqToBam {
    input {
        File trimmed_fastqs
        String sample_name
        String library_name
        String platform_unit
        String platform
    }

    Int cores = 1
    Int space_needed_gb = 10 + 2*round(size([trimmed_fastqs], "GB"))
    Int preemptible = 1
    Int maxRetries = 0
    runtime {
      cpu: cores
      docker: "quay.io/biocontainers/fgbio:1.3.0--0"
      memory: "6GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        /usr/local/bin/fgbio FastqToBam --input ~{trimmed_fastqs} --output output.bam --sample ~{sample_name} --library ~{library_name} --platform ~{platform} --platform-unit ~{platform_unit}
    >>>

    output {
        File bam = "output.bam"
    }
}
