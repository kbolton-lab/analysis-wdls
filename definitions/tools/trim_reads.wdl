version 1.0

task trimReads {
    input {
        File? fastq1
        File? fastq2
        String adapter_one
        String adapter_two
        Int? trimN = 3
    }

    Int cores = 1
    Int space_needed_gb = 10 + round(size([fastq1, fastq2], "GB"))
    Int preemptible = 1
    Int maxRetries = 0
    runtime {
      cpu: cores
      docker: "quay.io/biocontainers/cutadapt:3.5--py36hc5360cc_0"
      memory: "6GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        /usr/local/bin/cutadapt --interleaved -a ~{adapter_one} -A ~{adapter_two} ~{fastq1} ~{fastq2} | /usr/local/bin/cutadapt --interleaved -u ~{trimN} -u -~{trimN} -U ~{trimN} -U -~{trimN} -m 30 -o trimmed.fq.gz -
    >>>

    output {
        File trimmed_fastqs = "trimmed.fq.gz"
    }
}
