version 1.0

task annotateBam {
    input {
        File bam
        File fastq_with_umis
    }

    Int cores = 1
    Int space_needed_gb = 10 + 2*round(size([bam, fastq_with_umis], "GB"))
    Int preemptible = 1
    Int maxRetries = 0
    runtime {
      cpu: cores
      docker: "quay.io/biocontainers/fgbio:1.3.0--0"
      memory: "6GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        /usr/local/bin/fgbio AnnotateBamWithUmis --input ~{bam} --fastq ~{fastq_with_umis} --output umi_annotated.bam
    >>>

    output {
        File umi_bam = "umi_annotated.bam"
    }
}
