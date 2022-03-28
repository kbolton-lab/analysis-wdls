version 1.0

task fastQC {
    input {
        File bam
        File bam_bai
    }

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([bam, bam_bai], "GB")
    Int space_needed_gb = 1 + round(data_size)

    runtime {
        docker: "biocontainers/fastqc:v0.11.9_cv8"
        memory: "6GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String bamroot = basename(bam, ".bam")

    command <<<
        /usr/local/bin/fastqc ~{bam} -outdir $PWD
    >>>

    output {
        File fastqc = "~{bamroot}_fastqc.html"
    }
}

workflow wf {
    input {
        File bam
        File bam_bai
    }

    call fastQC {
        input:
        bam = bam,
        bam_bai = bam_bai
    }

    output {
        File fastqc = fastQC.fastqc
    }
}
