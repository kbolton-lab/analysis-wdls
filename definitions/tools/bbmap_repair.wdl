version 1.0

task bbmapRepair {
    input {
        File fastq1
        File fastq2
    }

    Int cores = 1
    Float data_size = size([fastq1, fastq2], "GB")
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/bbmap:38.92--he522d1c_0"
        memory: "6GB"
        cpu: cores
        bootDiskSizeGb: 10 + round(2*data_size)
        disks: "local-disk ~{10 + round(2*data_size)} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        repair.sh -Xmx6g repair=t overwrite=true interleaved=false outs=singletons.fq out1=R1.fixed.fastq.gz out2=R2.fixed.fastq.gz in1=~{fastq1} in2=~{fastq2}
    >>>

    output {
        Array[File] fastqs = ["R1.fixed.fastq.gz", "R2.fixed.fastq.gz"]
        File fastq1_repair = "R1.fixed.fastq.gz"
        File fastq2_repair = "R2.fixed.fastq.gz"
    }
}

workflow wf {
    input {
        File fastq1
        File fastq2
    }
    call bbmapRepair {
        input:
        fastq1 = fastq1,
        fastq2 = fastq2
    }
}
