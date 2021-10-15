version 1.0

task clipOverlap {
    input {
        File bam
        File reference
        File reference_fai
        File reference_dict
    }

    Int cores = 2
    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: "6GB"
        cpu: cores
        bootDiskSizeGb: 10 + round(5*data_size + reference_size)
        disks: "local-disk ~{10 + round(5*data_size + reference_size)} SSD"
    }

    command <<<
        /usr/local/bin/fgbio ClipBam --input ~{bam} --ref ~{reference} --clipping-mode Hard --clip-overlapping-reads true --output clipped.bam
    >>>

    output {
        File clipped_bam = "clipped.bam"
        File clipped_bam_bai = "clipped.bai"
    }
}

workflow wf {
    input {
        File bam
        File reference
        File reference_fai
        File reference_dict
    }

    call clipOverlap {
        input:
        bam = bam,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict
    }
}
