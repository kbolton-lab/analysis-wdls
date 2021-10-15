version 1.0

task realign {
    input {
        File bam
        File reference
        File reference_fai
        File reference_dict
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_pac
        File reference_sa
    }

    Int cores = 8
    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_amb, reference_ann, reference_bwt, reference_pac, reference_sa], "GB")

    runtime {
        docker: "mgibio/dna-alignment:1.0.0"
        memory: "16GB"
        cpu: cores
        bootDiskSizeGb: 10 + round(5*data_size + reference_size)
        disks: "local-disk ~{10 + round(5*data_size + reference_size)} SSD"
    }

    command <<<
        /bin/bash /usr/bin/umi_realignment.sh ~{bam} ~{reference} ~{cores}
    >>>

    output {
        File consensus_aligned_bam = "realigned.bam"
        File consensus_aligned_bam_bai = "realign.bam.bai"
    }
}

workflow wf {
    input {
        File bam
        File reference
        File reference_fai
        File reference_dict
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_pac
        File reference_sa
    }

    call realign {
        input:
        bam = bam,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        reference_amb = reference_amb,
        reference_ann = reference_ann,
        reference_bwt = reference_bwt,
        reference_pac = reference_pac,
        reference_sa = reference_sa
    }
}
