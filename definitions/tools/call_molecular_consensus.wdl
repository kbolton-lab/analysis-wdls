version 1.0

task callMolecularConsensus {
    input {
        File bam
    }

    Int cores = 2
    Int space_needed_gb = 10 + round(10*size(bam, "GB"))

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: "6GB"
        cpu: cores
        bootDiskSizeGb: 25
        disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        /usr/local/bin/fgbio CallMolecularConsensusReads --input ~{bam} --error-rate-pre-umi 45 --error-rate-post-umi 30 --min-input-base-quality 30 --min-reads 1 --output consensus_unaligned.bam
    >>>

    output {
        File consensus_bam = "consensus_unaligned.bam"
    }
}

workflow wf {
    input {
        File bam
    }

    call callMolecularConsensus {
        input:
        bam = bam
    }
}
