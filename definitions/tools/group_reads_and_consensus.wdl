version 1.0

task groupReadsAndConsensus {
    input {
        File bam
        Boolean umi_paired = true
    }

    Int cores = 1
    Int space_needed_gb = 5 + round(3*size(bam, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: "6GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -eo pipefail

        PAIRED=~{umi_paired}
        BAM=~{bam}

        if [ "$PAIRED" == true ]; then
            /usr/local/bin/fgbio GroupReadsByUmi --strategy paired --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input $BAM --output umi_grouped.bam
        else
            /usr/local/bin/fgbio GroupReadsByUmi --strategy adjacency --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input $BAM --output umi_grouped.bam
        fi
        /usr/local/bin/fgbio CallMolecularConsensusReads --input umi_grouped.bam --error-rate-pre-umi 45 --error-rate-post-umi 30 --min-input-base-quality 30 --min-reads 1 --output consensus_unaligned.bam
    >>>

    output {
        File consensus_bam = "consensus_unaligned.bam"
    }
}


workflow wf {
    input {
        File bam
        Boolean umi_paired = true
    }
    call groupReadsAndConsensus {
        input:
        bam = bam,
        umi_paired = umi_paired
    }
}
