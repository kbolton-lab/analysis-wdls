version 1.0

task groupReads {
    input {
        File bam
        Boolean umi_paired = true
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
        set -eo pipefail

        PAIRED=~{umi_paired}
        BAM=~{bam}

        if [ "$PAIRED" == true ]; then
            /usr/local/bin/fgbio GroupReadsByUmi --strategy paired --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input $BAM --output umi_grouped.bam
        else
            /usr/local/bin/fgbio GroupReadsByUmi --strategy adjacency --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input $BAM --output umi_grouped.bam
        fi
    >>>

    output {
        File grouped_bam = "umi_grouped.bam"
    }
}


workflow wf {
    input {
        File bam
        Boolean umi_paired = true
    }
    call groupReads {
        input:
        bam = bam,
        umi_paired = umi_paired
    }
}
