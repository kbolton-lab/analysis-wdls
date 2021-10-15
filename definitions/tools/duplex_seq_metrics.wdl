version 1.0

task duplexSeqMetrics {
    input {
        File bam
        File? target_intervals
        String description
        Boolean umi_paired = true
    }

    Int cores = 2
    Int space_needed_gb = 10 + round(10*size(bam, "GB"))

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: "6GB"
        cpu: cores
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        set -eo pipefail

        PAIRED=~{umi_paired}
        DESCRIPTION=~{description}
        INTERVALS=~{target_intervals}

        if [ "$PAIRED" == true ]; then
            if [[ -z "$DESCRIPTION" ]]; then
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input ~{bam} --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input ~{bam} --output duplex_seq.metrics --intervals ${INTERVALS}
                fi
            else
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input ~{bam} --description ${DESCRIPTION} --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input ~{bam} --description ${DESCRIPTION} --output duplex_seq.metrics --intervals ${INTERVALS}
                fi
            fi
        else
            echo "Sample not UMI Paired" > duplex_seq.metrics.txt
        fi
    >>>

    output {
        Array[File] duplex_seq_metrics = glob("duplex_seq.metrics.*")
    }
}

workflow wf {
    input {
        File bam
        File? target_intervals
        String description
        Boolean umi_paired = true
    }

    call duplexSeqMetrics {
        input:
        bam = bam,
        target_intervals = target_intervals,
        description = description,
        umi_paired = umi_paired
    }
}
