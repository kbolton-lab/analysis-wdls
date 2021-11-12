version 1.0

task filterClipAndCollectMetrics {
    input {
        File bam
        File reference
        File reference_fai
        File reference_dict
        Array[Int] min_reads = [1]
        Float? max_read_error_rate = 0.05
        Float? max_base_error_rate = 0.1
        Int min_base_quality = 1
        Float max_no_call_fraction = 0.5
        File? target_intervals
        String description
        Boolean umi_paired = true
    }

    Int cores = 1
    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: "6GB"
        cpu: cores
        bootDiskSizeGb: 5 + round(3*data_size + reference_size)
        disks: "local-disk ~{5 + round(3*data_size + reference_size)} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -eo pipefail

        /usr/local/bin/fgbio FilterConsensusReads --input ~{bam} --output consensus_filtered.bam --ref ~{reference} --min-reads ~{sep=" " min_reads} --max-read-error-rate ~{max_read_error_rate} --max-base-error-rate ~{max_base_error_rate} --min-base-quality ~{min_base_quality} --max-no-call-fraction ~{max_no_call_fraction}
        /usr/local/bin/fgbio ClipBam --input consensus_filtered.bam --ref ~{reference} --clipping-mode Hard --clip-overlapping-reads true --output clipped.bam

        PAIRED=~{umi_paired}
        DESCRIPTION=~{description}
        INTERVALS=~{target_intervals}

        if [ "$PAIRED" == true ]; then
            if [[ -z "$DESCRIPTION" ]]; then
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --output duplex_seq.metrics --intervals ${INTERVALS}
                fi
            else
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --description ${DESCRIPTION} --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --description ${DESCRIPTION} --output duplex_seq.metrics --intervals ${INTERVALS}
                fi
            fi
        else
            echo "Sample not UMI Paired" > duplex_seq.metrics.txt
        fi
    >>>

    output {
        File clipped_bam = "clipped.bam"
        File clipped_bam_bai = "clipped.bai"
        Array[File] duplex_seq_metrics = glob("duplex_seq.metrics.*")
    }
}

workflow wf {
    input {
        File bam
        File reference
        File reference_fai
        File reference_dict
        Array[Int] min_reads = [1]
        Float? max_read_error_rate = 0.05
        Float? max_base_error_rate = 0.1
        Int min_base_quality = 1
        Float max_no_call_fraction = 0.5
        File? target_intervals
        String description
        Boolean umi_paired = true

    }

    call filterClipAndCollectMetrics {
        input:
        bam = bam,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        min_reads = min_reads,
        max_read_error_rate = max_read_error_rate,
        max_base_error_rate = max_base_error_rate,
        min_base_quality = min_base_quality,
        max_no_call_fraction = max_no_call_fraction,
        target_intervals = target_intervals,
        description = description,
        umi_paired = umi_paired
    }
}
