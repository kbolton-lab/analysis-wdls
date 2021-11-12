version 1.0

task filterConsensus {
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
        bootDiskSizeGb: 10 + round(5*data_size + reference_size)
        disks: "local-disk ~{10 + round(5*data_size + reference_size)} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/local/bin/fgbio FilterConsensusReads --input ~{bam} --output consensus_filtered.bam --ref ~{reference} --min-reads ~{sep=" " min_reads} --max-read-error-rate ~{max_read_error_rate} --max-base-error-rate ~{max_base_error_rate} --min-base-quality ~{min_base_quality} --max-no-call-fraction ~{max_no_call_fraction}
    >>>

    output {
        File filtered_bam = "consensus_filtered.bam"
        File filtered_bam_bai = "consensus_filtered.bai"
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
    }

    call filterConsensus {
        input:
        bam = bam,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        min_reads = min_reads,
        max_read_error_rate = max_read_error_rate,
        max_base_error_rate = max_base_error_rate,
        min_base_quality = min_base_quality,
        max_no_call_fraction = max_no_call_fraction

    }
}
