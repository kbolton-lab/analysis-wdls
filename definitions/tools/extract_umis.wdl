version 1.0

task extractUmis {
    input {
        File bam
        Array[String] read_structure
        Boolean? umi_paired = true
    }

    Int cores = 1
    Int space_needed_gb = 10 + round(2*size(bam, "GB"))
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
        BAM="~{bam}"
        READ_STRUCUTRE="~{sep=" " read_structure}"

        if [ "$PAIRED" == true ]; then
            /usr/local/bin/fgbio ExtractUmisFromBam --molecular-index-tags ZA ZB --single-tag RX --input $BAM --read-structure $READ_STRUCUTRE --output umi_extracted.bam
        else
            /usr/local/bin/fgbio ExtractUmisFromBam --molecular-index-tags ZA --single-tag RX --input $BAM --read-structure $READ_STRUCUTRE --output umi_extracted.bam
        fi
    >>>

    output {
        File umi_extracted_bam = "umi_extracted.bam"
    }
}

workflow wf {
    input {
        File bam
        Array[String] read_structure
        Boolean? umi_paired = true
    }
    call extractUmis {
        input:
        bam = bam,
        read_structure = read_structure,
        umi_paired = umi_paired
    }
}
