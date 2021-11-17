version 1.0

task filterUmiLength {
    input {
        File? fastq1
        File? fastq2
        Int umi_length
    }

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([fastq1, fastq2], "GB")
    Int space_needed_gb = 10 + round(2*data_size)

    runtime {
        docker: "ubuntu:xenial"
        memory: "6GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -o pipefail
        set -o errexit
        set -o nounset

        FASTQ_ONE="~{fastq1}"
        FASTQ_TWO="~{fastq2}"
        UMI_LENGTH="~{umi_length}"

        zcat $FASTQ_ONE | awk -v regex="AACCGCCAGGAGT" -v umi_length="$UMI_LENGTH" 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; split(seq,a,regex); if (length(a[1]) == umi_length) {print header, seq, qheader, qseq}}' > R1_filtered.fastq
        gzip R1_filtered.fastq
        cp $FASTQ_TWO R2_filtered.fastq.gz
    >>>

    output {
        File fastq1_filtered = "R1_filtered.fastq.gz"
        File fastq2_filtered = "R2_filtered.fastq.gz"
    }
}

workflow wf {
    input {
        File? fastq1
        File? fastq2
        Int umi_length
    }
    call filterUmiLength {
        input:
        fastq1 = fastq1,
        fastq2 = fastq2,
        umi_length = umi_length
    }
}
