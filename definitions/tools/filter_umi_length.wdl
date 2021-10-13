version 1.0

task filter_umi_length {
    input {
        File fastq1
        File fastq2
        Int umi_length
    }

    Int cores = 8
    Float data_size = size([fastq1, fastq2], "GB")

    runtime {
        docker: "ubuntu:xenial"
        memory: "4GB"
        cpu: cores
        bootDiskSizeGb: 10 + round(5*data_size)
        disks: "local-disk ~{10 + round(5*data_size)} SSD"
    }

    String output_fastq_name = "R1_filtered"
    command <<<
        set -o pipefail
        set -o errexit
        set -o nounset

        FASTQ_ONE="~{fastq1}"
        FASTQ_TWO="~{fastq2}"
        UMI_LENGTH="~{umi_length}"

        zcat $FASTQ_ONE | awk -v regex="AACCGCCAGGAGT" -v umi_length="$UMI_LENGTH" 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; split(seq,a,regex); if (length(a[1]) == umi_length) {print header, seq, qheader, qseq}}' > ~{output_fastq_name}.fastq
        gzip ~{output_fastq_name}.fastq
        cp $FASTQ_TWO R2_filtered.fastq.gz
    >>>

    output {
        File fastq1 = "~{output_fastq_name}.fastq.gz"
        File fastq2 = "R2_filtered.fastq.gz"
    }
}

workflow wf {
    input {
        File fastq1
        File fastq2
        Int umi_length
    }
    call filter_umi_length {
        input:
        fastq1 = fastq1,
        fastq2 = fastq2,
        umi_length = umi_length
    }
}
