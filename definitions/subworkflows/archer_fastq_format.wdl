version 1.0

import "../types.wdl"

import "../tools/filter_umi_length.wdl" as ful
import "../tools/bbmap_repair.wdl" as br

workflow archer_fastq_format {
    input {
        SequenceData sequence
        Int umi_length=8
    }

    call ful.filter_umi_length as filter_umi_length {
        input:
        fastq1 = sequence.sequence.fastq1,
        fastq2 = sequence.sequence.fastq2,
        umi_length = umi_length
    }

    call br.bbmap_repair as repair {
        input:
        fastq1 = filter_umi_length.fastq1_filtered,
        fastq2 = filter_umi_length.fastq2_filtered
    }

    output {
        Array[File] fastqs = repair.fastqs
        File fastq1 = repair.fastq1_repair
        File fastq2 = repair.fastq2_repair
    }
}
