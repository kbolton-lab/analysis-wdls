version 1.0

import "../types.wdl"

import "../tools/filter_umi_length.wdl" as ful
import "../tools/bbmap_repair.wdl" as br

workflow archerFastqFormat {
    input {
        SequenceData sequence
        Int umi_length=8
    }

    call ful.filterUmiLength as filterUMI {
        input:
        fastq1 = sequence.sequence.fastq1,
        fastq2 = sequence.sequence.fastq2,
        umi_length = umi_length
    }

    call br.bbmapRepair as repair {
        input:
        fastq1 = filterUMI.fastq1_filtered,
        fastq2 = filterUMI.fastq2_filtered
    }

    output {
        Array[File] fastqs = repair.fastqs
        File fastq1 = repair.fastq1_repair
        File fastq2 = repair.fastq2_repair
    }
}
