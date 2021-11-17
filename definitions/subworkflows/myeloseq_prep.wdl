version 1.0

import "../types.wdl"

import "../tools/trim_reads.wdl" as tr
import "../tools/fgbio_fastq_to_bam.wdl" as fftb
import "../tools/fgbio_annotate_bam.wdl" as ab

workflow myeloseqPrep {
    input {
        SequenceData sequence
        String adapter_one
        String adapter_two
        String sample_name
        File fastq_with_umis
    }

    call tr.trimReads as trim_reads {
        input:
        fastq1 = sequence.sequence.fastq1,
        fastq2 = sequence.sequence.fastq2,
        adapter_one = adapter_one,
        adapter_two = adapter_two
    }

    call fftb.fgbioFastqToBam as fastq_to_bam {
        input:
        trimmed_fastqs = trim_reads.trimmed_fastqs,
        sample_name = sample_name,
        library_name = "Library",
        platform_unit = "Illumina",
        platform = "Illumina"
    }

    call ab.annotateBam as annotate_bam {
        input:
        bam = fastq_to_bam.bam,
        fastq_with_umis = fastq_with_umis
    }

    output {
        File bam = annotate_bam.umi_bam
    }
}
