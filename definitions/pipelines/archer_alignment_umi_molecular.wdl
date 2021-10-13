version 1.0

import "../types.wdl"

import "../subworkflows/archer_fastq_format.wdl" as fqf
import "../subworkflows/molecular_alignment.wdl" as ma

import "../tools/fastq_to_bam.wdl" as ftb
import "../tools/bam_to_cram.wdl" as btc
import "../tools/index_cram.wdl" as ic

workflow archer_alignment_umi_molecular {
    input {
        Array[SequenceData] sequence
        String sample_name
        Array[String] read_structure
        File reference
        File reference_fai
        File reference_dict
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_pac
        File reference_sa
        File? target_intervals
        Boolean? umi_paired = true
        Int? umi_length = 8
        Array[Int] min_reads = [1]
        Float? max_read_error_rate = 0.05
        Float? max_base_error_rate = 0.1
        Int min_base_quality = 1
        Float max_no_call_fraction = 0.5
    }

    scatter(seq_data in sequence) {
        call fqf.archer_fastq_format as format_fastq {
            input:
            sequence = seq_data,
            umi_length = umi_length
        }
    }

    call ftb.fastq_to_bam as fastq_to_bam {
        input:
        read1_fastq = format_fastq.fastq1,
        read2_fastq = format_fastq.fastq2,
        sample_name = sample_name,
        library_name = "Library",
        platform_unit = "Illumina",
        platform = "ArcherDX"
    }

    call ma.molecular_alignment as alignment_workflow {
        input:
        bam = fastq_to_bam.bam,
        sample_name = sample_name,
        read_structure = read_structure,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        reference_amb = reference_amb,
        reference_ann = reference_ann,
        reference_bwt = reference_bwt,
        reference_pac = reference_pac,
        reference_sa = reference_sa,
        target_intervals = target_intervals,
        umi_paired = umi_paired,
        min_reads = min_reads,
        max_read_error_rate = max_read_error_rate,
        max_base_error_rate = max_base_error_rate,
        min_base_quality = min_base_quality,
        max_no_call_fraction = max_no_call_fraction
    }

    call btc.bam_to_cram as bam_to_cram {
        input:
        bam = alignment_workflow.aligned_bam,
        bam_bai = alignment_workflow.aligned_bam_bai,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
    }

    call ic.index_cram as index_cram {
        input:
        cram = bam_to_cram.cram
    }

    output {
        File aligned_cram = index_cram.indexed_cram
        File aligned_cram_crai = index_cram_crai.indexed_cram_crai
        Array[File] adapter_histogram = alignment_workflow.adapter_histogram
        Array[File] duplex_seq_metrics = alignment_workflow.duplex_seq_metrics
    }
}
