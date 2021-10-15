version 1.0

import "umi_alignment.wdl" as ua

import "../tools/merge_bams.wdl" as mb
import "../tools/group_reads.wdl" as gr     #test
import "../tools/call_molecular_consensus.wdl" as cmc
import "../tools/realign.wdl" as r
import "../tools/filter_consensus.wdl" as fc
import "../tools/clip_overlap.wdl" as co
import "../tools/duplex_seq_metrics.wdl" as dsm  # test
import "../tools/index_bam.wdl" as ib

workflow molecularAlignment {
    input {
        Array[File] bam
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
        Array[Int] min_reads = [1]
        Float? max_read_error_rate = 0.05
        Float? max_base_error_rate = 0.1
        Int min_base_quality = 1
        Float max_no_call_fraction = 0.5
    }

    scatter(bam_data in bam) {
        call ua.umiAlignment as align {
            input:
            bam = bam_data,
            read_structure = read_structure,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            reference_amb = reference_amb,
            reference_ann = reference_ann,
            reference_bwt = reference_bwt,
            reference_pac = reference_pac,
            reference_sa = reference_sa,
            umi_paired = umi_paired
        }
    }

    call mb.mergeBams as merge {
        input:
        bams = align.aligned_bam
    }

    call gr.groupReads as group_reads_by_umi {
        input:
        bam = merge.merged_bam,
        umi_paired = umi_paired
    }

    call cmc.callMolecularConsensus as call_molecular_consensus {
        input:
        bam = group_reads_by_umi.grouped_bam
    }

    call r.realign as align_consensus {
        input:
        bam = call_molecular_consensus.consensus_bam,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        reference_amb = reference_amb,
        reference_ann = reference_ann,
        reference_bwt = reference_bwt,
        reference_pac = reference_pac,
        reference_sa = reference_sa
    }

    call fc.filterConsensus as filter_consensus {
        input:
        bam = align_consensus.consensus_aligned_bam,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        min_reads = min_reads,
        max_read_error_rate = max_read_error_rate,
        max_base_error_rate = max_base_error_rate,
        min_base_quality = min_base_quality,
        max_no_call_fraction = max_no_call_fraction
    }

    call co.clipOverlap as clip_overlap {
        input:
        bam = filter_consensus.filtered_bam,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict
    }

    call dsm.duplexSeqMetrics as collect_duplex_seq_metrics {
        input:
        bam = group_reads_by_umi.grouped_bam,
        target_intervals = target_intervals,
        description = sample_name,
        umi_paired = umi_paired
    }

    call ib.indexBam as index_bam {
        input:
        bam = clip_overlap.clipped_bam
    }

    output {
        File aligned_bam = index_bam.indexed_bam
        File aligned_bam_bai = index_bam.indexed_bam_bai
        Array[File] adapter_histogram = align.adapter_metrics
        Array[File] duplex_seq_metrics = collect_duplex_seq_metrics.duplex_seq_metrics
    }
}
