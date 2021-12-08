version 1.0

import "umi_alignment_no_extract.wdl" as uane

import "../tools/merge_bams.wdl" as mb
import "../tools/group_reads_and_consensus.wdl" as grac
import "../tools/realign.wdl" as r
import "../tools/filter_clip_and_collect_metrics.wdl" as fcacm
import "../tools/index_bam.wdl" as ib

workflow molecularAlignmentNoExtract {
    input {
        Array[File] bam
        String sample_name
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
        call uane.umiAlignmentNoExtract as align {
            input:
            bam = bam_data,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            reference_amb = reference_amb,
            reference_ann = reference_ann,
            reference_bwt = reference_bwt,
            reference_pac = reference_pac,
            reference_sa = reference_sa
        }
    }

    # I don't see why we need to spin up an instance to merge non-existent BAMs
    if (length(align.aligned_bam) > 1) {
        call mb.mergeBams as merge {
            input:
            bams = align.aligned_bam
        }
    }

    # Just take the BAM out of the array and move it to the next tool
    if (length(align.aligned_bam) == 1) {
        File single_bam = align.aligned_bam[0]
    }

    call grac.groupReadsAndConsensus as group_reads_by_umi_and_call_molecular_consensus {
        input:
        bam = select_first([merge.merged_bam, single_bam]),
        umi_paired = umi_paired,
    }

    call r.realign as align_consensus {
        input:
        bam = group_reads_by_umi_and_call_molecular_consensus.consensus_bam,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        reference_amb = reference_amb,
        reference_ann = reference_ann,
        reference_bwt = reference_bwt,
        reference_pac = reference_pac,
        reference_sa = reference_sa
    }

    call fcacm.filterClipAndCollectMetrics as filter_clip_collect_metrics {
        input:
        bam = align_consensus.consensus_aligned_bam,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        min_reads = min_reads,
        max_read_error_rate = max_read_error_rate,
        max_base_error_rate = max_base_error_rate,
        min_base_quality = min_base_quality,
        max_no_call_fraction = max_no_call_fraction,
        target_intervals = target_intervals,
        description = sample_name,
        umi_paired = umi_paired
    }

    call ib.indexBam as index_bam {
        input:
        bam = filter_clip_collect_metrics.clipped_bam
    }

    output {
        File aligned_bam = index_bam.indexed_bam
        File aligned_bam_bai = index_bam.indexed_bam_bai
        Array[File] adapter_histogram = align.adapter_metrics
        Array[File] duplex_seq_metrics = filter_clip_collect_metrics.duplex_seq_metrics
    }
}
