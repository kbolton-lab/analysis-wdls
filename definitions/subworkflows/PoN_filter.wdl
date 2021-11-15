version 1.0

import "../types.wdl"

import "../tools/msk_get_base_counts.wdl" as mgbc
import "../tools/normal_fisher.wdl" as nf
import "../tools/index_vcf.wdl" as iv
import "../tools/merge_vcf.wdl" as mv

workflow PoNFilter {
    input {
        File reference
        File reference_fai
        File reference_dict
        File caller_vcf
        File caller_vcf_tbi
        String? caller_prefix = "caller"
        File normal_bams_file
        Array[bam_and_bai] pon_bams
        Int? mapq = 5
        Int? baseq = 5
        String? pon_final_name = "pon.pileup"
        String? pon_pvalue = "0.05"
        Boolean arrayMode = true
    }

    if (arrayMode) {
        scatter (pon_bam in pon_bams) {
            call mgbc.mskGetBaseCounts as mskGetBaseCounts {
                input:
                reference = reference,
                reference_fai = reference_fai,
                reference_dict = reference_dict,
                normal_bam = pon_bam,
                pon_final_name = pon_final_name,
                vcf = caller_vcf,
                mapq = mapq,
                baseq = baseq
            }
        }

        call mv.mergeVcf as merge {
            input:
                vcfs = mskGetBaseCounts.pileup,
                vcf_tbis = mskGetBaseCounts.pileup_tbi
        }
    }
    if (!arrayMode) {
        call mgbc.mskGetBaseCountsWithFile as get_pileup_counts {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            normal_bams = normal_bams_file,
            pon_final_name = pon_final_name,
            vcf = caller_vcf,
            mapq = mapq,
            baseq = baseq
        }
    }

    call nf.normalFisher as call_R_fisher {
        input:
        vcf = caller_vcf,
        pon = select_first([merge.merged_vcf, get_pileup_counts.pileup]),
        pon_tbi = select_first([merge.merged_vcf_tbi, get_pileup_counts.pileup_tbi]),
        p_value = pon_pvalue,
        caller = caller_prefix
    }

    call iv.indexVcf as index_pon_vcf {
        input: vcf = call_R_fisher.pon_vcf
    }

    call iv.indexVcf as index_pon_filtered_vcf {
        input: vcf = call_R_fisher.pon_filtered_vcf
    }

    output {
        File processed_filtered_vcf = index_pon_filtered_vcf.indexed_vcf
        File processed_filtered_vcf_tbi = index_pon_filtered_vcf.indexed_vcf_tbi
        File pon_total_counts = select_first([merge.merged_vcf, get_pileup_counts.pileup])
        File pon_total_counts_tbi = select_first([merge.merged_vcf_tbi, get_pileup_counts.pileup_tbi])
    }
}
