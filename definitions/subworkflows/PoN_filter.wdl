version 1.0

import "../types.wdl"

import "../tools/msk_get_base_counts.wdl" as mgbc
import "../tools/normal_fisher.wdl" as nf
import "../tools/index_vcf.wdl" as iv
import "../tools/bcftools_merge.wdl" as bm

workflow PoNFilter {
    input {
        File reference
        File reference_fai
        File reference_dict
        File caller_vcf
        File normal_bams_file
        Array[bam_and_bai] pon_bams
        Int? mapq = 5
        Int? baseq = 5
        String? pon_final_name = "pon.pileup"
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

        call bm.bcftoolsMerge as merge {
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

    output {
        File pon_total_counts = select_first([merge.merged_vcf, get_pileup_counts.pileup])
        File pon_total_counts_tbi = select_first([merge.merged_vcf_tbi, get_pileup_counts.pileup_tbi])
    }
}
