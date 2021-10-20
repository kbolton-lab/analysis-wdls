version 1.0

import "../tools/bcftools_isec_complement.wdl" as bic
import "../tools/msk_get_base_counts.wdl" as mgbc
import "../tools/normal_fisher.wdl" as nf
import "../tools/index_vcf.wdl" as iv

workflow gnomadAndPoNFilter {
    input {
        File reference
        File reference_fai
        File reference_dict
        File caller_vcf
        File caller_vcf_tbi
        String? caller_prefix = "caller"
        File gnomAD_exclude_vcf
        File gnomAD_exclude_vcf_tbi
        File normal_bams
        Int? mapq = 5
        Int? baseq = 5
        String? pon_final_name = "pon.pileup"
        String? pon_pvalue = "0.05"
        Boolean? arrayMode = false
    }

    call bic.bcftoolsIsecComplement as isec_complement_gnomAD {
        input:
        vcf = caller_vcf,
        vcf_tbi = caller_vcf_tbi,
        exclude_vcf = gnomAD_exclude_vcf,
        exclude_vcf_tbi = gnomAD_exclude_vcf_tbi,
        output_vcf_name = caller_prefix + ".gnomAD_AF_filter.vcf",
        output_type = "v"
    }

    call mgbc.mskGetBaseCounts as get_pileup_counts {
        input:
            normal_bams = normal_bams,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            pon_final_name = pon_final_name,
            vcf = caller_vcf,
            baseq = baseq,
            mapq = mapq,
            arrayMode = arrayMode
    }

    call nf.normalFisher as call_R_fisher {
        input:
        vcf = isec_complement_gnomAD.complement_vcf,
        pon = get_pileup_counts.pileup,
        pon_tbi = get_pileup_counts.pileup_tbi,
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
        File processed_gnomAD_filtered_vcf = index_pon_vcf.indexed_vcf
        File processed_gnomAD_filtered_vcf_tbi = index_pon_vcf.indexed_vcf_tbi
        File processed_filtered_vcf = index_pon_filtered_vcf.indexed_vcf
        File processed_filtered_vcf_tbi = index_pon_filtered_vcf.indexed_vcf_tbi
    }
}
