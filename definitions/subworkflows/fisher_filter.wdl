version 1.0

import "../types.wdl"

import "../tools/normal_fisher.wdl" as nf
import "../tools/index_vcf.wdl" as iv

workflow fisherFilter {
    input {
        File caller_vcf
        String? caller_prefix = "caller"
        String? pon_pvalue = "0.05"
        File pileup_vcf
        File pileup_vcf_tbi
    }

    call nf.normalFisher as call_R_fisher {
        input:
        vcf = caller_vcf,
        pon = pileup_vcf,
        pon_tbi = pileup_vcf_tbi,
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
        File processed_unfiltered_vcf = index_pon_vcf.indexed_vcf
        File processed_unfiltered_vcf_tbi = index_pon_vcf.indexed_vcf_tbi
    }
}
