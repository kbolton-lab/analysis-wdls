version 1.0

import "../subworkflows/fisher_filter.wdl" as ff
import "../tools/annotate_vcf.wdl" as av

workflow annotateCaller {
    input {
        File vcf
        File vcf_tbi
        File fp_filter
        File fp_filter_tbi
        File pileup_file
        File pileup_file_tbi
        File vep
        File vep_tbi
        String caller_prefix
        String sample_name
        String? pon_pvalue = "0.05"
    }

    call ff.fisherFilter as fisher_test {
        input:
        caller_vcf = vcf,
        caller_prefix = caller_prefix,
        pon_pvalue = pon_pvalue,
        pileup_vcf = pileup_file,
        pileup_vcf_tbi = pileup_file_tbi
    }

    call av.annotateVcf as annotate_vcf {
        input:
        vcf = fisher_test.processed_filtered_vcf,
        vcf_tbi = fisher_test.processed_filtered_vcf_tbi,
        fp_filter = fp_filter,
        fp_filter_tbi = fp_filter_tbi,
        vep = vep,
        vep_tbi = vep_tbi,
        caller_prefix = caller_prefix,
        sample_name = sample_name
    }

    output {
        File pon_annotated_vcf = annotate_vcf.pon_annotated_vcf
        File pon_annotated_vcf_tbi = annotate_vcf.pon_annotated_vcf_tbi
        File final_annotated_vcf = annotate_vcf.final_annotated_vcf
        File final_annotated_vcf_tbi = annotate_vcf.final_annotated_vcf_tbi
    }
}
