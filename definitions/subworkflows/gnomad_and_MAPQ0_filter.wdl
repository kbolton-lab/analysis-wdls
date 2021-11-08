version 1.0

import "../tools/bcftools_isec_complement.wdl" as bic
import "../tools/mapq0.wdl" as mq

workflow gnomadAndMapQ0Filter {
    input {
        File caller_vcf
        File caller_vcf_tbi
        File bam
        File bam_bai
        String? caller_prefix = "caller"
        File gnomAD_exclude_vcf
        File gnomAD_exclude_vcf_tbi
        Float mapq0perc = 0.15
        
    }

    call bic.bcftoolsIsecComplement as gnomad {
        input:
            vcf = caller_vcf,
            vcf_tbi = caller_vcf_tbi,
            exclude_vcf = gnomAD_exclude_vcf,
            exclude_vcf_tbi = gnomAD_exclude_vcf_tbi,
            output_vcf_name = caller_prefix + ".gnomAD_AF_filter.vcf.gz",
            output_type = "z"
    }

    call mq.mapq0 as MQ0 {
        input:
            vcf = gnomad.complement_vcf,
            vcf_tbi = gnomad.complement_vcf,
            bam = bam,
            bam_bai = bam_bai,
            mapq0perc = mapq0perc,
            caller = caller_prefix,
            output_type = "z"
    }

    

    output {
        File mapq0_soft_filtered_vcf = MQ0.filtered_vcf
        File mapq0_soft_filtered_vcf_tbi = MQ0.filtered_vcf_tbi
    }
}
