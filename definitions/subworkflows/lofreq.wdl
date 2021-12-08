version 1.0

import "../subworkflows/fp_filter.wdl" as ff
import "../tools/index_vcf.wdl" as iv
import "../tools/lofreq_pass.wdl" as lp
#import "../tools/lofreq_call.wdl" as lc
import "../tools/lofreq_reformat.wdl" as lr

workflow lofreq {
    input {
        File reference
        File reference_fai
        File reference_dict
        File tumor_bam
        File tumor_bam_bai
        File? normal_bam
        File? normal_bam_bai
        File interval_bed
        String tumor_sample_name
        Float? min_var_freq = 0.05
        Boolean? tumor_only = false
    }

    call lp.lofreqPass as lofreqTask {
        input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = tumor_bam,
            tumor_bam_bai = tumor_bam_bai,
            normal_bam = normal_bam,
            normal_bam_bai = normal_bam_bai,
            interval_bed = interval_bed,
            tumor_only = tumor_only
    }

    call lr.lofreqReformat as reformat {
        input:
            vcf = lofreqTask.vcf,
            tumor_sample_name = tumor_sample_name
    }

    call ff.fpFilter {
      input:
          reference = reference,
          reference_fai = reference_fai,
          reference_dict = reference_dict,
          bam = tumor_bam,
          bam_bai = tumor_bam_bai,
          vcf = reformat.reformat_vcf,
          vcf_tbi = reformat.reformat_vcf_tbi,
          variant_caller = "lofreq",
          sample_name = tumor_sample_name,
          min_var_freq = min_var_freq
    }

    output {
      File unfiltered_vcf = fpFilter.unfiltered_vcf
      File unfiltered_vcf_tbi = fpFilter.unfiltered_vcf_tbi
      File filtered_vcf = fpFilter.filtered_vcf
      File filtered_vcf_tbi = fpFilter.filtered_vcf_tbi
    }
}
