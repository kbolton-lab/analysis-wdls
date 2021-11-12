version 1.0

import "../tools/lofreq_pass.wdl" as lp
#import "../tools/lofreq_call.wdl" as lc
import "../tools/lofreq_reformat.wdl" as lr

workflow lofreqNoFp {
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

    output {
      File unfiltered_vcf = reformat.reformat_vcf
      File unfiltered_vcf_tbi = reformat.reformat_vcf_tbi
    }
}
