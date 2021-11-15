version 1.0

import "../tools/lofreq_pass.wdl" as lp
#import "../tools/lofreq_call.wdl" as lc
import "../tools/lofreq_reformat.wdl" as lr
import "../tools/bcftools_norm.wdl" as bn
import "../tools/vcf_sanitize.wdl" as vs
import "../tools/vt_decompose.wdl" as vd

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

    call vs.vcfSanitize as sanitizeVcf {
      input: vcf=reformat.reformat_vcf
    }

    call bn.bcftoolsNorm as normalize {
        input:
        reference=reference,
        reference_fai=reference_fai,
        vcf=sanitizeVcf.sanitized_vcf,
        vcf_tbi=sanitizeVcf.sanitized_vcf_tbi
    }

    call vd.vtDecompose as decomposeVariants {
      input:
      vcf=normalize.normalized_vcf,
      vcf_tbi=normalize.normalized_vcf_tbi
    }

    output {
      File unfiltered_vcf = decomposeVariants.decomposed_vcf
      File unfiltered_vcf_tbi = decomposeVariants.decomposed_vcf_tbi
    }
}
