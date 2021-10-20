version 1.0

import "../subworkflows/fp_filter.wdl" as ff
import "../tools/index_vcf.wdl" as iv
import "../tools/merge_vcf.wdl" as mv
import "../tools/vardict.wdl" as v
import "../tools/split_interval_list_to_bed.wdl" as sil

workflow vardict {
  input {
    File reference
    File reference_fai
    File reference_dict

    File tumor_bam
    File tumor_bam_bai
    String tumor_sample_name
    # both or neither
    File? normal_bam
    File? normal_bam_bai
    String? normal_sample_name

    File interval_bed
    Int scatter_count
    Float? min_var_freq = 0.05

    Boolean? tumor_only = false
  }

#   call sil.splitIntervalListToBed {
#     input:
#     interval_list=interval_bed,
#     scatter_count=scatter_count
#   }

#   scatter(segment in splitIntervalListToBed.split_beds) {
#     call v.vardict as vardictTask {
#       input:
#       reference=reference,
#       reference_fai=reference_fai,
#       tumor_bam=tumor_bam,
#       tumor_bam_bai=tumor_bam_bai,
#       normal_bam=normal_bam,
#       normal_bam_bai=normal_bam_bai,
#       interval_bed=segment
#     }
#   }

#   call mv.mergeVcf {
#     input:
#     vcfs=mutectTask.vcf,
#     vcf_tbis=mutectTask.vcf_tbi
#   }

#   call iv.indexVcf {
#     input: vcf=mergeVcf.merged_vcf
#   }

    call v.vardict as vardict {
        input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = tumor_bam,
            tumor_bam_bai = tumor_bam_bai,
            normal_bam = normal_bam,
            normal_bam_bai = normal_bam_bai,
            interval_bed = interval_bed,
            min_var_freq = min_var_freq,
            tumor_sample_name = tumor_sample_name,
            normal_sample_name = normal_sample_name,
            tumor_only = tumor_only
    }

    call iv.indexVcf {
        input:
            vcf = vardict.vcf
    }

    call ff.fpFilter {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bam = tumor_bam,
            bam_bai = tumor_bam_bai,
            vcf = indexVcf.indexed_vcf,
            vcf_tbi = indexVcf.indexed_vcf_tbi,
            variant_caller = "vardict",
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
