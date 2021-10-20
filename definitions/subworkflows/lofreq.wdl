version 1.0

import "../subworkflows/fp_filter.wdl" as ff
import "../tools/index_vcf.wdl" as iv
import "../tools/merge_vcf.wdl" as mv
import "../tools/lofreq.wdl" as l
import "../tools/lofreq_reformat.wdl" as lr
import "../tools/split_interval_list.wdl" as sil

workflow lofreq {
    input {
        File reference
        File reference_fai
        File reference_dict
        File tumor_bam
        File tumor_bam_bai
        File? normal_bam
        File? normal_bam_bai
        File interval_list
        Int scatter_count
        String tumor_sample_name
        Float? min_var_freq = 0.05
        Boolean? tumor_only = false
    }

    call sil.splitIntervalList {
      input:
      interval_list = interval_list,
      scatter_count = scatter_count
    }

    scatter(segment in splitIntervalList.split_interval_lists) {
      call l.lofreq as lofreqTask {
        input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = tumor_bam,
            tumor_bam_bai = tumor_bam_bai,
            normal_bam = normal_bam,
            normal_bam_bai = normal_bam_bai,
            interval_list = segment,
            tumor_only = tumor_only
      }
    }

    call mv.mergeVcf {
      input:
          vcfs = lofreqTask.vcf,
          vcf_tbis = lofreqTask.vcf_tbi
    }

    call lr.lofreqReformat as reformat {
        input:
            vcf = mergeVcf.merged_vcf,
            tumor_sample_name = tumor_sample_name
    }

    call iv.indexVcf {
      input:
        vcf = reformat.reformat_vcf
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
