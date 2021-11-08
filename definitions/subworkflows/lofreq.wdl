version 1.0

import "../subworkflows/fp_filter.wdl" as ff
import "../tools/index_vcf.wdl" as iv
import "../tools/merge_vcf.wdl" as mv
import "../tools/lofreq_pass.wdl" as lp
import "../tools/lofreq_call.wdl" as lc
import "../tools/lofreq_intersect.wdl" as li
import "../tools/lofreq_reformat.wdl" as lr
import "../tools/split_interval_list_to_bed.wdl" as siltb

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

    call siltb.splitIntervalListToBed {
      input:
      interval_list = interval_list,
      scatter_count = scatter_count
    }

    scatter(segment in splitIntervalListToBed.split_beds) {
      call lp.lofreqPass as lofreqPassTask {
        input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = tumor_bam,
            tumor_bam_bai = tumor_bam_bai,
            normal_bam = normal_bam,
            normal_bam_bai = normal_bam_bai,
            interval_bed = segment,
            tumor_only = tumor_only
      }

      call lc.lofreqCall as lofreqCallTask {
        input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = tumor_bam,
            tumor_bam_bai = tumor_bam_bai,
            normal_bam = normal_bam,
            normal_bam_bai = normal_bam_bai,
            interval_bed = segment,
            min_vaf = min_var_freq,
            tumor_only = tumor_only
      }
    }

    call mv.mergeVcf as mergePass {
      input:
          vcfs = lofreqPassTask.vcf,
          vcf_tbis = lofreqPassTask.vcf_tbi
    }

    call mv.mergeVcf as mergeCalls {
      input:
          vcfs = lofreqCallTask.vcf,
          vcf_tbis = lofreqCallTask.vcf_tbi
    }

    call li.lofreqIntersect as lofreqFinal {
      input:
        call_vcf = mergeCalls.merged_vcf,
        call_vcf_tbi = mergeCalls.merged_vcf_tbi,
        pass_vcf = mergePass.merged_vcf,
        pass_vcf_tbi = mergePass.merged_vcf_tbi
    }

    call lr.lofreqReformat as reformat {
        input:
            vcf = lofreqFinal.intersect_vcf,
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
