version 1.0

import "../subworkflows/fp_filter.wdl" as ff
import "../tools/index_vcf.wdl" as iv
import "../tools/merge_vcf.wdl" as mv
import "../tools/mutect.wdl" as m
import "../tools/split_interval_list.wdl" as sil

workflow mutect {
  input {
    File reference
    File reference_fai
    File reference_dict

    File tumor_bam
    File tumor_bam_bai
    File? pon
    File? pon_tbi
    File? gnomad_file
    File? gnomad_file_tbi

    # both or neither
    Boolean? tumor_only = false
    File? normal_bam
    File? normal_bam_bai

    File interval_list
    Int scatter_count
    Boolean no_scatter = true
    String tumor_sample_name
    Float? min_var_freq = 0.05

    String variant_caller = "mutect"
  }

  if (no_scatter) {
      call m.mutect as mutectTaskNoScatter {
        input:
        reference=reference,
        reference_fai=reference_fai,
        reference_dict=reference_dict,
        tumor_bam=tumor_bam,
        tumor_bam_bai=tumor_bam_bai,
        normal_bam=normal_bam,
        normal_bam_bai=normal_bam_bai,
        interval_list=interval_list,
        tumor_only=tumor_only,
        pon = pon,
        pon_tbi = pon_tbi,
        gnomad = gnomad_file,
        gnomad_tbi = gnomad_file_tbi,
      }
  }

  if (!no_scatter) {
      call sil.splitIntervalList {
        input:
        interval_list=interval_list,
        scatter_count=scatter_count
      }

      scatter(segment in splitIntervalList.split_interval_lists) {
        call m.mutect as mutectTask {
          input:
          reference=reference,
          reference_fai=reference_fai,
          reference_dict=reference_dict,
          tumor_bam=tumor_bam,
          tumor_bam_bai=tumor_bam_bai,
          normal_bam=normal_bam,
          normal_bam_bai=normal_bam_bai,
          interval_list=segment,
          tumor_only=tumor_only
        }
      }

      call mv.mergeVcf {
        input:
        vcfs=mutectTask.vcf,
        vcf_tbis=mutectTask.vcf_tbi
      }
  }

  call ff.fpFilter {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=tumor_bam,
    bam_bai=tumor_bam_bai,
    vcf=select_first([mergeVcf.merged_vcf, mutectTaskNoScatter.vcf]),
    vcf_tbi=select_first([mergeVcf.merged_vcf_tbi, mutectTaskNoScatter.vcf_tbi]),
    variant_caller=variant_caller,
    sample_name=tumor_sample_name,
    min_var_freq=min_var_freq
  }

  output {
    File unfiltered_vcf = fpFilter.unfiltered_vcf
    File unfiltered_vcf_tbi = fpFilter.unfiltered_vcf_tbi
    File filtered_vcf = fpFilter.filtered_vcf
    File filtered_vcf_tbi = fpFilter.filtered_vcf_tbi
  }
}
