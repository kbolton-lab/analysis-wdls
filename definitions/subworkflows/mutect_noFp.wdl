version 1.0

import "../tools/merge_vcf.wdl" as mv
import "../tools/mutect.wdl" as m
import "../tools/split_interval_list.wdl" as sil

workflow mutectNoFp {
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
    Boolean no_scatter = true    # Recommended no scatter for Orientation to work properly
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

  output {
    File unfiltered_vcf = select_first([mergeVcf.merged_vcf, mutectTaskNoScatter.vcf])
    File unfiltered_vcf_tbi = select_first([mergeVcf.merged_vcf_tbi, mutectTaskNoScatter.vcf_tbi])
  }
}
