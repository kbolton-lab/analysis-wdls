version 1.0

import "../subworkflows/fp_filter.wdl" as ff
import "../tools/varscan.wdl" as v
import "../tools/intervals_to_bed.wdl" as itb
import "../tools/bcftools_filter_bcbio.wdl" as bfb
import "../tools/mapq0.wdl" as mq
import "../tools/split_interval_list_to_bed.wdl" as siltb
import "../tools/merge_vcf.wdl" as ca

workflow varscan {
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

    File interval_list
    Int scatter_count = 30
    Float? min_var_freq = 0.05

    Boolean? tumor_only = false
  }

    call siltb.splitIntervalListToBed {
      input:
      interval_list=interval_list,
      scatter_count=scatter_count
    }

  
    scatter(region_file in splitIntervalListToBed.split_beds) {
      call v.varscan as varscanTask {
        input:
          reference = reference,
          reference_fai = reference_fai,
          tumor_bam = tumor_bam,
          tumor_bam_bai = tumor_bam_bai,
          normal_bam = normal_bam,
          normal_bam_bai = normal_bam_bai,
          interval_bed = region_file,
          min_var_freq = min_var_freq,
          tumor_sample_name = tumor_sample_name,
          normal_sample_name = normal_sample_name,
          tumor_only = tumor_only
      }
    }

    call ca.mergeVcf {
      input: 
        vcfs=varscanTask.vcf,
        vcf_tbis=varscanTask.vcf_tbi
    }

    call ff.fpFilter {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bam = tumor_bam,
            bam_bai = tumor_bam_bai,
            vcf = mergeVcf.merged_vcf,
            vcf_tbi = mergeVcf.merged_vcf_tbi,
            variant_caller = "varscan",
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
