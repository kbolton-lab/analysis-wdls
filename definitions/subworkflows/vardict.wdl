version 1.0

import "../subworkflows/fp_filter.wdl" as ff
import "../tools/index_vcf.wdl" as iv
import "../tools/vardict.wdl" as v
import "../tools/intervals_to_bed.wdl" as itb
import "../tools/bcftools_filter_bcbio.wdl" as bfb

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

    File interval_list
    Int scatter_count
    Float? min_var_freq = 0.05

    String bcbio_filter_string = "((FMT/AF * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 2.0) || (FMT/DP < 10) || (INFO/QUAL < 45)))"
    Boolean? tumor_only = false
  }

    call itb.intervalsToBed as interval_list_to_bed {
        input:
            interval_list = interval_list
    }

    call v.vardict as vardict {
        input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = tumor_bam,
            tumor_bam_bai = tumor_bam_bai,
            normal_bam = normal_bam,
            normal_bam_bai = normal_bam_bai,
            interval_bed = interval_list_to_bed.interval_bed,
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

  call bfb.bcftoolsFilterBcbio as bcbio_filter {
      input:
        vcf = indexVcf.indexed_vcf,
        vcf = indexVcf.indexed_vcf_tbi,
        filter_string = bcbio_filter_string,
        filter_flag = "include",
        output_type = "z",
        output_vcf_name = "bcbiofilter.vcf.gz"
  }

  call iv.indexVcf as index_bcbio {
      input:
        vcf = bcbio_filter.filtered_vcf
  }

  output {
    File unfiltered_vcf = fpFilter.unfiltered_vcf
    File unfiltered_vcf_tbi = fpFilter.unfiltered_vcf_tbi
    File filtered_vcf = fpFilter.filtered_vcf
    File filtered_vcf_tbi = fpFilter.filtered_vcf_tbi
    File bcbio_filtered_vcf = index_bcbio.indexed_vcf
    File bcbio_filtered_vcf_tbi = index_bcbio.indexed_vcf_tbi
  }
}
