version 1.0

import "../tools/fp_filter.wdl" as ff
import "../tools/select_variants.wdl" as sv

workflow fpFilterNoNorm {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    File vcf
    File vcf_tbi
    String variant_caller
    String? sample_name
    Float? min_var_freq
  }

  call ff.fpFilter as firstFilter {
    input:
    reference=reference,
    reference_fai=reference_fai,
    bam=bam,
    vcf=vcf,
    sample_name=sample_name,
    min_var_freq=min_var_freq,
    output_vcf_basename = variant_caller + "_full"
  }

  call sv.selectVariants as hardFilter {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=firstFilter.filtered_vcf,
    vcf_tbi=firstFilter.filtered_vcf_tbi,
    exclude_filtered=true,
    output_vcf_basename = variant_caller + "_filtered"
  }

  output {
    File unfiltered_vcf = firstFilter.filtered_vcf
    File unfiltered_vcf_tbi = firstFilter.filtered_vcf_tbi
    File filtered_vcf = hardFilter.filtered_vcf
    File filtered_vcf_tbi = hardFilter.filtered_vcf_tbi
  }
}
