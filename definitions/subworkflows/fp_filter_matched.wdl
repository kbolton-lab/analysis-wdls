version 1.0

import "../tools/fp_filter.wdl" as ff
import "../tools/bcftools_norm.wdl" as bn
import "../tools/bcftools_extract_tumor.wdl" as e
import "../tools/vcf_sanitize.wdl" as vs
import "../tools/vt_decompose.wdl" as vd
import "../tools/select_variants.wdl" as sv

workflow fpFilter {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    File vcf
    File vcf_tbi
    String variant_caller
    String sample_name
    Float? min_var_freq
  }

  call vs.vcfSanitize as sanitizeVcf {
    input: vcf=vcf
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

  call ff.fpFilter as firstFilter {
    input:
    reference=reference,
    reference_fai=reference_fai,
    bam=bam,
    vcf=decomposeVariants.decomposed_vcf,
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

  call e.extractTumor as extractTumor {
    input:
        vcf=hardFilter.filtered_vcf,
        tumor_sample_name=sample_name,
        output_type = "z",
        output_vcf_basename = variant_caller + "_filtered_" + sample_name
  }
  

  output {
    File unfiltered_vcf = firstFilter.filtered_vcf
    File unfiltered_vcf_tbi = firstFilter.filtered_vcf_tbi
    File filtered_vcf = extractTumor.filtered_vcf
    File filtered_vcf_tbi = extractTumor.filtered_vcf_tbi
  }
}
