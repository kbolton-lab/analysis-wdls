version 1.0

import "../tools/vardict.wdl" as v
import "../tools/bcftools_filter_bcbio.wdl" as bfb
import "../tools/bcftools_norm.wdl" as bn
import "../tools/vcf_sanitize.wdl" as vs
import "../tools/vt_decompose.wdl" as vd

workflow vardictNoFp {
  input {
    File reference
    File reference_fai

    File tumor_bam
    File tumor_bam_bai
    String tumor_sample_name
    # both or neither
    File? normal_bam
    File? normal_bam_bai
    String? normal_sample_name

    File interval_bed
    Float? min_var_freq = 0.05

    ## does BQSR make all Vardict's FMT/NM = 0 ?????
    String bcbio_filter_string = "((FMT/AF * FMT/DP < 3) && ( FMT/MQ < 55.0 || FMT/DP < 10 || FMT/QUAL < 30 ))"
    # String bcbio_filter_string = "((FMT/AF * FMT/DP < 6) && ((FMT/MQ < 55.0 && FMT/NM > 1.0) || (FMT/MQ < 60.0 && FMT/NM > 2.0) || (FMT/DP < 10) || (FMT/QUAL < 45)))" # AF=15, MQ=24, NM=27, DP=8, QUAL=20
    Boolean? tumor_only = false
  }

  # vardict multi-threading is best
  # call itb.intervalsToBed as interval_list_to_bed {
  #     input:
  #         interval_list = interval_list
  # }

  call v.vardict as vardictTask {
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

  call bfb.bcftoolsFilterBcbio as bcbio_filter {
    input:
      vcf = vardictTask.vcf,
      vcf_tbi = vardictTask.vcf_tbi,
      filter_string = bcbio_filter_string,
      filter_flag = "include",
      output_type = "z",
      output_vcf_prefix = "vardict.bcbiofilter"
  }

  call vs.vcfSanitize as sanitizeVcf {
    input: vcf=bcbio_filter.filtered_vcf
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
    File unfiltered_vcf = vardictTask.vcf
    File unfiltered_vcf_tbi = vardictTask.vcf_tbi
    File bcbio_filtered_vcf = decomposeVariants.decomposed_vcf
    File bcbio_filtered_vcf_tbi = decomposeVariants.decomposed_vcf_tbi
  }
}
