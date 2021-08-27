version 1.0

import "../tools/gatk_haplotype_caller.wdl" as ghc

workflow gatkHaplotypecallerIterator {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bai
    String emit_reference_confidence  # enum [NONE, BP_RESOLUTION, GVCF]
    Array[String] gvcf_gq_bands
    Array[Array[String]] intervals
    File? dbsnp_vcf
    File? dbsnp_vcf_tbi
    File verify_bam_id_metrics
    Int? max_alternate_alleles
    Int? ploidy
    String? read_filter
    String? output_prefix
  }

  String pref = if defined(output_prefix) then select_first([output_prefix]) else ""
  # TODO: if only alphanumeric
  String base = if length(intervals) == 1 then intervals[0] else "output"
  String output_file_name = pref + "." + base + ".g.vcf.gz"

  scatter(interval_sublist in intervals) {
    call ghc.gatkHaplotypeCaller as haplotypeCaller {
      input:
      reference=reference,
      reference_fai=reference_fai,
      reference_dict=reference_dict,
      bam=bam,
      bai=bai,
      emit_reference_confidence=emit_reference_confidence,
      gvcf_gq_bands=gvcf_gq_bands,
      intervals=interval_sublist,
      dbsnp_vcf=dbsnp_vcf,
      dbsnp_vcf_tbi=dbsnp_vcf_tbi,
      verify_bam_id_metrics=verify_bam_id_metrics,
      max_alternate_alleles=max_alternate_alleles,
      ploidy=ploidy,
      read_filter=read_filter,
      output_prefix = output_prefix
    }
  }

  output {
    Array[File] gvcf = haplotypeCaller.gvcf
    Array[File] gvcf_tbi = haplotypeCaller.gvcf_tbi
  }
}
