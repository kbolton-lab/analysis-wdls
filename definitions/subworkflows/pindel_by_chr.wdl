version 1.0

import "../subworkflows/fp_filter.wdl" as ff
import "../subworkflows/pindel_cat.wdl" as pc
import "../tools/bgzip.wdl" as b
import "../tools/cat_all.wdl" as ca
import "../tools/index_vcf.wdl" as iv
import "../tools/pindel_somatic_filter.wdl" as psf
import "../tools/pindel_to_vcf.wdl" as ptv
import "../tools/remove_end_tags.wdl" as ret
import "../tools/split_interval_list_to_bed.wdl" as siltb
import "../tools/split_bed_to_chr.wdl" as siltc

workflow pindel {
  input {
    File reference
    File reference_fai
    File reference_dict
    File tumor_bam
    File tumor_bam_bai
    File? normal_bam
    File? normal_bam_bai
    File interval_bed
    String tumor_sample_name
    String? normal_sample_name
    Int insert_size = 400
    Int scatter_count = 50

    Boolean tumor_only = false
    String? ref_name = "GRCh38DH"
    String? ref_date = "20161216"
    Int? pindel_min_supporting_reads = 3
    Float? min_var_freq = 0.05
  }

  call siltc.splitBedToChr {
    input:
    interval_bed=interval_bed
  }

  scatter(region_file in splitBedToChr.split_chr) {
    call pc.pindelCat {
      input:
      reference=reference,
      reference_fai=reference_fai,
      reference_dict=reference_dict,
      tumor_bam=tumor_bam,
      tumor_bam_bai=tumor_bam_bai,
      normal_bam=normal_bam,
      normal_bam_bai=normal_bam_bai,
      region_file=region_file,
      insert_size=insert_size,
      tumor_sample_name=tumor_sample_name,
      normal_sample_name=normal_sample_name,
      tumor_only = tumor_only
    }
  }

  call ca.catAll {
    input: region_pindel_outs=pindelCat.per_region_pindel_out
  }

  if (tumor_only) {
      call ptv.pindelToVcf as pindel2vcf {
          input:
          reference=reference,
          reference_fai=reference_fai,
          reference_dict=reference_dict,
          pindel_output_summary=catAll.all_region_pindel_head,
          ref_name = ref_name,
          ref_date = ref_date,
          min_supporting_reads = pindel_min_supporting_reads
      }
  }

  if (!tumor_only) {
      call psf.pindelSomaticFilter as somaticFilter {
        input:
        reference=reference,
        reference_fai=reference_fai,
        reference_dict=reference_dict,
        pindel_output_summary=catAll.all_region_pindel_head
      }
  }

  call b.bgzip {
    input: file=select_first([somaticFilter.vcf, pindel2vcf.pindel_vcf])
  }

  call iv.indexVcf as index {
    input: vcf=bgzip.bgzipped_file
  }

  call ret.removeEndTags {
    input:
    vcf=index.indexed_vcf,
    vcf_tbi=index.indexed_vcf_tbi
  }

  call iv.indexVcf as reindex {
    input: vcf=removeEndTags.processed_vcf
  }

  call ff.fpFilter as filter {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=tumor_bam,
    bam_bai=tumor_bam_bai,
    vcf=reindex.indexed_vcf,
    vcf_tbi=reindex.indexed_vcf_tbi,
    variant_caller="pindel",
    sample_name=tumor_sample_name,
    min_var_freq = min_var_freq
  }

  output {
    File unfiltered_vcf = filter.unfiltered_vcf
    File unfiltered_vcf_tbi = filter.unfiltered_vcf_tbi
    File filtered_vcf = filter.filtered_vcf
    File filtered_vcf_tbi = filter.filtered_vcf_tbi
  }
}
