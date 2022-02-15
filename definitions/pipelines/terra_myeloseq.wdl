version 1.0

struct TrimmingOptions {
  File adapters
  Int min_overlap
}

struct LabelledFile {
  File file
  String label
}

struct bam_and_bai {
    File bam
    File bai
}

struct bam_and_bai_array {
    Array[File] bams
    Array[File] bais
}

# ---- vep_custom_annotation ----
struct Info {
  File file
  Array[File]? secondary_files
  String data_format  # enum, ['bed', 'gff', 'gtf', 'vcf', 'bigwig']
  String name
  Array[String]? vcf_fields
  Boolean? gnomad_filter
  Boolean check_existing
}

struct VepCustomAnnotation {
  Boolean force_report_coordinates
  String method  # enum, ['exact', 'overlap']
  Info annotation
}

workflow myeloseq {
    input {
        # Pipeline
        Int scatter_count = 20

        # Sequence and BAM Information
        File bam_file
        File? bam_file_bai
        Boolean? aligned = false
        String? tumor_name = "tumor"
        String tumor_sample_name
        File target_intervals
        File target_bed

        # Reference
        File reference
        File reference_fai
        File reference_dict
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_pac
        File reference_sa

        # FASTQ Preprocessing
        Boolean? umi_paired = true
        Array[Int] min_reads = [1]
        Float? max_read_error_rate = 0.05
        Float? max_base_error_rate = 0.1
        Int min_base_quality = 1
        Float max_no_call_fraction = 0.5

        # BQSR
        Array[File] bqsr_known_sites
        Array[File] bqsr_known_sites_tbi
        Array[String] bqsr_intervals

        # QC
        File bait_intervals
        Array[LabelledFile] per_base_intervals
        Array[LabelledFile] per_target_intervals
        Array[LabelledFile] summary_intervals
        File omni_vcf
        File omni_vcf_tbi
        String picard_metric_accumulation_level
        Int? qc_minimum_mapping_quality = 0
        Int? qc_minimum_base_quality = 0

        # Variant Calling
        Boolean? arrayMode = true       # Decide if you would rather use the File (--bam_fof) or the Array (--bam) does the same thing, just input type is different
        File pon_normal_bams_file       # on GCP, it's not possible to do File because the file paths are unaccessable for each VM instance, so you have to do ArrayMode
        Array[bam_and_bai] pon_bams     # This is just an array of Files... if you have the bam_fof it's easier just to use above and set this to empty array []
        Boolean tumor_only = true
        Float? af_threshold = 0.0001
        String? pon_pvalue = "1.59442843e-7"

        # Pindel
        Int pindel_insert_size = 400
        String? ref_name = "GRCh38DH"
        String? ref_date = "20161216"
        Int? pindel_min_supporting_reads = 3

        String bcbio_filter_string = "((FMT/AF * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 3.0) || (FMT/DP < 1600) || (INFO/QUAL < 27)))"

        # PoN2
        File mutect_pon2_file
        File mutect_pon2_file_tbi
        File lofreq_pon2_file
        File lofreq_pon2_file_tbi
        File vardict_pon2_file
        File vardict_pon2_file_tbi
        #File pindel_pon2_file
        #File pindel_pon2_file_tbi

        # R Stuff
        File bolton_bick_vars
        File mut2_bick
        File mut2_kelly
        File matches2
        File TSG_file
        File oncoKB_curated
        File pd_annotation_file
        File pan_myeloid
        File truncating
        File cosmic_dir_zip

        # VEP
        File vep_cache_dir_zip
        String vep_ensembl_assembly
        String vep_ensembl_version
        String vep_ensembl_species
        Array[String] vep_plugins = ["Frameshift", "Wildtype"]
        File? synonyms_file
        Boolean? annotate_coding_only = true
        Array[VepCustomAnnotation] vep_custom_annotations
        Array[String] variants_to_table_fields = ["CHROM","POS","ID","REF","ALT","set","AC","AF"]
        Array[String]? variants_to_table_genotype_fields = ["GT","AD"]
        Array[String]? vep_to_table_fields = ["HGVSc","HGVSp"]
        String vep_pick = "pick"
        Boolean everything = true

        # gnomAD
        Float filter_gnomADe_maximum_population_allele_frequency = 0.005
        File normalized_gnomad_exclude
        File normalized_gnomad_exclude_tbi
        String filter_flag = "include"
    }

    if (!aligned) {
        call markIlluminaAdapters as mark_illumina_adapters {
            input:
            bam = bam_file
        }

        call umiAlign as align {
            input:
            bam = mark_illumina_adapters.marked_bam,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            reference_amb = reference_amb,
            reference_ann = reference_ann,
            reference_bwt = reference_bwt,
            reference_pac = reference_pac,
            reference_sa = reference_sa
        }

        call groupReadsAndConsensus as group_reads_by_umi_and_call_molecular_consensus {
            input:
            bam = align.aligned_bam,
            umi_paired = umi_paired,
        }

        call realign as align_consensus {
            input:
            bam = group_reads_by_umi_and_call_molecular_consensus.consensus_bam,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            reference_amb = reference_amb,
            reference_ann = reference_ann,
            reference_bwt = reference_bwt,
            reference_pac = reference_pac,
            reference_sa = reference_sa
        }

        call filterClipAndCollectMetrics as filter_clip_collect_metrics {
            input:
            bam = align_consensus.consensus_aligned_bam,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            min_reads = min_reads,
            max_read_error_rate = max_read_error_rate,
            max_base_error_rate = max_base_error_rate,
            min_base_quality = min_base_quality,
            max_no_call_fraction = max_no_call_fraction,
            target_intervals = target_intervals,
            description = tumor_sample_name,
            umi_paired = umi_paired
        }

        call indexBam as index_align_bam {
            input:
            bam = filter_clip_collect_metrics.clipped_bam
        }
    }

    call bqsrApply as bqsr {
        input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        bam = select_first([index_align_bam.indexed_bam, bam_file]),
        bam_bai = select_first([index_align_bam.indexed_bam_bai, bam_file_bai]),
        intervals = bqsr_intervals,
        known_sites = bqsr_known_sites,
        known_sites_tbi = bqsr_known_sites_tbi
    }

    call indexBam as index_bam {
        input:
        bam = bqsr.bqsr_bam
    }

    call Metrics as collectAllMetrics {
        input:
        bam=index_bam.indexed_bam,
        bam_bai=index_bam.indexed_bam_bai,
        reference=reference,
        reference_fai=reference_fai,
        reference_dict=reference_dict,
        metric_accumulation_level=picard_metric_accumulation_level
    }

    call collectHsMetrics as collectRoiHsMetrics {
      input:
      bam=index_bam.indexed_bam,
      bam_bai=index_bam.indexed_bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      reference_dict=reference_dict,
      metric_accumulation_level="ALL_READS",
      bait_intervals=bait_intervals,
      target_intervals=target_intervals,
      per_target_coverage = false,
      per_base_coverage = false,
      output_prefix = "roi",
      minimum_mapping_quality=qc_minimum_mapping_quality,
      minimum_base_quality=qc_minimum_base_quality
    }

    scatter(interval in summary_intervals) {
      call collectHsMetrics as collectSummaryHsMetrics{
        input:
        bam=index_bam.indexed_bam,
        bam_bai=index_bam.indexed_bam_bai,
        reference=reference,
        reference_fai=reference_fai,
        reference_dict=reference_dict,
        bait_intervals=interval.file,
        target_intervals=interval.file,
        output_prefix="summary-~{interval.label}",
        minimum_mapping_quality=qc_minimum_mapping_quality,
        minimum_base_quality=qc_minimum_base_quality,
        metric_accumulation_level="ALL_READS",
        per_target_coverage=false,
        per_base_coverage=false
      }
    }

    scatter(interval in per_base_intervals) {
      call collectHsMetrics as collectPerBaseHsMetrics {
        input:
        bam=index_bam.indexed_bam,
        bam_bai=index_bam.indexed_bam_bai,
        reference=reference,
        reference_fai=reference_fai,
        reference_dict=reference_dict,
        bait_intervals=interval.file,
        target_intervals=interval.file,
        output_prefix="base-~{interval.label}",
        minimum_mapping_quality=qc_minimum_mapping_quality,
        minimum_base_quality=qc_minimum_base_quality,
        metric_accumulation_level="ALL_READS",
        per_target_coverage=false,
        per_base_coverage=true
      }
    }

    scatter(interval in per_target_intervals) {
      call collectHsMetrics as collectPerTargetHsMetrics{
        input:
        bam=index_bam.indexed_bam,
        bam_bai=index_bam.indexed_bam_bai,
        reference=reference,
        reference_fai=reference_fai,
        reference_dict=reference_dict,
        bait_intervals=interval.file,
        target_intervals=interval.file,
        output_prefix="target-~{interval.label}",
        minimum_mapping_quality=qc_minimum_mapping_quality,
        minimum_base_quality=qc_minimum_base_quality,
        metric_accumulation_level="ALL_READS",
        per_target_coverage=true,
        per_base_coverage=false
      }
    }

    call samtoolsFlagstat {
      input:
      bam=index_bam.indexed_bam,
      bam_bai=index_bam.indexed_bam_bai
    }

    call selectVariants {
      input:
      reference=reference,
      reference_fai=reference_fai,
      reference_dict=reference_dict,
      vcf=omni_vcf,
      vcf_tbi=omni_vcf_tbi,
      interval_list=target_intervals
    }

    call verifyBamId {
      input:
      bam=index_bam.indexed_bam,
      bam_bai=index_bam.indexed_bam_bai,
      vcf=selectVariants.filtered_vcf
    }

    call splitBedToChr as split_bed_to_chr {
        input:
        interval_bed = target_bed
    }

    scatter (chr_bed in split_bed_to_chr.split_chr) {
        # Mutect
        call mutectTumorOnly as mutectTask {
          input:
          reference=reference,
          reference_fai=reference_fai,
          reference_dict=reference_dict,
          tumor_bam=index_bam.indexed_bam,
          tumor_bam_bai=index_bam.indexed_bam_bai,
          interval_list=chr_bed
        }

        call vcfSanitize as mutectSanitizeVcf {
          input: vcf=mutectTask.vcf
        }

        call bcftoolsNorm as mutectNormalize {
            input:
            reference=reference,
            reference_fai=reference_fai,
            vcf=mutectSanitizeVcf.sanitized_vcf,
            vcf_tbi=mutectSanitizeVcf.sanitized_vcf_tbi
        }

        call vtDecompose as mutectDecomposeVariants {
          input:
          vcf=mutectNormalize.normalized_vcf,
          vcf_tbi=mutectNormalize.normalized_vcf_tbi
        }
        call bcftoolsIsecComplement as mutect_isec_complement_gnomAD {
            input:
            vcf = mutectDecomposeVariants.decomposed_vcf,
            vcf_tbi = mutectDecomposeVariants.decomposed_vcf_tbi,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            output_vcf_name = "mutect." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
            output_type = "z"
        }
        call pon2Percent as mutect_pon2 {
            input:
            vcf = mutect_isec_complement_gnomAD.complement_vcf,
            #vcf = mutectDecomposeVariants.decomposed_vcf,
            vcf2PON = mutect_pon2_file,
            vcf2PON_tbi = mutect_pon2_file_tbi,
            caller = "mutect",
            sample_name = tumor_sample_name
        }

        # Vardict
        call vardictTumorOnly as vardictTask {
          input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = index_bam.indexed_bam,
            tumor_bam_bai = index_bam.indexed_bam_bai,
            interval_bed = chr_bed,
            min_var_freq = af_threshold,
            tumor_sample_name = tumor_sample_name
        }

        call bcftoolsFilterBcbio as bcbio_filter {
          input:
            vcf = vardictTask.vcf,
            vcf_tbi = vardictTask.vcf_tbi,
            filter_string = bcbio_filter_string,
            filter_flag = "exclude",
            output_type = "z",
            output_vcf_prefix = "vardict.bcbiofilter"
        }

        call vcfSanitize as vardictSanitizeVcf {
          input: vcf=bcbio_filter.filtered_vcf
        }

        call bcftoolsNorm as vardictNormalize {
            input:
            reference=reference,
            reference_fai=reference_fai,
            vcf=vardictSanitizeVcf.sanitized_vcf,
            vcf_tbi=vardictSanitizeVcf.sanitized_vcf_tbi
        }

        call vtDecompose as vardictDecomposeVariants {
          input:
          vcf=vardictNormalize.normalized_vcf,
          vcf_tbi=vardictNormalize.normalized_vcf_tbi
        }
        call bcftoolsIsecComplement as vardict_isec_complement_gnomAD {
            input:
            vcf = vardictDecomposeVariants.decomposed_vcf,
            vcf_tbi = vardictDecomposeVariants.decomposed_vcf_tbi,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            output_vcf_name = "vardict." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
            output_type = "z"
        }
        call pon2Percent as vardict_pon2 {
            input:
            vcf = vardict_isec_complement_gnomAD.complement_vcf,
            #vcf = vardictDecomposeVariants.decomposed_vcf,
            vcf2PON = vardict_pon2_file,
            vcf2PON_tbi = vardict_pon2_file_tbi,
            caller = "vardict",
            sample_name = tumor_sample_name
        }

        # Lofreq
        call lofreqTumorOnly as lofreqTask {
            input:
                reference = reference,
                reference_fai = reference_fai,
                tumor_bam = index_bam.indexed_bam,
                tumor_bam_bai = index_bam.indexed_bam_bai,
                interval_bed = chr_bed
        }

        call lofreqReformat as reformat {
            input:
                vcf = lofreqTask.vcf,
                tumor_sample_name = tumor_sample_name
        }

        call vcfSanitize as lofreqSanitizeVcf {
          input: vcf=reformat.reformat_vcf
        }

        call bcftoolsNorm as lofreqNormalize {
            input:
            reference=reference,
            reference_fai=reference_fai,
            vcf=lofreqSanitizeVcf.sanitized_vcf,
            vcf_tbi=lofreqSanitizeVcf.sanitized_vcf_tbi
        }

        call vtDecompose as lofreqDecomposeVariants {
          input:
          vcf=lofreqNormalize.normalized_vcf,
          vcf_tbi=lofreqNormalize.normalized_vcf_tbi
        }
        call bcftoolsIsecComplement as lofreq_isec_complement_gnomAD {
            input:
            vcf = lofreqDecomposeVariants.decomposed_vcf,
            vcf_tbi = lofreqDecomposeVariants.decomposed_vcf_tbi,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            output_vcf_name = "lofreq." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
            output_type = "z"
        }
        call pon2Percent as lofreq_pon2 {
            input:
            vcf = lofreq_isec_complement_gnomAD.complement_vcf,
            #vcf = lofreqDecomposeVariants.decomposed_vcf,
            vcf2PON = lofreq_pon2_file,
            vcf2PON_tbi = lofreq_pon2_file_tbi,
            caller = "lofreq",
            sample_name = tumor_sample_name
        }

        # Pindel
        call pindelTumorOnly as pindelTask {
            input:
            reference=reference,
            reference_fai=reference_fai,
            reference_dict=reference_dict,
            tumor_bam = index_bam.indexed_bam,
            tumor_bam_bai = index_bam.indexed_bam_bai,
            region_file=chr_bed,
            insert_size=pindel_insert_size,
            tumor_sample_name=tumor_sample_name
        }
        call catOut as pindelCat {
          input: pindel_outs=[pindelTask.deletions, pindelTask.insertions, pindelTask.tandems, pindelTask.long_insertions, pindelTask.inversions]
        }
        call pindelToVcf as pindel2vcf {
            input:
            reference=reference,
            reference_fai=reference_fai,
            reference_dict=reference_dict,
            pindel_output_summary=pindelCat.pindel_out,
            ref_name = ref_name,
            ref_date = ref_date,
            min_supporting_reads = pindel_min_supporting_reads
        }
        call removeEndTags {
          input: vcf=pindel2vcf.pindel_vcf
        }
        call vcfSanitize as pindelSanitizeVcf {
          input: vcf=removeEndTags.processed_vcf
        }
        call bcftoolsNorm as pindelNormalize {
            input:
            reference=reference,
            reference_fai=reference_fai,
            vcf=pindelSanitizeVcf.sanitized_vcf,
            vcf_tbi=pindelSanitizeVcf.sanitized_vcf_tbi
        }
        call vtDecompose as pindelDecomposeVariants {
          input:
          vcf=pindelNormalize.normalized_vcf,
          vcf_tbi=pindelNormalize.normalized_vcf_tbi
        }

        scatter (caller_vcf in [mutect_pon2.annotated_vcf, vardict_pon2.annotated_vcf, lofreq_pon2.annotated_vcf]){
            call createFakeVcf as fake_vcf {
                input:
                vcf = caller_vcf,
                tumor_sample_name = tumor_sample_name
            }
        }

        call mergeVcf as mergeCallers {
            input:
            vcfs = fake_vcf.fake_vcf,
            vcf_tbis = fake_vcf.fake_vcf_tbi,
            merged_vcf_basename = "all_callers." + tumor_sample_name
        }

        call fpFilter as firstFilter {
          input:
          reference=reference,
          reference_fai=reference_fai,
          bam=index_bam.indexed_bam,
          vcf=mergeCallers.merged_vcf,
          sample_name=tumor_sample_name,
          min_var_freq=af_threshold,
          output_vcf_basename = "all_callers_full"
        }

        call selectVariants as hardFilter {
          input:
          reference=reference,
          reference_fai=reference_fai,
          reference_dict=reference_dict,
          vcf=firstFilter.filtered_vcf,
          vcf_tbi=firstFilter.filtered_vcf_tbi,
          exclude_filtered=true,
          output_vcf_basename = "all_callers_filtered"
        }

        scatter (pon_bam in pon_bams) {
            call mskGetBaseCounts as mskGetBaseCounts {
                input:
                reference = reference,
                reference_fai = reference_fai,
                reference_dict = reference_dict,
                normal_bam = pon_bam,
                pon_final_name = "all_callers." + tumor_sample_name + ".pon.pileup",
                vcf = mergeCallers.merged_vcf
            }
        }

        call bcftoolsMerge as pileup_merge {
            input:
                vcfs = mskGetBaseCounts.pileup,
                vcf_tbis = mskGetBaseCounts.pileup_tbi
        }

        scatter(custom_annotation in vep_custom_annotations) {
            call generateCustomString {
                input: custom_annotation = custom_annotation
            }
        }

        call vep {
          input:
              vcf=mergeCallers.merged_vcf,
              cache_dir_zip=vep_cache_dir_zip,
              reference=reference,
              reference_fai=reference_fai,
              reference_dict=reference_dict,
              plugins=vep_plugins,
              ensembl_assembly=vep_ensembl_assembly,
              ensembl_version= vep_ensembl_version,
              ensembl_species=vep_ensembl_species,
              synonyms_file=synonyms_file,
              custom_annotations = vep_custom_annotations,
              custom_annotation_string = generateCustomString.custom_string,
              custom_annotation_files = generateCustomString.custom_file,
              custom_annotation_files_tbi = generateCustomString.custom_file_tbi,
              coding_only=annotate_coding_only,
              everything=everything,
              pick=vep_pick
        }

        call normalFisher as mutect_call_R_fisher {
            input:
            vcf = mutect_pon2.annotated_vcf,
            pon = pileup_merge.merged_vcf,
            pon_tbi = pileup_merge.merged_vcf_tbi,
            p_value = pon_pvalue,
            caller = "mutect"
        }

        call indexVcf as mutect_index_pon_filtered_vcf {
            input: vcf = mutect_call_R_fisher.pon_filtered_vcf
        }

        call annotateVcf as mutect_annotate_vcf {
            input:
            vcf = mutect_index_pon_filtered_vcf.indexed_vcf,
            vcf_tbi = mutect_index_pon_filtered_vcf.indexed_vcf_tbi,
            fp_filter = firstFilter.filtered_vcf,
            fp_filter_tbi = firstFilter.filtered_vcf_tbi,
            vep = vep.annotated_vcf,
            vep_tbi = vep.annotated_vcf_tbi,
            caller_prefix = "mutect",
            sample_name = tumor_sample_name
        }

        call normalFisher as vardict_call_R_fisher {
            input:
            vcf = vardict_pon2.annotated_vcf,
            pon = pileup_merge.merged_vcf,
            pon_tbi = pileup_merge.merged_vcf_tbi,
            p_value = pon_pvalue,
            caller = "vardict"
        }

        call indexVcf as vardict_index_pon_filtered_vcf {
            input: vcf = vardict_call_R_fisher.pon_filtered_vcf
        }

        call annotateVcf as vardict_annotate_vcf {
            input:
            vcf = vardict_index_pon_filtered_vcf.indexed_vcf,
            vcf_tbi = vardict_index_pon_filtered_vcf.indexed_vcf_tbi,
            fp_filter = firstFilter.filtered_vcf,
            fp_filter_tbi = firstFilter.filtered_vcf_tbi,
            vep = vep.annotated_vcf,
            vep_tbi = vep.annotated_vcf_tbi,
            caller_prefix = "vardict",
            sample_name = tumor_sample_name
        }

        call normalFisher as lofreq_call_R_fisher {
            input:
            vcf = lofreq_pon2.annotated_vcf,
            pon = pileup_merge.merged_vcf,
            pon_tbi = pileup_merge.merged_vcf_tbi,
            p_value = pon_pvalue,
            caller = "lofreq"
        }

        call indexVcf as lofreq_index_pon_filtered_vcf {
            input: vcf = lofreq_call_R_fisher.pon_filtered_vcf
        }

        call annotateVcf as lofreq_annotate_vcf {
            input:
            vcf = lofreq_index_pon_filtered_vcf.indexed_vcf,
            vcf_tbi = lofreq_index_pon_filtered_vcf.indexed_vcf_tbi,
            fp_filter = firstFilter.filtered_vcf,
            fp_filter_tbi = firstFilter.filtered_vcf_tbi,
            vep = vep.annotated_vcf,
            vep_tbi = vep.annotated_vcf_tbi,
            caller_prefix = "lofreq",
            sample_name = tumor_sample_name
        }
    }

    call mergeVcf as merge_mutect_full {
        input:
            vcfs = mutectDecomposeVariants.decomposed_vcf,
            vcf_tbis = mutectDecomposeVariants.decomposed_vcf_tbi,
            merged_vcf_basename = "mutect_full." + tumor_sample_name
    }
    call mergeVcf as merge_mutect_pon {
        input:
            vcfs = mutect_annotate_vcf.pon_annotated_vcf,
            vcf_tbis = mutect_annotate_vcf.pon_annotated_vcf_tbi,
            merged_vcf_basename = "mutect." + tumor_sample_name + ".pileup.fisherPON"
    }
    call mergeVcf as merge_mutect_final {
        input:
            vcfs = mutect_annotate_vcf.final_annotated_vcf,
            vcf_tbis = mutect_annotate_vcf.final_annotated_vcf_tbi,
            merged_vcf_basename = "mutect." + tumor_sample_name + ".final.annotated"
    }

    call mergeVcf as merge_vardict_full {
        input:
            vcfs = vardictDecomposeVariants.decomposed_vcf,
            vcf_tbis = vardictDecomposeVariants.decomposed_vcf_tbi,
            merged_vcf_basename = "vardict_full." + tumor_sample_name
    }
    call mergeVcf as merge_vardict_pon {
        input:
            vcfs = vardict_annotate_vcf.pon_annotated_vcf,
            vcf_tbis = vardict_annotate_vcf.pon_annotated_vcf_tbi,
            merged_vcf_basename = "vardict." + tumor_sample_name + ".pileup.fisherPON"
    }
    call mergeVcf as merge_vardict_final {
        input:
            vcfs = vardict_annotate_vcf.final_annotated_vcf,
            vcf_tbis = vardict_annotate_vcf.final_annotated_vcf_tbi,
            merged_vcf_basename = "vardict." + tumor_sample_name + ".final.annotated"
    }

    call mergeVcf as merge_lofreq_full {
        input:
            vcfs = lofreqDecomposeVariants.decomposed_vcf,
            vcf_tbis = lofreqDecomposeVariants.decomposed_vcf_tbi,
            merged_vcf_basename = "lofreq_full." + tumor_sample_name
    }
    call mergeVcf as merge_lofreq_pon {
        input:
            vcfs = lofreq_annotate_vcf.pon_annotated_vcf,
            vcf_tbis = lofreq_annotate_vcf.pon_annotated_vcf_tbi,
            merged_vcf_basename = "lofreq." + tumor_sample_name + ".pileup.fisherPON"
    }
    call mergeVcf as merge_lofreq_final {
        input:
            vcfs = lofreq_annotate_vcf.final_annotated_vcf,
            vcf_tbis = lofreq_annotate_vcf.final_annotated_vcf_tbi,
            merged_vcf_basename = "lofreq." + tumor_sample_name + ".final.annotated"
    }

    call mergeVcf as merge_pindel_full {
        input:
            vcfs = pindelDecomposeVariants.decomposed_vcf,
            vcf_tbis = pindelDecomposeVariants.decomposed_vcf_tbi,
            merged_vcf_basename = "pindel_full." + tumor_sample_name
    }

    call mergeVcf as merge_pon {
        input:
            vcfs = pileup_merge.merged_vcf,
            vcf_tbis = pileup_merge.merged_vcf_tbi,
            merged_vcf_basename = tumor_sample_name + ".pon.total.counts"
    }

    call mergeVcf as merge_fp_filter {
        input:
            vcfs = firstFilter.filtered_vcf,
            vcf_tbis = firstFilter.filtered_vcf_tbi,
            merged_vcf_basename = tumor_sample_name + ".fpfilter"
    }

    call mergeVcf as merge_vep {
        input:
            vcfs = vep.annotated_vcf,
            vcf_tbis = vep.annotated_vcf_tbi,
            merged_vcf_basename = tumor_sample_name + ".vep"
    }

    call archerRAnnotate as annotateRMutect {
        input:
            vcf = merge_mutect_final.merged_vcf,
            caller = "mutect",
            bolton_bick_vars = bolton_bick_vars,
            mut2_bick = mut2_bick,
            mut2_kelly = mut2_kelly,
            matches2 = matches2,
            TSG_file = TSG_file,
            oncoKB_curated = oncoKB_curated,
            pd_annotation_file = pd_annotation_file,
            pan_myeloid = pan_myeloid,
            truncating = truncating,
            cosmic_dir_zip = cosmic_dir_zip
    }

    call archerRAnnotate as annotateRLofreq {
        input:
            vcf = merge_lofreq_final.merged_vcf,
            caller = "lofreq",
            bolton_bick_vars = bolton_bick_vars,
            mut2_bick = mut2_bick,
            mut2_kelly = mut2_kelly,
            matches2 = matches2,
            TSG_file = TSG_file,
            oncoKB_curated = oncoKB_curated,
            pd_annotation_file = pd_annotation_file,
            pan_myeloid = pan_myeloid,
            truncating = truncating,
            cosmic_dir_zip = cosmic_dir_zip
    }

    call archerRAnnotate as annotateRVardict {
        input:
            vcf = merge_vardict_final.merged_vcf,
            caller = "vardict",
            bolton_bick_vars = bolton_bick_vars,
            mut2_bick = mut2_bick,
            mut2_kelly = mut2_kelly,
            matches2 = matches2,
            TSG_file = TSG_file,
            oncoKB_curated = oncoKB_curated,
            pd_annotation_file = pd_annotation_file,
            pan_myeloid = pan_myeloid,
            truncating = truncating,
            cosmic_dir_zip = cosmic_dir_zip
    }

    call XGBModel as model {
        input:
            lofreq_tsv = annotateRLofreq.vcf_annotate_pd,
            mutect_tsv = annotateRMutect.vcf_annotate_pd,
            vardict_tsv  = annotateRVardict.vcf_annotate_pd,
            pindel_full_vcf = merge_pindel_full.merged_vcf,
            pon = merge_pon.merged_vcf,
            tumor_sample_name = tumor_sample_name
    }

    output {
        # Alignments
        File? aligned_bam = index_align_bam.indexed_bam
        File bqsr_bam = index_bam.indexed_bam

        # Tumor QC
        File tumor_insert_size_metrics = collectAllMetrics.insert_size_metrics
        File tumor_alignment_summary_metrics = collectAllMetrics.alignment_summary_metrics
        File tumor_hs_metrics = collectRoiHsMetrics.hs_metrics
        Array[File] tumor_per_target_coverage_metrics = select_all(collectPerTargetHsMetrics.per_target_coverage_metrics)
        Array[File] tumor_per_target_hs_metrics = collectPerTargetHsMetrics.hs_metrics
        Array[File] tumor_per_base_coverage_metrics = select_all(collectPerBaseHsMetrics.per_base_coverage_metrics)
        Array[File] tumor_per_base_hs_metrics = collectPerBaseHsMetrics.hs_metrics
        Array[File] tumor_summary_hs_metrics = collectSummaryHsMetrics.hs_metrics
        File tumor_flagstats = samtoolsFlagstat.flagstats
        File tumor_verify_bam_id_metrics = verifyBamId.verify_bam_id_metrics
        File tumor_verify_bam_id_depth = verifyBamId.verify_bam_id_depth

        # Mutect
        File mutect_full =  merge_mutect_full.merged_vcf                                # Raw Mutect Ouput
        File mutect_pon_annotated_vcf = merge_mutect_pon.merged_vcf                     # gnomAD Filtered + PoN Filtered + PoN2 Annotated
        File mutect_vep_annotated_vcf = merge_mutect_final.merged_vcf                   # gnomAD Filtered + PoN Filtered + PoN2 Annotated w/ VEP Annotation

        # Lofreq
        File lofreq_full = merge_lofreq_full.merged_vcf                                 # Raw Lofreq Ouput
        File lofreq_pon_annotated_vcf = merge_lofreq_pon.merged_vcf                     # gnomAD Filtered + PoN Filtered + PoN2 Annotated
        File lofreq_vep_annotated_vcf = merge_lofreq_final.merged_vcf                   # gnomAD Filtered + PoN Filtered + PoN2 Annotated w/ VEP Annotation

        # Vardict
        File vardict_full = merge_vardict_full.merged_vcf                               # Raw Vardict Ouput
        File vardict_pon_annotated_vcf = merge_vardict_pon.merged_vcf                   # gnomAD Filtered + PoN Filtered + PoN2 Annotated
        File vardict_vep_annotated_vcf = merge_vardict_final.merged_vcf                 # gnomAD Filtered + PoN Filtered + PoN2 Annotated w/ VEP Annotation

        # Pindel
        File pindel_full = merge_pindel_full.merged_vcf                               # Raw Vardict Ouput
        #File pindel_pon_annotated_vcf = merge_pindel_pon.merged_vcf                   # gnomAD Filtered + PoN Filtered + PoN2 Annotated
        #File pindel_vep_annotated_vcf = merge_pindel_final.merged_vcf                 # gnomAD Filtered + PoN Filtered + PoN2 Annotated w/ VEP Annotation

        File pon_total_counts = merge_pon.merged_vcf                                    # PoN Pileup Results
        File fpfilter_results = merge_fp_filter.merged_vcf
        File vep_results = merge_vep.merged_vcf

        #File gnomAD_exclude = get_gnomad_exclude.normalized_gnomad_exclude

        # R Things
        File mutect_annotate_pd = annotateRMutect.vcf_annotate_pd
        File lofreq_annotate_pd = annotateRLofreq.vcf_annotate_pd
        File vardict_annotate_pd = annotateRVardict.vcf_annotate_pd

        # Model
        File model_output = model.model_output
        File mutect_complex = model.mutect_complex
        File pindel_complex = model.pindel_complex
        File lofreq_complex = model.lofreq_complex
    }
}

task markIlluminaAdapters {
    input {
        File bam
    }

    Int cores = 1
    Int space_needed_gb = 10 + round(2*size(bam, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "mgibio/dna-alignment:1.0.0"
        memory: "6GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar MarkIlluminaAdapters INPUT=~{bam} OUTPUT=marked.bam METRICS=adapter_metrics.txt
    >>>

    output {
        File marked_bam = "marked.bam"
        File metrics = "adapter_metrics.txt"
    }
}

task umiAlign {
    input {
        File bam
        File reference
        File reference_fai
        File reference_dict
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_pac
        File reference_sa
    }

    Int cores = 8
    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_amb, reference_ann, reference_bwt, reference_pac, reference_sa], "GB")

    runtime {
      docker: "mgibio/dna-alignment:1.0.0"
      memory: "48GB"
      cpu: cores
      # 1 + just for a buffer
      # data_size*10 because bam uncompresses and streams to /dev/stdout and /dev/stdin, could have a couple flying at once
      bootDiskSizeGb: 10 + round(10*data_size + reference_size)
      disks: "local-disk ~{10 + round(10*data_size + reference_size)} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        /bin/bash /usr/bin/umi_alignment.sh ~{bam} ~{reference} ~{cores}
    >>>

    output {
        File aligned_bam = "aligned.bam"
        File aligned_bam_bai = "aligned.bai"
    }
}

task groupReadsAndConsensus {
    input {
        File bam
        Boolean umi_paired = true
    }

    Int cores = 1
    Int space_needed_gb = 10 + round(3*size(bam, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: "6GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -eo pipefail

        PAIRED=~{umi_paired}
        BAM=~{bam}

        if [ "$PAIRED" == true ]; then
            /usr/local/bin/fgbio GroupReadsByUmi --strategy paired --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input $BAM --output umi_grouped.bam
        else
            /usr/local/bin/fgbio GroupReadsByUmi --strategy adjacency --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input $BAM --output umi_grouped.bam
        fi
        /usr/local/bin/fgbio CallMolecularConsensusReads --input umi_grouped.bam --error-rate-pre-umi 45 --error-rate-post-umi 30 --min-input-base-quality 30 --min-reads 1 --output consensus_unaligned.bam
    >>>

    output {
        File consensus_bam = "consensus_unaligned.bam"
    }
}

task realign {
    input {
        File bam
        File reference
        File reference_fai
        File reference_dict
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_pac
        File reference_sa
    }

    Int cores = 8
    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_amb, reference_ann, reference_bwt, reference_pac, reference_sa], "GB")

    runtime {
        docker: "mgibio/dna-alignment:1.0.0"
        memory: "48GB"
        cpu: cores
        bootDiskSizeGb: 10 + round(10*data_size + reference_size)
        disks: "local-disk ~{10 + round(10*data_size + reference_size)} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /bin/bash /usr/bin/umi_realignment.sh ~{bam} ~{reference} ~{cores}
    >>>

    output {
        File consensus_aligned_bam = "realigned.bam"
        File consensus_aligned_bam_bai = "realigned.bai"
    }
}

task filterClipAndCollectMetrics {
    input {
        File bam
        File reference
        File reference_fai
        File reference_dict
        Array[Int] min_reads = [1]
        Float? max_read_error_rate = 0.05
        Float? max_base_error_rate = 0.1
        Int min_base_quality = 1
        Float max_no_call_fraction = 0.5
        File? target_intervals
        String description
        Boolean umi_paired = true
    }

    Int cores = 1
    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: "6GB"
        cpu: cores
        bootDiskSizeGb: 10 + round(3*data_size + reference_size)
        disks: "local-disk ~{10 + round(3*data_size + reference_size)} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -eo pipefail

        /usr/local/bin/fgbio FilterConsensusReads --input ~{bam} --output consensus_filtered.bam --ref ~{reference} --min-reads ~{sep=" " min_reads} --max-read-error-rate ~{max_read_error_rate} --max-base-error-rate ~{max_base_error_rate} --min-base-quality ~{min_base_quality} --max-no-call-fraction ~{max_no_call_fraction}
        /usr/local/bin/fgbio ClipBam --input consensus_filtered.bam --ref ~{reference} --clipping-mode Hard --clip-overlapping-reads true --output clipped.bam

        PAIRED=~{umi_paired}
        DESCRIPTION=~{description}
        INTERVALS=~{target_intervals}

        if [ "$PAIRED" == true ]; then
            if [[ -z "$DESCRIPTION" ]]; then
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --output duplex_seq.metrics --intervals ${INTERVALS}
                fi
            else
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --description ${DESCRIPTION} --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --description ${DESCRIPTION} --output duplex_seq.metrics --intervals ${INTERVALS}
                fi
            fi
        else
            echo "Sample not UMI Paired" > duplex_seq.metrics.txt
        fi
    >>>

    output {
        File clipped_bam = "clipped.bam"
        File clipped_bam_bai = "clipped.bai"
        Array[File] duplex_seq_metrics = glob("duplex_seq.metrics.*")
    }
}

task Metrics {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    String metric_accumulation_level
  }

  Int space_needed_gb = 10 + round(size([bam, bam_bai, reference, reference_fai, reference_dict],"GB"))
  Int preemptible = 1
  Int maxRetries = 0
  Int cores = 1
  runtime {
    memory: "6GB"
    docker: "broadinstitute/picard:2.23.6"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }
  String bamroot = basename(bam, ".bam")
  String summary_metrics = "~{bamroot}.AlignmentSummaryMetrics.txt"
  String size_metrics = "~{bamroot}.InsertSizeMetrics.txt"
  String size_histogram = "~{bamroot}.InsertSizeHistogram.pdf"
  command <<<
    /usr/bin/java -Xmx6g -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics INPUT=~{bam} OUTPUT=~{summary_metrics} REFERENCE_SEQUENCE=~{reference} METRIC_ACCUMULATION_LEVEL=~{metric_accumulation_level}
    /usr/bin/java -Xmx6g -jar /usr/picard/picard.jar CollectInsertSizeMetrics O=~{size_metrics} H=~{size_histogram} I=~{bam} REFERENCE_SEQUENCE=~{reference} METRIC_ACCUMULATION_LEVEL=~{metric_accumulation_level}
  >>>

  # TODO: how much space to allocate?
  output {
    File alignment_summary_metrics = summary_metrics
    File insert_size_histogram = size_histogram
    File insert_size_metrics = size_metrics
  }
}

task collectHsMetrics {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    String metric_accumulation_level

    File bait_intervals
    File target_intervals
    Boolean per_target_coverage = false
    Boolean per_base_coverage = false
    Int? minimum_base_quality
    Int? minimum_mapping_quality

    String output_prefix = "out"
  }

  Int space_needed_gb = 10 + round(size([bam, bam_bai, reference, reference_fai, reference_dict, bait_intervals, target_intervals], "GB"))
  Int preemptible = 1
  Int maxRetries = 0
  Int cores = 1
  runtime{
    memory: "6GB"
    docker: "broadinstitute/picard:2.23.6"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  String bamroot = basename(bam, ".bam")
  String hs_txt = "~{bamroot}.~{output_prefix}-HsMetrics.txt"
  String per_target_txt = "~{bamroot}.~{output_prefix}-PerTargetCoverage.txt"
  String per_base_txt = "~{bamroot}.~{output_prefix}-PerBaseCoverage.txt"
  command <<<
    /usr/bin/java -Xmx6g -jar /usr/picard/picard.jar CollectHsMetrics \
    O=~{hs_txt} \
    I=~{bam} \
    R=~{reference} \
    TARGET_INTERVALS=~{target_intervals} \
    METRIC_ACCUMULATION_LEVEL=~{metric_accumulation_level} \
    BAIT_INTERVALS=~{bait_intervals} \
    ~{if per_target_coverage then "PER_TARGET_COVERAGE=~{per_target_txt}" else ""} \
    ~{if per_base_coverage then "PER_BASE_COVERAGE=~{per_base_txt}" else ""} \
    ~{if defined(minimum_mapping_quality) then "MINIMUM_MAPPING_QUALITY=~{minimum_mapping_quality}" else ""} \
    ~{if defined(minimum_base_quality) then "MINIMUM_BASE_QUALITY=~{minimum_base_quality}" else ""}
  >>>

  # TODO: how much space to allocate?
  output {
    File hs_metrics = hs_txt
    File? per_target_coverage_metrics = per_target_txt
    File? per_base_coverage_metrics = per_base_txt
  }
}

task samtoolsFlagstat {
  input {
    File bam
    File bam_bai
  }

  Int space_needed_gb = 5 + round(size([bam, bam_bai], "GB")*2)
  Int preemptible = 1
  Int maxRetries = 0
  Int cores = 1
  runtime {
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  String outfile = basename(bam) + ".flagstat"
  command <<<
    /usr/local/bin/samtools flagstat ~{bam} > ~{outfile}
  >>>

  output {
    File flagstats = outfile
  }
}

task selectVariants {
  input {
    File reference
    File reference_fai
    File reference_dict

    File vcf
    File vcf_tbi

    File? interval_list
    Boolean exclude_filtered = false
    String output_vcf_basename = "select_variants"
    Array[String]? samples_to_include  # include genotypes from this sample

    # ENUM: one of ["INDEL", "SNP", "MIXED", "MNP", "SYMBOLIC", "NO_VARIATION"]
    String? select_type
  }

  Int space_needed_gb = 5 + round(size([vcf, vcf_tbi], "GB")*3 + size([reference, reference_fai, reference_dict, interval_list], "GB"))
  Int preemptible = 1
  Int maxRetries = 0
  Int cores = 1
  runtime {
    docker: "broadinstitute/gatk:4.2.0.0"
    memory: "6GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  String outfile = "~{output_vcf_basename}.vcf.gz"
  Array[String] samples = if defined(samples_to_include) then prefix("--sample-name ", select_first([samples_to_include])) else []
  command <<<
    /gatk/gatk --java-options -Xmx4g SelectVariants -O ~{outfile} \
    -R ~{reference} \
    --variant ~{vcf} \
    ~{if defined(interval_list) then "-L ~{interval_list}" else ""} \
    ~{if exclude_filtered then "--exclude-filtered ~{select_first([exclude_filtered])}" else ""} \
    ~{sep=" " samples} \
    ~{if defined(select_type) then "-select-type ~{select_type}" else ""}
  >>>

  output {
    File filtered_vcf = outfile
    File filtered_vcf_tbi = "~{outfile}.tbi"
  }
}

task verifyBamId {
  input {
    File vcf
    File bam
    File bam_bai
  }

  Int space_needed_gb = 5 + round(size([bam, bam_bai, vcf], "GB"))
  Int preemptible = 1
  Int maxRetries = 0
  Int cores = 1
  runtime {
    docker: "mgibio/verify_bam_id-cwl:1.1.3"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  String bamroot = basename(bam, ".bam")
  String outroot = "~{bamroot}.VerifyBamId"
  command <<<
    /usr/local/bin/verifyBamID --out ~{outroot} --vcf ~{vcf} --bam ~{bam} --bai ~{bam_bai}
  >>>

  output {
    File verify_bam_id_metrics = "~{outroot}.selfSM"
    File verify_bam_id_depth = "~{outroot}.depthSM"
  }
}

task bqsrApply {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bam_bai
    String output_name = "final"
    Array[File] known_sites
    Array[File] known_sites_tbi  # secondaryFiles...
    Array[String] intervals = ["chr1", "chr2", "chr3", "chr4", "chr5","chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
  }

  Int cores = 1
  Int space_needed_gb = 10 + round(size(known_sites, "GB") + size(known_sites_tbi, "GB") + size([reference, reference_fai, reference_dict], "GB") + size([bam, bam_bai], "GB") * 2)
  Int preemptible = 1
  Int maxRetries = 0
  runtime {
    cpu: cores
    docker: "broadinstitute/gatk:4.1.8.1"
    memory: "18GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    preemptible: preemptible
    maxRetries: maxRetries
  }

  command <<<
    /gatk/gatk --java-options -Xmx16g BaseRecalibrator -O bqsr.table ~{sep=" " prefix("-L ", intervals)} -R ~{reference} -I ~{bam} ~{sep=" " prefix("--known-sites ", known_sites)}
    /gatk/gatk --java-options -Xmx16g ApplyBQSR -O ~{output_name}.bam ~{sep=" " prefix("--static-quantized-quals ", [10, 20, 30])} -R ~{reference} -I ~{bam} -bqsr bqsr.table
  >>>

  output {
    File bqsr_bam = "~{output_name}.bam"
    File bqsr_bam_bai = "~{output_name}.bai"
  }
}

task indexBam {
  input { File bam }

  Int space_needed_gb = 10 + round(size(bam, "GB")*2)
  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    cpu: cores
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    memory: "6GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    preemptible: preemptible
    maxRetries: maxRetries
  }
  command <<<
    mv ~{bam} ~{basename(bam)}
    /usr/local/bin/samtools index ~{basename(bam)} ~{basename(bam)}.bai
  >>>
  output {
    File indexed_bam = basename(bam)
    File indexed_bam_bai = "~{basename(bam)}.bai"
  }
}

task mutectTumorOnly {
  input {
    File reference
    File reference_fai
    File reference_dict
    File? pon
    File? pon_tbi
    File? gnomad
    File? gnomad_tbi
    File tumor_bam
    File tumor_bam_bai

    File interval_list
  }

  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
  Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_list, "GB"))
  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    cpu: cores
    docker: "broadinstitute/gatk:4.2.0.0"
    memory: "32GB"
    bootDiskSizeGb: space_needed_gb
    disks: "local-disk ~{space_needed_gb} SSD"
    preemptible: preemptible
    maxRetries: maxRetries
  }

  String output_vcf = "mutect.filtered.vcf.gz"
  command <<<
    set -o pipefail
    set -o errexit

    /gatk/gatk Mutect2 --java-options "-Xmx20g" \
      -O mutect.vcf.gz \
      -R ~{reference} \
      -L ~{interval_list} \
      -I ~{tumor_bam} \
      ~{"--germline-resource " + gnomad} \
      ~{"-pon " + pon} \
      --read-index ~{tumor_bam_bai} \
      --f1r2-tar-gz mutect.f1r2.tar.gz \
      --max-reads-per-alignment-start 0 \

    /gatk/gatk LearnReadOrientationModel \
      -I mutect.f1r2.tar.gz \
      -O mutect.read-orientation-model.tar.gz

    /gatk/gatk FilterMutectCalls \
      -R ~{reference} \
      -V mutect.vcf.gz \
      --ob-priors mutect.read-orientation-model.tar.gz \
      -O ~{output_vcf} #Running FilterMutectCalls on the output vcf.
  >>>

  output {
    File vcf = output_vcf
    File vcf_tbi = output_vcf + ".tbi"
  }
}

task vcfSanitize {
  input {
    File vcf
  }

  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    memory: "6GB"
    cpu: cores
    docker: "mgibio/samtools-cwl:1.0.0"
    preemptible: preemptible
    maxRetries: maxRetries
  }

  # outbase should match in script but I don't want to risk changing it yet
  String outbase = basename(basename(vcf, ".gz"), ".vcf")
  command <<<
    set -eou pipefail

    # 1) removes lines containing non ACTGN bases, as they conflict with the VCF spec
    # and cause GATK to choke
    # 2) removes mutect-specific format tags containing underscores, which are likewise
    # illegal in the vcf spec
    base=`basename ~{vcf}`
    outbase=`echo $base | perl -pe 's/.vcf(.gz)?$//g'`
    echo "~{vcf}   $base    $outbase"
    if [[ "~{vcf}" =~ ".gz" ]];then
        #gzipped input
        gunzip -c "~{vcf}" | perl -a -F'\t' -ne 'print $_ if $_ =~ /^#/ || $F[3] !~ /[^ACTGNactgn]/' | sed -e "s/ALT_F1R2/ALTF1R2/g;s/ALT_F2R1/ALTF2R1/g;s/REF_F1R2/REFF1R2/g;s/REF_F2R1/REFF2R1/g" >$outbase.sanitized.vcf
    else
        #non-gzipped input
        cat "~{vcf}" | perl -a -F'\t' -ne 'print $_ if $_ =~ /^#/ || $F[3] !~ /[^ACTGNactgn]/' | sed -e "s/ALT_F1R2/ALTF1R2/g;s/ALT_F2R1/ALTF2R1/g;s/REF_F1R2/REFF1R2/g;s/REF_F2R1/REFF2R1/g" >$outbase.sanitized.vcf
    fi
    /opt/htslib/bin/bgzip $outbase.sanitized.vcf
    /usr/bin/tabix -p vcf $outbase.sanitized.vcf.gz
  >>>


  output {
    File sanitized_vcf = outbase + ".sanitized.vcf.gz"
    File sanitized_vcf_tbi = outbase + ".sanitized.vcf.gz.tbi"
  }
}

task bcftoolsNorm {
    input {
        File reference
        File reference_fai
        # File reference_dict

        File vcf
        File vcf_tbi
    }

    Int space_needed_gb = 10 + round(size([vcf, vcf_tbi], "GB") * 2 + size([reference, reference_fai], "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      memory: "6GB"
      docker: "kboltonlab/bst"
      disks: "local-disk ~{space_needed_gb} SSD"
      cpu: cores
      preemptible: preemptible
      maxRetries: maxRetries
    }
    command <<<
        /usr/local/bin/bcftools norm --multiallelics -any --output-type z --output bcftools_norm.vcf.gz ~{vcf} -f ~{reference}

        /usr/local/bin/tabix bcftools_norm.vcf.gz
    >>>

    output {
      File normalized_vcf = "bcftools_norm.vcf.gz"
      File normalized_vcf_tbi = "bcftools_norm.vcf.gz.tbi"
    }
}

task vtDecompose {
  input {
    File vcf
    File vcf_tbi
  }

  Int space_needed_gb = 5 + round(size([vcf, vcf_tbi], "GB")*2)
  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    memory: "6GB"
    docker: "kboltonlab/vt"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  command <<<
    /opt/vt/vt decompose -s -o decomposed.vcf.gz ~{vcf}

    /usr/bin/tabix decomposed.vcf.gz
  >>>

  output {
    File decomposed_vcf = "decomposed.vcf.gz"
    File decomposed_vcf_tbi = "decomposed.vcf.gz.tbi"
  }
}

task bcftoolsIsecComplement {
    input {
        File vcf
        File vcf_tbi
        File exclude_vcf
        File exclude_vcf_tbi
        String output_type = "z"
        String? output_vcf_name = "bcftools_isec.vcf"
    }

    Int space_needed_gb = 10 + 2*round(size([vcf, vcf_tbi, exclude_vcf, exclude_vcf_tbi], "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    runtime {
      cpu: cores
      docker: "kboltonlab/bst:latest"
      memory: "6GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        /usr/local/bin/bcftools isec -C -w1 ~{vcf} ~{exclude_vcf} --output-type ~{output_type} --output ~{output_vcf_name}.gz && /usr/local/bin/tabix ~{output_vcf_name}.gz

    >>>

    output {
        File complement_vcf = "~{output_vcf_name}.gz"
        File complement_vcf_tbi = "~{output_vcf_name}.gz.tbi"
    }
}

task vardictTumorOnly {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        String tumor_sample_name = "TUMOR"
        File interval_bed
        Float? min_var_freq = 0.005
    }

    Int cores = 16
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/vardictjava:1.0"
        memory: "96GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_name = "vardict"
    command <<<
        set -o pipefail
        set -o errexit



        /opt/VarDictJava/build/install/VarDict/bin/VarDict \
            -U -G ~{reference} \
            -f ~{min_var_freq} \
            -N ~{tumor_sample_name} \
            -b ~{tumor_bam} \
            -c 1 -S 2 -E 3 -g 4 ~{interval_bed} \
            -th ~{cores} | \
        /opt/VarDictJava/build/install/VarDict/bin/teststrandbias.R | \
        /opt/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl \
            -N "~{tumor_sample_name}" \
            -E \
            -f ~{min_var_freq} > ~{output_name}.vcf

        /usr/bin/bgzip ~{output_name}.vcf && /usr/bin/tabix ~{output_name}.vcf.gz
    >>>

    output {
        File vcf = "~{output_name}.vcf.gz"
        File vcf_tbi = "~{output_name}.vcf.gz.tbi"
    }
}

task lofreqTumorOnly {
    input {
      File reference
      File reference_fai

      File tumor_bam
      File tumor_bam_bai

      File interval_bed
      String? output_name = "lofreq.vcf"
    }

    Int cores = 4
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      docker: "kboltonlab/lofreq:latest"
      memory: "24GB"
      cpu: cores
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        set -o errexit
        set -o nounset

        /opt/lofreq/bin/lofreq indelqual --dindel -f ~{reference} -o output.indel.bam ~{tumor_bam}
        samtools index output.indel.bam
        /opt/lofreq/bin/lofreq call-parallel --pp-threads ~{cores} -A -B -f ~{reference} --call-indels --bed ~{interval_bed} -o ~{output_name} output.indel.bam --force-overwrite
        bgzip ~{output_name} && tabix ~{output_name}.gz
    >>>

    output {
        File vcf = "~{output_name}.gz"
        File vcf_tbi = "~{output_name}.gz.tbi"
    }
}

task lofreqReformat {
    input {
        File vcf
        String tumor_sample_name
    }

    Int space_needed_gb = 10 + 2*round(size(vcf, "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
      memory: "6GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        set -o errexit
        set -o nounset

        zcat ~{vcf} | grep "##" > lofreq.reformat.vcf
        echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">"  >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">"  >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">"  >> lofreq.reformat.vcf;
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{tumor_sample_name}" >> lofreq.reformat.vcf;
        zcat ~{vcf} | grep -v '#' | awk '{ n=split($8, semi, /;/); sample=""; format=""; for(i in semi){ split(semi[i], equ, /=/); if(i<=3){ if(i+1==4) sample=sample equ[2]; else sample=sample equ[2] ":"; if(i+1==4) format=format equ[1]; else format=format equ[1] ":";}}{print $0, "GT:"format, "0/1:"sample}}' OFS='\t' >> lofreq.reformat.vcf;
        bgzip lofreq.reformat.vcf && tabix lofreq.reformat.vcf.gz
    >>>

    output {
        File reformat_vcf = "lofreq.reformat.vcf.gz"
        File reformat_vcf_tbi = "lofreq.reformat.vcf.gz.tbi"
    }
}

task fpFilter {
  input {
    File reference
    File reference_fai
    File? reference_dict

    File bam
    File vcf

    String output_vcf_basename = "fpfilter"
    String sample_name = "TUMOR"
    Float? min_var_freq = 0.05
  }

  Int space_needed_gb = 10 + round(size(vcf, "GB")*2 + size([reference, reference_fai, reference_dict, bam], "GB"))
  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    memory: "6GB"
    bootDiskSizeGb: 25
    docker: "kboltonlab/fp_filter-wdl"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  String output_vcf = output_vcf_basename + ".vcf"
  command <<<
    /usr/bin/perl /usr/bin/fpfilter.pl --bam-readcount /usr/bin/bam-readcount --samtools /opt/samtools/bin/samtools --output ~{output_vcf} --reference ~{reference} --bam-file ~{bam} --vcf-file ~{vcf} --sample ~{sample_name} --min-var-freq ~{min_var_freq}

    /usr/bin/bgzip ~{output_vcf} && /usr/bin/tabix ~{output_vcf}.gz
  >>>

  output {
    File filtered_vcf = "~{output_vcf}.gz"
    File filtered_vcf_tbi = "~{output_vcf}.gz.tbi"
  }
}

task mskGetBaseCounts {
    input {
        File reference
        File reference_fai
        File reference_dict
        bam_and_bai normal_bam
        String? pon_final_name = "pon.pileup"
        File vcf
        Int? mapq = 5
        Int? baseq = 5
    }

    Int cores = 4
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float bam_size = size([normal_bam.bam, normal_bam.bai], "GB")
    Float vcf_size = size(vcf, "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + vcf_size)
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      docker: "kboltonlab/msk_getbasecounts:3.0"
      cpu: cores
      memory: "24GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        set -eou pipefail

        echo "REFERENCE: ~{reference_size}"
        echo "BAM: ~{bam_size}"
        echo "VCF: ~{vcf_size}"
        echo "SPACE_NEEDED: ~{space_needed_gb}"

        sample_name=$(samtools view -H ~{normal_bam.bam} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
        if [[ $(zgrep -v '#' ~{vcf} | wc -l) -lt 1 ]]; then
            echo "PRINT: ~{vcf}"
            printf "##fileformat=VCFv4.2\n" > ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth matching reference (REF) allele\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth matching alternate (ALT) allele\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=VF,Number=1,Type=Float,Description=\"Variant frequence (AD/DP)\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=DPP,Number=1,Type=Integer,Description=\"Depth on postitive strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=DPN,Number=1,Type=Integer,Description=\"Depth on negative strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=RDP,Number=1,Type=Integer,Description=\"Reference depth on postitive strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=RDN,Number=1,Type=Integer,Description=\"Reference depth on negative strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=ADP,Number=1,Type=Integer,Description=\"Alternate depth on postitive strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=ADN,Number=1,Type=Integer,Description=\"Alternate depth on negative strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=DPF,Number=1,Type=Integer,Description=\"Total fragment depth\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=RDF,Number=1,Type=Float,Description=\"Fragment depth matching reference (REF) allele\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=ADF,Number=1,Type=Float,Description=\"Fragment depth matching alternate (ALT) allele\">\n" >> ~{pon_final_name}.vcf
            printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample_name}\n" >> ~{pon_final_name}.vcf
        else
            echo "SCRIPT: ~{normal_bam.bam}"
            if [[ ~{vcf} == *.vcf.gz ]]; then
                bgzip -d ~{vcf}
                vcf_file=~{vcf}
                /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam ${sample_name}:~{normal_bam.bam} --vcf "${vcf_file%.*}" --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
            else
                /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam ${sample_name}:~{normal_bam.bam} --vcf ~{vcf} --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
            fi
        fi
        bgzip ~{pon_final_name}.vcf && tabix ~{pon_final_name}.vcf.gz
    >>>

    output {
        File pileup = "~{pon_final_name}.vcf.gz"
        File pileup_tbi = "~{pon_final_name}.vcf.gz.tbi"
    }
}

task vep {
  input {
    File vcf
    File cache_dir_zip
    File reference
    File reference_fai
    File reference_dict
    String ensembl_assembly
    String ensembl_version
    String ensembl_species
    Array[String] plugins
    Boolean coding_only = false
    Array[VepCustomAnnotation] custom_annotations = []
    Array[String]? custom_annotation_string = [""]
    Array[File]? custom_annotation_files = [""]
    Array[Array[File]?]? custom_annotation_files_tbi = [[""]]
    Boolean everything = true
    # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
    String pick = "flag_pick"
    String additional_args = "--pick_order canonical,rank,mane,ccds,appris,tsl,biotype,length --merged --buffer_size 1000 --af_gnomad"
    File? synonyms_file
  }

  Float cache_size = 3*size(cache_dir_zip, "GB")  # doubled to unzip
  Float vcf_size = 2*size(vcf, "GB")  # doubled for output vcf
  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Int space_needed_gb = 50 + round(reference_size + vcf_size + cache_size + size(synonyms_file, "GB"))
  runtime {
    memory: "64GB"
    bootDiskSizeGb: 30
    cpu: 4
    docker: "mgibio/vep_helper-cwl:vep_101.0_v2"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  String annotated_path = basename(basename(vcf, ".gz"), ".vcf") + "_annotated.vcf"
  String cache_dir = basename(cache_dir_zip, ".zip")

  command <<<
    custom_string="~{sep=" " custom_annotation_string}"
    if [[ -z "$custom_string" ]]; then
        echo ${custom_string} >> custom_string_validation.txt
    else
        for file_path in ~{sep=" " custom_annotation_files}; do
            custom_string=$(awk -v srch="<CUSTOM_FILE>" -v repl="$file_path" '!x{x=sub(srch,repl)}{print $0}' <<< $custom_string)
        done
        echo ${custom_string} >> custom_string_validation.txt
    fi

    #mkdir ~{cache_dir} && unzip -qq ~{cache_dir_zip} -d ~{cache_dir}
    unzip -qq ~{cache_dir_zip}

    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
    --format vcf \
    --vcf \
    --fork 4 \
    --terms SO \
    --transcript_version \
    --offline \
    --cache \
    --symbol \
    -o ~{annotated_path} \
    -i ~{vcf} \
    ~{if defined(synonyms_file) then "--synonyms ~{synonyms_file}" else ""} \
    --sift p \
    --polyphen p \
    ~{if coding_only then "--coding_only" else ""} \
    --~{pick} \
    --dir ~{cache_dir} \
    --fasta ~{reference} \
    ~{sep=" " prefix("--plugin ", plugins)}  \
    ~{if everything then "--everything" else ""} \
    --assembly ~{ensembl_assembly} \
    --cache_version ~{ensembl_version} \
    --species ~{ensembl_species} \
    ~{additional_args} \
    ${custom_string}

    bgzip ~{annotated_path} && tabix ~{annotated_path}.gz
  >>>

  output {
    File annotated_vcf = "~{annotated_path}.gz"
    File annotated_vcf_tbi = "~{annotated_path}.gz.tbi"
    File vep_summary = annotated_path + "_summary.html"
  }
}

task generateCustomString {
    input { VepCustomAnnotation custom_annotation }
    runtime { docker: "ubuntu:xenial" }
    command <<<
        /bin/echo '~{if custom_annotation.annotation.check_existing then "--check_existing" else ""} --custom <CUSTOM_FILE>,~{custom_annotation.annotation.name},~{custom_annotation.annotation.data_format},~{custom_annotation.method},~{if custom_annotation.force_report_coordinates then 1 else 0},~{sep="," custom_annotation.annotation.vcf_fields}'
    >>>
    output {
        String custom_string = read_string(stdout())
        File custom_file = custom_annotation.annotation.file
        Array[File]? custom_file_tbi = custom_annotation.annotation.secondary_files
    }
}

task pon2Percent {
    input {
        File vcf
        File vcf2PON
        File vcf2PON_tbi
        String caller = "caller"
        String sample_name = "tumor"
    }

    Int space_needed_gb = 10 + round(size([vcf, vcf2PON, vcf2PON_tbi],"GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/bst:latest"
      memory: "6GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        set -eou pipefail

        export name=~{caller}.~{sample_name}.pon2.annotated.vcf.gz

        printf "##INFO=<ID=PON_2AT2_percent,Number=1,Type=Integer,Description=\"If 2 PoN samples have variant at >=2 percent\">\n" > pon2.header;
        printf "##INFO=<ID=PON_NAT2_percent,Number=1,Type=Integer,Description=\"Number of samples with variant at >=2 percent\">\n" >> pon2.header;
        printf "##INFO=<ID=PON_MAX_VAF,Number=1,Type=Float,Description=\"The maximum VAF found in the PoN Samples\">\n" >> pon2.header;
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t1\t%INFO/NS\t%INFO/max_VAF\n" ~{vcf2PON} > normal2.txt
        bgzip -f normal2.txt
        tabix -f -s1 -b2 -e2 normal2.txt.gz
        bcftools annotate -a normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent,PON_NAT2_percent,PON_MAX_VAF ~{vcf} -Oz -o $name
        tabix $name
    >>>

    output {
        File annotated_vcf = "~{caller}.~{sample_name}.pon2.annotated.vcf.gz"
        File annotated_vcf_tbi = "~{caller}.~{sample_name}.pon2.annotated.vcf.gz.tbi"
    }
}

task splitBedToChr {
  input {
    File interval_bed
  }

  Int cores = 1
  Int space_needed_gb = 10 + round(size(interval_bed, "GB")*2)
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    cpu: cores
    memory: "6GB"
    docker: "kboltonlab/bst:latest"
    disks: "local-disk ~{space_needed_gb} SSD"
    preemptible: preemptible
    maxRetries: maxRetries
  }

  command <<<

    intervals=$(awk '{print $1}' ~{interval_bed} | uniq)
    for chr in ${intervals}; do
        grep -w $chr ~{interval_bed} > ~{basename(interval_bed, ".bed")}_${chr}.bed
    done
  >>>

  output {
    Array[File] split_chr = glob(basename(interval_bed, ".bed")+"_*.bed")
  }
}

task mergeVcf {
  input {
    Array[File] vcfs
    Array[File] vcf_tbis
    String merged_vcf_basename = "merged"
  }

  Int space_needed_gb = 10 + round(2*(size(vcfs, "GB") + size(vcf_tbis, "GB")))
  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    docker: "kboltonlab/bst:latest"
    memory: "6GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  String output_file = merged_vcf_basename + ".vcf.gz"
  command <<<
    /usr/local/bin/bcftools concat --allow-overlaps --remove-duplicates --output-type z -o ~{output_file} ~{sep=" " vcfs}

    /usr/local/bin/tabix ~{output_file}
  >>>

  output {
    File merged_vcf = output_file
    File merged_vcf_tbi = "~{output_file}.tbi"
  }
}

task createFakeVcf {
    input {
        File vcf
        String tumor_sample_name
    }

    Int cores = 1
    Int space_needed_gb = 10 + round(2*size(vcf, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
        memory: "6GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -o errexit
        set -o nounset

        echo -e "##fileformat=VCFv4.2" > fake.vcf
        echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> fake.vcf;
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{tumor_sample_name}" >> fake.vcf;
        zcat ~{vcf} | grep -v '#' | awk '{print $1, $2, $3, $4, $5, $6, "PASS\t.\tGT\t0/1"}' OFS='\t' >> fake.vcf;
        bgzip fake.vcf && tabix fake.vcf.gz
    >>>

    output {
        File fake_vcf = "fake.vcf.gz"
        File fake_vcf_tbi = "fake.vcf.gz.tbi"
    }
}

task bcftoolsFilterBcbio {
    input {
        File vcf
        File vcf_tbi
        String filter_flag = "exclude"
        String filter_string
        String? output_vcf_prefix = "bcftools_filter"
        String output_type = "z"
    }

    Int space_needed_gb = 10 + 2*round(size(vcf, "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/bst:latest"
      memory: "6GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    String ff = if filter_flag == "include" then "-i" else "-e"
    command <<<
        /usr/local/bin/bcftools filter ~{ff} "~{filter_string}" ~{vcf} --output-type ~{output_type} --output ~{output_vcf_prefix}.vcf.gz -s "BCBIO" -m+

        /usr/local/bin/tabix ~{output_vcf_prefix}.vcf.gz

    >>>

    output {
        File filtered_vcf = "~{output_vcf_prefix}.vcf.gz"
        File filtered_vcf_tbi = "~{output_vcf_prefix}.vcf.gz.tbi"
    }
}

task bcftoolsMerge {
  input {
    Array[File] vcfs
    Array[File] vcf_tbis
    String merged_vcf_basename = "merged"
  }

  Int space_needed_gb = 10 + round(2*(size(vcfs, "GB") + size(vcf_tbis, "GB")))
  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    docker: "kboltonlab/bst:latest"
    memory: "6GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  String output_file = merged_vcf_basename + ".vcf.gz"
  command <<<
    /usr/local/bin/bcftools merge --output-type z -o ~{output_file} ~{sep=" " vcfs}
    /usr/local/bin/tabix ~{output_file}
  >>>

  output {
    File merged_vcf = output_file
    File merged_vcf_tbi = "~{output_file}.tbi"
  }
}

task indexVcf {
  input {
    File vcf
  }

  Int space_needed_gb = 10 + round(3*size(vcf, "GB"))
  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    memory: "6GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  command <<<
    cp ~{vcf} ~{basename(vcf)}
    /usr/local/bin/tabix -p vcf ~{basename(vcf)}
  >>>
  output {
    File indexed_vcf = basename(vcf)
    File indexed_vcf_tbi = basename(vcf) + ".tbi"
  }
}

task normalFisher {
    input {
        File vcf
        File pon
        File pon_tbi
        String caller = "caller"
        String? p_value = "0.05"
    }


    Int space_needed_gb = 10 + round(size([vcf, pon, pon_tbi], "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/sam_bcftools_tabix_bgzip:1.0"
      memory: "6GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        set -eou pipefail

        if [[ "~{vcf}" == *.gz ]]; then
            name=$(basename ~{vcf} .vcf.gz)
        else
            name=$(basename ~{vcf} .vcf)
        fi

        bcftools +fill-tags -Oz -o RD.vcf.gz ~{pon} -- -t "PON_RefDepth=sum(RD)"
        bcftools +fill-tags -Oz -o RD_AD.vcf.gz RD.vcf.gz -- -t "PON_AltDepth=sum(AD)" && tabix RD_AD.vcf.gz

        printf "##INFO=<ID=PON_FISHER,Number=1,Type=Float,Description=\"P-value from Fisher's exact test with totals from PoN\">" > fisher.header;
        printf "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">" > sample.header;

        sample=`bcftools query -l ~{vcf}`
        bcftools view -H ~{vcf} | awk -v sampleID=$sample '{print $1, $2, $3, $4, $5, sampleID}' OFS='\t' > $sample.name;
        bgzip $sample.name;
        tabix $sample.name.gz -s1 -b2 -e2;
        bcftools annotate -a $sample.name.gz -h sample.header -c CHROM,POS,-,REF,ALT,SAMPLE ~{vcf} -Oz -o $name.sample.vcf.gz && tabix $name.sample.vcf.gz;

        ## Varscan has AD and RD instead of comma sep AD field
        ## you can't double quote for string match in bash like you can in zsh so need to make it a variable
        pat="[Vv]arscan"

        ## Lofreq has DP4 which splits into RefFwd, RefRev, AltFwd, AltRev
        patt="[Ll]ofreq"
        if [[ ~{caller} =~ $pat ]]
        then
            echo '
            #!/usr/bin/env Rscript

            args = commandArgs(trailingOnly=TRUE)

            if (length(args)==0) {
            stop("At least one argument must be supplied (input file).n", call.=FALSE)
            } else if (length(args)==1) {
            # default output file
            stop("Must supply (output file).n", call.=FALSE)
            }

            df = read.table(args[1], header=F, colClasses = c("character", "numeric", "character", "character", "numeric", "numeric", "numeric", "numeric"))

            if (length(colnames(df)) != 8) {
            stop("Must supply file with 8 columns: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%RD]\t[%AD]", call.=FALSE)
            }

            df$fisher.exact.pval <- apply(df, 1, function(x) {
            x <- as.numeric(x[-c(1,2,3,4)])
            if ((x[1]+x[2]==0) | (x[3]+x[4]==0)){
                return(0)
            } else if (x[2]==0 & x[1]!=0) {
                return(0)
            } else if ((x[1]==0 & x[2]!=0) | (x[3]==0 & x[4]!=0)) {
                return(1)
            } else if (x[2]/(x[1]+x[2]) >= x[4]/(x[3]+ x[4])) {
                return(1)
            } else {
                return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
            }
            })
            write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
            ' > fisherTestInput.R
            bcftools annotate -a RD_AD.vcf.gz -c PON_RefDepth,PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%RD]\t[%AD]\n' $name.sample.pileup.vcf.gz > $name.fisher.input;
        elif [[ ~{caller} =~ $patt ]]
        then
            echo '
            #!/usr/bin/env Rscript

            args = commandArgs(trailingOnly=TRUE)

            if (length(args)==0) {
            stop("At least one argument must be supplied (input file).n", call.=FALSE)
            } else if (length(args)==1) {
            # default output file
            stop("Must supply (output file).n", call.=FALSE)
            }

            df = read.table(args[1], header=F, colClasses = c("character", "numeric", "character", "character", "numeric", "numeric", "character"))

            #https://statisticsglobe.com/split-data-frame-variable-into-multiple-columns-in-r
            df <- cbind(df, data.frame(do.call("rbind", strsplit(as.character(df$V7), ",", fixed = TRUE))))[,-7]

            if (length(colnames(df)) < 8) {
            stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t%INFO/DP4", call.=FALSE)
            }

            df$fisher.exact.pval <- apply(df, 1, function(x) {
            x <- as.numeric(x[-c(1,2,3,4)])
            # Remember, Lofreq splits DP4 into RefFwd, RefRev and AltFwd, AltRev so technically ref = x[3] + x[4] and alt = x[5] + x[6]
            ref = x[3] + x[4]
            alt = x[5] + x[6]
            if ((x[1]+x[2]==0) | (ref+alt==0)){
                return(0)
            } else if (x[2]==0 & x[1]!=0) {
                return(0)
            } else if ((x[1]==0 & x[2]!=0) | (ref==0 & alt!=0)) {
                return(1)
            } else if (x[2]/(x[1]+x[2]) >= alt/(ref+alt)) {
                return(1)
            } else {
                return(fisher.test(matrix(c(x[1], x[2], ref, alt), ncol=2))$p.value)
            }
            })
            write.table(df[, -c(9:10)], file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
            ' > fisherTestInput.R
            bcftools annotate -a RD_AD.vcf.gz -c PON_RefDepth,PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t%INFO/DP4\n' $name.sample.pileup.vcf.gz > $name.fisher.input;
        else
            echo '
            #!/usr/bin/env Rscript

            args = commandArgs(trailingOnly=TRUE)

            if (length(args)==0) {
            stop("At least one argument must be supplied (input file).n", call.=FALSE)
            } else if (length(args)==1) {
            # default output file
            stop("Must supply (output file).n", call.=FALSE)
            }

            df = read.table(args[1], header=F, colClasses = c("character", "numeric", "character", "character", "numeric", "numeric", "character"))

            #https://statisticsglobe.com/split-data-frame-variable-into-multiple-columns-in-r
            df <- cbind(df, data.frame(do.call("rbind", strsplit(as.character(df$V7), ",", fixed = TRUE))))[,-7]
            if (length(colnames(df)) != 8) {
            stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%AD]", call.=FALSE)
            }

            df$fisher.exact.pval <- apply(df, 1, function(x) {
            x <- as.numeric(x[-c(1,2,3,4)])
            if ((x[1]+x[2]==0) | (x[3]+x[4]==0)){
                return(0)
            } else if (x[2]==0 & x[1]!=0) {
                return(0)
            } else if ((x[1]==0 & x[2]!=0) | (x[3]==0 & x[4]!=0)) {
                return(1)
            } else if (x[2]/(x[1]+x[2]) >= x[4]/(x[3]+x[4])) {
                return(1)
            } else {
                return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
            }
            })
            write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
            ' > fisherTestInput.R
            bcftools annotate -a RD_AD.vcf.gz -c PON_RefDepth,PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%AD]\n' $name.sample.pileup.vcf.gz > $name.fisher.input;

        fi
        chmod u+x fisherTestInput.R

        # Depending on how we split, we might have caller_vcf that doesn't have any variants called
        if [ -s $name.fisher.input ]; then
            LC_ALL=C.UTF-8 Rscript --vanilla ./fisherTestInput.R $name.fisher.input $name.fisher.output
            bgzip -f $name.fisher.output
            tabix -f -s1 -b2 -e2 $name.fisher.output.gz
            bcftools annotate -a $name.fisher.output.gz -h fisher.header -c CHROM,POS,REF,ALT,-,-,-,-,PON_FISHER $name.sample.pileup.vcf.gz -Oz -o $name.pileup.fisherPON.vcf.gz && tabix $name.pileup.fisherPON.vcf.gz
            bcftools filter -i "INFO/PON_FISHER<~{p_value}" $name.pileup.fisherPON.vcf.gz -Oz -o $name.filtered.pileup.fisherPON.vcf.gz && tabix $name.filtered.pileup.fisherPON.vcf.gz
        else
            bcftools annotate -h fisher.header $name.sample.pileup.vcf.gz -Oz -o $name.pileup.fisherPON.vcf.gz && tabix $name.pileup.fisherPON.vcf.gz
            bcftools filter -i "INFO/PON_FISHER<~{p_value}" $name.pileup.fisherPON.vcf.gz -Oz -o $name.filtered.pileup.fisherPON.vcf.gz && tabix $name.filtered.pileup.fisherPON.vcf.gz
        fi
    >>>

    output {
        File pon_vcf = select_first([basename(vcf, ".vcf.gz") + ".pileup.fisherPON.vcf.gz",basename(vcf, ".vcf") + ".pileup.fisherPON.vcf.gz"])
        File pon_vcf_tbi = select_first([basename(vcf, ".vcf.gz") + ".pileup.fisherPON.vcf.gz.tbi",basename(vcf, ".vcf") + ".pileup.fisherPON.vcf.gz.tbi"])
        File pon_filtered_vcf = select_first([basename(vcf, ".vcf.gz") + ".filtered.pileup.fisherPON.vcf.gz",basename(vcf, ".vcf") + ".filtered.pileup.fisherPON.vcf.gz"])
        File pon_filtered_vcf_tbi = select_first([basename(vcf, ".vcf.gz") + ".filtered.pileup.fisherPON.vcf.gz.tbi",basename(vcf, ".vcf") + ".filtered.pileup.fisherPON.vcf.gz.tbi"])
    }
}

task annotateVcf {
    input {
        File vcf
        File vcf_tbi
        File fp_filter
        File fp_filter_tbi
        File vep
        File vep_tbi
        String caller_prefix
        String sample_name
    }

    Int space_needed_gb = 10 + 2*round(size([vcf, vcf_tbi, fp_filter, fp_filter_tbi, vep], "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      memory: "6GB"
      docker: "kboltonlab/bst"
      disks: "local-disk ~{space_needed_gb} SSD"
      bootDiskSizeGb: space_needed_gb
      cpu: cores
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        zcat ~{fp_filter} | grep '##' | tail -n +4  > fp_filter.header;
        zcat ~{vep} | grep '##' | tail -n +3 > vep.header;

        bcftools annotate -a ~{fp_filter} -h fp_filter.header -c =FILTER ~{vcf} -Oz -o ~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz
        tabix ~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz
        bcftools annotate -a ~{vep} -h vep.header -c CSQ ~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz -Oz -o ~{caller_prefix}.~{sample_name}.final.annotated.vcf.gz
        tabix ~{caller_prefix}.~{sample_name}.final.annotated.vcf.gz
    >>>

    output {
        File final_annotated_vcf = "~{caller_prefix}.~{sample_name}.final.annotated.vcf.gz"
        File final_annotated_vcf_tbi = "~{caller_prefix}.~{sample_name}.final.annotated.vcf.gz.tbi"
        File pon_annotated_vcf = "~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz"
        File pon_annotated_vcf_tbi = "~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz.tbi"
    }
}

task fastqToBam {
    input {
        File fastq1
        File fastq2
        String sample_name
        String library_name
        String platform_unit
        String platform
    }

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = 10 + round(2*size([fastq1, fastq2], "GB"))

    runtime {
        docker: "mgibio/dna-alignment:1.0.0"
        memory: "6GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar FastqToSam FASTQ=~{fastq1} FASTQ2=~{fastq2} SAMPLE_NAME=~{sample_name} LIBRARY_NAME=~{library_name} PLATFORM_UNIT=~{platform_unit} PLATFORM=~{platform} OUTPUT=unaligned.bam
    >>>

    output {
        File bam = "unaligned.bam"
    }
}

# Pindel Stuff
task pindelTumorOnly {
    input {
        File reference
        File reference_fai
        File reference_dict
        File tumor_bam
        File tumor_bam_bai
        File region_file
        String tumor_sample_name
        String? chromosome
        Int insert_size = 400
    }

    Int cores = 4
    Int preemptible = 1
    Int maxRetries = 5
    Int space_needed_gb = 10 + round(size([reference, reference_fai, reference_dict, tumor_bam, tumor_bam_bai, region_file], "GB"))
    runtime {
        bootDiskSizeGb: 100
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        docker: "mgibio/cle:v1.4.2"
        memory: "24GB"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        echo -e "~{tumor_bam}\t~{insert_size}\t~{tumor_sample_name}" >> pindel.config

        /usr/bin/pindel -i pindel.config -w 30 -T ~{cores} -o all -f ~{reference} \
        ~{if defined(chromosome) then "-c ~{chromosome}" else ""} \
        ~{if defined(region_file) then "-j ~{region_file}" else ""}
    >>>

  # TODO: how much space to allocate?
    output {
        File deletions = "all_D"
        File insertions = "all_SI"
        File tandems = "all_TD"
        File long_insertions = "all_LI"
        File inversions = "all_INV"
    }
}

task catOut {
  input {
    Array[File] pindel_outs
  }

  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0
  Int space_needed_gb = 10 + round(size(pindel_outs, "GB")*2)
  runtime {
    memory: "6GB"
    docker: "ubuntu:xenial"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  command <<<
    /bin/cat ~{sep=" " pindel_outs} | /bin/grep "ChrID" /dev/stdin > pindel.head
  >>>

  output {
    File pindel_out = "pindel.head"
  }
}

task pindelToVcf {
    input {
        File pindel_output_summary
        File reference
        File reference_fai
        File reference_dict
        String? ref_name = "GRCh38DH"
        String? ref_date = "20161216"
        Int? min_supporting_reads = 3
        String? output_name = "pindel.vcf"
    }
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = 10 + round(size([reference, reference_fai, reference_dict, pindel_output_summary], "GB"))
    runtime {
      memory: "6GB"
      docker: "mgibio/cle:v1.3.1"
      disks: "local-disk ~{space_needed_gb} SSD"
      cpu: cores
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        /usr/bin/pindel2vcf -G -p ~{pindel_output_summary} -r ~{reference} -R ~{ref_name} -e ~{min_supporting_reads} -d ~{ref_date} -v ~{output_name}
    >>>

    output {
        File pindel_vcf = "~{output_name}"
    }
}

task removeEndTags {
  input {
    File vcf
  }
  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0
  Int space_needed_gb = 10 + round(size(vcf, "GB")*2)
  runtime {
    memory: "6GB"
    docker: "kboltonlab/bst:latest"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  String outfile = "pindel.noend.vcf.gz"
  command <<<
      /usr/local/bin/bgzip -c ~{vcf} > pindel.vcf.gz
      /usr/local/bin/tabix -p vcf pindel.vcf.gz
      /usr/local/bin/bcftools annotate -x INFO/END -Oz -o ~{outfile} pindel.vcf.gz
      /usr/local/bin/tabix -p vcf ~{outfile}
  >>>

  output {
    File processed_vcf = outfile
    File processed_vcf_tbi = "~{outfile}.tbi"
  }
}

task archerRAnnotate {
    input {
        File vcf
        String caller = "caller"
        File bolton_bick_vars
        File mut2_bick
        File mut2_kelly
        File matches2
        File TSG_file
        File oncoKB_curated
        File pd_annotation_file
        File pan_myeloid
        File truncating
        File cosmic_dir_zip
    }

    Float file_size = size([vcf, bolton_bick_vars, mut2_bick, mut2_kelly, matches2, TSG_file, oncoKB_curated, pd_annotation_file, pan_myeloid], "GB")
    Float cosmic_size = 3*size(cosmic_dir_zip, "GB")
    Int space_needed_gb = 10 + round(file_size + cosmic_size)
    Int cores = 8
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/r_docker_ichan:latest"
      memory: "6GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    String cosmic_dir = basename(cosmic_dir_zip, ".zip")

    command <<<
        set -eou pipefail

        if [[ "~{vcf}" == *.gz ]]; then
            name=$(basename ~{vcf} .vcf.gz)
        else
            name=$(basename ~{vcf} .vcf)
        fi

        unzip -qq ~{cosmic_dir_zip}

        LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/ArcherAnnotationScript.R --input ~{vcf} --out ${name}.tsv --caller ~{caller} \
        --bolton_bick_vars ~{bolton_bick_vars} \
        --mut2_bick ~{mut2_bick} \
        --mut2_kelly ~{mut2_kelly} \
        --matches2 ~{matches2} \
        --TSG_file ~{TSG_file} \
        --oncoKB_curated ~{oncoKB_curated} \
        --pd_annotation_file ~{pd_annotation_file} \
        --pan_myeloid ~{pan_myeloid} \
        --truncating ~{truncating} \
        --cosmic_dir ~{cosmic_dir}
    >>>

    output {
        File vcf_annotate_pd = basename(vcf, ".vcf.gz") + ".tsv"
    }
}

task XGBModel {
    input {
        File lofreq_tsv
        File mutect_tsv
        File vardict_tsv
        File pindel_full_vcf
        File pon
        String tumor_sample_name
    }

    Float file_size = size([lofreq_tsv, mutect_tsv, vardict_tsv, pindel_full_vcf, pon], "GB")
    Int space_needed_gb = 10 + round(file_size)
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/predict_xgb:latest"
      memory: "6GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        /opt/bin/xgbappcompiled/bin/xgbapp ~{lofreq_tsv} ~{mutect_tsv} ~{vardict_tsv} ~{pindel_full_vcf} ~{pon} ""
        echo "Model Finished..."
    >>>

    output {
        File model_output = "output_~{tumor_sample_name}.tsv.gz"
        File mutect_complex = "mutect_complex_~{tumor_sample_name}.tsv.gz"
        File pindel_complex = "pindel_complex_~{tumor_sample_name}.tsv.gz"
        File lofreq_complex = "lofreq_complex_~{tumor_sample_name}.tsv.gz"
    }
}
