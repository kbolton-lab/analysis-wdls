version 1.0

import "../types.wdl"

import "../subworkflows/archer_fastq_format.wdl" as fqf
import "../subworkflows/molecular_alignment.wdl" as ma
import "../subworkflows/qc_exome.wdl" as qe
import "../subworkflows/PoN_filter.wdl" as pf
import "../subworkflows/fp_filter_no_norm.wdl" as ffnn
import "../subworkflows/mutect_noFp.wdl" as m
import "../subworkflows/lofreq_noFp.wdl" as l
import "../subworkflows/vardict_noFp.wdl" as v
import "../subworkflows/annotate_caller.wdl" as ac

import "../tools/fastq_to_bam.wdl" as ftb
import "../tools/bqsr_apply.wdl" as ba
import "../tools/index_bam.wdl" as ib
import "../tools/bcftools_isec_complement.wdl" as bic
import "../tools/vep.wdl" as vep
import "../tools/pon2percent.wdl" as pp
import "../tools/split_bed_to_chr.wdl" as sbtc
import "../tools/merge_vcf.wdl" as mv
import "../tools/create_fake_vcf.wdl" as cfv
import "../tools/archer_R_annotate.wdl" as ara
import "../tools/xgb_model.wdl" as xm

workflow archerdx {
    input {
        # Pipeline
        Int scatter_count = 20

        # Sequence and BAM Information
        File bam
        File bam_bai
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
        String? pon_pvalue = "2.114164905e-6"

        # Pindel
        #Int pindel_insert_size = 400
        #String? ref_name = "GRCh38DH"
        #String? ref_date = "20161216"
        #Int? pindel_min_supporting_reads = 3

        String bcbio_filter_string = "((FMT/AF * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 3.0) || (FMT/DP < 6500) || (INFO/QUAL < 27)))"

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

        # gnomAD
        Float filter_gnomADe_maximum_population_allele_frequency = 0.005
        File normalized_gnomad_exclude
        File normalized_gnomad_exclude_tbi
        String filter_flag = "include"
    }


    call ba.bqsrApply as bqsr {
        input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        bam = bam,
        bam_bai = bam_bai,
        intervals = bqsr_intervals,
        known_sites = bqsr_known_sites,
        known_sites_tbi = bqsr_known_sites_tbi
    }

    call ib.indexBam as index_bam {
        input:
        bam = bqsr.bqsr_bam
    }

    call qe.qcExome as tumor_qc {
        input:
        bam = index_bam.indexed_bam,
        bam_bai = index_bam.indexed_bam_bai,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        bait_intervals = bait_intervals,
        target_intervals = target_intervals,
        per_base_intervals = per_base_intervals,
        per_target_intervals = per_target_intervals,
        summary_intervals = summary_intervals,
        omni_vcf = omni_vcf,
        omni_vcf_tbi = omni_vcf_tbi,
        picard_metric_accumulation_level = picard_metric_accumulation_level,
        minimum_mapping_quality = qc_minimum_mapping_quality,
        minimum_base_quality = qc_minimum_base_quality
    }

    call sbtc.splitBedToChr as split_bed_to_chr {
        input:
        interval_bed = target_bed
    }

    scatter (chr_bed in split_bed_to_chr.split_chr) {
        # Mutect
        call m.mutectNoFp as mutect {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            tumor_bam = index_bam.indexed_bam,
            tumor_bam_bai = index_bam.indexed_bam_bai,
            interval_list = chr_bed,
            scatter_count = scatter_count,
            tumor_only = tumor_only
        }
        #call bic.bcftoolsIsecComplement as mutect_isec_complement_gnomAD {
        #    input:
        #    vcf = mutect.unfiltered_vcf,
        #    vcf_tbi = mutect.unfiltered_vcf_tbi,
        #    exclude_vcf = normalized_gnomad_exclude,
        #    exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
        #    output_vcf_name = "mutect." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
        #    output_type = "z"
        #}
        call pp.pon2Percent as mutect_pon2 {
            input:
            #vcf = mutect_isec_complement_gnomAD.complement_vcf,
            vcf = mutect.unfiltered_vcf,
            vcf2PON = mutect_pon2_file,
            vcf2PON_tbi = mutect_pon2_file_tbi,
            caller = "mutect",
            sample_name = tumor_sample_name
        }

        # Vardict
        call v.vardictNoFp as vardict {
            input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = index_bam.indexed_bam,
            tumor_bam_bai = index_bam.indexed_bam_bai,
            interval_bed = chr_bed,
            bcbio_filter_string = bcbio_filter_string,
            tumor_sample_name = tumor_sample_name,
            min_var_freq = af_threshold,
            tumor_only = tumor_only
        }
        #call bic.bcftoolsIsecComplement as vardict_isec_complement_gnomAD {
        #    input:
        #    vcf = vardict.bcbio_filtered_vcf,
        #    vcf_tbi = vardict.bcbio_filtered_vcf_tbi,
        #    exclude_vcf = normalized_gnomad_exclude,
        #    exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
        #    output_vcf_name = "vardict." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
        #    output_type = "z"
        #}
        call pp.pon2Percent as vardict_pon2 {
            input:
            #vcf = vardict_isec_complement_gnomAD.complement_vcf,
            vcf = vardict.bcbio_filtered_vcf,
            vcf2PON = vardict_pon2_file,
            vcf2PON_tbi = vardict_pon2_file_tbi,
            caller = "vardict",
            sample_name = tumor_sample_name
        }

        # Lofreq
        call l.lofreqNoFp as lofreq {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            tumor_bam = index_bam.indexed_bam,
            tumor_bam_bai = index_bam.indexed_bam_bai,
            interval_bed = chr_bed,
            tumor_sample_name = tumor_sample_name,
            tumor_only = tumor_only
        }
        #call bic.bcftoolsIsecComplement as lofreq_isec_complement_gnomAD {
        #    input:
        #    vcf = lofreq.unfiltered_vcf,
        #    vcf_tbi = lofreq.unfiltered_vcf_tbi,
        #    exclude_vcf = normalized_gnomad_exclude,
        #    exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
        #    output_vcf_name = "lofreq." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
        #    output_type = "z"
        #}
        call pp.pon2Percent as lofreq_pon2 {
            input:
            #vcf = lofreq_isec_complement_gnomAD.complement_vcf,
            vcf = lofreq.unfiltered_vcf,
            vcf2PON = lofreq_pon2_file,
            vcf2PON_tbi = lofreq_pon2_file_tbi,
            caller = "lofreq",
            sample_name = tumor_sample_name
        }

        scatter (caller_vcf in [mutect_pon2.annotated_vcf, vardict_pon2.annotated_vcf, lofreq_pon2.annotated_vcf]){
            call cfv.createFakeVcf as fake_vcf {
                input:
                vcf = caller_vcf,
                tumor_sample_name = tumor_sample_name
            }
        }

        call mv.mergeVcf as mergeCallers {
            input:
            vcfs = fake_vcf.fake_vcf,
            vcf_tbis = fake_vcf.fake_vcf_tbi,
            merged_vcf_basename = "all_callers." + tumor_sample_name
        }

        call ffnn.fpFilterNoNorm as fpFilter {
            input:
            bam=index_bam.indexed_bam,
            bam_bai=index_bam.indexed_bam_bai,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            vcf = mergeCallers.merged_vcf,
            vcf_tbi = mergeCallers.merged_vcf_tbi,
            variant_caller="all_callers",
            sample_name=tumor_sample_name,
            min_var_freq=af_threshold
        }

        call pf.PoNFilter as PoN_filter {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            caller_vcf = mergeCallers.merged_vcf,
            normal_bams_file = pon_normal_bams_file,
            pon_bams = pon_bams,
            pon_final_name = "all_callers." + tumor_sample_name + ".pon.pileup",
            arrayMode = arrayMode
        }

        call vep.vep as vep {
            input:
            vcf = mergeCallers.merged_vcf,
            cache_dir_zip = vep_cache_dir_zip,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            plugins = vep_plugins,
            ensembl_assembly = vep_ensembl_assembly,
            ensembl_version = vep_ensembl_version,
            ensembl_species = vep_ensembl_species,
            synonyms_file = synonyms_file,
            custom_annotations = vep_custom_annotations,
            coding_only = annotate_coding_only,
            pick = vep_pick
        }

        call ac.annotateCaller as mutect_annotate_vcf {
            input:
            vcf = mutect_pon2.annotated_vcf,
            vcf_tbi = mutect_pon2.annotated_vcf_tbi,
            fp_filter = fpFilter.unfiltered_vcf,
            fp_filter_tbi = fpFilter.unfiltered_vcf_tbi,
            pileup_file = PoN_filter.pon_total_counts,
            pileup_file_tbi = PoN_filter.pon_total_counts_tbi,
            vep = vep.annotated_vcf,
            vep_tbi = vep.annotated_vcf_tbi,
            caller_prefix = "mutect",
            sample_name = tumor_sample_name,
            pon_pvalue = pon_pvalue
        }

        call ac.annotateCaller as vardict_annotate_vcf {
            input:
            vcf = vardict_pon2.annotated_vcf,
            vcf_tbi = vardict_pon2.annotated_vcf_tbi,
            fp_filter = fpFilter.unfiltered_vcf,
            fp_filter_tbi = fpFilter.unfiltered_vcf_tbi,
            pileup_file = PoN_filter.pon_total_counts,
            pileup_file_tbi = PoN_filter.pon_total_counts_tbi,
            vep = vep.annotated_vcf,
            vep_tbi = vep.annotated_vcf_tbi,
            caller_prefix = "vardict",
            sample_name = tumor_sample_name,
            pon_pvalue = pon_pvalue
        }

        call ac.annotateCaller as lofreq_annotate_vcf {
            input:
            vcf = lofreq_pon2.annotated_vcf,
            vcf_tbi = lofreq_pon2.annotated_vcf_tbi,
            fp_filter = fpFilter.unfiltered_vcf,
            fp_filter_tbi = fpFilter.unfiltered_vcf_tbi,
            pileup_file = PoN_filter.pon_total_counts,
            pileup_file_tbi = PoN_filter.pon_total_counts_tbi,
            vep = vep.annotated_vcf,
            vep_tbi = vep.annotated_vcf_tbi,
            caller_prefix = "lofreq",
            sample_name = tumor_sample_name,
            pon_pvalue = pon_pvalue
        }
    }

    call mv.mergeVcf as merge_mutect_full {
        input:
            vcfs = mutect.unfiltered_vcf,
            vcf_tbis = mutect.unfiltered_vcf_tbi,
            merged_vcf_basename = "mutect_full." + tumor_sample_name
    }
    call mv.mergeVcf as merge_mutect_pon {
        input:
            vcfs = mutect_annotate_vcf.pon_annotated_vcf,
            vcf_tbis = mutect_annotate_vcf.pon_annotated_vcf_tbi,
            merged_vcf_basename = "mutect." + tumor_sample_name + ".pileup.fisherPON"
    }
    call mv.mergeVcf as merge_mutect_final {
        input:
            vcfs = mutect_annotate_vcf.final_annotated_vcf,
            vcf_tbis = mutect_annotate_vcf.final_annotated_vcf_tbi,
            merged_vcf_basename = "mutect." + tumor_sample_name + ".final.annotated"
    }

    call mv.mergeVcf as merge_vardict_full {
        input:
            vcfs = vardict.unfiltered_vcf,
            vcf_tbis = vardict.unfiltered_vcf_tbi,
            merged_vcf_basename = "vardict_full." + tumor_sample_name
    }
    call mv.mergeVcf as merge_vardict_pon {
        input:
            vcfs = vardict_annotate_vcf.pon_annotated_vcf,
            vcf_tbis = vardict_annotate_vcf.pon_annotated_vcf_tbi,
            merged_vcf_basename = "vardict." + tumor_sample_name + ".pileup.fisherPON"
    }
    call mv.mergeVcf as merge_vardict_final {
        input:
            vcfs = vardict_annotate_vcf.final_annotated_vcf,
            vcf_tbis = vardict_annotate_vcf.final_annotated_vcf_tbi,
            merged_vcf_basename = "vardict." + tumor_sample_name + ".final.annotated"
    }

    call mv.mergeVcf as merge_lofreq_full {
        input:
            vcfs = lofreq.unfiltered_vcf,
            vcf_tbis = lofreq.unfiltered_vcf_tbi,
            merged_vcf_basename = "lofreq_full." + tumor_sample_name
    }
    call mv.mergeVcf as merge_lofreq_pon {
        input:
            vcfs = lofreq_annotate_vcf.pon_annotated_vcf,
            vcf_tbis = lofreq_annotate_vcf.pon_annotated_vcf_tbi,
            merged_vcf_basename = "lofreq." + tumor_sample_name + ".pileup.fisherPON"
    }
    call mv.mergeVcf as merge_lofreq_final {
        input:
            vcfs = lofreq_annotate_vcf.final_annotated_vcf,
            vcf_tbis = lofreq_annotate_vcf.final_annotated_vcf_tbi,
            merged_vcf_basename = "lofreq." + tumor_sample_name + ".final.annotated"
    }

    call mv.mergeVcf as merge_pon {
        input:
            vcfs = PoN_filter.pon_total_counts,
            vcf_tbis = PoN_filter.pon_total_counts_tbi,
            merged_vcf_basename = tumor_sample_name + ".pon.total.counts"
    }

    call mv.mergeVcf as merge_fp_filter {
        input:
            vcfs = fpFilter.unfiltered_vcf,
            vcf_tbis = fpFilter.unfiltered_vcf_tbi,
            merged_vcf_basename = tumor_sample_name + ".fpfilter"
    }

    call mv.mergeVcf as merge_vep {
        input:
            vcfs = vep.annotated_vcf,
            vcf_tbis = vep.annotated_vcf_tbi,
            merged_vcf_basename = tumor_sample_name + ".vep"
    }

    call ara.archerRAnnotate as annotateRMutect {
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
            cosmic_dir_zip = cosmic_dir_zip
    }

    call ara.archerRAnnotate as annotateRLofreq {
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
            cosmic_dir_zip = cosmic_dir_zip
    }

    call ara.archerRAnnotate as annotateRVardict {
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
            cosmic_dir_zip = cosmic_dir_zip
    }

    call xm.XGBModel as model {
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
        #File aligned_bam = alignment_workflow.aligned_bam
        File bqsr_bam = index_bam.indexed_bam

        # Tumor QC
        File tumor_insert_size_metrics = tumor_qc.insert_size_metrics
        File tumor_alignment_summary_metrics = tumor_qc.alignment_summary_metrics
        File tumor_hs_metrics = tumor_qc.hs_metrics
        Array[File] tumor_per_target_coverage_metrics = tumor_qc.per_target_coverage_metrics
        Array[File] tumor_per_target_hs_metrics = tumor_qc.per_target_hs_metrics
        Array[File] tumor_per_base_coverage_metrics = tumor_qc.per_base_coverage_metrics
        Array[File] tumor_per_base_hs_metrics = tumor_qc.per_base_hs_metrics
        Array[File] tumor_summary_hs_metrics = tumor_qc.summary_hs_metrics
        File tumor_flagstats = tumor_qc.flagstats
        File tumor_verify_bam_id_metrics = tumor_qc.verify_bam_id_metrics
        File tumor_verify_bam_id_depth = tumor_qc.verify_bam_id_depth

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
