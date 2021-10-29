version 1.0

import "../types.wdl"

import "../subworkflows/archer_fastq_format.wdl" as fqf
import "../subworkflows/molecular_alignment.wdl" as ma
import "../subworkflows/qc_exome.wdl" as qe
import "../subworkflows/gnomad_and_PoN_filter.wdl" as gapf
import "../subworkflows/mutect.wdl" as m
import "../subworkflows/lofreq.wdl" as l
import "../subworkflows/vardict.wdl" as v
import "../subworkflows/pindel.wdl" as p

import "../tools/fastq_to_bam.wdl" as ftb
import "../tools/bam_to_cram.wdl" as btc
import "../tools/index_cram.wdl" as ic
import "../tools/bqsr.wdl" as b
import "../tools/apply_bqsr.wdl" as ab
import "../tools/index_bam.wdl" as ib
import "../tools/interval_list_expand.wdl" as ile
import "../tools/vep.wdl" as vep
import "../tools/pon2percent.wdl" as pp

workflow archerdx {
    input {
        # Pipeline
        Int scatter_count = 20

        # Sequence and BAM Information
        # Array[SequenceData] sequence
        # Array[String] read_structure
        File tumor_bam
        File tumor_bai
        File? normal_bam
        File? normal_bai
        String? tumor_name = "tumor"
        String tumor_sample_name
        File? target_intervals

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
        Int? umi_length = 8
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
        Int target_interval_padding = 100
        Array[LabelledFile] per_base_intervals
        Array[LabelledFile] per_target_intervals
        Array[LabelledFile] summary_intervals
        File omni_vcf
        File omni_vcf_tbi
        String picard_metric_accumulation_level
        Int? qc_minimum_mapping_quality = 0
        Int? qc_minimum_base_quality = 0

        # Variant Calling
        Boolean? arrayMode = false       # Decide if you would raather use the File (--bam_fof) or the Array (--bam) does the same thing, just input type is different
        File pon_normal_bams
        Array[String] bams              # This is just an array of Files (as Strings)... if you have the bam_fof it's easier just to use above and set this to empty array []
        Boolean tumor_only = true
        Boolean mutect_artifact_detection_mode = false
        Float? mutect_max_alt_allele_in_normal_fraction
        Int? mutect_max_alt_alleles_in_normal_count
        Int pindel_insert_size = 400
        String? ref_name = "GRCh38DH"
        String? ref_date = "20161216"
        Int? pindel_min_supporting_reads = 3
        Float? af_threshold = 0.005
        String bcbio_filter_string = "((FMT/AF * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 2.0) || (FMT/DP < 10) || (INFO/QUAL < 45)))"
        String? pon_pvalue = "4.098606e-08"

        # PoN2
        File mutect_pon2_file
        File mutect_pon2_file_tbi
        File lofreq_pon2_file
        File lofreq_pon2_file_tbi
        File vardict_pon2_file
        File vardict_pon2_file_tbi
        File pindel_pon2_file
        File pindel_pon2_file_tbi

        # TODO: R Stuff
        # File impact_annotation
        # File topmed_annotation
        # File cosmic_annotation
        # File tsg_annotation
        # File oncoKB_annotation
        # File pd_table_annotation
        # File panmyeloid_annotation
        # File blacklist_annotation
        # File segemental_duplications_annotation
        # File simple_repeats_annotation
        # File repeat_masker_annotation

        # VEP
        File vep_cache_dir_zip
        String vep_ensembl_assembly
        String vep_ensembl_version
        String vep_ensembl_species
        Array[String] vep_plugins = ["Frameshift", "Wildtype"]
        File? synonyms_file
        Boolean? annotate_coding_only = true
        Array[VepCustomAnnotation] vep_custom_annotations
        File limit_variant_intervals
        Array[String] variants_to_table_fields = ["CHROM","POS","ID","REF","ALT","set","AC","AF"]
        Array[String]? variants_to_table_genotype_fields = ["GT","AD"]
        Array[String]? vep_to_table_fields = ["HGVSc","HGVSp"]
        String vep_pick = "flag_pick"

        # gnomAD
        Float filter_gnomADe_maximum_population_allele_frequency = 0.005
        File normalized_gnomad_exclude
        File normalized_gnomad_exclude_tbi
        String filter_flag = "include"
    }
    

    call b.bqsr as bqsr {
        input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        bam = alignment_workflow.aligned_bam,
        bam_bai = alignment_workflow.aligned_bam_bai,
        intervals = bqsr_intervals,
        known_sites = bqsr_known_sites,
        known_sites_tbi = bqsr_known_sites_tbi
    }

    call ab.applyBqsr as apply_bqsr {
        input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        bam = alignment_workflow.aligned_bam,
        bam_bai = alignment_workflow.aligned_bam_bai,
        bqsr_table = bqsr.bqsr_table
    }

    call ib.indexBam as index_bam {
        input:
        bam = apply_bqsr.bqsr_bam
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

    call ile.intervalListExpand as pad_target_intervals {
        input:
        interval_list = target_intervals,
        roi_padding = target_interval_padding
    }

    call m.mutect as mutect {
        input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        tumor_bam = index_bam.indexed_bam,
        tumor_bam_bai = index_bam.indexed_bam_bai,
        interval_list = target_intervals,
        scatter_count = scatter_count,
        tumor_sample_name = tumor_sample_name,
        min_var_freq = af_threshold,
        tumor_only = tumor_only
    }
    call gapf.gnomadAndPoNFilter as mutect_gnomad_pon_filters {
        input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        caller_vcf = mutect.unfiltered_vcf,
        caller_vcf_tbi = mutect.unfiltered_vcf_tbi,
        gnomAD_exclude_vcf = normalized_gnomad_exclude,
        gnomAD_exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
        caller_prefix = "mutect." + tumor_sample_name,
        arrayMode = arrayMode,
        normal_bams = pon_normal_bams,
        bams = bams,
        pon_final_name = "mutect." + tumor_sample_name + ".pon.pileup",
        pon_pvalue = pon_pvalue
    }
    call vep.vep as mutect_annotate_variants {
        input:
        vcf = mutect_gnomad_pon_filters.processed_filtered_vcf,
        cache_dir_zip = cache_dir_zip,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        plugins = plugins,
        ensembl_assembly = ensembl_assembly,
        ensembl_version = ensembl_version,
        ensembl_species = ensembl_species,
        synonyms_file = synonyms_file,
        custom_annotations = custom_annotations,
        coding_only = annotate_coding_only,
        pick = vep_pick
    }
    call pp.pon2Percent as mutect_pon2 {
        input:
        vcf = mutect_annotate_variants.annotated_vcf,
        vcf2PON = mutect_pon2_file,
        vcf2PON_tbi = mutect_pon2_file_tbi,
        caller = "mutect",
        sample_name = tumor_sample_name
    }

    # Vardict
    call v.vardict as vardict {
        input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        tumor_bam = index_bam.indexed_bam,
        tumor_bam_bai = index_bam.indexed_bam_bai,
        interval_list = target_intervals,
        scatter_count = scatter_count,
        tumor_sample_name = tumor_sample_name,
        min_var_freq = af_threshold,
        bcbio_filter_string = bcbio_filter_string,
        tumor_only = tumor_only
    }
    call gapf.gnomadAndPoNFilter as vardict_gnomad_pon_filters {
        input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        caller_vcf = vardict.unfiltered_vcf,
        caller_vcf_tbi = vardict.unfiltered_vcf_tbi,
        gnomAD_exclude_vcf = normalized_gnomad_exclude,
        gnomAD_exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
        caller_prefix = "vardict." + tumor_sample_name,
        arrayMode = arrayMode,
        normal_bams = pon_normal_bams,
        bams = bams,
        pon_final_name = "vardict." + tumor_sample_name + ".pon.pileup",
        pon_pvalue = pon_pvalue
    }
    call vep.vep as vardict_annotate_variants {
        input:
        vcf = vardict_gnomad_pon_filters.processed_filtered_vcf,
        cache_dir_zip = cache_dir_zip,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        plugins = plugins,
        ensembl_assembly = ensembl_assembly,
        ensembl_version = ensembl_version,
        ensembl_species = ensembl_species,
        synonyms_file = synonyms_file,
        custom_annotations = custom_annotations,
        coding_only = annotate_coding_only,
        pick = vep_pick
    }
    call pp.pon2Percent as vardict_pon2 {
        input:
        vcf = vardict_annotate_variants.annotated_vcf,
        vcf2PON = vardict_pon2_file,
        vcf2PON_tbi = vardict_pon2_file_tbi,
        caller = "vardict",
        sample_name = tumor_sample_name
    }

    call l.lofreq as lofreq {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            tumor_bam = index_bam.indexed_bam,
            tumor_bam_bai = index_bam.indexed_bam_bai,
            interval_list = target_intervals,
            scatter_count = scatter_count,
            tumor_sample_name = tumor_sample_name,
            min_var_freq = af_threshold,
            tumor_only = tumor_only
    }
    call gapf.gnomadAndPoNFilter as lofreq_gnomad_pon_filters {
        input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        caller_vcf = lofreq.unfiltered_vcf,
        caller_vcf_tbi = lofreq.unfiltered_vcf_tbi,
        gnomAD_exclude_vcf = normalized_gnomad_exclude,
        gnomAD_exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
        caller_prefix = "lofreq." + tumor_sample_name,
        arrayMode = arrayMode,
        normal_bams = pon_normal_bams,
        bams = bams,
        pon_final_name = "lofreq." + tumor_sample_name + ".pon.pileup",
        pon_pvalue = pon_pvalue
    }
    call vep.vep as lofreq_annotate_variants {
        input:
        vcf = lofreq_gnomad_pon_filters.processed_filtered_vcf,
        cache_dir_zip = cache_dir_zip,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        plugins = plugins,
        ensembl_assembly = ensembl_assembly,
        ensembl_version = ensembl_version,
        ensembl_species = ensembl_species,
        synonyms_file = synonyms_file,
        custom_annotations = custom_annotations,
        coding_only = annotate_coding_only,
        pick = vep_pick
    }
    call pp.pon2Percent as lofreq_pon2 {
        input:
        vcf = lofreq_annotate_variants.annotated_vcf,
        vcf2PON = lofreq_pon2_file,
        vcf2PON_tbi = lofreq_pon2_file_tbi,
        caller = "lofreq",
        sample_name = tumor_sample_name
    }

    call p.pindel as pindel {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            tumor_bam = index_bam.indexed_bam,
            tumor_bam_bai = index_bam.indexed_bam_bai,
            interval_list = target_intervals,
            scatter_count = scatter_count,
            insert_size = pindel_insert_size,
            tumor_sample_name = tumor_sample_name,
            ref_name = ref_name,
            ref_date = ref_date,
            pindel_min_supporting_reads = pindel_min_supporting_reads,
            min_var_freq = af_threshold,
            tumor_only = tumor_only
    }

    output {
        # Alignments
        File aligned_bam = alignment_workflow.aligned_bam
        File bqsr_bam = index_bam/indexed_bam

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
        File mutect_full = mutect.unfiltered_vcf                                                                # Raw Mutect Ouput
        File mutect_pon_annotated_unfiltered_vcf = mutect_gnomad_pon_filters.processed_gnomAD_filtered_vcf      # gnomAD Filtered w/ PoN Annotated
        File mutect_pon_annotated_filtered_vcf = mutect_pon2.annotated_vcf                                      # gnomAD Filtered + PoN Filtered + PoN2 Filtered w/ VEP Annotation
        File mutect_pon_total_counts = mutect_gnomad_pon_filters.pon_total_counts                               # PoN Pileup Results

        # Lofreq
        #File lofreq_full = lofreq.unfiltered_vcf                                                                # Raw Lofreq Ouput
        #File lofreq_pon_annotated_unfiltered_vcf = lofreq_gnomad_pon_filters.processed_gnomAD_filtered_vcf      # gnomAD Filtered w/ PoN Annotated
        #File lofreq_pon_annotated_filtered_vcf = lofreq_pon2.annotated_vcf                                      # gnomAD Filtered + PoN Filtered + PoN2 Filtered w/ VEP Annotation
        #File lofreq_pon_total_counts = lofreq_gnomad_pon_filters.pon_total_counts                               # PoN Pileup Results

        # Vardict
        #File vardict_full = vardict.unfiltered_vcf                                                                # Raw Vardict Ouput
        #File vardict_pon_annotated_unfiltered_vcf = vardict_gnomad_pon_filters.processed_gnomAD_filtered_vcf      # gnomAD Filtered w/ PoN Annotated
        #File vardict_pon_annotated_filtered_vcf = vardict_pon2.annotated_vcf                                      # gnomAD Filtered + PoN Filtered + PoN2 Filtered w/ VEP Annotation
        #File vardict_pon_total_counts = vardict_gnomad_pon_filters.pon_total_counts                               # PoN Pileup Results

        # Pindel
        #File pindel_full = pindel.unfiltered_vcf                                                                # Raw Pindel Ouput
        #File pindel_pon_annotated_unfiltered_vcf = pindel_gnomad_pon_filters.processed_gnomAD_filtered_vcf      # gnomAD Filtered w/ PoN Annotated
        #File pindel_pon_annotated_filtered_vcf = pindel_pon2.annotated_vcf                                      # gnomAD Filtered + PoN Filtered + PoN2 Filtered w/ VEP Annotation
        #File pindel_pon_total_counts = pindel_gnomad_pon_filters.pon_total_counts                               # PoN Pileup Results

        #File gnomAD_exclude = get_gnomad_exclude.normalized_gnomad_exclude

        # TODO: Maybe Implement R Things
        # File final_tsv = final_annotation.final_tsv
        # File column_check = final_annotation.column_check
    }
}

    # final_annotation:
    #     run: ../tools/annotate_CH_pd.cwl
    #     in:
    #         mutect_vcf: mutect_pon2/annotated_vcf
    #         varscan_vcf: varscan_pon2/annotated_vcf
    #         vardict_vcf: vardict_pon2/annotated_vcf
    #         pindel_vcf: pindel_pon2/annotated_vcf
    #         sample_name: tumor_sample_name
    #         impact_annotation: impact_annotation
    #         topmed_annotation: topmed_annotation
    #         cosmic_annotation: cosmic_annotation
    #         tsg_annotation: tsg_annotation
    #         oncoKB_annotation: oncoKB_annotation
    #         pd_table_annotation: pd_table_annotation
    #         panmyeloid_annotation: panmyeloid_annotation
    #         blacklist_annotation: blacklist_annotation
    #         segemental_duplications_annotation: segemental_duplications_annotation
    #         simple_repeats_annotation: simple_repeats_annotation
    #         repeat_masker_annotation: repeat_masker_annotation
    #     out:
    #         [final_tsv, column_check]
