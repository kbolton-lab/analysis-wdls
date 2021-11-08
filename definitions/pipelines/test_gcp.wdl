version 1.0

import "../types.wdl"

import "../subworkflows/gnomad_and_PoN_filter.wdl" as gapf
import "../subworkflows/lofreq.wdl" as l

import "../tools/vep.wdl" as vep
import "../tools/pon2percent.wdl" as pp

workflow test_gcp {
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

        # Variant Calling
        Boolean? arrayMode = true       # Decide if you would rather use the File (--bam_fof) or the Array (--bam) does the same thing, just input type is different
        File pon_normal_bams            # on GCP, it's not possible to do File because the file paths are unaccessable for each VM instance, so you have to do ArrayMode
        Array[File] bams                # This is just an array of Files... if you have the bam_fof it's easier just to use above and set this to empty array []
        Array[File] bams_bai
        Boolean tumor_only = true
        # Boolean mutect_artifact_detection_mode = false
        # Float? mutect_max_alt_allele_in_normal_fraction
        # Int? mutect_max_alt_alleles_in_normal_count
        # Int pindel_insert_size = 400
        # String? ref_name = "GRCh38DH"
        # String? ref_date = "20161216"
        # Int? pindel_min_supporting_reads = 3
        Float? af_threshold = 0.0001
        # String bcbio_filter_string = "((FMT/AF * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 2.0) || (FMT/DP < 10) || (INFO/QUAL < 45)))"
        String? pon_pvalue = "4.098606e-08"

        # PoN2
        /* File mutect_pon2_file
        File mutect_pon2_file_tbi */
        File lofreq_pon2_file
        File lofreq_pon2_file_tbi
        /* File vardict_pon2_file
        File vardict_pon2_file_tbi
        File pindel_pon2_file
        File pindel_pon2_file_tbi */

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

    call l.lofreq as lofreq {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            tumor_bam = bam,
            tumor_bam_bai = bam_bai,
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
        bams_bai = bams_bai,
        pon_final_name = "lofreq." + tumor_sample_name + ".pon.pileup",
        pon_pvalue = pon_pvalue
    }
    call vep.vep as lofreq_annotate_variants {
        input:
        vcf = lofreq_gnomad_pon_filters.processed_filtered_vcf,
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
    call pp.pon2Percent as lofreq_pon2 {
        input:
        vcf = lofreq_annotate_variants.annotated_vcf,
        vcf2PON = lofreq_pon2_file,
        vcf2PON_tbi = lofreq_pon2_file_tbi,
        caller = "lofreq",
        sample_name = tumor_sample_name
    }

    output {
        /* # Alignments
        File aligned_bam = alignment_workflow.aligned_bam
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
        File mutect_full = mutect.unfiltered_vcf                                                                # Raw Mutect Ouput
        File mutect_pon_annotated_unfiltered_vcf = mutect_gnomad_pon_filters.processed_gnomAD_filtered_vcf      # gnomAD Filtered w/ PoN Annotated
        File mutect_pon_annotated_filtered_vcf = mutect_pon2.annotated_vcf                                      # gnomAD Filtered + PoN Filtered + PoN2 Filtered w/ VEP Annotation
        File mutect_pon_total_counts = mutect_gnomad_pon_filters.pon_total_counts                               # PoN Pileup Results */

        # Lofreq
        File lofreq_full = lofreq.unfiltered_vcf                                                                # Raw Lofreq Ouput
        File lofreq_pon_annotated_unfiltered_vcf = lofreq_gnomad_pon_filters.processed_gnomAD_filtered_vcf      # gnomAD Filtered w/ PoN Annotated
        File lofreq_pon_annotated_filtered_vcf = lofreq_pon2.annotated_vcf                                      # gnomAD Filtered + PoN Filtered + PoN2 Filtered w/ VEP Annotation
        File lofreq_pon_total_counts = lofreq_gnomad_pon_filters.pon_total_counts                               # PoN Pileup Results

        /* # Vardict
        File vardict_full = vardict.unfiltered_vcf                                                                # Raw Vardict Ouput
        File vardict_pon_annotated_unfiltered_vcf = vardict_gnomad_pon_filters.processed_gnomAD_filtered_vcf      # gnomAD Filtered w/ PoN Annotated
        File vardict_pon_annotated_filtered_vcf = vardict_pon2.annotated_vcf                                      # gnomAD Filtered + PoN Filtered + PoN2 Filtered w/ VEP Annotation
        File vardict_pon_total_counts = vardict_gnomad_pon_filters.pon_total_counts                               # PoN Pileup Results

        # Pindel
        File pindel_full = pindel.unfiltered_vcf                                                                # Raw Pindel Ouput
        File pindel_pon_annotated_unfiltered_vcf = pindel_gnomad_pon_filters.processed_gnomAD_filtered_vcf      # gnomAD Filtered w/ PoN Annotated
        File pindel_pon_annotated_filtered_vcf = pindel_pon2.annotated_vcf                                      # gnomAD Filtered + PoN Filtered + PoN2 Filtered w/ VEP Annotation
        File pindel_pon_total_counts = pindel_gnomad_pon_filters.pon_total_counts                               # PoN Pileup Results */

    }
}
