version 1.0

import "../types.wdl"

import "../subworkflows/qc_exome.wdl" as qe
import "../subworkflows/gnomad_and_MAPQ0_filter.wdl" as gamf
import "../subworkflows/mutect.wdl" as m
import "../subworkflows/lofreq.wdl" as l
import "../subworkflows/vardict.wdl" as v
import "../subworkflows/varscan.wdl" as vs
import "../subworkflows/pindel.wdl" as p
import "../tools/bqsr.wdl" as b
import "../tools/apply_bqsr.wdl" as ab
import "../tools/index_bam.wdl" as ib
import "../tools/interval_list_expand.wdl" as ile
import "../tools/vep_brian.wdl" as vep

workflow somatic_tumor_only {
    input {
        # Pipeline
        Int scatter_count = 10

        # Sequence and BAM Information
        # Array[SequenceData] sequence
        # Array[String] read_structure
        File tumor_bam
        File tumor_bai
        File? normal_bam
        File? normal_bai
        String? tumor_name = "tumor"
        String tumor_sample_name
        File target_intervals
        File target_bed

        # Reference
        File reference
        File reference_fai
        File reference_dict
        File? reference_amb
        File? reference_ann
        File? reference_bwt
        File? reference_pac
        File? reference_sa

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
        Array[String]? bqsr_intervals

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
        File? pon_normal_bams
        Array[String]? bams              # This is just an array of Files (as Strings)... if you have the bam_fof it's easier just to use above and set this to empty array []
        Boolean tumor_only = true
        Boolean mutect_artifact_detection_mode = false
        Float? mutect_max_alt_allele_in_normal_fraction
        Int? mutect_max_alt_alleles_in_normal_count
        Int pindel_insert_size = 400
        String? ref_name = "GRCh38DH"
        String? ref_date = "20161216"
        Int? pindel_min_supporting_reads = 3
        Float? af_threshold = 0.005
        String bcbio_filter_string = "((FMT/AF * FMT/DP < 500) && ( INFO/MQ < 55.0 || FMT/DP < 800 || FMT/QUAL < 30 ))"
        String? pon_pvalue = "4.098606e-08"

        # PoN2
        File? mutect_pon2_file
        File? mutect_pon2_file_tbi
        File? lofreq_pon2_file
        File? lofreq_pon2_file_tbi
        File? vardict_pon2_file
        File? vardict_pon2_file_tbi
        File? pindel_pon2_file
        File? pindel_pon2_file_tbi

        

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
        String additional_args = "--sift p --polyphen p --pick --pick_order canonical,rank,mane,ccds,appris,tsl,biotype,length --everything 1 --merged --check_existing --buffer_size 1000 --af_gnomad"

        # gnomAD
        Float filter_gnomADe_maximum_population_allele_frequency = 0.005
        File? normalized_gnomad_file # option for mutect2
        File? normalized_gnomad_file_tbi # option for mutect2
        File normalized_gnomad_exclude # gnomad filter
        File normalized_gnomad_exclude_tbi # gnomad filter
        String filter_flag = "include"
        File? pon # option for mutect2
        File? pon_tbi # option for mutect2
    }
    

    call b.bqsr as bqsr {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bam = tumor_bam,
            bam_bai = tumor_bai,
            intervals = bqsr_intervals,
            known_sites = bqsr_known_sites,
            known_sites_tbi = bqsr_known_sites_tbi
    }

    call ab.applyBqsr as apply_bqsr {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bam = tumor_bam,
            bam_bai = tumor_bai,
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

    # Mutect PoN
    call m.mutect as mutect_pon {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            tumor_bam = index_bam.indexed_bam,
            tumor_bam_bai = index_bam.indexed_bam_bai,
            interval_list = target_intervals,
            pon = pon,
            pon_tbi = pon_tbi,
            gnomad_file = normalized_gnomad_file,
            gnomad_file_tbi = normalized_gnomad_file_tbi,
            scatter_count = scatter_count,
            tumor_sample_name = tumor_sample_name,
            min_var_freq = af_threshold,
            tumor_only = tumor_only,
            variant_caller = "mutect_pon"
    }

    call gamf.gnomadAndMapQ0Filter as mutect_pon_gnomad_mapq0_filters {
        input:
            caller_vcf = mutect_pon.unfiltered_vcf,
            caller_vcf_tbi = mutect_pon.unfiltered_vcf_tbi,
            bam = tumor_bam,
            bam_bai = tumor_bai,
            gnomAD_exclude_vcf = normalized_gnomad_exclude,
            gnomAD_exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            caller_prefix = "mutect_pon." + tumor_sample_name
    }

    call vep.vep as mutect_pon_annotate_variants {
        input:
            vcf = mutect_pon_gnomad_mapq0_filters.mapq0_soft_filtered_vcf,
            cache_dir_zip = vep_cache_dir_zip,
            reference = reference,
            reference_fai = reference_fai,
            plugins = vep_plugins,
            ensembl_assembly = vep_ensembl_assembly,
            ensembl_version = vep_ensembl_version,
            ensembl_species = vep_ensembl_species,
            synonyms_file = synonyms_file,
            custom_annotations = vep_custom_annotations,
            coding_only = annotate_coding_only,
            additional_args = additional_args
    }

     # Mutect no PoN
    call m.mutect as mutect {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            tumor_bam = index_bam.indexed_bam,
            tumor_bam_bai = index_bam.indexed_bam_bai,
            interval_list = target_intervals,
            gnomad_file = normalized_gnomad_file,
            gnomad_file_tbi = normalized_gnomad_file_tbi,
            scatter_count = scatter_count,
            tumor_sample_name = tumor_sample_name,
            min_var_freq = af_threshold,
            tumor_only = tumor_only,
            variant_caller = "mutect"
    }

    call gamf.gnomadAndMapQ0Filter as mutect_gnomad_mapq0_filters {
        input:
            caller_vcf = mutect.unfiltered_vcf,
            caller_vcf_tbi = mutect.unfiltered_vcf_tbi,
            bam = tumor_bam,
            bam_bai = tumor_bai,
            gnomAD_exclude_vcf = normalized_gnomad_exclude,
            gnomAD_exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            caller_prefix = "mutect." + tumor_sample_name
    }

    call vep.vep as mutect_annotate_variants {
        input:
            vcf = mutect_gnomad_mapq0_filters.mapq0_soft_filtered_vcf,
            cache_dir_zip = vep_cache_dir_zip,
            reference = reference,
            reference_fai = reference_fai,
            plugins = vep_plugins,
            ensembl_assembly = vep_ensembl_assembly,
            ensembl_version = vep_ensembl_version,
            ensembl_species = vep_ensembl_species,
            synonyms_file = synonyms_file,
            custom_annotations = vep_custom_annotations,
            coding_only = annotate_coding_only,
            additional_args = additional_args
    }


    # Vardict
    call v.vardict as vardict {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            tumor_bam = index_bam.indexed_bam,
            tumor_bam_bai = index_bam.indexed_bam_bai,
            interval_bed = target_bed,
            tumor_sample_name = tumor_sample_name,
            min_var_freq = af_threshold,
            bcbio_filter_string = bcbio_filter_string,
            tumor_only = tumor_only
    }

    call gamf.gnomadAndMapQ0Filter as vardict_gnomad_mapq0_filters {
        input:
            caller_vcf = vardict.bcbio_filtered_vcf,
            caller_vcf_tbi = vardict.bcbio_filtered_vcf_tbi,
            bam = tumor_bam,
            bam_bai = tumor_bai,
            gnomAD_exclude_vcf = normalized_gnomad_exclude,
            gnomAD_exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            caller_prefix = "vardict." + tumor_sample_name
    }
    
    call vep.vep as vardict_annotate_variants {
        input:
            vcf = vardict_gnomad_mapq0_filters.mapq0_soft_filtered_vcf,
            cache_dir_zip = vep_cache_dir_zip,
            reference = reference,
            reference_fai = reference_fai,
            plugins = vep_plugins,
            ensembl_assembly = vep_ensembl_assembly,
            ensembl_version = vep_ensembl_version,
            ensembl_species = vep_ensembl_species,
            synonyms_file = synonyms_file,
            custom_annotations = vep_custom_annotations,
            additional_args = additional_args
    }


    call l.lofreq as lofreq {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            tumor_bam = index_bam.indexed_bam,
            tumor_bam_bai = index_bam.indexed_bam_bai,
            interval_bed = target_bed,
            tumor_sample_name = tumor_sample_name,
            min_var_freq = af_threshold,
            tumor_only = tumor_only
    }

    call gamf.gnomadAndMapQ0Filter as lofreq_gnomad_mapq0_filters {
        input:
            caller_vcf = lofreq.unfiltered_vcf,
            caller_vcf_tbi = lofreq.unfiltered_vcf_tbi,
            bam = tumor_bam,
            bam_bai = tumor_bai,
            gnomAD_exclude_vcf = normalized_gnomad_exclude,
            gnomAD_exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            caller_prefix = "lofreq." + tumor_sample_name
    }

    call vep.vep as lofreq_annotate_variants {
        input:
            vcf = lofreq_gnomad_mapq0_filters.mapq0_soft_filtered_vcf,
            cache_dir_zip = vep_cache_dir_zip,
            reference = reference,
            reference_fai = reference_fai,
            plugins = vep_plugins,
            ensembl_assembly = vep_ensembl_assembly,
            ensembl_version = vep_ensembl_version,
            ensembl_species = vep_ensembl_species,
            synonyms_file = synonyms_file,
            custom_annotations = vep_custom_annotations,
            coding_only = annotate_coding_only,
            additional_args = additional_args
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

    call gamf.gnomadAndMapQ0Filter as pindel_gnomad_mapq0_filters {
        input:
            caller_vcf = pindel.unfiltered_vcf,
            caller_vcf_tbi = pindel.unfiltered_vcf_tbi,
            bam = tumor_bam,
            bam_bai = tumor_bai,
            gnomAD_exclude_vcf = normalized_gnomad_exclude,
            gnomAD_exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            caller_prefix = "pindel." + tumor_sample_name
    }

    call vep.vep as pindel_annotate_variants {
        input:
            vcf = pindel_gnomad_mapq0_filters.mapq0_soft_filtered_vcf,
            cache_dir_zip = vep_cache_dir_zip,
            reference = reference,
            reference_fai = reference_fai,
            plugins = vep_plugins,
            ensembl_assembly = vep_ensembl_assembly,
            ensembl_version = vep_ensembl_version,
            ensembl_species = vep_ensembl_species,
            synonyms_file = synonyms_file,
            custom_annotations = vep_custom_annotations,
            additional_args = additional_args
    }


   # Varscan2
    call vs.varscan as varscan {
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

    call gamf.gnomadAndMapQ0Filter as varscan_gnomad_mapq0_filters {
        input:
            caller_vcf = varscan.unfiltered_vcf,
            caller_vcf_tbi = varscan.unfiltered_vcf_tbi,
            bam = tumor_bam,
            bam_bai = tumor_bai,
            gnomAD_exclude_vcf = normalized_gnomad_exclude,
            gnomAD_exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            caller_prefix = "varscan." + tumor_sample_name
    }

    call vep.vep as varscan_annotate_variants {
        input:
            vcf = varscan_gnomad_mapq0_filters.mapq0_soft_filtered_vcf,
            cache_dir_zip = vep_cache_dir_zip,
            reference = reference,
            reference_fai = reference_fai,
            plugins = vep_plugins,
            ensembl_assembly = vep_ensembl_assembly,
            ensembl_version = vep_ensembl_version,
            ensembl_species = vep_ensembl_species,
            synonyms_file = synonyms_file,
            custom_annotations = vep_custom_annotations,
            coding_only = annotate_coding_only,
            additional_args = additional_args
    }


    output {
        # Alignments
        
        File bqsr_bam = index_bam.indexed_bam
        File bqsr_bam_bai = index_bam.indexed_bam_bai

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

        # Mutect PoN
        File mutect_pon_full = mutect_pon.unfiltered_vcf # mutect_full.vcf.gz
        File mutect_pon_full_tbi = mutect_pon.unfiltered_vcf_tbi # mutect_full.vcf.gz.tbi
        File mutect_pon_gnomad_mapq = mutect_pon_gnomad_mapq0_filters.mapq0_soft_filtered_vcf
        File mutect_pon_gnomad_mapq_tbi = mutect_pon_gnomad_mapq0_filters.mapq0_soft_filtered_vcf_tbi
        File mutect_pon_vep = mutect_pon_annotate_variants.annotated_vcf
        File mutect_pon_vep_tbi = mutect_pon_annotate_variants.annotated_vcf_tbi

        # Mutect no PoN
        File mutect_full = mutect.unfiltered_vcf # mutect_full.vcf.gz
        File mutect_full_tbi = mutect.unfiltered_vcf_tbi # mutect_full.vcf.gz.tbi
        File mutect_gnomad_mapq = mutect_gnomad_mapq0_filters.mapq0_soft_filtered_vcf
        File mutect_gnomad_mapq_tbi = mutect_gnomad_mapq0_filters.mapq0_soft_filtered_vcf_tbi
        File mutect_vep = mutect_annotate_variants.annotated_vcf
        File mutect_vep_tbi = mutect_annotate_variants.annotated_vcf_tbi

        # Vardict
        File vardict_full = vardict.unfiltered_vcf
        File vardict_full_tbi = vardict.unfiltered_vcf_tbi
        File vardict_gnomad_mapq = vardict_gnomad_mapq0_filters.mapq0_soft_filtered_vcf
        File vardict_gnomad_mapq_tbi = vardict_gnomad_mapq0_filters.mapq0_soft_filtered_vcf_tbi
        File vardict_vep = vardict_annotate_variants.annotated_vcf
        File vardict_vep_tbi = vardict_annotate_variants.annotated_vcf_tbi

        # Lofreq
        File lofreq_full = lofreq.unfiltered_vcf
        File lofreq_full_tbi = lofreq.unfiltered_vcf_tbi
        File lofreq_gnomad_mapq = lofreq_gnomad_mapq0_filters.mapq0_soft_filtered_vcf
        File lofreq_gnomad_mapq_tbi = lofreq_gnomad_mapq0_filters.mapq0_soft_filtered_vcf_tbi
        File lofreq_vep = lofreq_annotate_variants.annotated_vcf
        File lofreq_vep_tbi = lofreq_annotate_variants.annotated_vcf_tbi

        # Pindel
        File pindel_full = pindel.unfiltered_vcf
        File pindel_full_tbi = pindel.unfiltered_vcf_tbi
        File pindel_gnomad_mapq = pindel_gnomad_mapq0_filters.mapq0_soft_filtered_vcf
        File pindel_gnomad_mapq_tbi = pindel_gnomad_mapq0_filters.mapq0_soft_filtered_vcf_tbi
        File pindel_vep = pindel_annotate_variants.annotated_vcf
        File pindel_vep_tbi = pindel_annotate_variants.annotated_vcf_tbi

        # Varscan
        File varscan_full = varscan.unfiltered_vcf
        File varscan_full_tbi = varscan.unfiltered_vcf_tbi
        File varscan_gnomad_mapq = varscan_gnomad_mapq0_filters.mapq0_soft_filtered_vcf
        File varscan_gnomad_mapq_tbi = varscan_gnomad_mapq0_filters.mapq0_soft_filtered_vcf_tbi
        File varscan_vep = varscan_annotate_variants.annotated_vcf
        File varscan_vep_tbi = varscan_annotate_variants.annotated_vcf_tbi
    }
}
