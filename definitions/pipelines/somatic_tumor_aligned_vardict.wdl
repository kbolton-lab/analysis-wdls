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
    
    # Vardict
    call v.vardict as vardict {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            tumor_bam = tumor_bam,
            tumor_bam_bai = tumor_bai,
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
            reference_dict = reference_dict,
            plugins = vep_plugins,
            ensembl_assembly = vep_ensembl_assembly,
            ensembl_version = vep_ensembl_version,
            ensembl_species = vep_ensembl_species,
            synonyms_file = synonyms_file,
            custom_annotations = vep_custom_annotations,
            additional_args = additional_args
    }


    output {
        # Vardict
        File vardict_full = vardict.unfiltered_vcf
        File vardict_full_tbi = vardict.unfiltered_vcf_tbi
        File vardict_gnomad_mapq = vardict_gnomad_mapq0_filters.mapq0_soft_filtered_vcf
        File vardict_gnomad_mapq_tbi = vardict_gnomad_mapq0_filters.mapq0_soft_filtered_vcf_tbi
        File vardict_vep = vardict_annotate_variants.annotated_vcf
        File vardict_vep_tbi = vardict_annotate_variants.annotated_vcf_tbi
    }
}
