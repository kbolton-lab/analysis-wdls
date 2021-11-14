version 1.0

import "../types.wdl"

import "../subworkflows/archer_fastq_format.wdl" as fqf
import "../subworkflows/molecular_alignment.wdl" as ma
import "../subworkflows/qc_exome.wdl" as qe
import "../subworkflows/PoN_filter.wdl" as gapf
import "../subworkflows/fp_filter.wdl" as ff
import "../subworkflows/mutect_noFp.wdl" as m
import "../subworkflows/lofreq_noFp.wdl" as l
import "../subworkflows/vardict_noFp.wdl" as v

import "../tools/fastq_to_bam.wdl" as ftb
import "../tools/bqsr_apply.wdl" as ba
import "../tools/index_bam.wdl" as ib
import "../tools/bcftools_isec_complement.wdl" as bic
import "../tools/vep.wdl" as vep
import "../tools/pon2percent.wdl" as pp
import "../tools/split_bam_into_chr.wdl" as sbic
import "../tools/merge_vcf.wdl" as mv

workflow archerdx {
    input {
        # Pipeline
        Int scatter_count = 20

        # Sequence and BAM Information
        Array[SequenceData] sequence
        Array[String] read_structure
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
        Int pindel_insert_size = 400
        String? ref_name = "GRCh38DH"
        String? ref_date = "20161216"
        Int? pindel_min_supporting_reads = 3

        String bcbio_filter_string = "((FMT/AF * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 3.0) || (FMT/DP < 6500) || (INFO/QUAL < 27)))"

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

    scatter(seq_data in sequence) {
        call fqf.archerFastqFormat as format_fastq {
            input:
            sequence = seq_data,
            umi_length = umi_length
        }

        call ftb.fastqToBam as fastq_to_bam {
            input:
            fastq1 = format_fastq.fastq1,
            fastq2 = format_fastq.fastq2,
            sample_name = tumor_sample_name,
            library_name = "Library",
            platform_unit = "Illumina",
            platform = "ArcherDX"
        }
    }

    call ma.molecularAlignment as alignment_workflow {
        input:
        bam = fastq_to_bam.bam,
        sample_name = tumor_sample_name,
        read_structure = read_structure,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        reference_amb = reference_amb,
        reference_ann = reference_ann,
        reference_bwt = reference_bwt,
        reference_pac = reference_pac,
        reference_sa = reference_sa,
        target_intervals = target_intervals,
        umi_paired = umi_paired,
        min_reads = min_reads,
        max_read_error_rate = max_read_error_rate,
        max_base_error_rate = max_base_error_rate,
        min_base_quality = min_base_quality,
        max_no_call_fraction = max_no_call_fraction
    }

    call ba.bqsrApply as bqsr {
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

    call sbic.splitBamIntoChr as split_bam_into_chr {
        input:
        bam = index_bam.indexed_bam,
        interval_bed = target_bed
    }

    scatter (bam_chr in split_bam_into_chr.split_chr) {
        call ib.indexBam as index_chr_bam {
            input:
            bam = bam_chr
        }

        # Mutect
        call m.mutectNoFp as mutect {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            tumor_bam = index_chr_bam.indexed_bam,
            tumor_bam_bai = index_chr_bam.indexed_bam_bai,
            interval_list = target_intervals,
            scatter_count = scatter_count,
            tumor_only = tumor_only
        }
        call bic.bcftoolsIsecComplement as mutect_isec_complement_gnomAD {
            input:
            vcf = mutect.unfiltered_vcf,
            vcf_tbi = mutect.unfiltered_vcf_tbi,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            output_vcf_name = "mutect." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
            output_type = "z"
        }
        call pp.pon2Percent as mutect_pon2 {
            input:
            vcf = mutect_isec_complement_gnomAD.complement_vcf,
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
            interval_bed = target_bed,
            bcbio_filter_string = bcbio_filter_string,
            tumor_sample_name = tumor_sample_name,
            min_var_freq = af_threshold,
            tumor_only = tumor_only
        }
        call bic.bcftoolsIsecComplement as vardict_isec_complement_gnomAD {
            input:
            vcf = vardict.bcbio_filtered_vcf,
            vcf_tbi = vardict.bcbio_filtered_vcf_tbi,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            output_vcf_name = "vardict." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
            output_type = "z"
        }
        call pp.pon2Percent as vardict_pon2 {
            input:
            vcf = vardict_isec_complement_gnomAD.complement_vcf,
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
            interval_bed = target_bed,
            tumor_sample_name = tumor_sample_name,
            tumor_only = tumor_only
        }
        call bic.bcftoolsIsecComplement as lofreq_isec_complement_gnomAD {
            input:
            vcf = lofreq.unfiltered_vcf,
            vcf_tbi = lofreq.unfiltered_vcf_tbi,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            output_vcf_name = "lofreq." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
            output_type = "z"
        }
        call pp.pon2Percent as lofreq_pon2 {
            input:
            vcf = lofreq_isec_complement_gnomAD.complement_vcf,
            vcf2PON = lofreq_pon2_file,
            vcf2PON_tbi = lofreq_pon2_file_tbi,
            caller = "lofreq",
            sample_name = tumor_sample_name
        }

        # 1. Merge Callers Together and create a fake VCF
        # 2. Run fpfilter on fake VCF
        # 3. Run PoN on fake VCF
        # 4. Run VEP on fake VCF
        # 5. Annotate the Original Callers with the annotations
    }

    call mv.mergeVcf as merge_mutect_full {
        input:
            vcfs = mutect.unfiltered_vcf,
            vcf_tbis = mutect.unfiltered_vcf_tbi,
            merged_vcf_basename = "mutect_full." + tumor_sample_name
    }

    call mv.mergeVcf as merge_vardict_full {
        input:
            vcfs = vardict.unfiltered_vcf,
            vcf_tbis = vardict.unfiltered_vcf_tbi,
            merged_vcf_basename = "vardict_full." + tumor_sample_name
    }

    call mv.mergeVcf as merge_lofreq_full {
        input:
            vcfs = lofreq.unfiltered_vcf,
            vcf_tbis = lofreq.unfiltered_vcf_tbi,
            merged_vcf_basename = "lofreq_full." + tumor_sample_name
    }

    output {
        # Alignments
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
        File mutect_full =  merge_mutect_full.merged_vcf                                                          # Raw Mutect Ouput
        # File mutect_pon_annotated_unfiltered_vcf = mutect_gnomad_pon_filters.processed_gnomAD_filtered_vcf      # gnomAD Filtered w/ PoN Annotated
        # File mutect_pon_annotated_filtered_vcf = mutect_pon2.annotated_vcf                                      # gnomAD Filtered + PoN Filtered + PoN2 Filtered w/ VEP Annotation
        # File mutect_pon_total_counts = mutect_gnomad_pon_filters.pon_total_counts                               # PoN Pileup Results

        # Lofreq
        File lofreq_full = merge_lofreq_full.merged_vcf                                                           # Raw Lofreq Ouput
        # File lofreq_pon_annotated_unfiltered_vcf = lofreq_gnomad_pon_filters.processed_gnomAD_filtered_vcf      # gnomAD Filtered w/ PoN Annotated
        # File lofreq_pon_annotated_filtered_vcf = lofreq_pon2.annotated_vcf                                      # gnomAD Filtered + PoN Filtered + PoN2 Filtered w/ VEP Annotation
        # File lofreq_pon_total_counts = lofreq_gnomad_pon_filters.pon_total_counts                               # PoN Pileup Results

        # Vardict
        File vardict_full = merge_vardict_full.merged_vcf                                                            # Raw Vardict Ouput
        # File vardict_pon_annotated_unfiltered_vcf = vardict_gnomad_pon_filters.processed_gnomAD_filtered_vcf      # gnomAD Filtered w/ PoN Annotated
        # File vardict_pon_annotated_filtered_vcf = vardict_pon2.annotated_vcf                                      # gnomAD Filtered + PoN Filtered + PoN2 Filtered w/ VEP Annotation
        # File vardict_pon_total_counts = vardict_gnomad_pon_filters.pon_total_counts                               # PoN Pileup Results

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
