version 1.0

import "../types.wdl"

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

workflow manual {
    input {
        File input_bam
        File input_bam_bai

        String? tumor_name = "tumor"
        String tumor_sample_name
        File target_intervals
        File target_bed

        # Reference
        File reference
        File reference_fai
        File reference_dict

        # Variant Calling
        Boolean? arrayMode = true       # Decide if you would rather use the File (--bam_fof) or the Array (--bam) does the same thing, just input type is different
        File pon_normal_bams_file       # on GCP, it's not possible to do File because the file paths are unaccessable for each VM instance, so you have to do ArrayMode
        Array[bam_and_bai] pon_bams     # This is just an array of Files... if you have the bam_fof it's easier just to use above and set this to empty array []
        Boolean tumor_only = true
        Float? af_threshold = 0.0001
        String? pon_pvalue = "2.114164905e-6"

        String bcbio_filter_string = "((FMT/AF * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 3.0) || (FMT/DP < 6500) || (INFO/QUAL < 27)))"

        # PoN2
        File mutect_pon2_file
        File mutect_pon2_file_tbi
        File lofreq_pon2_file
        File lofreq_pon2_file_tbi
        File vardict_pon2_file
        File vardict_pon2_file_tbi

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
            tumor_bam = input_bam,
            tumor_bam_bai = input_bam_bai,
            interval_list = chr_bed,
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
            tumor_bam = input_bam,
            tumor_bam_bai = input_bam_bai,
            interval_bed = chr_bed,
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
            tumor_bam = input_bam,
            tumor_bam_bai = input_bam_bai,
            interval_bed = chr_bed,
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
            bam=input_bam,
            bam_bai=input_bam_bai,
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

    output {
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
    }
}
