version 1.0

import "../tools/kallisto.wdl" as k
import "../tools/transcript_to_gene.wdl" as ttg
import "../tools/samtools_sort.wdl" as ss
import "../tools/merge_bams.wdl" as mb
import "../tools/index_bam.wdl" as ib
import "../tools/mark_duplicates_and_sort.wdl" as mdas
import "../tools/stringtie.wdl" as st
import "../tools/generate_qc_metrics.wdl" as gqm
import "../tools/bam_to_bigwig.wdl" as btb
import "../subworkflows/bam_to_trimmed_fastq_and_hisat_alignments.wdl" as bttfaha

workflow rnaseq {
  input {
    File reference
    File reference_fai
    File reference_dict

    File reference_index
    File reference_index_1ht2
    File reference_index_2ht2
    File reference_index_3ht2
    File reference_index_4ht2
    File reference_index_5ht2
    File reference_index_6ht2
    File reference_index_7ht2
    File reference_index_8ht2

    File reference_annotation
    Array[String] read_group_id
    Array[File] instrument_data_bams
    Array[Array[String]] read_group_fields
    String? strand  # [first, second, unstranded]
    String sample_name

    File trimming_adapters
    String trimming_adapter_trim_end
    Int trimming_adapter_min_overlap
    Int trimming_max_uncalled
    Int trimming_min_readlength

    File kallisto_index
    File gene_transcript_lookup_table
    File refFlat
    File? ribosomal_intervals
  }

  scatter(bam in instrument_data_bams) {
    scatter(id in read_group_id) {
      scatter(fields in read_group_fields) {
        call bttfaha.bamToTrimmedFastqAndHisatAlignments {
          input:
          bam=bam,
          read_group_id=id,
          read_group_fields=fields,
          adapters=trimming_adapters,
          adapter_trim_end=trimming_adapter_trim_end,
          adapter_min_overlap=trimming_adapter_min_overlap,
          max_uncalled=trimming_max_uncalled,
          min_readlength=trimming_min_readlength,
          reference_index=reference_index,
          reference_index_1ht2=reference_index_1ht2,
          reference_index_2ht2=reference_index_2ht2,
          reference_index_3ht2=reference_index_3ht2,
          reference_index_4ht2=reference_index_4ht2,
          reference_index_5ht2=reference_index_5ht2,
          reference_index_6ht2=reference_index_6ht2,
          reference_index_7ht2=reference_index_7ht2,
          reference_index_8ht2=reference_index_8ht2,
          strand=strand
        }
      }
    }
  }

  call k.kallisto {
    input:
    kallisto_index=kallisto_index,
    strand=strand,
    fastqs=flatten(flatten(flatten(bamToTrimmedFastqAndHisatAlignments.fastqs)))
  }

  call ttg.transcriptToGene {
    input:
    transcript_table_h5=kallisto.expression_transcript_h5,
    gene_transcript_lookup_table=gene_transcript_lookup_table
  }

  # TODO: remove extra sort
  call mb.mergeBams as merge {
    input: bams=flatten(flatten(bamToTrimmedFastqAndHisatAlignments.aligned_bam))
  }

  call ss.samtoolsSort as positionSort {
    input: input_bam=merge.merged_bam
  }

  call ib.indexBam {
    input: bam=positionSort.sorted_bam
  }

  call mdas.markDuplicatesAndSort as markDup {
    input:
    bam=indexBam.indexed_bam,
    input_sort_order="coordinate"
  }

  call st.stringtie {
    input:
    bam=indexBam.indexed_bam,
    reference_annotation=reference_annotation,
    sample_name=sample_name,
    strand=strand
  }

  call gqm.generateQcMetrics {
    input:
    refFlat=refFlat,
    ribosomal_intervals=ribosomal_intervals,
    strand=strand,
    bam=indexBam.indexed_bam
  }

  call btb.bamToBigwig as cgpbigwigBamcoverage {
    input:
    bam=markDup.sorted_bam,
    bam_bai=markDup.sorted_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict
  }

  output {
    File final_bam = markDup.sorted_bam
    File final_bam_bai = markDup.sorted_bam_bai
    File stringtie_transcript_gtf = stringtie.transcript_gtf
    File stringtie_gene_expression_tsv = stringtie.gene_expression_tsv
    File transcript_abundance_tsv = kallisto.expression_transcript_table
    File transcript_abundance_h5 = kallisto.expression_transcript_h5
    File gene_abundance = transcriptToGene.gene_abundance
    File metrics = generateQcMetrics.metrics
    File? chart = generateQcMetrics.chart
    File fusion_evidence = kallisto.fusion_evidence
    File bamcoverage_bigwig = cgpbigwigBamcoverage.outfile
  }
}
