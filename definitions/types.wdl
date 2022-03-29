version 1.0

struct Sequence {
  File? bam
  File? fastq1
  File? fastq2
}
# assume either bam or fastqs defined
struct SequenceData {
  Sequence sequence
  String readgroup
}

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

# VEP can utilize other files for custom annotations outside of the normal available plugins
# This structure format is ported over from MGI's CWL Pipelines so it could potentially be better integrated
struct VepCustomAnnotation {
    Boolean check_existing
    File custom_file
    String name
    String data_format  # enum, ['bed', 'gff', 'gtf', 'vcf', 'bigwig']
    String method  # enum, ['exact', 'overlap']
    Boolean force_report_coordinates
    Array[String]? vcf_fields
    Array[File]? secondary_files
}

# The SpliceAI Plugin requires two files be provided, so rather than having 4 files passed individually
# it makes things cleaner to have a single structure hold all four.
struct VepSpliceAIPlugin {
    File? spliceAI_snv
    File? spliceAI_snv_tbi
    File? spliceAI_indel
    File? spliceAI_indel_tbi
}
