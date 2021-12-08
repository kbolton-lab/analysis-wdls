version 1.0

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

workflow wf {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    String metric_accumulation_level
  }

  call Metrics {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    metric_accumulation_level=metric_accumulation_level
  }
}
