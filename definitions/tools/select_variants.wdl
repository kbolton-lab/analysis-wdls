version 1.0

task selectVariants {
  input {
    File reference
    File reference_fai
    File reference_dict

    File vcf
    File vcf_tbi

    File? interval_list
    Boolean exclude_filtered = false
    String output_vcf_basename = "select_variants"
    Array[String]? samples_to_include  # include genotypes from this sample

    # ENUM: one of ["INDEL", "SNP", "MIXED", "MNP", "SYMBOLIC", "NO_VARIATION"]
    String? select_type
  }

  Int space_needed_gb = 5 + round(size([vcf, vcf_tbi], "GB")*3 + size([reference, reference_fai, reference_dict, interval_list], "GB"))
  Int preemptible = 1
  Int maxRetries = 0
  Int cores = 1
  runtime {
    docker: "broadinstitute/gatk:4.2.0.0"
    memory: "6GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  String outfile = "~{output_vcf_basename}.vcf.gz"
  Array[String] samples = if defined(samples_to_include) then prefix("--sample-name ", select_first([samples_to_include])) else []
  command <<<
    /gatk/gatk --java-options -Xmx4g SelectVariants -O ~{outfile} \
    -R ~{reference} \
    --variant ~{vcf} \
    ~{if defined(interval_list) then "-L ~{interval_list}" else ""} \
    ~{if exclude_filtered then "--exclude-filtered ~{select_first([exclude_filtered])}" else ""} \
    ~{sep=" " samples} \
    ~{if defined(select_type) then "-select-type ~{select_type}" else ""}
  >>>

  output {
    File filtered_vcf = outfile
    File filtered_vcf_tbi = "~{outfile}.tbi"
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict
    File vcf
    File vcf_tbi
    File? interval_list
    Boolean? exclude_filtered
    String output_vcf_basename = "select_variants"
    Array[String]? samples_to_include  # include genotypes from this sample
    # ENUM: one of ["INDEL", "SNP", "MIXED", "MNP", "SYMBOLIC", "NO_VARIATION"]
    String? select_type
  }

  call selectVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=vcf,
    vcf_tbi=vcf_tbi,
    interval_list=interval_list,
    exclude_filtered=exclude_filtered,
    output_vcf_basename=output_vcf_basename,
    samples_to_include=samples_to_include,
    select_type=select_type
  }
}
