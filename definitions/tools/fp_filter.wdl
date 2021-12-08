version 1.0

task fpFilter {
  input {
    File reference
    File reference_fai
    File? reference_dict

    File bam
    File vcf

    String output_vcf_basename = "fpfilter"
    String sample_name = "TUMOR"
    Float? min_var_freq = 0.05
  }

  Int space_needed_gb = 10 + round(size(vcf, "GB")*2 + size([reference, reference_fai, reference_dict, bam], "GB"))
  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    memory: "6GB"
    bootDiskSizeGb: 25
    docker: "kboltonlab/fp_filter-wdl:1.1"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  String output_vcf = output_vcf_basename + ".vcf"
  command <<<
    /usr/bin/perl /usr/bin/fpfilter.pl --bam-readcount /usr/bin/bam-readcount --samtools /usr/bin/samtools --output ~{output_vcf} --reference ~{reference} --bam-file ~{bam} --vcf-file ~{vcf} --sample ~{sample_name} --min-var-freq ~{min_var_freq}

    /usr/bin/bgzip ~{output_vcf} && /usr/bin/tabix ~{output_vcf}.gz
  >>>

  output {
    File filtered_vcf = "~{output_vcf}.gz"
    File filtered_vcf_tbi = "~{output_vcf}.gz.tbi"
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File? reference_dict

    File bam
    File vcf

    String sample_name = "TUMOR"
    Float? min_var_freq
    String output_vcf_basename = "fpfilter"
  }

  call fpFilter {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    vcf=vcf,
    sample_name=sample_name,
    min_var_freq=min_var_freq,
    output_vcf_basename=output_vcf_basename
  }
}
