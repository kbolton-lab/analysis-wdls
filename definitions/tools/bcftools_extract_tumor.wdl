version 1.0

task extractTumor {
    input {
        File vcf
        String tumor_sample_name
        String output_vcf_basename
        String output_type = "z"
    }

    Int space_needed_gb = 10 + round(size(vcf, "GB"))
    runtime {
      docker: "kboltonlab/bst:latest"
      memory: "4GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        /usr/local/bin/bcftools view -s ~{tumor_sample_name} ~{vcf} --output-type ~{output_type} --output ~{output_vcf_basename}.vcf.gz

        /usr/local/bin/tabix ~{output_vcf_basename}.vcf.gz

    >>>

    output {
        File filtered_vcf = "~{output_vcf_basename}.vcf.gz"
        File filtered_vcf_tbi = "~{output_vcf_basename}.vcf.gz.tbi"
    }
}
