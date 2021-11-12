version 1.0

task bcftoolsFilterBcbio {
    input {
        File vcf
        File vcf_tbi
        String filter_flag = "include"
        String filter_string
        String? output_vcf_prefix = "bcftools_filter"
        String output_type = "z"
    }

    Int space_needed_gb = 5 + 2*round(size(vcf, "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/bst:latest"
      memory: "4GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    String ff = if filter_flag == "include" then "-i" else "-e"
    command <<<
        /usr/local/bin/bcftools filter ~{ff} "~{filter_string}" ~{vcf} --output-type ~{output_type} --output ~{output_vcf_prefix}.vcf.gz -s "BCBIO" -m+

        /usr/local/bin/tabix ~{output_vcf_prefix}.vcf.gz

    >>>

    output {
        File filtered_vcf = "~{output_vcf_prefix}.vcf.gz"
        File filtered_vcf_tbi = "~{output_vcf_prefix}.vcf.gz.tbi"
    }
}
