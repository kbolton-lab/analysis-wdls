version 1.0

task bcftoolsFilterBcbio {
    input {
        File vcf
        File vcf_tbi
        String filter_flag = "include"
        String filter_string
        String? output_vcf_name = "bcftools_filter"
        String output_type = "z"
    }

    Int space_needed_gb = 5 + round(size(vcf, "GB"))
    runtime {
      docker: "kboltonlab/bst:latest"
      memory: "4GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
    }

    String ff = if filter_flag == "include" then "-i" else "-e"
    command <<<
        /usr/local/bin/bcftools filter ~{ff} "~{filter_string}" ~{vcf} --output-type ~{output_type} --output ~{output_vcf_name}.vcf.gz -s "BCBIO" -m+
        
        /usr/local/bin/tabix ~{output_vcf_name}.vcf.gz

    >>>

    output {
        File filtered_vcf = "~{output_vcf_name}.vcf.gz"
        File filtered_vcf_tbi = "~{output_vcf_name}.vcf.gz.tbi"
    }
}
