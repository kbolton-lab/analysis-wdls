version 1.0

task bcftoolsFilterBcbio {
    input {
        File vcf
        File vcf_tbi
        String filter_flag = "include"
        String filter_string
        String? output_vcf_name = "bcftools_filter.vcf.gz"
        String output_type = "z"
    }

    Int space_needed_gb = 5 + round(size(vcf, "GB"))
    runtime {
      docker: "mgibio/bcftools-cwl:1.9"
      memory: "4GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
    }

    String ff = if filter_flag == "include" then "-i" else "-e"
    command <<<
        /opt/bcftools/bin/bcftools filter ~{ff} ~{filter_string} ~{vcf} --output-type ~{output_type} --output ~{output_vcf_name} -s "BCBIO" -m+
    >>>

    output {
        File filtered_vcf = output_vcf_name
    }
}
