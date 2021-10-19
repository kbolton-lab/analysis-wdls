version 1.0

task bcftoolsIsecComplement {
    input {
        File vcf
        File vcf_tbi
        File exclude_vcf
        File exclude_vcf_tbi
        String output_type = "z"
        String? output_vcf_name = "bcftools_isec.vcf.gz"
    }

    Int space_needed_gb = 10 + round(size([vcf, vcf_tbi, exclude_vcf, exclude_vcf_tbi], "GB"))
    runtime {
      docker: "kboltonlab/bst:latest"
      memory: "4GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        /usr/local/bin/bcftools isec -C -w1 ~{vcf} ~{exclude_vcf} --output-type z --output ~{output_vcf_name}
    >>>

    output {
        File complement_vcf = ~{output_vcf_name}
    }
}
