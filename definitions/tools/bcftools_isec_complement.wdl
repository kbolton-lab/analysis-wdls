version 1.0

task bcftoolsIsecComplement {
    input {
        File vcf
        File vcf_tbi
        File exclude_vcf
        File exclude_vcf_tbi
        String output_type = "z"
        String? output_vcf_name = "bcftools_isec.vcf"
    }

    Int space_needed_gb = 10 + round(size([vcf, vcf_tbi, exclude_vcf, exclude_vcf_tbi], "GB"))
    runtime {
      docker: "kboltonlab/bst:latest"
      memory: "4GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        /usr/local/bin/bcftools isec -C -w1 ~{vcf} ~{exclude_vcf} --output-type ~{output_type} --output ~{output_vcf_name}.gz && /usr/local/bin/tabix ~{output_vcf_name}.gz

    >>>

    output {
        File complement_vcf = "~{output_vcf_name}.gz"
        File complement_vcf_tbi = "~{output_vcf_name}.gz.tbi"
    }
}

workflow wf {
    input {
        File vcf
        File vcf_tbi
        File exclude_vcf
        File exclude_vcf_tbi
        String output_type = "z"
        String? output_vcf_name = "bcftools_isec.vcf.gz"
    }

    call bcftoolsIsecComplement {
        input:
        vcf = vcf,
        vcf_tbi = vcf_tbi,
        exclude_vcf = exclude_vcf,
        exclude_vcf_tbi = exclude_vcf_tbi,
        output_type = output_type,
        output_vcf_name = output_vcf_name
    }
}
