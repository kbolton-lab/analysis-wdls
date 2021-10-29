version 1.0

task bcftoolsNorm {
    input {
        File reference
        File reference_fai
        # File reference_dict

        File vcf
        File vcf_tbi
    }

    Int space_needed_gb = 10 + round(size([vcf, vcf_tbi], "GB") + size([reference, reference_fai], "GB"))
    runtime {
      memory: "9GB"
      docker: "kboltonlab/bst"
      disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        /usr/local/bin/bcftools norm --multiallelics -any --output-type z --output bcftools_norm.vcf.gz ~{vcf} -f ~{reference}

        /usr/local/bin/tabix bcftools_norm.vcf.gz
    >>>

    output {
      File normalized_vcf = "bcftools_norm.vcf.gz"
      File normalized_vcf_tbi = "bcftools_norm.vcf.gz.tbi"
    }
}
