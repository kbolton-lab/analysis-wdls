version 1.0

task bcftoolsNorm {
    input {
        File reference
        File reference_fai
        # File reference_dict

        File vcf
        File vcf_tbi
    }

    Int space_needed_gb = 5 + round(size([vcf, vcf_tbi], "GB") * 2 + size([reference, reference_fai], "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      memory: "6GB"
      docker: "kboltonlab/bst"
      disks: "local-disk ~{space_needed_gb} SSD"
      cpu: cores
      preemptible: preemptible
      maxRetries: maxRetries
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
