version 1.0

task bcftoolsNorm {
    input {
        File reference
        File reference_fai
        File reference_dict

        File vcf
        File vcf_tbi
    }

    Int space_needed_gb = 10 + round(size([vcf, vcf_tbi], "GB") + size([reference, reference_fai, reference_dict], "GB"))
    runtime {
      memory: "9GB"
      docker: "mgibio/bcftools-cwl:1.9"
      disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        /opt/bcftools/bin/bcftools norm -any --multiallelics --output-type z --output bcftools_norm.vcf.gz ~{vcf} -f ~{reference}

    >>>

    output {
      File normalized_vcf = "bcftools_norm.vcf.gz"
    }
}
