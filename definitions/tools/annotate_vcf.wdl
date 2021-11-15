version 1.0

task annotateVcf {
    input {
        File vcf
        File vcf_tbi
        File fp_filter
        File fp_filter_tbi
        File pon_filter
        File pon_filter_tbi
        File vep
        String caller_prefix
        String sample_name
    }

    Int space_needed_gb = 10 + 2*round(size([vcf, vcf_tbi, fp_filter, fp_filter_tbi, pon_filter, pon_filter_tbi, vep], "GB"))
    Int cores = 16
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      memory: "96GB"
      docker: "kboltonlab/bst"
      disks: "local-disk ~{space_needed_gb} SSD"
      bootDiskSizeGb: space_needed_gb
      cpu: cores
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        bgzip ~{vep} && tabix ~{vep}.gz
        zcat ~{fp_filter} | tail -n +4  > fp_filter.header;
        zcat ~{pon_filter} | tail -n +3 > pon_filter.header;
        zcat ~{vep} | tail -n +3 > vep.header;

        bcftools annotate --threads 32 -a ~{fp_filter} -h fp_filter.header -c +FILTER ~{vcf} -Oz -o ~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz
        bcftools annotate --threads 32 -a ~{pon_filter} -h pon_filter.header -c PON_RefDepth,PON_AltDepth,PON_FISHER ~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz -Oz -o ~{caller_prefix}.~{sample_name}.fp_filter.pon.annotated.vcf.gz
        bcftools annotate --threads 32 -a ~{vep} -h vep.header -c CSQ ~{caller_prefix}.~{sample_name}.fp_filter.pon.annotated.vcf.gz -Oz -o ~{caller_prefix}.~{sample_name}.final.annotated.vcf.gz -Oz -o ~{caller_prefix}.~{sample_name}.final.annotated.vcf.gz

        tabix ~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz -Oz -o ~{caller_prefix}.~{sample_name}.fp_filter.pon.annotated.vcf.gz
        tabix ~{caller_prefix}.~{sample_name}.final.annotated.vcf.gz
    >>>

    output {
        File final_annotated_vcf = "~{caller_prefix}.~{sample_name}.final.annotated.vcf.gz"
        File final_annotated_vcf_tbi = "~{caller_prefix}.~{sample_name}.final.annotated.vcf.gz.tbi"
        File pon_annotated_vcf = "~{caller_prefix}.~{sample_name}.fp_filter.pon.annotated.vcf.gz"
        File pon_annotated_vcf_tbi = "~{caller_prefix}.~{sample_name}.fp_filter.pon.annotated.vcf.gz.tbi"
    }
}
