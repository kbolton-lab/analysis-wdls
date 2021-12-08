version 1.0

task pon2Percent {
    input {
        File vcf
        File vcf2PON
        File vcf2PON_tbi
        String caller = "caller"
        String sample_name = "tumor"
    }

    Int space_needed_gb = 10 + round(size([vcf, vcf2PON, vcf2PON_tbi],"GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/bst:latest"
      memory: "6GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        set -eou pipefail

        export name=~{caller}.~{sample_name}.pon2.annotated.vcf.gz

        printf "##INFO=<ID=PON_2AT2_percent,Number=1,Type=Integer,Description=\"If 2 PoN samples have variant at >=2 percent\">\n" > pon2.header;
        printf "##INFO=<ID=PON_NAT2_percent,Number=1,Type=Integer,Description=\"Number of samples with variant at >=2 percent\">\n" >> pon2.header;
        printf "##INFO=<ID=PON_MAX_VAF,Number=1,Type=Float,Description=\"The maximum VAF found in the PoN Samples\">\n" >> pon2.header;
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t1\t%INFO/NS\t%INFO/max_VAF\n" ~{vcf2PON} > normal2.txt
        bgzip -f normal2.txt
        tabix -f -s1 -b2 -e2 normal2.txt.gz
        bcftools annotate -a normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent,PON_NAT2_percent,PON_MAX_VAF ~{vcf} -Oz -o $name
        tabix $name
    >>>

    output {
        File annotated_vcf = "~{caller}.~{sample_name}.pon2.annotated.vcf.gz"
        File annotated_vcf_tbi = "~{caller}.~{sample_name}.pon2.annotated.vcf.gz.tbi"
    }
}

workflow wf {
    input {
        File vcf
        File vcf2PON
        File vcf2PON_tbi
        String caller = "caller"
        String sample_name = "tumor"
    }

    call pon2Percent {
        input:
            vcf = vcf,
            vcf2PON = vcf2PON,
            vcf2PON_tbi = vcf2PON_tbi,
            caller = caller,
            sample_name = sample_name
    }
}
