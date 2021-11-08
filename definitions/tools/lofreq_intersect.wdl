version 1.0

task lofreqIntersect {
    input {
        File call_vcf
        File pass_vcf
    }

    Int space_needed_gb = 10 + round(size([call_vcf, pass_vcf], "GB"))
    runtime {
      memory: "9GB"
      docker: "kboltonlab/bst"
      disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        printf "##FILTER=<ID=CALL,Description=\"A variant that was called by Lofreq's Caller without any filters\">" > lofreq.header;
        zcat ~{call_vcf} | sed 's/PASS/CALL/g' > call_to_pass.vcf
        bgzip call_to_pass.vcf && tabix call_to_pass.vcf.gz
        bcftools annotate --threads 32 -a ~{pass_vcf} -h lofreq.header -c FILTER call_to_pass.vcf.gz -Oz -o lofreq_intersect.vcf.gz
    >>>

    output {
        File intersect_vcf = "lofreq_intersect.vcf.gz"
    }
}

workflow wf {
    input {
        File call_vcf
        File pass_vcf
    }

    call lofreqIntersect {
        input:
        call_vcf = call_vcf,
        pass_vcf = pass_vcf
    }

    output {
        File intersect_vcf = lofreqIntersect.intersect_vcf
    }
}
