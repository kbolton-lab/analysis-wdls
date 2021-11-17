version 1.0

task createFakeVcf {
    input {
        File vcf
        String tumor_sample_name
    }

    Int cores = 1
    Int space_needed_gb = 10 + round(2*size(vcf, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
        memory: "6GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -o errexit
        set -o nounset

        echo -e "##fileformat=VCFv4.2" > fake.vcf
        echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> fake.vcf;
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{tumor_sample_name}" >> fake.vcf;
        zcat ~{vcf} | grep -v '#' | awk '{print $1, $2, $3, $4, $5, $6, "PASS\t.\tGT\t0/1"}' OFS='\t' >> fake.vcf;
        bgzip fake.vcf && tabix fake.vcf.gz
    >>>

    output {
        File fake_vcf = "fake.vcf.gz"
        File fake_vcf_tbi = "fake.vcf.gz.tbi"
    }
}


workflow wf {
    input {
        File vcf
        String tumor_sample_name
    }
    call createFakeVcf {
        input:
        vcf=vcf,
        tumor_sample_name=tumor_sample_name
    }
}
