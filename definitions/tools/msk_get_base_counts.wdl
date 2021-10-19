version 1.0

task mskGetBaseCounts {
    input {
        File reference
        File reference_fai
        File reference_dict
        File normal_bams
        String sample_name
        File vcf
        Int mapq
        Int baseq
    }

    Int cores = 4
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float bam_size = size(normal_bams, "GB")
    Float vcf_size = size(vcf, "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + vcf_size)
    runtime {
      docker: "kboltonlab/msk_getbasecounts:1.0"
      cpu: cores
      memory: "64GB"
      disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        set -eou pipefail

        if [[ ~{vcf} == *.vcf.gz ]]; then
            bgzip -d ~{vcf}
            /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam_fof ~{normal_bams} --vcf basename(~{vcf}, ".gz") --output ~{sample_name}.pileup.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
        else
            /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam_fof ~{normal_bams} --vcf ~{vcf} --output ~{sample_name}.pileup.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
        fi
        bgzip ~{sample_name}.pileup.vcf && tabix ~{sample_name}.pileup.vcf.gz
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%RD]\t[%AD]\n' ~{sample_name}.pileup.vcf.gz > ~{sample_name}.pileup.txt
    >>>

    output {
        File pileup = "~{sample_name}.pileup.vcf.gz"
        File pileup_counts = "~{sample_name}.pileup.txt"
    }
}

workflow wf {
    input {
        File reference
        File reference_fai
        File reference_dict
        File normal_bams
        String sample_name
        File vcf
        Int mapq
        Int baseq
    }

    call mskGetBaseCounts {
        input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        normal_bams = normal_bams,
        sample_name = sample_name,
        vcf = vcf,
        mapq = mapq,
        baseq = baseq
    }
}
