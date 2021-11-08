version 1.0

task lofreqTumorOnly {
    input {
      File reference
      File reference_fai

      File tumor_bam
      File tumor_bam_bai

      File interval_bed
      String? output_name = "lofreq"
      Float ? min_vaf = 0.0001
    }

    Int cores = 2
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))
    runtime {
      docker: "kboltonlab/lofreq:latest"
      memory: "12GB"
      cpu: cores
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        set -o errexit
        set -o nounset

        /opt/lofreq/bin/lofreq indelqual --dindel -f ~{reference} -o output.indel.bam ~{tumor_bam}
        /opt/lofreq/bin/lofreq call --no-default-filter -B -a 1 -b 1 -l ~{interval_bed} -f ~{reference} --call-indels -o lofreq.vcf output.indel.bam --force-overwrite
        /opt/lofreq/bin/lofreq filter -i lofreq.vcf -o ~{output_name}.filtered.vcf -v 5 -a ~{min_vaf} -A 0.9 --sb-incl-indels --print-all

        bgzip ~{output_name}.filtered.vcf && tabix ~{output_name}.filtered.vcf.gz
    >>>

    output {
        File vcf = "~{output_name}.filtered.vcf.gz"
        File vcf_tbi = "~{output_name}.filtered.vcf.gz.tbi"
    }
}

# THIS IS UNTESTED... Someone test please. Thanks
task lofreqNormal {
    input {
      File reference
      File reference_fai

      File tumor_bam
      File tumor_bam_bai

      File? normal_bam
      File? normal_bam_bai

      File interval_bed
      String? output_name = "lofreq"
    }

    Int cpu = 2
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))
    runtime {
      docker: "kboltonlab/lofreq:latest"
      memory: "12GB"
      cores: cpu
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        set -o errexit
        set -o nounset

        /opt/lofreq/bin/lofreq indelqual --dindel -f ~{reference} -o output.indel.bam ~{tumor_bam}
        /opt/lofreq/bin/lofreq somatic --call-indels -n ~{normal_bam} -t output.indel.bam -f ~{reference} -l ~{interval_bed} -o lofreq_ --threads 8
        tabix lofreq_somatic_final.snvs.vcf.gz
        tabix lofreq_somatic_final.indels.vcf.gz
        bcftools concat -a lofreq_somatic_final.snvs.vcf.gz lofreq_somatic_final.indels.vcf.gz > ~{output_name}.vcf
        bgzip ~{output_name}.vcf && tabix ~{output_name}.vcf.gz
    >>>

    output {
        File vcf = "~{output_name}.vcf.gz"
        File vcf_tbi = "~{output_name}.vcf.gz.tbi"
    }
}

workflow lofreqCall {
    input {
        File tumor_bam
        File tumor_bam_bai
        File? normal_bam
        File? normal_bam_bai
        File reference
        File reference_fai
        File interval_bed
        String? output_name = "lofreq"
        Float? min_vaf = 0.0001
        Boolean tumor_only = false
    }

    if (tumor_only) {
        call lofreqTumorOnly {
            input:
                tumor_bam = tumor_bam,
                tumor_bam_bai = tumor_bam_bai,
                reference = reference,
                reference_fai = reference_fai,
                interval_bed = interval_bed,
                output_name = output_name,
                min_vaf = min_vaf
        }
    }

    if (!tumor_only) {
        call lofreqNormal {
            input:
                tumor_bam = tumor_bam,
                tumor_bam_bai = tumor_bam_bai,
                normal_bam = normal_bam,
                normal_bam_bai = normal_bam_bai,
                reference = reference,
                reference_fai = reference_fai,
                interval_bed = interval_bed,
                output_name = output_name
        }
    }

    output {
      File vcf = select_first([lofreqNormal.vcf, lofreqTumorOnly.vcf])
      File vcf_tbi = select_first([lofreqNormal.vcf_tbi, lofreqTumorOnly.vcf_tbi])
    }
}
