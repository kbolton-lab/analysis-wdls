version 1.0

task vardictNormal {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        String tumor_sample_name = "TUMOR"
        File? normal_bam
        File? normal_bam_bai
        String? normal_sample_name = "NORMAL"
        File interval_bed
        Float? min_var_freq = 0.005
    }

    Int cores = 32
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))

    runtime {
        docker: "kboltonlab/vardictjava:1.0"
        memory: "128GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
    }

    String output_name = "vardict"
    command <<<
        set -o pipefail
        set -o errexit

        /opt/VarDictJava/build/install/VarDict/bin/VarDict \
            -U -G ~{reference} \
            -f ~{min_var_freq} \
            -N ~{tumor_sample_name} \
            -b "~{tumor_bam}|~{normal_bam}" \
            -c 1 -S 2 -E 3 -g 4 ~{interval_bed} \
            -th 64 | \
        /opt/VarDictJava/build/install/VarDict/bin/testsomatic.R | \
        /opt/VarDictJava/build/install/VarDict/bin/var2vcf_paired.pl \
            -N "~{tumor_sample_name}|~{normal_sample_name}" \
            -f ~{min_var_freq} > ~{output_name}.vcf

        /usr/bin/bgzip ~{output_name}.vcf && /usr/bin/tabix ~{output_name}.vcf.gz
        # extract tumor before bcbio filter or both normal and tumor will be canidates for the filter
        /usr/bin/bcftools view \
            -s ~{tumor_sample_name} \
            --threads 64 ~{output_name}.vcf.gz \
            -Oz -o ~{output_name}.sample.vcf.gz
    >>>

    output {
        File vcf = "~{output_name}.sample.vcf.gz"
        File vcf_tbi = "~{output_name}.sample.vcf.gz.tbi"
    }
}

task vardictTumorOnly {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        String tumor_sample_name = "TUMOR"
        File? normal_bam
        File? normal_bam_bai
        String? normal_sample_name = "NORMAL"
        File interval_bed
        Float? min_var_freq = 0.005
    }

    Int cores = 32
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))

    runtime {
        docker: "kboltonlab/vardictjava:1.0"
        memory: "128GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
    }

    String output_name = "vardict"
    command <<<
        set -o pipefail
        set -o errexit

        /opt/VarDictJava/build/install/VarDict/bin/VarDict \
            -U -G ~{reference} \
            -f ~{min_var_freq} \
            -N ~{tumor_sample_name} \
            -b "~{tumor_bam}" \
            -c 1 -S 2 -E 3 -g 4 ~{interval_bed} \
            -th 64 | \
        /opt/VarDictJava/build/install/VarDict/bin/teststrandbias.R | \
        /opt/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl \
            -N "~{tumor_sample_name}" \
            -E \
            -f ~{min_var_freq} > ~{output_name}.vcf

        /usr/bin/bgzip ~{output_name}.vcf && /usr/bin/tabix ~{output_name}.vcf.gz
    >>>

    output {
        File vcf = "~{output_name}.vcf.gz"
        File vcf_tbi = "~{output_name}.vcf.gz.tbi"
    }
}

workflow vardict {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        String tumor_sample_name = "TUMOR"
        File? normal_bam
        File? normal_bam_bai
        String? normal_sample_name = "NORMAL"
        File interval_bed
        Float? min_var_freq = 0.005
        Boolean tumor_only = false
    }

    if (tumor_only) {
        call vardictTumorOnly {
            input:
                reference = reference,
                reference_fai = reference_fai,
                tumor_bam = tumor_bam,
                tumor_bam_bai = tumor_bam_bai,
                normal_bam = normal_bam,
                normal_bam_bai = normal_bam_bai,
                interval_bed = interval_bed,
                min_var_freq = min_var_freq,
                tumor_sample_name = tumor_sample_name,
                normal_sample_name = normal_sample_name
      }
    }
    if (!tumor_only) {
        call vardictNormal {
            input:
                reference = reference,
                reference_fai = reference_fai,
                tumor_bam = tumor_bam,
                tumor_bam_bai = tumor_bam_bai,
                normal_bam = normal_bam,
                normal_bam_bai = normal_bam_bai,
                interval_bed = interval_bed,
                min_var_freq = min_var_freq,
                tumor_sample_name = tumor_sample_name,
                normal_sample_name = normal_sample_name
      }
    }

    output {
        File vcf = select_first([vardictNormal.vcf, vardictTumorOnly.vcf])
        File vcf_tbi = select_first([vardictNormal.vcf_tbi, vardictTumorOnly.vcf_tbi])
    }
}
