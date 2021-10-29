version 1.0

task varscanNormal {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        File? normal_bam
        File? normal_bam_bai
        File? interval_bed
        Int? strand_filter = 0
        Int? min_coverage = 8
        Float? min_var_freq = 0.005
        Float? p_value = 0.99
        String? tumor_sample_name = "TUMOR" # we don't use it in varscan somatic?
        String tumor_name_replace = "TUMOR"
        String? normal_sample_name = "NORMAL" # we don't use it in varscan somatic?
        String normal_name_replace = "NORMAL"
    }

    Int cores = 2
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))

    runtime {
        docker: "kboltonlab/varscan2:1.0"
        memory: "128GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
    }

    String output_name = "varscan"
    command <<<
        set -o pipefail
        set -o errexit
        export interval_bed="~{interval_bed}"

        if [ -z ${interval_bed+x} ]
        then
            #run without ROI
            java -jar /opt/varscan/VarScan.v2.4.2.jar somatic \
                <(/usr/bin/samtools mpileup --no-baq -f ~{reference} ~{normal_bam} ~{tumor_bam}) \
                "~{output_name}" \
                --strand-filter ~{strand_filter} \
                --min-coverage ~{min_coverage} \
                --min-var-freq ~{min_var_freq} \
                --p-value ~{p_value} \
                --mpileup 1 \
                --output-vcf
        else
            java -jar /opt/varscan/VarScan.v2.4.2.jar somatic \
                <(/usr/bin/samtools mpileup --no-baq -l ~{interval_bed} -f ~{reference} ~{normal_bam} ~{tumor_bam}) \
                "~{output_name}" \
                --strand-filter ~{strand_filter} \
                --min-coverage ~{min_coverage} \
                --min-var-freq ~{min_var_freq} \
                --p-value ~{p_value} \
                --mpileup 1 \
                --output-vcf
        fi
        
        /usr/bin/bgzip ~{output_name}.snp.vcf && /usr/bin/tabix ~{output_name}.snp.vcf.gz
        /usr/bin/bgzip ~{output_name}.indel.vcf && /usr/bin/tabix ~{output_name}.indel.vcf.gz
        /usr/bin/bcftools concat -a -D ~{output_name}.snp.vcf.gz ~{output_name}.indel.vcf.gz -Oz -o ~{output_name}.concat.vcf.gz && /usr/bin/tabix ~{output_name}.concat.vcf.gz

        printf "~{tumor_name_replace} ~{tumor_sample_name}\n~{normal_name_replace} ~{normal_sample_name}\n" > sample_update.txt
        /usr/bin/bcftools reheader ~{output_name}.concat.vcf.gz -s sample_update.txt -o ~{output_name}.vcf.gz && /usr/bin/tabix -p vcf ~{output_name}.vcf.gz


    >>>

    output {
        File vcf = "~{output_name}.vcf.gz"
        File vcf_tbi = "~{output_name}.vcf.gz.tbi"
    }
}


task varscanTumorOnly {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        File? normal_bam
        File? normal_bam_bai
        File? interval_bed
        Int? strand_filter = 0
        Int? min_coverage = 8
        Int? min_reads = 2
        Float? min_var_freq = 0.005
        Float? p_value = 0.99
        String? tumor_sample_name = "TUMOR"
    }

    Int cores = 2
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))

    runtime {
        docker: "kboltonlab/varscan2:1.0"
        memory: "16GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
    }

    String output_name = "varscan"
    command <<<
        set -o pipefail
        set -o errexit

        interval_bed="~{interval_bed}"
        echo ~{tumor_sample_name} > sample.list.txt


        if [ -z ${interval_bed+x} ]
        then
            #run without ROI
            java -jar /opt/varscan/VarScan.v2.4.2.jar mpileup2cns \
                <(/usr/bin/samtools mpileup --no-baq -f ~{reference} ~{tumor_bam}) \
                --strand-filter ~{strand_filter} \
                --min-coverage ~{min_coverage} \
                --min-var-freq ~{min_var_freq} \
                --min-reads2 ~{min_reads} \
                --p-value ~{p_value} \
                --mpileup 1 \
                --output-vcf \
                --variants \
                --vcf-sample-list sample.list.txt \
                > ~{output_name}.vcf
        else
            java -jar /opt/varscan/VarScan.jar mpileup2cns \
                <(/usr/bin/samtools mpileup --no-baq -l ~{interval_bed} -f ~{reference} ~{tumor_bam}) \
                $OUTPUT \
                --strand-filter ~{strand_filter} \
                --min-coverage ~{min_coverage} \
                --min-var-freq ~{min_var_freq} \
                --min-reads2 ~{min_reads} \
                --p-value ~{p_value} \
                --mpileup 1 \
                --variants \
                --output-vcf \
                --vcf-sample-list sample.list.txt \
                > ~{output_name}.vcf
        fi

        /usr/bin/bgzip ~{output_name}.vcf && /usr/bin/tabix ~{output_name}.vcf.gz
    >>>

    output {
        File vcf = "~{output_name}.vcf.gz"
        File vcf_tbi = "~{output_name}.vcf.gz.tbi"
    }
}

workflow varscan {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        String tumor_sample_name = "TUMOR"
        File? normal_bam
        File? normal_bam_bai
        String? normal_sample_name = "NORMAL"
        File? interval_bed
        Float? min_var_freq = 0.005
        Float? p_value = 0.99
        Int? strand_filter = 0
        Int? min_coverage = 8
        Int? min_reads = 2
        Boolean tumor_only = false
    }

    if (tumor_only) {
        call varscanTumorOnly {
            input:
                reference = reference,
                reference_fai = reference_fai,
                tumor_bam = tumor_bam,
                tumor_bam_bai = tumor_bam_bai,
                interval_bed = interval_bed,
                min_var_freq = min_var_freq,
                tumor_sample_name = tumor_sample_name

      }
    }
    if (!tumor_only) {
        call varscanNormal {
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
        File vcf = select_first([varscanNormal.vcf, varscanTumorOnly.vcf])
        File vcf_tbi = select_first([varscanNormal.vcf_tbi, varscanTumorOnly.vcf_tbi])
    }
}
