version 1.0

task mskGetBaseCountsWithFile {
    input {
        File reference
        File reference_fai
        File reference_dict
        File normal_bams
        String? pon_final_name = "pon.pileup"
        File vcf
        Int? mapq = 5
        Int? baseq = 5
    }

    Int cores = 4
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float bam_size = size(normal_bams, "GB")
    Float vcf_size = size(vcf, "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + vcf_size)
    runtime {
      docker: "kboltonlab/msk_getbasecounts:3.0"
      cpu: cores
      memory: "64GB"
      disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        set -eou pipefail

        if [[ ~{vcf} == *.vcf.gz ]]; then
            bgzip -d ~{vcf}
            vcf_file=~{vcf}
            /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam_fof ~{normal_bams} --vcf "${vcf_file%.*}" --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
        else
            /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam_fof ~{normal_bams} --vcf ~{vcf} --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
        fi
        bgzip ~{pon_final_name}.vcf && tabix ~{pon_final_name}.vcf.gz
    >>>

    output {
        File pileup = "~{pon_final_name}.vcf.gz"
        File pileup_tbi = "~{pon_final_name}.vcf.gz.tbi"
    }
}

task mskGetBaseCountsWithArray {
    input {
        File reference
        File reference_fai
        File reference_dict
        Array[File] normal_bams
        Array[File] normal_bams_bai
        String? pon_final_name = "pon.pileup"
        File vcf
        Int? mapq = 5
        Int? baseq = 5
    }

    Int cores = 4
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float bam_size = size(normal_bams, "GB")
    Float vcf_size = size(vcf, "GB")
    Int space_needed_gb = 10 + round(reference_size + 5*bam_size + vcf_size)
    runtime {
      docker: "kboltonlab/msk_getbasecounts:3.0"
      cpu: cores
      memory: "128GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      bootDiskSizeGb: space_needed_gb
    }

    command <<<
        set -eou pipefail

        echo "REFERENCE: ~{reference_size}"
        echo "BAM: ~{bam_size}"
        echo "VCF: ~{vcf_size}"
        echo "SPACE_NEEDED: ~{space_needed_gb}"

        bam_string=""
        for bam in ~{sep=" " normal_bams}; do
            sample_name=$(samtools view -H $bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
            bam_string="$bam_string --bam ${sample_name}:$bam"
        done

        if [[ ~{vcf} == *.vcf.gz ]]; then
            bgzip -d ~{vcf}
            vcf_file=~{vcf}
            /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} $bam_string --vcf "${vcf_file%.*}" --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
        else
            /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} $bam_string --vcf ~{vcf} --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
        fi
        bgzip ~{pon_final_name}.vcf && tabix ~{pon_final_name}.vcf.gz
    >>>

    output {
        File pileup = "~{pon_final_name}.vcf.gz"
        File pileup_tbi = "~{pon_final_name}.vcf.gz.tbi"
    }
}

workflow mskGetBaseCounts {
    input {
        File reference
        File reference_fai
        File reference_dict
        Boolean arrayMode = false
        File normal_bams
        Array[File] bams
        Array[File] bams_bai
        String? pon_final_name = "pon.pileup"
        File vcf
        Int? mapq = 5
        Int? baseq = 5
    }
    if (arrayMode) {
        call mskGetBaseCountsWithArray {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            normal_bams = bams,
            normal_bams_bai = bams_bai,
            pon_final_name = pon_final_name,
            vcf = vcf,
            mapq = mapq,
            baseq = baseq
        }
    }
    if (!arrayMode) {
        call mskGetBaseCountsWithFile {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            normal_bams = normal_bams,
            pon_final_name = pon_final_name,
            vcf = vcf,
            mapq = mapq,
            baseq = baseq
        }
    }

    output {
        File pileup = select_first([mskGetBaseCountsWithArray.pileup, mskGetBaseCountsWithFile.pileup])
        File pileup_tbi = select_first([mskGetBaseCountsWithArray.pileup_tbi, mskGetBaseCountsWithFile.pileup_tbi])
    }
}
