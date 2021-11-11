version 1.0

import "../types.wdl"
import "../tools/merge_vcf.wdl" as mv

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
    Int space_needed_gb = 100 + round(reference_size + 2*bam_size + vcf_size)
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

task mskGetBaseCounts {
    input {
        File reference
        File reference_fai
        File reference_dict
        bam_and_bai normal_bam
        String? pon_final_name = "pon.pileup"
        File vcf
        Int? mapq = 5
        Int? baseq = 5
    }

    Int cores = 4
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float bam_size = size([normal_bam.bam, normal_bam.bai], "GB")
    Float vcf_size = size(vcf, "GB")
    Int space_needed_gb = 5 + round(reference_size + 2*bam_size + vcf_size)
    runtime {
      docker: "kboltonlab/msk_getbasecounts:3.0"
      cpu: cores
      memory: "128GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      bootDiskSizeGb: space_needed_gb
    }

    File bam = normal_bam.bam

    command <<<
        set -eou pipefail

        echo "REFERENCE: ~{reference_size}"
        echo "BAM: ~{bam_size}"
        echo "VCF: ~{vcf_size}"
        echo "SPACE_NEEDED: ~{space_needed_gb}"

        sample_name=$(samtools view -H ~{bam} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)

        if [[ ~{vcf} == *.vcf.gz ]]; then
            bgzip -d ~{vcf}
            vcf_file=~{vcf}
            /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam ${sample_name}:~{bam} --vcf "${vcf_file%.*}" --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
        else
            /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam ${sample_name}:~{bam} --vcf ~{vcf} --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
        fi
        bgzip ~{pon_final_name}.vcf && tabix ~{pon_final_name}.vcf.gz
    >>>

    output {
        File pileup = "~{pon_final_name}.vcf.gz"
        File pileup_tbi = "~{pon_final_name}.vcf.gz.tbi"
    }
}

workflow wf {
    input {
        File reference
        File reference_fai
        File reference_dict
        Boolean arrayMode = false
        File normal_bams_file
        Array[bam_and_bai] pon_bams
        String? pon_final_name = "pon.pileup"
        File vcf
        Int? mapq = 5
        Int? baseq = 5
    }
    if (arrayMode) {
        scatter (pon_bam in pon_bams){
            call mskGetBaseCounts {
                input:
                reference = reference,
                reference_fai = reference_fai,
                reference_dict = reference_dict,
                normal_bam = pon_bam,
                pon_final_name = pon_final_name,
                vcf = vcf,
                mapq = mapq,
                baseq = baseq
            }
        }

        call mv.mergeVcf as merge {
            input:
                vcfs = mskGetBaseCounts.pileup,
                vcf_tbis = mskGetBaseCounts.pileup_tbi
        }
    }
    if (!arrayMode) {
        call mskGetBaseCountsWithFile {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            normal_bams = normal_bams_file,
            pon_final_name = pon_final_name,
            vcf = vcf,
            mapq = mapq,
            baseq = baseq
        }
    }

    output {
        File pileup = select_first([merge.merged_vcf, mskGetBaseCountsWithFile.pileup])
        File pileup_tbi = select_first([merge.merged_vcf_tbi, mskGetBaseCountsWithFile.pileup_tbi])
    }
}
