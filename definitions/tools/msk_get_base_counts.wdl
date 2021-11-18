version 1.0

import "../types.wdl"
import "../tools/bcftools_merge.wdl" as bm

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
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      docker: "kboltonlab/msk_getbasecounts:3.0"
      cpu: cores
      memory: "24GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
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
        bam_and_bai_array normal_bams
        String? pon_final_name = "pon.pileup"
        File vcf
        Int? mapq = 5
        Int? baseq = 5
    }

    Int cores = 4
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float bam_size = size(normal_bams.bams, "GB")
    Float vcf_size = size(vcf, "GB")
    Int space_needed_gb = 10 + round(reference_size + 6*bam_size*vcf_size)
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      docker: "kboltonlab/msk_getbasecounts:3.0"
      cpu: cores
      memory: "24GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        set -eou pipefail

        echo "REFERENCE: ~{reference_size}"
        echo "BAM: ~{bam_size}"
        echo "VCF: ~{vcf_size}"
        echo "SPACE_NEEDED: ~{space_needed_gb}"

        bam_string=""
        for bam in ~{sep=" " normal_bams.bams}; do
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
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + vcf_size)
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      docker: "kboltonlab/msk_getbasecounts:3.0"
      cpu: cores
      memory: "24GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        set -eou pipefail

        echo "REFERENCE: ~{reference_size}"
        echo "BAM: ~{bam_size}"
        echo "VCF: ~{vcf_size}"
        echo "SPACE_NEEDED: ~{space_needed_gb}"

        sample_name=$(samtools view -H ~{normal_bam.bam} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
        if [[ $(zgrep -v '#' ~{vcf} | wc -l) -lt 1 ]]; then
            echo "PRINT: ~{vcf}"
            printf "##fileformat=VCFv4.2\n" > ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth matching reference (REF) allele\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth matching alternate (ALT) allele\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=VF,Number=1,Type=Float,Description=\"Variant frequence (AD/DP)\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=DPP,Number=1,Type=Integer,Description=\"Depth on postitive strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=DPN,Number=1,Type=Integer,Description=\"Depth on negative strand\"\n>" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=RDP,Number=1,Type=Integer,Description=\"Reference depth on postitive strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=RDN,Number=1,Type=Integer,Description=\"Reference depth on negative strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=ADP,Number=1,Type=Integer,Description=\"Alternate depth on postitive strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=ADN,Number=1,Type=Integer,Description=\"Alternate depth on negative strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=DPF,Number=1,Type=Integer,Description=\"Total fragment depth\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=RDF,Number=1,Type=Float,Description=\"Fragment depth matching reference (REF) allele\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=ADF,Number=1,Type=Float,Description=\"Fragment depth matching alternate (ALT) allele\">\n" >> ~{pon_final_name}.vcf
            printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample_name}\n" >> ~{pon_final_name}.vcf
        else
            echo "SCRIPT: ~{normal_bam.bam}"
            if [[ ~{vcf} == *.vcf.gz ]]; then
                bgzip -d ~{vcf}
                vcf_file=~{vcf}
                /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam ${sample_name}:~{normal_bam.bam} --vcf "${vcf_file%.*}" --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
            else
                /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam ${sample_name}:~{normal_bam.bam} --vcf ~{vcf} --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
            fi
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
        Boolean multi = false
        File normal_bams_file
        Array[bam_and_bai] pon_bams
        Array[bam_and_bai_array] pon_bams_array
        String? pon_final_name = "pon.pileup"
        File vcf
        Int? mapq = 5
        Int? baseq = 5
    }
    if (arrayMode) {
        if (multi) {
            scatter (bam_set in pon_bams_array){
                call mskGetBaseCountsWithArray {
                    input:
                    reference = reference,
                    reference_fai = reference_fai,
                    reference_dict = reference_dict,
                    normal_bams = bam_set,
                    pon_final_name = pon_final_name,
                    vcf = vcf,
                    mapq = mapq,
                    baseq = baseq
                }
            }

            call bm.bcftoolsMerge as mergeMulti {
                input:
                    vcfs = mskGetBaseCountsWithArray.pileup,
                    vcf_tbis = mskGetBaseCountsWithArray.pileup_tbi
            }
        }

        if (!multi) {
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

            call bm.bcftoolsMerge as merge {
                input:
                    vcfs = mskGetBaseCounts.pileup,
                    vcf_tbis = mskGetBaseCounts.pileup_tbi
            }
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
        File pileup = select_first([merge.merged_vcf, mergeMulti.merged_vcf, mskGetBaseCountsWithFile.pileup])
        File pileup_tbi = select_first([merge.merged_vcf_tbi, mergeMulti.merged_vcf_tbi, mskGetBaseCountsWithFile.pileup_tbi])
    }
}
