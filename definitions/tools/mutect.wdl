version 1.0

task mutect {
  input {
    File reference
    File reference_fai
    File reference_dict

    File tumor_bam
    File tumor_bam_bai

    File? normal_bam
    File? normal_bam_bai

    File interval_list
  }

  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Float bam_size = size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB")
  Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_list, "GB"))
  runtime {
    docker: "broadinstitute/gatk:4.2.0.0"
    memory: "32GB"
    bootDiskSizeGb: space_needed_gb
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  String output_vcf = "mutect.filtered.vcf.gz"
  command <<<
    set -o pipefail
    set -o errexit

    if [ -z ~{normal_bam} ]; then
        /gatk/gatk Mutect2 --java-options "-Xmx20g" -O mutect.vcf.gz -R ~{reference} -L ~{interval_list} \
        -I ~{tumor_bam} --read-index ~{tumor_bam_bai} --max-reads-per-alignment-start 0
    else
        NORMAL=`samtools view -H ~{normal_bam} | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1`
        TUMOR=`samtools view -H ~{tumor_bam} | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1`

        /gatk/gatk Mutect2 --java-options "-Xmx20g" -O mutect.vcf.gz -R ~{reference} -L ~{interval_list} \
          -I ~{tumor_bam} --read-index ~{tumor_bam_bai} -tumor "$TUMOR" \
          -I ~{normal_bam} --read-index ~{normal_bam_bai} -normal "$NORMAL" \
          --max-reads-per-alignment-start 0

    /gatk/gatk FilterMutectCalls -R ~{reference} -V mutect.vcf.gz -O ~{output_vcf} #Running FilterMutectCalls on the output vcf.
  >>>

  output {
    File vcf = output_vcf
    File vcf_tbi = output_vcf + ".tbi"
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict
    File tumor_bam
    File tumor_bam_bai
    File? normal_bam
    File? normal_bam_bai
    File interval_list
  }

  call mutect {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    interval_list=interval_list
  }
}
