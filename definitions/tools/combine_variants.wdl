version 1.0

task combineVariants {
  input {
    File reference
    File reference_fai
    File reference_dict

    File mutect_vcf
    File mutect_vcf_tbi

    File varscan_vcf
    File varscan_vcf_tbi

    File strelka_vcf
    File strelka_vcf_tbi

    File pindel_vcf
    File pindel_vcf_tbi
  }

  Float ref_size = size([reference, reference_fai, reference_dict], "GB")
  Float mutect_size = size([mutect_vcf, mutect_vcf_tbi], "GB")
  Float varscan_size = size([varscan_vcf, varscan_vcf_tbi], "GB")
  Float strelka_size = size([strelka_vcf, strelka_vcf_tbi], "GB")
  Float pindel_size = size([pindel_vcf, pindel_vcf_tbi], "GB")
  Int space_needed_gb = 10 + round(ref_size + mutect_size + varscan_size + strelka_size + pindel_size)*2
  runtime {
    memory: "9GB"
    bootDiskSizeGb: 25
    docker: "mgibio/gatk-cwl:3.6.0"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  String outfile = "combined.vcf.gz"
  command <<<
    /usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T CombineVariants \
    -genotypeMergeOptions PRIORITIZE --rod_priority_list mutect,varscan,strelka,pindel -o ~{outfile} \
    -R ~{reference} \
    --variant:mutect ~{mutect_vcf} \
    --variant:varscan ~{varscan_vcf} \
    --variant:strelka ~{strelka_vcf} \
    --variant:pindel ~{pindel_vcf}
  >>>

  output {
    File combined_vcf = outfile
    File combined_vcf_tbi = outfile + ".tbi"
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict

    File mutect_vcf
    File mutect_vcf_tbi

    File varscan_vcf
    File varscan_vcf_tbi

    File strelka_vcf
    File strelka_vcf_tbi

    File pindel_vcf
    File pindel_vcf_tbi
  }

  call combineVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    mutect_vcf=mutect_vcf,
    mutect_vcf_tbi=mutect_vcf_tbi,
    varscan_vcf=varscan_vcf,
    varscan_vcf_tbi=varscan_vcf_tbi,
    strelka_vcf=strelka_vcf,
    strelka_vcf_tbi=strelka_vcf_tbi,
    pindel_vcf=pindel_vcf,
    pindel_vcf_tbi=pindel_vcf_tbi
  }
}
