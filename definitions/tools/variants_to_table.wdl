version 1.0

task variantsToTable {
  input {
    File reference
    File reference_fai
    File reference_dict
    File vcf
    File vcf_tbi
    Array[String] fields = ["CHROM", "POS", "ID", "REF", "ALT", "set"]
    Array[String] genotype_fields = ["GT", "AD", "DP", "AF"]
  }

  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Float vcf_size = size([vcf, vcf_tbi], "GB")
  Int space_needed_gb = 10 + round(vcf_size*2 + reference_size)
  runtime {
    memory: "6GB"
    bootDiskSizeGb: 25
    docker: "broadinstitute/gatk:4.1.8.1"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  String outfile = "variants.tsv"
  command <<<
    /gatk/gatk --java-options -Xmx4g VariantsToTable -O ~{outfile} \
    -R ~{reference} --variant ~{vcf} \
    ~{sep=" " prefix("-F ", fields)}
  >>>

  output {
    File variants_tsv = outfile
  }
}


workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict
    File vcf
    File vcf_tbi
    Array[String]? fields
    Array[String]? genotype_fields
  }

  call variantsToTable {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=vcf,
    vcf_tbi=vcf_tbi,
    fields=fields,
    genotype_fields=genotype_fields
  }
}
