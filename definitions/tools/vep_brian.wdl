version 1.0

import "../types.wdl"

task vep {
  input {
    File vcf
    File cache_dir_zip
    File reference
    File reference_fai
    String ensembl_assembly
    String ensembl_version
    String ensembl_species
    Array[String] plugins
    Boolean coding_only = false
    Array[VepCustomAnnotation] custom_annotations = []
    Boolean everything = true
    File? synonyms_file
    String? additional_args = "--sift p --polyphen p --pick --pick_order canonical,rank,mane,ccds,appris,tsl,biotype,length --everything 1 --merged --check_existing --buffer_size 1000 --af_gnomad"
  }

  Float cache_size = 2*size(cache_dir_zip, "GB")  # doubled to unzip
  Float vcf_size = 2*size(vcf, "GB")  # doubled for output vcf
  Float reference_size = size([reference, reference_fai], "GB")
  Int space_needed_gb = 10 + round(reference_size + vcf_size + cache_size + size(synonyms_file, "GB"))
  runtime {
    memory: "64GB"
    bootDiskSizeGb: 30
    cpu: 4
    docker: "kboltonlab/ic_vep"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  String annotated_path = basename(basename(vcf, ".gz"), ".vcf") + "_annotated.vcf"
  String cache_dir = basename(cache_dir_zip, ".zip")

  Int annotation_len = length(custom_annotations)  

  command <<<
     
    if [[ ~{annotation_len} -ge 1 ]]; then
      cus_annot=$(/usr/bin/python3 /opt/bin/theToast.py ~{write_json(custom_annotations)})
    else
      cus_annot=""
    fi
    echo $cus_annot
  
    unzip -qq ~{cache_dir_zip}

    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
    --format vcf \
    --vcf \
    --fork 4 \
    --terms SO \
    --transcript_version \
    --offline \
    --cache \
    --symbol \
    -o ~{annotated_path} \
    -i ~{vcf} \
    ~{if defined(synonyms_file) then "--synonyms ~{synonyms_file}" else ""} \
    ~{if coding_only then "--coding_only" else ""} \
    --dir ~{cache_dir} \
    --fasta ~{reference} \
    ~{sep=" " prefix("--plugin ", plugins)}  \
    --assembly ~{ensembl_assembly} \
    --cache_version ~{ensembl_version} \
    --species ~{ensembl_species} \
    ~{additional_args} $cus_annot

    /usr/local/bin/bgzip ~{annotated_path} && /usr/local/bin/tabix ~{annotated_path}.gz
  >>>

  output {
    File annotated_vcf = "~{annotated_path}.gz"
    File annotated_vcf_tbi = "~{annotated_path}.gz.tbi"
    File vep_summary = annotated_path + "_summary.html"
  }
}


workflow wf {
  input {
    File vcf
    File cache_dir_zip
    File reference
    File reference_fai
    Array[String] plugins
    String ensembl_assembly
    String ensembl_version
    String ensembl_species
    File? synonyms_file
    Array[VepCustomAnnotation] custom_annotations = []
    Boolean coding_only = false
    Boolean everything = true
    String additional_args = "--sift p --polyphen p --pick --pick_order canonical,mane,ccds,rank,appris,tsl,biotype,length --everything 1 --merged --check_existing --buffer_size 1000 --af_gnomad"
  }

  call vep {
    input:
    vcf=vcf,
    cache_dir_zip=cache_dir_zip,
    reference=reference,
    reference_fai=reference_fai,
    plugins=plugins,
    ensembl_assembly=ensembl_assembly,
    ensembl_version=ensembl_version,
    ensembl_species=ensembl_species,
    synonyms_file=synonyms_file,
    custom_annotations=custom_annotations,
    coding_only=coding_only,
    everything=everything,
    additional_args=additional_args
  }
}
