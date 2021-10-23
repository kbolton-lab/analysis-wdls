version 1.0

import "../types.wdl"

task vep {
  input {
    File vcf
    File cache_dir_zip
    File reference
    File reference_fai
    File reference_dict
    String ensembl_assembly
    String ensembl_version
    String ensembl_species
    Array[String] plugins
    Boolean coding_only = false
    #Array[VepCustomAnnotation] custom_annotations = []
    Array[String] custom_annotation_string =[""]
    Array[File] custom_annotation_files = [""]
    Boolean everything = true
    # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
    String pick = "flag_pick"
    File? synonyms_file
  }

  Float cache_size = 2*size(cache_dir_zip, "GB")  # doubled to unzip
  Float vcf_size = 2*size(vcf, "GB")  # doubled for output vcf
  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Int space_needed_gb = 10 + round(reference_size + vcf_size + cache_size + size(synonyms_file, "GB"))
  runtime {
    memory: "64GB"
    bootDiskSizeGb: 30
    cpu: 4
    docker: "mgibio/vep_helper-cwl:vep_101.0_v2"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  String annotated_path = basename(basename(vcf, ".gz"), ".vcf") + "_annotated.vcf"
  String cache_dir = basename(cache_dir_zip, ".zip")

  command <<<
    #mkdir ~{cache_dir} && unzip -qq ~{cache_dir_zip} -d ~{cache_dir}
    unzip -qq ~{cache_dir_zip}

    custom_files=()
    for file_path in ~{sep=" " custom_annotation_files}; do
        custom_files+=("$file_path")
    done

    custom_strings=()
    for custom_string in ~{sep=" " custom_annotation_string}; do
        custom_strings+=("$custom_string")
    done

    for i in $(seq 0 ~{length(custom_annotation_string)}); do
        vep_string="${vep_string} $(sed 's/<CUSTOM_FILE>/${custom_files[$i]}' <<< ${custom_strings[$i]})"
    done


    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
    --format vcf \
    --vcf \
    --fork 4 \
    --term SO \
    --transcript_version \
    --offline \
    --cache \
    --symbol \
    -o ~{annotated_path} \
    -i ~{vcf} \
    ~{if defined(synonyms_file) then "--synonyms ~{synonyms_file}" else ""} \
    --sift p \
    --polyphen p \
    --coding_only ~{coding_only} \
    --~{pick} \
    --dir ~{cache_dir} \
    --fasta ~{reference} \
    ~{sep=" " prefix("--plugin ", plugins)}  \
    ~{if everything then "--everything" else ""} \
    --assembly ~{ensembl_assembly} \
    --cache_version ~{ensembl_version} \
    --species ~{ensembl_species} \
    ${vep_string}
    #~{sep=" " custom_annotation_string}
  >>>

  output {
    File annotated_vcf = annotated_path
    File vep_summary = annotated_path + "_summary.html"
  }
}

task generateCustomString {
    input { VepCustomAnnotation custom_annotation }
    runtime { docker: "ubuntu:xenial" }
    command <<<
        /bin/echo '~{if custom_annotation.annotation.check_existing then "--check_existing" else ""} --custom <CUSTOM_FILE>,~{custom_annotation.annotation.name},~{custom_annotation.annotation.data_format},~{custom_annotation.method},~{if custom_annotation.force_report_coordinates then 1 else 0},~{sep="," custom_annotation.annotation.vcf_fields}'
    >>>
    output {
        String custom_string = read_string(stdout())
        File custom_file = custom_annotation.annotation.file
    }
}

workflow wf {
  input {
    File vcf
    File cache_dir_zip
    File reference
    File reference_fai
    File reference_dict
    Array[String] plugins
    String ensembl_assembly
    String ensembl_version
    String ensembl_species
    File? synonyms_file
    Array[VepCustomAnnotation] custom_annotations = []
    Boolean coding_only = false
    Boolean everything = true
    # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
    String pick = "flag_pick"
  }

  scatter(custom_annotation in custom_annotations) {
      call generateCustomString {
          input: custom_annotation = custom_annotation
      }
  }

  call vep {
    input:
    vcf=vcf,
    cache_dir_zip=cache_dir_zip,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    plugins=plugins,
    ensembl_assembly=ensembl_assembly,
    ensembl_version=ensembl_version,
    ensembl_species=ensembl_species,
    synonyms_file=synonyms_file,
    custom_annotation_string = generateCustomString.custom_string,
    custom_annotation_files = generateCustomString.custom_file,
    coding_only=coding_only,
    everything=everything,
    pick=pick
  }
}
