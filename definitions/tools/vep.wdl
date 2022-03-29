version 1.0

import "../types.wdl"

task vepTask {
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
        VepSpliceAIPlugin? spliceAI_files
        Boolean coding_only = false
        Array[VepCustomAnnotation] custom_annotations = []
        Boolean everything = true
        # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
        String pick = "flag_pick"
        String additional_args = "--pick_order canonical,rank,mane,ccds,appris,tsl,biotype,length --merged --buffer_size 1000 --af_gnomad"
        File? synonyms_file
    }

    Float cache_size = 3*size(cache_dir_zip, "GB")  # doubled to unzip
    Float vcf_size = 2*size(vcf, "GB")  # doubled for output vcf
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Int space_needed_gb = 50 + round(reference_size + vcf_size + cache_size + size(synonyms_file, "GB"))

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

    File test = spliceAI_files.spliceAI_snv

  command <<<
    if [[ ~{annotation_len} -ge 1 ]]; then
      custom_annotation=$(/usr/bin/python3 /opt/bin/jsonToVepString.py ~{write_json(custom_annotations)})
    else
      custom_annotation=""
    fi
    echo $custom_annotation

    #mkdir ~{cache_dir} && unzip -qq ~{cache_dir_zip} -d ~{cache_dir}
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
    --sift p \
    --polyphen p \
    ~{if coding_only then "--coding_only" else ""} \
    --~{pick} \
    --dir ~{cache_dir} \
    --fasta ~{reference} \
    ~{sep=" " prefix("--plugin ", plugins)}  \
    ~{if defined(spliceAI_files) then "--plugin SpliceAI,snv=~{test}" else ""} \
    ~{if everything then "--everything" else ""} \
    --assembly ~{ensembl_assembly} \
    --cache_version ~{ensembl_version} \
    --species ~{ensembl_species} \
    ~{additional_args} \
    ${custom_annotation}

    bgzip ~{annotated_path} && tabix ~{annotated_path}.gz
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
    File reference_dict
    Array[String] plugins
    VepSpliceAIPlugin? spliceAI_files
    String ensembl_assembly
    String ensembl_version
    String ensembl_species
    File? synonyms_file
    Array[VepCustomAnnotation] custom_annotations = []
    Boolean coding_only = false
    Boolean everything = true
    String pick = "pick"       # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
  }

  call vepTask {
      input:
          vcf = vcf,
          cache_dir_zip = cache_dir_zip,
          reference = reference,
          reference_fai = reference_fai,
          reference_dict = reference_dict,
          plugins = plugins,
          spliceAI_files = spliceAI_files,
          ensembl_assembly = ensembl_assembly,
          ensembl_version = ensembl_version,
          ensembl_species = ensembl_species,
          synonyms_file = synonyms_file,
          custom_annotations = custom_annotations,
          coding_only = coding_only,
          everything = everything,
          pick = pick
  }

  output {
      File annotated_vcf = vepTask.annotated_vcf
      File annotated_vcf_tbi = vepTask.annotated_vcf_tbi
      File vep_summary = vepTask.vep_summary
  }
}
