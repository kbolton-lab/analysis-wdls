version 1.0

task XGBModel {
    input {
        File lofreq_tsv
        File mutect_tsv
        File vardict_tsv
        File pindel_full_vcf
        File pon
        String tumor_sample_name
    }

    Float file_size = size([lofreq_tsv, mutect_tsv, vardict_tsv, pindel_full_vcf, pon], "GB")
    Int space_needed_gb = 10 + round(file_size)
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/xgb:latest"
      memory: "6GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        /opt/bin/xgbappcompiled/bin/xgbapp ~{lofreq_tsv} ~{mutect_tsv} ~{vardict_tsv} ~{pindel_full_vcf} ~{pon} ""
        echo "Model Finished..."
    >>>

    output {
        File model_output = "output_~{tumor_sample_name}.tsv.gz"
        File mutect_complex = "mutect_complex_~{tumor_sample_name}.tsv.gz"
        File pindel_complex = "pindel_complex_~{tumor_sample_name}.tsv.gz"
        File lofreq_complex = "lofreq_complex_~{tumor_sample_name}.tsv.gz"
    }
}
