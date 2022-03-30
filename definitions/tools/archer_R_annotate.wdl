version 1.0

task archerRAnnotate {
    input {
        File vcf
        String caller = "caller"
        File bolton_bick_vars
        File mut2_bick
        File mut2_kelly
        File matches2
        File TSG_file
        File oncoKB
        File oncoKB_curated
        File pd_annotation_file
        File truncating
        File cosmic_dir_zip
        String? pon_pvalue = "2.114164905e-6"

    }

    Float file_size = size([vcf, bolton_bick_vars, mut2_bick, mut2_kelly, matches2, TSG_file, oncoKB_curated, pd_annotation_file], "GB")
    Float cosmic_size = 3*size(cosmic_dir_zip, "GB")
    Int space_needed_gb = 10 + round(file_size + cosmic_size)
    Int cores = 8
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/r_docker_ichan:latest"
      memory: "6GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    String cosmic_dir = basename(cosmic_dir_zip, ".zip")

    command <<<
        set -eou pipefail

        if [[ "~{vcf}" == *.gz ]]; then
            name=$(basename ~{vcf} .vcf.gz)
        else
            name=$(basename ~{vcf} .vcf)
        fi

        unzip -qq ~{cosmic_dir_zip}

        LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/ArcherAnnotationScript.R --input ~{vcf} --out ${name} --caller ~{caller} \
        --bolton_bick_vars ~{bolton_bick_vars} \
        --mut2_bick ~{mut2_bick} \
        --mut2_kelly ~{mut2_kelly} \
        --matches2 ~{matches2} \
        --TSG_file ~{TSG_file} \
        --oncoKB_curated ~{oncoKB_curated} \
        --pd_annotation_file ~{pd_annotation_file} \
        --cosmic_dir ~{cosmic_dir} \
        --truncating ~{truncating} \
        --p_value ~{pon_pvalue}
    >>>

    output {
        File vcf_annotate_pd = basename(vcf, ".vcf.gz") + ".tsv"
    }
}

workflow wf {
    input {
        File vcf
        String caller = "caller"
        File bolton_bick_vars
        File mut2_bick
        File mut2_kelly
        File matches2
        File TSG_file
        File oncoKB
        File oncoKB_curated
        File pd_annotation_file
        File truncating
        File cosmic_dir_zip
    }

    call archerRAnnotate {
        input:
            vcf = vcf,
            caller = caller,
            bolton_bick_vars = bolton_bick_vars,
            mut2_bick = mut2_bick,
            mut2_kelly = mut2_kelly,
            matches2 = matches2,
            TSG_file = TSG_file,
            oncoKB = oncoKB,
            oncoKB_curated = oncoKB_curated,
            pd_annotation_file = pd_annotation_file,
            truncating = truncating,
            cosmic_dir_zip = cosmic_dir_zip
    }
}
