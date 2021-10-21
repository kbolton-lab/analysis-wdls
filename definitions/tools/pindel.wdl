version 1.0

task pindelNormal {
    input {
        File reference
        File reference_fai
        File reference_dict
        File tumor_bam
        File tumor_bam_bai
        File? normal_bam
        File? normal_bam_bai
        File region_file
        String tumor_sample_name
        String? normal_sample_name
        String? chromosome
        Int insert_size = 400
    }

    Int cores = 4
    Int space_needed_gb = 10 + round(size([reference, reference_fai, reference_dict, normal_bam, normal_bam_bai, tumor_bam, tumor_bam_bai, region_file], "GB"))
    runtime {
        bootDiskSizeGb: 100
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        docker: "mgibio/cle:v1.4.2"
        memory: "16GB"
    }

    command <<<
        echo -e "~{normal_bam}\t~{insert_size}\t~{normal_sample_name}" > pindel.config
        echo -e "~{tumor_bam}\t~{insert_size}\t~{tumor_sample_name}" >> pindel.config

        /usr/bin/pindel -i pindel.config -w 30 -T ~{cores} -o all -f ~{reference} \
        ~{if defined(chromosome) then "-c ~{chromosome}" else ""} \
        ~{if defined(region_file) then "-j ~{region_file}" else ""}
    >>>

  # TODO: how much space to allocate?
    output {
        File deletions = "all_D"
        File insertions = "all_SI"
        File tandems = "all_TD"
        File long_insertions = "all_LI"
        File inversions = "all_INV"
    }
}

task pindelTumorOnly {
    input {
        File reference
        File reference_fai
        File reference_dict
        File tumor_bam
        File tumor_bam_bai
        File region_file
        String tumor_sample_name
        String? chromosome
        Int insert_size = 400
    }

    Int cores = 4
    Int space_needed_gb = 10 + round(size([reference, reference_fai, reference_dict, tumor_bam, tumor_bam_bai, region_file], "GB"))
    runtime {
        bootDiskSizeGb: 100
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        docker: "mgibio/cle:v1.4.2"
        memory: "16GB"
    }

    command <<<
        echo -e "~{tumor_bam}\t~{insert_size}\t~{tumor_sample_name}" >> pindel.config

        /usr/bin/pindel -i pindel.config -w 30 -T ~{cores} -o all -f ~{reference} \
        ~{if defined(chromosome) then "-c ~{chromosome}" else ""} \
        ~{if defined(region_file) then "-j ~{region_file}" else ""}
    >>>

  # TODO: how much space to allocate?
    output {
        File deletions = "all_D"
        File insertions = "all_SI"
        File tandems = "all_TD"
        File long_insertions = "all_LI"
        File inversions = "all_INV"
    }
}

workflow pindel {
    input {
        File reference
        File reference_fai
        File reference_dict
        File tumor_bam
        File tumor_bam_bai
        File? normal_bam
        File? normal_bam_bai
        File region_file
        String tumor_sample_name
        String? normal_sample_name
        String? chromosome
        Int insert_size = 400
        Boolean tumor_only = false
    }

    if (tumor_only) {
        call pindelTumorOnly {
            input:
                reference=reference,
                reference_fai=reference_fai,
                reference_dict=reference_dict,
                tumor_bam=tumor_bam,
                tumor_bam_bai=tumor_bam_bai,
                region_file=region_file,
                tumor_sample_name=tumor_sample_name,
                chromosome=chromosome,
                insert_size=insert_size
        }
    }

    if (!tumor_only) {
        call pindelNormal {
            input:
                reference=reference,
                reference_fai=reference_fai,
                reference_dict=reference_dict,
                tumor_bam=tumor_bam,
                tumor_bam_bai=tumor_bam_bai,
                normal_bam=normal_bam,
                normal_bam_bai=normal_bam_bai,
                region_file=region_file,
                tumor_sample_name=tumor_sample_name,
                normal_sample_name=normal_sample_name,
                chromosome=chromosome,
                insert_size=insert_size
        }
    }

    output {
        File deletions = select_first([pindelNormal.deletions, pindelTumorOnly.deletions])
        File insertions = select_first([pindelNormal.insertions, pindelTumorOnly.insertions])
        File tandems = select_first([pindelNormal.tandems, pindelTumorOnly.tandems])
        File long_insertions = select_first([pindelNormal.long_insertions, pindelTumorOnly.long_insertions])
        File inversions = select_first([pindelNormal.inversions, pindelTumorOnly.inversions])
    }
}
