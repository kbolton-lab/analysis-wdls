version 1.0

task pindelToVcf {
    input {
        File pindel_output_summary
        File reference
        File reference_fai
        File reference_dict
        String? ref_name = "GRCh38DH"
        String? ref_date = "20161216"
        Int? min_supporting_reads = 3
        String? output_name = "pindel.vcf"
        String tumor_sample_name
    }

    Int space_needed_gb = 10 + round(size([reference, reference_fai, reference_dict, pindel_output_summary], "GB"))
    runtime {
      memory: "16GB"
      docker: "mgibio/cle:v1.3.1"
      disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        /usr/bin/pindel2vcf -G -p ~{pindel_output_summary} -r ~{reference} -R ~{ref_name} -e ~{min_supporting_reads} -d ~{ref_date} -v ~{output_name}
        # If pindel returns empty pindel.head file, need to account for empty file.
        is_empty=$(grep "~{tumor_sample_name}" ~{output_name})
        if [[ ${is_empty} == "" ]]; then
            grep "##" ~{output_name} > temp.vcf
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{tumor_sample_name}" >> temp.vcf
            mv temp.vcf ~{output_name}
        fi
    >>>

    output {
        File pindel_vcf = "~{output_name}"
    }
}
