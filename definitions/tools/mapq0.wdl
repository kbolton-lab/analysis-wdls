version 1.0

task mapq0 {
    input {
        File vcf
        File vcf_tbi
        File bam
        File bam_bai
        Float mapq0perc = 0.15
        String? output_vcf_prefix = "mapq0_filter"
        String caller = "caller"
        String output_type = "z"
    }


    Int cores = 2
    Float bam_size = size([bam, bam_bai], "GB")
    Int space_needed_gb = 10 + round(2*bam_size + size(vcf, "GB"))


    runtime {
        docker: "kboltonlab/bst:latest"
        memory: "12GB"
        cpu: "~{cores}"
        bootDiskSizeGb: "~{space_needed_gb}"
        disks: "local-disk ~{space_needed_gb} SSD"
    }

 
    command <<<
        set -ex -o pipefail
        #grab sites that don't already have the MQ0 field
        # rm mapq0counts
        zgrep -v "^#" ~{vcf} | grep -v "MQ0" | cut -f 1,2 | while read chr pos; do
            /usr/local/bin/samtools stats -d -@8 ~{bam} $chr:$pos-$pos > stats
            mapq0=$(grep "reads MQ0:" stats | cut -f3); printf "$chr\t$pos\t$mapq0\t" >> mapq0counts
            # sequences/reads
            grep -P "SN\tsequences:" stats | cut -f3 >> mapq0counts
        done
      

        printf "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Number of MAPQ == 0 reads covering this record\">\n##INFO=<ID=samtools_DP,Number=1,Type=Integer,Description=\"Samtools depth at this position\">\n" > MQ0.header;
        

        /usr/local/bin/bgzip -f mapq0counts
        /usr/local/bin/tabix mapq0counts.gz -s1 -b2 -e2;
        /usr/local/bin/bcftools annotate --threads 8 -a mapq0counts.gz -h MQ0.header -c CHROM,POS,MQ0,samtools_DP ~{vcf} | /usr/local/bin/bcftools filter -m+ -e "((INFO/MQ0) / (INFO/samtools_DP)) > ~{mapq0perc}" -s "MQ0" --threads 8 -O~{output_type} -o ~{caller}.~{output_vcf_prefix}.vcf.gz && /usr/local/bin/tabix ~{caller}.~{output_vcf_prefix}.vcf.gz

    >>>

    output {
        File filtered_vcf = "~{caller}.~{output_vcf_prefix}.vcf.gz"
        File filtered_vcf_tbi = "~{caller}.~{output_vcf_prefix}.vcf.gz.tbi"
    }
}

workflow wf {
    input {
        File vcf
        File vcf_tbi
        File bam
        File bam_bai
        Float mapq0perc = 0.15
        String? output_vcf_prefix = "mapq0_filter"
        String caller = "caller"
        String output_type = "z"
    }


    call mapq0 {
        input:
            vcf=vcf,
            vcf_tbi=vcf_tbi,
            bam=bam,
            bam_bai=bam_bai,
            mapq0perc=mapq0perc,
            output_vcf_prefix=output_vcf_prefix,
            caller=caller,
            output_type=output_type
    }

}
