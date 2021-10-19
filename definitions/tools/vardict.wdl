version 1.0

task vardict {
    input {
        File reference
        File reference_fai

        File tumor_bam
        File tumor_bam_bai

        File? normal_bam
        File? normal_bam_bai

        File interval_bed
        Float min_var_freq = 0.005
    }

    Int cores = 32
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))

    runtime {
        docker: "kboltonlab/vardictjava:1.0"
        memory: "128GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
    }

    String output_vcf = "vardict"
    command <<<
        set -o pipefail
        set -o errexit

        NORMAL_BAM=~{normal_bam}

        if [ -z ${NORMAL_BAM} ]; then
            TUMOR=$(samtools view -H ~{tumor_bam} | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)

            /opt/VarDictJava/build/install/VarDict/bin/VarDict \
                -U -G ~{reference} \
                -f ~{min_var_freq} \
                -N ${TUMOR} \
                -b "~{tumor_bam}" \
                -c 1 -S 2 -E 3 -g 4 ~{interval_bed} \
                -th 64 | \
            /opt/VarDictJava/build/install/VarDict/bin/teststrandbias.R | \
            /opt/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl \
                -N "$TUMOR" \
                -E \
                -f ~{min_var_freq} > ~{output_vcf}.vcf

            /usr/bin/bgzip ~{output_vcf}.vcf && /usr/bin/tabix ~{output_vcf}.vcf.gz

            /usr/bin/bcftools norm  \
                ~{output_vcf}.vcf.gz \
                -f $REF -m -any --threads 64 -Oz \
                -o ~{output_vcf}.normalized.vcf.gz && \
            /usr/bin/tabix ~{output_vcf}.normalized.vcf.gz 

        else
            NORMAL=$(samtools view -H ~{normal_bam} | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)
            TUMOR=$(samtools view -H ~{tumor_bam} | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)

            /opt/VarDictJava/build/install/VarDict/bin/VarDict \
                -U -G ~{reference} \
                -f ~{min_var_freq} \
                -N ${TUMOR} \
                -b "~{tumor_bam}|~{normal_bam}" \
                -c 1 -S 2 -E 3 -g 4 ~{interval_bed} \
                -th 64 | \
            /opt/VarDictJava/build/install/VarDict/bin/testsomatic.R | \
            /opt/VarDictJava/build/install/VarDict/bin/var2vcf_paired.pl \
                -N "$TUMOR|$NORMAL" \
                -f ~{min_var_freq} > ~{output_vcf}.vcf
            
            /usr/bin/bgzip ~{output_vcf}.vcf && /usr/bin/tabix ~{output_vcf}.vcf.gz
            # extract tumor before bcbio filter or both normal and tumor will be canidates for the filter
            /usr/bin/bcftools view \
                -s ${TUMOR} \
                --threads 64 ~{output_vcf}.vcf.gz | \
            /usr/bin/bcftools norm \
                -f $REF -m -any --threads 64 -Oz -o ~{output_vcf}.normalized.vcf.gz && \
            /usr/bin/tabix ~{output_vcf}.normalized.vcf.gz 
        fi
    >>>

    output {
        File vcf = "~{output_vcf}.normalized.vcf.gz"
        File vcf_tbi = "~{output_vcf}.normalized.vcf.gz.tbi"
    }
}

workflow wf {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        File? normal_bam
        File? normal_bam_bai
        File interval_bed
    }

    call vardict {
        input:
            reference=reference,
            reference_fai=reference_fai,
            tumor_bam=tumor_bam,
            tumor_bam_bai=tumor_bam_bai,
            normal_bam=normal_bam,
            normal_bam_bai=normal_bam_bai,
            interval_bed=interval_bed
  }
}
