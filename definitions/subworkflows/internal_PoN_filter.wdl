version 1.0

import "../types.wdl"

import "../tools/msk_get_base_counts.wdl" as mgbc
import "../tools/normal_fisher.wdl" as nf
import "../tools/index_vcf.wdl" as iv
import "../tools/bcftools_merge.wdl" as bm

workflow InternalPoNFilter {
    input {
        File reference
        File reference_fai
        File reference_dict
        File caller_vcf
        Array[bam_and_bai] pon_bams
        Int? mapq = 5
        Int? baseq = 5
        String? pon_final_name = "pon.pileup"
        String caller = "caller"
        String? p_value = "0.05"

    }

    scatter (pon_bam in pon_bams) {
        call mgbc.mskGetBaseCounts as mskGetBaseCounts {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            normal_bam = pon_bam,
            pon_final_name = pon_final_name,
            vcf = caller_vcf,
            mapq = mapq,
            baseq = baseq
        }
    }

    call bm.bcftoolsMerge as merge {
        input:
            vcfs = mskGetBaseCounts.pileup,
            vcf_tbis = mskGetBaseCounts.pileup_tbi
    }

    call normalFisher as fisher_test {
        input:
            vcf = caller_vcf,
            pon = merge.merged_vcf,
            pon_tbi = merge.merged_vcf_tbi,
            caller = caller,
            p_value = p_value
    }

    output {
        File pon_total_counts = merge.merged_vcf
        File pon_vcf = fisher_test.pon_vcf
        File pon_filtered_vcf = fisher_test.pon_filtered_vcf
    }
}

task normalFisher {
    input {
        File vcf
        File pon
        File pon_tbi
        String caller = "caller"
        String? p_value = "0.05"
    }


    Int space_needed_gb = 10 + round(size([vcf, pon, pon_tbi], "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/sam_bcftools_tabix_bgzip:1.0"
      memory: "6GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        set -eou pipefail

        if [[ "~{vcf}" == *.gz ]]; then
            name=$(basename ~{vcf} .vcf.gz)
        else
            name=$(basename ~{vcf} .vcf)
        fi

        bcftools +fill-tags -Oz -o RD.vcf.gz ~{pon} -- -t "Internal_PON_RefDepth=sum(RD)"
        bcftools +fill-tags -Oz -o RD_AD.vcf.gz RD.vcf.gz -- -t "Internal_PON_AltDepth=sum(AD)" && tabix RD_AD.vcf.gz

        printf "##INFO=<ID=Internal_PON_FISHER,Number=1,Type=Float,Description=\"P-value from Fisher's exact test with totals from PoN\">" > fisher.header;
        printf "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">" > sample.header;

        sample=`bcftools query -l ~{vcf}`
        bcftools view -H ~{vcf} | awk -v sampleID=$sample '{print $1, $2, $3, $4, $5, sampleID}' OFS='\t' > $sample.name;
        bgzip $sample.name;
        tabix $sample.name.gz -s1 -b2 -e2;
        bcftools annotate -a $sample.name.gz -h sample.header -c CHROM,POS,-,REF,ALT,SAMPLE ~{vcf} -Oz -o $name.sample.vcf.gz && tabix $name.sample.vcf.gz;

        ## Varscan has AD and RD instead of comma sep AD field
        ## you can't double quote for string match in bash like you can in zsh so need to make it a variable
        pat="[Vv]arscan"

        ## Lofreq has DP4 which splits into RefFwd, RefRev, AltFwd, AltRev
        patt="[Ll]ofreq"
        if [[ ~{caller} =~ $pat ]]
        then
            echo '
            #!/usr/bin/env Rscript

            args = commandArgs(trailingOnly=TRUE)

            if (length(args)==0) {
            stop("At least one argument must be supplied (input file).n", call.=FALSE)
            } else if (length(args)==1) {
            # default output file
            stop("Must supply (output file).n", call.=FALSE)
            }

            df = read.table(args[1], header=F, colClasses = c("character", "numeric", "character", "character", "numeric", "numeric", "numeric", "numeric"))

            if (length(colnames(df)) != 8) {
            stop("Must supply file with 8 columns: %CHROM\t%POS\t%REF\t%ALT\t%INFO/Internal_PON_RefDepth\t%INFO/Internal_PON_AltDepth\t[%RD]\t[%AD]", call.=FALSE)
            }

            df$fisher.exact.pval <- apply(df, 1, function(x) {
            x <- as.numeric(x[-c(1,2,3,4)])
            if ((x[1]+x[2]==0) | (x[3]+x[4]==0)){
                return(0)
            } else if (x[2]==0 & x[1]!=0) {
                return(0)
            } else if ((x[1]==0 & x[2]!=0) | (x[3]==0 & x[4]!=0)) {
                return(1)
            } else if (x[2]/(x[1]+x[2]) >= x[4]/(x[3]+ x[4])) {
                return(1)
            } else {
                return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
            }
            })
            write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
            ' > fisherTestInput.R
            bcftools annotate -a RD_AD.vcf.gz -c Internal_PON_RefDepth,Internal_PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/Internal_PON_RefDepth\t%INFO/Internal_PON_AltDepth\t[%RD]\t[%AD]\n' $name.sample.pileup.vcf.gz > $name.fisher.input;
        elif [[ ~{caller} =~ $patt ]]
        then
            echo '
            #!/usr/bin/env Rscript

            args = commandArgs(trailingOnly=TRUE)

            if (length(args)==0) {
            stop("At least one argument must be supplied (input file).n", call.=FALSE)
            } else if (length(args)==1) {
            # default output file
            stop("Must supply (output file).n", call.=FALSE)
            }

            df = read.table(args[1], header=F, colClasses = c("character", "numeric", "character", "character", "numeric", "numeric", "character"))

            #https://statisticsglobe.com/split-data-frame-variable-into-multiple-columns-in-r
            df <- cbind(df, data.frame(do.call("rbind", strsplit(as.character(df$V7), ",", fixed = TRUE))))[,-7]

            if (length(colnames(df)) < 8) {
            stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALT\t%INFO/Internal_PON_RefDepth\t%INFO/Internal_PON_AltDepth\t%INFO/DP4", call.=FALSE)
            }

            df$fisher.exact.pval <- apply(df, 1, function(x) {
            x <- as.numeric(x[-c(1,2,3,4)])
            # Remember, Lofreq splits DP4 into RefFwd, RefRev and AltFwd, AltRev so technically ref = x[3] + x[4] and alt = x[5] + x[6]
            ref = x[3] + x[4]
            alt = x[5] + x[6]
            if ((x[1]+x[2]==0) | (ref+alt==0)){
                return(0)
            } else if (x[2]==0 & x[1]!=0) {
                return(0)
            } else if ((x[1]==0 & x[2]!=0) | (ref==0 & alt!=0)) {
                return(1)
            } else if (x[2]/(x[1]+x[2]) >= alt/(ref+alt)) {
                return(1)
            } else {
                return(fisher.test(matrix(c(x[1], x[2], ref, alt), ncol=2))$p.value)
            }
            })
            write.table(df[, -c(9:10)], file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
            ' > fisherTestInput.R
            bcftools annotate -a RD_AD.vcf.gz -c Internal_PON_RefDepth,Internal_PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/Internal_PON_RefDepth\t%INFO/Internal_PON_AltDepth\t%INFO/DP4\n' $name.sample.pileup.vcf.gz > $name.fisher.input;
        else
            echo '
            #!/usr/bin/env Rscript

            args = commandArgs(trailingOnly=TRUE)

            if (length(args)==0) {
            stop("At least one argument must be supplied (input file).n", call.=FALSE)
            } else if (length(args)==1) {
            # default output file
            stop("Must supply (output file).n", call.=FALSE)
            }

            df = read.table(args[1], header=F, colClasses = c("character", "numeric", "character", "character", "numeric", "numeric", "character"))

            #https://statisticsglobe.com/split-data-frame-variable-into-multiple-columns-in-r
            df <- cbind(df, data.frame(do.call("rbind", strsplit(as.character(df$V7), ",", fixed = TRUE))))[,-7]
            if (length(colnames(df)) != 8) {
            stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALT\t%INFO/Internal_PON_RefDepth\t%INFO/Internal_PON_AltDepth\t[%AD]", call.=FALSE)
            }

            df$fisher.exact.pval <- apply(df, 1, function(x) {
            x <- as.numeric(x[-c(1,2,3,4)])
            if ((x[1]+x[2]==0) | (x[3]+x[4]==0)){
                return(0)
            } else if (x[2]==0 & x[1]!=0) {
                return(0)
            } else if ((x[1]==0 & x[2]!=0) | (x[3]==0 & x[4]!=0)) {
                return(1)
            } else if (x[2]/(x[1]+x[2]) >= x[4]/(x[3]+x[4])) {
                return(1)
            } else {
                return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
            }
            })
            write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
            ' > fisherTestInput.R
            bcftools annotate -a RD_AD.vcf.gz -c Internal_PON_RefDepth,Internal_PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/Internal_PON_RefDepth\t%INFO/Internal_PON_AltDepth\t[%AD]\n' $name.sample.pileup.vcf.gz > $name.fisher.input;

        fi
        chmod u+x fisherTestInput.R

        # Depending on how we split, we might have caller_vcf that doesn't have any variants called
        if [ -s $name.fisher.input ]; then
            LC_ALL=C.UTF-8 Rscript --vanilla ./fisherTestInput.R $name.fisher.input $name.fisher.output
            bgzip -f $name.fisher.output
            tabix -f -s1 -b2 -e2 $name.fisher.output.gz
            bcftools annotate -a $name.fisher.output.gz -h fisher.header -c CHROM,POS,REF,ALT,-,-,-,-,Internal_PON_FISHER $name.sample.pileup.vcf.gz -Oz -o $name.internalPON.vcf.gz && tabix $name.internalPON.vcf.gz
            bcftools filter -i "INFO/Internal_PON_FISHER<~{p_value}" $name.internalPON.vcf.gz -Oz -o $name.filtered.internalPON.vcf.gz && tabix $name.filtered.internalPON.vcf.gz
        else
            bcftools annotate -h fisher.header $name.sample.pileup.vcf.gz -Oz -o $name.internalPON.vcf.gz && tabix $name.internalPON.vcf.gz
            bcftools filter -i "INFO/Internal_PON_FISHER<~{p_value}" $name.internalPON.vcf.gz -Oz -o $name.filtered.internalPON.vcf.gz && tabix $name.filtered.internalPON.vcf.gz
        fi
    >>>

    output {
        File pon_vcf = select_first([basename(vcf, ".vcf.gz") + ".internalPON.vcf.gz",basename(vcf, ".vcf") + ".internalPON.vcf.gz"])
        File pon_vcf_tbi = select_first([basename(vcf, ".vcf.gz") + ".internalPON.vcf.gz.tbi",basename(vcf, ".vcf") + ".internalPON.vcf.gz.tbi"])
        File pon_filtered_vcf = select_first([basename(vcf, ".vcf.gz") + ".filtered.internalPON.vcf.gz",basename(vcf, ".vcf") + ".filtered.internalPON.vcf.gz"])
        File pon_filtered_vcf_tbi = select_first([basename(vcf, ".vcf.gz") + ".filtered.internalPON.vcf.gz.tbi",basename(vcf, ".vcf") + ".filtered.internalPON.vcf.gz.tbi"])
    }
}
