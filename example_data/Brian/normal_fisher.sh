#!/bin/bash
vcf=$1
pon=$2
caller=$3
p_val=$4

echo $vcf
echo $pon
echo $caller
echo $p_val


if [[ "$vcf" == *.gz ]]; then
    name=$(basename $vcf .vcf.gz)
    tabix -f $vcf
else
    name=$(basename $vcf .vcf)
    bgzip -f $vcf
    tabix -f $vcf.gz
    vcf=$vcf.gz
fi

if [[ "$pon" == *.gz ]]; then
    tabix -f $pon  
else
    bgzip -f $pon
    tabix -f $pon.gz
    pon=$pon.gz
fi


bcftools +fill-tags -Oz -o RD.vcf.gz $pon -- -t "PON_RefDepth=sum(RD)"
bcftools +fill-tags -Oz -o RD_AD.vcf.gz RD.vcf.gz -- -t "PON_AltDepth=sum(AD)" && tabix -f RD_AD.vcf.gz

printf "##INFO=<ID=PON_FISHER,Number=1,Type=Float,Description=\"P-value from Fisher's exact test with totals from PoN\">" > fisher.header;
printf "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">" > sample.header;

sample=$(bcftools query -l $vcf)
bcftools view -H $vcf | awk -v sampleID=$sample '{print $1, $2, $3, $4, $5, sampleID}' OFS='\t' > $sample.name;
bgzip -f $sample.name;
tabix -f $sample.name.gz -s1 -b2 -e2;
bcftools annotate -a $sample.name.gz -h sample.header -c CHROM,POS,-,REF,ALT,SAMPLE $vcf -Oz -o $name.sample.vcf.gz && tabix -f $name.sample.vcf.gz;

## Varscan has AD and RD instead of comma sep AD field
## you can't double quote for string match in bash like you can in zsh so need to make it a variable
pat="[Vv]arscan"

## Lofreq has DP4 which splits into RefFwd, RefRev, AltFwd, AltRev
patt="[Ll]ofreq"
if [[ ${caller} =~ $pat ]]
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
    stop("Must supply file with 8 columns: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%RD]\t[%AD]", call.=FALSE)
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
    bcftools annotate -a RD_AD.vcf.gz -c PON_RefDepth,PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%RD]\t[%AD]\n' $name.sample.pileup.vcf.gz > $name.fisher.input;
elif [[ ${caller} =~ $patt ]]
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
    stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t%INFO/DP4", call.=FALSE)
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
    bcftools annotate -a RD_AD.vcf.gz -c PON_RefDepth,PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t%INFO/DP4\n' $name.sample.pileup.vcf.gz > $name.fisher.input;
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
    stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%AD]", call.=FALSE)
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

    bcftools annotate -a RD_AD.vcf.gz -c PON_RefDepth,PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%AD]\n' $name.sample.pileup.vcf.gz > $name.fisher.input;

fi
chmod u+x fisherTestInput.R

# Depending on how we split, we might have caller_vcf that doesn't have any variants called
if [ -s $name.fisher.input ]; then
    LC_ALL=C.UTF-8 Rscript --vanilla ./fisherTestInput.R $name.fisher.input $name.fisher.output
    bgzip -f $name.fisher.output
    tabix -f -s1 -b2 -e2 $name.fisher.output.gz
    bcftools annotate -a $name.fisher.output.gz -h fisher.header -c CHROM,POS,REF,ALT,-,-,-,-,PON_FISHER $name.sample.pileup.vcf.gz -Oz -o $name.pileup.fisherPON.vcf.gz && tabix -f $name.pileup.fisherPON.vcf.gz
    bcftools filter -i "INFO/PON_FISHER < ${p_val}" -m+ -s "pon_fail" $name.pileup.fisherPON.vcf.gz -Oz -o $name.filtered.pileup.fisherPON.vcf.gz && tabix -f $name.filtered.pileup.fisherPON.vcf.gz
else
    bcftools annotate -h fisher.header $name.sample.pileup.vcf.gz -Oz -o $name.pileup.fisherPON.vcf.gz && tabix -f $name.pileup.fisherPON.vcf.gz
    bcftools filter -i "INFO/PON_FISHER < ${p_val}" $name.pileup.fisherPON.vcf.gz -Oz -o $name.$caller.filtered.pileup.fisherPON.vcf.gz && tabix -f $name.$caller.filtered.pileup.fisherPON.vcf.gz
fi
