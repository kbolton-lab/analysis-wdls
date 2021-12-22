#!/bin/bash
vcf=$1
echo "hello world"
pon=$2
caller=$3
p_val=$4

echo $p_val
touch file1.txt
pwd

if [[ "$vcf" == *.gz ]]; then
    name=$(basename $vcf .vcf.gz)
    bgzip -f $vcf
else
    name=$(basename $vcf .vcf)
fi
touch file2.txt
if [[ "$pon" == *.gz ]]; then
    bgzip -f $pon
fi

tabix -f $vcf.gz
tabix -f $pon.gz
bcftools +fill-tags -Oz -o RD.vcf.gz $pon.gz -- -t "PON_RefDepth=sum(RD)"
bcftools +fill-tags -Oz -o RD_AD.vcf.gz RD.vcf.gz -- -t "PON_AltDepth=sum(AD)" && tabix RD_AD.vcf.gz
touch file3.txt
printf "##INFO=<ID=PON_FISHER,Number=1,Type=Float,Description=\"P-value from Fisher's exact test with totals from PoN\">" > fisher.header;
printf "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">" > sample.header;

sample=$(bcftools query -l $vcf.gz)
bcftools view -H $vcf.gz | awk -v sampleID=$sample '{print $1, $2, $3, $4, $5, sampleID}' OFS='\t' > $sample.name;
bgzip -f $sample.name;
tabix $sample.name.gz -s1 -b2 -e2;
bcftools annotate -a $sample.name.gz -h sample.header -c CHROM,POS,-,REF,ALT,SAMPLE $vcf.gz -Oz -o $name.sample.vcf.gz && tabix $name.sample.vcf.gz;
touch file3.txt
## Varscan has AD and RD instead of comma sep AD field
## you can't double quote for string match in bash like you can in zsh so need to make it a variable
pat="[Vv]arscan"

## Lofreq has DP4 which splits into RefFwd, RefRev, AltFwd, AltRev
patt="[Ll]ofreq"