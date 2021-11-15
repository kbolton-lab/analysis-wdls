version 1.0

task vcfSanitize {
  input {
    File vcf
  }

  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    memory: "6GB"
    cpu: cores
    docker: "mgibio/samtools-cwl:1.0.0"
    preemptible: preemptible
    maxRetries: maxRetries
  }

  # outbase should match in script but I don't want to risk changing it yet
  String outbase = basename(basename(vcf, ".gz"), ".vcf")
  command <<<
    set -eou pipefail

    # 1) removes lines containing non ACTGN bases, as they conflict with the VCF spec
    # and cause GATK to choke
    # 2) removes mutect-specific format tags containing underscores, which are likewise
    # illegal in the vcf spec
    base=`basename ~{vcf}`
    outbase=`echo $base | perl -pe 's/.vcf(.gz)?$//g'`
    echo "~{vcf}   $base    $outbase"
    if [[ "~{vcf}" =~ ".gz" ]];then
        #gzipped input
        gunzip -c "~{vcf}" | perl -a -F'\t' -ne 'print $_ if $_ =~ /^#/ || $F[3] !~ /[^ACTGNactgn]/' | sed -e "s/ALT_F1R2/ALTF1R2/g;s/ALT_F2R1/ALTF2R1/g;s/REF_F1R2/REFF1R2/g;s/REF_F2R1/REFF2R1/g" >$outbase.sanitized.vcf
    else
        #non-gzipped input
        cat "~{vcf}" | perl -a -F'\t' -ne 'print $_ if $_ =~ /^#/ || $F[3] !~ /[^ACTGNactgn]/' | sed -e "s/ALT_F1R2/ALTF1R2/g;s/ALT_F2R1/ALTF2R1/g;s/REF_F1R2/REFF1R2/g;s/REF_F2R1/REFF2R1/g" >$outbase.sanitized.vcf
    fi
    /opt/htslib/bin/bgzip $outbase.sanitized.vcf
    /usr/bin/tabix -p vcf $outbase.sanitized.vcf.gz
  >>>


  output {
    File sanitized_vcf = outbase + ".sanitized.vcf.gz"
    File sanitized_vcf_tbi = outbase + ".sanitized.vcf.gz.tbi"
  }
}
