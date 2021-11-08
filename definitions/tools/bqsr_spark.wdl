version 1.0

task bqsr_spark {
  input {
    File reference
    File reference_fai   # secondaryFiles...
    File reference_dict  # secondaryFiles...

    File bam
    File bam_bai  # secondaryFiles...

    Array[File] known_sites
    Array[File] known_sites_tbi  # secondaryFiles...

  }

  Int cores = 16
  Int space_needed_gb = 10 + round(size(known_sites, "GB") + size(known_sites_tbi, "GB") + size([reference, reference_fai, reference_dict, bam, bam_bai], "GB"))
  
  runtime {
    docker: "broadinstitute/gatk:4.2.1.0"
    memory: "64GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} SSD"
  }

 
  String bam_name = basename(bam, ".bam")
  command <<<
    /gatk/gatk --java-options "-Xmx124G -XX:+UseParallelGC -XX:ParallelGCThreads=32" BQSRPipelineSpark \
        -R ~{reference} \
        -I ~{bam} \
        ~{sep=" " prefix("--known-sites ", known_sites)} \
        -O test.bqsr.bam --verbosity ERROR \
        -- --spark-runner LOCAL --spark-master local[16] \
        --conf spark.local.dir=/tmp/spark
  >>>

  output {
    File bqsr_bam = "test.bqsr.bam"
    File bqsr_bam_bai = "test.bqsr.bam.bai"
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict

    File bam
    File bam_bai

    Array[File] known_sites
    Array[File] known_sites_tbi
  }

  call bqsr_spark {
    input:
    reference = reference,
    reference_fai = reference_fai,
    reference_dict = reference_dict,
    bam = bam,
    bam_bai = bam_bai,
    known_sites = known_sites,
    known_sites_tbi = known_sites_tbi
  }
}
