version 1.0

task splitBamIntoChr {
  input {
    File bam
    File? interval_bed
    File? interval_list
    Array[String]? intervals = ["chr1", "chr2", "chr3", "chr4", "chr5","chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
  }

  Int cores = 1
  Int space_needed_gb = 10 + round(size(interval_list, "GB")*2)
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    cpu: cores
    memory: "6GB"
    docker: "kboltonlab/bst:latest"
    disks: "local-disk ~{space_needed_gb} SSD"
    preemptible: preemptible
    maxRetries: maxRetries
  }

  command <<<

    if ~{defined(interval_bed)}; then
        intervals=$(awk '{print $1}' ~{interval_bed} | uniq)
    elif ~{defined(interval_list)}; then
        intervals=$(grep -v '@' ~{interval_list} | awk '{print $1}' | uniq)
    else
        intervals=~{sep=' ' intervals}
    fi
    for chr in ${intervals}; do
        samtools view -b ~{bam} $chr > ~{basename(bam, ".bam")}_${chr}.bam
    done
  >>>

  output {
    Array[File] split_chr = glob(basename(bam, ".bam")+"_*.bam")
  }
}

workflow wf {
  input {
      File bam
      File? interval_bed
      File? interval_list
      Array[String]? intervals = ["chr1", "chr2", "chr3", "chr4", "chr5","chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
  }

  call splitBamIntoChr {
    input:
    bam = bam,
    interval_bed = interval_bed,
    interval_list=interval_list,
    intervals = intervals
  }
}
