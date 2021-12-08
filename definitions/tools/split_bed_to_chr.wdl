version 1.0

task splitBedToChr {
  input {
    File interval_bed
  }

  Int cores = 1
  Int space_needed_gb = 10 + round(size(interval_bed, "GB")*2)
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

    intervals=$(awk '{print $1}' ~{interval_bed} | uniq)
    for chr in ${intervals}; do
        grep -w $chr ~{interval_bed} > ~{basename(interval_bed, ".bed")}_${chr}.bed
    done
  >>>

  output {
    Array[File] split_chr = glob(basename(interval_bed, ".bed")+"_*.bed")
  }
}

workflow wf {
  input {
      File interval_bed
  }

  call splitBedToChr {
    input:
    interval_bed = interval_bed
  }
}
