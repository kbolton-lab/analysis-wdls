version 1.0

task bcftoolsMerge {
  input {
    Array[File] vcfs
    Array[File] vcf_tbis
    String merged_vcf_basename = "merged"
  }

  Int space_needed_gb = 5 + round(2*(size(vcfs, "GB") + size(vcf_tbis, "GB")))
  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0

  runtime {
    docker: "kboltonlab/bst:latest"
    memory: "6GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
  }

  String output_file = merged_vcf_basename + ".vcf.gz"
  command <<<
    /usr/local/bin/bcftools merge --output-type z -o ~{output_file} ~{sep=" " vcfs}
    /usr/local/bin/tabix ~{output_file}
  >>>

  output {
    File merged_vcf = output_file
    File merged_vcf_tbi = "~{output_file}.tbi"
  }
}
