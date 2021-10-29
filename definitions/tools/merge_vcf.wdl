version 1.0

task mergeVcf {
  input {
    Array[File] vcfs
    Array[File] vcf_tbis
    String merged_vcf_basename = "merged"
  }

  Int space_needed_gb = 10 + round(2*(size(vcfs, "GB") + size(vcf_tbis, "GB")))
  runtime {
    docker: "kboltonlab/bst:latest"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  String output_file = merged_vcf_basename + ".vcf.gz"
  command <<<
    /usr/local/bin/bcftools concat --allow-overlaps --remove-duplicates --output-type z -o ~{output_file} ~{sep=" " vcfs}
    
    /usr/local/bin/tabix ~{output_file}
  >>>

  output {
    File merged_vcf = output_file
    File merged_vcf_tbi = "~{output_file}.tbi"
  }
}
