version 1.0

task fastqToBam {
    input {
        File fastq1
        File fastq2
        String sample_name
        String library_name
        String platform_unit
        String platform
    }

    Int cores = 2
    Int space_needed_gb = 10 + round(10*size([fastq1, fastq2], "GB"))

    runtime {
        docker: "mgibio/dna-alignment:1.0.0"
        memory: "6GB"
        cpu: cores
        bootDiskSizeGb: 25
        disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar FastqToSam FASTQ=~{fastq1} FASTQ2=~{fastq2} SAMPLE_NAME=~{sample_name} LIBRARY_NAME=~{library_name} PLATFORM_UNIT=~{platform_unit} PLATFORM=~{platform} OUTPUT=unaligned.bam
    >>>

    output {
        File bam = "unaligned.bam"
    }
}

workflow wf {
    input {
        File fastq1
        File fastq2
        String sample_name
        String library_name
        String platform_unit
        String platform
    }
    call fastqToBam {
        input:
        fastq1 = fastq1,
        fastq2 = fastq2,
        sample_name = sample_name,
        library_name = library_name,
        platform_unit = platform_unit,
        platform = platform
    }
}
