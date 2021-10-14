version 1.0

task markIlluminaAdapters {
    input {
        File bam
    }

    Int cores = 2
    Int space_needed_gb = 10 + round(10*size(bam, "GB"))

    runtime {
        docker: "mgibio/dna-alignment:1.0.0"
        memory: "6GB"
        cpu: cores
        bootDiskSizeGb: 25
        disks: "local-disk ~{space_needed_gb} SSD"
    }

    command <<<
        /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar MarkIlluminaAdapters INPUT=~{bam} OUTPUT=marked.bam METRICS=adapter_metrics.txt
    >>>

    output {
        File marked_bam = "marked.bam"
        File metrics = "adapter_metrics.txt"
    }
}

workflow wf {
    input {
        File bam
    }
    call markIlluminaAdapters {
        input:
        bam = bam
    }
}
