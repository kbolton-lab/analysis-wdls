version 1.0

import "../tools/extract_umis.wdl" as eu #test
import "../tools/mark_illumina_adapters.wdl" as mia
import "../tools/umi_align.wdl" as ua

workflow umiAlignment {
    input {
        File bam
        Array[String] read_structure
        File reference
        File reference_fai
        File reference_dict
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_pac
        File reference_sa
        Boolean? umi_paired = true
    }

    call eu.extractUmis as extract_umis {
        input:
        bam = bam,
        read_structure = read_structure,
        umi_paired = umi_paired
    }

    call mia.markIlluminaAdapters as mark_illumina_adapters {
        input:
        bam = extract_umis.umi_extracted_bam
    }

    call ua.umiAlign as align {
        input:
        bam = mark_illumina_adapters.marked_bam,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        reference_amb = reference_amb,
        reference_ann = reference_ann,
        reference_bwt = reference_bwt,
        reference_pac = reference_pac,
        reference_sa = reference_sa
    }

    output {
        File aligned_bam = align.aligned_bam
        File aligned_bam_bai = align.aligned_bam_bai
        File adapter_metrics = mark_illumina_adapters.metrics
    }
}
