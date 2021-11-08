${(z)WOMTOOL} validate /Users/brian/test/WDL/analysis-wdls/definitions/tools/vardict.wdl

${(z)WOMTOOL}


${(z)CROMWELL} run \
    -o /Users/brian/test/WDL/analysis-wdls/example_data/Brian/options/vardict_options_home.json \
    -t wdl \
    -i /Users/brian/test/WDL/analysis-wdls/example_data/tools/vardict_tool_home.json \
    /Users/brian/test/WDL/analysis-wdls/definitions/tools/vardict.wdl

    ${(z)CROMWELL} run \
    -t wdl \
    -i /Users/brian/test/WDL/analysis-wdls/example_data/tools/norm_home.json \
    /Users/brian/test/WDL/analysis-wdls/definitions/tools/bcftools_norm.wdl

${(z)CROMWELL} run \
    -o /Users/brian/test/WDL/analysis-wdls/example_data/Brian/options/vep_options_home.json \
    -t wdl \
    -i /Users/brian/test/WDL/analysis-wdls/example_data/tools/vep_brian.json \
    /Users/brian/test/WDL/analysis-wdls/definitions/tools/vep_brian.wdl

${(z)CROMWELL} run \
    -o /Users/brian/test/WDL/analysis-wdls/example_data/Brian/options/mq0_options_home.json \
    -t wdl \
    -i /Users/brian/test/WDL/analysis-wdls/example_data/tools/mapq0_home.json \
    /Users/brian/test/WDL/analysis-wdls/definitions/tools/mapq0.wdl


${(z)CROMWELL} run \
    -o /Users/brian/test/WDL/analysis-wdls/example_data/Brian/options/mq0_options_home.json \
    -t wdl \
    -i /Users/brian/test/WDL/analysis-wdls/example_data/tools/mutect_tool_home.json \
    /Users/brian/test/WDL/analysis-wdls/definitions/tools/mutect.wdl




/opt/VarDictJava/build/install/VarDict/bin/VarDict \
            -U -G /Users/brian/Bolton/CWL_TESTS/GRCh38_full_analysis_set_plus_decoy_hla.chr22.fa \
            -f 0.005 \
            -N TUMOR \
            -b "/Users/brian/Bolton/CWL_TESTS/1006108_23153_0_0_chr22.bam|  \
            -c 1 -S 2 -E 3 -g 4 /Users/brian/Bolton/CWL_TESTS/xgen_plus_spikein.GRCh38.chr22.bed \
            -th 4

             -b "~{tumor_bam}|~{normal_bam}" \