bsub -G compute-bolton -g /bwileytest3 -oo vardict_test.log -eo vardict_error.log -q general -M 4G -R 'select[mem>4G] span[hosts=1] rusage[mem=4G]' -a "docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-34)" /usr/bin/java \
    -Dconfig.file=/storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/young/test/cromwell.storage1.config \
    -jar /opt/cromwell.jar run \
    -o /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/Brian/options/vardict_options.json \
    -t wdl \
    -i /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/tools/vardict_tool.json \
    /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/tools/vardict.wdl

/storage1/fs1/bolton/Active/data/presets/cromwell.config
/storage1/fs1/bolton/Active/data/presets/cromwell.config


/usr/bin/java \
    -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config \
    -jar /opt/cromwell.jar run \
    -o /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/Brian/options/vardict_options.json \
    -t wdl \
    -i /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/tools/vardict_tool.json \
    /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/tools/vardict.wdl



/opt/java/openjdk/bin/java \
    -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config \
    -jar /app/cromwell.jar run \
    -o /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/Brian/options/vardict_options.json \
    -t wdl \
    -i /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/tools/vardict_tool.json \
    /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/subworkflows/vardict.wdl


/opt/java/openjdk/bin/java \
    -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config \
    -jar /app/cromwell.jar run \
    -o /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/Brian/options/varscan_options.json \
    -t wdl \
    -i /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/tools/varscan_tool.json \
    /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/tools/varscan.wdl  

/opt/java/openjdk/bin/java \
    -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config \
    -jar /app/cromwell.jar run \
    -o /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/Brian/options/varscan_options.json \
    -t wdl \
    -i /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/tools/mapq0.json \
    /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/tools/mapq0.wdl

/usr/bin/java \
    -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config \
    -jar /opt/cromwell.jar run \
    -o /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/Brian/options/varscan_options.json \
    -t wdl \
    -i /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/subworkflows/gnomad_and_MAPQ0_filter.json \
    /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/subworkflows/gnomad_and_MAPQ0_filter.wdl

/opt/java/openjdk/bin/java \
    -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config \
    -jar /app/cromwell.jar run \
    -o /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/Brian/options/varscan_options.json \
    -t wdl \
    -i /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/subworkflows/varscan_sub.json \
    /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/subworkflows/varscan.wdl  

/opt/java/openjdk/bin/java \
    -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config \
    -jar /app/cromwell.jar run \
    -o /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/Brian/options/mutect_normal_options.json \
    -t wdl \
    -i /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/subworkflows/mutect_normal.json \
    /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/subworkflows/mutect_test.wdl  

/usr/bin/java \
    -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config \
    -jar /opt/cromwell.jar run \
    -o /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/Brian/options/vep_options.json \
    -t wdl \
    -i /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/tools/vep_brian.json \
    /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/tools/vep_brian.wdl

broadinstitute/cromwell:dev


/opt/java/openjdk/bin/java \
    -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config \
    -jar /app/cromwell.jar run \
    -o /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/Brian/options/vep_options.json \
    -t wdl \
    -i /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/tools/vep_brian.json \
    /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/tools/vep_brian.wdl


/opt/java/openjdk/bin/java \
    -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config \
    -jar /app/cromwell.jar run \
    -o /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/Brian/options/vep_options.json \
    -t wdl \
    -i /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/tools/vep_brian_no_custom.json \
    /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/tools/vep_brian.wdl







for sample in my_samples; do
    cd /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/results/$sample
    bsub -oo $sample.vardict.log -G compute-bolton -g /bwileytest -q general-interactive -M 8G -R 'select[mem>8G] rusage[mem=8G]' -a 'docker(broadinstitute/cromwell:dev)' /opt/java/openjdk/bin/java -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config -jar /app/cromwell.jar run -o twinstrand_template.json -t wdl -i /full/dir/$sample/$sample.input.json /wdl/definitions/pipelines/acherdx.wdl
done

for sample in $(ls -d * | grep DNA | tail -n+2); do
    cd /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/results/$sample
    cp somatic_tumor_only_template.json somatic_tumor_only_template_vardict.json
    sed -i "s|/storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/$sample.kraken_filtered.ends_trimmed.rehead.bam|/storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/results/$sample/final.bam|" somatic_tumor_only_template_vardict.json
    bsub -oo $sample.vardict.log -G compute-bolton -g /bwileytest -q general-interactive -M 8G -R 'select[mem>8G] rusage[mem=8G]' -a 'docker(broadinstitute/cromwell:dev)' /opt/java/openjdk/bin/java -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config -jar /app/cromwell.jar run -o twinstrand_template.json -t wdl -i somatic_tumor_only_template_vardict.json /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/pipelines/somatic_tumor_aligned_vardict.wdl
done

for sample in $(ls -d * | grep DNA); do
    cp $sample/vardict.$sample.mapq0_filter_annotated.vcf.gz vardict && gunzip -f vardict/vardict.$sample.mapq0_filter_annotated.vcf.gz
done


/opt/java/openjdk/bin/java -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config -jar /app/cromwell.jar run -t wdl -i /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/tools/fp_filter.json /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/tools/fp_filter.wdl

/usr/bin/perl /usr/bin/fpfilter.pl --bam-readcount /usr/bin/bam-readcount --samtools /opt/samtools/bin/samtools --output fpfilter.vcf --reference /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/GRCh38_full_analysis_set_plus_decoy_hla.fa --bam-file cromwell-executions/wf/03ae9e77-5d85-4ce0-9936-d4d609072384/call-fpFilter/inputs/208882098/1385418_23153_0_0.bqsr.bam --vcf-file cromwell-executions/wf/03ae9e77-5d85-4ce0-9936-d4d609072384/call-fpFilter/inputs/1948159467/1385418_23153_0_0.filtered.vcf.gz --sample UKB_5714275_236171340 --min-var-freq 0.05

/usr/bin/perl /usr/bin/fpfilter.pl --bam-readcount /usr/bin/bam-readcount --samtools /opt/samtools/bin/samtools --output fpfilter.vcf --reference /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/GRCh38_full_analysis_set_plus_decoy_hla.fa --bam-file cromwell-executions/wf/03ae9e77-5d85-4ce0-9936-d4d609072384/call-fpFilter/inputs/208882098/1385418_23153_0_0.bqsr.bam --vcf-file cromwell-executions/wf/03ae9e77-5d85-4ce0-9936-d4d609072384/call-fpFilter/inputs/1948159467/1385418_23153_0_0.filtered.vcf.gz --sample UKB_5714275_236171340 --min-var-freq 0.05


for sample in $(cat /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/results/samples.txt | tail -n+2); do
    cd /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/results/$sample
    bsub -oo $sample.pindel.log -G compute-bolton -g /bwileytest -q general -M 8G -R 'select[mem>8G] rusage[mem=8G]' -a 'docker(broadinstitute/cromwell:dev)' /opt/java/openjdk/bin/java -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config -jar /app/cromwell.jar run -o twinstrand_template.json -t wdl -i somatic_tumor_only_template.pindel.json /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/pipelines/somatic_tumor_aligned_pindel.wdl
done