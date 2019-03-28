# variant calling
These scripts were applied for variant calling and processing in the context of NAVIP: https://github.com/bpucker/NAVIP. Some scripts are only included for documentation purposes, while others were written in a generic way to facilitate re-use.



# GATK_variant_calling.py

This script is intended as documentation of the process. It is customized for best performance on the local compute cluster. Re-use would require adjustments to certain parts of the script.

python GATK_variant_calling.py \
--input_bam_file <PATH_TO_BAM_FILE> \
--ref_file <PATH_TO_REFERENCE_FILE> \
--directory <PATH_TO_DIRECTORY> \
--piccard <FULL_PATH_TO_PICCARD> \
--samtools <FULL_PATH_TO_SAMTOOLS> \
--gatk <FULL_PATH_TO_GATK> \
--varcallprepscript <FULL_PATH_TO variant_call_preparation.py> \
--varsortscript <FULL_PATH_TO sort_vcf_by_fasta.py> \
--bam_is_sorted (prevents sorting of bam file)



# variant_call_preparation.py

This script is used internaly to allow parallel processing of sequences in the reference data set.



# GATK1_BP.py

This script is intended as documentation of the process. It is customized for best performance on the local compute cluster. Re-use would require adjustments to certain parts of the script.

python GATK1_BP.py \
--input_bam_file <PATH_TO_BAM_FILE> \
--ref_file <PATH_TO_REFERENCE_FILE> \
--directory <PATH_TO_DIRECTORY> \
--gold_vcf <PATH_TO_GOLD_STANDARD_VCF_FILE> \
--piccard <FULL_PATH_TO_PICCARD> \
--samtools <FULL_PATH_TO_SAMTOOLS> \
--gatk <FULL_PATH_TO_GATK> \
--varcallprepscript <FULL_PATH_TO variant_call_preparation.py> \
--bam_is_sorted (prevents sorting of bam file)



# GATK2_BP.py

This script is intended as documentation of the process. It is customized for best performance on the local compute cluster. Re-use would require adjustments to certain parts of the script.

python GATK2_BP.py \
--ref_file <FULL_PATH_TO_REF_FILE> \
--vcf_dir <FULL_PATH_TO_DIR_WITH_PREPARED_VCFs> \
--out_dir <FULL_PATH_TO_OUTPUT_DIR> \
--gold_vcf <FULL_PATH_TO_GOLD_STANDARD_VCF> \
--piccard <FULL_PATH_TO_PICCARD> \
--samtools <FULL_PATH_TO_SAMTOOLS> \
--gatk <FULL_PATH_TO_GATK>



# VCF_combiner.py

python VCF_combiner.py \
--in <INPUT_DIRECTORY> \
--out <OUTPUT_VCF>



# sort_vcf_by_fasta.py

python sort_vcf_by_fasta.py \
--vcf <INPUT_VCF> \
--fasta <INPUT_FASTA_FILE> \
--output <OUTPUT_VCF_FILE>



# variant_validator.py

WARNING: number of sequences (chromosomes) should not exceed 9!


python variant_validator.py \
--assembly <FULL_PATH_TO_ASSEMBLY_FILE> \
--ref <FULL_PATH_TO_REFERENCE_FILE> \
--invcf <FULL_PATH_TO_INPUT_VCF_FILE> \
--flank <INT, size of query sequence> \
--outvcf <FULL_PATH_TO_OUTPUT_VCF> \
--chr <CHROMOSOME_TO_PROCESS> \
--outerr <FULL_PATH_TO_ERROR_OUTPUT_FILE>



# variant_validation_wrapper.py

python variant_validaton_wrapper.py \
--assembly <FULL_PATH_TO_ASSEMBLY_FILE> \
--ref <FULL_PATH_TO_REFERENCE_FILE> \
--vcf <FULL_PATH_TO_INPUT_VCF_FILE> \
--flank <INT, size of query sequence> \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY> \
--script <FULL_PATH_TO variant_validator.py>



# analyze_variant_set.py

python analyze_variant_set.py \
--vcf <FULL_PATH_TO_VCF_FILE (INPUT)> \
--fig  <FULL_PATH_TO_FIGURE_FILE (OUTPUT)> \
--report <FULL_PATH_TO_REPORT_FILE (OUTPUT)>

