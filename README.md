[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2616418.svg)](https://doi.org/10.5281/zenodo.2616418)

# variant calling
These scripts were applied for variant calling and processing in the context of NAVIP: https://github.com/bpucker/NAVIP. Some scripts are only included for documentation purposes, while others were written in a generic way to facilitate re-use. All scripts are written in Python v2.7.



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

This script combines the content of all VCF files detected in the provided input folder in a single VCF file.

python VCF_combiner.py \
--in <INPUT_DIRECTORY> \
--out <OUTPUT_VCF>



# sort_vcf_by_fasta.py

This script sorts a given VCF file based on the oder of sequences in a given FASTA file.

python sort_vcf_by_fasta.py \
--vcf <INPUT_VCF> \
--fasta <INPUT_FASTA_FILE> \
--output <OUTPUT_VCF_FILE>



# variant_validator.py

This script validates variants in a given VCF file by comparison against a high quality assembly. This assembly needs to be independent from the reads contributing to the analyzed variants.

WARNING: number of sequences (chromosomes) should not exceed 9!


python variant_validator.py \
--assembly <FULL_PATH_TO_ASSEMBLY_FILE> \
--ref <FULL_PATH_TO_REFERENCE_FILE> \
--invcf <FULL_PATH_TO_INPUT_VCF_FILE> \
--flank <INT, SIZE_OF_QUERY_SEQUENCE> \
--outvcf <FULL_PATH_TO_OUTPUT_VCF> \
--chr <CHROMOSOME_TO_PROCESS> \
--outerr <FULL_PATH_TO_ERROR_OUTPUT_FILE>



# variant_validation_wrapper.py

This script splits a given VCF file and allows parallel processing of variants in each sequence.

python variant_validation_wrapper.py \
--assembly <FULL_PATH_TO_ASSEMBLY_FILE> \
--ref <FULL_PATH_TO_REFERENCE_FILE> \
--vcf <FULL_PATH_TO_INPUT_VCF_FILE> \
--flank <INT, SIZE_OF_QUERY_SEQUENCE> \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY> \
--script <FULL_PATH_TO variant_validator.py>



# analyze_variant_set.py

This script calculates statistics and displays the genome-wide distribution of variants.

python analyze_variant_set.py \
--vcf <FULL_PATH_TO_VCF_FILE (INPUT)> \
--fig  <FULL_PATH_TO_FIGURE_FILE (OUTPUT)> \
--report <FULL_PATH_TO_REPORT_FILE (OUTPUT)>


# correct_VCF_format.py
Add a last column (FORMAT) to an existing VCF-like file to meet the VCF requirements.

python3 correct_VCF_format.py
--in <INPUT_VCF>
--out <OUTOUT_VCF>


# separate_SNVs_InDels.py
Separate SNVs and InDels from a VCF file by generating two separate new files.

python3 separate_SNVs_InDels.py
--in <INPUT_VCF_FILE>
--snvout <OUTPUT_SNV_FILE>
--indelout <OUTPUT_INDEL_FILE>


# compare_stop_gain_events.py
This script compares the stop_gain predictions of SnpEff and NAVIP.

python3 compare_stop_gain_events.py.py
--snpeffvcf <SnpEff_VCF_OUTPUT_FILE>
--navipvcf <NAVIP_VCF_OUTPUT_FILE>
--out <OUTPUT_FOLDER>


# dNdS_analysis.py

python dNdS_analysis.py
--in <NAVIP_OUTPUT_FILE>
--genes <GENE_IDs_FILE>
--out <OUTPUT_FOLDER>


# Reference (how to cite):
Baasner, J.-S., Howard, D., Pucker, B.(2019). Influence of neighboring small sequence variants on functional impact prediction. bioRxiv. doi:10.1101/596718
https://doi.org/10.1101/596718

