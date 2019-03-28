# variant_calling
The scripts were applied for variant calling and processing. Some scripts are only for documentation purposes, while others are written in a generic way to facilitate re-use.


# GATK_variant_calling.py

This script is intended as documentation of the process. It is customized for best performance on the local compute cluster. Re-use would require adjustments to certain parts of the script.

python GATK_variant_calling.py\n
				--input_bam_file <PATH_TO_BAM_FILE>\n
				--ref_file <PATH_TO_REFERENCE_FILE>\n
				--directory <PATH_TO_DIRECTORY>\n
				--piccard <FULL_PATH_TO_PICCARD>
				--samtools <FULL_PATH_TO_SAMTOOLS>
				--gatk <FULL_PATH_TO_GATK>
				--varcallprepscript <FULL_PATH_TO variant_call_preparation.py>
				--varsortscript <FULL_PATH_TO sort_vcf_by_fasta.py>
				--bam_is_sorted (prevents sorting of bam file)


# GATK1_BP.py

This script is intended as documentation of the process. It is customized for best performance on the local compute cluster. Re-use would require adjustments to certain parts of the script.

python GATK1_BP.py\n
				--input_bam_file <PATH_TO_BAM_FILE>\n
				--ref_file <PATH_TO_REFERENCE_FILE>\n
				--directory <PATH_TO_DIRECTORY>\n
				--gold_vcf <PATH_TO_GOLD_STANDARD_VCF_FILE>\n
				--piccard <FULL_PATH_TO_PICCARD>
				--samtools <FULL_PATH_TO_SAMTOOLS>
				--gatk <FULL_PATH_TO_GATK>
				--varcallprepscript <FULL_PATH_TO variant_call_preparation.py>
				--bam_is_sorted (prevents sorting of bam file)



# GATK2_BP.py

This script is intended as documentation of the process. It is customized for best performance on the local compute cluster. Re-use would require adjustments to certain parts of the script.

python GATK2_BP.py\n
	--ref_file <FULL_PATH_TO_REF_FILE>
	--vcf_dir <FULL_PATH_TO_DIR_WITH_PREPARED_VCFs>
	--out_dir <FULL_PATH_TO_OUTPUT_DIR>
	--gold_vcf <FULL_PATH_TO_GOLD_STANDARD_VCF>
	--piccard <FULL_PATH_TO_PICCARD>
	--samtools <FULL_PATH_TO_SAMTOOLS>
	--gatk <FULL_PATH_TO_GATK>


# VCF_combiner.py

python VCF_combiner.py
					--in <INPUT_DIRECTORY>
					--out <OUTPUT_VCF>
















