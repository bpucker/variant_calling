[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2616418.svg)](https://doi.org/10.5281/zenodo.2616418)

# variant calling
These scripts were applied for variant calling and processing in the context of NAVIP: https://github.com/bpucker/NAVIP. Some scripts are only included for documentation purposes, while others were written in a generic way to facilitate re-use. Scripts are written in Python v2.7 or Python v3.8.


# GATK_variant_calling.py

This script is intended as documentation of the process. It is customized for best performance on the local compute cluster. Re-use would require adjustments to certain parts of the script.


```
Usage
python GATK_variant_calling.py

Mandatory:
--input_bam_file     STR   Path to BAM file.
--ref_file           STR   Path to reference sequence file.
--directory          STR   Output folder
--piccard            STR   Full path to piccard tools.
--samtools           STR   Samtools path.
--gatk               STR   Path to GATK.
--varcallprepscript  STR   Path to variant_call_preparation.py.
--varsortscript      STR   Path to sort_vcf_by_fasta.py.

Optional:
--bam_is_sorted          (prevents sorting of bam file).

```

`--input_bam_file` specifies full path to BAM input file.

`--ref_file` specifies the full path to the reference genome sequence FASTA file.

`--directory` specifies the output folder.

`--piccard` specifies the full path to piccard tools.

`--samtools` specifies the full path to samtools.

`--gatk` specifies the full path to GATK.

`--varcallprepscript` specifies the full path to the Python script variant_call_preparation.py (see below).

`--varsortscript` specifies the full path to the Python script sort_vcf_by_fasta.py (see below).


# variant_call_preparation.py

This script is used internaly to allow parallel processing of sequences in the reference data set.



# GATK1_BP.py

This script is intended as documentation of the process. It is customized for best performance on the local compute cluster. Re-use would require adjustments to certain parts of the script.


```
Usage
python GATK1_BP.py

Mandatory:
--input_bam_file     STR   Path to BAM file.
--ref_file           STR   Path to reference sequence file.
--directory          STR   Output folder
--gold_vcf           STR   Path to gold standard VCF
--piccard            STR   Full path to piccard tools.
--samtools           STR   Samtools path.
--gatk               STR   Path to GATK.
--varcallprepscript  STR   Path to variant_call_preparation.py.

Optional:
--bam_is_sorted          (prevents sorting of bam file).

```

`--input_bam_file` specifies full path to BAM input file.

`--ref_file` specifies the full path to the reference genome sequence FASTA file.

`--directory` specifies the output folder.

`--gold_vcf` specifies the full path to the gold standard VCF.

`--piccard` specifies the full path to piccard tools.

`--samtools` specifies the full path to samtools.

`--gatk` specifies the full path to GATK.

`--varcallprepscript` specifies the full path to the Python script variant_call_preparation.py (see below).



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
This script performs an analysis of synonymous (dS) and non-synonymous (dN) variants in genes with premature stop codons.


python dNdS_analysis.py
--in <NAVIP_OUTPUT_FILE>
--genes <GENE_IDs_FILE>
--out <OUTPUT_FOLDER>

# compare_gene_exp_between_gene_groups.py
This scripts takes the average expression per gene and compares these values between two groups of genes.

python3 compare_gene_exp_between_gene_groups.py
--genes <SnpEff_VCF_OUTPUT_FILE>
--exp <NAVIP_VCF_OUTPUT_FILE>
--out <OUTPUT_FIGURE_FILE>
optional:
--gff <NAVIP_VCF_OUTPUT_FILE>


# Reference (how to cite):
Baasner, J.-S., Howard, D., Pucker, B.(2019). Influence of neighboring small sequence variants on functional impact prediction. bioRxiv. doi:10.1101/596718
https://doi.org/10.1101/596718

