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

```
Usage
python GATK2_BP.py

Mandatory:
--ref_file           STR   Path to reference sequence file.
--vcf_dir            STR   Path to VCF folder
--out_dir            STR   Path to output folder
--gold_vcf           STR   Path to gold standard VCF
--piccard            STR   Full path to piccard tools.
--samtools           STR   Samtools path.
--gatk               STR   Path to GATK.

Optional:
--bam_is_sorted          (prevents sorting of bam file).
```

`--ref_file` specifies the full path to the reference genome sequence FASTA file.

`--vcf_dir` specifies the folder containing the VCF files.

`--out_dir` specifies the output folder.

`--gold_vcf` specifies the full path to the gold standard VCF.

`--piccard` specifies the full path to piccard tools.

`--samtools` specifies the full path to samtools.

`--gatk` specifies the full path to GATK.


# VCF_combiner.py

This script combines the content of all VCF files detected in the provided input folder in a single VCF file.

```
Usage
python VCF_combiner.py

Mandatory:
--in   STR   Path to VCF input folder.
--out  STR   Path to output file.
```

`--in` specifies the path to the input VCF folder.

`--out` specifies the path to the output VCF file.


# sort_vcf_by_fasta.py

This script sorts a given VCF file based on the oder of sequences in a given FASTA file.

```
Usage
python sort_vcf_by_fasta.py

Mandatory:
--vcf     STR   Path to input VCF.
--fasta   STR   Path to input FASTA.
--output  STR   Path to output VCF.
```

`--vcf` specifies the VCF input file.

`--fasta` specifies the FASTA input file.

`--output` specifies the VCF output file.




# variant_validator.py

This script validates variants in a given VCF file by comparison against a high quality assembly. This assembly needs to be independent from the reads contributing to the analyzed variants.

WARNING: number of sequences (chromosomes) should not exceed 9!

```
Usage
python variant_validator.py

Mandatory:
--assembly   STR   Path to assembly file.
--ref        STR   Path to reference genome sequence file.
--invcf      STR   Path to input VCF file.
--flank      INT   Length of flanking sequences.
--outvcf     STR   Path to output VCF.
--chr        STR   Chromosome name.
--outerr     STR   Path to error output file.
```

`--assembly` specifies the full path to the assembly FASTA file.

`--ref` specifies the full path to the reference genome FASTA file.

`--invcf` specifies the full path to the input VCF file.

`--flank` specifies the size of the flanking sequences of variants to run the validation.

`--outvcf` specifies the full path to the output VCF.

`--chr` specifies the name of a chromsome to run the validation for one chromosome at a time.

`--outerr` specifies the full path to the error output file.



# variant_validation_wrapper.py

This script splits a given VCF file and allows parallel processing of variants in each sequence.

```
Usage
python variant_validation_wrapper.py

Mandatory:
--assembly  STR   Path to assembly file.
--ref       STR   Path to reference file.
--vcf       STR   Path to input VCF file.
--flank     INT   Length of flanking sequences.
--out       STR   Path to the output folder.
--script    STR   Path to variant_validator.py

```

`--assembly` specifies the full path to the assembly FASTA file.

`--ref` specifies the full path to the reference genome sequence FASTA file.

`--vcf` specifies the input VCF file.

`--flank` specifies the length of the variant flanking sequence used for validation.

`--out` specifies the output folder.

`--script` specifies the full path to the script variant_validator.py.


# analyze_variant_set.py

This script calculates statistics and displays the genome-wide distribution of variants.

```
Usage
python analyze_variant_set.py

Mandatory:
--vcf      STR   Path to input VCF file.
--fig      STR   Path to output figure.
--report   STR   Path to report file.
```

`--vcf` specifies the full path to the input VCF file.

`--fig` specifies the full path to the output figure file.

`--report` specifies the full path to the report file.


# correct_VCF_format.py
Add a last column (FORMAT) to an existing VCF-like file to meet the VCF requirements.

```
Usage
python3 correct_VCF_format.py

Mandatory:
--in   STR   Path to input VCF file.
--out  STR   Path to output VCF file.
```

`--in` specifies the full path to the input VCF file.

`--out` specifies the full path to the output VCF file.



# separate_SNVs_InDels.py
Separate SNVs and InDels from a VCF file by generating two separate new files.

```
Usage
python3 separate_SNVs_InDels.py

Mandatory:
--in        STR   Path to input VCF file.
--snvout    STR   Path to SNV output VCF file.
--indelout  STR   Path to InDel output VCF file.
```

`--in` specifies the full path to the input VCF file.

`--snvout` specifies the full path to the SNV output VCF file.

`--indelout` specifies the full path to the InDel output VCF file.




# compare_stop_gain_events.py
This script compares the stop_gain predictions of SnpEff and NAVIP.

```
Usage
python3 compare_stop_gain_events.py.py

Mandatory:
--snpeffvcf  STR   Path to SnpEff output file.
--navipvcf   STR   Path to NAVIP output file.
--out        STR   Path to output folder.
```

`--snpeffvcf` specifies the SnpEff output VCF file that is required as input for this script.

`--navipvcf` specifies the NAVIP output VCF file that is required as input for this script.

`--out` specifies the output folder.


# aa_ns_analysis.py
This script performs an analysis of synonymous (aa<sub>S</sub>) and non-synonymous (aa<sub>N</sub>) variants in genes with premature stop codons.

```
Usage
python aa_ns_analysis.py

Mandatory:
--in     STR   Path to NAVIP output file.
--genes  STR   Path to genes info file.
--out    STR   Path to output folder.
```

`--in` specifies the NAVIP output file as input for this script.

`--genes` specifies the gene info file that provides the IDs of genes with premature stop codons.

`--out` specifies the folder for all output files.



# compare_gene_exp_between_gene_groups.py
This scripts takes the average expression per gene and compares these values between two groups of genes.

```
Usage
python3 compare_gene_exp_between_gene_groups.py

Mandatory:
--genes  STR   Path to genes info file.
--exp    STR   Path to average expression file.
--out    STR   Path to output folder.

optional:
--gff    STR   Path to GFF file.
```

`--genes` specifies the gene info file that provides the IDs of genes with premature stop codons.

`--exp` specifies the path to a file with average gene expression. Gene IDs are in the first column, mean values in the second column, and median values in the third column.

`--out` specifies the folder for all output files.

`--gff` specifies the GFF3 file for background gene IDs.


# Reference (how to cite):
Baasner, J.-S., Howard, D., Pucker, B.(2019). Influence of neighboring small sequence variants on functional impact prediction. bioRxiv. doi:10.1101/596718
https://doi.org/10.1101/596718

