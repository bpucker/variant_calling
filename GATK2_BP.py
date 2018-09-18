### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###


__usage__ = """
	python GATK2_BP.py\n
	--ref_file <FULL_PATH_TO_REF_FILE>
	--vcf_dir <FULL_PATH_TO_DIR_WITH_PREPARED_VCFs>
	--out_dir <FULL_PATH_TO_OUTPUT_DIR>
	--gold_vcf <FULL_PATH_TO_GOLD_STANDARD_VCF>
	--piccard <FULL_PATH_TO_PICCARD>
	--samtools <FULL_PATH_TO_SAMTOOLS>
	--gatk <FULL_PATH_TO_GATK>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, re, sys, glob
from operator import itemgetter

# --- end of imports --- #

def index_ref_file( ref_fasta_file, piccard_tools, samtools ):
	"""! @brief construct indexed reference file """
	
	dictionary_file = '.'.join( ref_fasta_file.split('.')[:-1] ) + ".dict"
	cmd1 = [ 	"java -Xmx8g -jar ",
				piccard_tools,
				" CreateSequenceDictionary R=",
				ref_fasta_file,
				" O=",
				dictionary_file
			]
	
	cmd2 = [ 	samtools,
				" faidx ",
				ref_fasta_file,
			]
	
	cmd1 = "".join( cmd1 )
	cmd2 = "".join( cmd2 )
	
	os.popen( cmd1 )
	os.popen( cmd2 )


def load_vcf_data( vcf_file ):
	"""! @brief load all data from VCF file """
	
	data = []
	samples = []
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[ :len( "#CHROM" ) ] == "#CHROM":
				samples = line.strip().split('\t')[9:]
			elif line[0] != '#':				
				parts = line.strip().split('\t')
				data.append( { 'chr': parts[0], 'pos': int( parts[1] ), 'line': line } )
			line = f.readline()
	return data, samples


def construct_sorted_output_file( data, samples, outputfile ):
	"""! @brief write all data in sorted way (chr, pos) into output file """
	
	sorted_data = sorted( data, key=itemgetter( 'chr', 'pos', 'line' ) )
	with open( outputfile, "w" ) as out:
		out.write( "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join( samples ) + '\n' )
		for each in sorted_data:
			out.write( each['line'] )


def separate_and_filter_snps_and_indels( prefix, combined_raw_vcf_file, input_ref_file, GATK ):
	"""! @brief separation of SNPs and InDels and separate filtering of both data sets """
	
	snps_only_file = prefix + "raw_snps_only_file.g.vcf"
	indels_only_file = prefix + "raw_indels_only_file.g.vcf"
	
	clean_snps_file = prefix + "clean_snps_file.g.vcf"
	clean_indels_file = prefix + "clean_indels_file.g.vcf"
	
	clean_all_variants_file = prefix + "clean_all_variants_file.g.vcf"
	
	
	# --- extraction of SNPs --- #
	cmd = ''.join( [ 	"java -Xmx8g -jar ",
						GATK,
						" -T SelectVariants -R ",
						input_ref_file,
						" -V ",
						combined_raw_vcf_file,
						" -selectType SNP -o ",
						snps_only_file
					] )
	os.popen( cmd )
	
	# --- filtering of SNPs --- #
	cmd = ''.join( [ 	"java -Xmx6g -jar ",
						GATK,
						" -T VariantFiltration -R ",
						input_ref_file,
						" -V ",
						snps_only_file,
						' --filterExpression "QD < 2.0" --filterName "QD_filter"',
						' --filterExpression "FS > 60.0" --filterName "FS_filter"',
						' --filterExpression "MQ < 40.0" --filterName "MQ_filter"',
						" -o ",
						clean_snps_file
					] )
	os.popen( cmd )
	
	# --- extraction of InDels --- #
	cmd = ''.join( [ 	"java -Xmx8g -jar ",
						GATK,
						" -T SelectVariants -R ",
						input_ref_file,
						" -V ",
						combined_raw_vcf_file,
						" -selectType INDEL -o ",
						indels_only_file
					])
	os.popen( cmd )
	
	# --- filtering of InDels --- #
	cmd = ''.join( [ 	"java -Xmx6g -jar ",
						GATK,
						" -T VariantFiltration -R ",
						input_ref_file,
						" -V ",
						indels_only_file,
						' --filterExpression "QD < 2.0" --filterName "QD_filter"',
						' --filterExpression "FS > 200.0" --filterName "FS_filter"',
						' -o ',
						clean_indels_file
					] )
	os.popen( cmd )
	
	# --- combine all clean variants in one vcf file --- #
	snps, samples  = load_vcf_data( clean_snps_file )
	indels, samples = load_vcf_data( clean_indels_file )
	construct_sorted_output_file( snps+indels, samples, clean_all_variants_file )
	
	
	# --- print final report --- #
	print "\n\nFINAL REPORT\n\ncleaned SNPs are located here: " + clean_snps_file
	print "\ncleaned InDels are located here: " + clean_indels_file
	print "\nall cleaned variants are located here: " + clean_all_variants_file


def joint_genotyping( merged_vcf_file, all_vcf_files, GATK, input_ref_file ):
	"""! @brief merge all VCF files into one """
	
	cmd = [ "java -Xmx32g -jar " + GATK + " -T GenotypeGVCFs -R " + input_ref_file + " -o " + merged_vcf_file ]
	for vcf in all_vcf_files:
		cmd.append( " --variant " + vcf )
	print "".join( cmd ) 
	os.popen( "".join( cmd ) )


def run_SNP_VQSR( GATK, input_ref_file, merged_vcf_file, known_sites_vcf, output_dir ):
	"""! @brief run SNP recalibration calculation """
	
	recal_snp_file = output_dir + "SNP_recalibration.recal"
	snp_tranch_file = output_dir + "recalibrate_SNP.tranches"
	cmd = [ 	"java -Xmx32g -jar " + GATK,
					" -T VariantRecalibrator -R " + input_ref_file,
					" --input " + merged_vcf_file,
					" -resource:dbsnp,known=true,training=true,truth=true,prior=2.0 " + known_sites_vcf,
					" -an DP ",
					" -an QD ",
					" -an FS ",
					" -an SOR ",
					" -an MQ ",
					" -an MQRankSum ",
					" -an ReadPosRankSum ",
					#" -an InbreedingCoeff ",
					" -mode SNP ",
					" -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 ",
					" -recalFile " + recal_snp_file,
					" -tranchesFile " + snp_tranch_file,
					" -rscriptFile " + output_dir + "recalibrate_SNP_plots.R"
				]
	
	os.popen( "".join( cmd ) )
	return recal_snp_file, snp_tranch_file


def run_snp_filtering( GATK, input_ref_file, merged_vcf_file, snp_recal_file, snp_tranch_file, output_dir ):
	"""! @brief recalibration of SNPs """
	
	recalibrated_snps_raw_indels = output_dir + "recalibrated_snps_raw_indels.g.vcf"
	cmd = [ 	"java -Xmx32g -jar " + GATK,
					" -T ApplyRecalibration -R " + input_ref_file,
					" --input " + merged_vcf_file,
					" -mode SNP ",
					" --ts_filter_level 99.0 ",
					" -recalFile " + snp_recal_file,
					" -tranchesFile " + snp_tranch_file,
					" -o " + recalibrated_snps_raw_indels
				]
		
	os.popen( "".join( cmd ) )
	return recalibrated_snps_raw_indels


def run_indel_VQSR( GATK, input_ref_file, recalibrated_snps_raw_indels, known_sites_vcf, output_dir ):
	"""! @brief run InDel recalibration calculation """
	
	recal_indel_file = output_dir + "indel_recalibration.recal"
	indel_tranch_file = output_dir + "recalibrate_indel.tranches"
	cmd = [ 	"java -Xmx32g -jar " + GATK,
					" -T VariantRecalibrator -R " + input_ref_file,
					" --input " + recalibrated_snps_raw_indels,
					" -resource:dbsnp,known=true,training=true,truth=true,prior=2.0 " + known_sites_vcf,
					" -an QD ",
					" -an DP ",
					" -an FS ",
					" -an SOR ",
					" -an MQRankSum ",
					" -an ReadPosRankSum ",
					#" -an InbreedingCoeff ",
					" -mode INDEL ",
					" -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 ",
					" --maxGaussians 4 "
					" -recalFile " + recal_indel_file,
					" -tranchesFile " + indel_tranch_file,
					" -rscriptFile " + output_dir + "recalibrate_INDEL_plots.R"
				]	
				
	os.popen( "".join( cmd ) )
	return recal_indel_file, indel_tranch_file


def run_indel_filtering( GATK, input_ref_file, recalibrated_snps_raw_indels, indel_recal_file, indel_tranch_file, output_dir ):
	"""! @brief recalibration of InDels """
	
	recalibrated_snps_recalibrated_indels = output_dir + "recalibrated_snps_recalibrated_indels.g.vcf"
	cmd = [ 	"java -Xmx32g -jar " + GATK,
					" -T ApplyRecalibration -R " + input_ref_file,
					" --input " + recalibrated_snps_raw_indels,
					" -mode INDEL ",
					" --ts_filter_level 99.0 "
					" -recalFile " + indel_recal_file,
					" -tranchesFile " + indel_tranch_file,
					" -o " + recalibrated_snps_recalibrated_indels
				]
	os.popen( "".join( cmd ) )
	return recalibrated_snps_recalibrated_indels


def main( arguments ):
	
	GATK = arguments[ arguments.index( '--gatk' )+1 ]	#v3.8
	piccard_tools = arguments[ arguments.index( '--piccard' )+1 ]		#v2.5.0
	samtools = arguments[ arguments.index( '--samtools' )+1 ]
	
	input_ref_file = arguments[ arguments.index( '--ref_file' )+1 ]
	input_vcf_dir = arguments[ arguments.index( '--vcf_dir' )+1 ]
	output_dir = arguments[ arguments.index( '--out_dir' )+1 ]
	known_sites_vcf = arguments[ arguments.index( '--gold_vcf' )+1 ]
	
	all_vcf_files = []
	extensions = [ "*.g.vcf", "*/*.g.vcf", "*.vcf", "*/*.vcf", "*.g.VCF" ]
	for each in extensions:
		all_vcf_files += glob.glob( input_vcf_dir + each )
	
	# --- prepare reference file --- #
	index_ref_file( input_ref_file, piccard_tools, samtools )
	
	# --- do combined genotyping --- #
	merged_vcf_file = output_dir + "merged_variants.g.vcf"
	joint_genotyping( merged_vcf_file, all_vcf_files, GATK, input_ref_file )
	
	print "ok"
	# --- recalibrate variant quality scores (VQSR) --- #
	snp_recal_file, snp_tranch_file = run_SNP_VQSR( GATK, input_ref_file, merged_vcf_file, known_sites_vcf, output_dir )
	recalibrated_snps_raw_indels = run_snp_filtering( GATK, input_ref_file, merged_vcf_file, snp_recal_file, snp_tranch_file, output_dir )
	
	indel_recal_file, indel_tranch_file = run_indel_VQSR( GATK, input_ref_file, recalibrated_snps_raw_indels, known_sites_vcf, output_dir )
	recalibrated_snps_recalibrated_indels = run_indel_filtering( GATK, input_ref_file, recalibrated_snps_raw_indels, indel_recal_file, indel_tranch_file, output_dir )
	
	# --- separate and filter variants --- #
	separate_and_filter_snps_and_indels( output_dir, recalibrated_snps_recalibrated_indels, input_ref_file, GATK )


if __name__ == '__main__':
	
	if '--ref_file' in sys.argv and '--vcf_dir' in sys.argv and '--out_dir' in sys.argv and '--gold_vcf' in sys.argv and '--piccard' in sys.argv and '--samtools' in sys.argv and '--gatk' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	print "all done!"

