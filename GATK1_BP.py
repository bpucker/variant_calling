# Boas Pucker #
# bpucker@cebitec.uni-bielefeld.de #
# v0.7 #


import os, sys, datetime, time, re
from operator import itemgetter

# --- end of imports --- #

__usage__ = """ python GATK_variant_calling.py\n
				--input_bam_file <PATH_TO_BAM_FILE>\n
				--ref_file <PATH_TO_REFERENCE_FILE>\n
				--directory <PATH_TO_DIRECTORY>\n
				--gold_vcf <PATH_TO_GOLD_STANDARD_VCF_FILE>\n
				--piccard <FULL_PATH_TO_PICCARD>
				--samtools <FULL_PATH_TO_SAMTOOLS>
				--gatk <FULL_PATH_TO_GATK>
				--varcallprepscript <FULL_PATH_TO variant_call_preparation.py>
				--bam_is_sorted (prevents sorting of bam file)
				
				bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
			"""


def load_sequences( filename ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( filename ) as f:
		header = f.readline()[1:].strip().split(' ')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split(' ')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	print "number of identified sequences: " + str( len( sequences.keys() ) )
	return sequences


def sort_bam_file( prefix, input_bam_file, samtools_path, piccard_tools ):
	"""! @brief sort given bam file by position """
	
	# --- sort bam file --- #
	
	tmp_file = prefix + input_bam_file.split('/')[-1] + "_tmp.bam"
	
	cmd = samtools_path + " sort -@ 8 -m 8G -o " +  tmp_file + " " + input_bam_file
	os.popen( cmd )
	
	# --- adding read group IDs --- #
	output_file = prefix + input_bam_file.split('/')[-1] + "_sorted.bam"
	print "adding read group ... "
	cmd = [ "java -Xmx8g -jar ",
			piccard_tools,
			" AddOrReplaceReadGroups I=",
			tmp_file,
			" O=",
			output_file,
			" RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"
		]
	os.popen( "".join( cmd ) )
	
	return output_file


def index_bam_file( sorted_bam_file, piccard_tools ):
	"""! @brief index given BAM file """
	
	index_file = sorted_bam_file + ".bai"
	
	cmd = [ "java -Xmx8g -jar ",
			piccard_tools,
			" BuildBamIndex I=",
			sorted_bam_file,
			" O=",
			index_file
		 ]
	cmd = "".join( cmd )
	os.popen( cmd )
	return index_file


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


def prepare_files_for_variant_calling_on_cluster( prefix, input_bam_file, headers, raw_bam_files, input_ref_file, samtools, piccard_tools, GATK, para_jobs, variant_call_preparator ):
	"""! @brief do all preparation in one file per chromosome on cluster """
	
	IDs_to_check = []
	files_to_delete = []
	
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, header in enumerate( headers ):
		
		ID = "PP" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + '.sh'
		files_to_delete.append( sh_file )
		out_file = prefix + ID + '.out'
		files_to_delete.append( out_file )
		err_file = prefix + ID + '.report'
		
		cmd = ''.join( [ 	"python " + variant_call_preparator,
							" --header " + header,
							' --original_bam_file ' + input_bam_file,
							" --bam_file ",
							raw_bam_files[ idx ],
							" --ref_file ",
							input_ref_file,
							" --samtools " + samtools,
							" --GATK " + GATK,
							" --piccard_tools " + piccard_tools
						] )
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + '"',
																"| qsub -cwd",
																"-N",
																ID,
																"-l vf=30G",
																"-l arch=lx-amd64",
																"-l idle=1",
																"-P fair_share",
																"-o",
																out_file,
																"-e",
																err_file
															] ) + '\n'
							)
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "PP" + batch_ID + "_\d{4}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 10 )
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "PP" + batch_ID + "_\d{4}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 10 )
	
	# --- delete temp files --- #
	time.sleep( 60 )
	for each in files_to_delete:
		try:
			os.remove( each )
		except:
			"file cannot be deleted, because it is missing"


def recalibration_round( GATK, input_ref_file, duplicated_marked_bam_files, known_sites_vcf, prefix, para_jobs ):
	"""! @brief calculate base recallibration """
	
	IDs_to_check = []
	files_to_delete = []
	recal_files = []
	
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, bam_file in enumerate( duplicated_marked_bam_files ):		
		ID = "RC" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + '.sh'
		files_to_delete.append( sh_file )
		out_file = prefix + ID + '.out'
		files_to_delete.append( out_file )
		err_file = prefix + ID + '.report'
		
		recal_file = bam_file + "recal_data.table"
		recal_files.append( recal_file )
		
		cmd = "java -Xmx32g -jar " + GATK + " -T BaseRecalibrator -R " + input_ref_file + " -I " + bam_file + " -knownSites " + known_sites_vcf + " -nct 4 -U ALLOW_N_CIGAR_READS " + " -o " + recal_file
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + '"',
																"| qsub -cwd",
																"-N",
																ID,
																"-l vf=10G",
																"-l arch=lx-amd64",
																"-l idle=1",
																'-pe multislot 4',
																"-P fair_share",
																"-o",
																out_file,
																"-e",
																err_file
															] ) + '\n'
							)
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "RC" + batch_ID + "_\d{4}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 10 )
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "RC" + batch_ID + "_\d{4}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 10 )
	
	# --- delete temp files --- #
	for each in files_to_delete:
		try:
			os.remove( each )
		except:
			"file cannot be deleted, because it is missing"
	
	return recal_files


def base_recalibration( GATK, input_ref_file, duplicated_marked_bam_files, recal_files, prefix, para_jobs ):
	"""! @brief runs actual recalibration """
	
	IDs_to_check = []
	files_to_delete = []
	resulting_bam_files = []
	
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, bam_file in enumerate( duplicated_marked_bam_files ):		
		ID = "PR" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + '.sh'
		files_to_delete.append( sh_file )
		out_file = prefix + ID + '.out'
		files_to_delete.append( out_file )
		err_file = prefix + ID + '.report'
		
		outputfile = bam_file + "recalibrated.bam"
		resulting_bam_files.append( outputfile )
		
		cmd = "java -Xmx32g -jar " + GATK + " -T PrintReads -R " + input_ref_file + " -I " + bam_file + " -BQSR " + recal_files[ idx ] + " -U ALLOW_N_CIGAR_READS " + " -o " + outputfile
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + '"',
																"| qsub -cwd",
																"-N",
																ID,
																"-l vf=35G",
																"-l arch=lx-amd64",
																"-l idle=1",
																"-P fair_share",
																"-o",
																out_file,
																"-e",
																err_file
															] ) + '\n'
							)
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "PR" + batch_ID + "_\d{4}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 10 )
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "PR" + batch_ID + "_\d{4}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 10 )
	
	# --- delete temp files --- #
	for each in files_to_delete:
		try:
			os.remove( each )
		except:
			"file cannot be deleted, because it is missing"
	
	return resulting_bam_files


def run_realignment_on_cluster( prefix, bam_files, input_ref_file, target_lists, GATK, para_jobs ):
	"""! @brief write shell script and submit realignment job to cluster """
	
	IDs_to_check = []
	files_to_delete = []
	resulting_bam_files = []
	
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, bam_file in enumerate( bam_files ):
		output_file = bam_file + "_InDels_realigned.bam"
		resulting_bam_files.append( output_file )
		
		ID = "GR" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + '.sh'
		files_to_delete.append( sh_file )
		out_file = prefix + ID + '.out'
		files_to_delete.append( out_file )
		err_file = prefix + ID + '.report'
		
		cmd = "java -Xmx32g -jar " + GATK + " -T IndelRealigner -R " + input_ref_file + " -I " + bam_file + " -targetIntervals " + target_lists[ idx ] + " -U ALLOW_N_CIGAR_READS " + " -o " + output_file
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + '"',
																"| qsub -cwd",
																"-N",
																ID,
																"-l vf=35G",
																"-l arch=lx-amd64",
																"-l idle=1",
																"-P fair_share",
																"-o",
																out_file,
																"-e",
																err_file
															] ) + '\n'
							)
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "GR" + batch_ID + "_\d{4}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 10 )
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "GR" + batch_ID + "_\d{4}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 10 )
	
	# --- delete temp files --- #
	for each in files_to_delete:
		try:
			os.remove( each )
		except:
			"file cannot be deleted, because it is missing"
	
	return resulting_bam_files


def run_haplotype_caller_on_cluster( prefix, bam_files, input_ref_file, GATK, para_jobs ):
	"""! @brief run variant detection itself """
	
	IDs_to_check = []
	files_to_delete = []
	resulting_vcf_files = []
	
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, bam_file in enumerate( bam_files ):
		output_file = bam_file + "_raw_variants.g.vcf"
		resulting_vcf_files.append( output_file )
		
		ID = "HC" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + '.sh'
		files_to_delete.append( sh_file )
		out_file = prefix + ID + '.out'
		files_to_delete.append( out_file )
		err_file = prefix + ID + '.report'
		
		cmd = "java -Xmx32g -jar " + GATK + " -T HaplotypeCaller -R " + input_ref_file + " -I " + bam_file + " --genotyping_mode DISCOVERY --emitRefConfidence GVCF -stand_call_conf 10 -nct 4 -U ALLOW_N_CIGAR_READS " + " -o " + output_file
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + '"',
																"| qsub -cwd",
																"-N",
																ID,
																"-l vf=10G",
																"-l arch=lx-amd64",
																"-l idle=1",
																'-pe multislot 4',
																"-P fair_share",
																"-o",
																out_file,
																"-e",
																err_file
															] ) + '\n'
							)
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "HC" + batch_ID + "_\d{4}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 1 )
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "HC" + batch_ID + "_\d{4}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 1 )
	
	# --- delete temp files --- #
	for each in files_to_delete:
		try:
			os.remove( each )
		except:
			"file cannot be deleted, because it is missing"
	
	return resulting_vcf_files


def combine_vcf_files( prefix, vcf_files ):
	"""! @brief combines VCF files of all chromosomes for final filtering """
	
	combined_raw_vcf_file = prefix + "combined_raw_vcf_file.g.vcf"
	
	# --- collect header lines of all VCF files --- #
	fileformat_line = ""
	filter_lines = []
	format_lines = []
	GATK_command_lines = []
	info_lines = []
	contig_lines = []
	header = ""
	
	with open( vcf_files[0], "r" ) as f:
		line = f.readline()
		while line:
			if '##fileformat=' in line:
				fileformat_line = line
			elif '##FILTER=' in line:
				filter_lines.append( line )
			elif '##FORMAT=' in line:
				format_lines.append( line )
			elif '##GATKCommandLine.HaplotypeCaller=' in line:
				GATK_command_lines.append( line )
			elif '##INFO=' in line:
				info_lines.append( line )
			elif '##contig=' in line:
				contig_lines.append( line )
			elif '#CHROM' in line:
				header = '\t'.join( line.split('\t')[:-1] ) + '\tDATA_SET\n'
			
			if line[0] != '#':
				break
			line = f.readline()
	
	for vcf_file in vcf_files[1:]:
		with open( vcf_file, "r" ) as f:
			line = f.readline()
			while line:
				if '##GATKCommandLine.HaplotypeCaller=' in line:
					GATK_command_lines.append( line )
				
				if line[0] != '#':
					break
				line = f.readline()
	
	# --- write header lines into output VCF file --- #
	with open( combined_raw_vcf_file, "w" ) as out:
		out.write( fileformat_line )
		out.write( ''.join( filter_lines ) )
		out.write( ''.join( format_lines ) )
		out.write( ''.join( GATK_command_lines ) )
		out.write( ''.join( info_lines ) )
		out.write( ''.join( contig_lines ) )
		out.write( header )
		
		# --- write variant content into combined VCF file --- #
		number_of_raw_variants = 0
		for vcf_file in vcf_files:
			with open( vcf_file, "r" ) as f:
				line = f.readline()
				while line:
					if line[0] != '#':
						out.write( line )
						number_of_raw_variants += 1
					line = f.readline()
	print "number of raw variants: " + str( number_of_raw_variants )
	return combined_raw_vcf_file


def main( arguments):
	"""! @brief calls all functions involved in the GATK variant detection """
	
	
	### NEEDS TO BE SET ###
	input_bam_file = arguments[ arguments.index( '--input_bam_file' )+1 ]
	input_ref_file = arguments[ arguments.index( '--ref_file' )+1 ]
	known_sites_vcf = arguments[ arguments.index( '--gold_vcf' )+1 ]
	
	
	prefix = arguments[ arguments.index( '--directory' )+1 ]
	if prefix[-1] != '/':
		prefix += "/"
	
	if '--para_jobs' in arguments:
		para_jobs = int( arguments[ arguments.index( '--para_jobs' )+1 ] )
	else:
		para_jobs=100
	
	
	# --- setting tool path names --- #
	piccard_tools = arguments[ arguments.index( '--piccard' )+1 ]		#v2.5.0
	samtools = arguments[ arguments.index( '--samtools' )+1 ]
	GATK = arguments[ arguments.index( '--gatk' )+1 ]	#v3.8
	
	variant_call_preparator = arguments[ arguments.index( '--varcallprepscript' )+1 ] #variant_call_preparation.py
	
	# --- construct directory for all files --- #
	if os.path.exists( prefix ):
		sys.exit( "given directory already exists - nothing should be overwritten" )
	else:
		os.makedirs( prefix )
	
	# -- transfering and indexing reference file --- #
	print "transferring reference file ..."
	cmd = "cp " + input_ref_file + " " + prefix
	os.popen( cmd )
	cmd = "cp " + known_sites_vcf + " " + prefix
	os.popen( cmd )
	known_sites_vcf = prefix + known_sites_vcf.split('/')[-1]
	input_ref_file = prefix + input_ref_file.split('/')[-1]
	index_ref_file( input_ref_file, piccard_tools, samtools )
	
	
	# --- sorting and indexing BAM file prior to splitting --- #
	if not '--bam_is_sorted' in arguments:
		print "transferring BAM file ... "
		cmd = "cp " + input_bam_file + " " + prefix
		os.popen( cmd )
		print "sorting BAM file ... "
		sorted_bam_file = prefix + input_bam_file.split('/')[-1] + "_sorted.bam"
		sort_bam_file( prefix, input_bam_file, samtools, piccard_tools )
	else:
		sorted_bam_file = input_bam_file
	index_bam_file( sorted_bam_file, piccard_tools )
	
	
	# --- preparation of splitting and processing separate BAM files on cluster --- #
	headers = load_sequences( input_ref_file ).keys()
	raw_bam_files = []
	duplicated_marked_bam_files = []
	target_lists = []
	
	for header in headers:	#everything is done for a single chromosome
		bam_file = prefix + header + ".bam"
		
		raw_bam_files.append( bam_file )
		duplicated_marked_bam_files.append( bam_file + "_RG_added.bam_duplicates_marked.bam" )
		target_lists.append( bam_file + "_RG_added.bam_duplicates_marked.bam_realignment_targets.list" )
	
	
	# --- prepare files --- #
	prepare_files_for_variant_calling_on_cluster( prefix, sorted_bam_file, headers, raw_bam_files, input_ref_file, samtools, piccard_tools, GATK, para_jobs, variant_call_preparator )
	
	
	# --- base recalibration --- #
	first_recal_files = recalibration_round( GATK, input_ref_file, duplicated_marked_bam_files, known_sites_vcf, prefix, para_jobs )
	second_recal_files = recalibration_round( GATK, input_ref_file, duplicated_marked_bam_files, known_sites_vcf, prefix, para_jobs )
	recalibrated_bam_files = base_recalibration( GATK, input_ref_file, duplicated_marked_bam_files, second_recal_files, prefix, para_jobs )
	
	
	# --- do realignment of reads around InDels --- #
	resulting_bam_files = run_realignment_on_cluster( prefix, recalibrated_bam_files, input_ref_file, target_lists, GATK, para_jobs )
	
	
	# --- identify varaints via HaplotypeCaller --- #
	raw_vcf_files = run_haplotype_caller_on_cluster( prefix, resulting_bam_files, input_ref_file, GATK, para_jobs )
	
	
	# --- combine all VCF files for final filtering --- #
	combined_raw_vcf_file = combine_vcf_files( prefix, raw_vcf_files )
	
	
	# --- sort variants in combined vcf file --- #
	sorted_combined_raw_vcf_file = combined_raw_vcf_file + "_sorted.g.vcf"
	
	dict_file_for_sorting = ".".join( input_ref_file.split( '.' )[:-1] ) + '.dict'
	cmd = "java -Xmx32g -jar " + piccard_tools + " SortVcf I=" + combined_raw_vcf_file + " O=" + sorted_combined_raw_vcf_file + " SEQUENCE_DICTIONARY=" + dict_file_for_sorting
	os.popen( cmd )
	

if __name__ == '__main__':
	
	
	if '--input_bam_file' in sys.argv and '--ref_file' in sys.argv and '--gold_vcf' in sys.argv and '--directory' in sys.argv and '--piccard' in sys.argv and '--samtools' in sys.argv and '--gatk' in sys.argv and '--varcallprepscript' in sys.argv:
		main( sys.argv)
	else:
		sys.exit( __usage__ )
	
	print "all done!"
