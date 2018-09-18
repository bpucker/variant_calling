import os, sys

# --- end of imports --- #

def index_bam_file( sorted_bam_file, piccard_tools, java ):
	"""! @brief index given BAM file """
	
	index_file = sorted_bam_file + ".bai"
	
	cmd = [ java + " -Xmx8g -jar ",
			piccard_tools,
			" BuildBamIndex I=",
			sorted_bam_file,
			" O=",
			index_file
		 ]
	cmd = "".join( cmd )
	os.popen( cmd )
	return index_file


def add_read_groups( input_bam_file, output_bam_file, piccard_tools, java ):
	
	cmd = "".join( [ java + " -Xmx8g -jar ",
					piccard_tools,
					" AddOrReplaceReadGroups ",
					" I= ",
					input_bam_file,
					" O=",
					output_bam_file,
					" RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"
					])
	os.popen( cmd )


def mark_duplicates( bam_file, piccard_tools, java ):
	"""! @brief mark PCR duplicated reads """
	
	duplicates_marked_bam_file = bam_file + "_duplicates_marked.bam"
	metrics_file = bam_file + "_metrics.txt"
	
	cmd = [ java + " -Xmx8G -jar ",
			piccard_tools,
			" MarkDuplicates INPUT=",
			bam_file,
			" OUTPUT=", 
			duplicates_marked_bam_file,
			" METRICS_FILE=",
			metrics_file,
			"OPTICAL_DUPLICATE_PIXEL_DISTANCES=2500 CREATE_INDEX=true TMP_DIR=/tmp"
			]
	cmd = "".join( cmd )
	os.popen( cmd )
	return duplicates_marked_bam_file


def construct_target_list_for_realigner( duplicates_marked_bam_file, ref_fasta_file, GATK, header, java ):
	"""! @brief sorted+indexed BAM file and indexed reference file are used for the construction of a target list for the realignment around indels """
	
	target_list_file = duplicates_marked_bam_file + "_realignment_targets.list"
	
	cmd = [ java + " -Xmx8g -jar ",
			GATK,
			" -T RealignerTargetCreator -R ",	#--fix_misencoded_quality_scores --allow_potentially_misencoded_quality_scores 
			ref_fasta_file,
			" -I ",
			duplicates_marked_bam_file,
			" -L " + header,
			" -U ALLOW_N_CIGAR_READS ",
			" -o ",
			target_list_file
			]
	cmd = "".join( cmd )
	os.popen( cmd )
	return target_list_file


def main( arguments ):
		
		# --- required inputs --- #
		original_bam_file = arguments[ arguments.index( '--original_bam_file' )+1 ]
		ref_file = arguments[ arguments.index( '--ref_file' )+1 ]
		bam_file = arguments[ arguments.index( '--bam_file' )+1 ]
		header = arguments[ arguments.index( '--header' )+1 ]
		samtools = arguments[ arguments.index( '--samtools' )+1 ]
		GATK = arguments[ arguments.index( '--GATK' )+1 ]
		piccard_tools = arguments[ arguments.index( '--piccard_tools' )+1 ]
		
		java = "java"	#should be replaced by full path to java version
		
		# --- construction of bam file --- #
		cmd = samtools + " view " + original_bam_file + " " + header + " -b > " + bam_file
		os.popen( cmd )
		
		
		# --- indexing BAM file --- #
		index_bam_file( bam_file, piccard_tools, java )
		
		
		# --- adding read group --- #
		rg_bam_file = bam_file + "_RG_added.bam"
		add_read_groups( bam_file, rg_bam_file, piccard_tools, java )
		
		
		# --- indexing BAM file --- #
		index_bam_file( rg_bam_file, piccard_tools, java )
		
		
		## --- mark PCR duplicates --- #
		mark_duplicates( rg_bam_file, piccard_tools, java )
		duplicates_marked_bam_file = rg_bam_file + "_duplicates_marked.bam"
		
		
		## --- indexing BAM file --- #
		index_bam_file( duplicates_marked_bam_file, piccard_tools, java )
		
		
		# --- construct realignmet target list --- #
		construct_target_list_for_realigner( duplicates_marked_bam_file, ref_file, GATK, header, java )


if __name__ == '__main__':
	
	main( sys.argv )
