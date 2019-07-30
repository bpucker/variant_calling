### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python mapping_validator.py\n
					--bam <BAM_FILE>
					--fasta <FASTA_FILE>
					--out <OUTPUT_DIRECTORY>
					
					--bam_is_sorted <PREVENTS_EXTRA_SORTING_OF_BAM_FILE>
					
					feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
					"""

__cite__ = """ based on Pucker & Brockington, 2018: https://doi.org/10.1186/s12864-018-5360-z """


import os, sys

# --- end of imports --- #


def load_seq_lengths( fasta_file ):
	"""! @brief load sequence lengths from given FASTA file """
	
	seq_lengths = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip().split(" ")[0]
		seq = 0
		line = f.readline()
		while line:
			if line[0] == '>':
					seq_lengths.update( { header: seq } )
					header = line.strip()[1:].split(" ")[0]
					seq = 0
			else:
				seq += len( line.strip() )
			line = f.readline()
		seq_lengths.update( { header: seq } )
	return seq_lengths


def validate_cov( cov_file, seq_lengths, result_file, window_size ):
	"""! @brief check coverage file """
	
	valid = True
	
	with open( result_file, "w" ) as out:
		with open( cov_file, "r" ) as f:
			line = f.readline()
			chromosome = line.split('\t')[0]
			covs = []
			while line:
				parts = line.strip().split('\t')
				if parts[0] != chromosome:
					if len( covs ) != seq_lengths[ chromosome ]:
						out.write( "ERROR: coverage values not covering chromosome length - " + chromosome + "\n" )
						valid = False
					chunks = [ covs [ i:i + window_size ] for i in xrange( 0, len( covs ), window_size ) ]
					for idx, chunk in enumerate( chunks ):
						if len( chunk ) >= 1:
							try:
								mean = sum( chunk ) / len( chunk )
								out.write( str( mean ) + "\t" + str( len( chunk ) ) + '\n' )
								if mean == 0:
									out.write( "ERROR: coverage is zero - " + chromosome + " - block idx: " + str( idx ) + "\n" )
									valid = False
							except ZeroDivisionError:
								out.write( "ZeroDivisionError: " + chromosome + " - " + str( idx ) + "\n" )
					covs = []
					chromosome = parts[0]
				covs.append( float( parts[2] ) )
				line = f.readline()
			if len( covs ) != seq_lengths[ chromosome ]:
				out.write( "ERROR: coverage values not covering chromosome length - " + chromosome + "\n" )
				valid = False
			chunks = [ covs [ i:i + window_size ] for i in xrange( 0, len( covs ), window_size ) ]
			for idx, chunk in enumerate( chunks ):
				if len( chunk ) >= 1:
					try:
						mean = sum( chunk ) / len( chunk )
						out.write( str( mean ) + "\t" + str( len( chunk ) ) + '\n' )
						if mean == 0:
							out.write( "ERROR: coverage is zero - " + chromosome + " - block idx: " + str( idx ) + "\n" )
							valid = False
					except ZeroDivisionError:
						out.write( "ZeroDivisionError: " + chromosome + " - " + str( idx ) + "\n" )
			out.write( "FINAL STATUS: valid? >> " + str( valid ) + '\n'  )
	return valid


def main( arguments ):
	"""! @brief run everything """
	
	bam_file = arguments[ arguments.index( '--bam' )+1 ]
	fasta_file = arguments[ arguments.index( '--fasta' )+1 ]
	prefix = arguments[ arguments.index( '--out' )+1 ]
	
	if not prefix[-1] == "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	samtools = "samtools"
	bedtools = "genomeCoverageBed"
	
	window_size = 500000
	
	# --- sorting of BAM file if necessary --- #
	if '--bam_is_sorted' in arguments:
		sorted_bam_file = bam_file
	else:
		#print "sorting BAM file ..."
		sorted_bam_file = output_file + "_sorted.bam"
		cmd = samtools + " sort -m 5000000000 --threads 8 " + bam_file + " > " + sorted_bam_file
		os.popen( cmd )
	
	
	# --- calculate read coverage depth per position --- #
	#print "calculating coverage per position ...."
	cov_file = prefix + sorted_bam_file.split('/')[-1].split('.bam')[0] + ".cov"
	cmd = bedtools + " -d -split -ibam " + sorted_bam_file + " > " + cov_file
	os.popen( cmd )
	
	# --- check completeness of file --- #
	seq_lengths = load_seq_lengths( fasta_file )
	
	result_file = cov_file.replace( ".cov", ".results" )
	status = validate_cov( cov_file, seq_lengths, result_file, window_size )
	if not status:
		print "ERROR detected in " + sorted_bam_file
	else:
		print "OK!"


if __name__ == '__main__':
	
	if '--bam' in sys.argv and '--out' in sys.argv and '--fasta' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
	print "all done!"
