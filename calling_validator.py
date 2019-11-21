### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.35 ###

__usage__ = """
					python calling_validator.py\n
					--vcf <VCF_FILE>
					--fasta <FASTA_FILE>
					--out <OUTPUT_DIRECTORY>
					
					feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
					"""

import sys, os

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


def validate_vcf( vcf_file, seq_lengths, result_file, window_size ):
	"""! @brief check coverage file """
	
	valid = True
	
	vcf_chromosomes = {}
	
	with open( result_file, "w" ) as out:
		with open( vcf_file, "r" ) as f:
			line = f.readline()
			chromosome = False
			variant_pos = []
			while line:
				if line[0] != '#':
					if not chromosome:
						chromosome = line.split('\t')[0]
						vcf_chromosomes.update( { line.split('\t')[0]: None } )
					parts = line.strip().split('\t')
					if parts[0] != chromosome:
						chunks = {}
						for pos in variant_pos:
							val = pos / window_size
							try:
								chunks[ val ] += 1
							except KeyError:
								chunks.update( { val: 0 } )
						try:
							for i in range( max( chunks.keys() ) ):
								try:
									if chunks[ i ] > 0:
										out.write( str( chunks[ i ]  ) + "\t" + str( window_size ) + '\n' )
									else:
										out.write( "ERROR: no variants - " + chromosome + " - block idx: " + str( i ) + "\n" )
										valid = False
								except KeyError:
									out.write( "ERROR: no variants - " + chromosome + " - block idx: " + str( i ) + "\n" )
									valid = False
						except ValueError:
							out.write( "ERROR: no variants - " + chromosome + " - block idx: ?\n" )
							valid = False
						variant_pos = []
						chromosome = parts[0]
						vcf_chromosomes.update( { parts[0]: None } )
					variant_pos.append( int( parts[1] ) )
				line = f.readline()
			chunks = {}
			for pos in variant_pos:
				val = pos / window_size
				try:
					chunks[ val ] += 1
				except KeyError:
					chunks.update( { val: 0 } )
			for i in range( max( chunks.keys() ) ):
				try:
					if chunks[ i ] > 0:
						out.write( str( chunks[ i ]  ) + "\t" + str( window_size ) + '\n' )
					else:
						out.write( "ERROR: coverage is zero - " + chromosome + " - block idx: " + str( i ) + "\n" )
						valid = False
				except KeyError:
					out.write( "ERROR: coverage is zero - " + chromosome + " - block idx: " + str( i ) + "\n" )
					valid = False
			
			# --- check if all sequences are present --- #
			for key in seq_lengths.keys():
				try:
					vcf_chromosomes[ key ]
				except:
					out.write( "ERROR: missing chromosome - " + key + "\n" )
					valid = False
			
			out.write( "FINAL STATUS: valid? >> " + str( valid ) + '\n'  )
	return valid


def main( arguments ):
	"""! @brief run everything """
	
	vcf_file = arguments[ arguments.index( '--vcf' )+1 ]
	fasta_file = arguments[ arguments.index( '--fasta' )+1 ]
	prefix = arguments[ arguments.index( '--out' )+1 ]
	
	if not prefix[-1] == "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	window_size = 500000
		
	# --- check completeness of file --- #
	seq_lengths = load_seq_lengths( fasta_file )
	
	result_file = prefix + vcf_file.split('/')[-1].lower().replace( ".vcf", ".results" )
	status = validate_vcf( vcf_file, seq_lengths, result_file, window_size )
	if not status:
		print "ERROR detected in " + vcf_file
	else:
		print "OK!"


if __name__ == '__main__':
	
	if '--vcf' in sys.argv and '--out' in sys.argv and '--fasta' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
	print "all done!"
