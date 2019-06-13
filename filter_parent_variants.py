### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.15 ###

__usage__ = """
					python filter_parent_variants.py
					--vcf <FULL_PATH_TO_VCF_FILE>
					--cov <FULL_PATH_TO_COVERAGE_FILE>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					--min_cov_cov <MINIMAL_COV_OF_POSITION_IN_COV>
					--max_cov_cov <MAXIMAL_COV_OF_POSITION_IN_COV>
					--min_cov_vcf <MINIMAL_COV_OF_VARIANT_IN_VCF>
					--max_cov_vcf <MAXIMAL_COV_OF_VARIANT_IN_VCF>
					--black_vcf <FULL_PATH_TO_BLACKLIST_VCF>
					"""


import sys, os

# --- end of imports --- #


def load_cov( cov_file ):
	"""! @brief load all information from coverage file """
	
	cov = {}
	with open( cov_file, "r" ) as f:
		line = f.readline()
		header = line.split('\t')[0]
		tmp = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != header:
				cov.update( { header: tmp } )
				header = parts[0]
				tmp = []
			tmp.append( float( parts[-1] ) )
			line = f.readline()
		cov.update( { header: tmp } )
	return cov


def load_variants( black_vcf ):
	"""! @brief load all variant positions from given VCF """
	
	variants = {}
	
	with open( black_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				variants.update( { parts[0] + '_%_' + ( parts[1].zfill( 9 ) ): None } )
			line = f.readline()
	return variants


def main( arguments ):
	"""! @brief run everything """
	
	input_vcf = arguments[ arguments.index( '--vcf' ) + 1 ]
	cov_file = arguments[ arguments.index( '--cov' ) + 1 ]
	output_vcf = arguments[ arguments.index( '--out' ) + 1 ]
	
	if '--min_cov_cov' in arguments:
		min_cov_cutoff = int( arguments[ arguments.index( '--min_cov_cov' ) + 1 ] )
	else:
		min_cov_cutoff = 10
	
	if '--max_cov_cov' in arguments:
		max_cov_cutoff = int( arguments[ arguments.index( '--max_cov_cov' ) + 1 ] )
	else:
		max_cov_cutoff = 100
	
		
	if '--min_cov_vcf' in arguments:
		vcf_min_cov = int( arguments[ arguments.index( '--min_cov_vcf' ) + 1 ] )
	else:
		vcf_min_cov = 10
	
	if '--max_cov_vcf' in arguments:
		vcf_max_cov = int( arguments[ arguments.index( '--max_cov_vcf' ) + 1 ] )
	else:
		vcf_max_cov = 100
	
	if '--black_vcf' in arguments:
		black_vcf = arguments[ arguments.index( '--black_vcf' ) + 1 ]
		black = load_variants( black_vcf )
	else:
		black = {}
	
	cov = load_cov( cov_file )
	
	with open( output_vcf, "w" ) as out:
		with open( input_vcf, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] == "#":
					out.write( line )
				else:
					parts = line.strip().split('\t')
					if parts[-1][:3] == "1/1":
						dp = int( parts[-1].split(':')[2] )
						if vcf_min_cov < dp < vcf_max_cov:
							if min_cov_cutoff < cov[ parts[0] ][ int( parts[1] ) ] < max_cov_cutoff:
								try:
									black[ parts[0] + '_%_' + ( parts[1].zfill( 9 ) ) ]
								except KeyError:
									out.write( line )
				line = f.readline()


if '--vcf' in sys.argv and '--cov' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
