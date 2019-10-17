### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python vcf_cleaner.py
					--in <INPUT_VCF>
					--out <OUTPUT_VCF>
					
					Filter criteria:
					1) PASS
					2) max InDel / MNP length
					3) min read cutoff / max read cutoff (coverage)
					4) biallelic
					"""

import sys

# --- end of imports --- #


def main( arguments ):
	"""! @brief run everything """
	
	input_vcf = arguments[ arguments.index('--in')+1 ]
	output_vcf = arguments[ arguments.index('--out')+1 ]
	
	max_length = 1
	min_reads = 20
	max_reads = 500
	
	filter_counter = 0
	counter = 0
	low_cov_counter = 0
	high_cov_counter = 0
	remaining_counter = 0
	multiallelic_cunter = 0
	with open( output_vcf, "w" ) as out:
		with open( input_vcf, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] == "#":
					out.write( line )
				else:
					parts = line.strip().split('\t')
					if parts[ 6 ] == "PASS":
						if max( [ len( parts[3] ), len( parts[4] ) ] ) <= max_length:	#length cutoff
							try:
								x, y = map( float, parts[-1].split(':')[1].split(',')[:2] )
								if min_reads < x+y < max_reads:
									if not ',' in parts[4]:
										out.write( line )
										remaining_counter += 1
									else:
										multiallelic_cunter += 1
								else:
									if x+y < min_reads:
										low_cov_counter += 1
									else:
										high_cov_counter += 1
							except IndexError:	#no genotyping information available for this line
								pass
						else:
							counter += 1
					else:
						filter_counter += 1
				line = f.readline()
	
	print "number of variants removed due to filter: " + str( filter_counter )
	print "InDel too large:" + str( counter )
	print "coverage too low: " + str( low_cov_counter )
	print "coverage too high: " + str( high_cov_counter )
	print "number of multiallelic variants: " + str( multiallelic_cunter )
	print "number of remaining variants: " + str( remaining_counter )

if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
