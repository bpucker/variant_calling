### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.1 ###

__usage__ = """
					python3 correct_VCF_format.py
					--in <INPUT_VCF_FILE>
					--out <OUTPUT_VCF_FILE>
					"""

import os, sys

# --- end of imports --- #

def main( arguments ):
	"""! @brief run everything """

	input_vcf = arguments[ arguments.index('--in')+1 ]
	output_vcf = arguments[ arguments.index('--out')+1 ]

	with open( output_vcf, "w" ) as out:
		with open( input_vcf, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] == "#":
					if line[1] == "#":
						out.write( line )
					else:
						out.write( "\t".join( [ "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT" ] ) + "\n" )
				else:
					parts = line.strip().split('\t')
					parts.append( "." )
					out.write( "\t".join( parts )+"\n" )
				line = f.readline()


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
