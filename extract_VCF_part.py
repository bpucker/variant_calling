### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1  ###


__usage__ = """
					python extract_VCF_part.py
					--in <FULL_PATH_TO_INPUT_FILE>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					--chr <STR, name of chromosome of interest>
					--start <INT, start of region of interest>
					--end <INT, end of region of interest>
					"""

import sys

# --- end of imports --- #

def main( parameters ):
	"""! @brief run everything """

	input_vcf = parameters[ parameters.index( '--in' ) + 1 ]
	output_vcf = parameters[ parameters.index( '--out' ) + 1 ]


	chromosome = parameters[ parameters.index( '--chr' ) + 1 ]

	start = int( parameters[ parameters.index( '--start' ) + 1 ] )
	end = int( parameters[ parameters.index( '--end' ) + 1 ] )

	counter = 0
	with open( output_vcf, "w" ) as out:
		with open( input_vcf, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] == "#":
					out.write( line )
				else:
					parts = line.strip().split('\t')
					if parts[0] == chromosome:
						if start <= int( parts[1] ) <= end:
							out.write( line )
							counter += 1
				line = f.readline()
	print "number of extracted lines: " + str( counter )


if '--in' in sys.argv and '--out' in sys.argv and '--chr' in sys.argv and '--start' in sys.argv and '--end' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
