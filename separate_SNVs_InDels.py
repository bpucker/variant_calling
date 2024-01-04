### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.1 ###

__usage__ = """
					python3 separate_SNVs_InDels.py
					--in <INPUT_VCF_FILE>
					--snvout <OUTPUT_SNV_FILE>
					--indelout <OUTPUT_INDEL_FILE>
					"""

import os, sys

# --- end of imports --- #

def main( arguments ):
	"""! @brief run everything """
	
	input_vcf = arguments[ arguments.index('--in')+1 ]
	snv_output_vcf = arguments[ arguments.index('--snvout')+1 ]
	indel_output_vcf = arguments[ arguments.index('--indelout')+1 ]

	snv_counter = 0
	indel_counter = 0
	with open( snv_output_vcf, "w" ) as snv_out:
		with open( indel_output_vcf, "w" ) as indel_out:
			with open( input_vcf, "r" ) as f:
				line = f.readline()
				while line:
					if line[0] == "#":
						snv_out.write( line )
						indel_out.write( line )
					else:
						parts = line.strip().split('\t')
						if len( parts[3] ) == 1 and len( parts[4] ) == 1:
							snv_out.write( line )
							snv_counter += 1
						else:
							indel_out.write( line )
							indel_counter += 1
					line = f.readline()
	
	sys.stdout.write( "Total variant number: " + str( snv_counter+indel_counter ) + "\n" )
	sys.stdout.write( "SNVs: " + str( snv_counter ) + "\n" )
	sys.stdout.write( "InDels: " + str( indel_counter ) + "\n" )
	sys.stdout.flush()


if '--in' in sys.argv and '--snvout' in sys.argv and '--indelout' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
