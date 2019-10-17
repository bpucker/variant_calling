### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
				python filter_VCF_by_gold_standard.py
				--vcf <VCF (INPUT)>
				--gold <GOLD_VCF (INPUT)>
				--out <VCF (OUTPUT)>
				"""

import sys
import matplotlib.pyplot as plt

# --- end of imports --- #

def load_all_variants( vcf ):
	"""! @brief load all variant positions """
	
	variants = {}
	with open( vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split( '\t' )
				variants.update( { parts[0] + "_%_" + parts[1]: None } )
			line = f.readline()
	return variants


def main( arguments ):
	"""! @brief runs everything """
	
	input_vcf = arguments[ arguments.index('--vcf')+1 ]
	gold_vcf = arguments[ arguments.index('--gold')+1 ]
	output_vcf = arguments[ arguments.index('--out')+1 ]
		
	gold = load_all_variants( gold_vcf )
	
	counter = 0
	with open( output_vcf, "w" ) as out:
		with open(  input_vcf, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					try:
						gold[ parts[0] + "_%_" + parts[1] ]
						out.write( line )
						counter += 1
					except KeyError:
						pass
				else:
					out.write( line )
				line = f.readline()
	print "number of remaining variants: " + str( counter )

if '--vcf' in sys.argv and '--gold' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
