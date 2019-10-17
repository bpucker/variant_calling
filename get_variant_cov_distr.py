### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python get_variant_cov_distr.py
					--vcf <VCF_FILE(INPUT)>
					--fig <FIGURE_FILE(OUTPUT)>
					"""

import sys
import matplotlib.pyplot as plt

# --- end of imports --- #

def main( arguments ):
	
	input_vcf = arguments[ arguments.index('--vcf')+1 ]
	fig_file = arguments[ arguments.index('--fig')+1 ]

	variant_covs = []

	with open(  input_vcf, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if parts[6] == "PASS":
				try:
					variant_covs.append( int( parts[-1].split(':')[2] ) )
				except IndexError:
					pass
			line = f.readline()


	fig, ax = plt.subplots()
	ax.hist( variant_covs, bins=1000, color="lime" )
	ax.set_xlim( 0, int( 5* ( sum( variant_covs ) / float( len( variant_covs ) ) ) ) )
	ax.set_xlabel( "coverage" )
	ax.set_ylabel( "number of variant positions" )

	fig.savefig( fig_file, dpi=300 )


if '--vcf' in sys.argv and '--fig' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
