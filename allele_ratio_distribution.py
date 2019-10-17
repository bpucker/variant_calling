### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
				python allele_ratio_distribution.py
				--vcf <VCF>
				--fig <FIGURE_FILE>
				"""

import sys



import matplotlib.pyplot as plt

# --- end of imports --- #

def main( arguments ):
	"""! @brief runs everything """
	
	input_vcf = arguments[ arguments.index('--vcf')+1 ]
	fig_file = arguments[ arguments.index('--fig')+1 ]

	ratios = []
	
	with open(  input_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				try:
					x, y = map( float, parts[-1].split(':')[1].split(',') )
					ratios.append( y / (x+y) )
				except IndexError:
					print line
			line = f.readline()

	print min( ratios )
	print max( ratios )
	
	fig, ax = plt.subplots()
	ax.hist( ratios, bins=100, color="lime" )
	ax.set_xlim( 0, 1 )
	ax.set_xlabel( "frequence of alternative allele" )
	ax.set_ylabel( "number of variants" )

	fig.savefig( fig_file, dpi=300 )


if '--vcf' in sys.argv and '--fig' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
