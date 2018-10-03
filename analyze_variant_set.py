### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python analyze_variant_set.py
					--vcf <FULL_PATH_TO_VCF_FILE (INPUT)>
					--fig  <FULL_PATH_TO_FIGURE_FILE (OUTPUT)>
					--report <FULL_PATH_TO_REPORT_FILE (OUTPUT)>
					"""

import matplotlib.pyplot as plt
import numpy as np
import sys

# --- end of imports --- #


def load_variant_infos( input_vcf ):
	"""! @brief load variant infos """
	
	snps = 0
	mnps = 0
	insertion = 0
	deletion = 0
	
	variants = {}
	
	with open( input_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if len( parts[3] ) == 1:
					if len( parts[4] ) == len( parts[3] ):
						snps += 1
					elif len( parts[3] ) > len( parts[4] ):
						deletion += 1
					elif len( parts[4] ) > len( parts[3] ):
						insertion += 1
					else:
						print "ERROR: alternative allele is missing/gap?!"
						print line
				elif len( parts[3] ) == len( parts[4] ):
					mnps += 1
				elif len( parts[3] ) > len( parts[4] ):
					deletion += 1
				elif len( parts[3] ) < len( parts[4] ):
					insertion += 1
				else:
					print "ERROR: no classification possible!"
				try:
					variants[ parts[0] ].append( int( parts[1] ) )
				except KeyError:
					variants.update( { parts[0]: [ int( parts[1] ) ] } )
			line =f.readline()
	return variants, snps, mnps, insertion, deletion


def get_variant_dist_distribution( variants, fig_file ):
	"""! @brief get variant distribution """
	
	distances = []
	for values in variants.values():
		for idx, pos in enumerate( sorted( values ) ):
			if idx == len( values )-1:
				pass
			else:
				distances.append( values[idx+1]-pos )
	
	avg_mean = np.mean( distances )
	avg_median = np.median( distances )
	
	fig, ax = plt.subplots()
	
	minimal = 0.0
	maximal = 20000
	ax.set_yscale('log')
	ax.hist( distances, bins=1000, color="green", range=( minimal, maximal ), log=True )
		
	ax.set_xlabel( "distance of variants [bp]" )
	ax.set_ylabel( "number of variants" )
	
	ax.set_xlim( minimal, maximal )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	
	fig.savefig( fig_file, dpi=900 )
	
	return avg_mean, avg_median


def main( arguments ):
	"""! @brief run all parts """
	
	input_vcf = arguments[ arguments.index( '--vcf' )+1 ]
	fig_file = arguments[ arguments.index( '--fig' )+1 ]
	report_file = arguments[ arguments.index( '--report' )+1 ]
	
	variants, snps, mnps, insertion, deletion = load_variant_infos( input_vcf )
	
	dist_mean, dist_median = get_variant_dist_distribution( variants, fig_file )
	total = snps+mnps+insertion+deletion
	with open( report_file, "w" ) as out:
		out.write( "SNPs: " + str( snps ) + '\t' + str( 100.0*snps / total )+"%" + '\n' )
		out.write( "MNPs: " + str( mnps ) + '\t' + str( 100.0*mnps / total ) +"%" + '\n' )
		out.write( "insertions: " + str( insertion ) + '\t' + str( 100.0*insertion / total )+"%"  + '\n' )
		out.write( "deletions: " + str( deletion ) + '\t' + str( 100.0*deletion / total ) +"%" + '\n' )
		out.write( "total: " + str( total ) + '\n' )
		out.write( "average variant distance (mean): " + str( dist_mean ) + '\n' )
		out.write( "average variant distance (median): " + str( dist_median ) + '\n' )


if __name__ == '__main__':
	if '--vcf' in sys.argv and '--fig' in sys.argv and '--report' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
