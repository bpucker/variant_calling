### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python genome_wide_distribution_of_variants.py
					--vcf <FULL_PATH_TO_INPUT_VCF>
					--out <FULL_PATH_TO_OUTPUT_FOLDER>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""


import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys, os

# --- end of imports --- #


def load_variants_from_vcf( vcf_file ):
	"""! @brief loads the variant informaiton from a SnpEff output VCF file """
	
	snps_per_chr = [ [], [], [], [], [] ]
	indels_per_chr = [ [], [], [], [], [] ]
	
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				try:
					parts = line.strip().split('\t')
					if len( parts[3] ) == len( parts[4] ):	#SNP/MNP
						if len( parts[3] ) == 1:	#SNP
							snps_per_chr[ int( parts[0][-1] ) - 1 ].append( int( parts[1] ) )
					elif len( parts[3] ) != len( parts[4] ) :	#InDel
						indels_per_chr[ int( parts[0][-1] ) - 1 ].append( int( parts[1] ) )
				except:
					pass	#print line
			line = f.readline()
	
	return snps_per_chr, indels_per_chr


def construct_plot( snps_per_chr, indels_per_chr, chr_lengths, result_file, result_table, resolution, y_max_lim ):
	"""! @brief construct variant over Col-0 genome distribution plot """
	
	chromosome_length_units = [ 60, 40, 50, 40, 60 ]	#hard coded values!!

	f = plt.figure( figsize=(10,15) )
	gs = gridspec.GridSpec( 5, 1 )
	x_lim = [ 30, 20, 25, 20, 30 ]
	
	with open( result_table, "w" ) as out:
		for idx, chr_length in enumerate( chr_lengths ):
			chr_name = 'Chr' + str( idx + 1 )
			current_grid = gridspec.GridSpecFromSubplotSpec( 1, 60, subplot_spec = gs[ idx ] )
			
			# ---- identify positions of variants --- #
			upper_lim = resolution
			lower_lim = 0
			
			snp_data = []
			indel_data = []
			while True:
				if upper_lim >= chr_length:
					break
				else:
					snp_tmp = []
					indel_tmp = []
					for SNP in snps_per_chr[ idx ]:
						if SNP <= upper_lim and SNP > lower_lim:
							snp_tmp.append( 'X' )
					for indel in indels_per_chr[ idx ]:
						if indel <= upper_lim and indel > lower_lim:
							indel_tmp.append( 'X' )
					snp_data.append( len( snp_tmp ) )
					indel_data.append( len( indel_tmp ) )
				upper_lim += resolution
				lower_lim += resolution
			
			print "length of snp_data: " + str( len( snp_data ) )
			
			# --- improving x-axis --- #
			ax_a = plt.subplot( current_grid[ 0: chromosome_length_units[ idx ] ] )	#axes[ idx ]	#( len( snp_data ) - 1 )
			ax_a.set_xlim( 0, x_lim[ idx ] )	#len( snp_data ) / 2.0
			ax_a.set_ylim( 0, y_max_lim )
			
			# --- plotting SNP and InDel distribution --- #
			for i, snps in enumerate( snp_data ):
				if i == 0:	#add labels only once!
					ax_a.plot( ( ((i*0.5)+0.2), ((i*0.5)+0.2) ), ( 0, snps ), "-", color="black", label="SNVs" )
					ax_a.plot( ( ((i*0.5)+0.3), ((i*0.5)+0.3) ), ( 0, indel_data[ i ] ), "-", color="red", label="InDels" )
				else:
					ax_a.plot( ( ((i*0.5)+0.2), ((i*0.5)+0.2) ), ( 0, snps ), "-", color="black" )
					ax_a.plot( ( ((i*0.5)+0.3), ((i*0.5)+0.3) ), ( 0, indel_data[ i ] ), "-", color="red" )
			
			# --- writing data into output table --- #
			out.write( 'Chr' + str( idx+1 ) + "SNVs:\t" + '\t'.join( map( str, snp_data ) ) + '\n' )
			out.write( 'Chr' + str( idx+1 ) + "InDels:\t" + '\t'.join( map( str, indel_data ) ) + '\n' )
			
			
			## --- improving ticks of y-axis a --- #
			max_yticks = 3
			yloc = plt.MaxNLocator(max_yticks)
			ax_a.yaxis.set_major_locator( yloc )
			
			
			current_labels = ax_a.get_xticks()
			labels = []
			for each in current_labels:
				labels.append( each  )
			ax_a.set_xticks([])
			ax_a.set_xticks( labels )
			
			ax_a.set_title( chr_name )
			ax_a.set_xlabel( "Mbp" )
			ax_a.set_ylabel( "number of variants" )	#"counts", "SNPs (black), InDel(red)"
			ax_a.legend( prop={'size':10} )
		
	gs.update( hspace=0.75 )
	#plt.show()
	f.subplots_adjust( left=0.09, right=0.985, top=0.98, bottom=0.03 )
	f.savefig( result_file, dpi=600 )
	plt.close('all')


def main( arguments ):
	"""! @brief run everything """
	
	vcf_file = arguments[ arguments.index( '--vcf' )+1 ]
	
	output_dir = arguments[ arguments.index( '--out' )+1 ]
	
	if output_dir[-1] != '/':
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	result_file = output_dir + "genome_wide_small_variants.png"
	result_table = output_dir + "genome_wide_small_variants.txt"
	
	snps_per_chr, indels_per_chr = load_variants_from_vcf( vcf_file )
	
	print "number of SNVs: " + str( len( [ x for each in snps_per_chr for x in each ] ) )
	print "number of InDels: " + str( len( [ x for each in indels_per_chr for x in each ]) )
	
	chr_lengths = [ 30427671, 19698289, 23459830, 18585056, 26975502 ]
	
	
	# --- generate statistics --- #
	for idx, snps in enumerate( snps_per_chr ):
		print str( len( snps ) ) + "(1 SNP in " + str( int( float( chr_lengths[ idx ] ) / len( snps ) ) ) + "bp)"
		print str( len( indels_per_chr[ idx ] ) ) + "(1 InDel in " + str( int( float( chr_lengths[ idx ] ) / len( indels_per_chr[ idx ] ) ) ) + "bp)"
	
	
	construct_plot( snps_per_chr, indels_per_chr, chr_lengths, result_file, result_table, resolution = 500000, y_max_lim=10500 )


if '--vcf' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
