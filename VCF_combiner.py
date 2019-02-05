### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###

__usage__ = """
					python VCF_combiner.py
					--in <INPUT_DIRECTORY>
					--out <OUTPUT_VCF>
					"""

import re, glob

# --- end of imports --- #

def load_vcf_content( vcf ):
	"""! @brief load content of VCF file """
	
	variants = {}
	with open( vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if not ',' in parts[4]:
					if parts[6] == "PASS":
						variants.update( { parts[0] + "&" + parts[1].zfill(8) + "&.&" + parts[3] + "&" + parts[4]: parts[-1] } )
			line = f.readline()
	return variants


def main( arguments ):
	"""! @brief run all parts """
	
	inputd_dir = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	cutoff = 5
	
	filenames = glob.glob( inputd_dir + "*.vcf" )
	
		
	# --- loading data --- #
	all_keys = {}
	data = {}
	for filename in filenames:
		ID = filename.split('/')[-1].split('.')[0]
		variants = load_vcf_content( filename )
		for key in variants.keys():
			try:
				all_keys[ key ]
			except KeyError:
				all_keys.update( { key: None } )
		data.update( { ID: variants } )
	
	# --- generating output file --- #
	sorted_IDs = sorted( data.keys() )
	with open( output_file, "w" ) as out:
		out.write( "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\t" + "\t".join( sorted_IDs ) + '\n' )
		for key in sorted( all_keys.keys() ):
			new_line = key.split('&') + [ ".", "PASS", ".", "GT:AD:DP:GQ:PL" ]
			new_line[1] = str( int( new_line[1]  ) )
			counter = 0
			for ID in sorted_IDs:
				try:
					new_line.append( data[ ID ][ key ] )
					counter += 1
				except KeyError:
					new_line.append( "./.:0,0:0:.:0,0,0" )
			if counter > cutoff:
				out.write( "\t".join( new_line ) + "\n" )

if '--in' in sys.argv and '--out' in sys.argv:
	main()
else:
	sys.exit( __usage__ )
