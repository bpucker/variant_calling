### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.3 ###

__usage__ = """
	python combine_single_VCFs.py
	--in <FULL_PATH_TO_VCF_FILE_FOLDER>
	--out <FULL_PATH_TO_OUTPUT_FILE>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import glob, sys

# --- end of imports --- #


def load_variants( vcf ):
	"""! @brief load all variants from given VCF """
	
	variants = {}
	with open( vcf, "r" ) as f:
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts[3] ) == 1 and len( parts[4] ) == 1:
				if parts[ 6 ] == "PASS":
					variants.update( { parts[0]+'_%_'+(parts[1].zfill(8)): { 'ref': parts[3], 'alt': parts[4], 'info': parts[-1], 'filter': True } } )
				else:
					variants.update( { parts[0]+'_%_'+(parts[1].zfill(8)): { 'ref': parts[3], 'alt': parts[4], 'info': parts[-1], 'filter': False } } )
			line = f.readline()
	return variants


def main( arguments ):
	"""! @brief run everything """
	
	input_folder = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	# --- load variants from all VCF files into separate dictionaries --- #
	vcfs = sorted( glob.glob( input_folder + "*.vcf" ) )
	IDs = []
	
	global_variants = {}
	for vcf in vcfs:
		ID = vcf.split('/')[-1].split('.')[0]
		IDs.append( ID )
		variants = load_variants( vcf )
		global_variants.update( { ID: variants } )
	
	# --- merge variant positions of all VCF files --- #
	merged_keys = {}
	for each in global_variants.values():
		for key in each.keys():
			if each[ key ]['filter']:
				try:
					merged_keys[ key ]
				except KeyError:
					merged_keys.update( { key: None } )
	merged_keys = sorted( merged_keys.keys() )
	print "number of keys: " + str( len( merged_keys ) )
	
	# --- generate output file --- #
	with open( output_file, "w" ) as out:
		out.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join( IDs ) + '\n' )
		for key in merged_keys:
			ref = False
			alt = False
			for ID in IDs:
				try:
					if not ref:
						ref = global_variants[ ID ][ key ]['ref']
					if not alt:
						alt = global_variants[ ID ][ key ]['alt']
					else:
						if global_variants[ ID ][ key ]['alt'] != alt:
							ref = False
							alt = False
							break
				except KeyError:
					pass
			if ref and alt:
				new_line = [ key.split('_%_')[0], str( int( key.split('_%_')[-1] ) ), ".", ref, alt, ".", "PASS", ".", "GT:AD:DP:GQ:PL" ]
				#str( int() ) to remove any leading 0s
				for ID in IDs:
					try:
						new_line.append( global_variants[ ID ][ key ]['info'] )
					except KeyError:
						new_line.append( "./.:0,0:0:.:0,0,0" )
				out.write( "\t".join( new_line ) + '\n' )


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
