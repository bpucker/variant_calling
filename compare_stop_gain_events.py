### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.2 ###

#WARNING: This script was written to handle VCFs containing only SNVs (no InDels). Inclusion of InDels could lead to complex cases of stop loss events that might not be handled properly.

__usage__ = """
					python3 compare_stop_gain_events.py.py
					--snpeffvcf <SnpEff_VCF_OUTPUT_FILE>
					--navipvcf <NAVIP_VCF_OUTPUT_FILE>
					--out <OUTPUT_FOLDER>
					"""

import os, sys, re
import matplotlib.pyplot as plt
from operator import itemgetter

# --- end of imports --- #

def load_snpeff_stop_gained_cases( input_vcf_file ):
	"""! @brief load SnpEff stop_gained cases with two SNVs per codon """
	
	all_snpeff_variants = {}	#load all variants for later
	stop_gained_cases = {}
	cds_pos_of_first_snpeff_stop_per_gene = {}	#collect the CDS position of the first stop variant in each gene (can replace gene orientation information)
	with open( input_vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				all_snpeff_variants.update( { parts[0] + "_%_" + parts[1]: parts } )
				if "stop_gained" in parts[7]:
					if "," in parts[7]:
						subparts = parts[7].split(',')
						for subpart in subparts:
							if "stop_gained" in subpart:
								parts[7] = subpart
					stop_gained_cases.update( { parts[0] + "_%_" + parts[1]: parts } )
			line = f.readline()
	
	# --- filter stop gains by gene to report only first stop gain per gene --- #
	stop_gains_per_gene = {}
	for variant in list( stop_gained_cases.values() ):
		gene = variant[7].split('|')[3]
		try:
			stop_gains_per_gene[ gene ].append( {'pos': int( variant[1] ), 'chr': variant[0], 'cdspos': int( variant[7].split('|')[12].split('/')[0] ), 'content': variant } )
		except KeyError:
			stop_gains_per_gene.update( { gene: [ {'pos': int( variant[1] ), 'chr': variant[0], 'cdspos': int( variant[7].split('|')[12].split('/')[0] ), 'content': variant } ] } )
	final_stop_gained_cases = {}
	for gene in list( stop_gains_per_gene.keys() ):
		stops = stop_gains_per_gene[ gene ]
		if len( stops ) > 1:	#find the first predicted premature stop codon and ignore the following ones
			stops = sorted( stops, key=itemgetter('cdspos') )
		final_stop_gained_cases.update( { stops[0]['chr'] + "_%_" + str( stops[0]['pos'] ): stops[0]['content'] } )
		cds_pos_of_first_snpeff_stop_per_gene.update( { gene: int( stops[0]['content'][7].split('|')[12].split('/')[0] ) } )
	return final_stop_gained_cases, all_snpeff_variants, cds_pos_of_first_snpeff_stop_per_gene


def identify_matches( snpeff_stop_gained_cases, navip_vcf_file, verbose, cds_pos_of_first_snpeff_stop_per_gene ):
	"""! @brief identify matches between SnpEff and NAVIP predictions for potential stop_gained cases with two SNVs per codon """
	
	matches = []
	black_list = {}	#count each genomic position only once (reduces redundancy through alternative transcripts)
	navip_exclusive_premature_stops = []
	with open( navip_vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				try:
					black_list[ parts[0] + "_%_" + parts[1] ]
				except KeyError:	#this position was not analyzed previously for a different transcript isoform
					try:
						hit = snpeff_stop_gained_cases[ parts[0] + "_%_" + parts[1] ]	#check if this position matches a SnpEff stop_gain
						matches.append( [ hit, parts ] )
					except KeyError:
						if 'Stop gained' in line:	#check if this position might be a stop_gain only detected by NAVIP
							geneID = line.split('\t')[7].split('|')[3]
							cds_position = int( line.split('\t')[7].split('|')[12].split('/')[0] )
							try:
								if cds_position < cds_pos_of_first_snpeff_stop_per_gene[ geneID ]:	#SnpEff also predicted a stop_gain, but downstream
									navip_exclusive_premature_stops.append( line )
							except KeyError:	#SnpEff did not predict a stop_gain for this gene
								navip_exclusive_premature_stops.append( line )
					black_list.update( { parts[0] + "_%_" + parts[1]: None } )
			line = f.readline()
	
	# --- filter stop gains by gene to report only first stop gain per gene --- #
	stop_gains_per_gene = {}
	for variant in navip_exclusive_premature_stops:
		gene = variant.split('\t')[7].split('|')[3]
		try:
			stop_gains_per_gene[ gene ].append( {'pos': int( variant.split('\t')[1] ), 'chr': variant.split('\t')[0], 'cdspos': int( variant.split('\t')[7].split('|')[12].split('/')[0] ), 'content': variant } )
		except KeyError:
			try:
				stop_gains_per_gene.update( { gene: [ {'pos': int( variant.split('\t')[1] ), 'chr': variant.split('\t')[0], 'cdspos': int( variant.split('\t')[7].split('|')[12].split('/')[0] ), 'content': variant } ] } )
			except IndexError:
				sys.stdout.write( "ERROR: cds position not detected: " + variant[7] + "\n" )
				sys.stdout.flush()
	final_navip_exclusive_premature_stops = []
	for gene in list( stop_gains_per_gene.keys() ):
		stops = stop_gains_per_gene[ gene ]
		if len( stops ) > 1:	#find the first predicted premature stop codon and ignore the following ones
			stops = sorted( stops, key=itemgetter('cdspos') )
		final_navip_exclusive_premature_stops.append( stops[0]['content'] )
	return matches, final_navip_exclusive_premature_stops


def generate_summary_file( summary_file, navip_snpeff_matches, verbose ):
	"""! @brief generate summary file """
	
	# --- collect substitution results --- #
	substitution_results = []
	with open( summary_file, "w" ) as out:
		out.write( "#Chr\tPos\tRef\tAlt\tSnpEff\tNAVIP\n" )
		for each in navip_snpeff_matches:
			chromosome = each[0][0]
			position = each[0][1]
			ref_allele = each[0][3]
			alt_allele = each[0][4]
			snpeff_anno = each[0][7].split('|')[10][2:]
			navip_anno = each[1][7].split('|')[10][2:]
			if snpeff_anno != navip_anno:
				if verbose:
					sys.stdout.write( each[0][7] + "\n" )
					sys.stdout.write( each[1][7] + "\n\n" )
					sys.stdout.flush()
			pos = re.findall( "[0-9]+", navip_anno )[0]
			result = navip_anno.split( pos )[1]
			substitution_results.append( result )
			out.write( "\t".join( [ chromosome, position, ref_allele, alt_allele, snpeff_anno, navip_anno ] ) + "\n" )
	
	# --- summarize substitution results --- #
	counts = {}
	for aa in sorted( list( set( substitution_results ) ) ):
		count = substitution_results.count( aa )
		sys.stdout.write( aa + ": " + str( count ) + "\n" )
		sys.stdout.flush()
		counts.update( { aa: count } )
	return counts


def generate_aa_sub_figure( amino_acid_substitution_figure_file, data ):
	"""! @brief generate a barplot showing the different amino acid substitutions """
	
	values = []
	names = []
	for aa in sorted( list( data.keys() ) ):
		if aa != "*":	#ignore stop codons
			names.append( aa )
			values.append( data[ aa ] )
	
	plt.figure()
	bars = plt.bar( range( len( values ) ), values, color='#1F77B4' )	#color matching NAVIP cInDel output
	
	plt.xlabel('amino acid substitution instead of premature stop codon')
	plt.ylabel('number of events')

	plt.xticks(range(len(names)), names)	#set names of bars

	for bar in bars:
		yval = bar.get_height()
		plt.text(bar.get_x() + bar.get_width()/2, yval + 0.1, round(yval, 2), ha='center', va='bottom')

	plt.tight_layout()  # Optional: Adjust layout to prevent overlap
	plt.savefig( amino_acid_substitution_figure_file, dpi=300 )


def output_navip_exclusive_stop_gains( navip_exclusive_stop_gain_file, navip_exclusive_premature_stops ):
	"""! @brief write NAVIP exclusive predicted premature stop codons into output file """
	
	with open( navip_exclusive_stop_gain_file, "w" ) as out:
		out.write( "\t".join( [ "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT" ] ) + "\n" )
		for each in navip_exclusive_premature_stops:
			parts = each.strip().split('\t')
			if len( parts ) == 7:	#account for NAVIP issue in old version
				parts.append( "." )
			out.write( "\t".join( parts ) + "\n" )


def write_navip_stop_gain_genes_to_file( navip_exclusive_premature_stops, navip_snpeff_matches, navip_stop_codon_genes_file ):
	"""! @brief write all genes with premature stop codons predicted by NAVIP into output file """
	
	navip_stop_gain_genes = []
	for each in navip_exclusive_premature_stops:
		navip_stop_gain_genes.append( each.split('\t')[7].split('|')[3] )	#extraction of gene ID from NAVIP annotation string
	for each in navip_snpeff_matches:
		navip_stop_gain_genes.append( each[1][7].split('|')[3] )	#NAVIP string has index 1 in match list; annotation string has index 7 in parts list; gene ID has index 3 in annotation string
	genes = sorted( list( set( navip_stop_gain_genes ) ) )	#make all entries unique
	with open( navip_stop_codon_genes_file, "w" ) as out:
		out.write( "\n".join( genes ) + "\n" )
	return genes


def write_snpeff_stop_gain_genes_to_file( navip_snpeff_matches, snpeff_stop_codon_genes_file ):
	"""! @brief write all genes with premature stop codons predicted by SnpEff into output file """
	
	snpeff_stop_gain_genes = []
	for each in navip_snpeff_matches:
		snpeff_stop_gain_genes.append( each[0][7].split('|')[3] )	#SnpEff string has index 1 in match list; annotation string has index 7 in parts list; gene ID has index 3 in annotation string
	genes = sorted( list( set( snpeff_stop_gain_genes ) ) )	#make all entries unique
	with open( snpeff_stop_codon_genes_file, "w" ) as out:
		out.write( "\n".join( genes ) + "\n" )
	return genes


def generate_detail_summary_file( complex_comparison_file, navip_exclusive_premature_stops, navip_snpeff_matches, all_snpeff_variants ):
	"""! @brief summarize information about premature stop codons in one file """
	
	data_per_gene = {}
	# --- process NAVIP-specific stop gains --- #
	for each in navip_exclusive_premature_stops:
		chromosome = each.split('\t')[0]
		position = int( each.split('\t')[1] )
		ref_allele = each.split('\t')[3]
		alt_allele = each.split('\t')[4]
		gene = each.split('\t')[7].split('|')[3]	#extraction of gene ID from NAVIP annotation string
		
		spart = all_snpeff_variants[ chromosome + "_%_" + str( position ) ]
		
		scds_change = spart[7].split('|')[9]
		spep_change = spart[7].split('|')[10]
		scds_pos = spart[7].split('|')[12]
		snpep_pos = spart[7].split('|')[13]
		
		snpeff_string = "|".join( [ scds_change, spep_change, scds_pos, snpep_pos ] )
		
		ncds_change = each.split('\t')[7].split('|')[9]
		npep_change = each.split('\t')[7].split('|')[10]
		ncds_pos = each.split('\t')[7].split('|')[12]
		npep_pos = each.split('\t')[7].split('|')[13]
		
		navip_string = "|".join( [ ncds_change, npep_change, ncds_pos, npep_pos ] )
		
		data_per_gene.update( { gene: [ { 'chr': chromosome, 'pos': position, 'ref': ref_allele, 'alt': alt_allele, 'navip': navip_string, 'snpeff': snpeff_string } ] } )
	
	# --- process all other stop gains --- #
	for match in navip_snpeff_matches:
		spart = match[0]	#SnpEff part
		npart = match[1]	#NAVIP part
		
		chromosome = spart[0]
		position = int( spart[1] )
		ref_allele = spart[3]
		alt_allele = spart[4]
		gene = spart[7].split('|')[3]
		
		scds_change = spart[7].split('|')[9]
		spep_change = spart[7].split('|')[10]
		scds_pos = spart[7].split('|')[12]
		spep_pos = spart[7].split('|')[13]
		snpeff_string = "|".join( [ scds_change, spep_change, scds_pos, spep_pos ] )
		
		ncds_change = npart[7].split('|')[9]
		npep_change = npart[7].split('|')[10]
		ncds_pos = npart[7].split('|')[12]
		npep_pos = npart[7].split('|')[13]
		navip_string = "|".join( [ ncds_change, npep_change, ncds_pos, npep_pos ] )
		
		try:
			data_per_gene[ gene ].append( { 'chr': chromosome, 'pos': position, 'ref': ref_allele, 'alt': alt_allele, 'navip': navip_string, 'snpeff': snpeff_string } )
		except KeyError:
			data_per_gene.update( { gene: [ { 'chr': chromosome, 'pos': position, 'ref': ref_allele, 'alt': alt_allele, 'navip': navip_string, 'snpeff': snpeff_string } ] } )
	
	# --- write data into output file --- #
	with open( complex_comparison_file, "w" ) as out:
		out.write( "#NAVIP/SnpEff strings contain: CDS change | PEP change | position in CDS | position in PEP\n" )
		out.write( "\t".join( [ "#Chr", "Pos", "Gene", "Ref", "Alt", "NAVIP", "SnpEff" ] ) + "\n" )
		for gene in sorted( list( data_per_gene.keys() ) ):
			entries = data_per_gene[ gene ]
			for entry in sorted( entries, key=itemgetter('pos') ):
				new_line = [ entry['chr'], str( entry['pos'] ), gene, entry['ref'], entry['alt'], entry['navip'], entry['snpeff'] ]
				out.write( "\t".join( new_line ) + "\n" )


def main( arguments ):
	"""! @brief run everything """
	
	input_vcf_file = arguments[ arguments.index('--snpeffvcf')+1 ]
	navip_vcf_file = arguments[ arguments.index('--navipvcf')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	verbose = False
	
	if output_folder[-1] != "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	snpeff_stop_gained_cases, all_snpeff_variants, cds_pos_of_first_snpeff_stop_per_gene = load_snpeff_stop_gained_cases( input_vcf_file )
	#all_snpeff_variants = dictionary with 'chr + _%_ + pos' as key and list of parts as value
	#cds_pos_of_first_snpeff_stop_per_gene = dictionary with CDS position of SnpEff stop gain assigned to gene ID
	sys.stdout.write( "Number of stop_gained cases in SnpEff results: " + str( len( list( snpeff_stop_gained_cases.keys() ) ) ) + "\n" )
	sys.stdout.flush()
	
	# --- identify matching NAVIP predictions --- #
	navip_snpeff_matches, navip_exclusive_premature_stops = identify_matches( snpeff_stop_gained_cases, navip_vcf_file, verbose, cds_pos_of_first_snpeff_stop_per_gene )
	sys.stdout.write( "Number of premature stop codons only detectable by NAVIP: " + str( len( navip_exclusive_premature_stops ) ) + "\n" )
	sys.stdout.flush()
	
	# --- generate summary output file --- #
	summary_file = output_folder + "summary_based_on_snpeff_stop_gains.txt"
	data = generate_summary_file( summary_file, navip_snpeff_matches, verbose )
	
	# --- generate figure --- #
	amino_acid_substitution_figure_file = output_folder + "aa_sub_figure.png"
	generate_aa_sub_figure( amino_acid_substitution_figure_file, data )
	
	# --- generate output file with NAVIP exclusive premature stop codons --- #
	navip_exclusive_stop_gain_file = output_folder + "navip_exclusive_stop_gains.vcf"
	output_navip_exclusive_stop_gains( navip_exclusive_stop_gain_file, navip_exclusive_premature_stops )
	
	# --- generate NAVIP stop codon gene info file --- #
	navip_stop_codon_genes_file = output_folder + "navip_stop_codon_genes.txt"
	navip_stop_gain_genes = write_navip_stop_gain_genes_to_file( navip_exclusive_premature_stops, navip_snpeff_matches, navip_stop_codon_genes_file )
	
	# --- generate SnpEff stop codon gene info file --- #
	snpeff_stop_codon_genes_file = output_folder + "snpeff_stop_codon_genes.txt"
	snpeff_stop_gain_genes = write_snpeff_stop_gain_genes_to_file( navip_snpeff_matches, snpeff_stop_codon_genes_file )
	
	# --- generate a complex summary file with infos contributed by both tools --- #
	complex_comparison_file = output_folder + "all_genes_with_premature_stop_codons_and_details.txt"
	generate_detail_summary_file( complex_comparison_file, navip_exclusive_premature_stops, navip_snpeff_matches, all_snpeff_variants )


if '--snpeffvcf' in sys.argv and '--navipvcf' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
