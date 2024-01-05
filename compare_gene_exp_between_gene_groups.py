### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.1 ###


__usage__ = """
					python3 compare_gene_exp_between_gene_groups.py
					--genes <SnpEff_VCF_OUTPUT_FILE>
					--exp <NAVIP_VCF_OUTPUT_FILE>
					--out <OUTPUT_FIGURE_FILE>
					
					optional:
					--gff <NAVIP_VCF_OUTPUT_FILE>
					"""

import os, sys, re
from operator import itemgetter
import plotly.io as pio
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import mannwhitneyu
import statistics

# --- end of imports --- #

def load_gene_IDs( genes_info_file ):
	"""! @brief load gene IDs from input file """
	
	genes = {}
	with open( genes_info_file, "r" ) as f:
		lines = f.read().strip().split('\n')
		for line in lines:
			genes.update( { line: None } )
	return genes


def load_gene_IDs_from_gff_file( gff_file, IDtag="ID" ):
	"""! @brief load gene IDs from given GFF3 file """
	
	gene_IDs = []
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "gene":
					gene = parts[-1].split( IDtag + '=')[1]
					if ";" in gene:
						gene = gene.split(';')[0]
					gene_IDs.append( gene )
			line = f.readline()
	return gene_IDs


def load_exp_per_gene( exp_file ):
	"""! @brief load average expression value per gene """
	
	avg_exp_per_gene = {}
	with open( exp_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			geneID = parts[0]
			if "." in geneID:	#if representative transcript is used
				geneID = geneID.split('.')[0]
			avg_exp_per_gene.update( { geneID: float( parts[1] ) } )
			line = f.readline()
	return avg_exp_per_gene


def plot_exp_of_groups( avg_exp_per_gene, candidate_genes, figure_file, supplied_genes, display_cutoff ):
	"""! @brief plot dN/dS ratios of defined genes against all other genes """
	
	group1_values = []
	group2_values = []
	if len( supplied_genes ) > 0:
		for gene in supplied_genes:
			try:
				candidate_genes[ gene ]
				try:
					group1_values.append( min( [ avg_exp_per_gene[ gene ], display_cutoff ] ) )
				except KeyError:
					group1_values.append( 0 )	#if gene ID not in exp data, gene does not have any expression
			except KeyError:
				try:
					group2_values.append( min( [ avg_exp_per_gene[ gene ], display_cutoff ] ) )
				except KeyError:
					group2_values.append( 0 )	#if gene ID not in exp data, gene does not have any expression
	else:
		for gene in list( avg_exp_per_gene.keys() ):
			try:
				candidate_genes[ gene ]
				group1_values.append( min( [ avg_exp_per_gene[ gene ], display_cutoff ] ) )
			except KeyError:
				group2_values.append( min( [ avg_exp_per_gene[ gene ], display_cutoff ] ) )
	
	group1_values_df = pd.DataFrame( group1_values, columns=["group1"] )
	group2_values_df = pd.DataFrame( group2_values, columns=["group2"] )
	
	ng1 = len( group1_values )
	ng2 = len( group2_values )
	
	fig = go.Figure()
	
	legend1 = 'genes with stop_gained\n(n=' + str( ng1 ) + ")"
	fig.add_trace(go.Violin( y=group1_values_df['group1'], name=legend1, showlegend=False, meanline=dict(visible=True, color='black') ) )
	
	legend2 = 'other genes\n(n=' + str( ng2 ) + ")"
	fig.add_trace(go.Violin( y=group2_values_df['group2'], name=legend2, showlegend=False, meanline=dict(visible=True, color='black') ) )
	
	fig.update_layout(violingap=0, violinmode='overlay')
	fig.update_yaxes(range=[0, display_cutoff*1.1])
	fig.update_layout(yaxis_title="average expression per gene")
	
	#fig.update_layout(xaxis=dict(tickfont=dict(size=20)),  #change x-tick fontsize
    #              yaxis=dict(titlefont=dict(size=20)))  #change y-axis labele fontsize

	pio.write_image( fig, figure_file, scale=10 )
	return group1_values, group2_values


def main( arguments ):
	"""! @brief run everything """
	
	genes_file = arguments[ arguments.index('--genes')+1 ]
	exp_file = arguments[ arguments.index('--exp')+1 ]
	output_figure_file = arguments[ arguments.index('--out')+1 ]
	
	display_cutoff=50
	
	if '--gff' in arguments:
		gff_file = arguments[ arguments.index('--gff')+1 ]
		supplied_genes = load_gene_IDs_from_gff_file( gff_file, IDtag="ID" )
	else:
		supplied_genes = []
	
	# --- load data --- #
	candidate_gene_IDs = load_gene_IDs( genes_file )
	avg_exp_per_gene = load_exp_per_gene( exp_file )
	
	# --- generate plot --- #
	g1, g2 = plot_exp_of_groups( avg_exp_per_gene, candidate_gene_IDs, output_figure_file, supplied_genes, display_cutoff )
	
	# -- run statistical test --- #
	u_value, p_value = mannwhitneyu( g1, g2, alternative='two-sided' )
	print( "U: " + str( u_value ) + "\tp-value: " + str( p_value ) )
	print( "group1: mean=" + str( statistics.mean( g1 ) ) + "\tmedian=" + str( statistics.median( g1 ) ) )
	print( "group2: mean=" + str( statistics.mean( g2 ) ) + "\tmedian=" + str( statistics.median( g2 ) ) )


if '--genes' in sys.argv and '--exp' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
