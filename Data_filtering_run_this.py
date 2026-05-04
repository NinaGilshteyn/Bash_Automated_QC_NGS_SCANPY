# Written by Nina Gilshteyn, M.S.

from filter_mod import filtercells
from NEW_QC_mod import get_2D_QC_viz
from NEW_QC_mod import get_3D_QC_viz
import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(description = 'input filename')
parser.add_argument('filename', help = 'dfname in directory')
args = parser.parse_args()

### read in file name in directory and edit file-name string
QCfile = args.filename
df_file = QCfile[16:]
df_file = df_file 
print(df_file)
outputname = df_file[:-4]

# set parameters based on file 
if df_file == 'K6_parsebio_lungfat__downsampled_from_84677_number_cells.csv':
    ### set Gene Cutoffs 
    gene_cutoff_bigger_number = 4500
    gene_cutoff_smaller_number = 0
    gene_cutoff_only_number = 'NA'

    ### set count cutoffs
    count_cutoff_bigger_number = 20000
    count_cutoff_smaller_number = 0
    count_cutoff_only_number = 'NA'

    ###   set mitocounts
    mito_cutoff_bigger_number = 'NA'
    mito_cutoff_smaller_number = 'NA'
    mito_cutoff_only_number = 30

    ### determine if cutoffs are inclusive or exclusive based on documentation
    gene_filter_switch = 'genes_exclusive'
    count_filter_switch = 'transcripts_exclusive'
    mito_filter_swtich = 'mito_exclusive'

elif df_file == 'K1_parsebio_lungfat__downsampled_from_133849_number_cells.csv':
    
    ### set Gene Cutoffs 
    gene_cutoff_bigger_number = 5000
    gene_cutoff_smaller_number = 0
    gene_cutoff_only_number = 'NA'

    ### set count cutoffs
    count_cutoff_bigger_number = 20000
    count_cutoff_smaller_number = 0
    count_cutoff_only_number = 'NA'

    ###   set mitocounts
    mito_cutoff_bigger_number = 'NA'
    mito_cutoff_smaller_number = 'NA'
    mito_cutoff_only_number = 30  # many have 2000 expressed genes and 42% mito so above 43 is likely the dead cells 

    ### determine if cutoffs are inclusive or exclusive based on documentation
    gene_filter_switch = 'genes_exclusive'
    count_filter_switch = 'transcripts_exclusive'
    mito_filter_swtich = 'mito_exclusive'


### read in files and format them 

df1 = pd.read_csv(df_file, header=0, index_col=0)
QC_df1 = pd.read_csv(QCfile, header = 0 , index_col = 0)
if list(df1.index) == list(QC_df1.index): 
    df1 = pd.concat([QC_df1, df1], axis =1)
else: 
    print('indices are not equal. rewrite program to merge')

#### input files and parameters into function 

### for Gene number 
df2 = filtercells(x=df1, filter_switch=gene_filter_switch, bigger_number=gene_cutoff_bigger_number, smaller_number=gene_cutoff_smaller_number, only_number=gene_cutoff_only_number,number_of_filtration_step=1, dfname_description_of_what_filtered=outputname, display_switch=-1)
### for transcript number 
df3 = filtercells(x=df2, filter_switch=count_filter_switch, bigger_number=count_cutoff_bigger_number, smaller_number=count_cutoff_smaller_number, only_number=count_cutoff_only_number,number_of_filtration_step=2, dfname_description_of_what_filtered=outputname, display_switch=-1)
### for  mito percent 
df4 = filtercells(x=df3, filter_switch=mito_filter_swtich, bigger_number=mito_cutoff_bigger_number, smaller_number=mito_cutoff_smaller_number, only_number=mito_cutoff_only_number,number_of_filtration_step=3, dfname_description_of_what_filtered=outputname, display_switch=-1)

#filtercells(x, filter_switch, bigger_number, smaller_number, only_number,number_of_filtration_step, dfname_description_of_what_filtered, display_switch)
####################### output results 
output1 = 'POST_QC_purified_gex_with_QC_stats_' + df_file
df4.to_csv(output1)

print(df4.head(4)) 
a = df4['n_genes_by_counts']
b = df4['total_counts']
c = df4['pct_counts_mt']
updatedQC = pd.concat([a,b,c], axis =1)
updatedQC.columns = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
updatedQC.index = df4.index
output3 = 'POST_QC_QC_stats_from_Purified_gex' + df_file
updatedQC.to_csv(output3)


df4 = df4.drop('n_genes_by_counts', axis=1)
df4 = df4.drop('total_counts', axis=1)
df4 = df4.drop('pct_counts_mt', axis=1)
#df4 = df4.drop('barcodes', axis=1)
#df4 = df4.T

output2 = 'POST_QC_purified_gex_in_standard_format_'+df_file
df4.to_csv(output2)


###################### visualized new QC


COLORS = ['lightcoral', 'mediumseagreen', 'lightsteelblue']
counter = 0
for col_name in updatedQC:
    col_data = updatedQC[col_name]
    print('counter:', counter)
    COLOR = COLORS[counter]
    binz = int(len(updatedQC)**0.5)
    a = sns.displot(col_data, bins = binz, log_scale = (False,True), color = COLOR, edgecolor = 'k', height=6, aspect=1.2)
    a.savefig(f'POST_QC_{col_name}_{outputname}_histogram_log_scale.png')
    a.savefig(f'POST_QC_{col_name}_{outputname}_histogram_log_scale.eps')
    plt.close(a.figure) 
    a = sns.displot(col_data, bins = binz, log_scale = (False,False), color = COLOR, edgecolor = 'k', height=6, aspect=1.2)
    a.savefig(f'POST_QC_{col_name}_{outputname}_histogram_linear_scale.png')
    a.savefig(f'POST_QC_{col_name}_{outputname}_histogram_linear_scale.eps')

    counter += 1
    plt.close(a.figure)

POSTQC_outputname = 'POST_QC_' + outputname
get_3D_QC_viz(x_vals=updatedQC['n_genes_by_counts'], y_vals=updatedQC['total_counts'], z_vals=updatedQC['pct_counts_mt'], name_of_dataset=POSTQC_outputname)
get_2D_QC_viz(x_vals=updatedQC['n_genes_by_counts'], y_vals=updatedQC['total_counts'], z_vals=updatedQC['pct_counts_mt'], name_of_dataset=POSTQC_outputname)
