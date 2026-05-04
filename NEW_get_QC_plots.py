# Written by Nina Gilshteyn, M.S.

import pandas as pd
from QCpltGenerator_mod_for_get_QC_plots import QCplt  # old QC function

### my new QC library ####

from NEW_QC_mod import get_hists_and_input_to_multi_dimension_viz_function
from NEW_QC_mod import get_3D_QC_viz
from NEW_QC_mod import get_2D_QC_viz
##########################

import argparse


parser = argparse.ArgumentParser(description = 'input filename')
parser.add_argument('filename', help = 'dfname in directory')
args = parser.parse_args()


file = args.filename

print(file)

df = pd.read_csv(file, header = 0, index_col=0)

df=df.T
print('double check if col names are barcodes')
print(df.head(3))

#def QCplt(df, dfname):

DFNAME = file[:-4]
print(DFNAME)
DFNAME_2 = DFNAME + '_' + 'old_QC_mod_double_check'
QCplt(df,dfname=DFNAME_2)


#####  stuff from new library 
df=df.T
print('double check if col names are genes')
print(df.head(3))

QCstats = get_hists_and_input_to_multi_dimension_viz_function(data=df,dfname=DFNAME)

#"n_genes_by_counts","total_counts","pct_counts_mt"

get_3D_QC_viz(x_vals= QCstats["n_genes_by_counts"], y_vals= QCstats["total_counts"], z_vals= QCstats["pct_counts_mt"], name_of_dataset=DFNAME)
get_2D_QC_viz(x_vals= QCstats["n_genes_by_counts"], y_vals= QCstats["total_counts"], z_vals= QCstats["pct_counts_mt"], name_of_dataset=DFNAME)













