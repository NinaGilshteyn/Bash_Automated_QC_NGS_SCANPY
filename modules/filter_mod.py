# Written by Nina Gilshteyn, M.S.

import pandas as pd
import numpy as np
import seaborn as sns
#from QCpltGenerator_mod_for_get_QC_plots
#number of genes,net transcripts
import matplotlib.pyplot as plt

def filtercells(x, filter_switch, bigger_number, smaller_number, only_number,number_of_filtration_step, dfname_description_of_what_filtered, display_switch):

### option to display impure stats 

    if display_switch > 0:
        binz = int(len(x)**0.5)
        #"n_genes_by_counts","total_counts","pct_counts_mt"

        plt.title('number of expressed genes')
        plt.suptitle('unpurified data')        
        d = sns.displot(x['n_genes_by_counts'], bins = binz)
        plt.show()
        plt.close()
        d = sns.displot(x['total_counts'],bins = binz)
        plt.title('transcripts')   
        plt.suptitle('unpurified data')
        plt.show()
        plt.close()
        d = sns.displot(x['pct_counts_mt'],bins = binz)
        plt.title('percent counts that are mitochonrial DNA')   
        plt.suptitle('unpurified data')
        plt.show()
    else: 
        print('Filtering Data :D')

#'genes_exclusive', 'transcripts_exclusive', 'mito_exclusive','genes_lower_inclusive', 'transcripts_lower_inclusive', 'genes_upper_inclusive', 'transcripts_upper_inclusive', 

#### identify indices of residual cells with the following logic statement
# "n_genes_by_counts","total_counts","pct_counts_mt"
#################
    if filter_switch == 'genes_exclusive':
        dff = x.query(f'n_genes_by_counts > {bigger_number} | n_genes_by_counts < {smaller_number}') 
        print('bounds')
        #print(lowerbound)
        print('indicides that meet conditions')
        print(dff)
        # these are the indexs that meet the conditions and are being filtered out 
        logic_statment = f'{smaller_number}<n_genes_by_counts<{bigger_number}'
##################
    elif filter_switch == 'transcripts_exclusive':
        dff = x.query(f'total_counts > {bigger_number} | total_counts < {smaller_number}')
        logic_statment = f'{smaller_number}<total_counts<{bigger_number}'
        print('indicies that meet conditions')
        print(dff)
##################
    elif filter_switch == 'mito_exclusive':
        dff = x.query(f'pct_counts_mt > {only_number}')
        logic_statment = f'pct_counts_mt<{only_number}'
##################
    elif filter_switch == 'mito_exclusive2':
        dff = x.query(f'pct_counts_mt < {only_number}')
        logic_statment = f'pct_counts_mt<{only_number}'

##################
    elif filter_switch == 'mito_inclusive':
        dff = x.query(f'pct_counts_mt >= {bigger_number}')
        logic_statment = f'pct_counts_mt<={bigger_number}'
##################
    elif filter_switch == 'genes_lower_inclusive':
        dff = x.query(f'n_genes_by_counts >= {bigger_number} | n_genes_by_counts < {smaller_number}') 
        # these are the indexs that meet the conditions and are being filtered out 
        logic_statment = f'{smaller_number}<n_genes_by_counts<{bigger_number}'
##################
    elif filter_switch == 'transcripts_lower_inclusive':
        dff = x.query(f'total_counts >= {bigger_number} | total_counts < {smaller_number}')
        logic_statment = f'{smaller_number}<net_transcripts<{bigger_number}'
        print('indicies that meet conditions')
        print(dff)
        print('bounds')
        print(bigger_number)
        print('input')
        print(x)

##################
    elif filter_switch == 'genes_upper_inclusive':
        dff = x.query(f'n_genes_by_counts > {bigger_number} | n_genes_by_counts <= {smaller_number}') 
        # these are the indexs that meet the conditions and are being filtered out 
        logic_statment = f'{smaller_number}<n_genes_by_counts<{bigger_number}'
##################
    elif filter_switch == 'transcripts_upper_inclusive':
        dff = x.query(f'total_counts > {bigger_number} | total_counts <= {smaller_number}')
        logic_statment = f'{smaller_number}<net_transcripts<={bigger_number}'
##################
    elif filter_switch == 'genes_inclusive':
        dff = x.query(f'n_genes_by_counts >= {bigger_number} | n_genes_by_counts <= {smaller_number}') 
        # these are the indexs that meet the conditions and are being filtered out 
        logic_statment = f'{smaller_number}<n_genes_by_counts<={bigger_number}'
##################
    elif filter_switch == 'transcripts_inclusive':
        dff = x.query(f'total_counts >= {bigger_number} | total_counts <= {smaller_number}')
        logic_statment = f'{smaller_number}<net_transcripts<={bigger_number}'
##################
    elif filter_switch == 'genes_specific':
        dff = x.query(f'n_genes_by_counts == {only_number}')
        logic_statment = f'n_genes_by_counts!={only_number}'
        print(logic_statment)
##################
    elif filter_switch == 'transcripts_specific':
        dff = x.query(f'total_counts == {only_number}')
        logic_statment = f'total_counts!={only_number}'
        print(logic_statment)

#### write which residual cells were removed to disk 
    
    print('returning filtered data and writing removed cells to disk')
    dff.to_csv(f'{dfname_description_of_what_filtered}__when_{logic_statment}_residual_cells_after_filtering_step_number_{number_of_filtration_step}.csv')
    
#### identify indicies of pure cells as different than the chosen residual cells 

    df = x.loc[~x.index.isin(dff.index)]  ### Find's the index that is not in the filtrate index the list 
   
    #print('Cells that do not match this criteria', dff)
    # .loc allows to acccess group row and col labels
    # X.index will be the row index
    # isin will look for each row index that is in x and df and return a true in a boolean array
    # but the tilde in front means that it is not. So it will filter the x index rows that are not in df in$

#### visualize QC stats of purified cells and write to disk

    binz = int(len(df)**0.5) # N is amount of cells
    a = sns.displot(df['n_genes_by_counts'], bins = binz)
    plt.xlabel('number of expressed genes')
    plt.title(f'{dfname_description_of_what_filtered} filter step: {number_of_filtration_step}')
    plt.suptitle(f'{logic_statment}')
    a.savefig(f'{dfname_description_of_what_filtered}_number_expressed_genes_hist_when_{logic_statment}_after_filtering_step_{number_of_filtration_step}.png')
    #plt.show()
    plt.close()
    b = sns.displot(df['total_counts'], bins = binz)
    plt.xlabel('net transcripts')
    plt.title(f'{dfname_description_of_what_filtered} filter step: {number_of_filtration_step}')
    plt.suptitle(f'{logic_statment}')
    b.savefig(f'{dfname_description_of_what_filtered}_net_transcripts__hist_when_{logic_statment}_after_filtering_step_{number_of_filtration_step}.png')
    #plt.show()
    b = sns.displot(df['pct_counts_mt'], bins = binz)
    plt.xlabel('net transcripts')
    plt.title(f'{dfname_description_of_what_filtered} filter step: {number_of_filtration_step}')
    plt.suptitle(f'{logic_statment}')
    b.savefig(f'{dfname_description_of_what_filtered}_net_transcripts_hist_when__{logic_statment}_after_filtering_step_{number_of_filtration_step}.png')
    #plt.show()

    return df



### test on dummy to filter cell3 and cell4 

#a = pd.DataFrame({'Cell1':[1,1,1,5],'Cell2':[2,2,2,1], 'Cell3':[0,0,0,0], 'Cell4':[5,5,5,0]}, index = ['gene1','gene2', 'gene3', 'gene4'])
#a = a.T
#print(a)

#QC_data = pd.DataFrame({'n_genes_by_counts':[4,4,0,3], 'total_counts':[8,7,0,15], 'pct_counts_mt':[30,20,30,10]}, index=['cell1', 'cell2','cell3','cell4'])

#QC_data.index= a.index

#dfff = pd.concat([QC_data, a], axis = 1)
#print(dfff)
 
# here is the usage : (x, filter_switch, bigger_number, smaller_number, only_number,number_of_filtration_step, dfname_description_of_what_filtered, display_switch)

# let's filter cell 3 

#b = filtercells(x=dfff, filter_switch='genes_specific', bigger_number='NA', smaller_number='NA', only_number =0, number_of_filtration_step=1, dfname_description_of_what_filtered='dummy_empty_cell3', display_switch=-1)


#print('removed cell3?OB') 
#print(b)

# let's filter out cell4

#c = filtercells(x=b, filter_switch='transcripts_lower_inclusive', bigger_number=14, smaller_number=0, only_number ='NA', number_of_filtration_step=2, dfname_description_of_what_filtered='dummy_doublet', display_switch=-1)

#print('removed cell4?')
#print(c)



