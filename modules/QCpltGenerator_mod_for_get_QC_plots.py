# Written by Nina Gilshteyn, M.S.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#print("usage: df, "dfname", "fxn name: QCplt")
#from hist_loop_notfancy import basic_hist_generator
import seaborn as sns 

def QCplt(df, dfname):
    print('max sure the input columns are cells, not genes')
    print(df.head(3))
    plt.figure()
    genes = np.sum(df != 0 , axis = 0) #for each column it says how many times the row does not equal 0 so when cells are columns this is how many rexpress gene in each cell and how  many total counts for the cell 
    counts = np.sum(df, axis = 0)  # sum of each columns 
    genes = pd.DataFrame(genes)
    counts = pd.DataFrame(counts)
    #genes = genes.reset_index(drop=True)
    #counts = counts.reset_index(drop=True)
    #fuck.columns = ['number_of_genes']
    #shit.columns = ['net_transcripts']
     
    print(counts.head(4))
    print(genes.head(4))
 
    cell_barcodes = pd.DataFrame(df.columns)
    #print('head of input', df.head(3))
    print(cell_barcodes.head(4))

    cell_barcodes.index = genes.index
    counts.index = genes.index
    print('isolated cell barcodes', cell_barcodes.head(4))
    QCdata = pd.concat([cell_barcodes,genes,counts],axis = 1) #number genes on the left transcrupots right
    print(QCdata.head(3))
    #number_of_empty_cells = len(QCdata) - len(newQCdata)

    QCdata.columns = ['cellbarcodes','number_of_genes', 'net_transcripts']
    QCdata.to_csv(f'{dfname}_newQCdata_all_cells.csv')
    newQCdata=QCdata[QCdata['number_of_genes']!=0]
    number_of_empty_cells = len(QCdata) - len(newQCdata)

    if number_of_empty_cells != 0:
        newQCdata.to_csv(f'{dfname}_newQCdata_nozero_cellswithout_{number_of_empty_cells}_number_empty_cells.csv')
    else:
        print('number of empty cells is 0')
        print('number of empty cells:', number_of_empty_cells)
    #number_of_empty_cells = len(QCdata) - len(newQCdata)
    N = len(newQCdata)
    a = plt.scatter(newQCdata['number_of_genes'], newQCdata['net_transcripts'])
    plt.xlabel('number of genes')
    plt.ylabel('number of transcripts')
    plt.title(f'{dfname} N = {N}')
    #plt.show()
    if number_of_empty_cells != 0:
        a = a.get_figure()
        a.savefig(f'{dfname}_QC_without_{number_of_empty_cells}_number_empty_cells.png')
        plt.close(a.figure)
    
    else:
        a = a.get_figure()
        a.savefig(f'{dfname}_QC_plot.png')
        plt.close(a.figure)

    binz = int(len(QCdata)**0.5)
    a = sns.displot(newQCdata['number_of_genes'], bins=binz)
    plt.xlabel('number of expressed genes')
    plt.ylabel(f'number of cells N = {N}')
    plt.title(f'{dfname}')

    if number_of_empty_cells != 0:
        a.savefig(f'{dfname}_QC_hist_number_of_genes_with_{number_of_empty_cells}_numberofemptycells_removed.png')
        plt.close(a.figure)
    else:
        a.savefig(f'{dfname}_QC_hist_number_expressed_genes.png')
        plt.close(a.figure)

    
    a = sns.displot(newQCdata['net_transcripts'], bins=binz)
    plt.xlabel('net transcripts')
    plt.ylabel(f'number of cells N ={N}')
    plt.title(f'{dfname}')
    if number_of_empty_cells != 0:
        a.savefig(f'{dfname}_QC_hist_net_transcripts_with_{number_of_empty_cells}_numberofemptycells_removed.png')
        plt.close(a.figure)
    else:
        a.savefig(f'{dfname}_QC_hist_net_transcripts.png')
        plt.close(a.figure)



#A = pd.DataFrame({'cell1':[0,0,0,0],'cell2':[4,0,0,0], 'cell3':[0,0,0,0], 'cell4':[0,2,2,1]}, index = ['gene1','gene2','gene3','gene4'])

#QCplt(A, '11111111testQCpltfxn')
#print(A)
