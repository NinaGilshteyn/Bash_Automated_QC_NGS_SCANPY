# Written by Nina Gilshteyn, M.S.

import plotly.graph_objects as go
#load required libraries
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns


def get_hists_and_input_to_multi_dimension_viz_function(data,dfname):
# verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.settings.verbosity = 3
    print('make sure columns are gene names')
    adata = sc.AnnData(data)
    print('adata')
    print(adata)
    # annotate mitochondrial genes as 'mt' and calculate qc metrics
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    stats_output = pd.concat([adata.obs["n_genes_by_counts"],adata.obs["total_counts"],adata.obs["pct_counts_mt"]], axis =1)
    stats_output.columns = ["n_genes_by_counts","total_counts","pct_counts_mt"]
    stats_output.to_csv(f'QC_stats_output_{dfname}.csv')

    
    print('genes by counts')
    print(adata.obs["n_genes_by_counts"])
    print('total counts')
    print(adata.obs["total_counts"])
    print('pct_counts_mt')
    print(adata.obs["pct_counts_mt"])

    print('stats_output')
    print(stats_output)
    binz = int(len(stats_output)**0.5)
    COLORS = ['hotpink', '#39FF14', 'aqua']
    counter = 0
    for col_name in stats_output:
        col_data = stats_output[col_name]
        print('counter:', counter)
        COLOR = COLORS[counter]
        a = sns.displot(col_data, bins = binz, log_scale = (False,True), color = COLOR, edgecolor = 'k', height=6, aspect=1.2)
        a.savefig(f'{col_name}_{dfname}_histogram_log_scale.png')
        plt.close(a.figure)
        a = sns.displot(col_data, bins = binz, log_scale = (False,False), color = COLOR, edgecolor = 'k', height=6, aspect=1.2)
        a.savefig(f'{col_name}_{dfname}_histogram_linear_scale.png')
        counter += 1
        plt.close(a.figure)
    return stats_output 

def get_3D_QC_viz(x_vals, y_vals, z_vals, name_of_dataset):   
    # Create figure
    print('make sure x vals are gene expression, yvals are tot count, and z is mito')
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=x_vals,
        y=y_vals,
        z=z_vals,
        mode='markers+text',
        marker=dict(
            color = z_vals,
            colorscale='Rainbow',         # Use rainbow color scale
            colorbar=dict(title=f'Mito Percent:{name_of_dataset}',titlefont=dict(color='white'),tickfont=dict(color='white')),
            cmin=min(z_vals),
            cmax=max(z_vals),
            sizemode='diameter',
            sizeref=2,
            sizemin=2
        ),
    ))

    # Customize layout for dark grey background
    fig.update_layout(
        scene=dict(
            xaxis=dict(title = 'Number Expressed Genes',backgroundcolor='rgb(40, 40, 40)', color='white', gridcolor='gray'),
            yaxis=dict(title = 'Total RNA', backgroundcolor='rgb(40, 40, 40)', color='white', gridcolor='gray'),
            zaxis=dict(title = 'percent mitochondrial counts', backgroundcolor='rgb(40, 40, 40)', color='white', gridcolor='gray'),
        ),
        paper_bgcolor='rgb(30, 30, 30)',
        plot_bgcolor='rgb(30, 30, 30)',
    )

    # Save plots
    fig.write_html(f"QC_3D_plot_{name_of_dataset}.html")
    fig.write_image(f"QC_3D_plot_gene_{name_of_dataset}.png")


def get_2D_QC_viz(x_vals, y_vals, z_vals, name_of_dataset):  
    plt.scatter(x_vals, y_vals, c=z_vals, cmap='viridis', vmin=0)  # vmin ensures color starts at 0
    plt.colorbar(label='Percent Mitochondrial Counts')  # Optional colorbar
    plt.xlabel('Number of Expressed Genes')
    plt.ylabel('Total Counts')
    plt.title(f'{name_of_dataset}')
    plt.savefig(f'QC_2D_scatter_Colored_by_mito_percent_{name_of_dataset}.png')
    plt.close()
