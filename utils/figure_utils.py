import scanpy as sc

import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.style as style
import matplotlib

import plotly.figure_factory as ff
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

import scipy
import scipy.stats as stats
from scipy.stats import spearmanr
from scipy.stats import zscore
from scipy.stats import pearsonr
from sklearn.metrics import pairwise_distances

# other imports
import tqdm

# utils function for querying pairwise results
def get_pairwise_results(ad_data, cluster_id, filename=None,
                         heatmap_color=[[0,'rgb(250,250,250)'], [1,'rgb(102,0,204)']], heatmap_min_color=0, heatmap_max_color=4, plotly_color_template='plotly_white'):
    '''
        [Summary]
            Creates a plotly stacked plot consisting of a histogram and a heatmap with shared x-axes to view
            pairwise differential expression results from an anndata object.
        [Parameters]
            ad_data
            cluster_id
            filename
            heatmap_color
            heatmap_min_color
            heatmap_max_color
            plotly_color_template
        [Return]
    '''
    # Get pairwise results
    cluster_data = ad_data.varm['pairwise_cluster_count'].loc[:, str(cluster_id)].sort_values(ascending=False).to_frame().copy()
    cluster_data = cluster_data[cluster_data.iloc[:, 0] > 0]
    exp_mat = ad_data.varm['cluster_means'].loc[cluster_data.index]

    # Histogram Plots
    Hist_Counts = go.Bar(
        x=cluster_data.index,
        y=cluster_data.iloc[:, 0],
        marker=dict(
            color='slateblue'
        ),
        opacity=0.6
    )

    # Heatmap Plot
    Matrix_Heatmap = go.Heatmap(
        z=exp_mat.T,
        x=exp_mat.index,
        y=exp_mat.columns,
        colorscale=heatmap_color,
        zmin=heatmap_min_color,
        zmax=heatmap_max_color,
        colorbar={
            'len': 0.68,
            'y': 0.40
        }
    )

    # Put the plot together
    subplot = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.02, row_heights=[0.3, 0.9])
    plots = [Hist_Counts, Matrix_Heatmap]

    # Add plotly graph objects onto the subplots
    for row in range(len(plots)):
        subplot.add_trace(
            plots[row],
            row=row + 1,
            col=1
        )

    # Customize/update layouts
    subplot.update_layout(
        template=plotly_color_template,
        height=600,
        width=1000,
        title=dict(
            text=f'Cluster {cluster_id} | Pairwise Analysis Results -- Cluster Enriched Genes'
        ),
        xaxis=dict(
            showgrid=False
        ),
        yaxis=dict(
            showgrid=False,
            title='Pairwise Comparison<br>Gene Counts'
        ),
        xaxis2=dict(
            showgrid=False,
        ),
        yaxis2=dict(
            showgrid=False,
            title='Mean Gene Expression<br>Across Leiden Clusters'
        ),
        xaxis3=dict(
            showgrid=False,
            title='Cluster Enriched Genes'
        ),
        yaxis3=dict(
            visible=False,
            showgrid=False
        ),
    )

    # show the plots
    subplot.show()
    
    if filename is not None:
            # Save as HTML and PNG files using the base filename
            html_filename = f'{filename}.html'
            png_filename = f'{filename}.png'
            
            subplot.write_html(html_filename)
            subplot.write_image(png_filename)
            
# utils function -- computing dendrogram -- figure 5
def compute_linkage(data_obsm, data_obs_index, data_obs_label, feature_names, 
                    dist_metric='cosine', linkage_dist_metric='euclidean', linkage_method='average'):
    # params: data, feature_names, obs_label, metric, linkage_method
    '''
        [Summary]
            Computes a linkages to be used for constructing dendrograms
        [Parameters]
            data_obsm           : The features and values to be used to compute linkages e.g. gene expression values, PCs
            data_obs_index      : The designated name or index for each of the sample rows in data_obsm
            data_obs_label      : An added categorical label column to -- used for grouping the rows together
            feature_names       : The column names or names of the features
            dist_metric         : The distance metric to use e.g. cosine/euclidean etc. see sklearn.metrics.pairwise_distances documentation
            linkage_dist_metric : The distance metric used for linkage method see scipy.cluster.hierarchy.linkage documentation 
            linkage_method      : Linkage method used to compute linkages e.g. single, average, complete etc.
        [Returns]
            Returns the computed linkages, group_labels to be used for dendrogram as well as the 
            paired distance matrix for each group (unordered/unclustered matrix)
    '''
    # build the data matrix and add associated group labels
    features_matrix = pd.DataFrame(data_obsm, index=data_obs_index, columns=feature_names)
    features_matrix.loc[:, 'clust_labels'] = data_obs_label
    
    # compute the mean of the values per label in the data set and create a new matrix
    features_mean = features_matrix.groupby('clust_labels').mean()
    
    # compute the pairwise distances of the groups in the data
    features_dist = pd.DataFrame(
        pairwise_distances(features_mean, metric=dist_metric), 
        index=features_mean.index, columns=features_mean.index)
    
    # compute the linkages 
    compute_linkage = scipy.cluster.hierarchy.linkage(features_dist, method=linkage_method, metric=linkage_dist_metric)
    
    return compute_linkage, features_mean.index, features_dist

# custom function to add -- Scipy/Seaborn Hierarchical Clustering to anndata by replacing .uns() sc.tl.dendrogram metadata results
# the dendrogram_info keys within the computed sc.tl.dendrogram are the values that are computed when scipy.cluster.hierarchical.dendrogram is computed
def create_dendrogram(data, groupby, use_rep, pairwise_dist_metric='cosine', linkage_metric='euclidean', linkage_method='average'):
    '''
        Returns metadata information that can be used to override scanpy computed
        dendrograms and or create dendrogram from scratch
    '''
    # create a feature matrix
    feature_matrix = pd.DataFrame(data.obsm[use_rep], index=data.obs_names, columns=[i for i in range(1, data.obsm[use_rep].shape[1] + 1)])
    feature_matrix.loc[:,'group_labels'] = data.obs[groupby]
    feature_matrix_mean = feature_matrix.groupby('group_labels').mean() 
    
    # compute the pairwise distancees between the two groups
    pairwise_dist = pd.DataFrame(
        pairwise_distances(feature_matrix_mean, metric=pairwise_dist_metric),
        index=feature_matrix_mean.index,
        columns=feature_matrix_mean.index
    )
    
    # compute the linakges using the pairwise distances
    linkage = scipy.cluster.hierarchy.linkage(
        pairwise_dist, method=linkage_method, metric=linkage_metric)
    
    # using the linkage results, compute the dendrogram
    labels=feature_matrix_mean.index
    dend = scipy.cluster.hierarchy.dendrogram(linkage, labels=labels, no_plot=True)
    
    # extract infromation and pack into a dictionary
    packed_dict = {
        'linkage': linkage,
        'groupby': [groupby],
        'use_rep': use_rep,
        'cor_method': pairwise_dist_metric,
        'linkage_method': linkage_method,
        'categories_ordered':dend['ivl'],
        'categories_idx_ordered': dend['leaves'],
        'dendrogram_info': dend,
        'correlation_matrix':pairwise_dist.iloc[dend['leaves'],dend['leaves']]
    }
    
    return packed_dict

# function used in figure 4 -- to plot transcription factor and transporters
def filter_gene_expression(ad_data, cluster_labels, target_genes, percent_group_threshold=0.4, min_clusters=1, return_mean_exp=False):
    """
        Given an anndata.var column name, filters and compute the mean expression 
        of selected genes per cluster in an anndata object if return means is set to true.
        This is the same filtering set up used in the sheath/socket analyses.

        Parameters:
        - ad_data (anndata.AnnData): An Anndata object containing gene expression data.
        - cluster_labels (str): The column name in ad_data.obs that contains cluster labels.
        - target_genes (str): The column name in ad_data.var that specifies the target genes to analyze.
        - percent_group_threshold (float, optional, default: 0.4): The threshold for gene expression as a percentage.
          Genes with expression above this threshold in a cluster are considered.
        - min_clusters (int, optional, default: 1): The minimum number of clusters in which a gene must exceed the threshold to be selected.
        - return_mean_exp (bool, optional, default: False): If True, the function returns a DataFrame containing the mean expression
          of selected genes per cluster. If False, it returns a list of selected gene names.

        Returns:
        - If return_mean_exp is True:
            expression_means (pandas.DataFrame): A DataFrame containing the mean expression of selected genes per cluster.
        - If return_mean_exp is False:
            gene_sets (list): A list of selected gene names.

        This function filters the gene expression data for the specified target genes, calculates the fraction of cells expressing
        each gene per cluster, selects genes that meet the specified threshold and are expressed in a minimum number of clusters,
        and either computes the mean expression of these selected genes per cluster or returns a list of selected gene names.
    """
    # Filter the data for the specified target genes
    ad_data_sub = ad_data[:, ad_data.var[target_genes]]

    # Initialize a DataFrame to store the fraction of gene expression per cluster
    fraction_gene_exp = pd.DataFrame(0.0, index=ad_data_sub.obs[cluster_labels].values.categories, columns=ad_data_sub.var_names)

    # Calculate the fraction of cells expressing each gene per cluster
    for gene in tqdm.tqdm(ad_data_sub.var_names, total=len(ad_data_sub.var_names), desc='Getting Fractions'):
        gene_expression = np.ravel(ad_data_sub[:, gene].X.todense())
        cluster_counts = ad_data_sub.obs[cluster_labels].value_counts()
        gene_counts = ad_data_sub.obs[cluster_labels][gene_expression > 0].value_counts()
        fraction_gene_exp.loc[:, gene] = list(gene_counts / cluster_counts)

    # Sort columns by the number of clusters with expression above the threshold -- sortubg if the genes!
    fraction_gene_exp = fraction_gene_exp.loc[:, (fraction_gene_exp >= percent_group_threshold).sum().sort_values(ascending=False).index]

    # Select genes with expression above the threshold in at least min_clusters clusters
    gene_sets = fraction_gene_exp.columns[(fraction_gene_exp >= percent_group_threshold).sum() >= min_clusters]
    
    if return_mean_exp:
        # Compute mean expression of selected genes per cluster
        expression_means = pd.DataFrame(ad_data_sub.X.toarray().copy(), columns=ad_data_sub.var_names, index=ad_data_sub.obs_names)
        expression_means.loc[:, cluster_labels] = ad_data_sub.obs.loc[:, cluster_labels].copy()
        expression_means = expression_means.groupby(cluster_labels).mean().loc[:, gene_sets]

        return expression_means
    else:
        return gene_sets


# utils functions for ttest of gene expression between herm and males -- adds annotation
# modified version of: https://stackoverflow.com/questions/67505252/plotly-box-p-value-significant-annotation
def add_p_value_annotation(fig, array_columns, subplot=None, _format=dict(interline=0.07, text_height=1.07, color='black',line_placement=0.2)):
    ''' Adds notations giving the p-value between two box plot data (t-test two-sided comparison)
    
    Parameters:
    ----------
    fig: figure
        plotly boxplot figure
    array_columns: np.array
        array of which columns to compare 
        e.g.: [[0,1], [1,2]] compares column 0 with 1 and 1 with 2
    subplot: None or int
        specifies if the figures has subplots and what subplot to add the notation to
    _format: dict
        format characteristics for the lines

    Returns:
    -------
    fig: figure
        figure with the added notation
    '''
    # Specify in what y_range to plot for each pair of columns
    y_range = np.zeros([len(array_columns), 2])
    for i in range(len(array_columns)):
        y_range[i] = [1.01+i*_format['interline'], 1.02+i*_format['interline']]

    # Get values from figure
    fig_dict = fig.to_dict()

    # Get indices if working with subplots
    if subplot:
        if subplot == 1:
            subplot_str = ''
        else:
            subplot_str =str(subplot)
        indices = [] #Change the box index to the indices of the data for that subplot
        for index, data in enumerate(fig_dict['data']):
            #print(index, data['xaxis'], 'x' + subplot_str)
            if data['xaxis'] == 'x' + subplot_str:
                indices = np.append(indices, index)
        indices = [int(i) for i in indices]
        print((indices))
    else:
        subplot_str = ''

    # Print the p-values
    for index, column_pair in enumerate(array_columns):
        if subplot:
            data_pair = [indices[column_pair[0]], indices[column_pair[1]]]
        else:
            data_pair = column_pair

        # Mare sure it is selecting the data and subplot you want
        #print('0:', fig_dict['data'][data_pair[0]]['name'], fig_dict['data'][data_pair[0]]['xaxis'])
        #print('1:', fig_dict['data'][data_pair[1]]['name'], fig_dict['data'][data_pair[1]]['xaxis'])

        # Get the p-value
        pvalue = stats.ttest_ind(
            fig_dict['data'][data_pair[0]]['y'],
            fig_dict['data'][data_pair[1]]['y'],
            equal_var=False,
        )[1]
        # if pvalue >= 0.05:
        #     symbol = 'ns'
        # elif pvalue >= 0.01: 
        #     symbol = '*'
        # elif pvalue >= 0.001:
        #     symbol = '**'
        # elif pvalue >= 0.0001:
        #     symbol = '***'
        # else:
        #     symbol = '****'
            
# should probably re run this iw th the following:
        if pvalue >= 0.05: 
            symbol = 'ns' 
        elif pvalue < 0.05 and pvalue >= 0.01: 
            symbol = '*' 
        elif pvalue < 0.01 and pvalue >= 0.001: 
            symbol = '**' 
        elif pvalue < 0.001 and pvalue >= 0.0001: 
            symbol = '***' 
        else: 
            symbol = '****'

        # Vertical line
        fig.add_shape(type="line",
            xref="x"+subplot_str, yref="y"+subplot_str+" domain",
            x0=column_pair[0], y0=y_range[index][0] + _format['line_placement'], 
            x1=column_pair[0], y1=y_range[index][1] + _format['line_placement'],
            line=dict(color=_format['color'], width=2,)
        )
        # Horizontal line
        fig.add_shape(type="line",
            xref="x"+subplot_str, yref="y"+subplot_str+" domain",
            x0=column_pair[0], y0=y_range[index][1] + _format['line_placement'], 
            x1=column_pair[1], y1=y_range[index][1] + _format['line_placement'],
            line=dict(color=_format['color'], width=2,)
        )
        # Vertical line
        fig.add_shape(type="line",
            xref="x"+subplot_str, yref="y"+subplot_str+" domain",
            x0=column_pair[1], y0=y_range[index][0] + _format['line_placement'], 
            x1=column_pair[1], y1=y_range[index][1] + _format['line_placement'],
            line=dict(color=_format['color'], width=2,)
        )
        ## add text at the correct x, y coordinates
        ## for bars, there is a direct mapping from the bar number to 0, 1, 2...
        fig.add_annotation(dict(font=dict(color=_format['color'],size=20),
            x=(column_pair[0] + column_pair[1])/2,
            y=y_range[index][1]*_format['text_height'] + _format['line_placement'],
            showarrow=False,
            text=symbol,
            textangle=0,
            xref="x"+subplot_str,
            yref="y"+subplot_str+" domain"
        ))
    return fig

# utils function for ttest of gene expression between herm and males -- creates, render and saves the plot
# dependent on add_p_value_annotation
def sexes_gene_expression_comparison(ad_data, gene_name, herm_color='blue', male_color='red', custom_form={'interline': 0.07, 'text_height': 1.07, 'color': 'gray', 'line_placement': -0.1}, filename=None):
    # Create a DataFrame for gene expression and sex information
    gene_expression = pd.DataFrame(ad_data[:, gene_name].X.toarray(),
                                    index=ad_data.obs_names,
                                    columns=[gene_name])
    gene_expression.loc[:, 'sex'] = ad_data.obs['sex'].values.copy()

    
    
    # Plot the distributions
    color_map = {'Hermaphrodite': herm_color, 'Male': male_color}
    max_val = gene_expression.loc[:, gene_name].max()
    fig = px.violin(data_frame=gene_expression, x='sex', y=gene_name,
                    color='sex', color_discrete_map=color_map, points=False)
    
    # plot aesthetics
    fig.update_layout(
        width=700,
        height=900,
        plot_bgcolor='rgba(0,0,0,0)',
        yaxis=dict(title=f'<i>{gene_name}</i> expression'),
        showlegend=True,
        title=dict(text=f'<i><b>{gene_name}</b></i> Expression', font=dict(size=25)),
        font=dict(family='Arial', color='black'),
        margin=dict(pad=20),
    )

    fig.update_traces(box_visible=False, meanline_visible=True,
                      opacity=0.9, marker_line_color='rgba(0,0,0,0)', line=dict(width=0),
                      marker_line_width=0.5, showlegend=False, width=0.98, marker_size=4)
    
    fig.update_xaxes(showgrid=False, tickangle=90, tickfont=dict(family='Arial', size=18),
                     showline=True, gridcolor='lightgray', linecolor='black')
    fig.update_yaxes(showgrid=False, tickfont=dict(family='Arial', size=20),
                     showline=True, gridcolor='lightgray', range=[-0.001, max_val + 0.5], linecolor='black')

    # Perform t-test and add p-value annotation
    herm_values = gene_expression[gene_expression['sex'] == 'Hermaphrodite'][gene_name]
    male_values = gene_expression[gene_expression['sex'] == 'Male'][gene_name]
    # _, p_value = stats.ttest_ind(herm_values, male_values)

    fig = add_p_value_annotation(fig=fig, array_columns=[[0, 1]], _format=custom_form)

    fig.show()
    if filename:
        fig.write_image(f'{filename}.png')
        fig.write_html(f'{filename}.html')
