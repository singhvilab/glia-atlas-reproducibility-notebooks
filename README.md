# worm-glia-atlas
This is a repository with tutorials notebooks for the worm glia scRNA-seq atlas available at https://wormglia.org.

## Environment Setup
------------------------------------------
To run the tutorial, clone this repository and follow the environment set up below.
```
git clone https://github.com/settylab/worm-glia-atlas.git
```

##### Option I: Creating a `conda` environment using the provided yaml file:
```
envName=worm-glia-atlas

conda env create -n "$envName" --file envs/environment.yaml

conda activate "$envName"
```
##### Option II: Create an environment and install the required packages using `pip`:
```
envName=<your-environment-name>

conda create -n "$envName" python=3.8.10 pip=21.1.3

conda activate "$envName"

pip install -r envs/requirements.txt
```

### Important Dependencies
The following packages need to be installed for running the notebooks, which can also be installed by following the environment setup above:
1. `scanpy`: https://scanpy.readthedocs.io/en/stable/installation.html
4. `sklearn`: https://scikit-learn.org/stable/install.html
2. `plotly`: https://plotly.com/python/getting-started/
3. `tqdm` : https://github.com/tqdm/tqdm#installation


## Pairwise Differential Expression Analysis
------------------------------------------
The tutorial notebook for <b><i>pairwise differential expression analysis</i></b> is available <b>[here](https://github.com/settylab/worm-glia-atlas/blob/main/notebooks/pairwise-differential-results.ipynb)</b>.
Pairwise differential analysis performs pairwise differential analysis to identify cluster enriched genes rather than one-vs-all approach. 

#### <b>Inputs</b>
The input is an `anndata` object with normalized, log-transformed data and an `obs` variable containing information about clusters for pairwise comparison.
1. `anndata.obs[LEIDEN_NAME]`: Where `LEIDEN_NAME` is the name of the column in `anndata.obs` field containing the groups to be used for the anlaysis (leiden clusters in this analysis).

#### <b>Outputs</b>
The `anndata` object is updated with the following information
1. `anndata.varm['pairwise_cluster_count']`: Gene X Cluster matrix indicating how many comparisons the gene is differential in.
2. `anndata.varm['cluster_means']`: Gene X Cluster matrix of mean expression of gene per cluster.
3. HTML files of the pairwise analyses results can saved using the `plot_pairwise_results()` function by specifying a path to the `save` parameter.

## Feature Ranking Analysis 
------------------------------------------
The tutorial notebooks for <b><i>Sheath/Socket</i> & <i>Pan-Glia</i> marker analysis</b> is available <b>here</b>. A `logistic regression` model is trained and employed for `binary classification` of cells using gene expression. 

Subsequently, a ranking of the learned features within the model is then performed with the objective being to rank features that are highly informative and correlated with the specified target classes or cell type. 

These analyses can be readily extended to other datasets by providing the appropriate inputs, as outlined below.

#### <b>Inputs</b>
The key inputs for this analyses is the `anndata` object with normalized & log-transformed counts, imputed gene expression values as well as the following anndata attribute fields below:

- `anndata.obs[CLASS_LABELS]`: Where `CLASS_LABELS` is the name of the column in `anndata.obs` containing the ground truth labels for each cells in the anndata object.
- `anndata.obs[CLUSTER_LABELS]`: Where `CLUSTER_LABELS` is the name of the column in `anndata.obs` containing the cluster labels for each cells in the anndata object.
- `anndata.var[USE_GENES]`: Where `USE_GENES` is the name of the column in `anndata.var` containing boolean values that specifies whether a gene is to be used for analysis or ignored (default is `highly_variable` genes columns). 
- `anndata.layers[USE_LAYER]`: Where `USE_LAYER` is a key in `anndata.layers` dictionary corresponding to the imputed Cell X Gene count matrix. If not specified, will use the normalized and log-transformed counts as values for the constructed feature matrix & feature ranking analysis.

#### <b>Outputs</b>
The output of the analysis is a trained logistic regression model and an updated anndata object as follows:

- `anndata.obs['data_splits']`: A new column in `anndata.obs` containing labels indicating whether each cell belongs to the training, validation, or test dataset after the feature matrix and target vector are split accordingly.
- `anndata.uns['model_selection_metrics']`: A DataFrame object is stored in `anndata.uns`, containing mean accuracy scores of trained regularized models on the training, validation, and test datasets.
- `anndata.uns['<target_class>_marker_results']`: A dictionary object is stored in `anndata.uns`, containing results of the feature ranking analysis specific to a designated target class.
- `anndata.uns['<target_class>_probEst_Summary']`: A DataFrame object is stored in `anndata.uns`, containing the mean probability estimates for each cluster belonging to the specified target class.
- `anndata.uns['<target_class>_AUROCC_Summary']`: A DataFrame object is stored in `anndata.uns`, containing summary information about the AUROCC (Area Under the Receiver Operating Characteristic Curve) scores for each cluster belonging to the specified target class.

## Citations
------------------------------------------
Worm glia atlas manuscript is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.03.21.533668v1). Please cite our paper if you use these tutorials for your analyses:

```
@article {Purice2023.03.21.533668,
	author = {Maria D. Purice and Elgene J.A. Quitevis and R. Sean Manning and Liza J. Severs and Nina-Tuyen Tran and Violet Sorrentino and Manu Setty and Aakanksha Singhvi},
	title = {Molecular heterogeneity of C. elegans glia across sexes},
	elocation-id = {2023.03.21.533668},
	year = {2023},
	doi = {10.1101/2023.03.21.533668},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2023/03/24/2023.03.21.533668},
	eprint = {https://www.biorxiv.org/content/early/2023/03/24/2023.03.21.533668.full.pdf},
	journal = {bioRxiv}
}

```
