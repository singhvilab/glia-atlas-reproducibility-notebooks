# worm-glia-atlas-notebooks
This is a repository containing figure reproducibility notebooks for the worm glia scRNA-seq atlas available at https://wormglia.org.

## Environment Setup
------------------------------------------
To run the notebooks, clone this repository and follow the environment set up below.
```
git clone git@github.com:eL-Gene/worm-glia-atlas-notebooks.git
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
