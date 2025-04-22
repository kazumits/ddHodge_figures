# Figures of ddHodge 

This repository provides the Jupyter notebooks that reproduce the results in [our preprint of ddHodge](https://www.biorxiv.org/content/10.1101/2025.04.16.649050v1).

The Julia package of ddHodge is available at [ddHodge.jl](https://github.com/kazumits/ddHodge.jl).

## Preface

Most of the code written here was devoted to making plots.

In addition to reproducing our results, these helper tools can also be used for simulation and data analysis using Julia.

* [tools/Dynamical2D.jl] : A Julia module for setting up 2D dynamical systems, such as calculating analytical div, grad and curl in curved space.
* [tools/H5ADHelper.jl] : A Julia module to add a new column to the `obs` DataFrame of the loaded [anndata](https://anndata.readthedocs.io/en/stable/index.html) file.

Please refer to the examples below for specific usage:

* [tools/Dynamical2D.jl] : [notebook/surfsim_fig.ipynb](notebook/surfsim_fig.ipynb)
* [tools/H5ADHelper.jl] : [notebook/vasa_fig.ipynb](notebook/vasa_fig.ipynb)

The basic usage of `H5ADHelper.jl` is also explained in [ddHodge.jl](https://github.com/kazumits/ddHodge.jl).

## Figures

These figures were produced using ddHodge ([version 0.7.0](https://github.com/kazumits/ddHodge.jl/releases/tag/v0.7.0)).

### Figure 1

Showing an example of artificial dynamics in 3D and performing velocity embedding:

* [notebook/draw3D_fig.ipynb](notebook/draw3D_fig.ipynb)

### Figure 2

Validation of ddHodge performance using an artificial 2D gradient system and toggle-switch model:

* [notebook/sim2D_fig.ipynb](notebook/sim2D_fig.ipynb)
* [notebook/dynsim_fig.ipynb](notebook/dynsim_fig.ipynb)

### Figure 3

Validation of ddHodge in higher dimensional dynamics, such as artificial dynamics on a sphere and the Lorenz system:

* [notebook/surfsim_fig.ipynb](notebook/surfsim_fig.ipynb)
* [notebook/sim3D_fig.ipynb](notebook/sim3D_fig.ipynb)

### Figure 4-5

Using ddHodge to analyse scRNA-seq data:

* [notebook/vasa_fig.ipynb](notebook/vasa_fig.ipynb)
* [notebook/fucci_fig.ipynb](notebook/fucci_fig.ipynb)

The RNA velocities were calculated as follows:

* [notebook/scvelo_fucci.ipynb](notebook/scvelo_fucci.ipynb)
* [notebook/scvelo_vasa.ipynb](notebook/scvelo_vasa.ipynb)

The RNA velocity was calculated using scvelo (version 0.3.2).

## Data availability

For Figures 4-5, some of the public data are needed to reproduce our results. The small size (< 20MB) of the pre-processed data has been uploaded to the [data](data/) directory.

The datasets not included in this repository are availabel at:

* Spliced/unspliced read count of VASA-seq: [GSE176588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176588)
* List of histone genes (HistoneGenes.tsv): [VASAseq_2022](https://github.com/hemberg-lab/VASAseq_2022)
* FUCCI dataset (RNAData/a.loom) in: [input.zip](https://drive.google.com/file/d/1yYvkZVXz2EnoWBrLYi5_hu4WVlnmuFWd/view) (Note: the file was renamed to `fucci_RNAData.loom` in our analysis)

For Figures 1-3, the notebooks are self-contained. No additional files are required.
