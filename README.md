MuSiC2
======================

`MuSiC2` is an iterative algorithm aiming to improve cell type deconvolution for bulk RNA-seq data using scRNA-seq data as reference when the bulk data are generated from samples with multiple clinical conditions where at least one condition is different from the scRNA-seq reference. The key idea of MuSiC2 is that, when the bulk samples and single-cell reference samples are from different clinical conditions, the majority of genes shall still share similar cell-type-specific gene expression pattern regardless of clinical conditions. By removing genes with cell-type-specific differential expression (DE) between samples with different clinical conditions from the single-cell reference, MuSiC2 holds the potential to yield more accurate cell type proportion estimates.

<p align="center"> 
<img src="./Figure 1.jpg" width="700">
</p>

How to cite `MuSiC2`
-------------------
Please cite the following publication:

> *MuSiC2: cell type deconvolution for multi-condition bulk RNA-seq data*<br />
> <small>J. Fan, Y. Lyu, Q. Zhang, X. Wang, R. Xiao, M. Li<br /></small>

Installation
------------

``` r
# install devtools if necessary
if (!"devtools" %in% rownames(installed.packages())) {
  install.packages('devtools')
}
# install the MuSiC2 package
if (!"MuSiC2" %in% rownames(installed.packages())) {
  devtools::install_github('Jiaxin-Fan/MuSiC2')
}
# load
library(MuSiC2)
```

More Information
-----------------
Please see [Tutorial](https://jiaxin-fan.github.io/MuSiC2/articles/introduction.html).

