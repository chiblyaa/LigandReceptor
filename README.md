
<!-- README.md is generated from README.Rmd. Please edit that file -->
# LigandReceptor

<!-- badges: start -->
<!-- badges: end -->
The goal of LigandReceptor is to simplify the process of identifying potential ligand-receptor pairs in scRNAseq data

## Installation

You can install LigandReceptor from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("chiblyaa/LigandReceptor")
```

## Example

To use this function, first obtain a list of genes expressed in a scRNAseq dataset using the SEURAT function FindAllMarkers():

``` r
## basic example code
library(Seurat)
library(LigandReceptor)

seuratDEGS <- LigandReceptor::test_dataset #test_dataset obtained using FindAllMarkers from SEURAT
ncells = 23 #number of unique cell identities in SEURAT object
celltypelabels = as.vector(unique(test_dataset$cluster)) # create vector of cell identities
colors = topo.colors(ncells) # create vector of colors of length = ncells
LRdatabase <- LigandReceptor::LRdatabase # Load ligand-receptor pair reference database

## This will generate a table with ligand-receptor pair for the specified cell types in celltypelabels
LR.pairs <- LigandReceptorPairsTable(ncells, celltypelabels, seuratDEGS, LRdatabase)

# This will generate the chord plot associated with the table
#PairsPlot <- function(filename, ncells, celltypelabels, cellcolors, seuratDEGS, LRdatabase, subsetgenes, from = celltypelabels, to = celltypelabels)
```

The generated table looks like this:

``` r
head(LR.pairs, 10)
#>                 from      to value
#> 1            End bud End bud     1
#> 2        Krt19+ duct End bud     7
#> 3         Basal duct End bud     0
#> 4      Myoepithelial End bud     2
#> 5  Bpifa2+ Proacinar End bud     6
#> 6    Smgc+ Proacinar End bud     0
#> 7      Mitotic cells End bud     1
#> 8      Serous acinar End bud     5
#> 9  Seromucous acinar End bud     5
#> 10         Gstt1+ ID End bud     4
#>                                                                                     pairs
#> 1                                                                             Calm2_pde1c
#> 2  Tnc_sdc4, Lamc2_itgb4, Lamb3_itgb4, Lama5_itgb4, Hsp90aa1_cftr, Hbegf_cd9, Calm2_pde1c
#> 3                                                                                        
#> 4                                                              Hsp90b1_asgr1, Calm1_pde1c
#> 5            Tgm2_sdc4, Sema6a_plxna4, Lamc1_itgb4, Lamb1_itgb4, Cxcl12_sdc4, Calm1_pde1c
#> 6                                                                                    <NA>
#> 7                                                                             Calm2_pde1c
#> 8                            Spint1_st14, Hsp90aa1_cftr, Hp_asgr1, Hbegf_cd9, Calm2_pde1c
#> 9                        Lamc1_itgb4, Hsp90b1_asgr1, Hbegf_cd9, Dusp18_itgb4, Calm3_pde1c
#> 10                                      Spint1_st14, Hsp90aa1_cftr, Hp_asgr1, Calm1_pde1c
```

The corresponding plot would look like this:
