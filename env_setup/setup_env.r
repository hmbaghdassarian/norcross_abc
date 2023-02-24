library(BiocManager)
library(devtools)
library(remotes)

# monocle
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'HDF5Array',
#                        'terra', 'ggrastr'), 
#                      version = c('0.44.0', '0.24.0', '1.20.0', 
#                                 '3.54.1', '1.1-31', '0.36.0', '1.20.0', 
#                                 '1.28.0', '1.14.1', '1.26.0', 
#                                 '1.5-21', '1.0.1'),
#                      update = F)

# packages = c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'HDF5Array',
#                        'terra', 'ggrastr')
# versions = c('0.44.0', '0.24.0', '1.20.0', 
#                                 '3.54.1', '1.1-31', '0.36.0', '1.20.0', 
#                                 '1.28.0', '1.14.1', '1.26.0', 
#                                 '1.5-21', '1.0.1')
# for (i in seq_along(packages)){
#     package<-packages[[i]]
#     version<-versions[[i]]
#     BiocManager::install(package, version=version, update = F)
#  }

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'), update = F)

devtools::install_github('cole-trapnell-lab/monocle3@b545460966874948eb11a57a225594a107f1694d', upgrade=F) 

# seurat wrappers
remotes::install_github('satijalab/seurat-wrappers@d28512f804d5fe05e6d68900ca9221020d52cf1d', upgrade=F) 

# liana
remotes::install_github('saezlab/liana@0167d373d428403940df27ffa389977b755eec8a', upgrade=F)