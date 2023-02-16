library(BiocManager)
library(devtools)
library(remotes)

# monocle
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
devtools::install_github('cole-trapnell-lab/monocle3') # @b545460966874948eb11a57a225594a107f1694d

# seurat wrappers
remotes::install_github('satijalab/seurat-wrappers') # @d28512f804d5fe05e6d68900ca9221020d52cf1d