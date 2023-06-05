suppressPackageStartupMessages({
    suppressWarnings({
        library(hdWGCNA, quietly = T)
        library(WGCNA, quietly = T)   
    })
})

set.seed(12345)
WGCNA::enableWGCNAThreads(nThreads = 30)

get.metacells<-function(so, ...){

    so <- hdWGCNA::SetupForWGCNA(
      seurat_obj = so,
      wgcna_name = "wgcna.metacells"
    )

    so<-hdWGCNA::MetacellsByGroups(
      seurat_obj = so,
      reduction = "pca",
      assay = 'RNA',
      slot = 'counts',
      wgcna_name = 'wgcna.metacells', 
      ...
    )

    so<-NormalizeMetacells(so)
    so<-so@misc$wgcna.metacells$wgcna_metacell_obj
    
    return(so)
}