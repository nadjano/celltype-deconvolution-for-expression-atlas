library('Seurat')
library('plyr')
library('dplyr')
library('tools') 
library('tidyr')
#library('Biobase')
library('tidyverse')
library(stringr)


get_300_cells_per_celltype <- function(seurat){

    #Split cells by type
    cellSplits <- SplitObject(seurat, split.by = "cellType")
    cellSplits_red = cellSplits

    for (i in 1:length(levels(seurat$cellType))){
        if (length(rownames(cellSplits[[i]])) < 300){
            n = length(rownames(cellSplits[[i]]))
            }
        else {
            n = 300
        }
        cell.list <- WhichCells(cellSplits[[i]], downsample = n)
        cellSplits_red[i] <- seurat[, cell.list]
        }  

    seurat = merge(cellSplits_red[[1]], y = cellSplits_red[2:length(levels(seurat$cellType))] )    
    return(seurat)
}

remove_rare_celltypes = function(seurat){
  row_count <- as.data.frame(x = table(seurat$cellType))
  
  non_rare = row_count[row_count$Freq >100 ,]
  
  #sc$cellType_freq = mapvalues(df_cellType$sc.cellType, from = row_count$Var1, to = row_count$Freq)
  
  seurat <- subset(seurat, subset = cellType %in% paste(non_rare$Var1))
  return(seurat)
  

}
split_sc_per_tissue <- function(seurat, tissue, output){
    
    
    if (!file.exists(output)) {
        seurat = readRDS(seurat)

        cellSplits <- SplitObject(seurat, split.by = "tissue")
                
            
        for (i in 1:length(levels(seurat$tissue))){
        
        
            sc = cellSplits[[i]]
            sc$tissue = droplevels(x = sc$tissue)
            sc$cellType = droplevels(x = sc$cellType)
            this_tissue = levels(unique((sc$tissue)))
            if (this_tissue ==  tissue){
                sc = get_300_cells_per_celltype(sc)
                #sc = equal_celltype_per_sample(sc)
                sc = remove_rare_celltypes(sc)
                saveRDS(sc, output)

            }
         
    }

    } else {
                print(paste(output, "already exists. Skipping..."))
            }
}



args <- commandArgs(trailingOnly = TRUE)
tissue = args[1]
sc_ref = args[2]
output = args[3]


split_sc_per_tissue(sc_ref, tissue, output)


   