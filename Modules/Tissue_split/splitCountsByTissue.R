library(dplyr)
#install.packages('stringr')
#library(stringr)


split_counts_per_tissue <- function(bulk_counts, design_file){

    # get all tissues from sdrf file that are present in experiment
    tissues = unique(design_file$Factor.Value.organism.part.)
   
    tissue_dict = c()
    for (tissue in tissues){
       
        # get runs that belong to tissue
        tissue_runs = design_file[design_file$Factor.Value.organism.part. == tissue, ]
        tissue_runs = unique(tissue_runs$Run)
        # extract column from count file for each tissue
        bulk_counts_new = bulk_counts[,colnames(bulk_counts) %in%  tissue_runs]
        # save counts as .rds file
        saveRDS(bulk_counts_new,  gsub(" ", "", paste0('Tissue_splits/' ,exp_name, '-',
                tissue , '-fpkms.rds')))
        tissue_dict[tissue] = list(tissue_runs)

    }
   
    #return(tissue_dict)
}


args <- commandArgs(trailingOnly = TRUE)
counts = args[1]
sdrf = args[2]
exp_name = args[3]
output = args[4]


counts = read.csv(counts, sep = "\t",check.names = FALSE )
# remove duplicated columns
counts = counts %>% 
  select(unique(colnames(.)))

rownames(counts) = counts[,1]
counts = counts[,2:ncol(counts)]

print(head(counts))

sdrf = read.csv(sdrf, sep = "\t")
print(sdrf)

split_counts_per_tissue(counts, sdrf)


