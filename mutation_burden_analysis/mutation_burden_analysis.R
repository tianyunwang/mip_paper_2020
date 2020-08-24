###################################################################
# Mutation burden analyses for ultra-rare LGD and MIS30 variants
###################################################################

# set up working directory
setwd("/path/to/directory")

# read in my_dataset
my_dataset <- read.table("MIP125_genes_ultra_rare_LGD_MIS30_ExAC10X_0.9_noFalsePositive.txt", header = T)

# one-sided Fisher's exact test
my_dataset$Fisher_LGDp <- rep("NA",nrow(my_dataset))
my_dataset$Fisher_MIS30p <- rep("NA",nrow(my_dataset))

for (i in 1:nrow(my_dataset)){
  LGD_Matrix <- matrix(c(my_dataset[i,]$mip_case_lof,my_dataset[i,]$exac_lof,my_dataset[i,]$mip_case - my_dataset[i,]$mip_case_lof, 45376 - my_dataset[i,]$exac_lof),ncol = 2)
  my_dataset[i,]$Fisher_LGDp <- fisher.test(LGD_Matrix,alternative = "greater")$p.value
  
  MIS30_Matrix <- matrix(c(my_dataset[i,]$mip_case_mis30,my_dataset[i,]$exac_mis30,my_dataset[i,]$mip_case - my_dataset[i,]$mip_case_mis30, 45376 - my_dataset[i,]$exac_mis30),ncol = 2)    
  my_dataset[i,]$Fisher_MIS30p <- fisher.test(MIS30_Matrix,alternative = "greater")$p.value
}

my_dataset$Fisher_LGDq <- rep("NA",nrow(my_dataset))
my_dataset$Fisher_MIS30q <- rep("NA",nrow(my_dataset))

# multiple tests correction for the p-value by the number of genes in test (FDR significance).
my_dataset$Fisher_LGDq <- p.adjust(my_dataset$Fisher_LGDp, method = "BH", n = 125)
my_dataset$Fisher_MIS30q <- p.adjust(my_dataset$Fisher_MIS30p, method = "BH", n = 125)

# write output to file
write.csv(my_dataset, file = "MIP125_genes_ultra_rare_LGD_MIS30_ExAC10X_0.9_noFalsePositive_mutation_burden.csv", quote = F)


