######################################################################
# de novo enrichment analysis - CH model
# CH model was firstly described in O'Roak et al 2012 (PMID:23160955),
# and further developed and applied in Wang et al 2016 (PMID:27824329) 
# and Coe et al 2019 (PMID:30559488). Refer those for more details.
######################################################################

# set up working directory
setwd("/path/to/directory")

# load packages; install first if not
# install.packages("tidyverse")
library(tidyverse)

# load denovop function # FDR correction default as FALSE
source("denovopcallable_default_FALSE.R")

# The CH model gene-level null p-values for LGD and all MISsense
priors <- read.table("CH_Model_AutUS_R_SNPEFF_null_pvalue_LGD_MIS.txt")

# The CH model gene-level null p-values for LGD and MIS30 only
caddpriors <- read.table("CH_Model_AutUS_R_SNPEFF_null_pvalue_LGD_MIS30.txt")

# sample size
denovo_db_NDD_Size <- 10927
SPARK_ASD_Size <- 6499

# run chmodel - LGD + MIS
LGD_MIS <- read.table("125genes_10927NDD_6499ASD_denovo_LGD_MIS.txt", row.names= 1)
variants <- round((denovo_db_NDD_Size + SPARK_ASD_Size) * 1.8, 0)
CH_model_LGD_MIS <- denovop(genes=LGD_MIS, variants=variants, priors=priors)

# run chmodel - LGD + MIS30
LGD_MIS30 <- read.table("125genes_10927NDD_6499ASD_denovo_LGD_MIS30.txt", row.names= 1)
variants <- round((denovo_db_NDD_Size + SPARK_ASD_Size) * 1.8, 0)
CH_model_LGD_MIS30 <- denovop(genes=LGD_MIS30, variants=variants, priors=caddpriors)

# MERGE P VALUES: LGD, MIS, MIS30, ALT
CH_model_LGD_MIS$gene <- row.names(CH_model_LGD_MIS)
CH_model_LGD_MIS30$gene <- row.names(CH_model_LGD_MIS30)

CH_model_combine <- merge(CH_model_LGD_MIS, CH_model_LGD_MIS30,by="gene", all=TRUE)
colnames(CH_model_combine) <- c("gene","dnLGD","dnMIS","dnLGD_pValue","dnMIS_pValue","dnALT_pValue",
                                "dnLGD.2","dnMIS30","dnLGD_pValue.2","dnMIS30_pValue","dnLGD_dnMIS30_pValue")
CH_model_combine_output <- CH_model_combine[,c(1:6,8,10)]

## Adjust P-Values For Multiple Comparisons - Benjamini & Hochberg (1995) ("BH" or its alias "fdr")
## NOTE!! The q-values in paper were corrected for all genes with DNM in the 10927 NDD trios (aka. denovo-db) and 6499 ASD trios (aka. SPARK), 
## in this example code, DNM counts only provided for the 125 genes, so the q-values achieved here may slightly different as what in paper.
CH_model_combine_output$dnLGD_qValue <- p.adjust(CH_model_combine_output$dnLGD_pValue,method = "BH", n = 18946)
CH_model_combine_output$dnMIS_qValue <- p.adjust(CH_model_combine_output$dnMIS_pValue,method = "BH", n = 18946)
CH_model_combine_output$dnMIS30_qValue <- p.adjust(CH_model_combine_output$dnMIS30_pValue,method = "BH", n = 18946)
CH_model_combine_output$dnALT_qValue <- p.adjust(CH_model_combine_output$dnALT_pValue,method = "BH", n = 18946)

# write to output
CH_model_combine_output <- CH_model_combine_output %>% select(gene, dnLGD, dnMIS, dnMIS30, 
                                                              dnLGD_pValue, dnLGD_qValue, dnMIS_pValue, dnMIS_qValue, 
                                                              dnMIS30_pValue, dnMIS30_qValue, dnALT_pValue, dnALT_qValue)

write.csv(CH_model_combine_output,"125genes_10927NDD_6499ASD_DNM_CH_model_p_q.csv",row.names = F)

