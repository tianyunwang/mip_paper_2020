#################################################################
# de novo enrichment analysis - denovolyzeR
#################################################################

# set up working directory
setwd("/path/to/directory")

# load packages (install first if not)
# install.packages("denovolyzeR")
library(denovolyzeR)
# install.packages("tidyverse")
library(tidyverse)

# sample size
denovo_db_NDD_Size <- 10927
SPARK_ASD_Size <- 6499

# read in DNM counts
dnm_count <- read.delim("125genes_10927NDD_6499ASD_denovo_LGD_MIS.txt",as.is=TRUE)
colnames(dnm_count) <- c("gene","LOF","MIS")

# Reformat dataframe to denovolyzer format
## LOF SECTION
dnm_table_LOF <- data.frame(gene="",class="")

for(i in 1:nrow(dnm_count)){
  slice_1 <- rep(dnm_count$gene[i],dnm_count$LOF[i])
  slice_2 <- rep("lof",dnm_count$LOF[i])
  slice_combo <- data.frame(gene=slice_1,class=slice_2)
  dnm_table_LOF <- rbind(slice_combo,dnm_table_LOF)
}

dnm_table_LOF <- (dnm_table_LOF[dnm_table_LOF$gene %in% dnm_count$gene,])


## MIS SECTION
dnm_table_MIS <- data.frame(gene="",class="")

for(i in 1:nrow(dnm_count)){
  slice_1 <- rep(dnm_count$gene[i],dnm_count$MIS[i])
  slice_2 <- rep("mis",dnm_count$MIS[i])
  slice_combo <- data.frame(gene=slice_1,class=slice_2)
  dnm_table_MIS <- rbind(slice_combo,dnm_table_MIS)
}

dnm_table_MIS <- (dnm_table_MIS[dnm_table_MIS$gene %in% dnm_count$gene,])

# merge table lof + mis
dnm_table <- rbind(dnm_table_LOF, dnm_table_MIS)

# Run denovolyzebygene
denovolyzeR_stats <- denovolyzeByGene(genes=dnm_table$gene,
                                          classes=dnm_table$class,
                                          nsamples=(denovo_db_NDD_Size + SPARK_ASD_Size),
                                          includeClasses = c("lof","prot","mis"))

## Adjust P-Values For Multiple Comparisons - Benjamini & Hochberg (1995) ("BH" or its alias "fdr")
## NOTE!! The q-values in paper were corrected for all genes with DNM in the 10927 NDD trios (aka. denovo-db) and 6499 ASD trios (aka. SPARK), 
## in this example code, DNM counts only provided for the 125 genes, so the q-values achieved here may slightly different as what in paper.
denovolyzeR_stats$lof_qValue <- p.adjust(denovolyzeR_stats$lof_pValue, method = "BH", n = 19618)
denovolyzeR_stats$mis_qValue <- p.adjust(denovolyzeR_stats$mis_pValue, method = "BH", n = 19618)
denovolyzeR_stats$prot_qValue <- p.adjust(denovolyzeR_stats$prot_pValue, method = "BH", n = 19618)

# write to output
denovolyzeR_stats <- denovolyzeR_stats %>% select(gene, lof_observed, mis_observed, prot_observed,
                                                  lof_pValue, lof_qValue, mis_pValue, mis_qValue, 
                                                  prot_pValue, prot_qValue)

write.csv(denovolyzeR_stats,"125genes_10927NDD_6499ASD_DNM_denovolyzeR_p_q.csv",row.names = F)
