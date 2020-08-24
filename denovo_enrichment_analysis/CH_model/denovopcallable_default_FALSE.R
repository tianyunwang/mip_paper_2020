denovop <- function(genes,variants,priors,correct=FALSE,genecount=18946){
  # Function replicating cALTruncp in R with exact binomial probabilities
  # Inputs:
  #   genefile: Filename for input data. format has three coluMISn with headers 
  #             (GeneName  MIS	LGD)
  #   variants: Number of variants in simulation. Either number of denovo variants
  #             in study (exome) or denovo rate * number of samples (targeted)
  #   correct:  OPTIONAL Default FALSE, TRUE for FDR q-values in place of p
  #   genecount: OPTIONAL Total number of genes assayed (used for FDR correction)
  
  # total variants:  Enter Number of Denovo variants in study (i.e. 1.8 x number of samples)
  Var <- variants
  
  LGDp <- rep(NaN,dim(genes)[1])
  MISp <- rep(NaN,dim(genes)[1])
  ALTp <- rep(NaN,dim(genes)[1])
  
  # do the following for each gene
  for (i in 1:dim(genes)[1]){
    genepriors <- priors[priors[,4] == as.character(row.names(genes)[i]),]
    if (dim(genepriors)[1] > 1){
      genepriors <- genepriors[1,]
    }
    if (dim(genepriors)[1] > 0){
      #binom test for MIS
      MISp[i] <- binom.test(as.numeric(genes$MIS[i]), Var, p=as.numeric(genepriors[5]),alternative="greater")$p.value
      #binom test for LGD
      LGDp[i] <- binom.test(as.numeric(genes$LGD[i]), Var, p=as.numeric(genepriors[6]),alternative="greater")$p.value
      #Protein ALTering p
      ALTp[i] <- binom.test(as.numeric(genes$LGD[i])+as.numeric(genes$MIS[i]), Var, p=as.numeric(genepriors[5])+as.numeric(genepriors[6]),alternative="greater")$p.value
    }
  }
  if (correct){
    genes$LGDp <- p.adjust(LGDp,method="BH",n=genecount)
    genes$MISp <- p.adjust(MISp,method="BH",n=genecount)
    genes$ALTp <- p.adjust(ALTp,method="BH",n=genecount)
    
  }else{
    genes$LGDp <- LGDp
    genes$MISp <- MISp
    genes$ALTp <- ALTp
  }
  return(genes)
}