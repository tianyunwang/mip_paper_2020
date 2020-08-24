###################################################################
# plot Figure 2a,2c
###################################################################

# set up work directory
setwd("/path/to/directory")

# read in my_dataset
my_dataset <- read.delim(file = "MIP125_genes_log10qvalue.txt", header=T)

# plot Figure 2a - mutation burden analysis
plot(as.numeric(as.character(my_dataset$fisher_LGD)),as.numeric(as.character(my_dataset$fisher_MIS30)),
     xlab="LGD (-log10(q))", ylab="MIS30 (-log10(q))", col="black",pch=20, main ="mutation burden analysis")
segments(x0 = -log10(0.05), y0 = -10, x1 = -log10(0.05), y1 = 100, col = "black", lty = 2)
segments(x0 = -10, y0 = -log10(0.05), x1 = 100, y1 = -log10(0.05), col = "black", lty = 2)

MIS30 <- my_dataset[my_dataset$fisher_MIS30 > (-log10(0.05)) & my_dataset$fisher_LGD < (-log10(0.05)),]
LGD <- my_dataset[my_dataset$fisher_LGD > (-log10(0.05)) & my_dataset$fisher_MIS30 < (-log10(0.05)),]
ALL <- my_dataset[my_dataset$fisher_MIS30 < (-log10(0.05)) & my_dataset$fisher_LGD < (-log10(0.05)),]

points(as.numeric(as.character(MIS30$fisher_LGD)),as.numeric(as.character(MIS30$fisher_MIS30)), col="blue",pch=20)
points(as.numeric(as.character(LGD$fisher_LGD)),as.numeric(as.character(LGD$fisher_MIS30)), col="red",pch=20)
points(as.numeric(as.character(ALL$fisher_LGD)),as.numeric(as.character(ALL$fisher_MIS30)), col="grey",pch=20)

# add gene name - optional
text(as.numeric(as.character(MIS30$fisher_MIS30))~as.numeric(as.character(MIS30$fisher_LGD)),
     labels=MIS30$Label,data=MIS30, cex=1.5, font=3,col="blue", pos=4)

text(as.numeric(as.character(LGD$fisher_MIS30))~as.numeric(as.character(LGD$fisher_LGD)),
     labels=LGD$Label,data=LGD, cex=1.5, font=3,col="red", pos=4)

text(as.numeric(as.character(ALL$fisher_MIS30))~as.numeric(as.character(ALL$fisher_LGD)),
     labels=ALL$Label,data=ALL, cex=1.5, font=3,col="black", pos=2)


# plot Figure 2c - all genes
plot(as.numeric(as.character(my_dataset$ch_dnLGD)),as.numeric(as.character(my_dataset$dy_dnLGD)),
     xlab="CH model (-log10(q))", ylab="denovolyzeR (-log10(q))", col="red",pch=20, main ="de novo enrichment analysis")
points(as.numeric(as.character(my_dataset$ch_dnMIS)),as.numeric(as.character(my_dataset$dy_dnMIS)), col="blue",pch=20)
points(as.numeric(as.character(my_dataset$ch_dnALT)),as.numeric(as.character(my_dataset$dy_dnALT)), col="black",pch=20)

segments(x0 = -log10(0.05), y0 = -10, x1 = -log10(0.05), y1 = 100, col = "black", lty = 2)
segments(x0 = -10, y0 = -log10(0.05), x1 = 100, y1 = -log10(0.05), col = "black", lty = 2)

segments(x0 = -5, y0 = 15, x1 = 15, y1 = 15, col = "black", lty = 3)
segments(x0 = 15, y0 = -5, x1 = 15, y1 = 15, col = "black", lty = 3)

# add gene name - optional
dnMIS <- my_dataset[my_dataset$ch_dnMIS>15 | my_dataset$dy_dnMIS>15,]
text(as.numeric(as.character(dnMIS$dy_dnMIS))~as.numeric(as.character(dnMIS$ch_dnMIS)),
     labels=dnMIS$Gene,data=dnMIS, cex=1, font=3,col="blue", pos=4)
dnLGD <- my_dataset[my_dataset$ch_dnLGD>15 | my_dataset$dy_dnLGD>15,]
text(as.numeric(as.character(dnLGD$dy_dnLGD))~as.numeric(as.character(dnLGD$ch_dnLGD)),
     labels=dnLGD$Gene,data=dnLGD, cex=1, font=3,col="red", pos=2)
dnALT <- my_dataset[my_dataset$ch_dnALT>15 | my_dataset$dy_dnALT>15,]
text(as.numeric(as.character(dnALT$dy_dnALT))~as.numeric(as.character(dnALT$ch_dnALT)),
     labels=dnALT$Gene,data=dnALT, cex=1, font=3,col="black", pos=4)

# plot Figure 2c - subset
dnMIS2 <- my_dataset[my_dataset$ch_dnMIS<15 & my_dataset$dy_dnMIS<15,]
dnLGD2 <- my_dataset[my_dataset$ch_dnLGD<15 & my_dataset$dy_dnLGD<15,]
dnALT2 <- my_dataset[my_dataset$ch_dnALT<15 & my_dataset$dy_dnALT<15,]

dnMIS3 <- dnMIS2[dnMIS2$ch_dnMIS > (-log10(0.05)) | dnMIS2$dy_dnMIS > (-log10(0.05)),]
dnLGD3 <- dnLGD2[dnLGD2$ch_dnLGD > (-log10(0.05)) | dnLGD2$dy_dnLGD > (-log10(0.05)),]
dnALT3 <- dnALT2[dnALT2$ch_dnALT > (-log10(0.05)) | dnALT2$dy_dnALT > (-log10(0.05)),]


plot(as.numeric(as.character(dnLGD2$ch_dnLGD)),as.numeric(as.character(dnLGD2$dy_dnLGD)),
     xlab="CH model (-log10(q))", ylab="denovolyzeR (-log10(q))", col="red",pch=20, main ="de novo enrichment analysis")
points(as.numeric(as.character(dnMIS2$ch_dnMIS)),as.numeric(as.character(dnMIS2$dy_dnMIS)), col="blue",pch=20)
points(as.numeric(as.character(dnALT2$ch_dnALT)),as.numeric(as.character(dnALT2$dy_dnALT)), col="black",pch=20)

segments(x0 = -log10(0.05), y0 = -10, x1 = -log10(0.05), y1 = 100, col = "black", lty = 2)
segments(x0 = -10, y0 = -log10(0.05), x1 = 100, y1 = -log10(0.05), col = "black", lty = 2)

# add gene name - optional
text(as.numeric(as.character(dnMIS3$dy_dnMIS))~as.numeric(as.character(dnMIS3$ch_dnMIS)),
     labels=dnMIS3$Gene,data=dnMIS3, cex=0.7, font=3,col="blue", pos=4)
text(as.numeric(as.character(dnLGD3$dy_dnLGD))~as.numeric(as.character(dnLGD3$ch_dnLGD)),
     labels=dnLGD3$Gene,data=dnLGD3, cex=0.7, font=3,col="red", pos=3)
text(as.numeric(as.character(dnALT3$dy_dnALT))~as.numeric(as.character(dnALT3$ch_dnALT)),
     labels=dnALT3$Gene,data=dnALT3, cex=0.7,  font=3, col="black", pos=2)

