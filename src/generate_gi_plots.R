# Author: Bryn Reinstadler
# Date: 2016 - 2017
# Purpose: To prepare visualizations of genetic interaction data for 
#    the genes of interest in a dual-knock-out combinatorial CRISPR screen.


# Libraries & Packages
library(ape)
library(ggplot2)
library(ggrepel)
library(beeswarm)
library(limma)
library(reshape)
library(openxlsx)
library(pheatmap)
library(gtools)
library(grid)
library(gridExtra)
library(RColorBrewer)

source("./calculate_pi_score.R")

data = readRDS("./processed_data/clean.b12.rds")
genelist <- as.character(c(read.csv("./raw_data/genelist.csv",stringsAsFactors = F,
                                    header=T))[[1]])

## big heatmap of heatmaps
# data generated by CLUSTAL
dat <- ape::read.tree(text="(((((((((((SLC25A5:0.04270,SLC25A6:0.04270):0.02278,SLC25A4:0.06548)
                      :0.12499,SLC25A31:0.19046):0.82875,((((SLC25A23:
                      0.26773,SLC25A24:0.26773):0.02515,SLC25A25:0.29288):0.16066,SLC25A41:0.45354):0.36086,
                      (SLC25A42:0.64008,SLC25A16:0.64008):0.17431):0.20482):0.03447,SLC25A43:1.05368):0.09449,SLC25A19:1.14817):
                      0.07886,((((SLC25A10:0.59087,SLC25A11:0.59087):0.19913,(((UCP2:0.16588,UCP3:0.16588):
                      0.14387,UCP1:0.30974):0.39190,((SLC25A30:0.11662,SLC25A14:0.11662):0.51504,SLC25A27:0.63167):
                      0.06998):0.08836):0.11386,(SLC25A34:0.40355,SLC25A35:0.40355):0.50031):0.30018,(((((((SLC25A33:
                      0.27244,SLC25A36:0.27244):0.61790,(SLC25A28:0.22486,SLC25A37:0.22486):
                      0.66548):0.11140,SLC25A32:1.00174):0.09422,SLC25A17:1.09596):0.02526,((SLC25A39:
                      0.40392,SLC25A40:0.40392):0.60015,SLC25A38:1.00407):0.11715):0.02007,SLC25A44:1.14129):
                      0.00434,((((((SLC25A13:0.13750,SLC25A12:0.13750):0.43722,(SLC25A22:0.24465,SLC25A18:0.24465):
                      0.33007):0.32447,SLC25A21:0.89920):0.11377,SLC25A26:1.01297):0.04462,SLC25A1:1.05759):
                      0.03918,((SLC25A15:0.07205,SLC25A2:0.07205):0.86324,((((SLC25A45:0.47247,SLC25A48:
                      0.47247):0.06614,SLC25A47:0.53861):0.07890,SLC25A29:0.61750):0.14624,SLC25A20:0.76374):
                      0.17155):0.16148):0.04886):0.05842):0.02299):0.37789,SLC25A3:1.60492):0.02483,((SLC25A51:
                      0.02361,SLC25A52:0.02361):0.89682,SLC25A53:0.92043):0.70932):0.51056,(MTCH1:
                      0.43482,MTCH2:0.43482):1.70549):0.12256,SLC25A46:2.26287);")
dend <- chronos(dat)
plot(dend)
hc <- as.hclust(dend)
slc.dendro <- as.dendrogram(hc)

#plot the dendrogram
pdf("./figs/dendrogram.pdf",height=4,width=11)
plot(slc.dendro)
dev.off()

my_palette <- colorRampPalette(c("cyan", "black", "yellow"), 
                               interpolate="spline")(n = 100)

lfc.data <- data

sko.lfcs <- data[xor(data$U6.negctrl == T, data$H1.negctrl == T),]
sko.lfcs$gene <- sko.lfcs$U6.gene
sko.lfcs$gene[which(sko.lfcs$U6.negctrl == T)] <- sko.lfcs$H1.gene[which(sko.lfcs$U6.negctrl == T)]

for(condition in c("Glc15", "Gal15", "Anti15", "Pyr15")) {
  # final.gi.df <- pi_score(lfc.data, condition,
  #                         genelist,
  #                         sig.test = "nonparametric",
  #                         symmetrize=T)
  # saveRDS(final.gi.df,paste0("./processed_data/final_gi_df_",condition,".rds"))
  final.gi.df <- readRDS(paste0("./processed_data/final_gi_df_",condition,".rds"))
  
  wrong.names.piscores <- final.gi.df[c(dat$tip.label,"SLC16A1","SLC16A7","SLC16A11","SLC22A4","SLC30A6",
                                        "SLC30A9","SLC37A4","MCL1","BCL2L1","EEF2"), 
                                      c(dat$tip.label,"SLC16A1","SLC16A7","SLC16A11","SLC22A4","SLC30A6",
                                        "SLC30A9","SLC37A4","MCL1","BCL2L1","EEF2"),1] 
  # wrong.names.piscores <- (wrong.names.piscores + t(wrong.names.piscores))/2
  diag(wrong.names.piscores) <- NA
  
  df.barplot <- aggregate(. ~ gene, data=sko.lfcs[,c("gene",condition)], FUN=median)
  rownames(df.barplot) <- df.barplot$gene
  colnames(df.barplot) <- c("gene", "avg.lfc")
  df.barplot <- df.barplot[c(dat$tip.label,"SLC16A1","SLC16A7","SLC16A11","SLC22A4","SLC30A6",
                             "SLC30A9","SLC37A4","MCL1","BCL2L1","EEF2"),]
  highest.abs <- max(abs(df.barplot[,2]))
  df.barplot$SKO_LFC <- cut(df.barplot[,2],breaks=seq(-highest.abs, highest.abs,l=11),
                            include.lowest=T)
  annotation_row=data.frame(SKO_LogFC = df.barplot[,"SKO_LFC",drop=F])
  SKO_LFC <- colorRampPalette(c("red","white","blue"),interpolate="spline")(10)
  names(SKO_LFC) <- rev(levels(df.barplot$SKO_LFC))
  barplot_colors_list = list(SKO_LFC = SKO_LFC)
  
  pdf(paste0("./figs/heatmap_piscore_",condition,".pdf"),
      width=11, height=11)
  pheatmap(wrong.names.piscores, cluster_rows=F, cluster_cols=F,
           breaks=c(seq(min(-0.301,c(wrong.names.piscores),na.rm=T), -0.002, length.out=49), 
                    seq(-0.0019999, 0.0019999, length.out=50),
                    seq(0.002, max(c(0.301,wrong.names.piscores),na.rm=T), length.out=49)
           ),
           col=c(my_palette[1:50], rep("black", 50), my_palette[51:100]),
           main=paste0(condition," B1+2"),
           annotation_row=annotation_row,
           annotation_colors=barplot_colors_list,
           cellwidth=9, cellheight=9,na_col="grey25",
           border_color = NA)
  dev.off()
}
