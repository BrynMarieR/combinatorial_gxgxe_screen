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

data = readRDS("./processed_data/clean.b12.rds")
genelist <- as.character(c(read.csv("./raw_data/genelist.csv",stringsAsFactors = F,
                                    header=T))[[1]])

plot_specific_dko_one_way <- function(data, geneA, geneB, condition="Glc15",
                                      return.plot =F, ylims=c(-3.5,2)) {
  sko.subset <- data[(data$U6.negctrl + data$H1.negctrl == 1),]
  sko.subset$promoter <- c("H1","U6")[as.factor(
    sko.subset$H1.negctrl == 1)]
  sko.subset$gene <- sko.subset$H1.gene
  sko.subset[which(sko.subset$promoter == "U6"),"gene"] <- 
    sko.subset[which(sko.subset$promoter == "U6"),"U6.gene"]
  
  ctrl.subset <- data[(data$U6.negctrl + data$H1.negctrl == 2),]
  ctrl.subset$gene <- "ctrl"
  ctrl.subset$promoter <- "ctrl"
  
  sko.A <- sko.subset[which(sko.subset$gene == geneA),]
  sko.B <- sko.subset[which(sko.subset$gene == geneB),]
  
  dko.subset.ab <- data[c(which(data$U6_H1 == paste0(geneA,"_",geneB))),]
  dko.subset.ba <- data[c(which(data$U6_H1 == paste0(geneB,"_",geneA))),]
  dko.subset <- rbind(dko.subset.ab, dko.subset.ba)
  
  nrows.df <- nrow(ctrl.subset) + nrow(sko.A) + nrow(sko.B) + nrow(dko.subset)
  mega.df <- data.frame(matrix(rep("", 2*nrows.df),nrow=nrows.df),
                        stringsAsFactors=F)
  
  mega.df[1:nrow(ctrl.subset),1] <- "ctrl_ctrl"
  mega.df[1:nrow(ctrl.subset),2] <- ctrl.subset[,condition]
  index <- nrow(ctrl.subset) + 1
  
  mega.df[index:(index + nrow(sko.A) - 1),1] <- paste0("sko_",geneA)
  mega.df[index:(index + nrow(sko.A) - 1),2] <- sko.A[,condition]
  index <- index + nrow(sko.A)
  
  mega.df[index:(index + nrow(sko.B) - 1),1] <- paste0("sko_",geneB)
  mega.df[index:(index + nrow(sko.B) - 1),2] <- sko.B[,condition]
  index <- index + nrow(sko.B)
  
  mega.df[index:(index + nrow(dko.subset) - 1),1] <- paste0("dko_",geneA,"_",geneB)
  mega.df[index:(index + nrow(dko.subset) - 1),2] <- dko.subset[,condition]
  
  if(geneA == geneB) {
    mega.df <- mega.df[!duplicated(mega.df),]
  }
  
  plot.df <- data.frame(
    ko.status = mega.df[,1],
    logFC = mega.df[,2],
    simple.ko.status = as.factor(substr(mega.df$X1, 1, 3))
  )
  plot.df$logFC <- as.numeric(as.character(plot.df$logFC))
  
  plot.df$ko.status <- factor(plot.df$ko.status, ordered=T,
                              levels=c("ctrl_ctrl",
                                       paste0("sko_",geneA),
                                       paste0("sko_",geneB),
                                       paste0("dko_",geneA,"_",geneB)
                              ))
  
  calc.exp.df <- aggregate(logFC ~ ko.status, data=plot.df, median)
  expected <- sum(calc.exp.df[1:3,2])
  
  g <- ggplot(plot.df, aes(x=ko.status,y=logFC,fill=ko.status)) + 
    geom_violin() + theme_classic() + 
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.5, aes(color=simple.ko.status)) +
    ylim(ylims) + 
    ylab("LFC") + scale_x_discrete(labels=c("ctrl", paste0(geneA," KO"), 
                                            paste0(geneB, " KO"), "DKO")) + 
    theme(legend.position="none",axis.text.x=element_text(size=12,angle=45,vjust=0.5)) +
    scale_fill_manual(values=c("#999999","#ffa0b1","#00BFC4","#C77CFF")) + 
    scale_color_manual(values=rep("black",4)) + 
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.5) +
    geom_segment(color="#E7570F", size=1, aes(x = 3.8, y = expected, xend=4.2, yend=expected))
  if(return.plot) {
    return(g + ggtitle(condition) + xlab(""))
  } else {
    print(g + ggtitle(paste0(geneA," + ",geneB)) + xlab(condition))
  }
}

# SLC25A36+SLC25A1
# SLC25A36+SLC25A19
# SLC25A36+SLC25A20
# SLC25A36+SLC25A25
# SLC25A36+SLC25A3
# SLC25A36+SLC25A35

## gxgxe violin plots by condition
gene_pairs <- matrix(c(
                       "SLC25A33","SLC25A36",
                       "SLC25A1","SLC25A36",
                       # "SLC25A3","SLC25A36",
                       # "SLC25A19","SLC25A36",
                       # "SLC25A20","SLC25A36",
                       # "SLC25A25","SLC25A36",
                       # "SLC25A35","SLC25A36",
                       "SLC25A5","SLC25A6",
                       "MTCH1","MTCH2",
                       "SLC25A16","SLC25A42",
                       "SLC25A28","SLC25A37",
                       # "SLC25A13","SLC25A45",
                       # "SLC25A13","SLC25A33",
                       "SLC25A37","SLC25A39",
                       # "SLC25A11","SLC25A31",
                       # "SLC25A14","SLC25A30",
                       # "MTCH2","SLC25A1",
                       # "SLC25A43","SLC25A1",
                       # "MTCH2","SLC25A2",
                       "MTCH2","SLC25A36",                       
                       "MTCH2","SLC25A20",
                       "SLC25A11","SLC25A19",
                       # "SLC25A32","SLC25A26",
                       # "MTCH2", "SLC25A19"
                       "BCL2L1","MCL1",
                       "UCP1","SLC25A31"
                       ), ncol=2,byrow=T)
pdf("./figs/violin_gxg_by_condition.pdf", 
    # family="Helvetica", this is the default
    height=3.8,width=6.5)
for(i in 1:nrow(gene_pairs)) {
  gene_pair = gene_pairs[i,]
  gene_pair_data_subset <- data[which(data$U6_H1 == paste0(gene_pair[1],"_", gene_pair[2]) | 
                                        data$U6_H1 == paste0(gene_pair[2],"_", gene_pair[1]) | 
                                        (data$U6.gene == gene_pair[1] & data$H1.negctrl==T) |
                                        (data$U6.gene == gene_pair[2] & data$H1.negctrl==T) |
                                        (data$H1.gene == gene_pair[1] & data$U6.negctrl==T) |
                                        (data$H1.gene == gene_pair[2] & data$U6.negctrl==T) | 
                                        (data$H1.negctrl == T & data$U6.negctrl == T) ),]
  ylims = c( min(gene_pair_data_subset[,c("Glc15","Gal15","Anti15","Pyr15")]) - 0.15, 
             max(gene_pair_data_subset[,c("Glc15","Gal15","Anti15","Pyr15")]) + 0.15)
  pglc <- plot_specific_dko_one_way(data, gene_pair[1], gene_pair[2], "Glc15",
                                    return.plot=T, ylims)
  pglc <- pglc + theme(aspect.ratio = 3/2)
  pgal <- plot_specific_dko_one_way(data, gene_pair[1], gene_pair[2], "Gal15",
                                    return.plot=T, ylims)
  pgal <- pgal + theme(axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                       aspect.ratio = 3/2)
  panti <- plot_specific_dko_one_way(data, gene_pair[1], gene_pair[2], "Anti15",
                                     return.plot=T, ylims)
  panti <- panti + theme(axis.title.y = element_blank(),
                       axis.text.y = element_blank(),
                       aspect.ratio = 3/2)
  ppyr <- plot_specific_dko_one_way(data, gene_pair[1], gene_pair[2], "Pyr15",
                                    return.plot=T, ylims)
  ppyr <- ppyr + theme(axis.title.y = element_blank(),
                       axis.text.y = element_blank(),
                       aspect.ratio = 3/2)

  grid.arrange(cbind(ggplotGrob(pglc), ggplotGrob(pgal), ggplotGrob(panti), ggplotGrob(ppyr), size = "last"),
               top=textGrob(paste0(gene_pairs[i,1]," + ",gene_pairs[i,2]),
                            gp = gpar(fontsize = 24)),
               nrow=1)
}
dev.off()





## paralog pi scores
paralog_pairs <- matrix(c(
  "SLC25A5","SLC25A6",
  "MTCH1","MTCH2",
  "SLC25A28","SLC25A37",
  "SLC25A16","SLC25A42",
  "SLC25A37","SLC25A39",
  "SLC25A11","SLC25A19",  
  "MTCH2","SLC25A20",
  "SLC25A1","SLC25A36",
  "MTCH2","SLC25A36", 
  "BCL2L1","MCL1",
  "UCP1","SLC25A31"
), ncol=2,byrow=T)
paralog_data <- data.frame(matrix(rep(0, nrow(paralog_pairs)*4), nrow=nrow(paralog_pairs)))
colnames(paralog_data) <- c("Glc15", "Gal15", "Anti15", "Pyr15")
rownames(paralog_data) <- apply(paralog_pairs, 1, paste0,collapse="+")
paralog_data_noz <- paralog_data
for(condition in c("Glc15", "Gal15", "Anti15", "Pyr15")) {
  final.gi.df <- readRDS(paste0("./processed_data/final_gi_df_",condition,".rds"))
  
  # calculate z scores of pi scores
  final.gi.df.z <- (final.gi.df[,,1] - mean(final.gi.df[,,1])) / sd(final.gi.df[,,1])
  # ranks
  # by rows
  #final.gi.df.z <- matrix(rank(final.gi.df[1:63,1:63,1]), nrow=63,byrow=T)
  #rownames(final.gi.df.z) <- rownames(final.gi.df[,,])
  #colnames(final.gi.df.z) <- colnames(final.gi.df[,,])
  paralog_data[,condition] <- apply(paralog_pairs, 1, function(s) {
    final.gi.df.z[s[1],s[2]]
  })
  paralog_data_noz[,condition] <- apply(paralog_pairs, 1, function(s) {
    final.gi.df[s[1],s[2],1]
  })
}
pdf("./figs/paralog_heatmap.pdf")
pheatmap(paralog_data,breaks=c(-12.5,seq(-5,5,l=98), 5.00001),
         col=colorRampPalette(c("cyan", "black", "yellow"), 
                              interpolate="spline")(n = 100),
         cluster_rows=F,cluster_cols=F,cellwidth=15,cellheight=15)
dev.off()

write.xlsx(paralog_data, "paralog_data_zscores.xlsx",row.names=T)
write.xlsx(paralog_data_noz, "paralog_data_raw.xlsx",row.names=T)
