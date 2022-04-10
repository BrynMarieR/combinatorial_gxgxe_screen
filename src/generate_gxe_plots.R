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
## heatmap for single gxe

df <- data.frame()
par(mfrow=c(1,3))
tmp.data <- data
for(condition in c("Glc15", "Gal15","Anti15","Pyr15")) {
  data.name = "B12"
  tmp.data$sub <- tmp.data[,condition] - tmp.data[,"pDNA"]
  hist(tmp.data$sub,main=paste0(data.name,": ", condition))
  tmp.data <- tmp.data[which(tmp.data$U6.negctrl + tmp.data$H1.negctrl >= 1),]
  data.sub.agg.u6 <- aggregate(sub ~ U6.gene, mean, data=tmp.data)
  data.sub.agg.h1 <- aggregate(sub ~ H1.gene, mean, data=tmp.data)
  
  tmp.df <- data.frame(matrix(ncol = 0, 
                              nrow = nrow(tmp.data)*2))
  tmp.df$gene <- c(tmp.data$U6.gene, tmp.data$H1.gene)
  tmp.df$condition <- rep(condition, nrow(tmp.data)*2)
  tmp.df$batch <- rep(data.name, nrow(tmp.data)*2)
  tmp.df$promoter <- c(rep("U6.gene", nrow(tmp.data)),
                       rep("H1.gene", nrow(tmp.data)))
  tmp.df$zscores <- rep(tmp.data$sub - mean(tmp.data[which(tmp.data$U6.negctrl + tmp.data$H1.negctrl == 2),"sub"]) / 
                          sd(tmp.data[which(tmp.data$U6.negctrl + tmp.data$H1.negctrl == 2),"sub"]),2)
  tmp.df$sub <- c(rank(tmp.data$sub), rank(tmp.data$sub))
  
  
  df <- rbind(df, tmp.df)
}

#df <- df[!grepl("OR", df$gene),]
df$condition <- factor(df$condition, ordered=T,
                       levels=c("Glc15","Gal15","Pyr15","Anti15"))
df.median <- aggregate(zscores ~ gene + condition, data=df, FUN = median)
df.sd <- aggregate(zscores ~ gene + condition, data=df, FUN = sd)
df.median.by.cond <- data.frame(gene = df.median$gene[1:length(unique(df.median$gene))],
                              Glucose = df.median[which(df.median$condition == "Glc15"),"zscores"],
                              Galactose = df.median[which(df.median$condition == "Gal15"),"zscores"],
                              Antimycin = df.median[which(df.median$condition == "Anti15"),"zscores"],
                              Pyruvate = df.median[which(df.median$condition == "Pyr15"),"zscores"])
rownames(df.median.by.cond) <- df.median.by.cond$gene
df.median.by.cond$gene <- NULL

df.median.by.cond <- rbind(df.median.by.cond, c( median(df.median.by.cond[grepl("OR", rownames(df.median.by.cond)), "Glucose"]),
                                             median(df.median.by.cond[grepl("OR", rownames(df.median.by.cond)), "Galactose"]),
                                             median(df.median.by.cond[grepl("OR", rownames(df.median.by.cond)), "Antimycin"]),
                                             median(df.median.by.cond[grepl("OR", rownames(df.median.by.cond)), "Pyruvate"])))
rownames(df.median.by.cond) <- c(rownames(df.median.by.cond)[1:(nrow(df.median.by.cond)-1)],"Control")
df.median.by.cond <- df.median.by.cond[!grepl("OR", rownames(df.median.by.cond)),]

pdf("./figs/gxe_heatmap.pdf",height=10.5,width=2.3)
# pheatmap(df.median.by.cond[order(df.median.by.cond$Glucose),], cluster_rows = F, cluster_cols=F,
#          breaks=seq(-1.6,1.6,l=150),
#          color=c(colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)[1:50],
#                  rep("white",50), colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)[50:100]),
#          cellwidth=10, cellheight=10)
pheatmap(df.median.by.cond[order(df.median.by.cond$Glucose,decreasing=T),], cluster_rows = F, cluster_cols=F,
         breaks=seq(-1.6,1.6,l=115),
         color=c(colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)[1:50],
         rep("white",15), colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)[50:100]),
         cellwidth=10, cellheight=10)
dev.off()

write.xlsx(df.median.by.cond, "df_median_by_cond.xlsx", row.names=T)

genes.of.interest = c("SLC25A19","SLC25A32")

pdf("./figs/violin_single_gene_by_condition_zscores.pdf", 
    # family="Helvetica", this is the default
    height=2.3,width=3.8)
for(gene in genes.of.interest) {
  columns.of.interest=c("Glc15","Gal15",
                        "Anti15","Pyr15")
  
  gene.subset <- df[which(df$gene == gene),]
  
  ctrl.subset <- df[grepl("OR",df$gene),]
  ctrl.subset$gene <- "ctrl"
  ctrl.subset$promoter <- "ctrl"
  
  total.subset <- rbind(gene.subset,ctrl.subset)
  total.subset$condition <- factor(total.subset$condition, ordered=T,
                         levels=columns.of.interest)
  total.subset$condition_ctrl <- paste0(total.subset$condition, "_", c("ko", "ctrl")[as.numeric(grepl("ctrl", total.subset$promoter)) + 1])
  total.subset$condition_ctrl <- factor(total.subset$condition_ctrl, ordered=T,
                              levels=as.vector(sapply(columns.of.interest, function(s) c(paste0(s, "_ctrl"),c(paste0(s, "_ko"))))))
  total.subset$ctrl <- c("ko", "ctrl")[as.numeric(grepl("ctrl", total.subset$promoter)) + 1]
  # variable = condition
  g <- ggplot(total.subset, aes(x=condition_ctrl,y=zscores,group = condition_ctrl,fill=ctrl)) + 
    geom_violin() + theme_classic() + xlab("") + ylab("z-scores") +
    scale_fill_manual(values=c("ctrl"="#999999","ko"="#ffa0b1")) +
    theme(axis.text.x=element_text(size=8, angle=45, vjust=0.7)) + 
    ggtitle(paste0(gene, "xCtrl"))
  print(g)
}
dev.off()



## control vs single gene by condition
single_ko_violins <- function(data,genes.of.interest=NA,
                              columns.of.interest=c("Glc15","Gal15",
                                                    "Anti15","Pyr15"),
                              ctrl.side.by.side=F, use.ctrl = T) {
  if(any(is.na(genes.of.interest))) {
    genes.of.interest <- as.character(c(read.csv("data/genelist.csv",stringsAsFactors = F,
                                                 header=F))[[1]])
  }
  
  sko.subset <- data[(data$U6.negctrl + data$H1.negctrl == 1),]
  sko.subset$promoter <- c("H1","U6")[as.factor(
    sko.subset$H1.negctrl == 1)]
  sko.subset$gene <- sko.subset$H1.gene
  sko.subset[which(sko.subset$promoter == "U6"),"gene"] <- 
    sko.subset[which(sko.subset$promoter == "U6"),"U6.gene"]
  
  ctrl.subset <- data[(data$U6.negctrl + data$H1.negctrl == 2),]
  ctrl.subset$gene <- "ctrl"
  ctrl.subset$promoter <- "ctrl"
  
  for(gene in genes.of.interest) {
    gene.subset <- sko.subset[which(sko.subset$H1.gene == gene |
                                      sko.subset$U6.gene == gene),]
    total.subset <- rbind(gene.subset,ctrl.subset)
    df2 <- melt(total.subset[,c("promoter",columns.of.interest)])
    df2$promoter <- factor(df2$promoter,ordered=T,
                          levels=c("U6","H1","ctrl"))
    df2$condition <- factor(df2$variable, ordered=T,
                           levels=columns.of.interest)
    df2$condition_ctrl <- paste0(df2$condition, "_", c("ko", "ctrl")[as.numeric(grepl("ctrl", df2$promoter)) + 1])
    df2$condition_ctrl <- factor(df2$condition_ctrl, ordered=T,
                                levels=as.vector(sapply(columns.of.interest, function(s) c(paste0(s, "_ctrl"),c(paste0(s, "_ko"))))))
    df2$ctrl <- c("ko", "ctrl")[as.numeric(grepl("ctrl", df2$promoter)) + 1]
    # variable = condition
    g <- ggplot(df2, aes(x=condition_ctrl,y=value,group = condition_ctrl,fill=ctrl)) + 
      geom_violin() + theme_classic() + xlab("") + ylab("LFC") +
      scale_fill_manual(values=c("ctrl"="#999999","ko"="#ffa0b1")) +
      theme(axis.text.x=element_text(size=8, angle=45, vjust=0.7)) + 
      ggtitle(paste0(gene, "xCtrl"))
    print(g)
  }
}

pdf("./figs/violin_single_gene_by_condition.pdf", 
    # family="Helvetica", this is the default
    height=2.3,width=3.8)
for(gene in genelist) {
  single_ko_violins(data, gene)
}
dev.off()



#### gxe density plots, removing the OR
df.no.or <- df #df[!grepl("OR", df$gene),]
df.no.or$gene <- gsub("OR(.*)+", "OR", df.no.or$gene)
df.no.or$condition <- as.character(df.no.or$condition)
# rename glucose to be AAglucose so it is the base factor
df.no.or$condition[which(df.no.or$condition == "Glc15")] = "aaGlc15"
fit <- lm(zscores ~ gene + gene:condition, data=df.no.or)
# write.table(summary(fit)$coefficients, file="results/fit_coefficients_singleko_condition_interaction.tsv",
#             sep="\t")
pvals = summary(fit)$coefficients[,4]
pvals.adj = p.adjust(pvals, method="fdr")

plot.new()
plot(-log10(pvals.adj) ~ summary(fit)$coefficients[,1])
points(-log10(pvals.adj)[grepl(":", names(pvals.adj))] ~ summary(fit)$coefficients[,1][grepl(":", names(pvals.adj))],
       pch=19,col="red")
head(rev(sort(-log10(pvals.adj[grepl(":", names(pvals.adj))]))), 50)

openxlsx::write.xlsx(summary(fit)$coefficients[which(summary(fit)$coefficients[,4] < 0.05 & 
                                  grepl(":", rownames(summary(fit)$coefficients), fixed=T) & 
                                  abs(summary(fit)$coefficients[,1]) > 0.25),,drop=F],
            file="top_hits_pval_0.05_effect_size_0.25.xlsx",row.names=T)


gs <- list()
for(n in 1:length(unique(df.no.or$gene))) {
  gene <- unique(df.no.or$gene)[n]
  gs[[gene]] <- ggplot(df.no.or[which(df.no.or$gene == gene),], 
                       aes(x=zscores, group=condition,
                           color=condition,fill=condition)) + 
    geom_density(alpha=0.3) + ggtitle(gene) + theme_classic()
  #print(gs[[n]])
}

pdf("./figs/density_plots.pdf", height=6,width=8)
for(i in seq(1, length(unique(df.no.or$gene)), by=4)) {
  if(i + 3 > length(unique(df.no.or$gene))) {
    grid.arrange(gs[[unique(df.no.or$gene)[i]]], gs[[unique(df.no.or$gene)[i+1]]],
                 gs[[unique(df.no.or$gene)[i+2]]],
                 nrow=2,
                 padding=10)
  } else {
    grid.arrange(gs[[unique(df.no.or$gene)[i]]], gs[[unique(df.no.or$gene)[i+1]]],
                 gs[[unique(df.no.or$gene)[i+2]]], gs[[unique(df.no.or$gene)[i+3]]],
                 nrow=2,
                 padding=10)
  }
}
dev.off()
