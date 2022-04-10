# Author: Bryn Reinstadler
# Date: 2016 - 2017
# Purpose: To prepare visualizations of quality control metrics for 
#    the genes of interest in a dual-knock-out combinatorial CRISPR screen.


# Libraries & Packages

library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
library(limma)
library(reshape)
library(openxlsx)
library(pheatmap)
library(RColorBrewer)
library(gtools)

data = readRDS("./processed_data/clean.b12.rds")

cor.test(data[which(data$batch == 1),"Glc15"], data[which(data$batch == 2),"Glc15"])

## correlation between guides of different genes
u6.genes <- aggregate(Glc15 ~ U6.gene, FUN=median,data=data[which(data$H1.negctrl == TRUE),])
h1.genes <- aggregate(Glc15 ~ H1.gene, FUN=median,data=data[which(data$U6.negctrl == TRUE),])
tmp.df.corr.guides <- data.frame(Sp = u6.genes$Glc15,
                                 Sa = h1.genes$Glc15,
                                 gene = h1.genes$H1.gene)
pdf("./figs/correlation_sp_sa.pdf", 
    # family="Helvetica", this is the default
    height=6.5,width=6.5, useDingbats=F)
g <- ggplot(tmp.df.corr.guides, aes(x=Sp, y=Sa)) + 
  coord_fixed(ratio=1, xlim=c(-3, 1.1), ylim=c(-3,1.1)) + 
  theme(aspect.ratio=1) + geom_point() + theme_classic() +
  xlab("Sp median logFC") + ylab("Sa median logFC") + ggtitle("Correlation between Sp and Sa guides per gene") + 
  geom_abline(intercept=0, slope=1, col="red") + 
  geom_text_repel(data=subset(tmp.df.corr.guides, abs(tmp.df.corr.guides[,1] - tmp.df.corr.guides[,2]) > 1),
                   label=subset(tmp.df.corr.guides, abs(tmp.df.corr.guides[,1] - tmp.df.corr.guides[,2]) > 1)$gene)
print(g)
dev.off()

cor.test(tmp.df.corr.guides$Sp, tmp.df.corr.guides$Sa)

## Correlation of targeting w different controls
# for each position separately
sko.only.h1 <- data[which(data$H1.negctrl == 1 & data$U6.negctrl == 0),]
sko.only.u6 <- data[which(data$H1.negctrl == 0 & data$U6.negctrl == 1),]
sko.only.h1.agg <- aggregate(Glc15 ~ U6.sequence + H1.gene, data=sko.only.h1, FUN=median)
sko.only.u6.agg <- aggregate(Glc15 ~ H1.sequence + U6.gene, data=sko.only.u6, FUN=median)

mat.sko.only.h1.agg <- reshape(sko.only.h1.agg, direction="wide", timevar="H1.gene", idvar="U6.sequence")
rownames(mat.sko.only.h1.agg) <- mat.sko.only.h1.agg$U6.sequence
mat.sko.only.h1.agg$U6.sequence <- NULL
# remove two bad guides that had been filtered
mat.sko.only.h1.agg$`Glc15.OR6B3-chr02` <- 
  mat.sko.only.h1.agg$`Glc15.OR7E37P-chr13` <- NULL

pheatmap(cor(mat.sko.only.h1.agg),cluster_rows=F,cluster_cols=F,
         breaks=seq(-1,1,length.out=11),color=rev(brewer.pal(11,"RdBu")))
h1.mean <- mean(cor(mat.sko.only.h1.agg), na.rm=T)
h1.sd <- sd(cor(mat.sko.only.h1.agg), na.rm=T)

mat.sko.only.u6.agg <- reshape(sko.only.u6.agg, direction="wide", 
                               timevar="U6.gene", idvar="H1.sequence")
rownames(mat.sko.only.u6.agg) <- mat.sko.only.u6.agg$H1.sequence
mat.sko.only.u6.agg$H1.sequence <- NULL
pheatmap(cor(mat.sko.only.u6.agg),cluster_rows=F,cluster_cols=F,
         breaks=seq(-1,1,length.out=11),color=rev(brewer.pal(11,"RdBu")))
u6.mean <- mean(cor(mat.sko.only.u6.agg), na.rm=T)
u6.sd <- sd(cor(mat.sko.only.u6.agg), na.rm=T)

t.test(as.vector(cor(mat.sko.only.u6.agg)),as.vector(cor(mat.sko.only.h1.agg)))

pdf("./figs/correlation_guides_diff_ctrls.pdf", 
    # family="Helvetica", this is the default
    height=3.8,width=3.8, useDingbats=F)
par(mfrow=c(1,1))
plot(c(h1.mean, u6.mean) ~ c(1,1.5), xlim=c(0.75,1.75),
     ylim=c(-0.2, 1),pch=19,col="red",
     xaxt="n",xlab="Promoter",ylab="Pearson Correlation",
     main="Pearson correlation between \n targeting guides with different controls")
abline(h=0, lty=2)
axis(side=1, at=c(1,1.5), labels=c("Sa","Sp"),tick=T)
segments(x0=c(1,1.5),y0=c(h1.mean - h1.sd, u6.mean-u6.sd),
         x1=c(1,1.5),y1=c(h1.mean + h1.sd, u6.mean+u6.sd))
dev.off()

#### A51 and A46 guide issues
dat.a51.sko <- data[(data$U6.gene == "SLC25A51" & data$H1.negctrl) |
                      (data$H1.gene == "SLC25A51" & data$U6.negctrl),]
dat.a51.sko$promoter <- c("Sp")
dat.a51.sko$promoter[which(dat.a51.sko$U6.gene == "SLC25A51")] = "Sa"

a51.df <- melt(dat.a51.sko[,c("promoter","Glc15","Gal15","Anti15","Pyr15")])
colnames(a51.df) <- c("promoter","condition","value")
a51.df$promoter_condition <- paste0(a51.df$promoter, a51.df$condition)
a51.df$promoter_condition <- factor(a51.df$promoter_condition, 
                                    levels=c("SpGlc15","SaGlc15",
                                             "SpGal15","SaGal15",
                                             "SpAnti15","SaAnti15",
                                             "SpPyr15","SaPyr15"))
a51.df$promoter <- factor(a51.df$promoter,levels=c("Sp","Sa"))

g <- ggplot(a51.df, aes(x=promoter_condition,y=value,group=promoter_condition,color=promoter)) + 
  geom_quasirandom() + xlab("") + ylab("LFC") + theme_classic() + 
  scale_color_manual(values=c("#007fad","#00BFC4")) +
  theme(axis.text.x=element_text(size=8, angle=45, vjust=0.7)) + 
  ggtitle("SLC25A51xCtrl")
print(g)
ggsave("./figs/slc25a51_beeswarms.pdf", g,height=5, width=7, useDingbats=F)

dat.a46.sko <- data[(data$U6.gene == "EEF2" & data$H1.negctrl) |
                      (data$H1.gene == "EEF2" & data$U6.negctrl),]
dat.a46.sko$promoter <- c("Sp")
dat.a46.sko$promoter[which(dat.a46.sko$H1.gene == "EEF2")] = "Sa"

a46.h1.df <- melt(dat.a46.sko[which(dat.a46.sko$promoter == "Sa"),c("H1.sequence","promoter","Glc15","Gal15","Anti15","Pyr15")])
colnames(a46.h1.df) <- c("guide","promoter","variable","value")
a46.u6.df <- melt(dat.a46.sko[which(dat.a46.sko$promoter == "Sp"),c("U6.sequence","promoter","Glc15","Gal15","Anti15","Pyr15")])
colnames(a46.u6.df) <- c("guide","promoter","variable","value")
a46.df <- rbind(a46.h1.df,a46.u6.df)
colnames(a46.df) <- c("guide","promoter","condition","value")

a46.df$promoter_condition <- paste0(a46.df$promoter, a46.df$condition)
a46.df$promoter_condition <- factor(a46.df$promoter_condition, 
                                    levels=c("SpGlc15","SaGlc15",
                                             "SpGal15","SaGal15",
                                             "SpAnti15","SaAnti15",
                                             "SpPyr15","SaPyr15"))
a46.df$promoter <- factor(a46.df$promoter,levels=c("Sp","Sa"))

g <- ggplot(a46.df, aes(x=promoter_condition,y=value,group=promoter_condition,color=promoter)) + 
  geom_quasirandom() + theme_classic() + xlab("") + ylab("LFC") +
  scale_color_manual(values=c("#007fad","#00BFC4")) +
  theme(axis.text.x=element_text(size=8, angle=45, vjust=0.7)) + 
  ggtitle("EEF2xCtrl")
print(g)
ggsave("./figs/eef2_beeswarms.pdf", g,height=5, width=7, useDingbats=F)
