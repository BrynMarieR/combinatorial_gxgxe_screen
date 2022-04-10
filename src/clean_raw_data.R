# Author: Bryn Reinstadler
# Date: 2016 - 2017
# Purpose: To clean and prepare raw data (incl exploratory data analysis) for 
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

# Follows methods of "Combinatorial G x G x E CRISPR screen discovers 
# SLC25A39 in mitochondrial glutathione transport 
# linking iron homeostasis to OXPHOS"

raw.b1 <- openxlsx::read.xlsx("./raw_data/TableS1_Screen_Raw_Data.xlsx",
                              sheet=2)
raw.b2 <- openxlsx::read.xlsx("./raw_data/TableS1_Screen_Raw_Data.xlsx",
                              sheet=3)

colnames(raw.b2) <- gsub("p0","d0",colnames(raw.b2))
raw.b2$pDNA <- raw.b1$pDNA  # the two are identical and pDNA was not resubmitted

# ordering matches from b1 to b2
identical(raw.b2$H1.well, raw.b1$H1.well)

raw.b1$U6_H1 <- paste0(raw.b1$U6.gene, "_", raw.b1$H1.gene)
raw.b2$U6_H1 <- paste0(raw.b2$U6.gene, "_", raw.b2$H1.gene)

raw.b1$H1.negctrl <- grepl("OR|cutting", raw.b1$H1.gene)
raw.b1$U6.negctrl <- grepl("OR|cutting", raw.b1$U6.gene)
raw.b2$H1.negctrl <- grepl("OR|cutting", raw.b2$H1.gene)
raw.b2$U6.negctrl <- grepl("OR|cutting", raw.b2$U6.gene)

# read counts
g1 <- ggplot(raw.b1, aes(x=pDNA)) + geom_density() + 
  geom_vline(xintercept=500, color="red") + theme_classic() + 
  xlab("Read Count")
ggsave("rc_hist.pdf", g1, useDingbats=F)


hist(raw.b1$pDNA,breaks=1000, xlab="Read Counts",main="Histogram of Read Counts")
abline(v=500,col="red")

# pDNA is the same in both batches; remove guides with median rep < 500 in initial pDNA
u6.raw.b1.agg <- aggregate(pDNA ~ U6.sequence,data=raw.b1,FUN=median)
u6.raw.b1.agg[which(u6.raw.b1.agg$pDNA < 500),]
h1.raw.b1.agg <- aggregate(pDNA ~ H1.sequence,data=raw.b1,FUN=median)
h1.raw.b1.agg[which(h1.raw.b1.agg$pDNA < 500),]

bad.guides.from.pdna <- c(u6.raw.b1.agg[which(u6.raw.b1.agg$pDNA < 500),"U6.sequence"],
                          h1.raw.b1.agg[which(h1.raw.b1.agg$pDNA < 500),"H1.sequence"])

# we are interested in four conditions: Gal, Anti, Pyr, Glc 15 day (repeated between batches)
process_raw <- function(raw_data, 
                        columns.of.interest = c("pDNA","d0","Glc15","Gal15",
                                                "Anti15","Pyr15"),
                        color.special.genes=NA) {
  dataset.freqs <- raw_data
  for(x in columns.of.interest) {
    dataset.freqs[,x] <- raw_data[,x] / sum(raw_data[,x])
  }
  
  dataset.log2freqs <- dataset.freqs
  for(x in columns.of.interest) {
    dataset.log2freqs[,x] <- log2(dataset.freqs[,x]*1000000 + 1) # cannot have 0s
  }
  
  fc.data.raw <- dataset.log2freqs
  for(x in columns.of.interest) {
    fc.data.raw[,x] <- dataset.log2freqs[,x] - dataset.log2freqs[,"pDNA"]
  }
  return(fc.data.raw)
}

columns.of.interest = c("pDNA","d0","Glc15","Gal15",
                        "Anti15","Pyr15")

fc.data.b1 <- process_raw(raw.b1,columns.of.interest)
fc.data.b2 <- process_raw(raw.b2,columns.of.interest) #,color.special.genes=c("OR11H1","OR4F21"))

fc.data.b1[,columns.of.interest] <- sapply(columns.of.interest, function(y) {
  as.numeric(fc.data.b1[,y])
})
fc.data.b2[,columns.of.interest] <- sapply(columns.of.interest, function(y) {
  as.numeric(fc.data.b2[,y])
})

# filter uninteresting columns
fc.data.b1 <- fc.data.b1[,c("U6.gene","U6.well","H1.gene","H1.well",
                            "U6.sequence","H1.sequence",
                            "U6_H1",
                            columns.of.interest,
                            "U6.negctrl","H1.negctrl")]
fc.data.b2 <- fc.data.b2[,c("U6.gene","U6.well","H1.gene","H1.well",
                            "U6.sequence","H1.sequence",
                            "U6_H1",
                            columns.of.interest,
                            "U6.negctrl","H1.negctrl")]

# check basic correlations
for(col in columns.of.interest[-which(columns.of.interest == "pDNA")]) {
  tmp.df <- data.frame(b1 = fc.data.b1[,col],
                       b2 = fc.data.b2[,col])
  g <- ggplot(tmp.df, aes(x=b1, y=b2)) + geom_point(alpha=0.05,col="black") + 
    xlab("Batch 1 logFC") + ylab("Batch 2 logFC") + 
    theme_classic() + ggtitle(paste0(col,": Batch 1 and 2 Correlations")) + 
    stat_density_2d(aes(fill = ..level..), geom="polygon")
  ggsave(paste0("./figs/correlation_b12_",col,".pdf"), g, useDingbats=F)
}




# visualisations to confirm need for knock out
qc_control_double_knockouts <- function(data,columns.of.interest) {
  aggregate.ctrls.guide.level.u6 <- aggregate(. ~ U6.sequence + U6.gene, data[
    (data$H1.negctrl + data$U6.negctrl == 2),
    c("U6.sequence", "U6.gene",columns.of.interest
    )], FUN=median)
  cutoff <- max(abs(range(aggregate.ctrls.guide.level.u6[,columns.of.interest])))
  breakslist = seq(-cutoff, cutoff, length.out=100)
  colorpalette <- colorRampPalette(c("blue","white","red"),interpolate="spline")(length(breakslist))
  pheatmap(aggregate.ctrls.guide.level.u6[,3:ncol(aggregate.ctrls.guide.level.u6)],
           labels_row = paste0(#aggregate.ctrls.guide.level.u6$U6.sequence,"_",
                               aggregate.ctrls.guide.level.u6$U6.gene),
           main="Sp controls",cluster_cols = F,
           cellwidth = 10,cellheight = 10,
           breaks=breakslist,color=colorpalette)
  
  aggregate.ctrls.guide.level.h1 <- aggregate(. ~ H1.sequence + H1.gene, data[
    (data$H1.negctrl + data$U6.negctrl == 2),
    c("H1.sequence", "H1.gene", columns.of.interest
    )], FUN=median)
  cutoff <- max(abs(range(aggregate.ctrls.guide.level.h1[,columns.of.interest])))
  breakslist = seq(-cutoff, cutoff, length.out=100)
  pheatmap(aggregate.ctrls.guide.level.h1[,3:ncol(aggregate.ctrls.guide.level.h1)],
           labels_row = paste0(#aggregate.ctrls.guide.level.h1$H1.sequence,"_",
                               aggregate.ctrls.guide.level.h1$H1.gene),
           main="Sa controls",cluster_cols = F,
           cellwidth = 10,cellheight = 10,
           breaks=breakslist,color=colorpalette)
}

pdf("./figs/bad_OR_heatmaps.pdf")
qc_control_double_knockouts(rbind(fc.data.b1, fc.data.b2),columns.of.interest) 
dev.off()

diagnostic_histograms <- function(dataset,columns.of.interest) {
  melted.ds <- melt(dataset[,columns.of.interest])
  melted.ds$variable <- factor(melted.ds$variable,ordered=T,
                               levels=c(
                                 "pDNA","d0","Glc15","Gal15","Pyr15","Anti15"
                               ))
  print(ggplot(melted.ds, aes(x=value, colour=variable, fill=variable)) + geom_density(alpha=0.2)
        + theme_classic() + xlab("LogFC") + ylab("Density") + xlim(c(-4,4)) + 
          scale_color_manual(values=c("pDNA"="red","d0"="orange",
                                      "Glc15"="#440154FF", "Gal15"="#31688EFF",
                                      "Anti15"="#FDE725FF","Pyr15"="#35B779FF")) + 
          scale_fill_manual(values=c("pDNA"="red","d0"="orange",
                                      "Glc15"="#440154FF", "Gal15"="#31688EFF",
                                      "Anti15"="#FDE725FF","Pyr15"="#35B779FF")))
}

# unrepresented after filtering < 500, same in b2 since pDNA is the same
setdiff(unique(fc.data.b1[,"U6_H1"]), unique(fc.data.b1[which(raw.b1$pDNA >= 500),"U6_H1"]))
# [1] "UCP2_OR7E37P-chr13"        "OR7G3-chr19_OR7E37P-chr13"

# read count filtering is based on empirical histograms
# decide not to use non-cutting controls
filtered.fc.data.b1 <- fc.data.b1[!(grepl("OR11H1|OR4F21|cutting",fc.data.b1$U6_H1)) & 
                                    (raw.b1$pDNA >= 500),]
filtered.fc.data.b2 <- fc.data.b2[!(grepl("OR11H1|OR4F21|cutting",fc.data.b2$U6_H1)) & 
                                    (raw.b2$pDNA >= 500),]

comb.b12 <- rbind(fc.data.b1, fc.data.b2)
diagnostic_histograms(comb.b12[grepl("cutting", comb.b12$U6.gene) &
                                   grepl("cutting", comb.b12$H1.gene),
                                 columns.of.interest])
diagnostic_histograms(comb.b12[grepl("OR", comb.b12$U6.gene) &
                                   grepl("OR", comb.b12$H1.gene),
                                 columns.of.interest])


plot_guide_corr_heatmaps <- function(data, genes.of.interest=NA,
                                     columns.of.interest=c("pDNA","d0","Glc15","Gal15",
                                                           "Anti15","Pyr15"),
                                     return.poor.guides=F) {
  if(any(is.na(genes.of.interest))) {
    genes.of.interest <- as.character(c(read.csv("./raw_data/genelist.csv",stringsAsFactors = F,
                                                 header=T))[[1]])
  }
  poor.guides <- rep(NA, nrow(data))
  poor.guides.index <- 1
  par(mfrow=c(1,1))
  for (gene in genes.of.interest) {
    subset.u6 <- data[(data$U6.gene == gene) & 
                        (data$H1.negctrl == 1),]
    subset.h1 <- data[(data$H1.gene == gene) & 
                        (data$U6.negctrl == 1),]
    total.comb <- rbind(subset.u6, subset.h1)
    
    if(dim(subset.h1)[1] == 0 | dim(subset.u6)[1] == 0) {
      print(gene)
      print("guides from one promoter missing entirely")
      filtered.cormat.means <- c(filtered.cormat.means, NA)
      prefilt.means <- c(prefilt.means, NA)
      worst.two.means <- c(worst.two.means, NA)
    } else {
      guide.agg.u6 <- aggregate(. ~ U6.sequence, 
                                subset.u6[,c("U6.sequence",columns.of.interest)], FUN=mean)
      rownames(guide.agg.u6) <- paste0("U6_",guide.agg.u6$U6.sequence)
      guide.agg.u6$U6.sequence <- NULL
      
      guide.agg.h1 <- aggregate(. ~ H1.sequence, 
                                subset.h1[,c("H1.sequence",columns.of.interest)], mean)
      rownames(guide.agg.h1) <- paste0("H1_",guide.agg.h1$H1.sequence)
      guide.agg.h1$H1.sequence <- NULL
      
      
      
      if(length(columns.of.interest) < 3) { # we can't do proper correlations; just do normal heatmap
        max.lfc <- max(abs(unlist(rbind(guide.agg.u6, guide.agg.h1)))) + 0.5 # gives tolerance for things near zero
        breakslist = seq(-round(max.lfc), round(max.lfc), length.out=11)
        colorpalette <- colorRampPalette(c("blue","white","red"), interpolate="spline")(length(breakslist))
        pheatmap(rbind(guide.agg.u6, guide.agg.h1),main=gene,
                 color=colorpalette, breaks=breakslist,
                 cellheight = 20, cellwidth=20)
        poor.guides.for.gene <- c()
      } else {
        breakslist = seq(-1, 1, by=0.1)
        colorpalette <- colorRampPalette(c("blue","white","red"), interpolate="spline")(length(breakslist))
        pheatmap(cor(t(rbind(guide.agg.u6, guide.agg.h1))),main=gene,
                 color=colorpalette, breaks=breakslist,
                 cellheight = 20, cellwidth=20)
        poor.guides.for.gene <- 
          names(which(apply(cor(t(rbind(guide.agg.u6, guide.agg.h1))), 1, sum) < 0))
      }
      
      
      if(length(poor.guides.for.gene) > 0 ) {
        poor.guides[poor.guides.index:(poor.guides.index + length(poor.guides.for.gene) - 1)] <- 
          poor.guides.for.gene
        poor.guides.index <- poor.guides.index + length(poor.guides.for.gene)
      }
    }
  }
  if(return.poor.guides & length(columns.of.interest) >= 3) {
    return(gsub("^(.*)+_","",poor.guides[!is.na(poor.guides)]))
  } else if (return.poor.guides & (length(columns.of.interest) < 3)) {
    print("With so few columns, it is impossible to find correlation between guides.")
  }
}

plot_guide_corr_heatmaps(filtered.fc.data.b1)
plot_guide_corr_heatmaps(filtered.fc.data.b2)


b1.bad.guides <- plot_guide_corr_heatmaps(filtered.fc.data.b1, return.poor.guides = T)
b2.bad.guides <- plot_guide_corr_heatmaps(filtered.fc.data.b2, return.poor.guides = T)

length(b1.bad.guides)
# [1] 40
length(b2.bad.guides)
# [1] 39
length(intersect(b1.bad.guides, b2.bad.guides))
# [1] 29
merge.guides <- unique(c(b1.bad.guides, b2.bad.guides, bad.guides.from.pdna))

saveRDS(filtered.fc.data.b1, "./processed_data/unfiltered.b1.rds")
saveRDS(filtered.fc.data.b2, "./processed_data/unfiltered.b2.rds")

# remove guides that were bad in either batch just for safety
filtered.guides.fc.data.b1 <- filtered.fc.data.b1[-which(filtered.fc.data.b1$U6.sequence %in% merge.guides | 
                                                           filtered.fc.data.b1$H1.sequence %in% merge.guides),]
filtered.guides.fc.data.b2 <- filtered.fc.data.b2[-which(filtered.fc.data.b2$U6.sequence %in% merge.guides | 
                                                           filtered.fc.data.b2$H1.sequence %in% merge.guides),]

plot_guide_corr_heatmaps(filtered.guides.fc.data.b1)
plot_guide_corr_heatmaps(filtered.guides.fc.data.b2)

saveRDS(filtered.guides.fc.data.b1,"./processed_data/clean.b1.rds")
saveRDS(filtered.guides.fc.data.b2,"./processed_data/clean.b2.rds")

filtered.guides.fc.data.b1$batch <- c(1)
filtered.guides.fc.data.b2$batch <- c(2)
combined.filtered.guides.fc.data <- rbind(filtered.guides.fc.data.b1, 
                                          filtered.guides.fc.data.b2)

saveRDS(combined.filtered.guides.fc.data,"./processed_data/clean.b12.rds")

diagnostic_histograms(combined.filtered.guides.fc.data[grepl("OR", combined.filtered.guides.fc.data$U6.gene) &
                                 grepl("OR", combined.filtered.guides.fc.data$H1.gene),
                               columns.of.interest])
