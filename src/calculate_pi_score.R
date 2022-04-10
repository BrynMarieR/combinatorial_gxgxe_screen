# Author: Bryn Reinstadler
# Date: 2016 - 2017
# Purpose: To calculate the pi score (interaction score; epistasis score) for 
#    the genes of interest in a dual-knock-out combinatorial CRISPR screen.


# Libraries & Packages

library(pheatmap)
library(RColorBrewer)
library(gplots)


# Helper functions

makeRects <- function(tfMat,border){
  cAbove = expand.grid(1:dim(tfMat)[1],1:dim(tfMat)[2])[tfMat,]
  xl=cAbove[,1]-0.49
  yb=cAbove[,2]-0.49
  xr=cAbove[,1]+0.49
  yt=cAbove[,2]+0.49
  rect(xl,yb,xr,yt,border=border,lwd=2)
}

pi_score <- function(lfc.data, condition, genelist, sig.test = "parametric", symmetrize=F) {
  # "lfc.data": dataframe with the following columns:
  # U6.gene H1.gene U6.sequence H1.sequence condition U6.negctrl H1.negctrl batch
  #
  # The columns are mostly self-explanatory, except note that the condition column 
  # should be the logfc from the pDNA level and the U6.negctrl and H1.negctrl
  # should be 0 or 1 depending on whether the U6 gene or the H1 gene are a
  # negative control in your dataset. (1 indicates a negative control)
  # 
  # "condition": string which is the same as the name for your column with the LogFC in them.
  #
  # "genelist": vector that is used to build the output matrix. Typically doesn't include control genes,
  # although the controls are used for estimation of the interaction term.
  #
  # The method was originally developed as in this paper:
  # http://www.nature.com.ezp-prod1.hul.harvard.edu/articles/nmeth.1581
  # d_ijk = w + m_i + m'_j + pi_ij + e_ijk
  # d_ijk is the log of the k-th measurement
  # of the double-knock-out of genes i and j,
  # w is the log of the quant phenotype in 
  # undisturbed cells, and pi_ij is the pairwise
  # interaction term.  m and m' are the single knockout phenotypes.
  
  # see also the code repository from the above paper: 
  # http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.696.294&rep=rep1&type=pdf
  
  if(!("batch" %in% colnames(lfc.data))) {
    print("WARNING: no 'batch' column detected...")
    print("Assuming only one batch (MAY CAUSE ERRORS IF >1 BATCH!)...")
    lfc.data$batch = c(1)
  }
  
  pi.mat.by.guide <- list(max(lfc.data$batch))
  for(b in 1:max(lfc.data$batch)) {
    # for each batch, put all the data into a single wide matrix called "z"
    tmp.lfc <- lfc.data[which(lfc.data$batch == b),]
    t.z <- tmp.lfc[,c("U6.sequence","H1.sequence",condition)]
    z = stats::reshape(t.z, direction="wide", 
                       idvar="U6.sequence", timevar="H1.sequence")
    rownames(z) <- z$U6.sequence
    z$U6.sequence <- NULL
    colnames(z) <- gsub(paste0(condition,"."),"",colnames(z))
    u6.neg.rows <- which(rownames(z) %in% 
                           tmp.lfc[which(lfc.data$U6.negctrl == 1),"U6.sequence"])
    h1.neg.cols <- which(colnames(z) %in% 
                           tmp.lfc[which(lfc.data$H1.negctrl == 1),"H1.sequence"])
    
    # keep track of the number of rows and columns in z
    nr <- nrow(z)
    nc <- ncol(z)
    
    # this is a vestige of the old code and is kept because it may be useful for future
    # implementations. We only have one "template plate" (=="TP") but this 
    # vector may be used to control for within-batch differences
    TP <- rep(1,nr)
    
    # the "query" negative controls "QueryNeg" are the h1.neg columns
    QueryNeg <- h1.neg.cols
    # the target/template negative controls "TemplateNeg" are the u6.neg rows
    TemplateNeg <- u6.neg.rows
    
    # initalize variables for Tukey's median polish procedure
    t <- 0
    r <- numeric(nr) # initializes empty vector of length nr
    c <- numeric(nc)
    oldsum <- 0
    
    maxiter <- 100 # can be tweaked if desired, but I have never had any issues w convergence
    na.rm = TRUE
    eps = 1e-4 # convergence criteria; can be tweaked if desired but I have never had any issues 
    for (iter in 1L:maxiter) {
      # this is the row delta
      rdelta <- apply(z, 1L, median, na.rm = na.rm)
      z <- z - matrix(rdelta, nrow = nr, ncol = nc)
      r <- r + rdelta
      
      # this is the column delta
      cdelta <- apply(z, 2L, function(s) {
        tapply(s, TP, median, na.rm = na.rm)
      })
      # move to sapply because there's only one column delta
      # and only one "template plate"
      z <- z - matrix(cdelta, nrow=nr, ncol=nc,byrow=T)
      c <- c + cdelta
      if (is.null(QueryNeg)) {
        delta = median(c, na.rm = na.rm)
        if (!is.finite(delta)) {
          delta = 0
        }
      } else {
        delta <- mean(c[QueryNeg], na.rm = TRUE)
        if (!is.finite(delta)) {
          delta = median(c, na.rm = na.rm)
          if (!is.finite(delta)) {
            delta = 0
          }
        }
      }
      c <- c - delta
      t <- t + delta
      if (is.null(TemplateNeg)) {
        delta = median(r, na.rm = na.rm)
        if (!is.finite(delta)) {
          delta = 0
        }
      } else {
        delta <- mean(r[TemplateNeg], na.rm = TRUE)
        if (!is.finite(delta)) {
          delta = median(r, na.rm = na.rm)
          if (!is.finite(delta)) {
            delta = 0
          }
        }
      }
      r <- r - delta
      t <- t + delta
      newsum <- sum(abs(z), na.rm = na.rm)
      converged <- newsum == 0 || abs(newsum - oldsum) < eps * 
        newsum
      if (converged) { break; }
      oldsum <- newsum
    }
    if (!converged) {
      warning(gettextf("median polish procedure did not converge in %d iterations", 
                       maxiter), domain = NA)
    } else {
      print("Successful convergence!")
    }
    
    pi.mat.by.guide[[b]] <- z
  }
  # collect the mean of all of the guide-by-guide scores in the "GI data frame" = gi.df
  # only do so for genes in the "genelist" which typically doesn't include control genes,
  # although the controls were used for estimation of the interaction term.
  gi.df <- data.frame(matrix(
    rep(NA, length(genelist)*length(genelist)),
    nrow = length(genelist), ncol=length(genelist)
  ))
  gi.array <- array(data = c(gi.df, gi.df), dim=c(length(genelist), length(genelist), 2))
  
  rownames(gi.array) <- genelist
  colnames(gi.array) <- genelist
  
  for(generow in rownames(gi.array)) {
    for(genecol in colnames(gi.array)) {
      generow.guides = c(unique(lfc.data[which(lfc.data$U6.gene == generow), "U6.sequence"]),
                         unique(lfc.data[which(lfc.data$H1.gene == generow), "H1.sequence"]))
      genecol.guides = c(unique(lfc.data[which(lfc.data$H1.gene == genecol), "H1.sequence"]),
                         unique(lfc.data[which(lfc.data$U6.gene == genecol), "U6.sequence"]))
      
      values = c()
      values.batches = c()
      if(!symmetrize) {
        for(b in 1:max(lfc.data$batch)) {
          new.values <- unlist(c(pi.mat.by.guide[[b]][
                       generow.guides[which(generow.guides %in% rownames(pi.mat.by.guide[[b]]))], 
                       genecol.guides[which(genecol.guides %in% colnames(pi.mat.by.guide[[b]]))]]))
          values <- c(values, new.values)
          values.batches <- c(values.batches, rep(b, length(new.values)))
        }
      } else {
        for(b in 1:max(lfc.data$batch)) {
          new.values <- unlist(c(pi.mat.by.guide[[b]][
            generow.guides[which(generow.guides %in% rownames(pi.mat.by.guide[[b]]))], 
            genecol.guides[which(genecol.guides %in% colnames(pi.mat.by.guide[[b]]))]]))
          new.values <- c(new.values, unlist(
            c(pi.mat.by.guide[[b]][
              genecol.guides[which(genecol.guides %in% rownames(pi.mat.by.guide[[b]]))], 
              generow.guides[which(generow.guides %in% colnames(pi.mat.by.guide[[b]]))]])
          ))
          values <- c(values, new.values)
          values.batches <- c(values.batches, rep(b, length(new.values)))
        }
      }
      
      # currently doesn't adjust for batch, but info is stored in 'values.batches'
      gi.array[generow, genecol, 1] <- median(values,na.rm=T)
      if(sig.test == "parametric") {
        gi.array[generow, genecol, 2] <- t.test(values)$p.value
      } else {
        gi.array[generow, genecol, 2] <- wilcox.test(values)$p.value
      }
    }
  }
  
  gi.array <- apply(gi.array, c(1,2,3), as.numeric)
  
  if(!symmetrize) {
    gi.array[,,2] <- p.adjust(gi.array[,,2], "BH")
  } else {
    # did many fewer tests, only penalize those tests
    tmp <- gi.array[,,2]
    gi.array[,,2][lower.tri(gi.array[,,2],diag=T)] <- 
      p.adjust(tmp[lower.tri(tmp,diag=T)], "BH")
    gi.array[,,2][upper.tri(gi.array[,,2],diag=T)] <- 
      p.adjust(tmp[upper.tri(tmp,diag=T)], "BH")
  }
  
  return(gi.array)
}

# working minimal example: uncomment to run
# # set your own working directory
# lfc.data <- readRDS("./processed_data/clean.b12.rds")
# 
# final.gi.df <- pi_score(lfc.data, "Normal", 
#                         unique(lfc.data[which(lfc.data$U6.negctrl == 0),"U6.gene"]),
#                         sig.test = "parametric",
#                         symmetrize=T)
# 
# # Highlight FDR < 0.1
# sig <- final.gi.df[,,2] < 0.1
# heatmap.2(final.gi.df[,,1], Rowv=F, Colv=F, dendrogram = "none",
#           col = colorRampPalette(c("blue","white","red"))(50),
#           add.expr={makeRects(sig[,rev(1:ncol(sig))],"black")}, 
#           main="Heatmap of pi-scores, FDR < 0.1",
#           trace="none")