#!/bin/env R

# Determine the optimal order of leaves for the clustering of this tree.

library(Matrix)
library(reshape)
library(gplots)

determine_order <- function(D, method, outdir){

  hc <- hclust(D, method=method);
  write.csv(hc$order, file=sprintf('%s/network_order.%s.tsv', outdir, method), row.names=FALSE, col.names=FALSE, quote=FALSE)

}

SM <- read.csv(file="analysis/diffusion_input_network_sparse.tsv",head=FALSE,sep="\t")
SM_r <- sparseMatrix(SM[,1], SM[,2], x=SM[,3])
M = 1 - as.matrix(SM_r)
D = as.dist(M)

# http://www.inside-r.org/packages/cran/cba/docs/order.optimal

methods = c("ward.D", "single", "complete", "average", "median", "centroid");

apply(methods, 1, function(method){ determine_order(D, method, 'analysis') })

method='single';

O = read.csv(sprintf('analysis/network_order.%s.tsv', method, head=TRUE, sep='\t'))
S = read.csv('analysis/all_diffused_scores.tsv', head=TRUE, sep='\t')

O <- transform(O, x=as.numeric(x))

colnames(S) <- c("beta", 'name', 'score');
S <- transform(S, score=as.numeric(score))


heatmap_ready <- as.matrix(cast(S, beta ~ name))
HM_norm      <- t(apply(heatmap_ready, 1, function(x)(x-min(x))/(max(x)-min(x))))

jpeg('test.png', width=3000, height=3000)
heatmap(HM_norm, dendogram="none", Colv=O$x, Rowv=1:dim(heatmap_ready)[1]);
dev.off()


