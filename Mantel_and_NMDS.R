# for now, the dist files are produced by Diversity-EcoLog.R and Phylo_new respectively
# check these scripts for details


# MANTEL TEST AND NMDS PLOTS #

# this script takes: 
# --> the distances matrices from the OTU tables and the Ecolog plates 
# --> the sample metadata

# this script does:
# --> Mantel tests 
# --> NMDS plots
# --> Mantel plots

# this script exports:
# --> plots

# load packages
###############
library(ggplot2)
###############

# import files
###############
#weighted unifrac distance based on rarefied OTU table
load("unifrac_dist.RData")

#weighted jaccard distance based on rarefied OTU table
load("OTUdist.RData")

#euclidian distance Ecolog carbnon source utilization (Y/N)
load("NLSdist.RData")

#euclidian distance Ecolog carbon uptake rate
load("EcoLYNdist.RData")

# Importing metadata
ID<-read.table("gbgID",header=TRUE, stringsAsFactors=F)
###############

# find common samples in all matrices
Common_samples <- Reduce(intersect, list( labels( unifrac_dist), labels( OTUdist), 
                        labels( NLSdist), labels( EcoLYNdist)))

# subset all distance matrices to common samples
dist.list <- list( unifrac_dist, OTUdist, NLSdist, EcoLYNdist)

# function to exclude samples and sort remaining samples
exclude_samples <- function(x) {
  mat <- as.matrix(x)
  W <- which( dimnames( mat)[[ 1 ]] %in% Common_samples)
  mat <- mat[W,W]
  Common.sorted <- sort( Common_samples)
  mat <- mat[ Common.sorted, Common.sorted]
  mat <- as.dist( mat)
  return(mat)
}

#exclude samples
dist.list <- lapply(dist.list, exclude_samples)

#reassign dist objects 
unifrac_dist  <- dist.list[[1]]
OTUdist <- dist.list[[2]]
NLSdist <- dist.list[[3]]
EcoLYNdist <- dist.list[[4]]






