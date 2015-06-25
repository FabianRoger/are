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
library(reshape2)
library(gridExtra)
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

##### subset all distance matrices to common samples ###########################

dist.list <- list( unifrac_dist = unifrac_dist, OTUdist = OTUdist, 
                   NLSdist = NLSdist, EcoLYNdist = EcoLYNdist)

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

####### create dataframe that labels the dist combinations by Lake #############

Comb <- expand.grid( ID$gbgID, ID$gbgID)
Comb <- paste( Comb[ ,1], Comb[ ,2], sep = "|")
Comb <- data.frame( pairs = Comb, Lake = NA, DAT = NA, LakeDat = NA ,
                    stringsAsFactors = F)

# assign pairs within each Lake
for (i in levels( as.factor( ID$Lake))) {
  LakeComb <- expand.grid( ID[ ID$Lake == i, ]$gbgID, ID[ ID$Lake == i, ]$gbgID)
  LakeComb <- paste( LakeComb[ ,1], LakeComb[ ,2], sep = "|")
  Comb[ Comb$pairs %in% LakeComb, ]$Lake <- i
  }

# assign pairs within each Date
for (i in levels( as.factor( ID$DAT))) {
  DATComb <- expand.grid( ID[ ID$DAT == i, ]$gbgID, ID[ ID$DAT == i, ]$gbgID)
  DATComb <- paste( DATComb[ ,1], DATComb[ ,2], sep = "|")
  Comb[ Comb$pairs %in% DATComb, ]$DAT <- i
}

# assign pairs within each DatexLake combination
for ( i in levels( as.factor( ID$DAT))) {
  for ( n in levels( as.factor( ID$Lake))) {
    DATLakeComb <- expand.grid( ID[ ID$DAT == i & ID$Lake == n, ]$gbgID,
                                ID[ ID$DAT == i & ID$Lake == n, ]$gbgID)
    DATLakeComb <- paste( DATLakeComb[ ,1], DATLakeComb[ ,2], sep = "|")
    Comb[ Comb$pairs %in% DATLakeComb, ]$LakeDat <- paste( i, n, sep="x")
  }
}

# assign unassigned combinations as "full"
Comb[ is.na( Comb)] <- "full"

####### summerize all distance matrices to one dataframe with assignments ######


for (i in 1:4) {
  dist.name <- names(dist.list)[i]
  dist.df <- data.frame(melt(as.matrix(dist.list[[i]])))
  dist.df$pairs <- paste(dist.df$Var1, dist.df$Var2, sep = "|")
  dist.df <- dist.df[ , c( "value", "pairs")]
  colnames( dist.df)[ 1] <- dist.name
  assign(paste(dist.name,"DF",sep="."), dist.df)
}


# join dataframes
distanceDF <- join( join( join( NLSdist.DF, OTUdist.DF), unifrac_dist.DF), 
                    EcoLYNdist.DF)

#add pair assignments
distanceDF <- join( distanceDF, Comb, type = "left")

############## plot distance scatterplots ######################################

##### unifrac_dist  ~ EcoLYNdist #####

# exclude among Lake distances and joined 0 distances  
distanceDF_Lake <- distanceDF[ which(apply(distanceDF[,c( "unifrac_dist", "EcoLYNdist")], 1, sum) > 0),]
distanceDF_Lake <- distanceDF_Lake[distanceDF_Lake$Lake != "full", ]

G2 <- ggplot( distanceDF_Lake, aes( x = unifrac_dist, y = EcoLYNdist, 
                              colour = Lake, shape = Lake))+
  geom_point( alpha = 0.3)+
  facet_wrap(~Lake)+
  stat_smooth( method = "lm", se = F, size = 1)+
  theme_bw( base_size = 15)+
  theme( legend.position = "none")+
  scale_colour_manual( values = c( "orange", "darkred", "darkgreen", "darkblue")) +
  labs( x="", y= "")

##### unifrac_dist  ~ NLSdist #####

# exclude among Lake distances and joined 0 distances  
distanceDF_Lake <- distanceDF[ which(apply(distanceDF[,c( "unifrac_dist", "NLSdist")], 1, sum) > 0),]
distanceDF_Lake <- distanceDF_Lake[distanceDF_Lake$Lake != "full", ]

G4 <- ggplot( distanceDF_Lake, aes( x = unifrac_dist, y = NLSdist, 
                              colour = Lake, shape = Lake))+
  geom_point( alpha = 0.3)+
  facet_wrap(~Lake)+
  stat_smooth( method = "lm", se = F, size = 1)+
  theme_bw( base_size = 15)+
  theme( legend.position = "none")+
  scale_colour_manual( values = c( "orange", "darkred", "darkgreen", "darkblue")) + 
  labs( x = "phylogenetic composition \n (weighted unifrac)",
        y = "")

##### OTUdist  ~ EcoLYNdist #####

# exclude among Lake distances and joined 0 distances  
distanceDF_Lake <- distanceDF[ which(apply(distanceDF[,c( "OTUdist", "EcoLYNdist")], 1, sum) > 0),]
distanceDF_Lake <- distanceDF_Lake[distanceDF_Lake$Lake != "full", ]

G1 <- ggplot( distanceDF_Lake, aes( x = OTUdist, y = EcoLYNdist, 
                              colour = Lake, shape = Lake))+
  geom_point( alpha = 0.3)+
  facet_wrap(~Lake)+
  stat_smooth( method = "lm", se = F, size = 1)+
  theme_bw( base_size = 15)+
  theme( legend.position = "none")+
  scale_colour_manual( values = c( "orange", "darkred", "darkgreen", "darkblue")) + 
  labs( y = "Functional distance \n (carbon sources utilization)",
        x = "")

##### OTUdist  ~ NLSdist #####

# exclude among Lake distances and joined 0 distances  
distanceDF_Lake <- distanceDF[ which(apply(distanceDF[,c( "OTUdist", "NLSdist")], 1, sum) > 0),]
distanceDF_Lake <- distanceDF_Lake[distanceDF_Lake$Lake != "full", ]

G3 <- ggplot( distanceDF_Lake, aes( x = OTUdist, y = NLSdist, 
                              colour = Lake, shape = Lake))+
  geom_point( alpha = 0.3)+
  facet_wrap(~Lake)+
  stat_smooth( method = "lm", se = F, size = 1)+
  theme_bw( base_size = 15)+
  theme( legend.position = "none")+
  scale_colour_manual( values = c( "orange", "darkred", "darkgreen", "darkblue")) +
  labs( x = "species composition \n (presence/absence)",
        y = "Functional distance \n (average carbon uptake rate)")


grid.arrange( G1, G2, G3, G4, main=textGrob(
  "community composition vs functional distance",gp=gpar(fontsize=20)))








