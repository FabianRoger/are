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
library(vegan)
###############

# import files
###############
#weighted unifrac distance based on rarefied OTU table
load("unifrac_dist.RData")

#weighted jaccard distance based on rarefied OTU table
load("OTUdist.RData")

#euclidian distance Ecolog carbon source utilization (Y/N)
load("NLSdist.RData")

#euclidian distance Ecolog carbon uptake rate
load("EcoLYNdist.RData")

# Importing metadata
ID<-read.table("gbgID",header=TRUE, stringsAsFactors=F)

# change factor level for Dates
ID$DAT <- factor( ID$DAT, levels = c( "14Jun", "28Jun", "12Jul"))

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



############## calculate Mantel test ###########################################

# function to calculate mantel test between tow matrices in a give Lake

MantelLakes <- function(mat1,mat2,Lake,method,perm) {
  LakeID <- ID[ID$Lake == Lake,]$gbgID
  W1<- which(labels(mat1) %in% LakeID)
  W2 <- which( labels( mat2) %in% LakeID)
  dist1 <- as.dist( as.matrix( mat1)[ W1, W1])
  dist2 <- as.dist( as.matrix( mat2 )[ W2, W2])
  RES <- mantel(dist1,dist2,method=method, permutations = perm)
  return(RES)
}



Comp <- combn(names(dist.list),2)[,2:5]

RES <- data.frame(comp =  rep(paste(Comp[ 1, ], Comp[ 2, ], sep = "_" ),each = 4),
                  Lake = rep( levels( as.factor( ID$Lake)), 4), r = NA, p = NA,
                  r2 = NA )



for ( comp in 1:4 ) {
  for (L in levels( as.factor( ID$Lake))) {
    Mant <-  MantelLakes(mat1 = dist.list[[which(names(dist.list) == Comp[1,comp])]],
                         mat2 = dist.list[[which(names(dist.list) == Comp[2,comp])]],
                         Lake = L,
                         method = "pear",
                         perm = 9999)
    comP <- paste(Comp[1,comp],Comp[2,comp],sep="_") 
    print(comP)
    RES[RES$Lake == L & RES$comp == comP,]$r <- Mant$statistic
    RES[RES$Lake == L & RES$comp == comP,]$r2<- Mant$statistic^2
    RES[RES$Lake == L & RES$comp == comP,]$p <- Mant$signif
    
  }
}

RES[,3:5] <- apply(RES[,3:5], 2, signif, digits = 2 )

RES$comp <- factor(RES$comp, labels = c(levels(as.factor(RES$comp))[c(1,3,2,4)]))

ggplot(RES, aes(x=Lake, y=r2, fill=Lake, label = p))+
  geom_bar(stat="identity")+
  geom_text(aes(y = r2+0.01))+
  facet_wrap(~comp)+
  theme_bw(base_size=15)+
  scale_fill_manual(values=c("orange","darkred","darkgreen","darkblue"))





############## plot distance scatterplots ######################################

# create dataframe for plotting rectangles

DF.Rect <- data.frame( xmin = 0.3, xmax = 0.5, ymin = 0.75, ymax = 0.95)

##### OTUdist  ~ EcoLYNdist #####

# exclude among Lake distances and joined 0 distances  
distanceDF_Lake <- distanceDF[ which(apply(distanceDF[,c( "OTUdist", "EcoLYNdist")], 1, sum) > 0),]
distanceDF_Lake <- distanceDF_Lake[distanceDF_Lake$Lake != "full", ]

G1 <- ggplot( distanceDF_Lake, aes( x = OTUdist, y = EcoLYNdist))+
  geom_point( aes(alpha = 0.3,  colour = Lake, shape = Lake))+
  geom_rect( data = DF.Rect, aes(x=NULL,y=NULL, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
             fill = "white", alpha = 0.8) +
  geom_text(data = RES[13:16, ], aes( x = 0.4, y = 0.9, label = paste("r^2 : ", r2, sep = "")),
            parse = TRUE , colour = "black", size = 3) +
  geom_text(data = RES[13:16, ], aes( x = 0.4, y = 0.8, label = paste("p : ", p, sep = "")),
            parse = TRUE , colour = "black", size = 3) +
  facet_wrap(~Lake)+
  stat_smooth( method = "lm", se = F, size = 1, colour = "black")+
  theme_bw( base_size = 15)+
  theme( legend.position = "none")+
  scale_colour_manual( values = c( "orange", "darkred", "darkgreen", "darkblue")) +
  labs( x="", y= "functional distance\n(carbon source utilization)" )

##### unifrac_dist  ~ EcoLYNdist #####

# exclude among Lake distances and joined 0 distances  
distanceDF_Lake <- distanceDF[ which(apply(distanceDF[,c( "unifrac_dist", "EcoLYNdist")], 1, sum) > 0),]
distanceDF_Lake <- distanceDF_Lake[distanceDF_Lake$Lake != "full", ]

G2 <- ggplot( distanceDF_Lake, aes( x = unifrac_dist, y = EcoLYNdist ))+
  geom_point( aes(alpha = 0.3,  colour = Lake, shape = Lake))+
  geom_rect( data = DF.Rect, aes(x=NULL,y=NULL, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
             fill = "white", alpha = 0.8) +
  geom_text(data = RES[9:12, ], aes( x = 0.4, y = 0.9, label = paste("r^2 : ", r2, sep = "")),
            parse = TRUE , colour = "black", size = 3) +
  geom_text(data = RES[9:12, ], aes( x = 0.4, y = 0.8, label = paste("p : ", p, sep = "")),
            parse = TRUE , colour = "black", size = 3) +
  facet_wrap(~Lake)+
  stat_smooth( method = "lm", se = F, size = 1, colour = "black")+
  theme_bw( base_size = 15)+
  theme( legend.position = "none")+
  scale_colour_manual( values = c( "orange", "darkred", "darkgreen", "darkblue")) +
  labs( x="", y= "")

##### OTUdist  ~ NLSdist #####

# exclude among Lake distances and joined 0 distances  
distanceDF_Lake <- distanceDF[ which(apply(distanceDF[,c( "OTUdist", "NLSdist")], 1, sum) > 0),]
distanceDF_Lake <- distanceDF_Lake[distanceDF_Lake$Lake != "full", ]


G3 <- ggplot( distanceDF_Lake, aes( x = OTUdist, y = NLSdist))+
  geom_point( aes(alpha = 0.3,  colour = Lake, shape = Lake))+ 
  geom_rect( data = DF.Rect, aes(x=NULL,y=NULL, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
             fill = "white", alpha = 0.8) +
  geom_text(data = RES[5:8, ], aes( x = 0.4, y = 0.9, label = paste("r^2 : ", r2, sep = "")),
            parse = TRUE , colour = "black", size = 3) +
  geom_text(data = RES[5:8, ], aes( x = 0.4, y = 0.8, label = paste("p : ", p, sep = "")),
            parse = TRUE , colour = "black", size = 3) +
  facet_wrap(~Lake)+
  stat_smooth( method = "lm", se = F, size = 1, colour = "black")+
  theme_bw( base_size = 15)+
  theme( legend.position = "none")+
  scale_colour_manual( values = c( "orange", "darkred", "darkgreen", "darkblue")) +
  labs( x="species composition\n(presence/absence)", y= "functional distance\n(average carbon uptake rate)")


##### unifrac_dist  ~ NLSdist #####

# exclude among Lake distances and joined 0 distances  
distanceDF_Lake <- distanceDF[ which(apply(distanceDF[,c( "unifrac_dist", "NLSdist")], 1, sum) > 0),]
distanceDF_Lake <- distanceDF_Lake[distanceDF_Lake$Lake != "full", ]

G4 <- ggplot( distanceDF_Lake, aes( x = unifrac_dist, y = NLSdist, )) +
  geom_point( aes(alpha = 0.3,  colour = Lake, shape = Lake))+
  geom_rect( data = DF.Rect, aes(x=NULL,y=NULL, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
             fill = "white", alpha = 0.8) +
  geom_text(data = RES[1:4, ], aes( x = 0.4, y = 0.9, label = paste("r^2 : ", r2, sep = "")),
            parse = TRUE , colour = "black", size = 3) +
  geom_text(data = RES[1:4, ], aes( x = 0.4, y = 0.8, label = paste("p : ", p, sep = "")),
            parse = TRUE , colour = "black", size = 3) +
  facet_wrap(~Lake)+
  stat_smooth( method = "lm", se = F, size = 1, colour = "black")+
  theme_bw( base_size = 15)+
  theme( legend.position = "none")+
  scale_colour_manual( values = c( "orange", "darkred", "darkgreen", "darkblue")) +
  labs( x="phylogenetic composition\n(weighted unifrac)", y= "")


Figure_5 <- arrangeGrob( G1, G2, G3, G4, main=textGrob(
  "community composition vs functional distance",gp=gpar(fontsize=20)))

ggsave( filename = "Figure_5.pdf", plot = Figure_5, width = 10, height = 10)


############## plot NMDS plots #################################################


########## unifrac ##########

fit <- metaMDS(unifrac_dist, k=2) # k is the number of dim

# extract points for plotting
fitp<-data.frame(fit$points)

#join metadata
fitp$gbgID<-rownames(fitp)

fitp<-join(fitp,ID)

# for 2d plot ()

find_hull <- function(df) df[chull(df$MDS1, df$MDS2), ]
hulls <- ddply(fitp[,c("MDS1","MDS2","Lake","DAT")], .(Lake,DAT), find_hull)

G_nmds_1 <- ggplot(fitp,aes(x=MDS1,y=MDS2))+
  geom_point(data=fitp,aes(colour=Lake,shape=Lake,size=4))+
  geom_text(data=fitp,aes(colour=Lake,label=DIL,hjust=-0.7,size=4))+
  geom_polygon(data=hulls,aes(alpha=0.8,fill=Lake))+
  facet_wrap(~DAT)+
  theme_bw(base_size=15)+
  scale_colour_manual(values=c("orange","darkred","darkgreen","darkblue"))+
  scale_fill_manual(values=c("orange","darkred","darkgreen","darkblue"))+
  theme(legend.position="none")+
  guides(size=F,alpha=F)+
  labs(title=paste("nmds of OTU data \n weighted unifrac; stress = ",signif(fit$stress,2)))


######### OTUdist  #########

fit <- metaMDS(OTUdist, k=2) # k is the number of dim

# extract points for plotting
fitp<-data.frame(fit$points)

# flipping y axis
fitp$MDS2 <- -1 * fitp$MDS2

#join metadata
fitp$gbgID<-rownames(fitp)

fitp<-join(fitp,ID)

# for 2d plot ()

find_hull <- function(df) df[chull(df$MDS1, df$MDS2), ]
hulls <- ddply(fitp[,c("MDS1","MDS2","Lake","DAT")], .(Lake,DAT), find_hull)

G_nmds_2 <- ggplot(fitp,aes(x=MDS1,y=MDS2))+
  geom_point(data=fitp,aes(colour=Lake,shape=Lake,size=4))+
  geom_text(data=fitp,aes(colour=Lake,label=DIL,hjust=-0.7,size=4))+
  geom_polygon(data=hulls,aes(alpha=0.8,fill=Lake))+
  facet_wrap(~DAT)+
  theme_bw(base_size=15)+
  scale_colour_manual(values=c("orange","darkred","darkgreen","darkblue"))+
  scale_fill_manual(values=c("orange","darkred","darkgreen","darkblue"))+
  theme(legend.position="none")+
  guides(size=F,alpha=F)+
  labs(title=paste("nmds of OTU data \n bray-curtis presence/absence; stress = ",signif(fit$stress,2)))


######### EcoLYNdist  #########

fit <- metaMDS(EcoLYNdist, k=2) # k is the number of dim

# extract points for plotting
fitp<-data.frame(fit$points)

#join metadata
fitp$gbgID<-rownames(fitp)

fitp<-join(fitp,ID)

# for 2d plot ()

find_hull <- function(df) df[chull(df$MDS1, df$MDS2), ]
hulls <- ddply(fitp[,c("MDS1","MDS2","Lake","DAT")], .(Lake,DAT), find_hull)

G_nmds_3 <- ggplot(fitp,aes(x=MDS1,y=MDS2))+
  geom_point(data=fitp,aes(colour=Lake,shape=Lake,size=4))+
  geom_text(data=fitp,aes(colour=Lake,label=DIL,hjust=-0.7,size=4))+
  geom_polygon(data=hulls,aes(alpha=0.8,fill=Lake))+
  facet_wrap(~DAT)+
  theme_bw(base_size=15)+
  scale_colour_manual(values=c("orange","darkred","darkgreen","darkblue"))+
  scale_fill_manual(values=c("orange","darkred","darkgreen","darkblue"))+
  theme(legend.position="none")+
  guides(size=F,alpha=F)+
  labs(title=paste("nmds of carbon source utilization pattern \n euclidian; stress = ",signif(fit$stress,2)))


######### NLSdist  #########

fit <- metaMDS(NLSdist, k=2) # k is the number of dim

# extract points for plotting
fitp<-data.frame(fit$points)

#join metadata
fitp$gbgID<-rownames(fitp)

fitp<-join(fitp,ID)

# for 2d plot ()

find_hull <- function(df) df[chull(df$MDS1, df$MDS2), ]
hulls <- ddply(fitp[,c("MDS1","MDS2","Lake","DAT")], .(Lake,DAT), find_hull)

G_nmds_4 <- ggplot(fitp,aes(x=MDS1,y=MDS2))+
  geom_point(data=fitp,aes(colour=Lake,shape=Lake,size=4))+
  geom_text(data=fitp,aes(colour=Lake,label=DIL,hjust=-0.7,size=4))+
  geom_polygon(data=hulls,aes(alpha=0.8,fill=Lake))+
  facet_wrap(~DAT)+
  theme_bw(base_size=15)+
  scale_colour_manual(values=c("orange","darkred","darkgreen","darkblue"))+
  scale_fill_manual(values=c("orange","darkred","darkgreen","darkblue"))+
  theme(legend.position="none")+
  guides(size=F,alpha=F)+
  labs(title=paste("nmds of carbon source uptake rates \n euclidian, stress = ",signif(fit$stress,2)))


grid.arrange(G_nmds_1,G_nmds_2,G_nmds_3,G_nmds_4)

########################### heatmaps ###########################################
setwd("~/Documents/01_PhD/01_Research/02_rare_Biosphere/R scripts")

EcoLYN<-as.matrix(read.table("EcoLYN.txt",sep="\t"))
NLSr<-as.matrix(read.table("NLSr.txt",sep="\t"))

EcoLYN <- EcoLYN*100

ID$DAT <- factor(ID$DAT, levels = c("14Jun", "28Jun", "12Jul"))
ID$DIL <- factor(ID$DIL, levels = c("0","1","2","3","4","5","6","7","8","9","10","S"))
rownames(ID) <- ID$barcode

phylo_YN <- phyloseq( otu_table( NLSr, taxa_are_rows = F), 
                      sample_data(ID))

# sample order
IDphy <- ID[ID$barcode %in% sample_names(phylo_YN),]
sample.order <- as.character(IDphy[with(IDphy, order(Lake,DIL,DAT)),]$barcode)

NLSrplot <- plot_heatmap(phylo_YN, sample.order = sample.order, 
                         sample.label="Lake")

NLSrplot +
  facet_wrap(DAT ~ Lake, scales = "free_x")+
  theme_bw(base_size=15)



