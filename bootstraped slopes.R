setwd("~/Documents/01_PhD/01_Research/02_rare_Biosphere/R scripts")

library(ggplot2)
library(plyr)
library(GGally)
library(gridExtra)
library(boot)

ID<-read.table("gbgID",header=TRUE, stringsAsFactors=F)
DIV<-read.table("Diversity_5000_100p.txt",stringsAsFactors=F)
PD <- read.table("PhyloDiv.txt", stringsAsFactors=F)


DIV<-join(ID,DIV)
DIV <- join(DIV,PD)


####### variable selection #########

#check correlation among diversity metrics
#ggpairs(DIV[,c(7:16)])

# Chao, Ace, S, PD, PSR have all pairwise correltations > 0.92 --> only S is kept
cor(DIV[,c(7,8,12,14,16)],use="pairwise.complete.obs")

#check corrleation among remaining metrics
#ggpairs(DIV[,c(9:12,13,15)])

# Hill1 is highly correlated to Hill2 ( r= 0.98), Hill2 is dropped. 
# Pielous J is tightly linked to both HIll1 and HIll2 but the relationship is not linear. Both, J and Hill1 are kept

#check corrleation among remaining metrics
#ggpairs(DIV[,c(9,11,12,13,15)])

# the highest corelation is among Richness and PSV (0.836). PSV is also sensitive to subsampling so it may be dropped later
# for know, HIll1, Pielou's J, Richness, PSV and PSE are kept
DIV  <- DIV[,c(1:6,9,11,12,13,15)]

##### in accordance with the consideration above and the considerations about PSV and J (see manuscript) I only keep PSE, S and Hill1
DIV <- DIV[,-c(8,10)]

#check corrleation among remaining metrics
ggpairs(DIV[,7:9])

# average Diversiyt
avDIV<-ddply(DIV[,-3],.(Lake,DIL), numcolwise(mean,na.rm=T))


############# import Data   #############

############# Biomass #############

## full Dataset ##
BM<-read.table("BMfull.txt",stringsAsFactors=F)

#### maximum Biomass ###

maxBM<-ddply(BM, .(Lake,DIL), function(x) mean(x[with(x, order(-Cells)),][1:5,]$Cells))
colnames(maxBM)<-c("Lake","DIL","Cells")
maxBM <- join(maxBM,avDIV)

### exclude undiluted treatment form data  and exclude NA###
maxBM <- na.omit(maxBM[maxBM$DIL !="0",])

#### avBM at sampling data (sampling date and 2 days prior to sampling date) ###

# subset for sampling dates +-1
avBM <- BM[BM$DAT %in% c(610,612,614,624,626,628,708,710,712),]

# average over 3 days
avBM[avBM$DAT %in% c(610,612,614),]$DAT <- "14Jun"
avBM[avBM$DAT %in% c(624,626,628),]$DAT <- "28Jun"
avBM[avBM$DAT %in% c(708,710,712),]$DAT <- "12Jul"

avBM <- ddply(avBM,.(DAT,Lake,DIL), summarize, Cells = mean(Cells, na.rm=T))

#join Diversity data
avBM <- join(avBM,DIV)

### exclude undiluted treatment form data  and exclude NA###
avBM <- na.omit(avBM[avBM$DIL !="0",])



################# Stability_detr.txt ###############

CV <- read.table("Stability_detr.txt", sep="\t", header=T, stringsAsFactors = F)

CV <- join(avDIV, CV[,c(1,2,13)])

### exclude undiluted treatment form data and exclude NA ###
CV <- na.omit(CV[CV$DIL !="0",])

################# EcoLog.txt ###############

EP <- read.table("EcoPlates.txt", sep="\t", header=T, stringsAsFactors = F)

# join DIV
EP <- join(DIV, EP)

### exclude undiluted treatment form data and exclude NA###
EP <- na.omit(EP[EP$DIL !="0",])

### calculate average EcoLog ###
avEP <- ddply(EP, .(Lake,DIL), summarise, Hill1=mean(Hill1, na.rm=T), S=mean(S, na.rm=T), 
              PSEs = mean(PSEs, na.rm=T), NC = mean(NC, na.rm=T),
              mean.r = mean(mean.r, na.rm=T))


################# Nutrients ################

NUT <- read.table("Nut.txt", sep="\t", header=T, stringsAsFactors=F)
NUT <- join(NUT,ID[ID$DAT=="12Jul",c(3,4,6)])
NUT <- join(avDIV, NUT)

### exclude undiluted treatment form data ###
NUT <- NUT[NUT$DIL !="0",]

############ bootstrap all slopes, p values and R, squares for all 5 diversity metrics and 4 Lakes #############


# function to extract slope (of scaled regression), rsquare and p-value
# Y is response varibale, V predictor V

extract_Slope <- function(x,Y,V) {
    LM <- lm(scale(x[,Y]) ~ scale(x[,V]), data=x)
    c(coef(LM)[[2]])
}

# function to bootstrap with
boot_Slope <- function(d,i,Y,V) {
  extract_Slope(na.omit(d[i,]),Y,V)
}


########### parameter for bootstrapping #########


Lakes <- c("Botan", "Delsjön","Lerum","Surtesjön")
DIV_ind  <- c("Hill1","S","PSEs")
DAT  <- c( "14Jun", "28Jun", "12Jul")

# number of bootstrappes
R <- 10000

################ maximum Biomass ###########
Start  <- Sys.time()
# Data frames to store Results in

SlopesNH4 <- SlopesN23 <- SlopesEPmr <- SlopesEPNC <- 
  SlopesCV  <- Slopes_maxBM <- data.frame(Lake = rep(Lakes, each=3), DIV = rep(DIV_ind,4),
                                          medSlope = NA, meanSlope = NA, minSl = NA, maxSl = NA)
# Slopes_maxBM

for (L in Lakes) {
  for (V in DIV_ind) {
    BOOT <- boot(data= maxBM[maxBM$Lake == L,] , boot_Slope , R=R, V=V, Y="Cells")
    Slopes_maxBM [Slopes_maxBM$Lake == L & Slopes_maxBM$DIV == V ,]$medSlope <- median(BOOT$t[,1]) # median slope
    Slopes_maxBM [Slopes_maxBM$Lake == L & Slopes_maxBM$DIV == V ,]$meanSlope <- mean(BOOT$t[,1]) # mean slope
    Slopes_maxBM [Slopes_maxBM$Lake == L & Slopes_maxBM$DIV == V ,]$minSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
    Slopes_maxBM [Slopes_maxBM$Lake == L & Slopes_maxBM$DIV == V ,]$maxSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # max slope
  }
}

################ average Biomass at sampling date ###########

Slopes_avEPmr <- Slopes_avEPNC <- Slopes_avBM  <-  data.frame(Lake = rep(Lakes, each=9), DAT = rep(rep(DAT,each=3),4),
                                                           DIV = rep(DIV_ind,12), medSlope = NA, meanSlope = NA,
                                                           minSl = NA, maxSl = NA)

for (L in Lakes) {
  for (D in DAT) {
    for (V in DIV_ind) {
      BOOT <- boot(data= avBM[avBM$Lake == L & avBM$DAT == D,] , boot_Slope , R=R, V=V, Y="Cells")
      Slopes_avBM [Slopes_avBM$Lake == L & Slopes_avBM$DAT == D & Slopes_avBM$DIV == V ,]$medSlope <- median(BOOT$t[,1]) # median slope
      Slopes_avBM [Slopes_avBM$Lake == L & Slopes_avBM$DAT == D & Slopes_avBM$DIV == V ,]$meanSlope <- mean(BOOT$t[,1]) # mean slope
      Slopes_avBM [Slopes_avBM$Lake == L & Slopes_avBM$DAT == D & Slopes_avBM$DIV == V ,]$minSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
      Slopes_avBM [Slopes_avBM$Lake == L & Slopes_avBM$DAT == D & Slopes_avBM$DIV == V ,]$maxSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # av slope
    }
  }
}


################ Stability - CV ###########

# SlopesCV

for (L in Lakes) {
  for (V in DIV_ind) {
    BOOT <- boot(data= CV[CV$Lake == L,] , boot_Slope , R=R, V=V, Y="Stability")
    SlopesCV [SlopesCV$Lake == L & SlopesCV$DIV == V ,]$medSlope <- median(BOOT$t[,1]) # median slope
    SlopesCV [SlopesCV$Lake == L & SlopesCV$DIV == V ,]$meanSlope <- mean(BOOT$t[,1]) # mean slope
    SlopesCV [SlopesCV$Lake == L & SlopesCV$DIV == V ,]$minSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
    SlopesCV [SlopesCV$Lake == L & SlopesCV$DIV == V ,]$maxSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # max slope
  }
}

################# EcoLog.txt ###############

########## average NC ##########

#SlopesEPNC

for (L in Lakes) {
  for (V in DIV_ind) {
    BOOT <- boot(data= avEP[avEP$Lake == L,] , boot_Slope , R=R, V=V, Y="NC")
    SlopesEPNC [SlopesEPNC$Lake == L & SlopesEPNC$DIV == V ,]$medSlope <- median(BOOT$t[,1]) # median slope
    SlopesEPNC [SlopesEPNC$Lake == L & SlopesEPNC$DIV == V ,]$meanSlope <- mean(BOOT$t[,1]) # mean slope
    SlopesEPNC [SlopesEPNC$Lake == L & SlopesEPNC$DIV == V ,]$minSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
    SlopesEPNC [SlopesEPNC$Lake == L & SlopesEPNC$DIV == V ,]$maxSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # max slope
  }
}

########## Sampling NC ##########

#Slopes_avEPNC

for (L in Lakes) {
  for (D in DAT) {
    for (V in DIV_ind) {
      BOOT <- boot(data= EP[EP$Lake == L & EP$DAT == D,] , boot_Slope , R=R, V=V, Y="NC")
      Slopes_avEPNC [Slopes_avEPNC$Lake == L & Slopes_avEPNC$DAT == D & Slopes_avEPNC$DIV == V ,]$medSlope <- median(BOOT$t[,1]) # median slope
      Slopes_avEPNC [Slopes_avEPNC$Lake == L & Slopes_avEPNC$DAT == D & Slopes_avEPNC$DIV == V ,]$meanSlope <- mean(BOOT$t[,1]) # mean slope
      Slopes_avEPNC [Slopes_avEPNC$Lake == L & Slopes_avEPNC$DAT == D & Slopes_avEPNC$DIV == V ,]$minSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
      Slopes_avEPNC [Slopes_avEPNC$Lake == L & Slopes_avEPNC$DAT == D & Slopes_avEPNC$DIV == V ,]$maxSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # av slope
    }
  }
}


########## average uptake rate ##########

#SlopesEPmr

for (L in Lakes) {
  for (V in DIV_ind) {
    BOOT <- boot(data= avEP[avEP$Lake == L,] , boot_Slope , R=R, V=V, Y="mean.r")
    SlopesEPmr [SlopesEPmr$Lake == L & SlopesEPmr$DIV == V ,]$medSlope <- median(BOOT$t[,1]) # median slope
    SlopesEPmr [SlopesEPmr$Lake == L & SlopesEPmr$DIV == V ,]$meanSlope <- mean(BOOT$t[,1]) # mean slope
    SlopesEPmr [SlopesEPmr$Lake == L & SlopesEPmr$DIV == V ,]$minSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
    SlopesEPmr [SlopesEPmr$Lake == L & SlopesEPmr$DIV == V ,]$maxSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # max slope
  }
}

########## sampling day uptake rate ##########

#Slopes_avEPmr

for (L in Lakes) {
  for (D in DAT) {
    for (V in DIV_ind) {
      BOOT <- boot(data= EP[EP$Lake == L & EP$DAT == D,] , boot_Slope , R=R, V=V, Y="mean.r")
      Slopes_avEPmr [Slopes_avEPmr$Lake == L & Slopes_avEPmr$DAT == D & Slopes_avEPmr$DIV == V ,]$medSlope <- median(BOOT$t[,1]) # median slope
      Slopes_avEPmr [Slopes_avEPmr$Lake == L & Slopes_avEPmr$DAT == D & Slopes_avEPmr$DIV == V ,]$meanSlope <- mean(BOOT$t[,1]) # mean slope
      Slopes_avEPmr [Slopes_avEPmr$Lake == L & Slopes_avEPmr$DAT == D & Slopes_avEPmr$DIV == V ,]$minSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
      Slopes_avEPmr [Slopes_avEPmr$Lake == L & Slopes_avEPmr$DAT == D & Slopes_avEPmr$DIV == V ,]$maxSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # av slope
    }
  }
}

################# Nutrients ################


################# NO23###########

#SlopesN23

for (L in Lakes) {
  for (V in DIV_ind) {
    BOOT <- boot(data= NUT[NUT$Lake == L,] , boot_Slope , R=R, V=V, Y="NO23")
    SlopesN23 [SlopesN23$Lake == L & SlopesN23$DIV == V ,]$medSlope <- median(BOOT$t[,1]) # median slope
    SlopesN23 [SlopesN23$Lake == L & SlopesN23$DIV == V ,]$meanSlope <- mean(BOOT$t[,1]) # mean slope
    SlopesN23 [SlopesN23$Lake == L & SlopesN23$DIV == V ,]$minSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
    SlopesN23 [SlopesN23$Lake == L & SlopesN23$DIV == V ,]$maxSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # max slope
  }
}


################# NH4 ###########


#SlopesNH4

for (L in Lakes) {
  for (V in DIV_ind) {
    BOOT <- boot(data= NUT[NUT$Lake == L,] , boot_Slope , R=R, V=V, Y="NH4")
    SlopesNH4 [SlopesNH4$Lake == L & SlopesNH4$DIV == V ,]$medSlope <- median(BOOT$t[,1]) # median slope
    SlopesNH4 [SlopesNH4$Lake == L & SlopesNH4$DIV == V ,]$meanSlope <- mean(BOOT$t[,1]) # mean slope
    SlopesNH4 [SlopesNH4$Lake == L & SlopesNH4$DIV == V ,]$minSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
    SlopesNH4 [SlopesNH4$Lake == L & SlopesNH4$DIV == V ,]$maxSl <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # max slope
  }
}


End  <- Sys.time()
End-Start

##############################################################################################################################

G1 <- ggplot(Slopes_maxBM , aes(x=medSlope, y=DIV, colour=Lake, xmin=minSl, xmax=maxSl))+
  geom_point(size=3)+
  facet_wrap(~Lake,nrow=1)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1),breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="maximum Biomass",x="Slope")

G1_2 <- ggplot(Slopes_avBM , aes(x=medSlope, y=DIV, colour=Lake, xmin=minSl, xmax=maxSl))+
  geom_point(size=3)+
  facet_wrap(~Lake*DAT,nrow=4)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  # scale_x_continuous(limits = c(-1, 1),breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="maximum Biomass",x="Slope")

G2 <- ggplot(SlopesEPNC, aes(x=medSlope, y=DIV, colour=Lake, xmin=minSl, xmax=maxSl))+
  geom_point(size=3)+
  facet_wrap(~Lake,nrow=1)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1), breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="Number of C sources", x="Slope")

G2_2 <- ggplot(Slopes_avEPNC , aes(x=medSlope, y=DIV, colour=Lake, xmin=minSl, xmax=maxSl))+
  geom_point(size=3)+
  facet_wrap(~Lake*DAT,nrow=4)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  # scale_x_continuous(limits = c(-1, 1),breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="Number of C sources",x="Slope")

G3 <- ggplot(SlopesEPmr, aes(x=medSlope, y=DIV, colour=Lake, xmin=minSl, xmax=maxSl))+
  geom_point(size=3)+
  facet_wrap(~Lake,nrow=1)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1), breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="mean carbon uptake rate", x="Slope")

G3_2 <- ggplot(Slopes_avEPmr , aes(x=medSlope, y=DIV, colour=Lake, xmin=minSl, xmax=maxSl))+
  geom_point(size=3)+
  facet_wrap(~Lake*DAT,nrow=4)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  # scale_x_continuous(limits = c(-1, 1),breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="mean carbon uptake rate",x="Slope")


G4 <- ggplot(SlopesCV, aes(x=medSlope, y=DIV, colour=Lake, xmin=minSl, xmax=maxSl))+
  geom_point(size=3)+
  facet_wrap(~Lake,nrow=1)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1),breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="Stability",x="Slope")


G5 <- ggplot(SlopesNH4, aes(x=medSlope, y=DIV, colour=Lake, xmin=minSl, xmax=maxSl))+
  geom_point(size=3)+
  facet_wrap(~Lake,nrow=1)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1), breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="nutrient depletion (NH4)", x="Slope")

G6 <- ggplot(SlopesN23, aes(x=medSlope, y=DIV, colour=Lake, xmin=minSl, xmax=maxSl))+
  geom_point(size=3)+
  facet_wrap(~Lake,nrow=1)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1), breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="nutrient depletion (NO2 + NO3)", x="Slope")

grid.arrange(G1,G2,G3,G4,G5,G6)
grid.arrange(G1_2,G2_2,G3_2)

#save.image("bootstrap_10000.RData")





######################################################################################################################################
######## correlation of slopes ############
RES <- c("Slopes","SlopesCV","SlopesEPmr","SlopesEPNC",
         "SlopesN23","SlopesNH4")

Slopes_combined <- Slopes[0,]
Slopes_combined$EF <- character()

for (i in 1:length(RES)) {
  Temp <- get(RES[i])
  Temp$EF <- RES[i]
  Slopes_combined <- rbind(Slopes_combined,Temp)
  
}

ggplot(Slopes_combined, aes(x=DIV,y=meanSlope,colour=DIV))+
  geom_boxplot()

L_EF <- by(Slopes_combined,Slopes_combined$EF,list)

L_EF_pairs <- data.frame(maxBM=L_EF[[1]]$meanSlope)

res <- c("stability","mean_r","N_C","Nut_N23","Nut_N4")
for (i in 2:6) {
  L_EF_pairs$t <- L_EF[[i]]$meanSlope
  colnames(L_EF_pairs)[i] <- res[i-1]
}

ggpairs(L_EF_pairs)
cor(L_EF_pairs, method="sp")

########
L_EF <- by(Slopes_combined,Slopes_combined$Lake,list)

L_EF_pairs <- data.frame(Botan=L_EF[[1]]$meanSlope)

res <- c("Delsjön","Lerum","Surtesjön")

for (i in 2:4) {
  L_EF_pairs$t <- L_EF[[i]]$meanSlope
  colnames(L_EF_pairs)[i] <- res[i-1]
}

ggpairs(L_EF_pairs)
cor(L_EF_pairs, method="sp")

########
L_EF <- by(Slopes_combined,Slopes_combined$DIV,list)

L_EF_pairs <- data.frame(Hill1=L_EF[[1]]$meanSlope)

res <- c("S","PSE")

for (i in 2:3) {
  L_EF_pairs$t <- L_EF[[i]]$meanSlope
  colnames(L_EF_pairs)[i] <- res[i-1]
}

ggpairs(L_EF_pairs)
cor(L_EF_pairs, method="sp")


################################################################
##### new clustering with OTU data and distance matrix ########
####### McMurdie, Paul J, and Susan Holmes. 2014. “Waste Not, Want Not: Why Rarefying Microbiome Data Is Inadmissible.”
####### PLoS Computational Biology

library("edgeR")

OTU <- read.table("OTU_woS_onlyBAC.txt")
OTUm <- as.matrix(OTU[,-1])
dimnames(OTUm) <- list(OTU[,1],colnames(OTU[,-1]))

# group by Lakes
OTU.L <- DGEList(OTUm, group=ID[match(colnames(OTUm),ID$gbgID),]$Lake)

# The method used in the edgeR vignette is to keep only those genes (OTUs) 
# that have at least 1 read per million in at least 3 samples

OTU.L <- OTU.L[rowSums(1e+06 * OTU.L$counts/expandAsMatrix(OTU.L$samples$lib.size, 
                                                           dim(OTU.L)) > 1) >= 3, ]
dim(OTU.L) # 263 OTUs excluded from dataset

OTU.L <- calcNormFactors( OTU.L )

plotMDS(OTU.L, col=as.numeric(as.factor(
  ID[match(colnames(OTU.L$counts),ID$gbgID),]$Lake))+1,
  labels = ID[match(colnames(OTU.L$counts),ID$gbgID),]$DIL,
  method ="bcv")

# estimate data.set wide dispersion, than estimate sample dispersion

OTU.L <- estimateCommonDisp( OTU.L )
OTU.L <- estimateTagwiseDisp( OTU.L , prior.df = 36 )


meanVarPlot <- plotMeanVar( OTU.L , show.raw.vars=TRUE ,
                            show.tagwise.vars=TRUE ,
                            show.binned.common.disp.vars=FALSE ,
                            show.ave.raw.vars=FALSE ,
                            dispersion.method = "qcml" , NBline = TRUE ,
                            nbins = 100 ,
                            pch = 16 ,
                            xlab ="Mean Expression (Log10 Scale)" ,
                            ylab = "Variance (Log10 Scale)" ,
                            main = "Mean-Variance Plot" )

exT <- exactTest( OTU.L, pair = c(2,4) )
topTags( exT , n = 20)


