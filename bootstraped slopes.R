setwd("~/Documents/01_PhD/01_Research/02_rare_Biosphere/R scripts")

library(ggplot2)
library(plyr)
library(GGally)
library(gridExtra)
library(boot)
library(reshape2)

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
ggpairs(DIV[,7:9], diag = list(continuous = "bar"))

# average Diversiyt
avDIV<-ddply(DIV[,-3],.(Lake,DIL), numcolwise(mean,na.rm=T))

############ plot diversity gradient ##########################################

#reorder factor level
DIV$DAT <- factor( DIV$DAT, levels = c( "14Jun", "28Jun", "12Jul"))
DIV$DIL <- factor( DIV$DIL, levels = c( "0", "1", "2", "3", "4", "5",
                                        "6", "7", "8", "9", "10", "S"))
# Hill1
G_Hill <- ggplot(DIV, aes( x = DIL, y = Hill1, colour = Lake, shape = Lake, group = 1))+
  geom_point()+
  facet_wrap( ~DAT * Lake)+
  theme_bw(base_size = 8)+
  stat_smooth( method = "lm", se = F, colour = "grey", alpha = 0.5)+
  scale_y_log10(breaks=c(2,4,8,16,32,64,128,256,512))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs( x = "dillution factor", y = "Hill1" , 
        title = "effective number of species (Hill1) ")+
  theme( legend.position = "none")

# S
G_S <- ggplot(DIV, aes( x = DIL, y = S, colour = Lake, shape = Lake, group = 1))+
  geom_point()+
  facet_wrap( ~DAT * Lake)+
  theme_bw(base_size = 8)+
  stat_smooth( method = "lm", se = F, colour = "grey", alpha = 0.5)+
  scale_y_log10(breaks=c(2,4,8,16,32,64,128,256,512))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs( x = "dillution factor", y = "S", 
        title = "OTU richness (S) ")+
  theme( legend.position = "none")

# PSE
G_PSE <- ggplot(DIV, aes( x = DIL, y = PSEs, colour = Lake, shape = Lake, group = 1))+
  geom_point()+
  facet_wrap( ~DAT * Lake)+
  theme_bw(base_size = 8)+
  stat_smooth( method = "lm", se = F, colour = "grey", alpha = 0.5)+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs( x = "dillution factor", y = "PSE", 
        title = "phylogenetic species eveness (PSE)")+
  theme( legend.position = "none")

G_div <- arrangeGrob(G_S, G_Hill, G_PSE, nrow = 2)
ggsave( filename = "Figure_1.pdf", plot = G_div, width = 8, height = 8)


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

# calculate total dissolven inorganic nitrogen (DIN)
NUT$DIN <- NUT$NO23 + NUT$NH4

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

SlopesDIN <- SlopesEPmr <- SlopesEPNC <- 
  SlopesCV  <- Slopes_maxBM <- data.frame(
    Lake = rep(Lakes, each=3), DIV = rep(DIV_ind,4),medSlope = NA,
    meanSlope = NA, minSl99 = NA, maxSl99 = NA, minSl95 = NA, maxSl95 = NA,
    minSl90 = NA, maxSl90 = NA)

# Slopes_maxBM

for (L in Lakes) {
  for (V in DIV_ind) {
    BOOT <- boot(data= maxBM[maxBM$Lake == L,] , boot_Slope , R=R, V=V, Y="Cells")
    Slopes_maxBM [Slopes_maxBM$Lake == L & Slopes_maxBM$DIV == V ,]$medSlope <- median(BOOT$t[,1]) # median slope
    Slopes_maxBM [Slopes_maxBM$Lake == L & Slopes_maxBM$DIV == V ,]$meanSlope <- mean(BOOT$t[,1]) # mean slope
    Slopes_maxBM [Slopes_maxBM$Lake == L & Slopes_maxBM$DIV == V ,]$minSl99 <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
    Slopes_maxBM [Slopes_maxBM$Lake == L & Slopes_maxBM$DIV == V ,]$maxSl99 <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # max slope
    Slopes_maxBM [Slopes_maxBM$Lake == L & Slopes_maxBM$DIV == V ,]$minSl95 <- boot.ci(BOOT, conf=0.95, type="bca", index=1)$bca[[4]] # min slope
    Slopes_maxBM [Slopes_maxBM$Lake == L & Slopes_maxBM$DIV == V ,]$maxSl95 <- boot.ci(BOOT, conf=0.95, type="bca", index=1)$bca[[5]] # max slope
    Slopes_maxBM [Slopes_maxBM$Lake == L & Slopes_maxBM$DIV == V ,]$minSl90 <- boot.ci(BOOT, conf=0.90, type="bca", index=1)$bca[[4]] # min slope
    Slopes_maxBM [Slopes_maxBM$Lake == L & Slopes_maxBM$DIV == V ,]$maxSl90 <- boot.ci(BOOT, conf=0.90, type="bca", index=1)$bca[[5]] # max slope
  }
}


################ Stability - CV ###########

# SlopesCV

for (L in Lakes) {
  for (V in DIV_ind) {
    BOOT <- boot(data= CV[CV$Lake == L,] , boot_Slope , R=R, V=V, Y="Stability")
    SlopesCV [SlopesCV$Lake == L & SlopesCV$DIV == V ,]$medSlope <- median(BOOT$t[,1]) # median slope
    SlopesCV [SlopesCV$Lake == L & SlopesCV$DIV == V ,]$meanSlope <- mean(BOOT$t[,1]) # mean slope
    SlopesCV [SlopesCV$Lake == L & SlopesCV$DIV == V ,]$minSl99 <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
    SlopesCV [SlopesCV$Lake == L & SlopesCV$DIV == V ,]$maxSl99 <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # max slope
    SlopesCV [SlopesCV$Lake == L & SlopesCV$DIV == V ,]$minSl95 <- boot.ci(BOOT, conf=0.95, type="bca", index=1)$bca[[4]] # min slope
    SlopesCV [SlopesCV$Lake == L & SlopesCV$DIV == V ,]$maxSl95 <- boot.ci(BOOT, conf=0.95, type="bca", index=1)$bca[[5]] # max slope
    SlopesCV [SlopesCV$Lake == L & SlopesCV$DIV == V ,]$minSl90 <- boot.ci(BOOT, conf=0.90, type="bca", index=1)$bca[[4]] # min slope
    SlopesCV [SlopesCV$Lake == L & SlopesCV$DIV == V ,]$maxSl90 <- boot.ci(BOOT, conf=0.90, type="bca", index=1)$bca[[5]] # max slope
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
    SlopesEPNC [SlopesEPNC$Lake == L & SlopesEPNC$DIV == V ,]$minSl99 <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
    SlopesEPNC [SlopesEPNC$Lake == L & SlopesEPNC$DIV == V ,]$maxSl99 <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # max slope
    SlopesEPNC [SlopesEPNC$Lake == L & SlopesEPNC$DIV == V ,]$minSl95 <- boot.ci(BOOT, conf=0.95, type="bca", index=1)$bca[[4]] # min slope
    SlopesEPNC [SlopesEPNC$Lake == L & SlopesEPNC$DIV == V ,]$maxSl95 <- boot.ci(BOOT, conf=0.95, type="bca", index=1)$bca[[5]] # max slope
    SlopesEPNC [SlopesEPNC$Lake == L & SlopesEPNC$DIV == V ,]$minSl90 <- boot.ci(BOOT, conf=0.90, type="bca", index=1)$bca[[4]] # min slope
    SlopesEPNC [SlopesEPNC$Lake == L & SlopesEPNC$DIV == V ,]$maxSl90 <- boot.ci(BOOT, conf=0.90, type="bca", index=1)$bca[[5]] # max slope
  }
}


########## average uptake rate ##########

#SlopesEPmr

for (L in Lakes) {
  for (V in DIV_ind) {
    BOOT <- boot(data= avEP[avEP$Lake == L,] , boot_Slope , R=R, V=V, Y="mean.r")
    SlopesEPmr [SlopesEPmr$Lake == L & SlopesEPmr$DIV == V ,]$medSlope <- median(BOOT$t[,1]) # median slope
    SlopesEPmr [SlopesEPmr$Lake == L & SlopesEPmr$DIV == V ,]$meanSlope <- mean(BOOT$t[,1]) # mean slope
    SlopesEPmr [SlopesEPmr$Lake == L & SlopesEPmr$DIV == V ,]$minSl99 <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
    SlopesEPmr [SlopesEPmr$Lake == L & SlopesEPmr$DIV == V ,]$maxSl99 <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # max slope
    SlopesEPmr [SlopesEPmr$Lake == L & SlopesEPmr$DIV == V ,]$minSl95 <- boot.ci(BOOT, conf=0.95, type="bca", index=1)$bca[[4]] # min slope
    SlopesEPmr [SlopesEPmr$Lake == L & SlopesEPmr$DIV == V ,]$maxSl95 <- boot.ci(BOOT, conf=0.95, type="bca", index=1)$bca[[5]] # max slope
    SlopesEPmr [SlopesEPmr$Lake == L & SlopesEPmr$DIV == V ,]$minSl90 <- boot.ci(BOOT, conf=0.90, type="bca", index=1)$bca[[4]] # min slope
    SlopesEPmr [SlopesEPmr$Lake == L & SlopesEPmr$DIV == V ,]$maxSl90 <- boot.ci(BOOT, conf=0.90, type="bca", index=1)$bca[[5]] # max slope
  }
}



################# Nutrients ################


################# DIN ###########

#SlopesDIN

for (L in Lakes) {
  for (V in DIV_ind) {
    BOOT <- boot(data= NUT[NUT$Lake == L,] , boot_Slope , R=R, V=V, Y="NO23")
    SlopesDIN [SlopesDIN$Lake == L & SlopesDIN$DIV == V ,]$medSlope <- median(BOOT$t[,1]) # median slope
    SlopesDIN [SlopesDIN$Lake == L & SlopesDIN$DIV == V ,]$meanSlope <- mean(BOOT$t[,1]) # mean slope
    SlopesDIN [SlopesDIN$Lake == L & SlopesDIN$DIV == V ,]$minSl99 <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[4]] # min slope
    SlopesDIN [SlopesDIN$Lake == L & SlopesDIN$DIV == V ,]$maxSl99 <- boot.ci(BOOT, conf=0.99, type="bca", index=1)$bca[[5]] # max slope
    SlopesDIN [SlopesDIN$Lake == L & SlopesDIN$DIV == V ,]$minSl95 <- boot.ci(BOOT, conf=0.95, type="bca", index=1)$bca[[4]] # min slope
    SlopesDIN [SlopesDIN$Lake == L & SlopesDIN$DIV == V ,]$maxSl95 <- boot.ci(BOOT, conf=0.95, type="bca", index=1)$bca[[5]] # max slope
    SlopesDIN [SlopesDIN$Lake == L & SlopesDIN$DIV == V ,]$minSl90 <- boot.ci(BOOT, conf=0.90, type="bca", index=1)$bca[[4]] # min slope
    SlopesDIN [SlopesDIN$Lake == L & SlopesDIN$DIV == V ,]$maxSl90 <- boot.ci(BOOT, conf=0.90, type="bca", index=1)$bca[[5]] # max slope
  }
}



End  <- Sys.time()
End-Start

##############################################################################################################################

G1 <- ggplot(Slopes_maxBM , aes(x=medSlope, y=DIV, colour=Lake))+
  geom_point(size=3)+
  facet_wrap(~Lake,nrow=1)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0.3, size = 1, aes(xmin=minSl99, xmax=maxSl99))+
  geom_errorbarh(height=0.25, aes(xmin=minSl95, xmax=maxSl95))+
  geom_errorbarh(height=0.2, aes(xmin=minSl90, xmax=maxSl90))+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1),breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="maximum Biomass",x="Slope")


G2 <- ggplot(SlopesEPNC, aes(x=medSlope, y=DIV, colour=Lake))+
  geom_point(size=3)+
  facet_wrap(~Lake,nrow=1)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0.3, size = 1, aes(xmin=minSl99, xmax=maxSl99))+
  geom_errorbarh(height=0.25, aes(xmin=minSl95, xmax=maxSl95))+
  geom_errorbarh(height=0.2, aes(xmin=minSl90, xmax=maxSl90))+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1), breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="Number of C sources", x="Slope")


G3 <- ggplot(SlopesEPmr, aes(x=medSlope, y=DIV, colour=Lake))+
  geom_point(size=3)+
  facet_wrap(~Lake,nrow=1)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0.3, size = 1, aes(xmin=minSl99, xmax=maxSl99))+
  geom_errorbarh(height=0.25, aes(xmin=minSl95, xmax=maxSl95))+
  geom_errorbarh(height=0.2, aes(xmin=minSl90, xmax=maxSl90))+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1), breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="mean carbon uptake rate", x="Slope")


G4 <- ggplot(SlopesCV, aes(x=medSlope, y=DIV, colour=Lake))+
  geom_point(size=3)+
  facet_wrap(~Lake,nrow=1)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0.3, size = 1, aes(xmin=minSl99, xmax=maxSl99))+
  geom_errorbarh(height=0.25, aes(xmin=minSl95, xmax=maxSl95))+
  geom_errorbarh(height=0.2, aes(xmin=minSl90, xmax=maxSl90))+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1),breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="Stability",x="Slope")


G5 <- ggplot(SlopesDIN, aes(x=medSlope, y=DIV, colour=Lake))+
  geom_point(size=3)+
  facet_wrap(~Lake,nrow=1)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_errorbarh(height=0.3, size = 1, aes(xmin=minSl99, xmax=maxSl99))+
  geom_errorbarh(height=0.25, aes(xmin=minSl95, xmax=maxSl95))+
  geom_errorbarh(height=0.2, aes(xmin=minSl90, xmax=maxSl90))+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1), breaks=c(-1,0,1))+
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  labs(title="nutrient depletion (DIN)", x="Slope")

grid.arrange(G1,G2,G3,G4,G5)
Figure_3 <- arrangeGrob(G1,G2,G3,G4,G5)

ggsave("Figure_3.pdf", Figure_3, width = 8, height = 8)

save.image("bootstrap_10000_3CI.RData")

##### without runing code, load image #####

load("bootstrap_10000_new.RData")

######################## alternative test with Spearman's rank correlation #####

RankDIN <- RankEPmr <- RankEPNC <- 
  RankCV  <- Rank_maxBM <- data.frame(
    Lake = rep(Lakes, each=3), DIV = rep(DIV_ind,4),rho = NA,
    pval = NA)


################ maximum Biomass ###########

for (L in Lakes) {
  for (V in DIV_ind) {
    COR <- cor.test(maxBM[maxBM$Lake == L,"Cells"], maxBM[maxBM$Lake == L,V],
                    method = "spearman")
    Rank_maxBM [Rank_maxBM$Lake == L & Rank_maxBM$DIV == V ,]$rho <- round(COR$estimate, 2)
    Rank_maxBM [Rank_maxBM$Lake == L & Rank_maxBM$DIV == V ,]$pval <- round(COR$p.value, 3)
    
  }
}

########## average uptake rate ##########

for (L in Lakes) {
  for (V in DIV_ind) {
    COR <- cor.test(avEP[avEP$Lake == L, "mean.r"], avEP[avEP$Lake == L, V ],
                    method = "spearman")
    RankEPmr [RankEPmr$Lake == L & RankEPmr$DIV == V ,]$rho <- round(COR$estimate, 2)
    RankEPmr [RankEPmr$Lake == L & RankEPmr$DIV == V ,]$pval <- round(COR$p.value, 3)
    
  }
}

########## average NC ##########

for (L in Lakes) {
  for (V in DIV_ind) {
    COR <- cor.test(avEP[avEP$Lake == L, "NC"], avEP[avEP$Lake == L, V ],
                    method = "spearman")
    RankEPNC [RankEPNC$Lake == L & RankEPNC$DIV == V ,]$rho <- round(COR$estimate, 2)
    RankEPNC [RankEPNC$Lake == L & RankEPNC$DIV == V ,]$pval <- round(COR$p.value, 3)
    
  }
}

################ Stability - CV ###########

for (L in Lakes) {
  for (V in DIV_ind) {
    COR <- cor.test(CV[CV$Lake == L, "Stability"], CV[CV$Lake == L, V ],
                    method = "spearman")
    RankCV [RankCV$Lake == L & RankCV$DIV == V ,]$rho <- round(COR$estimate, 2)
    RankCV [RankCV$Lake == L & RankCV$DIV == V ,]$pval <- round(COR$p.value, 3)
    
  }
}


################# Nutrients ################

for (L in Lakes) {
  for (V in DIV_ind) {
    COR <- cor.test(NUT[NUT$Lake == L, "DIN"], NUT[NUT$Lake == L, V ],
                    method = "spearman")
    RankDIN [RankDIN$Lake == L & RankDIN$DIV == V ,]$rho <- round(COR$estimate, 2)
    RankDIN [RankDIN$Lake == L & RankDIN$DIV == V ,]$pval <- round(COR$p.value, 3)
    
  }
}


####### colour code p-values in graph ###
RankDIN$sig  <- 0
RankDIN[ RankDIN$pval <= 0.05, ]$sig  <- 1 # no significant correlations

RankEPmr$sig  <- 0
RankEPmr[ RankEPmr$pval <= 0.05, ]$sig  <- 1

RankEPNC$sig  <- 0
RankEPNC[ RankEPNC$pval <= 0.05, ]$sig  <- 1
  
RankCV$sig  <- 0
RankCV[ RankCV$pval <= 0.05, ]$sig  <- 1

Rank_maxBM$sig  <- 0
Rank_maxBM[ Rank_maxBM$pval <= 0.05, ]$sig  <- 1



G1R <- ggplot(Rank_maxBM, aes(x=rho, y=DIV, fill = Lake, label = pval, colour = as.factor(sig)))+
  geom_vline(xintercept = 0, linetype="dashed", colour = "lightgrey")+
  geom_point(size = 3, shape = 25)+
  geom_text(size=3, vjust = -1)+
  facet_wrap(~Lake,nrow=1)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1), breaks=c(-1, 0, 1))+
  scale_fill_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  scale_colour_manual(values=c("black", "red"))+
  labs(title="maximum Biomass", x="Spearman's rho")

G3R <- ggplot(RankEPmr, aes(x=rho, y=DIV, fill = Lake, label = pval, colour = as.factor(sig)))+
  geom_vline(xintercept = 0, linetype="dashed", colour = "lightgrey")+
  geom_point(size = 3, shape = 25)+
  geom_text(size=3, vjust = -1)+
  facet_wrap(~Lake,nrow=1)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1), breaks=c(-1, 0, 1))+
  scale_fill_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  scale_colour_manual(values=c("black", "red"))+
  labs(title="mean carbon uptake rate", x="Spearman's rho")

G2R <- ggplot(RankEPNC, aes(x=rho, y=DIV, fill = Lake, label = pval, colour = as.factor(sig)))+
  geom_vline(xintercept = 0, linetype="dashed", colour = "lightgrey")+
  geom_point(size = 3, shape = 25)+
  geom_text(size=3, vjust = -1)+
  facet_wrap(~Lake,nrow=1)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1), breaks=c(-1, 0, 1))+
  scale_fill_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  scale_colour_manual(values=c("black", "red"))+
  labs(title="Number of C sources", x="Spearman's rho")

G4R <- ggplot(RankCV, aes(x=rho, y=DIV, fill = Lake, label = pval, colour = as.factor(sig)))+
  geom_vline(xintercept = 0, linetype="dashed", colour = "lightgrey")+
  geom_point(size = 3, shape = 25)+
  geom_text(size=3, vjust = -1)+
  facet_wrap(~Lake,nrow=1)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1), breaks=c(-1, 0, 1))+
  scale_fill_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  scale_colour_manual(values=c("black", "red"))+
  labs(title="Stability", x="Spearman's rho")

G5R <- ggplot(RankDIN, aes(x=rho, y=DIV, fill = Lake, label = pval, colour = as.factor(sig)))+
  geom_vline(xintercept = 0, linetype="dashed", colour = "lightgrey")+
  geom_point(size = 3, shape = 25)+
  geom_text(size=3, vjust = -1)+
  facet_wrap(~Lake,nrow=1)+
  theme_bw(base_size=15)+
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-1, 1), breaks=c(-1, 0, 1))+
  scale_fill_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  scale_colour_manual(values=c("black", "red"))+
  labs(title="nutrient depletion (DIN)", x="Spearman's rho")

grid.arrange(G1R,G2R,G3R,G4R,G5R)

Figure_S_4 <- arrangeGrob( G1R,G2R,G3R,G4R,G5R)

ggsave( filename = "Figure_S_4.pdf", plot = Figure_S_4, width = 10, height = 10)










################### with alternative EF and by DATE where possible ##############




################# Biomass ###############

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

################# EcoLog.txt ###############

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

##############################################################################################################################

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


grid.arrange(G1_2,G2_2,G3_2)


################################################################################
#                             scatterplots                                     #
################################################################################

###### maxBM ########

#reshape to long format

maxBM_l <- melt( maxBM, id.vars = c("Lake","DIL","Cells"))
maxBM_l$variable <- factor( maxBM_l$variable, levels = c( "S", "PSEs", "Hill1"))

G1_s <- ggplot( maxBM_l, aes( x = value, y = Cells, colour= Lake))+
  geom_point( ) +
  facet_wrap( ~ variable * Lake, scales = "free_x") +
  stat_smooth( method = "lm", se = FALSE, linetype = "dashed") +
  theme_bw( base_size = 15)+
  theme(legend.position = "none")+
  labs( title = "maximum Biomass", y = "Cells", x = "Diversity") +
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))




###### Stability ########

#reshape to long format
CV_l <- melt( CV, id.vars = c("Lake","DIL","Stability"))
CV_l$variable <- factor( CV_l$variable, levels = c( "S", "PSEs", "Hill1"))

G2_s <- ggplot( CV_l, aes( x = value, y = Stability, colour= Lake))+
  geom_point( ) +
  facet_wrap( ~ variable * Lake, scales = "free_x") +
  stat_smooth( method = "lm", se = FALSE, linetype = "dashed") +
  theme_bw( base_size = 15)+
  theme(legend.position = "none")+
  labs( title = "Stability", y = "Stability", x = "Diversity") +
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))

###### Ecolog #########

#reshape to long format
EP_l <- melt( EP[, c("Lake", "Hill1", "S", "PSEs", "NC", "mean.r")],
              id.vars = c("Lake","NC", "mean.r"))
EP_l$variable <- factor( EP_l$variable, levels = c( "S", "PSEs", "Hill1"))

G3_s <- ggplot( EP_l, aes( x = value, y = NC, colour= Lake))+
  geom_point( ) +
  facet_wrap( ~ variable * Lake, scales = "free_x") +
  stat_smooth( method = "lm", se = FALSE, linetype = "dashed") +
  theme_bw( base_size = 15)+
  theme(legend.position = "none")+
  labs( title = "Number of C sources", y = "NC", x = "Diversity") +
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))


G4_s <- ggplot( EP_l, aes( x = value, y = mean.r, colour= Lake))+
  geom_point( ) +
  facet_wrap( ~ variable * Lake, scales = "free_x") +
  stat_smooth( method = "lm", se = FALSE, linetype = "dashed") +
  theme_bw( base_size = 15)+
  theme(legend.position = "none")+
  labs( title = "mean carbon uptake rate", y = "average r", x = "Diversity") +
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))


###### Nutrients #########

#reshape to long format
NUT_l <- melt( NUT[, c("Lake", "Hill1", "S", "PSEs", "DIN")],
              id.vars = c("Lake", "DIN"))
NUT_l$variable <- factor( NUT_l$variable, levels = c( "S", "PSEs", "Hill1"))

G5_s <- ggplot( NUT_l, aes( x = value, y = DIN, colour= Lake))+
  geom_point( ) +
  facet_wrap( ~ variable * Lake, scales = "free_x") +
  stat_smooth( method = "lm", se = FALSE, linetype = "dashed") +
  theme_bw( base_size = 15)+
  theme(legend.position = "none")+
  labs( title = "Nutrient depletion", y = "dissolved inorganic Nitrogen",
        x = "Diversity") +
  scale_colour_manual(values=c("orange", "darkred","darkgreen","darkblue"))



grid.arrange(G1_s,G3_s,G4_s,G2_s,G5_s)

Figure_4 <- arrangeGrob(G1_s,G3_s,G4_s,G2_s,G5_s)
ggsave("Figure_4.pdf", Figure_4, width = 18, height = 25)





####################### correlation of slopes ##################################

# vector with DFs with Slopes

RES <- c("Slopes_maxBM","SlopesCV","SlopesEPmr","SlopesEPNC",
         "SlopesDIN")

#join into common DF
Slopes_combined <- Slopes_maxBM[0,]
Slopes_combined$EF <- character()

for (i in 1:length(RES)) {
  Temp <- get(RES[i])
  Temp$EF <- RES[i]
  Slopes_combined <- rbind(Slopes_combined,Temp)
  
}

# change factor level of diversity
Slopes_combined$DIV <- factor(Slopes_combined$DIV, levels = c("S", "PSEs", "Hill1"))

# invert sign for NUT
Slopes_combined[ Slopes_combined$EF == "SlopesDIN",][,3:6]  <- apply( 
  Slopes_combined[ Slopes_combined$EF == "SlopesDIN",][,3:6], 2, function(x) -1*x)


ggplot(Slopes_combined, aes(x=DIV,y=medSlope,shape = EF, fill = Lake,
                            ymin=minSl, ymax=maxSl))+
  geom_point( size = 4, position = position_jitter(w = 0.1, h = 0), alpha=0.8 )+
  facet_wrap(~ Lake)+
  scale_fill_manual(values=c("orange", "darkred","darkgreen","darkblue"))+
  geom_hline(yintersect = 0, linetype = "dashed")+
  labs( x = "diversity metric", y = "median slopes", 
        title = "median slope of the diversity ~ EF relationship") +
  theme_bw(base_size = 15) + 
  guides(colour ="none", fill = "none")+
  theme(legend.position="bottom")+
  scale_shape_manual( values = c(21:25),
    name="Ecosystem Function",
                    breaks=c(levels(as.factor(Slopes_combined$EF))),
                    labels=c("maximum\nBiomass", "Stability", "Nutrient\ndeplition",
                             "number of\nC sources", "C source\nuptake rate"))

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


