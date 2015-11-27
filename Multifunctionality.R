library(multifunc)

setwd("~/Documents/01_PhD/01_Research/02_rare_Biosphere/R scripts")

ID<-read.table("gbgID",header=TRUE, stringsAsFactors=F)
DIV<-read.table("Diversity_5000_100p.txt",stringsAsFactors=F)
PD <- read.table("PhyloDiv.txt", stringsAsFactors=F)


DIV<-join(ID,DIV) 
DIV <- join(DIV,PD)
DIV <- DIV[ ,c("gbgID", "Lake", "DAT", "DIL", "Hill1", "S", "PSEs")]

# average Diversiyt
avDIV<-ddply(DIV[,which( ! colnames (DIV) %in% c("gbgID","DAT"))],
             .(Lake,DIL), numcolwise(mean,na.rm=T))

############# Biomass #############

## full Dataset ##
BM<-read.table("BMfull.txt",stringsAsFactors=F)

#### maximum Biomass ###

maxBM<-ddply(BM, .(Lake,DIL), function(x) mean(x[with(x, order(-Cells)),][1:5,]$Cells))
colnames(maxBM)<-c("Lake","DIL","Cells")

### exclude undiluted treatment form data  and exclude NA###
maxBM <- na.omit(maxBM[maxBM$DIL !="0",])

################# Stability_detr.txt ###############

CV <- read.table("Stability_detr.txt", sep="\t", header=T, stringsAsFactors = F)

CV <- CV[ , c("Lake","DIL","Stability")]

### exclude undiluted treatment form data and exclude NA ###
CV <- na.omit(CV[CV$DIL !="0",])

################# EcoLog.txt ###############

EP <- read.table("EcoPlates.txt", sep="\t", header=T, stringsAsFactors = F)

# join DIV
EP <- join(ID, EP)

### exclude undiluted treatment form data and exclude NA###
EP <- na.omit(EP[EP$DIL !="0",])

### calculate average EcoLog ###
avEP <- ddply(EP, .(Lake,DIL), summarise, NC = mean(NC, na.rm=T),
              mean.r = mean(mean.r, na.rm=T))


################# Nutrients ################

NUT <- read.table("Nut.txt", sep="\t", header=T, stringsAsFactors=F)
NUT <- join(NUT, ID[ID$DAT=="12Jul",])
NUT <- join(avDIV, NUT)

# calculate total dissolven inorganic nitrogen (DIN)
NUT$DIN <- NUT$NO23 + NUT$NH4

### exclude undiluted treatment form data ###
NUT <- NUT[NUT$DIL !="0",]

NUT <- NUT[, c("Lake","DIL","DIN")]

######## join all data frames #########

MultiDF <- Reduce(function(x, y) merge(x, y, all=TRUE), 
       list(NUT, avEP, CV, maxBM, avDIV[avDIV$DIL != "0",]))

########################  calculate Multifunc #################################

# load and subset data

MultiDF$DIN <- (max(MultiDF$DIN)-(MultiDF$DIN))
  
allVars <- qw(DIN, mean.r, NC, Stability, Cells)
#allVars <- qw(DIN, mean.r, NC, Cells)



varIdx <- which(names(MultiDF) %in% allVars)

vars <- whichVars(MultiDF, allVars) 

mydataThresh<-getFuncsMaxed(MultiDF, vars, threshmin=0.05, threshmax=0.99,
                            prepend=c("Hill1"), maxN =3)

mydataLinearSlopes<-getCoefTab(funcMaxed ~ Hill1, data=mydataThresh,
                               coefVar="Hill1", fun = lm) 

#CN <- colnames(mydataLinearSlopes)
colnames(mydataLinearSlopes) <- CN


ggplot(mydataLinearSlopes, aes(x=thresholds)) +
  geom_ribbon(fill="grey50", aes(x=thresholds*100, ymin=Estimate-
                                   1.96*mydataLinearSlopes[["Std. Error"]],
                                 ymax=Estimate+1.96*mydataLinearSlopes[["Std. Error"]])) +
  geom_point(aes(x=thresholds*100, y=Estimate)) +
  labs(x= "\nThreshold (%)", title = "Hill1") +
  stat_abline(intercept=0, slope=0, lwd=1, linetype=2) +
  theme_bw(base_size=14)


gcPlot<-subset(mydataThresh, mydataThresh$thresholds %in% qw(0.2, 0.4, 0.6, 0.8))
gcPlot$percent<-paste(100*gcPlot$thresholds, "%", sep="")


# plot 4 selected threshold values
ggplot(gcPlot, aes(x=Hill1, y=funcMaxed))+
  geom_point()+
  facet_wrap(~percent)+
  stat_smooth(method="lm", family = quasipoisson(link="identity"),colour="red", lwd=1.2 )+
  labs(y = expression("Number of Functions" >= Threshold), x = ("Hill1"))+
  theme_bw(base_size=15)

################################################################################


EP <- (read.table("EcoLb_NLS_07-12"))

# mark replicates on plate
EP$Well <- rep(c(1:32), 3)

# mark NA as 0
EP[is.na(EP)] <- 0

# take mean over replicates
EP <- ddply(EP, .(Plate, Well), summarize, r = median(r))

# exclude well 1 (blank)
EP <- EP[EP$Well != "1",]

#rename Plate
colnames(EP) <- c("BOT", "Well", "r")

# cast dataframe in wide foramt
EP <- dcast(EP, BOT~Well)

#join ID and diveristy data
EP <- join(EP, ID[ID$DAT == "12Jul", c("BOT","Lake","DIL")])

EP <- join(EP, DIV[DIV$DAT == "12Jul", c("Lake", "DIL", "Hill1", "S", "PSEs")])

#exclude dilution 0
EP <- EP[EP$DIL != "0", ]

########################  calculate Multifunc #################################

allVars <- qw(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
              26,27,28,29,30,31,32)

varIdx <- which(names(EP) %in% allVars)

vars <- whichVars(EP, allVars) 

mydataThresh<-getFuncsMaxed(EP, vars, threshmin=0.05, threshmax=0.99,
                            prepend=c("Hill1"), maxN =3)

mydataLinearSlopes<-getCoefTab(funcMaxed ~ Hill1, data=mydataThresh,
                               coefVar="Hill1", fun = lm) 

#CN <- colnames(mydataLinearSlopes)
colnames(mydataLinearSlopes) <- CN


ggplot(mydataLinearSlopes, aes(x=thresholds)) +
  geom_ribbon(fill="grey50", aes(x=thresholds*100, ymin=Estimate-
                                   1.96*mydataLinearSlopes[["Std. Error"]],
                                 ymax=Estimate+1.96*mydataLinearSlopes[["Std. Error"]])) +
  geom_point(aes(x=thresholds*100, y=Estimate)) +
  labs(x= "\nThreshold (%)", title = "Hill1") +
  stat_abline(intercept=0, slope=0, lwd=1, linetype=2) +
  theme_bw(base_size=14)


gcPlot<-subset(mydataThresh, mydataThresh$thresholds %in% qw(0.2, 0.4, 0.6, 0.8))
gcPlot$percent<-paste(100*gcPlot$thresholds, "%", sep="")


# plot 4 selected threshold values
ggplot(gcPlot, aes(x=Hill1, y=funcMaxed))+
  geom_point()+
  facet_wrap(~percent)+
  stat_smooth(method="lm", family = quasipoisson(link="identity"),colour="red", lwd=1.2 )+
  labs(y = expression("Number of Functions" >= Threshold), x = ("Hill1"))+
  theme_bw(base_size=15)










