###############
# Script to make a pie chart of the Lietrature overview summerized in the google spreadsheets

# source data are from here:
# https://docs.google.com/spreadsheets/d/18ssd0m6cfpi70IRSRYkUvHUoKJbpylEZgpRaS825dbw/edit#gid=0

# based on the categorization here:
# https://docs.google.com/spreadsheets/d/1fR_uzB27T-mWSeMaeCHYJ9kUve5AJ8lMZRqbpTX-vOg/edit#gid=0
##############

library(ggplot2)
library("reshape2")

#############

# import Data
LitDat <- read.table("Literatire_overview.txt", sep = "\t", header = T)

# code missing values  as 0
LitDat[ is.na(LitDat)] <- 0

# reshape data
LitDatR <- melt(LitDat, id.var = c("EF", "total"))

# sort variable
LitDatR$variable <- factor(LitDatR$variable, levels = c("none", "pos", "neg", "amb"))

# sort factor level after total
LitDatR <- LitDatR[ with( LitDatR, order( total, EF, variable)), ]
LitDatR$EF <- factor(LitDatR$EF, levels = unique(LitDatR$EF), labels = c("enhancing \n plant productivity",
                                                                         "invasion resistance",
                                                                         "extracellular enzyme \n multifunctionality",
                                                                         "temporal stability \n of biomass",
                                                                         "bacterial activity",
                                                                         "nitrogen cycling",
                                                                         "resilience",
                                                                         "resistance",
                                                                         "degradation \n of carbon substrates",
                                                                         "yield",
                                                                         "total"))
                                                                

# add poisition variabel for labels
group_midpoints = function(g) {
  cumsums = c(0, cumsum(g$value))
  diffs = diff(cumsums)
  pos = head(cumsums, -1) + (0.5 * diffs)
  return(data.frame(EF=g$EF, variable = g$variable, value= g$value, pos=pos ))
}

LitDatR = ddply(LitDatR, .(EF), group_midpoints)

LitDatR[LitDatR == 0] <- NA

Lit_graph <- ggplot(LitDatR[LitDatR$EF != "total", ], aes( x = EF, y = value, fill = variable, order = variable))+
  geom_bar(stat = "identity", position = "stack", width = 0.5)+
  geom_text( aes(label = value, y = pos), size = 4)+
  coord_flip()+
  labs( y = "number of relationships in literature", x = "",
        title = "bacterial biodiveristy - ecosystem functioning experiments \n using dillution to extinction",
        fill = NULL)+ 
  scale_fill_manual(values = c( "orange", "darkgreen", "red", "grey"),
                    labels = c("no\nrelationship", "positive\nrelationship", "negative\nrealtionship", "ambiguous\ndata"))+
  theme_bw( base_size = 10)+
  theme(legend.position = "bottom")

ggsave("Figure_5.pdf", Lit_graph, width = 8, height = 6)

####### simulation broad function ########

# "B" <- broad function
# "S"  <- specialized function
# N  <- number of species

COM <- c(rep("A", 90), rep("B", 10))

for (i in 1:100) { print(table(sample(COM,10))) }

