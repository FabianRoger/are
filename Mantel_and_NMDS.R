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
################################################################################

library(ggplot2)

################################################################################

# import files
################################################################################

#unifrac_dist
load("unifrac_dist.RData")



################################################################################
