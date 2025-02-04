


library(ape)
library(phangorn)
library(phytools)


#tree <- read.tree(file="./treepl_examl_T400F.astral.rooted.tre")
tree <- read.tree(file="./treepl_castles_T400F.wastral.rooted.tre")

plot(tree, cex=0.1, show.tip.label=FALSE); axisPhylo()


# Delete outgroup

Ingroup <- Descendants(tree, node=getMRCA(tree, tip=c('Sapayoa_ae_L2329', 'Tyranus_savana_L48429')))[[1]]

tree2 <- keep.tip(tree, tip=Ingroup)

plot(tree2, cex=0.1, show.tip.label=FALSE); axisPhylo()


# Check ultrametricity 
is.ultrametric(tree2)

# Make the tree exactly unrametric:
#tree2 <- nnls.tree(cophenetic(tree2), tree, rooted=TRUE)

tree2 <- force.ultrametric(tree2)

is.ultrametric(tree2)




######################
### CoMET analysis ###
######################

library(TESS)

# May et al. 2016
#"treats the number of specation-rate shifts, extinction-rate shifts, and mass-extinction events as random variables, and estimates their joint posterior distribution."

# CoMET with Empirical Hyperpriors and no mass extinctions #

# tess.analysis(tree2,
#               empiricalHyperPriors = TRUE,
#               MRCA = TRUE,
#               estimateNumberMassExtinctions = FALSE,
#               MAX_ITERATIONS = 1000000,
#               dir = "CoMETtreepl_concat_T400F.astral")
#               
# output <- tess.process.output(dir="CoMETtreepl_concat_T400F.astral", tree= tree2)


tess.analysis(tree2,
              empiricalHyperPriors = TRUE,
              MRCA = TRUE,
              estimateNumberMassExtinctions = FALSE,
              MAX_ITERATIONS = 1000000,
              dir = "CoMETtreepl_castles_T400F.wastral")
              
output <- tess.process.output(dir="CoMETtreepl_castles_T400F.wastral", tree= tree2)


# DIAGNOSTICS #

effectiveSize(output$numSpeciationCategories)
effectiveSize(output$numExtinctionCategories)

layout(matrix(1:4,nrow=2,ncol=2,byrow=TRUE))
tess.plot.singlechain.diagnostics(output, parameters=c("speciation rates", "extinction rates"), las=2)


# PLOTS

layout.mat <- matrix(1:4,nrow=2,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.output(output, fig.types = c("speciation rates", "speciation shift times","extinction rates", "extinction shift times"), las=2)


tess.plot.output(output, fig.types = "net-diversification rates", las=2, lwd=3)
