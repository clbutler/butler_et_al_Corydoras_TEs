#########################
#Phylogentic paper
#Phylogeny
#Author: Chris Butler
###########################

#RESET
rm(list = ls())
par(mfrow=c(1,1))

#LOAD LIBRARIES 
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)

####### mtDNAPhylogenetic Tree##########
library('phytools')
library('ape')
library('devtools')
library('l1ou')
library('gdata')




#set working directory to the danio folder
setwd('/Users/chrisbutler/Documents/chapter4/transcriptomes_NCBI/repeatmasker_output /Danio/RMtrips/')

myfiles = list.files(pattern=".csv", full.names=TRUE) 
read_csv_filename <- function(filename){
  ret <- read.csv(filename, header = TRUE, stringsAsFactors=FALSE)
  ret$Source <- filename #EDIT
  ret
}
RMfiles = lapply(myfiles, read_csv_filename)
combined_RMfiles <-do.call(rbind, RMfiles)

setwd('/Users/chrisbutler/Documents/chapter4//transcriptomes_NCBI/')
transcript_lengths <- read.csv('NCBI_transcriptomelengths.csv', stringsAsFactors=FALSE)


combined_RMfiles$Species[grepl("Sample56_", combined_RMfiles$Source)] <- "C. maculifer"
combined_RMfiles$Species[grepl("sample4.", combined_RMfiles$Source)] <- "A. fuscoguttatus"
combined_RMfiles$Species[grepl("sample39.", combined_RMfiles$Source)] <- "S. prionotus"
combined_RMfiles$Species[grepl("sample63.", combined_RMfiles$Source)] <- "C. hastatus"
combined_RMfiles$Species[grepl("sample28.", combined_RMfiles$Source)] <- "C. elegans"
combined_RMfiles$Species[grepl("sample31.", combined_RMfiles$Source)] <- "C. paleatus"
combined_RMfiles$Species[combined_RMfiles$Source == "./noasterisk_cleaned6_sample2.fasta.out_RM_TRIPS.csv"] <- "C. aeneus"
combined_RMfiles$Species[grepl("sample15.", combined_RMfiles$Source)] <- "C. haraldschultzei"
combined_RMfiles$Species[grepl("sample57.", combined_RMfiles$Source)] <- "C. araguaiaensis"
combined_RMfiles$Species[grepl("pCoco", combined_RMfiles$Source)] <- "I. punctatus"


combined_RMfiles <- combined_RMfiles %>% dplyr::group_by(Species) %>% mutate(., totalTElength = sum(mergedfraglength))

phylo <- merge(combined_RMfiles, transcript_lengths, by = c("Species"))
phylo <- mutate(phylo, percTEs = (totalTElength/Length )*100)

phylo_summary <- phylo %>% dplyr::select(Species, percTEs)
phylo_summary <- unique(phylo_summary)

#plot tree
#import tree
setwd('/users/chrisbutler/Documents/trees')
      tree <- read.newick('/users/chrisbutler/Documents/trees/raxml_radtree_edit.tre') #Load Corydoras newick trees
species <- tree$tip.label

tree <- chronos(tree)
is.ultrametric(tree)
#subset nodes I want
prunes <- c("Mega1_1[&!rotate=true]","Aen1_1","Fow1_1","Apo1_1", "Skro1_1", "Eleg1_1", "Pyg1_1", "Nat1_1", "Ara1_1", "Imi1_1")
pruned.tree <- drop.tip(tree, setdiff(tree$tip.label, prunes))                                             
plot(pruned.tree)
nodelabels()
#rename tip labels
pruned.tree$tip.label[pruned.tree$tip.label == "Skro1_1"] <- "S. prionotus"
pruned.tree$tip.label[pruned.tree$tip.label == "Aen1_1"] <- "C. aeneus"
pruned.tree$tip.label[pruned.tree$tip.label == "Nat1_1"] <- "C. paleatus"
pruned.tree$tip.label[pruned.tree$tip.label ==  "Eleg1_1"] <- "C. elegans"
pruned.tree$tip.label[pruned.tree$tip.label ==  "Fow1_1"] <- "C. maculifer"
pruned.tree$tip.label[pruned.tree$tip.label ==  "Imi1_1"] <- "C. haraldschultzei"
pruned.tree$tip.label[pruned.tree$tip.label ==  "Ara1_1"] <- "C. araguaiaensis"
pruned.tree$tip.label[pruned.tree$tip.label ==  "Pyg1_1"] <- "C. hastatus"
pruned.tree$tip.label[pruned.tree$tip.label ==  "Apo1_1"] <- "A. fuscoguttatus"
pruned.tree$tip.label[pruned.tree$tip.label ==  "Mega1_1[&!rotate=true]"] <- "I. punctatus"

pruned.tree <- rotate(pruned.tree, 12)
pruned.tree <- rotate(pruned.tree, 13)
pruned.tree <- rotate(pruned.tree, 14)
pruned.tree <- rotate(pruned.tree, 15)

#create vector with masked % per species
maskedvec <- phylo_summary$percTEs
names(maskedvec) <- phylo_summary$Species
maskedvec

# contmap plot masked per species - blue
maskedtree <-contMap(pruned.tree, maskedvec, plot = FALSE, sig = 1)
maskedtree <- setMap(maskedtree,colors=c("white","white", "lightblue1", "lightblue1", "dodgerblue3", "navyblue", 'navyblue' ))
plot(maskedtree,leg.txt = "TE Abundance %", lwd = 10, fsize = 0.8, ftype = "b", legend = 52)



#bootstrap phylogenetic shifts
L1ou_masked <- adjust_data(pruned.tree, maskedvec)
eMODEL <- estimate_shift_configuration(L1ou_masked$tree, L1ou_masked$Y)
result <- l1ou_bootstrap_support(eMODEL, nItrs = 100) #bootstrapping with replacement 100 times
e.l <- round(result$detection.rate*100, digits=1) #obtain bootstrap support values
# to avoid annotating edges with support at or below 10%
e.l <- ifelse(e.l>10, paste0(e.l,"%"), NA)
#plot
plot(eMODEL, edge.label=e.l, edge.label.ann=TRUE, label.offset=0.02)
