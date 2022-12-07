#
#relationship between genome size and TE abundance
#
library(ape)
library(geiger)
library(nlme)
library(phytools)

######GENOME SIZE VS TE COPIES#######
setwd('/Users/chrisbutler/Documents/chapter4/abundancevsgenome_size/')
gensize <- read.csv('genome_sizes.csv') #import data on genome sizes 
gensize <- gensize %>% dplyr::select(Species, Hap.C.value..pg., Lineage) #select species, c values and lineage
gensize$Species <- as.character(gensize$Species)
gensize <- gensize %>% dplyr::filter(!(Hap.C.value..pg. =="")) #remove species with no entries
gensize$Hap.C.value..pg. <- as.character(gensize$Hap.C.value..pg.) 
gensize$Hap.C.value..pg. <- as.numeric(gensize$Hap.C.value..pg.)
gensize <- subset(gensize, Lineage != 'Lineage 0') #remove outgroups
gensize <- gensize %>% group_by(Lineage) %>% mutate(., LineageMeanGenomeSize = mean(Hap.C.value..pg.)) #average genome size per lineage
gensize <- gensize %>% dplyr::select(Lineage, LineageMeanGenomeSize)
gensize <- unique(gensize)
gensize$Lineage <-gsub('Lineage ', '', gensize$Lineage) #tidy column

#insert TE abundance
setwd('/Users/chrisbutler/Documents/chapter4/abundancevsgenome_size/')
TEabundance <- read.csv('TEabundances_perc.csv')
TEabundance <- merge(TEabundance, gensize, by = c('Lineage'))

TEabundance$Lineage <- as.factor(TEabundance$Lineage)
TEabundance$species <- TEabundance$Species

TEabundance <- TEabundance %>% tibble::column_to_rownames(var="Species")

#now get phylogeny information
#plot tree
#import tree
setwd('/users/chrisbutler/Documents/trees')
tree <- read.newick('/users/chrisbutler/Documents/trees/raxml_radtree_edit.tre') #Load Corydoras newick trees
species <- tree$tip.label

tree <- force.ultrametric(tree, method=c("extend"))
is.ultrametric(tree)
#subset nodes I want
prunes <- c("Aen1_1","Fow1_1","Apo1_1", "Skro1_1", "Eleg1_1", "Pyg1_1", "Nat1_1", "Ara1_1", "Imi1_1")
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


pruned.tree <- rotate(pruned.tree, 10)
pruned.tree <- rotate(pruned.tree, 12)
pruned.tree <- rotate(pruned.tree, 11)
pruned.tree <- rotate(pruned.tree, 13)

obj<-name.check(pruned.tree, TEabundance)
obj

#When Lambda = 0
pglsModel_zero <- gls(TE_EDTA ~ LineageMeanGenomeSize, correlation = corPagel(value = 0, phy = pruned.tree,
                                                                              fixed = TRUE, form = ~species), data = TEabundance) 
summary(pglsModel_zero)



#When Lambda = 1
pglsModel_one<- gls(TE_EDTA ~ LineageMeanGenomeSize, correlation = corPagel(1, phy = pruned.tree, form = ~species, fixed = TRUE),
                      data = TEabundance)
summary(pglsModel_one)

#When Lambda = ML
pglsModel_ML <- gls(TE_EDTA ~ LineageMeanGenomeSize,correlation=corPagel(1, phy = pruned.tree,form = ~species, fixed = FALSE), method = "ML", data = TEabundance)
summary(pglsModel_ML)




#http://blog.phytools.org/2012/11/fitting-model-in-phylogenetic.html - with reference to Liam



ggplot(TEabundance, aes(x = LineageMeanGenomeSize, y = TE_EDTA)) + geom_point(size = 3) + geom_abline(intercept = 5.13, slope = 0.55, colour = "orange", size = 1.5, alpha = .6) +
  xlab("\nMean Genome Size per Lineage / pg") +
  ylab('Expressed TE Abundance \ %\n') + theme_bw() +
  theme(text = element_text(size = 20)) +
  geom_text(aes(label=Lineage),hjust=-0.7, vjust=-0.9, size = 6) + geom_abline(intercept = 5.10, slope = 0.51, colour = "darkgreen", size = 1.5, alpha = .6) + 
  geom_abline(intercept = 5.37, slope = 0.66, colour = '#7570b3', size = 1.5, alpha = .6) + ylim(min = 4.0,max = 8.0)




  

#setwd('/Users/chrisbutler/Documents/chapter4/abundancevsgenome_size/')
#ggsave('abundancevsgenomesize-withlineage.pdf')
