#Horizontal Transfer Analysis #####


##Reset
rm(list = ls())
par(mfrow=c(1,1))

##Load libraries##

library(stringr)
library(ggplot2)
library(reshape2)
library(seqinr)
library(splitstackshape)  
library(tidyr)
library(dplyr)
library(forcats)
library(rstudioapi)

#load custom functions
read_csv_filename <- function(filename){
  ret <- read.csv(filename, fill = TRUE, row.names=NULL, stringsAsFactors = FALSE)
  ret$Source <- filename #EDIT
  ret
}

current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )



###################################
#STEP 1 - input RepeatMasker files#
###################################

setwd('../RepeatMasker_outputs/Cory_TE/RM_trips/')

myfiles = list.files(pattern=".csv", full.names=TRUE) 

RMfiles = lapply(myfiles, read_csv_filename)
transcriptome <-do.call(rbind, RMfiles)
transcriptome <- filter(transcriptome, Source != "./noasterisk_IpCoco_mRNA.fa.out_RM_TRIPS.csv") #remove channel catfish

noisonosim <- transcriptome



######Adding TE type and class column#####
noisonosim$TE_type[grepl("TcMar", noisonosim$matching_class)] <- "TIR TC1-Mariner-like"


#filter just Mariner
noisonosim <- noisonosim %>% filter(TE_type == 'TIR TC1-Mariner-like')

###################################
#STEP 2 - add length of all transcripts and only keep those > 80% consists of TE
###################################

read_csv_filename <- function(filename){
  ret <- read.csv(filename, fill = TRUE, row.names=NULL, stringsAsFactors = FALSE, header = FALSE, sep = "\t")
  ret$Source <- filename #EDIT
  ret
}


current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
setwd('./transcript_lengths/')
myfiles = list.files(pattern="*.txt", full.names=TRUE) #myfiles is now a list of every headers csv file


all_sequencenames = lapply(myfiles, read_csv_filename) #read_csv_filenames of the myfiles list
combined_sequencenames <- do.call(rbind, all_sequencenames) #combined all the csv files into one dataframe
combined_sequencenames$Source <- gsub('.*\\_', '', combined_sequencenames$Source) #clean up file name column
combined_sequencenames$Source <- gsub('.fasta.gz.txt', '', combined_sequencenames$Source) #clean up file name column
combined_sequencenames <- rename(combined_sequencenames, length = V2)
combined_sequencenames <- rename(combined_sequencenames, qry_id = V1)

noisonosim$Source <- gsub('.fasta.out_RM_TRIPS.csv', '', noisonosim$Source)
noisonosim$Source <- gsub('.*\\_', '', noisonosim$Source)
noisonosim$Source <- gsub('S', 's', noisonosim$Source)
noisonosim$Source[(noisonosim$Source == 'final')] <- 'sample56'

noisonosim <- merge(noisonosim, combined_sequencenames, by = c('Source', 'qry_id'), all.x = TRUE )


noisonosim <- mutate(noisonosim, perc_TE_transcript = (mergedfraglength/length)*100)

test_80 <- filter(noisonosim, perc_TE_transcript >= 80) #if extracting Mariner elements over 80% of length of transcript - 8,010 out of 56,492 - ~14% of elements

###################################
#STEP 3 - write bed file of mariner elements
###################################

current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
setwd('./bed_files')
x <- split(test_80, test_80$Source)

filter_list <- c("sample56", "sample15", "sample2", "sample28", 'sample31', 'sample39', 'sample4', 'sample57', 'sample63')

for (i in filter_list) {
  x[[i]] <- dplyr::select(x[[i]],'qry_id', 'merged_qrystart', 'merged_qryend', 'repeat_id')
}

for (i in filter_list) {
  write.table(x[[i]], file = paste0("mariner", i ,".bed"), col.names = FALSE, sep = '\t', quote = FALSE, row.names = FALSE)
}

###########
#STEP 4 - run blast searches against Mariner elements, 
###########

#With the bed files found in the bed_files output you can run two shell scripts 
#1) getfasta_script.sh - this turns bed into a fasta - adjust for species needed (https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html)
#2) blast_repbase.sh - this runs a blast search against extracted elements and RepBase - adjust for species needed

#########################################
#STEP 5 - import and clean blast outfmts
#########################################
read_csv_filename <- function(filename){
  ret <- read.csv(filename, fill = TRUE, row.names=NULL, stringsAsFactors = FALSE, header = FALSE, sep = "\t")
  ret$Source <- filename #EDIT
  ret
}

#create a data frame with repeatname and species it originates from
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
repeat_ids <- read.csv('repeatmasker_headers.csv', header = FALSE)
repeat_ids <- repeat_ids %>% dplyr::select(V1, V2)
repeat_ids$V2 <- gsub('@', '', repeat_ids$V2)
repeat_ids <- rename(repeat_ids, repeat_id = V1)
repeat_ids <- rename(repeat_ids, Species = V2)
repeat_ids$repeat_id <- gsub('>', '', repeat_ids$repeat_id)

#import blast files
setwd('./blast_outputs/')
file_list <- list.files(pattern = ".outfmt")
combined_blasthits <- lapply(file_list, read_csv_filename)
combined_blasthits <- do.call(rbind, combined_blasthits)


#select only the best matches
combined_blasthits <- rename(combined_blasthits, qry_id = V1)
combined_blasthits <- rename(combined_blasthits, E_value = V11)
combined_blasthits <- rename(combined_blasthits, repeat_id = V2)
combined_blasthits$E_value <- as.numeric(combined_blasthits$E_value)
combined_blasthits <- merge(combined_blasthits, repeat_ids, by = "repeat_id", all.x = TRUE)

#remove any with hits against non mariner elenents 
combined_blasthits <- filter(combined_blasthits, str_detect(repeat_id, "TcMar"))

#Turn into taxonomic groups
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Danio')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'salar')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Xenopus')] <- 'Xenopus sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Esox')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Rana')] <- 'Rana sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Takifugu')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Anolis')] <- 'Anolis sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Gasterosteus')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Cichlidae')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Gasterosteus')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Cichlidae')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Petromyzon')] <- 'Petromyzon sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Oryzias')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Rhizopus')] <- 'Rhizopus sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Drosophila')] <- 'Drosophila sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Lethenteron')] <- 'Lethenteron sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Solenopsis_invicta')] <- 'Solenopsis sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Cyprinidae')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Pleuronectes')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Rhodnius')] <- 'Rhodnius sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Locusta')] <- 'Locusta sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Lepeophtheirus')] <- 'Lepeophtheirus sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Schmidtea')] <- 'Schmidtea sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Belone_belone')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Neogobius')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Oncorhynchus_nerka')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Harpegnathos_saltator')] <- 'Harpegnathos sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Anopheles_albimanus')] <- 'Anopheles sp.'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Perca_fluviatili')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Scophthalmus_maximus')] <- 'Teleost'
combined_blasthits$taxonomic_group[str_detect(combined_blasthits$Species, 'Primate')] <- 'Primate sp.'

combined_blasthits_best <- combined_blasthits %>% group_by(qry_id, Source) %>% slice(which.min(E_value))
combined_blasthits_best <- unique(combined_blasthits_best)
#there are 1,314 corydoras Mariner TEs with hits against the RefBase library 

###############
#STEP 6 - extract best teleost hit for each Mariner element
################

#extract best teleost hit
teleost_hit <- combined_blasthits %>% filter(., taxonomic_group == 'Teleost')
teleost_hit_best <- teleost_hit %>% group_by(qry_id, Source) %>% slice(which.min(E_value))
teleost_hit_best <- teleost_hit_best %>%  dplyr::select(qry_id, Source, V3)
teleost_hit_best <- rename(teleost_hit_best, Teleost_hit_perc = V3)


#merge best hit and best teleost hit
combined_blasthits_best <- merge(combined_blasthits_best, teleost_hit_best, by = c("qry_id", "Source"), all.x = TRUE)
combined_blasthits_best$Teleost_hit_perc <- replace_na(combined_blasthits_best$Teleost_hit_perc, 0)
combined_blasthits_best <- combined_blasthits_best %>% mutate(difference_in_perc = V3 - Teleost_hit_perc)

#if difference between best non-teleost hit and best teleost hit is less than 2% than keep teleost hit
combined_blasthits_best$altered_taxonomic_group[combined_blasthits_best$difference_in_perc < 2] <- 'Teleost'
combined_blasthits_best$altered_taxonomic_group <- ifelse(is.na(combined_blasthits_best$altered_taxonomic_group), combined_blasthits_best$taxonomic_group, combined_blasthits_best$altered_taxonomic_group)

#add Species_Cory
combined_blasthits_best$Species_Cory[combined_blasthits_best$Source == "HT_DN_mariner_2.outfmt6"] <- "C. aeneus"
combined_blasthits_best$Species_Cory[combined_blasthits_best$Source == "HT_DN_mariner_4.outfmt6"] <- "A. fuscoguttatus"
combined_blasthits_best$Species_Cory[combined_blasthits_best$Source == "HT_DN_mariner_31.outfmt6"] <- "C. paleatus"
combined_blasthits_best$Species_Cory[combined_blasthits_best$Source == "HT_DN_mariner_39.outfmt6"] <- "S. prionotus"
combined_blasthits_best$Species_Cory[combined_blasthits_best$Source == "HT_DN_mariner_28.outfmt6"] <- "C. elegans"
combined_blasthits_best$Species_Cory[combined_blasthits_best$Source == "HT_DN_mariner_15.outfmt6"] <- "C. haraldschultzei"
combined_blasthits_best$Species_Cory[combined_blasthits_best$Source == "HT_DN_mariner_57.outfmt6"] <- "C. araguaiaensis"
combined_blasthits_best$Species_Cory[combined_blasthits_best$Source == "HT_DN_mariner_63.outfmt6"] <- "C. hastatus"
combined_blasthits_best$Species_Cory[combined_blasthits_best$Source == "HT_DN_mariner_56.outfmt6"] <- "C. maculifer"

combined_blasthits_best <- combined_blasthits_best %>% ungroup() %>% add_count(name = "all_transposons") 
combined_blasthits_best <- combined_blasthits_best %>% group_by(altered_taxonomic_group) %>% add_count(name = "all_transposonspertaxa") 


combined_blasthits_best <- combined_blasthits_best %>% ungroup() %>% group_by(Species_Cory) %>% add_count(name = "total_transposons")
combined_blasthits_best <- combined_blasthits_best %>% ungroup() %>% group_by(Species_Cory, altered_taxonomic_group) %>% add_count(name = "HitsperTaxa")

total_barplot <- combined_blasthits_best %>% ungroup() %>%  dplyr::select(altered_taxonomic_group, all_transposons, all_transposonspertaxa)
total_barplot <- total_barplot %>% mutate(., perc_per_taxa = (all_transposonspertaxa/all_transposons)*100)
total_barplot <- total_barplot %>% dplyr::select(altered_taxonomic_group, perc_per_taxa)
total_barplot <- unique(total_barplot)
total_barplot$libtype <- 'Corydoras C115 custom/denovo'

###########
#STEP 7 - save output and plotting for figure 6A
###########
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
#write.csv(total_barplot, 'broaddenovo_HT.csv')

ggplot(total_barplot, aes(x = fct_reorder(altered_taxonomic_group, perc_per_taxa), y = perc_per_taxa)) + geom_bar(stat = "identity", colour = "darkblue", fill = "cornflowerblue") + theme_bw() +
  ylab("\n% of Total RepBase Hits") + theme(text = element_text(size = 25))  + 
  theme(axis.text.x=element_text(angle = 290)) + xlab('Closest Taxonomic Match\n') + coord_flip() + scale_y_continuous(limits = c(0,80)) 


##########
#STEP 7 - analysis and plotting for figure 6b
##########

#keep those elements which have not come from teleosts
combined_blasthits_plot <- filter(combined_blasthits_best, altered_taxonomic_group != 'Teleost')


#plotting time
combined_blasthits_plot <- mutate(combined_blasthits_plot, HT_prop_perspecies = (HitsperTaxa/total_transposons)*100)
combined_blasthits_plot <- combined_blasthits_plot %>% dplyr::select(Species_Cory, HT_prop_perspecies)
combined_blasthits_plot <- unique(combined_blasthits_plot)


#add Species_Cory lineage
combined_blasthits_plot$Lineage[combined_blasthits_plot$Species_Cory == "C. maculifer"] <- "1"
combined_blasthits_plot$Lineage[combined_blasthits_plot$Species_Cory == "A. fuscoguttatus"] <- "2"
combined_blasthits_plot$Lineage[combined_blasthits_plot$Species_Cory == "S. prionotus"] <- "3"
combined_blasthits_plot$Lineage[combined_blasthits_plot$Species_Cory == "C. hastatus"] <- "4"
combined_blasthits_plot$Lineage[combined_blasthits_plot$Species_Cory == "C. elegans"] <- "5"
combined_blasthits_plot$Lineage[combined_blasthits_plot$Species_Cory == "C. paleatus"] <- "8.5"
combined_blasthits_plot$Lineage[combined_blasthits_plot$Species_Cory == "C. aeneus"] <- "7"
combined_blasthits_plot$Lineage[combined_blasthits_plot$Species_Cory == "C. haraldschultzei"] <- "8"
combined_blasthits_plot$Lineage[combined_blasthits_plot$Species_Cory == "C. araguaiaensis"] <- "9"

combined_blasthits_plot$Lineage <- as.numeric(combined_blasthits_plot$Lineage)
combined_blasthits_plot <- rename(combined_blasthits_plot, Taxa = altered_taxonomic_group)

ggplot(combined_blasthits_plot, aes(x = fct_reorder(Species_Cory, -Lineage), y = HT_prop_perspecies, fill = Taxa)) + geom_bar(stat = "identity", colour = "white", alpha = .6) +
  theme_bw() + xlab("Corydoradinae Species\n") + ylab("\n Potential Horizontally Transferred Mariner Elements / %") +   theme(text = element_text(size = 20)) +
  scale_fill_brewer(palette = "Paired") + coord_flip() 

combined_blasthits_plot %>% group_by(Species_Cory) %>% summarise(., total_HT = sum(HT_prop_perspecies))

##########
#STEP 8 - additional phylogenetic shift analysis
##########

#does the % of horizontally transferred elements changes across the Cory phylogeny 

HTT <- combined_blasthits_best
HTT <- filter(HTT, altered_taxonomic_group != 'Teleost')
HTT <- cSplit(HTT, 'qry_id', sep = ":") 
HTT <- cSplit(HTT, 'qry_id_2', sep = "-") 
HTT <- mutate(HTT, length_TE = (qry_id_2_2 - qry_id_2_1)+1)
HTT <- HTT %>% group_by(Species_Cory) %>% mutate(., total_HT_length_species = sum(length_TE) )
HTT <- HTT %>% dplyr::select(Species_Cory, total_HT_length_species)
HTT <- unique(HTT)

setwd('../')
transcript_lengths <- read.csv('NCBI_transcriptomelengths.csv', stringsAsFactors=FALSE)
transcript_lengths <- rename(transcript_lengths, Species_Cory = Species)

HTT <- merge(HTT, transcript_lengths, by = "Species_Cory", all.x = TRUE)
HTT <- mutate(HTT, total_HT_length_prop = (total_HT_length_species/Length)*100)

####### Phylogenetic Tree##########
library('phytools')
library('ape')
library('devtools')
library('l1ou')
library('gdata')
library('tidyverse')

#plot tree
#import tree
setwd('Nuclear_tree/')
tree <- read.newick('/users/chrisbutler/Documents/trees/raxml_radtree_edit.tre') #Load Corydoras newick trees
species <- tree$tip.label

tree <-  chronos(tree)
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
pruned.tree$tip.label[pruned.tree$tip.label ==  "Mega1_1[&!rotate=true]"] <- "I. punctatus"

pruned.tree <- rotate(pruned.tree, 10)
pruned.tree <- rotate(pruned.tree, 11)
pruned.tree <- rotate(pruned.tree, 12)
pruned.tree <- rotate(pruned.tree, 13)

#create vector with masked % per species
maskedvec <- HTT$total_HT_length_prop
names(maskedvec) <- HTT$Species_Cory
maskedvec

# contmap plot masked per species - blue
maskedtree <-contMap(pruned.tree, maskedvec, plot = FALSE, sig = 1)
maskedtree <- setMap(maskedtree,colors=c("white","white", "lightblue1", "lightblue1", "dodgerblue3", "navyblue", 'navyblue' ))
plot(maskedtree,leg.txt = "TE Abundance %", lwd = 10, fsize = 0.8, ftype = "b", legend = 52)

#bootstrap phylogenetic shifts
L1ou_masked <- adjust_data(pruned.tree, maskedvec)
eMODEL <- estimate_shift_configuration(L1ou_masked$tree, L1ou_masked$Y)
result <- l1ou_bootstrap_support(eMODEL, nItrs = 100) #bootstrapping with replacement 1,000 times
e.l <- round(result$detection.rate*1000, digits=1) #obtain bootstrap support values
# to avoid annotating edges with support at or below 10%
e.l <- ifelse(e.l>10, paste0(e.l,"%"), NA)
#plot
plot(eMODEL, edge.label=e.l, edge.label.ann=TRUE, label.offset=0.02)

#conclusion - as a proportion of total transcriptome length there is no evidence that a significant increase in horizontally transferred Mariner elements has occurred 

#########################
#STEP 9 - additional pooling interesting horizontally transferred elements.
##########################

HT_elements <- combined_blasthits_best %>% dplyr::select(qry_id, repeat_id, altered_taxonomic_group, Species_Cory)
HT_elements <- rename(HT_elements, blastn_hit = repeat_id)
HT_elements$qry_id <- gsub('\\:.*', '', HT_elements$qry_id)

#combined bed output files (ie every mariner element coodinates - 80% length of transcript ) into a single data frame
combined_bed <-do.call("rbind", lapply(x, as.data.frame))
combined_bed$Source <- rownames(combined_bed)
combined_bed$Source <- gsub('\\..*','',combined_bed$Source)

combined_bed$Species_Cory[grepl("sample56", combined_bed$Source)] <- "C. maculifer"
combined_bed$Species_Cory[grepl("sample4", combined_bed$Source)] <- "A. fuscoguttatus"
combined_bed$Species_Cory[grepl("sample39", combined_bed$Source)] <- "S. prionotus"
combined_bed$Species_Cory[grepl("sample63", combined_bed$Source)] <- "C. hastatus"
combined_bed$Species_Cory[grepl("sample28", combined_bed$Source)] <- "C. elegans"
combined_bed$Species_Cory[grepl("sample31", combined_bed$Source)] <- "C. paleatus"
combined_bed$Species_Cory[grepl("sample15", combined_bed$Source)] <- "C. haraldschultzei"
combined_bed$Species_Cory[grepl("sample57", combined_bed$Source)] <- "C. araguaiaensis"
combined_bed$Species_Cory[ (combined_bed$Source == "sample2")] <- "C. aeneus"

combined_bed <- merge(combined_bed, HT_elements, by = c("Species_Cory", "qry_id"), all.x = TRUE)

#what elements are interesting?
interesting <- combined_bed
interesting <- interesting %>% filter(!is.na(blastn_hit))
interesting <- interesting %>% group_by(repeat_id, altered_taxonomic_group) %>% add_count(name = "number_of_origins") #how many hits against a non teleost
interesting <- filter(interesting, altered_taxonomic_group != "Teleost")
interesting <- interesting %>% dplyr::select(Species_Cory, repeat_id, altered_taxonomic_group, number_of_origins)
interesting <- unique(interesting)

#HTT TEs of interest - in terms of the number of cory species they appear in and abundance are
#TE_00003492
#TE_00002042
#TE_00002955
#TE_00003910
#TE_00001723 - Mariner17






