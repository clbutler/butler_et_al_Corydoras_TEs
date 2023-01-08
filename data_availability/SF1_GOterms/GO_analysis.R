#########################
#GO term enrichment 
###########################

#RESET
rm(list = ls())
par(mfrow=c(1,1))

#LOAD LIBRARIES 
library(dplyr)
library(ggplot2)
library(splitstackshape)
library(reshape2)
library(forcats)
library(rstudioapi)
library(stringr)


#
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )
setwd('./trinnotate_NCBI/')

#load input files
import <- list.files(pattern="*.csv", full.names = TRUE)
read_csv_filename <- function(filename){
  ret <- read.csv(filename, header = TRUE, stringsAsFactors=FALSE)
  ret$Source <- filename #EDIT
  ret
}
RM <- lapply(import, read_csv_filename) 
names(RM) <- import
combined_RM <-do.call("rbind", lapply(RM, as.data.frame))

# clean source column
combined_RM$Source <- gsub('.fasta_trinotate.csv', '', combined_RM$Source)
combined_RM$Source <- gsub('.*_', '', combined_RM$Source)

combined_RM$Source[(combined_RM$Source == "final")] <- "sample56"

#remove transcripts with zero protein blast hit
combined_RM <- filter(combined_RM, prot_id != '.')
GOterms <- combined_RM

#load TE RM_trips
setwd('../../RepeatMasker_outputs/Cory_TE/RM_trips/')

myfiles = list.files(pattern=".csv", full.names=TRUE) 
read_csv_filename <- function(filename){
  ret <- read.csv(filename, header = TRUE, stringsAsFactors=FALSE)
  ret$Source <- filename #EDIT
  ret
}
RMfiles = lapply(myfiles, read_csv_filename)
combined_RMfiles <-do.call(rbind, RMfiles)

combined_RMfiles$Source[(combined_RMfiles$Source == "./noasterisk_denovo_Sample56_NCBI_3_final.fasta.out_RM_TRIPS.csv")] <- "sample56"

combined_RMfiles$Source <- gsub('.fasta.out_RM_TRIPS.csv', '', combined_RMfiles$Source)
combined_RMfiles$Source <- gsub('/noasterisk_', '', combined_RMfiles$Source)
combined_RMfiles$Source <- gsub('.*_', '', combined_RMfiles$Source)

#merge
GOterms <- rename(GOterms, qry_id = transcript_id)


GOterms <- cSplit(GOterms, "gene_ontology_blast", "`") 

GOterms <- GOterms %>% dplyr::select(qry_id, Source, gene_ontology_blast_001, gene_ontology_blast_002, gene_ontology_blast_003, gene_ontology_blast_004, gene_ontology_blast_005, gene_ontology_blast_006, gene_ontology_blast_007, gene_ontology_blast_008, gene_ontology_blast_009, gene_ontology_blast_010, gene_ontology_blast_011, gene_ontology_blast_012, gene_ontology_blast_013, gene_ontology_blast_014, gene_ontology_blast_015, gene_ontology_blast_016, gene_ontology_blast_017, gene_ontology_blast_018, gene_ontology_blast_019, gene_ontology_blast_020, gene_ontology_blast_021)
GOterms <- melt(GOterms, id = c('qry_id', 'Source'))
GOterms <- filter(GOterms, str_detect(value, 'biological')) #only interested in biological GO terms
GOterms$value <- gsub('\\^', '`', GOterms$value)
GOterms <- cSplit(GOterms, "value", "`")


#now to get plotting

GOterms <- GOterms %>% group_by(Source, value_3) %>% add_count(name = "species_GO_occur") #get total "GO" terms to use as 'reference'

GOterms <- merge(GOterms, combined_RMfiles, by = c('qry_id', 'Source'), all.x = TRUE) #add TE data

#now only select transcripts with TEs in
GOterms <- GOterms %>% filter(!is.na(repeat_id))

GOterms <- GOterms %>% group_by(Source, value_3) %>% add_count(name = "TE_GO_occur") #get "GO" terms associated with TEs

GOterms <- mutate(GOterms, GO_TEenrich = TE_GO_occur/species_GO_occur )

#filter GO terms that are present very rarely 
GOterms <- filter(GOterms, species_GO_occur > 5)

#filter TE based terms - reflects the TE ORF
filter_list <- c("transposition, DNA-mediated",  
"transposition, RNA-mediated",
"reverse transcription involved in RNA-mediated transposition")  

GOterms <- filter(GOterms, !value_3 %in% filter_list)


GOterms <- GOterms %>% dplyr::select(Source, value_3, GO_TEenrich)
GOterms <- unique(GOterms)

GOterms <- GOterms %>% ungroup() %>% group_by(Source) %>% top_n(., n = 10) #top 10 most enriched GO terms per species

##############plot#################
GOterms$Species[grepl("sample56", GOterms$Source)] <- "C. maculifer"
GOterms$Species[grepl("sample4", GOterms$Source)] <- "A. fuscoguttatus"
GOterms$Species[grepl("sample39", GOterms$Source)] <- "S. prionotus"
GOterms$Species[grepl("sample63", GOterms$Source)] <- "C. hastatus"
GOterms$Species[grepl("sample28", GOterms$Source)] <- "C. elegans"
GOterms$Species[grepl("sample31", GOterms$Source)] <- "C. paleatus"
GOterms$Species[grepl("sample15", GOterms$Source)] <- "C. haraldschultzei"
GOterms$Species[grepl("sample57", GOterms$Source)] <- "C. araguaiaensis"
GOterms$Species[ (GOterms$Source == "sample2")] <- "C. aeneus"

GOterms$present <- "Yes"

GOterms <- GOterms %>% group_by(value_3) %>% add_count(name = "occurences")

GOterms$Species <- factor(GOterms$Species, levels = c('C. maculifer', "A. fuscoguttatus", "S. prionotus", 'C. hastatus', 'C. elegans', 'C. aeneus', 'C. haraldschultzei','C. paleatus',  'C. araguaiaensis'))

GOterms_histogram <- GOterms %>% dplyr::select(value_3, occurences)
GOterms_histogram <- unique(GOterms_histogram)

ggplot(GOterms_histogram, aes(x = occurences)) + geom_histogram(position = "identity") + xlab("\n Frequency of enriched GO terms across different Corydoras species") +
  ylab('Count\n') + theme_light() + theme(text = element_text(size = 20)) 

GOterms_histogram <- GOterms_histogram %>% group_by(occurences) %>% add_count(name = 'frequency')
GOterms_histogram <- GOterms_histogram %>% dplyr::select(occurences, frequency)
GOterms_histogram <- unique(GOterms_histogram)


ggplot(GOterms_histogram, aes(x = occurences, y = frequency)) + geom_bar(stat = "identity") + xlab("\n Frequency of enriched GO terms across different Corydoras species") +
  ylab('Count\n') + theme_light() + theme(text = element_text(size = 20)) 


ggplot(GOterms, aes(x = Species, y = fct_reorder(value_3, occurences))) +
  geom_point(aes(fill = present), colour = "black", pch = 22, size = 3) + xlab ('\nSpecies') + theme_light() + ylab("Biological GO Term\n")  +
  theme(axis.text.x=element_text(angle = 290)) + theme(text = element_text(size = 10)) + theme(legend.position = 'NONE')


