#########################
#Figure 2 Bubbleplot
###########################

#RESET
rm(list = ls())
par(mfrow=c(1,1))

#LOAD LIBRARIES 
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(rstudioapi)



######
#transcriptome
######

current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

#set working directory 
setwd('../RepeatMasker_outputs/Cory_TE/RM_trips/')

myfiles = list.files(pattern=".csv", full.names=TRUE) 
read_csv_filename <- function(filename){
  ret <- read.csv(filename, header = TRUE, stringsAsFactors=FALSE)
  ret$Source <- filename #EDIT
  ret
}
RMfiles = lapply(myfiles, read_csv_filename)
transcriptome <-do.call(rbind, RMfiles)

setwd(dirname(current_path ))
transcript_lengths <- read.csv('../NCBI_transcriptomelengths.csv', stringsAsFactors=FALSE)

unique(transcriptome$matching_class)
transcriptome$type[grepl("Gypsy", transcriptome$matching_class)] <- "Gypsy"
transcriptome$type[grepl("Copia", transcriptome$matching_class)] <- "Copia"
transcriptome$type[grepl("LINE", transcriptome$matching_class)] <- "LINE"
transcriptome$type[grepl("TcMar", transcriptome$matching_class)] <- "Mariner"
transcriptome$type[grepl("BEL", transcriptome$matching_class)] <- "BEL"
transcriptome$type[grepl("Mutator", transcriptome$matching_class)] <- "Mutator"
transcriptome$type[grepl("ERV", transcriptome$matching_class)] <- "ERV"
transcriptome$type[grepl("hAT", transcriptome$matching_class)] <- "hAT"
transcriptome$type[grepl("Harbinger", transcriptome$matching_class)] <- "Harbinger"
transcriptome$type[grepl("CACTA", transcriptome$matching_class)] <- "CACTA"
transcriptome$type[grepl("Helitron", transcriptome$matching_class)] <- "Helitron"
transcriptome$type[grepl("SINE", transcriptome$matching_class)] <- "SINE"
transcriptome$type[grepl("PiggyBac", transcriptome$matching_class)] <- "PiggyBac"
transcriptome$type[grepl("DIRS", transcriptome$matching_class)] <- "DIRS"
transcriptome$type[grepl("PLE", transcriptome$matching_class)] <- "Penelope"
transcriptome$type[grepl("P_nMITE", transcriptome$matching_class)] <- "P"

transcriptome$Class <- transcriptome$matching_class
transcriptome$Class[grepl("BEL", transcriptome$Class)] <- "LTR"
transcriptome$Class[grepl("Copia", transcriptome$Class)] <- "LTR"
transcriptome$Class[grepl("ERV", transcriptome$Class)] <- "LTR"
transcriptome$Class[grepl("Gypsy", transcriptome$Class)] <- "LTR"
transcriptome$Class[grepl("DIRS", transcriptome$Class)] <- "LTR"
transcriptome$Class[grepl("LINE", transcriptome$Class)] <- "LINE"
transcriptome$Class[grepl("PLE", transcriptome$Class)] <- "PLE"
transcriptome$Class[grepl("SINE", transcriptome$Class)] <- "SINE"
transcriptome$Class[grepl("CACTA", transcriptome$Class)] <- "DNA"
transcriptome$Class[grepl("Harbinger", transcriptome$Class)] <- "DNA"
transcriptome$Class[grepl("hAT", transcriptome$Class)] <- "DNA"
transcriptome$Class[grepl("Merlin", transcriptome$Class)] <- "DNA"
transcriptome$Class[grepl("Mutator", transcriptome$Class)] <- "DNA"
transcriptome$Class[grepl("DNA_P", transcriptome$Class)] <- "DNA"
transcriptome$Class[grepl("PiggyBac", transcriptome$Class)] <- "DNA"
transcriptome$Class[grepl("TcMar", transcriptome$Class)] <- "DNA"
transcriptome$Class[grepl("Helitron", transcriptome$Class)] <- "Helitron"
transcriptome$Class[grepl("EnSpm", transcriptome$Class)] <- "DNA"
transcriptome$Class[grepl("Dada", transcriptome$Class)] <- "DNA"
transcriptome$Class[grepl("Maverick", transcriptome$Class)] <- "Maverick"
transcriptome$Class[grepl("DNA", transcriptome$Class)] <- "DNA"


transcriptome$Species[grepl("Sample56_", transcriptome$Source)] <- "C. maculifer"
transcriptome$Species[grepl("sample4.", transcriptome$Source)] <- "A. fuscoguttatus"
transcriptome$Species[grepl("sample39.", transcriptome$Source)] <- "S. prionotus"
transcriptome$Species[grepl("sample63.", transcriptome$Source)] <- "C. hastatus"
transcriptome$Species[grepl("sample28.", transcriptome$Source)] <- "C. elegans"
transcriptome$Species[grepl("sample31.", transcriptome$Source)] <- "C. paleatus"
transcriptome$Species[transcriptome$Source == "./noasterisk_cleaned6_sample2.fasta.out_RM_TRIPS.csv"] <- "C. aeneus"
transcriptome$Species[grepl("sample15.", transcriptome$Source)] <- "C. haraldschultzei"
transcriptome$Species[grepl("sample57.", transcriptome$Source)] <- "C. araguaiaensis"
transcriptome$Species[grepl("pCoco", transcriptome$Source)] <- "I. punctatus"

transcriptome <- merge(transcriptome, transcript_lengths, by = "Species")

transcriptome <- transcriptome %>% group_by(Species, type) %>% mutate(., type_length = sum(mergedfraglength))
table <- transcriptome %>% group_by(Species, type) %>% summarise(., type_perc_length = (type_length/Length)*100)
table <- unique(table)

table <- tidyr::spread(table, type, type_perc_length)

table$Species <- factor(table$Species, levels = c('I. punctatus', 'C. maculifer', 'A. fuscoguttatus', 'S. prionotus', 'C. hastatus', 'C. elegans', 'C. aeneus', 'C. haraldschultzei', 'C. paleatus', 'C. araguaiaensis'))

table <- table %>% dplyr::select('BEL', 'Copia', 'ERV', 'Gypsy', 'DIRS', 'LINE', 'Penelope', 'SINE', 'CACTA', 'Harbinger', 'hAT', 'Mariner', 'Mutator', 'PiggyBac', 'Helitron')


#write.table(table, "summary.table.txt" )

###### Bubbleplots

###### bubble plots ########
bubble <- table

bubble <- tidyr::gather(bubble, TE_type, Abundance, BEL:Helitron)


#Create new column grouping TE families by Curcio & Derbyshire mechanism of transposition

bubble$Mechanism[grepl("SINE", bubble$TE_type)] <- "Non-autonomous Retrotransposon"
bubble$Mechanism[grepl("DIRS", bubble$TE_type)] <- "Y Retrotransposon"
bubble$Mechanism[grepl("Penelope", bubble$TE_type)] <- "Target Prime Retrotransposon"
bubble$Mechanism[grepl("LINE", bubble$TE_type)] <- "Target Prime Retrotransposon"
bubble$Mechanism[grepl("BEL", bubble$TE_type)] <- "LTR Retrotransposon"
bubble$Mechanism[grepl("Copia", bubble$TE_type)] <- "LTR Retrotransposon"
bubble$Mechanism[grepl("ERV", bubble$TE_type)] <- "LTR Retrotransposon"
bubble$Mechanism[grepl("Gypsy", bubble$TE_type)] <- "LTR Retrotransposon"
bubble$Mechanism[grepl("Helitron", bubble$TE_type)] <- "Rolling Circle DNA transposon"
bubble$Mechanism <- bubble$Mechanism %>% tidyr::replace_na("Terminal Inverted Repeat DNA transposon") #Everything else is a TIR DNA transposon

#Create new column grouping mechanisms into classic class I retrotransposons and class II DNA transposons
bubble$Class[grepl("Retrotransposon", bubble$Mechanism)] <- "Class I Retrotransposon"
bubble$Class[grepl("DNA ", bubble$Mechanism)] <- "Class II DNA Transposon"

#Rename again to abbreviations
bubble$Mechanism[grepl('LTR', bubble$Mechanism)] <- 'DDE'
bubble$Mechanism[grepl('Terminal', bubble$Mechanism)] <- 'DDE'
bubble$Mechanism[grepl('Target', bubble$Mechanism)] <- 'TP'
bubble$Mechanism[grepl('Non-autonomous', bubble$Mechanism)] <- 'NA'
bubble$Mechanism[grepl('Self', bubble$Mechanism)] <- 'SS'
bubble$Mechanism[grepl('Y', bubble$Mechanism)] <- 'Y'
bubble$Mechanism[grepl('Rolling', bubble$Mechanism)] <- 'Y2'



#turn copies into a discrete scale
bubble$grouped_copies[bubble$Abundance < 0.1] <- '< 0.1'
bubble$grouped_copies[bubble$Abundance >=0.1 & bubble$Abundance < 0.2] <- '0.1-0.2'
bubble$grouped_copies[bubble$Abundance >=0.2 & bubble$Abundance < 0.3] <- '0.2-0.3'
bubble$grouped_copies[bubble$Abundance >=0.3 & bubble$Abundance < 0.4] <- '0.3-0.4'
bubble$grouped_copies[bubble$Abundance >=0.4 & bubble$Abundance < 0.5] <- '0.4-0.5'
bubble$grouped_copies[bubble$Abundance >=0.5 & bubble$Abundance < 1.0] <- '0.5-1.0'
bubble$grouped_copies[bubble$Abundance >=1.0 & bubble$Abundance < 1.5] <- '1.0-1.5'
bubble$grouped_copies[bubble$Abundance >=1.5 & bubble$Abundance < 2.0] <- '1.5-2.0'
bubble$grouped_copies[bubble$Abundance >=2.0] <- '> 2.0'


#Change order of factors.
bubble$grouped_copies <- factor(bubble$grouped_copies, levels = c("< 0.1","0.1-0.2", "0.2-0.3", "0.3-0.4", "0.5-1.0", "0.4-0.5", "1.0-1.5", "1.5-2.0", "60.0 - 80.0", "> 2.0"))
bubble$TE_type <- factor(bubble$TE_type, levels = c('BEL', 'Copia', 'ERV', 'Gypsy', 'DIRS', 'LINE', 'Penelope', 'SINE', 'CACTA', 'Harbinger', 'hAT', 'Mariner', 'Mutator', 'PiggyBac', 'Helitron'))
bubble$Mechanism <- factor(bubble$Mechanism, levels = c("DDE", "Y", "TP" , "Non-autonomous", "SS", "Y2"))


setwd('/Users/chrisbutler/Documents/chapter4/summary_table/')

ggplot(bubble, aes(x = Species, y = TE_type)) +
geom_point(aes(fill = grouped_copies), colour = "black", pch = 21, size = 7) + xlab ('\nSpecies') + theme_light() + ylab("TE Family\n")  + 
  scale_fill_brewer(palette = "Reds") +
  labs(fill = 'TE Abundance %') + facet_grid(Mechanism ~ Class, scale = 'free', space = "free_y") +
theme(axis.text.x=element_text(angle = 290)) + theme(text = element_text(size = 17))


  