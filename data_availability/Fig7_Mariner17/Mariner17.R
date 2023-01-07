#######
#Figure 7 - Mariner17
#######

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


################
#extract all Mariner17 copies - not just 'full length' transcripts and create .bed files
################

current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

#set working directory 
setwd('../RepeatMasker_outputs/Cory_TE/RM_trips/')
#load custom functions
read_csv_filename <- function(filename){
  ret <- read.csv(filename, fill = TRUE, row.names=NULL, stringsAsFactors = FALSE)
  ret$Source <- filename #EDIT
  ret
}

myfiles = list.files(pattern=".csv", full.names=TRUE) 

RMfiles = lapply(myfiles, read_csv_filename)
transcriptome <-do.call(rbind, RMfiles)
transcriptome <- filter(transcriptome, Source != "./noasterisk_IpCoco_mRNA.fa.out_RM_TRIPS.csv") #remove channel catfish

noisonosim <- transcriptome

######Adding TE type and class column#####
noisonosim$TE_type[grepl("TE_00001723", noisonosim$repeat_id)] <- "Mariner-17_DR"

noisonosim$Source <- gsub('.fasta.out_RM_TRIPS.csv', '', noisonosim$Source)
noisonosim$Source <- gsub('.*\\_', '', noisonosim$Source)
noisonosim$Source <- gsub('S', 's', noisonosim$Source)
noisonosim$Source[(noisonosim$Source == 'final')] <- 'sample56'



#filter just Mariner
Mariner17 <- noisonosim %>% filter(TE_type == 'Mariner-17_DR')


filter_list <- c("sample56", "sample15", "sample2", "sample28", 'sample31', 'sample39', 'sample4', 'sample57', 'sample63')

#create vector with masked % per species
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )
setwd('./Mariner17_blast/')
x <- split(Mariner17, Mariner17$Source)

for (i in filter_list) {
  x[[i]] <- dplyr::select(x[[i]],'qry_id', 'merged_qrystart', 'merged_qryend', 'repeat_id')
}



for (i in filter_list) {
  write.table(x[[i]], file = paste0("TE_00001723", i ,".bed"), col.names = FALSE, sep = '\t', quote = FALSE, row.names = FALSE)
}

################
#extract all Mariner17 copies - do the same but for the danio searching
################

current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

#set working directory 
setwd('../RepeatMasker_outputs/Danio/RMtrips/')
#load custom functions
read_csv_filename <- function(filename){
  ret <- read.csv(filename, fill = TRUE, row.names=NULL, stringsAsFactors = FALSE)
  ret$Source <- filename #EDIT
  ret
}

myfiles = list.files(pattern=".csv", full.names=TRUE) 

RMfiles = lapply(myfiles, read_csv_filename)
transcriptome <-do.call(rbind, RMfiles)
transcriptome <- filter(transcriptome, Source != "./noasterisk_IpCoco_mRNA.fa.out_RM_TRIPS.csv") #remove channel catfish

noisonosim <- transcriptome

######Adding TE type and class column#####


noisonosim$Source <- gsub('.fasta.out_RM_TRIPS.csv', '', noisonosim$Source)
noisonosim$Source <- gsub('.*\\_', '', noisonosim$Source)
noisonosim$Source <- gsub('S', 's', noisonosim$Source)
noisonosim$Source[(noisonosim$Source == 'final')] <- 'sample56'



#filter just Mariner
Mariner17 <- noisonosim %>% filter(repeat_id == 'Mariner-17_DR')


filter_list <- c("sample56", "sample15", "sample2", "sample28", 'sample31', 'sample39', 'sample4', 'sample57', 'sample63')

#create vector with masked % per species
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )
setwd('./Mariner17_blast/')
x <- split(Mariner17, Mariner17$Source)

for (i in filter_list) {
  x[[i]] <- dplyr::select(x[[i]],'qry_id', 'merged_qrystart', 'merged_qryend', 'repeat_id')
}

for (i in filter_list) {
  write.table(x[[i]], file = paste0("mariner17_DR", i ,".bed"), col.names = FALSE, sep = '\t', quote = FALSE, row.names = FALSE)
}

#######
#SEQUENCE DIVERGENCE PLOTS - FIGURE 7
#######
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )
setwd('./Xenopus/')

#TC1-5 element within Xenopus 

xenopus <- read.csv('noasterisk_GCF_000004195.4_UCB_Xtro_10.0_rna.fa.out_RM_TRIPS.csv')
xenopus <- filter(xenopus, repeat_id == 'Tc1-5_Xt')
mean(xenopus$perc_div)
std <- function(x) sd(x) / sqrt(length(x)) # Create own function
std(xenopus$perc_div)
xenopus$species <- 'Xenopus'

xenopus <- xenopus %>% dplyr::select(repeat_id, perc_div, species)


#set working directory to the HT rerun folder
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )
setwd('./HT_RMtrips/')


myfiles = list.files(pattern=".csv", full.names=TRUE) 
read_csv_filename <- function(filename){
  ret <- read.csv(filename, header = TRUE, stringsAsFactors=FALSE)
  ret$Source <- filename #EDIT
  ret
}
RMfiles = lapply(myfiles, read_csv_filename)
combined_RMfiles <-do.call(rbind, RMfiles)

cory_tc1_5 <- filter(combined_RMfiles, repeat_id == "Mariner-17_Cory")
sep_species <- rename(cory_tc1_5, species = Source)

cory_tc1_5$species <- 'Corydoradinae'
cory_tc1_5 <- cory_tc1_5 %>% dplyr::select(repeat_id, perc_div, species)

single_element <- rbind(xenopus, cory_tc1_5)

single_element %>% group_by(species) %>% summarise(., mean = mean(perc_div))
single_element %>% group_by(species) %>% summarise(., std = std(perc_div))

ggplot(single_element, aes(x = perc_div, fill = species)) + geom_density(alpha = 0.5) + xlab('\n% Divergence from Consensus') + ylab('Density\n ') + 
  theme_bw() + theme(text = element_text(size = 25)) + labs(fill = "Species")  + scale_x_continuous(limits = c(0, 40))

stats <- single_element
stats <- single_element %>% group_by(species) %>% mutate(., meandiv = mean(perc_div))


#stats
#ggplot(stats, aes(x = perc_div, fill = species)) + geom_histogram() + facet_wrap(~species, scales = "free_y")
wilcox.test(perc_div ~ species, single_element)

#seperate Cory species
ggplot(sep_species, aes(x = perc_div)) + geom_density() + facet_wrap(~species, scales = "free_y") #true of all Cory species - not driven by one


#compare against other mariner elements
all_mariner <- combined_RMfiles
all_mariner$sequencetype <- "All Mariner"
all_mariner$sequencetype[(all_mariner$repeat_id == "Mariner-17_Cory")] <- "Mariner-17_Cory"
all_mariner <- filter(all_mariner, repeat_id != 'Tc1-5_Xt')
all_mariner <- filter(all_mariner, repeat_id != 'Mariner-12_EL')

all_mariner <- all_mariner %>% dplyr::select('Source', 'perc_div', 'sequencetype')


ggplot(all_mariner, aes(x = perc_div, fill = sequencetype)) + geom_density(alpha = 0.5) + xlab('\n% Divergence from Consensus') + ylab('Density\n ') + 
  theme_bw() + theme(text = element_text(size = 25)) + labs(fill = "TE") 



#stats
all_mariner <- all_mariner %>% group_by(sequencetype) %>%  mutate(., mean = mean(perc_div))
all_mariner <- all_mariner %>% group_by(sequencetype) %>%  mutate(., SD = std(perc_div))
unique(all_mariner$mean)
unique(all_mariner$SD)
wilcox.test(perc_div ~ sequencetype, all_mariner)



