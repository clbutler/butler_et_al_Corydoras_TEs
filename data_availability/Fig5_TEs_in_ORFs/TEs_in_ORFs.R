#########################
#Figure 5
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
library(ape)
library(phytools)
library(l1ou)
library(rstudioapi)


#add ORF type - needed to only look at complete CDS 
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

setwd('./orf_type/')
import <- list.files(pattern = "*transdecoder.txt", full.names = TRUE)

ORFs <- lapply(import, read.table)
names(ORFs) <- import
combined_ORFs <- do.call("rbind", lapply(ORFs, as.data.frame))
combined_ORFs <- tibble::rownames_to_column(combined_ORFs, "Sample")
combined_ORFs$Sample <- sub('\\_.*', '', combined_ORFs$Sample)
combined_ORFs$Sample <- sub('/', '', combined_ORFs$Sample)
combined_ORFs$Sample <- sub('.', '', combined_ORFs$Sample)


combined_ORFs$V2<- sub('.*\\~', '', combined_ORFs$V2)
combined_ORFs <- rename(combined_ORFs, prot_id = V2 )
combined_ORFs$Transcript <- combined_ORFs$prot_id
combined_ORFs$Transcript <- gsub("\\..*", "", combined_ORFs$Transcript)
combined_ORFs$Sample <- gsub('s', 'S', combined_ORFs$Sample)

#highest log likelihood per ORF
combined_ORFs$logscore <- combined_ORFs$V6
combined_ORFs$logscore <- sub(".*\\=", "", combined_ORFs$logscore)
combined_ORFs$logscore <- as.numeric(combined_ORFs$logscore)
combined_ORFs <- combined_ORFs %>% group_by(Sample, Transcript) %>% filter(., logscore == max(logscore))


combined_ORFs <- filter(combined_ORFs, V4 == "type:complete") #only keep complete transcripts



combined_ORFs$Species[combined_ORFs$Sample == "Sample15"] <- "C. haraldschultzei"
combined_ORFs$Species[combined_ORFs$Sample == "Sample2"] <- "C. aeneus"
combined_ORFs$Species[combined_ORFs$Sample == "Sample28"] <- "C. elegans"
combined_ORFs$Species[combined_ORFs$Sample == "Sample31"] <- "C. paleatus"
combined_ORFs$Species[combined_ORFs$Sample == "Sample39"] <- "S. prionotus"
combined_ORFs$Species[combined_ORFs$Sample == "Sample4"] <- "A. fuscoguttatus"
combined_ORFs$Species[combined_ORFs$Sample == "Sample56"] <- "C. maculifer"
combined_ORFs$Species[combined_ORFs$Sample == "Sample57"] <- "C. araguaiaensis"
combined_ORFs$Species[combined_ORFs$Sample == "Sample63"] <- "C. hastatus"
#############################################################################

#add repeat data 
setwd(dirname(current_path ))
setwd('../RepeatMasker_outputs/Cory_TE/RM_trips/')

myfiles = list.files(pattern=".csv", full.names=TRUE) 
read_csv_filename <- function(filename){
  ret <- read.csv(filename, header = TRUE, stringsAsFactors=FALSE)
  ret$Source <- filename #EDIT
  ret
}
RMfiles = lapply(myfiles, read_csv_filename)
combined_RMfiles <-do.call(rbind, RMfiles)


combined_RMfiles$Species[grepl("Sample56_", combined_RMfiles$Source)] <- "C. maculifer"
combined_RMfiles$Species[grepl("sample4.", combined_RMfiles$Source)] <- "A. fuscoguttatus"
combined_RMfiles$Species[grepl("sample39.", combined_RMfiles$Source)] <- "S. prionotus"
combined_RMfiles$Species[grepl("sample63.", combined_RMfiles$Source)] <- "C. hastatus"
combined_RMfiles$Species[grepl("sample28.", combined_RMfiles$Source)] <- "C. elegans"
combined_RMfiles$Species[grepl("sample31.", combined_RMfiles$Source)] <- "C. paleatus"
combined_RMfiles$Species[combined_RMfiles$Source == "./noasterisk_cleaned6_sample2.fasta.out_RM_TRIPS.csv"] <- "C. aeneus"
combined_RMfiles$Species[grepl("sample15.", combined_RMfiles$Source)] <- "C. haraldschultzei"
combined_RMfiles$Species[grepl("sample57.", combined_RMfiles$Source)] <- "C. araguaiaensis"


combined_RMfiles <- rename(combined_RMfiles, Transcript = qry_id)

############################################
#How many TEs per2k overall

TEcount <- combined_RMfiles %>% group_by(Species) %>% add_count(name = "number_TEs_species")

setwd(dirname(current_path ))
transcript_lengths <- read.csv('../NCBI_transcriptomelengths.csv', stringsAsFactors=FALSE)

TEcount <- merge(TEcount, transcript_lengths, by = c("Species"))
TEcount <- mutate(TEcount, TEcountper2k = (number_TEs_species/Length)*2000)
TEcount <- TEcount %>% dplyr::select(Species, TEcountper2k) %>% unique()

############################################

combined_RMfiles <- merge(combined_ORFs, combined_RMfiles, by = c("Transcript", "Species"), all.x = TRUE)


#TEs in ORF

combined_RMfiles$ORF1 <- combined_RMfiles$V7
combined_RMfiles$ORF1 <- sub(".*\\:", "", combined_RMfiles$ORF1)
combined_RMfiles$ORF1 <- gsub("\\(.*", "", combined_RMfiles$ORF1)
combined_RMfiles$ORF2 <- combined_RMfiles$ORF1
combined_RMfiles$ORF1 <- gsub("\\-.*", "", combined_RMfiles$ORF1)
combined_RMfiles$ORF2 <- sub(".*\\-", "", combined_RMfiles$ORF2)

combined_RMfiles$ORF1 <- as.numeric(combined_RMfiles$ORF1)
combined_RMfiles$ORF2 <- as.numeric(combined_RMfiles$ORF2)

combined_RMfiles <- mutate(combined_RMfiles, ORF_length = (ORF2 - ORF1)+1)

combined_RMfiles$TEinORF[(!grepl('TE', combined_RMfiles$repeat_id))] <- 0

TE_ORF <- combined_RMfiles %>% dplyr::select(Transcript, Species, repeat_id, merged_qrystart, merged_qryend, ORF1, ORF2, ORF_length, TEinORF)

TE_ORF$overlap[(TE_ORF$merged_qrystart > TE_ORF$ORF1 & TE_ORF$merged_qryend < TE_ORF$ORF2)] <- 1
TE_ORF$overlap[is.na(TE_ORF$overlap)] <- 0

TE_ORF <- TE_ORF %>% dplyr::select(Transcript, Species, overlap, ORF_length)
TE_ORF <- TE_ORF %>% group_by(Species, Transcript) %>% mutate(., TEsinORF = sum(overlap))
TE_ORF <- TE_ORF %>% dplyr::select(Transcript, Species, ORF_length, TEsinORF)
TE_ORF <- unique(TE_ORF)

TE_ORF <- mutate(TE_ORF, TE_ORF_per2k = (TEsinORF/ORF_length)*2000)

#per species
TE_ORF_Species <- TE_ORF %>% group_by(Species) %>% mutate(., mean_TEsinORF_per2k = mean( TE_ORF_per2k))
TE_ORF_Species <- TE_ORF_Species %>% dplyr::select(Species, mean_TEsinORF_per2k) %>% unique()
TE_ORF_Species$Species <- factor(TE_ORF_Species$Species, levels = c('C. maculifer', 'A. fuscoguttatus', 'S. prionotus', 'C. hastatus', 'C. elegans', 'C. aeneus', 'C. haraldschultzei', 'C. paleatus', 'C. araguaiaensis'))


#based on haplotype data
#add ploidy level 

TE_ORF$ploidy[(TE_ORF$Species == "C. haraldschultzei")] <- "Polyploid"
TE_ORF$ploidy[(TE_ORF$Species == "C. aeneus")] <- "Polyploid"
TE_ORF$ploidy[(TE_ORF$Species == "C. elegans")] <- "Polyploid"
TE_ORF$ploidy[(TE_ORF$Species == "C. hastatus")] <- "Polyploid"
TE_ORF$ploidy[(TE_ORF$Species == "S. prionotus")] <- "Polyploid"
TE_ORF$ploidy[(TE_ORF$Species ==  "A. fuscoguttatus")] <- "Polyploid"

TE_ORF$ploidy[(TE_ORF$Species ==  "C. paleatus")] <- "Recent Polyploid"
TE_ORF$ploidy[(TE_ORF$Species ==  "C. araguaiaensis")] <- "Recent Polyploid"

TE_ORF$ploidy[(TE_ORF$Species ==  "C. maculifer")] <- "Non-Polyploid"


TE_ORF <- TE_ORF %>% group_by(ploidy) %>% mutate(., mean_TEsinORF_per2k = mean( TE_ORF_per2k))
std <- function(x) sd(x) / sqrt(length(x))
TE_ORF <- TE_ORF %>% group_by(ploidy) %>% mutate(., std_TEsinORF_per2k = std( TE_ORF_per2k))

#######stats######

kruskal.test(TE_ORF_per2k ~ ploidy, data = TE_ORF)
library('FSA')
dunnTest(TE_ORF_per2k ~ ploidy , data = TE_ORF, method = "bonferroni")

####################
#plotting
####################

TE_ORF_plot <- TE_ORF %>% dplyr::select(ploidy, mean_TEsinORF_per2k, std_TEsinORF_per2k)
TE_ORF_plot <- unique(TE_ORF_plot)

ggplot(TE_ORF_plot, aes(x = ploidy, y = mean_TEsinORF_per2k)) + geom_bar(stat = "identity", colour = "darkblue", fill = "cornflowerblue") + xlab ('\nPloidy Level')  +
         ylab("TE insertions into ORF per 2kb \n") + theme_bw() + 
   geom_errorbar(aes(ymin=mean_TEsinORF_per2k-std_TEsinORF_per2k, ymax=mean_TEsinORF_per2k+std_TEsinORF_per2k), width=.2) +
  theme(text = element_text(size = 25)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0.0, 0.1)


####################same again but standardise by total TEs#############   
TE_ORF_stand <- merge(TE_ORF, TEcount, by = "Species", all.x = TRUE)
TE_ORF_stand <- mutate(TE_ORF_stand, stand_TEs_ORF = (TE_ORF_per2k/TEcountper2k))

TE_ORF_stand <- TE_ORF_stand %>% group_by(ploidy) %>%  mutate(., mean_stand = mean(stand_TEs_ORF))    
TE_ORF_stand <- TE_ORF_stand %>% group_by(ploidy) %>%  mutate(., std_stand = std(stand_TEs_ORF))   

TE_ORF_stand_perspecies <- TE_ORF_stand %>% group_by(Species) %>%  mutate(., mean_stand = mean(stand_TEs_ORF))    
TE_ORF_stand_perspecies <- TE_ORF_stand_perspecies %>% dplyr::select(Species, mean_stand) %>% unique()

kruskal.test(stand_TEs_ORF ~ ploidy, data = TE_ORF_stand)
dunnTest(stand_TEs_ORF ~ ploidy , data = TE_ORF_stand, method = "bonferroni")



TE_ORF_stand_plot <- TE_ORF_stand %>% dplyr::select(ploidy, mean_stand, std_stand)
TE_ORF_stand_plot <- unique(TE_ORF_stand_plot)

#ggplot(TE_ORF_stand_plot, aes(x = ploidy, y = mean_stand)) + geom_bar(stat = "identity", colour = "darkblue", fill = "cornflowerblue") + xlab ('\nPolidy Level')  +
  #ylab("TE insertions into ORF_stand per 2kb (standardised)\n") + theme_bw() + 
 
  # theme(text = element_text(size = 25)) +   geom_errorbar(aes(ymin= mean_stand - std_stand, ymax= mean_stand + std_stand), width=.2)

##########################
#based on snp data
#add ploidy level 

TE_ORF$ploidy[(TE_ORF$Species == "C. haraldschultzei")] <- "Polyploid"
TE_ORF$ploidy[(TE_ORF$Species == "C. aeneus")] <- "Polyploid"


TE_ORF$ploidy[(TE_ORF$Species ==  "C. paleatus")] <- "Recent Polyploid"
TE_ORF$ploidy[(TE_ORF$Species ==  "C. araguaiaensis")] <- "Recent Polyploid"

TE_ORF$ploidy[(TE_ORF$Species ==  "C. maculifer")] <- "Non-Polyploid"
TE_ORF$ploidy[(TE_ORF$Species == "C. elegans")] <- "Non-Polyploid"
TE_ORF$ploidy[(TE_ORF$Species == "C. hastatus")] <- "Non-Polyploid"
TE_ORF$ploidy[(TE_ORF$Species == "S. prionotus")] <- "Non-Polyploid"
TE_ORF$ploidy[(TE_ORF$Species ==  "A. fuscoguttatus")] <- "Non-Polyploid"


TE_ORF <- TE_ORF %>% group_by(ploidy) %>% mutate(., mean_TEsinORF_per2k = mean( TE_ORF_per2k))
std <- function(x) sd(x) / sqrt(length(x))
TE_ORF <- TE_ORF %>% group_by(ploidy) %>% mutate(., std_TEsinORF_per2k = std( TE_ORF_per2k))

#######stats######

kruskal.test(TE_ORF_per2k ~ ploidy, data = TE_ORF)
library('FSA')
dunnTest(TE_ORF_per2k ~ ploidy , data = TE_ORF, method = "bonferroni")

####################
#plotting
####################

TE_ORF_plot <- TE_ORF %>% dplyr::select(ploidy, mean_TEsinORF_per2k, std_TEsinORF_per2k)
TE_ORF_plot <- unique(TE_ORF_plot)

ggplot(TE_ORF_plot, aes(x = ploidy, y = mean_TEsinORF_per2k)) + geom_bar(stat = "identity", colour = "darkblue", fill = "cornflowerblue") + xlab ('\nPloidy Level')  +
  ylab("TE insertions into ORF per 2kb \n") + theme_bw() + 
  geom_errorbar(aes(ymin=mean_TEsinORF_per2k-std_TEsinORF_per2k, ymax=mean_TEsinORF_per2k+std_TEsinORF_per2k), width=.2) +
  theme(text = element_text(size = 25)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0.0, 0.1) 


####################same again but standardise by total TEs#############   
TE_ORF_stand <- merge(TE_ORF, TEcount, by = "Species", all.x = TRUE)
TE_ORF_stand <- mutate(TE_ORF_stand, stand_TEs_ORF = (TE_ORF_per2k/TEcountper2k))





TE_ORF_stand <- TE_ORF_stand %>% group_by(ploidy) %>%  mutate(., mean_stand = mean(stand_TEs_ORF))    
TE_ORF_stand <- TE_ORF_stand %>% group_by(ploidy) %>%  mutate(., std_stand = std(stand_TEs_ORF))   

TE_ORF_stand_perspecies <- TE_ORF_stand %>% group_by(Species) %>%  mutate(., mean_stand = mean(stand_TEs_ORF))    
TE_ORF_stand_perspecies <- TE_ORF_stand_perspecies %>% dplyr::select(Species, mean_stand) %>% unique()

kruskal.test(stand_TEs_ORF ~ ploidy, data = TE_ORF_stand)
dunnTest(stand_TEs_ORF ~ ploidy , data = TE_ORF_stand, method = "bonferroni")



#TE_ORF_stand_plot <- TE_ORF_stand %>% dplyr::select(ploidy, mean_stand, std_stand)
#TE_ORF_stand_plot <- unique(TE_ORF_stand_plot)

#ggplot(TE_ORF_stand_plot, aes(x = ploidy, y = mean_stand)) + geom_bar(stat = "identity", colour = "darkblue", fill = "cornflowerblue") + xlab ('\nPolidy Level')  +
#ylab("TE insertions into ORF_stand per 2kb (standardised)\n") + theme_bw() + 

# theme(text = element_text(size = 25)) +   geom_errorbar(aes(ymin= mean_stand - std_stand, ymax= mean_stand + std_stand), width=.2)
