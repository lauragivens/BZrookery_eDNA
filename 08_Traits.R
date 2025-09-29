library(taxize)
library(traits)
library(rfishbase)
library(dplyr)
library(stringr)
#library(rredlist)

dir_data <- '/Users/lauragivens/Desktop/R/BZrookery_eDNA/Rdata'
dir_man <- "/Volumes/Fuji/Mangroves"
dir_results <- "/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/cutadapt/results"
#dir_results <- "/Users/lauragivens/Desktop/MangroveResults"

#setwd(dir_results)

taxa_df <- readRDS(paste0(dir_results,"/tax.troph.combined.wordsize.df.rds"))

#two databases that we can pull from with rfishbase
#fishbase=marine fish species
#sealifebase=species consumed by marine fish/mammals

# ~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~ # fish (fishbase) # ~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~
# get chordate species  
fish.names <- taxa_df %>% 
  filter(Phylum=="Chordata" & !is.na(Species)) %>% 
  .$Species %>% unique() 
length(fish.names)

# select only valid species in fishbase   
valid.fish.names <- fish.names %>% rfishbase::validate_names(server="fishbase") 
valid.fish.names <- valid.fish.names[!is.na(valid.fish.names)]

# check lengths  
length(fish.names)
length(valid.fish.names)

# make new list of species that did not pass 
invalid.fish.names <- fish.names[!fish.names%in%valid.fish.names] 
length(invalid.fish.names)

# look for synonyms  
syn.fish.names <- rfishbase::synonyms(invalid.fish.names,server="fishbase")
View(syn.fish.names)

# get the ecology info for those species  
fish.troph <- rfishbase::ecology(species_list=valid.fish.names,
                                 server="fishbase"
)

#check if all species were returned
nrow(fish.troph)
length(valid.fish.names)

#add the ones that weren't to invalid.fish.names
invalid.fish.names <- c(invalid.fish.names,valid.fish.names[!valid.fish.names%in%fish.troph$Species])

# make a variable with the column names  
fishbase.colnames <- colnames(fish.troph)
fishbase.colnames[str_detect(fishbase.colnames,"Troph")] #which ones provide trophic data?  

# this only gets us the fish though, and only species that are valid in fishbase  
# what about non-fish taxa and anything that isn't annotated to species level?  
# ~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~ # non-fish (sealifebase) # ~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~
# get  species  
names <- taxa_df %>% 
  filter(Phylum!="Chordata" & !is.na(Species)) %>% 
  .$Species %>% unique() 
length(names)

# add the species that didn't pass the previous version in case they're available in sealifebase  
all.names <- c(invalid.fish.names,names)
length(all.names)

# select only valid species in fishbase   
valid.nonfish.names <- all.names %>% rfishbase::validate_names(server="sealifebase")
valid.nonfish.names <- valid.nonfish.names[!is.na(valid.nonfish.names)]

length(valid.nonfish.names)
invalid.nonfish.names <- all.names[!all.names%in%valid.nonfish.names]
length(invalid.nonfish.names)

# get the ecology info for those species  
nonfish.troph <- rfishbase::ecology(species_list=valid.nonfish.names,
                                    server="sealifebase"
)

# make a variable with the column names  
sealifebase.colnames <- colnames(nonfish.troph)
sealifebase.colnames[str_detect(sealifebase.colnames,"Troph")] #which ones provide trophic data?  

# ~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~ # combine # ~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~

# standardize the column names  
colnames(fish.troph) <- str_to_title(colnames(fish.troph))
colnames(nonfish.troph) <- str_to_title(colnames(nonfish.troph))
# are they the same? 
identical(colnames(fish.troph),colnames(nonfish.troph)) #no 

# combine the two dataframes with bind_rows, because that will fill any missing values with NA  
# include 'db' as an ID column so we can tell later what database the results came from  
traits.df <- bind_rows(fish.troph, nonfish.troph,.id="db") %>% filter(!duplicated(Species))

# merge with the taxa table  # merge with the taxa table  species_plantarum_binomials
taxa_traits_df <- merge((taxa_df %>% rownames_to_column("seqid")),traits.df,by="Species",all.x=TRUE)

taxa_traits_df <- relocate(taxa_traits_df,"Species",.after="Genus") 
#rownames(taxa_traits_df) <- taxa_traits_df$seqid
#taxa_traits_sub <- select(taxa_traits_df, -c(1:14))

# check numbers  
length(valid.fish.names)
length(valid.nonfish.names)

# is the size of the combined traits df from fishbase and sealifebase the same size as the two incoming?  
nrow(fish.troph)
nrow(nonfish.troph)
nrow(traits.df)

# is the size of the combined traits df the same size as the original taxa df?  
nrow(taxa_traits_df)
nrow(taxa_df)


# ~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~ # additional information # ~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~

# fisheries statistics  
# IUCN information (library(rredlist))

# ~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(taxa_traits_df,paste0(dir_results,"/taxa_traits_df.rds"))
#saveRDS(taxa_traits_sub,paste0(dir_results,"/taxa_traits_sub.rds"))
write.csv(taxa_traits_df,paste0(dir_results,"/taxa_traits_df.csv"))
#write.csv(taxa_traits_sub,paste0(dir_results,"/taxa_traits_sub.csv"))

save.image(paste0(dir_data,'/08_Traits.RData'))
