library(phyloseq)
library(vegan)
library(tidyverse)

dir_data <- '/Users/lauragivens/Desktop/R/BZrookery_eDNA/Rdata'
dir_man <- "/Volumes/Fuji/Mangroves"
#dir_results <- "/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/cutadapt/results"
dir_results <- "/Users/lauragivens/Desktop/MangroveResults"
setwd(dir_results) 


bold_output <- readxl::read_xlsx("dada2-uniqueseqs_identification_result.xlsx") 
bold_output[bold_output=='no-match'] <- NA
bold_selected <- bold_output %>% mutate(.,seqid=word(id,sep=";"),.before=1) %>% as.data.frame()
rownames(bold_selected) <- bold_selected$seqid

taxa_nt_sub <- readRDS(paste0(dir_results,"/taxtable.v2.nt.rds"))
taxa_nt_troph <- readRDS(paste0(dir_results,"/taxtable.v2.nt.wtroph.rds"))

bold_sub <- bold_selected %>% filter(!is.na(Phylum) & Phylum!="IncompleteTaxonomy")

View(bold_sub)

bold.tomerge <- (bold_sub[c(1,3:9,11)]  %>% 
                   mutate(db=1))
nt.tomerge <- (taxa_nt_sub[2:7] %>% 
                 rownames_to_column(var="seqid") %>%
                 mutate(db=2,pct_identity=NA,records=NA))

combined_tax <- rbind(bold.tomerge,
                      nt.tomerge
                      ) %>% mutate(db=case_when(db==1 ~ "BOLD",
                                                db==2 ~ "nt"))

# select duplicates
# combined_tax %>% 
#  group_by(seqid) %>% filter(n() == 2) %>%
#  View()

check.tax <- combined_tax %>% 
  group_by(seqid) %>% filter(n() == 2) %>% 
  arrange(seqid) %>%
  # mutate(across(everything(),!duplicated())) %>% 
  mutate(across(where(is.character), ~ duplicated(.x), .names = "{.col}_dup"),
         across(where(is.logical), ~ case_when( (.x==TRUE) ~ .x))
  ) 

#check.tax %>% 
#  group_by(seqid) %>%
#  summarize(across(where(is.character),~ unique(.x))) %>%
#  View()
#combined_tax %>% 
#  group_by(seqid) %>% 
#  filter(n()==2) %>% 
#  summarise(across(where(is.character), ~ na.omit(.x))) %>% View()

check.tax %>% 
  group_by(seqid) %>% filter(n() == 2) %>% 
  arrange(seqid) %>%
  mutate_at() #case_when dup == TRUE, .default=NA? 
  View()

check.tax %>%   
  group_by(seqid) %>% filter(n() == 2) %>% 
  arrange(seqid) %>%
  #mutate(Species = map_chr(Species, ~toString(sort(str_split(.x, " ")[[1]])))) %>%
  #https://stackoverflow.com/questions/60322712/merge-rows-containing-similar-strings-using-dplyr
  summarize(Phylum=first(Phylum),
            Class=first(Class),
            Order=first(Order),
            Family=first(Family),
            Genus=first(Genus),
            Species=first(Species)
            ) %>%
  View()


#this will combine rows that are the same but does nothing for rows where Species != Species
check.tax %>%   
  group_by(seqid,Phylum,Class,Order,Family,Genus,Species) %>% 
  summarize(Phylum=first(Phylum),
            Class=first(Class),
            Order=first(Order),
            Family=first(Family),
            Genus=first(Genus),
            Species=first(Species)
  ) %>%
  View()

#check.tax %>%   
#  group_by(seqid,Phylum,Class,Order,Family,Genus,Species) %>% 
#  summarize_all(na.omit) %>%
#  View()
#check.tax %>%   
#  group_by(seqid) %>% 
#  summarize_all(na.omit) %>%
#  View()
check.tax %>%   
  mutate(Class=case_when(Class=="Actinopteri" ~ "Actinopterygii",.default=Class)) %>%
  group_by(seqid) %>% 
  summarise_all(funs(toString(unique(na.omit(.)))))  %>% #unique for duplicated like col4
  #mutate_all(funs(na.locf(., na.rm = FALSE, fromLast = FALSE)))%>%filter(row_number()==n()) %>%
  #summarise_all(funs(toString(na.omit(.)))) %>%
  View()

# collapses rows and results in commas for any disagreeing cells
check.tax %>%   
  mutate(Class=case_when(Class=="Actinopteri" ~ "Actinopterygii",.default=Class)) %>% #replacing old version of name
  group_by(seqid) %>% 
  summarise_all(list(~toString(unique(na.omit(.))))) %>% #unique for cells that are the same / duplicates
  mutate(across(c(Phylum,Class,Order,Family,Genus), #check the higher level taxa
                ~ case_when(str_detect(.,",") & !str_detect(Species,",") ~ word(.,sep=","), 
                            #if the cell has a comma in it
                            #indicating that the two cells that merged did not agree
                            #check if Species has a comma
                            #if not, just keep the first value 
                            #(because the difference is likely because of naming convention changes)
                            .default=.x)
                )
         ) %>% 
  mutate_all(., list(~na_if(.,""))) %>%
  #filter(.,grepl(",",Species)) %>% View() 
  #look at the dataframe with commas in the Species, what does that look like
  #was going to say anything with commas should get replaced with NA 
  #but that also affects cells where naming convention appears to be the only difference
  #e.g., Class name is not the same but Family and Genus are 
  View()


check.tax %>%   
  group_by(seqid) %>% 
  summarise_all(list(~toString(unique(na.omit(.))))) %>% #unique for cells that are the same / duplicates
  mutate(across(c(Phylum,Class,Order,Family,Genus), #check the higher level taxa
                ~ case_when(str_detect(.,",") & !str_detect(Species,",") ~ word(.,sep=","), 
                            .default=.x)
  )
  ) %>% 
  mutate_all(., list(~na_if(.,""))) %>%
  mutate(Species=case_when(str_detect(Species,",") ~ NA,
                                      .default=Species),
         Genus=case_when(str_detect(Genus,",") & is.na(Species) ~ NA, 
                         #if Genus has a comma and Species is NA 
                         #Genus is NA
                         str_detect(Genus,",") & !is.na(Species) ~ word(Genus,sep=","),
                         #if Genus has a comma and Species has a value
                         #keep the first value of Genus
                         .default=Genus),
         Family=case_when(str_detect(Family,",") & is.na(Genus) ~ NA, 
                         #if Family has a comma and Genus is NA 
                         #Family is NA
                         str_detect(Family,",") & !is.na(Genus) ~ word(Family,sep=","),
                         #if Family has a comma and Genus has a value
                         #keep the first value of Family
                         .default=Family),
         Order=case_when(str_detect(Order,",") & is.na(Family) ~ NA, 
                         #if Order has a comma and Family is NA 
                         #Order is NA
                         str_detect(Order,",") & !is.na(Family) ~ word(Order,sep=","),
                         #if Order has a comma and Family has a value
                         #keep the first value of Order
                         .default=Order),
         Class=case_when(str_detect(Class,",") & is.na(Order) ~ NA, 
                         #if Class has a comma and Order is NA 
                         #Class is NA
                         str_detect(Class,",") & !is.na(Order) ~ word(Class,sep=","),
                         #if Class has a comma and Order has a value
                         #keep the first value of Class
                         .default=Class),
         Phylum=case_when(str_detect(Phylum,",") & is.na(Class) ~ NA, 
                         #if Phylum has a comma and Class is NA 
                         #Phylum is NA
                         str_detect(Phylum,",") & !is.na(Class) ~ word(Phylum,sep=","),
                         #if Phylum has a comma and Class has a value
                         #keep the first value of Phylum
                         .default=Phylum)
         ) %>%
  View()


combined_tax %>%   
  group_by(seqid) %>% 
  summarise_all(list(~toString(unique(na.omit(.))))) %>% #unique for cells that are the same / duplicates
  mutate(across(c(Phylum,Class,Order,Family,Genus), #check the higher level taxa
                ~ case_when(str_detect(.,",") & !str_detect(Species,",") ~ word(.,sep=","), 
                            .default=.x)
  )
  ) %>% 
  mutate_all(., list(~na_if(.,""))) %>%
  mutate(Species=case_when(str_detect(Species,",") ~ NA,
                           .default=Species),
         Genus=case_when(str_detect(Genus,",") & is.na(Species) ~ NA, 
                         #if Genus has a comma and Species is NA 
                         #Genus is NA
                         str_detect(Genus,",") & !is.na(Species) ~ word(Genus,sep=","),
                         #if Genus has a comma and Species has a value
                         #keep the first value of Genus
                         .default=Genus),
         Family=case_when(str_detect(Family,",") & is.na(Genus) ~ NA, 
                          #if Family has a comma and Genus is NA 
                          #Family is NA
                          str_detect(Family,",") & !is.na(Genus) ~ word(Family,sep=","),
                          #if Family has a comma and Genus has a value
                          #keep the first value of Family
                          .default=Family),
         Order=case_when(str_detect(Order,",") & is.na(Family) ~ NA, 
                         #if Order has a comma and Family is NA 
                         #Order is NA
                         str_detect(Order,",") & !is.na(Family) ~ word(Order,sep=","),
                         #if Order has a comma and Family has a value
                         #keep the first value of Order
                         .default=Order),
         Class=case_when(str_detect(Class,",") & is.na(Order) ~ NA, 
                         #if Class has a comma and Order is NA 
                         #Class is NA
                         str_detect(Class,",") & !is.na(Order) ~ word(Class,sep=","),
                         #if Class has a comma and Order has a value
                         #keep the first value of Class
                         .default=Class),
         Phylum=case_when(str_detect(Phylum,",") & is.na(Class) ~ NA, 
                          #if Phylum has a comma and Class is NA 
                          #Phylum is NA
                          str_detect(Phylum,",") & !is.na(Class) ~ word(Phylum,sep=","),
                          #if Phylum has a comma and Class has a value
                          #keep the first value of Phylum
                          .default=Phylum)
  ) %>%
  View()

tax.combined.df <- combined_tax %>%   
  group_by(seqid) %>% 
  summarise_all(list(~toString(unique(na.omit(.))))) %>% #unique for cells that are the same / duplicates
  mutate(across(c(Phylum,Class,Order,Family,Genus), #check the higher level taxa
                ~ case_when(str_detect(.,",") & !str_detect(Species,",") ~ word(.,sep=","), 
                            .default=.x)
  )
  ) %>% 
  mutate_all(., list(~na_if(.,""))) %>%
  mutate(Species=case_when(str_detect(Species,",") ~ NA,
                           .default=Species),
         Genus=case_when(str_detect(Genus,",") & is.na(Species) ~ NA, 
                         #if Genus has a comma and Species is NA 
                         #Genus is NA
                         str_detect(Genus,",") & !is.na(Species) ~ word(Genus,sep=","),
                         #if Genus has a comma and Species has a value
                         #keep the first value of Genus
                         .default=Genus),
         Family=case_when(str_detect(Family,",") & is.na(Genus) ~ NA, 
                          #if Family has a comma and Genus is NA 
                          #Family is NA
                          str_detect(Family,",") & !is.na(Genus) ~ word(Family,sep=","),
                          #if Family has a comma and Genus has a value
                          #keep the first value of Family
                          .default=Family),
         Order=case_when(str_detect(Order,",") & is.na(Family) ~ NA, 
                         #if Order has a comma and Family is NA 
                         #Order is NA
                         str_detect(Order,",") & !is.na(Family) ~ word(Order,sep=","),
                         #if Order has a comma and Family has a value
                         #keep the first value of Order
                         .default=Order),
         Class=case_when(str_detect(Class,",") & is.na(Order) ~ NA, 
                         #if Class has a comma and Order is NA 
                         #Class is NA
                         str_detect(Class,",") & !is.na(Order) ~ word(Class,sep=","),
                         #if Class has a comma and Order has a value
                         #keep the first value of Class
                         .default=Class),
         Phylum=case_when(str_detect(Phylum,",") & is.na(Class) ~ NA, 
                          #if Phylum has a comma and Class is NA 
                          #Phylum is NA
                          str_detect(Phylum,",") & !is.na(Class) ~ word(Phylum,sep=","),
                          #if Phylum has a comma and Class has a value
                          #keep the first value of Phylum
                          .default=Phylum)
  ) 

#how many unique Chordate species  
tax.combined.df %>% filter(Phylum=="Chordata") %>% .$Species %>% unique() %>% length()
#are any seqids repeated? 
summary(duplicated(tax.combined.df$seqid))

# Now we have some weird values from the nt tax database that didn't fill out correctly
# they leave an NA and fill the rest of the taxonomy 
# how do we fill that in if we have the information  

# TBD  


# --------------------------------------------------------------------------------- #
# conform tax table  
rownames(tax.combined.df) <- tax.combined.df$seqid

# upload asv table  
curated_lulu <- readRDS(paste0(dir_results,'/lulu-clustertable.rds'))
curated_asv <- curated_lulu$curated_table
cols_asv <- word(colnames(curated_asv),sep = "_",end=3) 
colnames(curated_asv) <- cols_asv

# upload metadata
samplelist <- read.csv("/Users/lauragivens/Desktop/R/BZrookery_eDNA/Metadata_BZrookery.csv")
rownames(samplelist) <- samplelist$SampleName_Long
samplelist$SampleName_Long <- NULL

# assemble ps object
otu <- otu_table((curated_asv),taxa_are_rows = TRUE)
meta <- sample_data(samplelist)
tax <- tax_table(as.matrix(tax.combined.df))

ps <- phyloseq(otu,meta,tax)
ps

##########################################################################################################
# Save 
saveRDS(ps,paste0(dir_results,"/ps.combined.rds"))

write.csv(tax.combined.df,paste0(dir_results,"/tax.combined.df.csv"))
saveRDS(tax.combined.df,paste0(dir_results,"/tax.combined.df.rds"))
save.image(paste0(dir_data,'/CombineBOLDnt.RData'))





  
