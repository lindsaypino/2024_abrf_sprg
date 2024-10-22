
# setup -------------------------------------------------------------------


library(tidyverse)
library(rio)
library(janitor)
source("R/00_functions.R")


# data import -------------------------------------------------------------

combined_fasta <- import_fasta_as_df(fname = "../fasta/Combined_proteomes.fasta")

pep <- import("dia_44min_ratio/QuantReport.elib.peptides.txt") %>% clean_names()


# global functions --------------------------------------------------------


digest_tryp <- function(sequence){
  
  tryCatch(
    digest_aa_sequence(sequence, enzyme = "trypsin", missed = 1)$peptide,
    error = function(e) NA,
    warning = function(w) NA
  )
}




# In silico digestion of FASTA --------------------------------------------


# Digestion settings
min_pep_len <- 6
proton <- 1.007276466


# cleaning up organism column in FASTA dataframe
# Trout proteins are labeled with weird annotation.
combined_fasta <- combined_fasta %>% 
  mutate(Organism = ifelse(!combined_fasta$Organism %in% c("HUMAN", "BOVIN"), "TROUT", Organism))


# Digesting the combined proteomes
# filtering for a minimum peptide length
combined_prot <- combined_fasta %>% 
  mutate(peptide = map(ProteinSequence, digest_tryp)) %>% # digesting proteomes...THIS TAKES A LONG TIME!
  unnest(peptide) %>% 
  select(-ProteinSequence) %>% 
  mutate(pep_len = nchar(peptide)) %>% # Counting peptide length
  filter(pep_len >= min_pep_len) # Filtering based on peptide length


# Grouping peptides based on which ORGANISM they were identified in
# THIS ALSO TAKES A LONG TIME!!!!
peptide_group <- combined_prot %>% 
  group_by(peptide, pep_len) %>% 
  summarize(ProteinAccession = paste(ProteinAccession, collapse = ";"),
            ProteinDescription = paste(ProteinDescription, collapse = ";"),
            GeneName = paste(GeneName, collapse = ";"),
            Organism = paste(Organism, collapse = ";")) %>% 
  ungroup()

# Summarizing organism info
# These are needed for filtering and plotting
peptide_group$organism_short <- unlist(lapply(peptide_group$Organism, function(x){
  paste(unique(unlist(strsplit(x, split = ";"))), collapse = ";")
}))

# Counting how many orgranism each peptide is derived
peptide_group <- peptide_group %>% 
  mutate(organism_n = str_count(organism_short, pattern = ";") + 1)

# plotting how peptides break down
ggplot(peptide_group) +
  geom_bar(aes(x = organism_n))

# plotting how peptides break down
ggplot(peptide_group) +
  geom_bar(aes(x = organism_short))



# data export -------------------------------------------------------------


export(peptide_group, file = "Combined_proteome_in_silico_digestion.csv")



# Venn Diagram ------------------------------------------------------------

library(venneuler)
library(eulerr)

plot(euler(c(HUMAN=0.202239262, 
                 BOVIN=0.203133213, 
                 TROUT=0.474668545, 
                 "HUMAN&BOVIN"=0.092548043, 
                 "HUMAN&TROUT"=0.002512678, 
                 "BOVIN&TROUT"=0.001895723 ,
                 "HUMAN&BOVIN&TROUT"=0.023002536)))



venn_weights <- table(peptide_group$organism_short)

names(venn_weights)[grep(";", names(venn_weights))] <- gsub(";", "&", names(venn_weights)[grep(";", names(venn_weights))])

venn_weights

plot(euler(c(HUMAN=1005370, 
             BOVIN=1009814, 
             TROUT=2359668, 
             "HUMAN&BOVIN"=460074, 
             "HUMAN&TROUT"=12491, 
             "BOVIN&TROUT"=9424 ,
             "HUMAN&BOVIN&TROUT"=114350)))


venn_weights
venneuler(venn_weights)

venn_data <- peptide_group %>% 
  select(peptide, organism_short) %>% 
  mutate(var = 1) %>% 
  mutate(organism_short = strsplit(organism_short, split = ";")) %>% 
  unnest(organism_short)

venn_data <- venn_data %>% 
  spread(organism_short, var)

venn_data[is.na(venn_data)] <- 0

venn_data <- venn_data %>% 
  column_to_rownames("peptide")

plot(venneuler(venn_data))
