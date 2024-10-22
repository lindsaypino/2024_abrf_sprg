
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


# data export -------------------------------------------------------------


export(peptide_group, file = "Combined_proteome_in_silico_digestion.csv")

