
# setup -------------------------------------------------------------------

library(tidyverse)
library(rio)
# library(GEOquery) # for using gunzip()
library(ampir)
library(UpSetR)
source("../Github/r-mass-spectrometry-tools/bin/digest_aa_sequence.R")
source("../Github/r-mass-spectrometry-tools/bin/Calculate_monoisotopic_mass.R")


# importing files ---------------------------------------------------------

 
# gz_files <- dir("resources/", full.names = TRUE)[grep(".faa.gz$", dir("resources/"))]
# 
# for(i in 1:length(gz_files)){
#   gunzip(gz_files[i])
# }


data <- lapply(dir("resources/", full.names = TRUE), function(x){
  
  # Reads in .faa files as dataframe
  read_faa(x) %>% 
    mutate(organism = x)
  
})

data <- bind_rows(data)

data$organism <- gsub("resources/", "", data$organism)

data$organism <- unlist(lapply(data$organism, function(x){
  unlist(strsplit(x, split = "_"))[1]
}))

unique(data$organism)

# global functions --------------------------------------------------------


# Digesting protein
digest_tryp <- function(sequence){
  
  tryCatch(
    digest_aa_sequence(sequence, enzyme = "trypsin", missed = 0)$peptide,
    error = function(e) NA,
    warning = function(w) NA
  )
}


# digesting proteomes -----------------------------------------------------


min_pep_len <- 7
max_pep_lin <- 40
proton <- 1.007276466

# Digesting all proteomes
# filtering for a minimum peptide length
system.time(
  proteome <- data %>% 
    mutate(peptide = map(seq_aa, digest_tryp)) %>% # digesting proteome TAKES A LONG TIME!
    unnest(peptide) %>% 
    select(-seq_aa) %>% 
    mutate(pep_len = nchar(peptide)) %>% # Counting peptide length
    filter(pep_len >= min_pep_len,
           pep_len <= max_pep_lin) # Filtering for peptide size
) # 3.25 min to finish on pc

195.38/60

# Calculating monoisotopic mass
# All lysines have acetyl modification
proteome <- proteome %>% 
  mutate(monoisotopic_mass = calculate_monoisotopic_mass(peptide, IAA = TRUE, AcK = FALSE))

# Calculating most likely charge state based on K|R|H + 1 (for n-terminus)
proteome$z <- sapply(proteome$peptide, function(x){
  length(unlist(gregexpr("K|R|H", x))) + 1 # For N-terminus
})

# calculating precursor m/z
# assuming z=2; z=3
proteome <- proteome %>% 
  filter(z %in% c(2,3,4)) %>% 
  mutate(precursor_mz = (monoisotopic_mass + z * proton) / z)

ggplot(proteome) +
  geom_density(aes(x = precursor_mz, color = organism)) +
  scale_x_log10()

ggplot(proteome) +
  geom_bar(aes(x = z, fill = organism), position = "stack")





# upset plot --------------------------------------------------------------


data_wide <- proteome %>% 
  select(organism, peptide) %>% 
  distinct() %>% 
  mutate(present = 1) %>% 
  spread(organism, present)


upset_data <- data_wide %>% 
  column_to_rownames("peptide")

upset_data[!is.na(upset_data)] <- 1
upset_data[is.na(upset_data)] <- 0

upset_data <- upset_data %>% 
  rownames_to_column("peptide")

names(upset_data)

png(filename = "peptide_comparison_across_organisms_upset_plot.png", width = 1024, height = 768)

upset(upset_data, 
      nsets = 25, 
      nintersects = 23,
      sets = c("Human", "Cow", "Trout", "AppleLeaves", "RiceFlour", "SpinachLeaves", "TomatoLeaves"),
      keep.order = TRUE, 
      set.metadata = NULL, 
      intersections = NULL,
      matrix.color = "gray23", 
      main.bar.color = "gray23",
      mainbar.y.label = "Intersection Size", 
      mainbar.y.max = NULL,
      sets.bar.color = "gray23", 
      sets.x.label = "Set Size",
      point.size = 2, 
      line.size = 0.7, 
      mb.ratio = c(0.7, 0.3),
      expression = NULL, 
      att.pos = NULL, 
      # att.color = main.bar.color,
      order.by = c("freq"), 
      decreasing = T,
      show.numbers = "yes", 
      number.angles = 0, 
      group.by = "degree",
      cutoff = NULL, 
      queries = NULL, 
      query.legend = "none",
      shade.color = "gray88", 
      shade.alpha = 0.25, 
      matrix.dot.alpha = 0.5,
      empty.intersections = NULL, 
      color.pal = 1, 
      boxplot.summary = NULL,
      attribute.plots = NULL, 
      scale.intersections = "identity",
      scale.sets = "identity", 
      text.scale = 1.75, 
      set_size.angles = 0,
      set_size.show = FALSE, 
      set_size.numbers_size = NULL,
      set_size.scale_max = NULL)

dev.off()

# data export -------------------------------------------------------------


export(proteome, file = "data/refseq_digested_proteomes.csv")
