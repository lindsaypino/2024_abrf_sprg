
# setup -------------------------------------------------------------------

library(tidyverse)
library(rio)
library(janitor)
library(patchwork)
library(ggrepel)
source("R/00_functions.R")
library(venneuler)



# data import -------------------------------------------------------------

# combined_fasta <- import_fasta_as_df(fname = "../fasta/Combined_proteomes.fasta")

pep <- import("dia_44min_ratio/QuantReport.elib.peptides.txt") %>% clean_names()

peptide_group <- import("Combined_proteome_in_silico_digestion.csv")



ground_truth <- data.frame(
  organism_short = c("TROUT", "HUMAN", "BOVIN"),
  A = c(50, 45, 5),
  B = c(50, 20, 30),
  C = c(50, 3, 47)
)

ground_truth <- ground_truth %>% 
  mutate(A_v_B = log2(A) - log2(B),
         A_v_C = log2(A) - log2(C),
         B_v_C = log2(B) - log2(C))

ground_truth <- ground_truth %>% 
  gather(comparison, log2fc, contains("_v_")) %>% 
  mutate(comparison = str_replace_all(comparison, "_v_", " v "))



# Formatting peptide dataframe --------------------------------------------

# Cleaning up peptide column
# removing modification
pep <- pep %>% 
  mutate(peptide = str_remove_all(peptide, "\\[\\+57.021464\\]"))

# Cleaning up variable names
# Merging with fasta dataframe
pep_df <- pep %>% 
  rename(A = x220602_abrf_ratio_a_mz_ml,
         B = x220602_abrf_ratio_b_mz_ml,
         C = x220602_abrf_ratio_c_mz_ml) %>% 
  select(-protein) %>% 
  right_join(peptide_group, .)

# transposing data to long format
# transforming to Log2
pep_condition_long <- pep_df %>% 
  gather(condition, abundance, A, B, C) %>% 
  mutate(abundance = ifelse(abundance == 0, NA, log2(abundance)))

# Visualizing the peptide abundances across the 3 conditions
# Consider equalizing the medians, which would impact the final ratios. 
# the medians are pretty close for now.
ggplot(pep_condition_long) +
  geom_boxplot(aes(x = condition, y = abundance, color = condition)) +
  geom_hline(yintercept = median(pep_condition_long$abundance, na.rm = TRUE))

# Calculating the ratios
pep_log2fc <- pep_condition_long %>% 
  spread(condition, abundance) %>% 
  mutate(`A v B` = A-B,
         `A v C` = A-C,
         `B v C` = B-C) %>% 
  gather(comparison, log2fc, contains(" v "))

# removing trypsin peptides because they probably came from the Worthington Trypsin
pep_log2fc <- pep_log2fc %>% 
  filter(ProteinDescription != "Cationic Trypsin")



# Protein summarization ---------------------------------------------------


###################
#### Functions ####
###################

# Summarizing Protein abundance
protein_summary_medpolish <- function(data){
  
  # selecting the data
  tmp <- data %>%
    select(peptide, condition, abundance) %>%
    # unite(id, group, sample_id, run_order, sep = "_") %>%
    # unite(feature, PeptideSeq, PeptideModSeq, PrecursorCharge, sep = "_") %>%
    mutate(abundance = 2^abundance) %>%
    spread(peptide, abundance) %>%
    column_to_rownames("condition")
  
  # Tukey Median polish
  meddata  <-  medpolish(tmp,
                         na.rm=TRUE,
                         trace.iter = FALSE)
  
  # Extracting the results
  tmpresult <- data.frame(abundance = (meddata$overall + meddata$row)) %>%
    rownames_to_column("id") %>%
    mutate(abundance = log2(abundance))
  
}

# Counting the number of peptides per protein
count_peptide_per_group <- function(x){
  length(unique(x$peptide))
}

###################
## End Functions ##
###################


proteinquants <- pep_condition_long %>%
  group_by(ProteinAccession, ProteinDescription, GeneName, Organism, organism_short, organism_n) %>%
  nest() %>%
  mutate(peptide_n = map_dbl(data, count_peptide_per_group)) %>%
  mutate(abundance = map(data, protein_summary_medpolish)) %>%
  select(-data) %>%
  unnest(abundance)

prot_log2fc <- proteinquants %>% 
  spread(id, abundance) %>% 
  mutate(`A v B` = A-B,
         `A v C` = A-C,
         `B v C` = B-C) %>% 
  gather(comparison, log2fc, contains(" v "))


prot_wide <- proteinquants %>%
  spread(id, abundance)


# plots -------------------------------------------------------------------


## Navarro plot using unique peptides
group_median <- pep_log2fc %>% 
  filter(organism_n == 1) %>% 
  group_by(organism_short, comparison) %>% 
  summarize(log2fc = median(log2fc, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(label = round(2^log2fc, digits = 3))

# left side plot
nav_1 <- ggplot(pep_log2fc %>% filter(organism_n == 1, B > 5),
                aes(x = B, y = log2fc, color = organism_short)) +
  geom_point(size = 1, shape = 21) +
  geom_hline(data = group_median, 
             aes(yintercept = log2fc, color = organism_short)) +
  geom_hline(data = ground_truth,
             aes(yintercept = log2fc, color = organism_short), linetype = "dashed") +
  geom_text(data = group_median, 
                  aes(x = 30, label = label), vjust = 1) +
  facet_grid(vars(comparison)) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 22) +
  theme(strip.text = element_blank()) +
  guides(color = "none") +
  labs(y = expression(Log[2]~fold~change),
       x = expression(Log[2]~B))

# Right side plot
nav_2 <- ggplot(pep_log2fc %>% filter(organism_n == 1, B > 5),
                aes(x = organism_short, y = log2fc, color = organism_short)) +
  geom_boxplot(size = 0.75) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 22) +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(vars(comparison)) +
  labs(y = NULL, 
       x = NULL,
       color = "Organism")

# Stitching plots together
nav_1 + nav_2 + 
  plot_layout(guides = "collect", widths = c(2,1))

# Saving
ggsave(filename = "figures/Peptide_lfqbench_figure.png", width = 9.5, height = 5.5, units = "in")
ggsave(filename = "figures/Peptide_lfqbench_figure.svg", width = 9.5, height = 5.5, units = "in")


pep_df %>% 
  filter(!is.na(organism_short)) %>% 
  ggplot() +
  geom_bar(aes(x = organism_short)) +
  coord_flip()

pep_venn <- pep_df %>% 
  select(peptide, organism_short) %>% 
  filter(!is.na(organism_short)) %>% 
  mutate(val = 1) %>% 
  spread(organism_short, val)



plot(venneuler(c(HUMAN=0.208871157, BOVIN=0.195728702, TROUT=0.293006337, 
                 "HUMAN&BOVIN"=0.178537902, "HUMAN&TROUT"=0.007157944, 
                 "BOVIN&TROUT"=0.004693734 ,"HUMAN&BOVIN&TROUT"=0.112004224)))


venn_weights <- table(pep_venn$organism_short)/nrow(pep_venn)

names(venn_weights)[grep(";", names(venn_weights))] <- gsub(";", "&", names(venn_weights)[grep(";", names(venn_weights))])





# Plot 2 ------------------------------------------------------------------


group_median_2 <- pep_log2fc %>% 
  group_by(organism_short, comparison) %>% 
  summarize(log2fc = median(log2fc, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(label = round(2^log2fc, digits = 3))


nav_3 <- ggplot(pep_log2fc %>% filter(!is.na(Organism), B > 5),
                aes(x = B, y = log2fc, color = organism_short)) +
  geom_point(size = 0.75) +
  geom_hline(data = group_median_2, 
             aes(yintercept = log2fc, color = organism_short)) +
  facet_grid(vars(comparison)) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 18) +
  theme(strip.text = element_blank()) +
  guides(color = "none") +
  labs(y = expression(Log[2]~fold~change),
       x = expression(Log[2]~B))

nav_4 <- ggplot(pep_log2fc %>% filter(!is.na(Organism), B > 5),
                aes(x = organism_short, y = log2fc, color = organism_short)) +
  geom_boxplot(size = 0.75) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 18) +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(vars(comparison)) +
  labs(y = NULL, 
       x = NULL,
       color = "Organism")

nav_3 + nav_4 + 
  plot_layout(guides = "collect", widths = c(2,1))

ggsave(filename = "raw/lfqbench_figure_2.png", width = 12, height = 8)


# Counting peptides per organism group
ggplot(pep_df %>% 
         filter(!is.na(organism_short))) +
  geom_bar(aes(x = organism_short)) +
  coord_flip() +
  labs(x = NULL)

ggsave(filename = "figures/Peptides_per_organism.png")



# Protein plots -----------------------------------------------------------

## Navarro plot using unique peptides
prot_group_median <- prot_log2fc %>% 
  filter(organism_n == 1) %>% 
  group_by(organism_short, comparison) %>% 
  summarize(log2fc = median(log2fc, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(label = round(2^log2fc, digits = 3))

# left side plot
prot_1 <- ggplot(prot_log2fc %>% filter(organism_n == 1, B > 5),
                aes(x = B, y = log2fc, color = organism_short)) +
  geom_point(size = 0.75) +
  geom_hline(data = prot_group_median, 
             aes(yintercept = log2fc, color = organism_short)) +
  # geom_smooth(method = "loess") +
  geom_text_repel(data = prot_group_median, 
                  aes(x = 30, label = label), position = position_nudge(0.01), hjust = 1) +
  facet_grid(vars(comparison)) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 18) +
  theme(strip.text = element_blank()) +
  guides(color = "none") +
  labs(y = expression(Log[2]~fold~change),
       x = expression(Log[2]~B))

# Right side plot
prot_2 <- ggplot(prot_log2fc %>% filter(organism_n == 1, B > 5),
                aes(x = organism_short, y = log2fc, color = organism_short)) +
  geom_boxplot(size = 0.75) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 18) +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(vars(comparison)) +
  labs(y = NULL, 
       x = NULL,
       color = "Organism")

# Stitching plots together
prot_1 + prot_2 + 
  plot_layout(guides = "collect", widths = c(2,1)) +
  plot_annotation(title = "Protein Group LFQbench")

# Saving
ggsave(filename = "figures/Protein_lfqbench_figure.png", width = 12, height = 8)


# # Figure 2
# 
# prot_group_median_2 <- pep_log2fc %>% 
#   group_by(organism_short, comparison) %>% 
#   summarize(log2fc = median(log2fc, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   mutate(label = round(2^log2fc, digits = 3))
# 
# 
# prot_3 <- ggplot(prot_log2fc %>% filter(!is.na(Organism), B > 5),
#                 aes(x = B, y = log2fc, color = organism_short)) +
#   geom_point(size = 0.75) +
#   geom_hline(data = prot_group_median_2, 
#              aes(yintercept = log2fc, color = organism_short)) +
#   facet_grid(vars(comparison)) +
#   scale_color_brewer(palette = "Dark2") +
#   theme_classic(base_size = 18) +
#   theme(strip.text = element_blank()) +
#   guides(color = "none") +
#   labs(y = expression(Log[2]~fold~change),
#        x = expression(Log[2]~B))
# 
# prot_4 <- ggplot(prot_log2fc %>% filter(!is.na(Organism), B > 5),
#                 aes(x = organism_short, y = log2fc, color = organism_short)) +
#   geom_boxplot(size = 0.75) +
#   scale_color_brewer(palette = "Dark2") +
#   theme_classic(base_size = 18) +
#   theme(axis.text = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks = element_blank()) +
#   facet_grid(vars(comparison)) +
#   labs(y = NULL, 
#        x = NULL,
#        color = "Organism")
# 
# prot_3 + prot_4 + 
#   plot_layout(guides = "collect", widths = c(2,1)) +
#   plot_annotation(title = "Protein Group LFQbench")
# 
# # Saving
# ggsave(filename = "raw/lfqbench_2_protein_figure.png", width = 12, height = 8)
# 
# 
# # Counting peptides per organism group
# ggplot(prot_wide %>% 
#          filter(!is.na(organism_short))) +
#   geom_bar(aes(x = organism_short)) +
#   coord_flip() +
#   labs(x = "Organism group",
#        title = "Protein groups")
# ggsave(filename = "raw/Protein_group_per_organism_group_count.png")


# data export -------------------------------------------------------------

export(pep_df, file = "final_data/peptide_data.csv")
export(prot_wide, file = "final_data/protein_data.csv")
