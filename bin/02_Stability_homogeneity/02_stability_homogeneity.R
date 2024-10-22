
# setup -------------------------------------------------------------------

# library(RCurl)
library(tidyverse)
library(janitor)
library(rio)
library(pheatmap)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(patchwork)


# data import -------------------------------------------------------------

# dataURL <- getURL("ftp://massive.ucsd.edu/2022%20Multi-Species%20Standard%20Study/42/EncyclopeDIA/sPRG_stability_homogeneity_QuantReport.elib.peptides.txt", userpwd="sprg:Pr0teomics")
# temp <- read_table(dataURL)

sprg_peptide <- import("42/raw/sPRG_stability_homogeneity_QuantReport.elib.peptides.txt")
in_silico_pep <- import("Combined_proteome_in_silico_digestion.csv")


# data cleanup ------------------------------------------------------------

# Cleaning up the Quant table
sprg <- sprg_peptide %>% 
  gather(FileName, Intensity, contains(".mzML")) %>% 
  mutate(FileName = str_remove_all(FileName, ".mzML"),
         abundance = ifelse(Intensity < 1, 0, log2(Intensity)),
         Peptide = str_remove_all(Peptide, "\\[\\+57.021464\\]"))


# Making an annotation file
annot <- sprg %>% 
  select(FileName) %>% 
  distinct() %>% 
  mutate(ratio = str_sub(FileName, start = -1L, end = -1L),
         experiment = case_when(
           grepl("230109_RT_", FileName) ~ "Homogeneity",
           TRUE ~ "Stability"),
         condition = case_when(
           grepl("RT_01", FileName) ~ "01",
           grepl("RT_06", FileName) ~ "06",
           grepl("RT_11", FileName) ~ "11",
           grepl("neg20", FileName) ~ "neg20",
           grepl("exploris480", FileName) ~ "initial")
  )


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
  gather(comparison, log2fc, contains("_v_"))


# Normalization -----------------------------------------------------------


group_median <- tapply(sprg$abundance, sprg$FileName, median, na.rm = TRUE)
global_median <- median(sprg$abundance, na.rm = TRUE)


sprg <- sprg %>% 
  mutate(abundance_norm = ifelse(abundance < 1, 0, abundance - group_median[.$FileName] + global_median))


# spreading to wide format
sprg_wide <- sprg %>% 
  filter(abundance > 1) %>% 
  select(-Intensity, -abundance, -numFragments) %>% 
  spread(FileName, abundance_norm)




# Filtering peptides
# Merging with annotation
# Merging with In silico digest of combined proteomes
sprg_ratio <- sprg %>% 
  filter(abundance > 0) %>% 
  full_join(annot, .) %>% 
  select(-Protein, -FileName, -Intensity, -abundance) %>% 
  rename(peptide = Peptide) %>% 
  right_join(in_silico_pep, .)

# spreading samples
sprg_ratio <- sprg_ratio %>% 
  spread(ratio, abundance_norm) 

# calculating log2fc for each channel combo
sprg_log2fc <- sprg_ratio %>% 
  mutate(A_v_B = A - B,
         A_v_C = A - C,
         B_v_C = B - C) %>% 
  gather(comparison, log2fc, contains("_v_"))

# 
# # Not sure how to make this work
# sprg_log2fc_wide <- sprg_log2fc %>%
#   filter(!is.na(pep_len),
#          organism_n == 1) %>%
#   select(-A, -B, -C, -Organism) %>%
#   select(peptide, pep_len, organism_short, experiment, condition, comparison, log2fc) %>%
#   unite(id, experiment, condition, comparison, sep = "_") %>%
#   spread(id, log2fc)
# 
# names(sprg_log2fc_wide)
# 
# sprg_log2fc_wide %>%
#   ggplot(aes(x = Homogeneity_01_A_v_C, y = Homogeneity_06_A_v_C, color = organism_short)) +
#   geom_point(alpha = 0.1) +
#   scale_color_brewer(palette = "Dark2") 
# 
# sprg_log2fc_wide %>%
#   ggplot(aes(x = Homogeneity_01_A_v_C, y = Homogeneity_11_A_v_C, color = organism_short)) +
#   geom_point(alpha = 0.1) +
#   scale_color_brewer(palette = "Dark2") 
# 
# sprg_log2fc_wide %>%
#   ggplot(aes(x = Homogeneity_06_A_v_C, y = Homogeneity_11_A_v_C, color = organism_short)) +
#   geom_point(alpha = 0.1) +
#   scale_color_brewer(palette = "Dark2") 
# 
#  
 
# GGally::ggscatmat(sprg_peptide, columns = 4:8, 
#                   # color = "organism_short", 
#                   alpha = 0.01) +
#   scale_color_brewer(palette = "Dark2")
# 
# 
# ab_wide <- sprg_log2fc %>%
#   filter(!is.na(pep_len),
#          organism_n == 1,
#          comparison == "A_v_B") %>%
#   select(-A, -B, -C, -Organism) %>%
#   select(peptide, pep_len, organism_short, experiment, condition, comparison, log2fc) %>%
#   unite(id, experiment, condition, comparison, sep = "_") %>%
#   spread(id, log2fc)
# 
# GGally::ggscatmat(ab_wide, columns = 4:8, 
#                   # color = "organism_short", 
#                   alpha = 0.01) +
#   scale_color_brewer(palette = "Dark2")
# 
# 
# ac_wide <- sprg_log2fc %>%
#   filter(!is.na(pep_len),
#          organism_n == 1,
#          comparison == "A_v_C") %>%
#   select(-A, -B, -C, -Organism) %>%
#   select(peptide, pep_len, organism_short, experiment, condition, comparison, log2fc) %>%
#   unite(id, experiment, condition, comparison, sep = "_") %>%
#   spread(id, log2fc)
# 
# GGally::ggscatmat(ac_wide, columns = 4:8, 
#                   color = "organism_short",
#                   alpha = 0.01) +
#   scale_color_brewer(palette = "Dark2")
# 
# 
# bc_wide <- sprg_log2fc %>%
#   filter(!is.na(pep_len),
#          organism_n == 1,
#          comparison == "B_v_C") %>%
#   select(-A, -B, -C, -Organism) %>%
#   select(peptide, pep_len, organism_short, experiment, condition, comparison, log2fc) %>%
#   unite(id, experiment, condition, comparison, sep = "_") %>%
#   spread(id, log2fc)
# 
# GGally::ggscatmat(bc_wide, columns = 4:8, 
#                   color = "organism_short",
#                   alpha = 0.01) +
#   scale_color_brewer(palette = "Dark2")
# 


# plots -------------------------------------------------------------------

peptide_count <- length(unique(sprg$Peptide))
length(unique(sprg$FileName))
paste0(peptide_count, " peptides quantified.")

# Histogram
sprg %>% 
  filter(Intensity > 1) %>% 
  ggplot() +
  geom_bar(aes(x = FileName)) +
  geom_hline(yintercept = peptide_count) +
  annotate("text", y = peptide_count * 1.05, x = 15/2, label = paste0(peptide_count, " peptides quantified"), size = 5) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Peptides (1% FDR)",
       x = "FileName")

ggsave(filename = "Stability_homogeneity_barplot.png")

sprg_ratio %>% 
  filter(!is.na(organism_short)) %>% 
  # filter(organism_n == 1) %>% 
  unite(id, experiment, condition, sep = "_", remove = FALSE) %>% 
  ggplot() +
  geom_bar(aes(x = id, fill= organism_short)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL,
       fill = "Organism(s)")
ggsave(filename = "peptide_barplot_organism.png")

sprg_ratio %>% 
  filter(organism_n == 1) %>% 
  ggplot(aes(x = B, y = A-C, color = organism_short)) +
  geom_point(alpha = 0.1) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(limits = c(10, 35)) +
  facet_wrap(vars(condition)) +
  theme_minimal(base_size = 18) +
  coord_fixed() +
  labs(y = "Log2 Fold Change (A - C)",
       x = "B (Log2)",
       color = "Organism")

ggsave(filename = "LFQBench_figure_stability_homogeneity.png", width = 10, height = 6)


sprg_ratio %>% 
  filter(organism_n == 1) %>%
  ggplot(aes(x = A-C, fill = organism_short)) +
  geom_histogram(bins = 50, color = "black", linewidth = 0.01) +
  scale_fill_brewer(palette = "Dark2") +
  # scale_x_continuous(limits = c(10, 35)) +
  facet_wrap(vars(condition)) +
  theme_minimal(base_size = 18) +
  labs(x = "Log2 Fold Change (A - C)",
       color = "Organism")

# Density plot
sprg_ratio %>% 
  filter(organism_n == 1) %>%
  ggplot(aes(x = A-C, color = organism_short)) +
  geom_density() +
  scale_fill_brewer(palette = "Dark2") +
  # scale_x_continuous(limits = c(10, 35)) +
  facet_wrap(vars(condition), ncol = 1) +
  theme_minimal(base_size = 18) +
  labs(x = "Log2 Fold Change (A - C)",
       color = "Organism")


sprg_log2fc %>% 
  filter(organism_n == 1,
         comparison == "A_v_C") %>%
  ggplot(aes(x = log2fc, fill = organism_short)) +
  geom_histogram(bins = 50, color = "black") +
  scale_fill_brewer(palette = "Dark2") +
  # scale_x_continuous(limits = c(10, 35)) +
  facet_grid(vars(condition), vars(comparison)) +
  theme_minimal(base_size = 18) +
  labs(x = "Log2 Fold Change",
       color = "Organism")


# LFQBench for all samples
sprg_log2fc %>% 
  filter(organism_n == 1,
         comparison == "A_v_C") %>%
  ggplot() + 
  geom_point(aes(x = B, y = log2fc, color = organism_short), alpha = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(cols = vars(condition)) +
  theme_bw(base_size = 14) +
  # coord_fixed() +
  labs(y = "Log2 fold change (A v C)",
       x = "Log2 B")


# sprg_ratio %>% 
#   filter(organism_n == 1,
#          !is.na(B)) %>% 
#   ggplot(aes(x = cut_width(B, 2), y = A-C, color = organism_short)) +
#   geom_boxplot(alpha = 0.1, outlier.alpha = 0) +
#   scale_color_brewer(palette = "Dark2") +
#   # scale_x_continuous(limits = c(10, 35)) +
#   facet_wrap(vars(condition)) +
#   theme_minimal(base_size = 18) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(y = "Log2 Fold Change (A - C)",
#        x = "B (Log2)",
#        color = "Organism")
# 



sprg_log2fc %>% 
  filter(organism_n == 1) %>% 
  ggplot(aes(x = log2fc, color = condition)) +
  geom_density() +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(rows = vars(organism_short),
             cols = vars(comparison)) +
  theme_classic(base_size = 18)
ggsave(filename = "Density_plot_log2FC.png")

# Log2 median normalization
ggplot(sprg) +
  geom_boxplot(aes(x = FileName, y = abundance_norm)) +
  geom_hline(yintercept = median(sprg$abundance, na.rm = TRUE)) +
  coord_flip() +
  theme_classic(base_size = 18)
ggsave(filename = "abundance_boxplot_original.png")



# Correlation -------------------------------------------------------------




unique(sprg$FileName)
# spreading to wide format
data_cor <- sprg_ratio %>% 
  filter(abundance > 1) %>% 
  filter(grepl("230109_", FileName)) %>%
  mutate(FileName = str_remove_all(FileName, "230109_")) %>% 
  select(-Intensity, -abundance, -numFragments) %>% 
  spread(FileName, abundance_norm)


mat_abundance <- cor(select_if(data_cor, is.numeric), 
                     use = "pairwise.complete.obs", 
                     method = "spearman")


# Data frame with column annotations.
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(mat)

# List with colors for each annotation.
mat_colors <- list(group = brewer.pal(3, "Set1"))
names(mat_colors$group) <- unique(col_groups)


pheatmap(mat_abundance,
         cellwidth = 15, 
         cellheight = 15, 
         scale = "none",
         color = viridis(45),
         # color = brewer.pal(50, "Spectral"),
         fontsize = 8)


data_cor %>% 
  ggplot() +
  geom_point(aes(x = RT_11A, y = RT_06A, color = ), alpha = 0.2) +
  coord_fixed()

test <- data_cor %>% 
  rename(peptide = Peptide) %>% 
  right_join(., in_silico_pep)

test %>% 
  # filter(organism_n == 1) %>% 
  ggplot() +
  geom_point(aes(x = RT_11A, y = RT_06A, color = organism_short), alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_brewer(palette = "Dark2") +
  coord_fixed()



# LFQ Bench figure --------------------------------------------------------


# plots -------------------------------------------------------------------


## Navarro plot using unique peptides
group_median <- sprg_log2fc %>% 
  filter(organism_n == 1) %>% 
  group_by(organism_short, comparison, experiment, condition) %>% 
  summarize(log2fc = median(log2fc, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(label = round(2^log2fc, digits = 3))




# left side plot
plot_1 <- sprg_log2fc %>% 
  filter(organism_n == 1, 
         condition == "01",
         comparison == "A_v_C",
         B > 5
         ) %>% 
  ggplot(aes(x = B, y = log2fc, color = organism_short)) +
  geom_point(size = 0.75, shape = 21) +
  geom_hline(data = group_median %>% filter(condition == "01", comparison == "A_v_C"), 
             aes(yintercept = log2fc, color = organism_short)) +
  geom_hline(data = ground_truth %>% filter(comparison == "A_v_C"),
             aes(yintercept = log2fc, color = organism_short), linetype = "dashed") +
  geom_text(data = group_median %>% filter(condition == "01", comparison == "A_v_C"), 
                  aes(x = 32, label = label), vjust = 0) +
  facet_grid(vars(comparison)) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(-10, 10)) +
  theme_minimal(base_size = 18) +
  theme(strip.text = element_blank()) +
  guides(color = "none") +
  labs(y = expression(Log[2]~fold~change~(A - C)),
       x = NULL,
       subtitle = "Batch 01")

# Right side plot
plot_6 <- sprg_log2fc %>% 
  filter(organism_n == 1, 
         condition == "06",
         comparison == "A_v_C",
         B > 5
  ) %>% 
  ggplot(aes(x = B, y = log2fc, color = organism_short)) +
  geom_point(size = 0.75, shape = 21) +
  geom_hline(data = group_median %>% filter(condition == "06", comparison == "A_v_C"), 
             aes(yintercept = log2fc, color = organism_short)) +
  geom_hline(data = ground_truth %>% filter(comparison == "A_v_C"),
             aes(yintercept = log2fc, color = organism_short), linetype = "dashed") +
  geom_text(data = group_median %>% filter(condition == "06", comparison == "A_v_C"), 
                  aes(x = 32, label = label), vjust = 0) +
  facet_grid(vars(comparison)) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(-10, 10)) +
  theme_minimal(base_size = 18) +
  theme(strip.text = element_blank(),
        axis.text.y = element_blank()) +
  guides(color = "none") +
  labs(y = NULL,
       x = NULL,
       subtitle = "Batch 06")

# Right side plot
plot_11 <- sprg_log2fc %>% 
  filter(organism_n == 1, 
         condition == "11",
         comparison == "A_v_C",
         B > 5
  ) %>% 
  ggplot(aes(x = B, y = log2fc, color = organism_short)) +
  geom_point(size = 0.75, shape = 21) +
  geom_hline(data = group_median %>% filter(condition == "11", comparison == "A_v_C"), 
             aes(yintercept = log2fc, color = organism_short)) +
  geom_hline(data = ground_truth %>% filter(comparison == "A_v_C"),
             aes(yintercept = log2fc, color = organism_short), linetype = "dashed") +
  geom_text(data = group_median %>% filter(condition == "11", comparison == "A_v_C"), 
                  aes(x = 32, label = label), vjust = 0) +
  facet_grid(vars(comparison)) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(-10, 10)) +
  theme_minimal(base_size = 18) +
  theme(strip.text = element_blank(),
        axis.text.y = element_blank()) +
  guides(color = "none") +
  labs(y = NULL,
       x = NULL,
       subtitle = "Batch 11")

# Right side plot
plot_neg20 <- sprg_log2fc %>% 
  filter(organism_n == 1, 
         condition == "neg20",
         comparison == "A_v_C",
         B > 5
  ) %>% 
  ggplot(aes(x = B, y = log2fc, color = organism_short)) +
  geom_point(size = 0.75, shape = 21) +
  geom_hline(data = group_median %>% filter(condition == "neg20", comparison == "A_v_C"), 
             aes(yintercept = log2fc, color = organism_short)) +
  geom_hline(data = ground_truth %>% filter(comparison == "A_v_C"),
             aes(yintercept = log2fc, color = organism_short), linetype = "dashed") +
  geom_text(data = group_median %>% filter(condition == "neg20", comparison == "A_v_C"), 
                  aes(x = 32, label = label), vjust = 0) +
  facet_grid(vars(comparison)) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(-10, 10)) +
  theme_minimal(base_size = 18) +
  theme(strip.text = element_blank(),
        axis.text.y = element_blank()) +
  guides(color = "none") +
  labs(y = NULL, 
       x = NULL, 
       subtitle = "Batch 06")


# Stitching plots together
plot_1 | plot_11 | plot_6 | plot_neg20 +
  plot_layout(guides = "collect")

# ggsave(filename = "figures/LFQBench_figure_stability_homogeneity.png", width = 10, height = 3)
ggsave(filename = "figures/LFQBench_figure_stability_homogeneity_2.png", width = 10, height = 3)



# Right side plot
sprg_log2fc %>% 
  filter(organism_n == 1, 
         condition == "neg20",
         comparison == "A_v_C",
         B > 5
  ) %>% 
  ggplot(aes(x = B, y = log2fc, color = organism_short)) +
  geom_point(size = 0.75) +
  geom_hline(data = group_median %>% filter(condition == "neg20",
                                            comparison == "A_v_C"), 
             aes(yintercept = log2fc, color = organism_short)) +
  geom_hline(data = ground_truth %>% filter(comparison == "A_v_C"),
             aes(yintercept = log2fc, color = organism_short), linetype = "dashed") +
  geom_text(data = group_median %>% filter(condition == "neg20",
                                           comparison == "A_v_C"), 
            aes(x = 32, label = label), 
            vjust = 0) +
  facet_grid(vars(comparison)) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(-10, 10)) +
  theme_minimal(base_size = 18) +
  theme(strip.text = element_blank(),
        axis.text.y = element_blank()) +
  guides(color = "none") +
  labs(y = NULL,
       x = NULL,
       subtitle = "Batch 06")

