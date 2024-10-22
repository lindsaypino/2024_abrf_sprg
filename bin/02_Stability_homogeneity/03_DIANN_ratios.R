
# libraries ---------------------------------------------------------------

library(tidyverse)
library(janitor)
library(rio)


# data import -------------------------------------------------------------

pep <- import("DIANN/SPRG_DIANN_analysis.tsv")

names(pep)

pep %>% 
  sample_n(1E4) %>% 
  ggplot() +
  geom_point(aes(x = Precursor.Quantity+10, y = Precursor.Normalised+10), alpha = 0.1) +
  facet_wrap(~Run) +
  scale_y_log10() +
  scale_x_log10()

table(pep$Run)
