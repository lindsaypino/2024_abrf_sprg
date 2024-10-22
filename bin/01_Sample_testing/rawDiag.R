

remotes::install_github('loosolab/MARMoSET', host='github.molgen.mpg.de/api/v3')

install.packages('http://fgcz-ms.uzh.ch/~cpanse/rawDiag_0.0.41.tar.gz', repo=NULL)

library(tidyverse)
library(rawDiag)
library(MARMoSET)




a <- read.raw("raw/dia_44min_ratio/220602_ABRF_ratio_A.raw", rawDiag = FALSE)
b <- read.raw("raw/dia_44min_ratio/220602_ABRF_ratio_B.raw", rawDiag = FALSE)
c <- read.raw("raw/dia_44min_ratio/220602_ABRF_ratio_C.raw", rawDiag = FALSE)
raw <- read.raw("../20210428_Ex480_Bartolomei_placenta/raw/dia/20210428-Ex480-Bartolomei-placenta_prot_1122-10_T_1_168.raw", rawDiag = FALSE)


data <- bind_rows(a, b, c)

data %>% 
  ggplot(aes(x = filename, y = log2(TIC))) +
  geom_boxplot()

raw %>% 
  # filter(MSOrder == "Ms") %>% 
  ggplot() +
  geom_line(aes(x = StartTime, y = AGCFill))


raw %>% 
  # filter(MSOrder == "Ms") %>% 
  ggplot() +
  geom_boxplot(aes(x = cut_interval(StartTime, 44), y = AGCFill))


ggplot(data) +
  geom_boxplot(aes(x = cut_interval(StartTime, 44), y = AGCFill)) +
  facet_grid(vars(filename), vars(MSOrder))

ggplot(data) +
  geom_boxplot(aes(x = cut_interval(StartTime, 44), y = IonInjectionTimems)) +
  facet_grid(vars(filename), vars(MSOrder))




names(data)


 
json <- generate_json(path_data = "raw/dia_44min_ratio/220602_ABRF_ratio_A.raw")


flat_json <- flatten_json(json = json)

create_term_match_table(instrument_list = c('Thermo EASY-nLC', 'Q Exactive - Orbitrap_MS'),
                        origin_key = 'jpr_guidelines_ms')

instrument_names(json, 1)

MARMoSET::testfile_json
