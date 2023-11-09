
library(tidyverse)
library(rvest)
library(rredlist)

source("Functions.R")

## BOTW taxonomy
BOTW_taxo <- data.table::fread("Data/Taxonomy/BOTW_Taxonomy_clean.csv", header = TRUE)

## 11188 recognized sp
sp_list <- BOTW_taxo %>% filter(Recognized == "R") %>% select(cn, sn)  %>%
  rename("Species" = "sn", "CommonName" = "cn")


## scrape birdlife for hab and forest dependency (24s for 10 species, ~7hrs for all)
BL_scr_res <- scrape_birdlife(sp_list)

## scrape IUCN for status and pop (29s for 10 species, ~9hrs for all)
IUCN_scr_res <- retrieve_info(sp_list)

