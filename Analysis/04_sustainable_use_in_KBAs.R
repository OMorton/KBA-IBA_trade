library(tidyverse)
library(tidyterra)
library(terra)

options(scipen = 999)

sp_AOH <- read.csv("Data/KBA_AOH/AOH_list.csv") # 12143 AOHs
sp_IDs <- read.csv("Data/KBA_AOH/SISID_list.csv") %>% mutate(SIS.ID = as.character(SIS.ID)) # 10666 sp
KBAs <- read.csv("Data/KBA_AOH/KBA_list.csv") # 15861 KBAs

KBA_AOH <- read.csv("Data/KBA_AOH/output_all_correct_summary_AOH.csv") 

full_traded_taxonomy <- read.csv("Data/Trade/Traded_AOH_taxonomy_NOV24upd.csv") %>% 
  select(!c(multi_allNR, multi_sname, oldname_multi)) %>% mutate(AOH_SIS = as.character(AOH_SIS))#4258 sp
#Taxonomy <- read.csv("Data/Taxonomy/BL_Taxonomy_full.csv", na.strings=c("", " ","NA")) # 33212 sp + ssp (R and NR)
Trade_dat <- read.csv("Data/Trade/Hughes_et_al_trade_data.csv") # 15686 traded birds and mammals
#bird_list <- read.csv("D:/Data/Spatial/Bird_AOH_2022/Birds_list_AOH.csv")

Trade_dat %>% group_by(taxa, Traded) %>% tally() # 4265
full_traded_taxonomy %>% group_by(Traded) %>% tally() #4356 (6274 no traded)
## 10517 sp
## 10666 species to start with. Then 146 are removed because their AOH maps have no suitable habitat, 
##then three more are removed for other small errors. So that gives us the 10517
KBA_AOH_sum <- KBA_AOH %>%
  group_by(SIS.ID, Family, Order) %>% 
  summarise(area_inKBA = sum(as.numeric(Area_suitable_habitat_KBAs.ha.)),
            AOH_area.ha. = sum(AOH_area.ha.)) %>%
  mutate(SIS.ID = as.character(SIS.ID),
         perc_inKBA = area_inKBA/AOH_area.ha. *100) %>%
  left_join(full_traded_taxonomy, by = c("SIS.ID" = "AOH_SIS")) %>%
  ungroup() %>%
  mutate(TradedF = ifelse(is.na(Traded), "Not.Tr", "Tr"),
         AOH_area.ha.log = log10(AOH_area.ha.),
         AOH_area.ha.log.z = (AOH_area.ha.log - mean(AOH_area.ha.log) / sd(AOH_area.ha.log)))


## fix missing orders
KBA_AOH_sum_fix <- KBA_AOH_sum %>%
  mutate(Order = case_when(Sname == "Threnetes niger" ~ "CAPRIMULGIFORMES",
                           Sname == "Threnetes leucurus" ~ "CAPRIMULGIFORMES",
                           Sname == "Chalybura buffonii" ~ "CAPRIMULGIFORMES",
                           Sname == "Chalybura urochrysia" ~ "CAPRIMULGIFORMES",
                           Sname == "Eugenes fulgens" ~ "CAPRIMULGIFORMES",
                           Sname == "Macropygia tenuirostris" ~ "COLUMBIFORMES",
                           Sname == "Macropygia emiliana" ~ "COLUMBIFORMES",
                           Sname == "Gallinago paraguaiae" ~ "CHARADRIIFORMES",
                           Sname == "Accipiter fasciatus" ~ "ACCIPITRIFORMES",
                           Sname == "Accipiter cirrocephalus" ~ "ACCIPITRIFORMES",
                           Sname == "Pitta elegans" ~ "PASSERIFORMES",
                           Sname == "Zosterops chloris" ~ "PASSERIFORMES",
                           Sname == "Locustella castanea" ~ "PASSERIFORMES",
                           Sname == "Acrocephalus stentoreus" ~ "PASSERIFORMES",
                           Sname == "Acrocephalus australis" ~ "PASSERIFORMES",
                           Sname == "Phylloscopus sarasinorum" ~ "PASSERIFORMES",
                           Sname == "Carpodacus sibiricus" ~ "PASSERIFORMES",
                           Sname == "Rhynchospiza strigiceps" ~ "PASSERIFORMES",
                           Sname == "Melopyrrha portoricensis" ~ "PASSERIFORMES",
                           Sname == "Actenoides bougainvillei" ~ "CORACIIFORMES",
                           Sname == "Actenoides excelsus" ~ "CORACIIFORMES",
                           Sname == "Rhea pennata" ~ "STRUTHIONIFORMES",
                           Sname == "Rhea tarapacensis" ~ "STRUTHIONIFORMES",
                           Sname == "Cyornis banyumas" ~ "PASSERIFORMES",
                           Sname == "Otus magicus" ~ "STRIGIFORMES",
                           Sname == "Chlorostilbon poortmani" ~ "CAPRIMULGIFORMES",
                           Sname == "Crypturellus cinnamomeus" ~ "STRUTHIONIFORMES",
                           Sname == "Crypturellus occidentalis" ~ "STRUTHIONIFORMES",
                           Sname == "Bubo africanus" ~ "STRIGIFORMES",
                           Sname == "Glaucidium brasilianum" ~ "STRIGIFORMES",
                           Sname == "Glaucidium tucumanum" ~ "STRIGIFORMES",
                           Sname == "Otus manadensis" ~ "STRIGIFORMES",
                           Sname == "Tunchiornis ochraceiceps" ~ "PASSERIFORMES",
                           Sname == "Dicaeum sanguinolentum" ~ "PASSERIFORMES",
                           Sname == "Passerella schistacea" ~ "PASSERIFORMES",
                           Sname == "Passerella megarhyncha" ~ "PASSERIFORMES",
                           Sname == "Bradypterus baboecala" ~ "PASSERIFORMES",
                           Sname == "Bradypterus centralis" ~ "PASSERIFORMES",
                           Sname == "Gracupica contra" ~ "PASSERIFORMES",
                           Sname == "Alaudala rufescens" ~ "PASSERIFORMES",
                           Sname == "Oenanthe lugens" ~ "PASSERIFORMES",
                           .default = Order))

unique(KBA_AOH_sum_fix$Order)
KBA_AOH_sum_fix %>% group_by(Traded) %>% tally()
KBA_AOH_sum_fix %>% group_by(Order, Traded) %>% tally()
## 4039 (differs to the 4083 as some AOHs contain 0 AOH habitat)

### Sustainable and Unsustainable use in KBAs ####

## Qs
## NT or Thr threatend by trade
## LC with no decline traded
## NT or Thr with no major threat and no decline
## Second and third keep as a joint analysis.

## Get threat taxonomy data for traded species
taxonomy_wth_threat <- read.csv("Data/IUCN/Full_taxonomy_wth_threat.csv")
mini_thr_taxonomy <- taxonomy_wth_threat %>% select(AOH_SIS, Sname, Recognized, Traded,
                                                    BL_name, Scheffers_name, IUCN2023_v8name,
                                                    redlistCategory, populationTrend,
                                                    timing_sc, scope_sc, severity_sc, score_upd,
                                                    Trade_a_threat, Trade_not_threat, Trade_not_threat_nodecline,
                                                    Trade_thr_state)

## Match to traded-AOH data (n = 4033)
tr_thr_AOH <- KBA_AOH_sum_fix %>% filter(Traded == 1) %>% 
  left_join(mini_thr_taxonomy, by = c("Sname", "Traded", "BL_name", "Scheffers_name")) %>%
  ## remove the 6 taxa now lumped to other taxa
  filter(!is.na(redlistCategory))

## 1. Are threateend species less likely to be threatened by BRU inside KBAs
thr_sp_thr_by_trade <- tr_thr_AOH %>% filter(!redlistCategory %in%
                                               c("Data Deficient", "Least Concern")) %>%
  mutate(Trade_a_threat = case_when(Trade_a_threat == "Threat" ~ 1, 
                                    is.na(Trade_a_threat) ~ 0))
thr_sp_thr_by_trade %>% group_by(Trade_a_threat) %>% tally()
## 8 Orders removed (EURYPYGIFORMES, PODICIPEDIFORMES, CATHARTIFORMES,
## SULIFORMES, MESITORNITHIFORMES, MUSOPHAGIFORMES, TROGONIFORMES, PHOENICOPTERIFORMES)
## totalling 19 species to be removed from coefs (n = 874 species)
f <- thr_sp_thr_by_trade %>% group_by(Order) %>% tally()
thr_sp_thr_by_trade2 <- thr_sp_thr_by_trade %>% 
  #group_by(Order) %>% filter(n()>5) %>%
  mutate(AOH_area.ha.log = log10(AOH_area.ha.),
         AOH_area.ha.log.z = (AOH_area.ha.log - mean(AOH_area.ha.log) / sd(AOH_area.ha.log)))

Thrsp_Thr_mod <- brm(Trade_a_threat ~ AOH_area.ha.log.z + perc_inKBA + (1+perc_inKBA|Order), 
                     family = bernoulli(),
                     data = thr_sp_thr_by_trade2,
                     file = "Outputs/Models/ThrSp_ThrByTrade_KBA_perc_allorders.rds",
                     cores = 4, iter = 1000, chains = 4)

## 2. Are LC species more likely to be sustainably used with increase area in KBAs
all_nodecline <- tr_thr_AOH %>% 
  mutate(no_decline = case_when(Trade_not_threat_nodecline == "no decline" ~ 1, 
                                is.na(Trade_not_threat_nodecline) ~ 0))
all_nodecline %>% group_by(no_decline) %>% tally()
f <- all_nodecline %>% group_by(Order) %>% tally()
## GAVIIFORMES, LEPTOSOMIFORMES, OPISTHOCOMIFORMES,
## PHOENICOPTERIFORMES, EURYPYGIFORMES, CARIAMIFORMES,
## MESITORNITHIFORMES, COLIIFORMES
## 20 species  to be removed from coefs(n = 4013)
all_nodecline2 <- all_nodecline %>% 
  #group_by(Order) %>% filter(n()>5) %>%
  mutate(AOH_area.ha.log = log10(AOH_area.ha.),
         AOH_area.ha.log.z = (AOH_area.ha.log - mean(AOH_area.ha.log) / sd(AOH_area.ha.log)))

Allsp_NoDec_mod <- brm(no_decline ~ AOH_area.ha.log.z + perc_inKBA + (1+perc_inKBA|Order), 
                       family = bernoulli(),
                       data = all_nodecline2,
                       file = "Outputs/Models/AllSp_nodecline_KBA_perc_allorders.rds",
                       cores = 4, iter = 1000, chains = 4)


## plotting
newdat_thr <- data.frame(AOH_area.ha.log.z = median(thr_sp_thr_by_trade2$AOH_area.ha.log.z ),
                         perc_inKBA = seq(from = 0, to = 100, length.out = 101))

thrsp.draws <- add_epred_draws(object = Thrsp_Thr_mod, newdata = newdat_thr, re_formula = NA)

thrsp.draws.sum <- thrsp.draws %>% group_by(perc_inKBA) %>% median_hdci(.epred, .width = .9)
thrsp.draws.sum7 <- thrsp.draws %>% group_by(perc_inKBA) %>% median_hdci(.epred, .width = .7)

thrsp_thr_by_trade_plt <- ggplot(thr_sp_thr_by_trade2, aes(perc_inKBA, Trade_a_threat)) +
  geom_point(alpha = .2) +
  geom_ribbon(data = thrsp.draws.sum, aes(x = perc_inKBA, y = .epred, ymin = .lower, ymax = .upper),
              fill = "#b2182b") +
  geom_ribbon(data = thrsp.draws.sum7, aes(x = perc_inKBA, y = .epred, ymin = .lower, ymax = .upper),
              fill = "#d6604d") +
  geom_line(data = thrsp.draws.sum, aes(perc_inKBA, .epred), linewidth = 1) +
  xlab("% of AOH in KBA") +
  ylab("Probabilty of use \n being a threat") +
  theme_minimal(base_size = 10)

thr.coefs.draws <- as.data.frame(coef(Thrsp_Thr_mod, summary = FALSE)$Order[,,"perc_inKBA"]) %>%
  pivot_longer(everything(), names_to = "Order", values_to = "coef")

thr.coefs.draws.sum <- thr.coefs.draws %>% 
  filter(!Order %in% c("AVIIFORMES", "LEPTOSOMIFORMES", "OPISTHOCOMIFORMES",
                       "PHOENICOPTERIFORMES", "EURYPYGIFORMES", "CARIAMIFORMES",
                       "MESITORNITHIFORMES", "COLIIFORMES")) %>% 
  group_by(Order) %>% median_hdci(coef, .width = .9)
thr.coefs.draws.sum7 <- thr.coefs.draws %>% 
  filter(!Order %in% c("AVIIFORMES", "LEPTOSOMIFORMES", "OPISTHOCOMIFORMES",
                       "PHOENICOPTERIFORMES", "EURYPYGIFORMES", "CARIAMIFORMES",
                       "MESITORNITHIFORMES", "COLIIFORMES")) %>% 
  group_by(Order) %>% median_hdci(coef, .width = .7)

thrsp_thr_by_trade_coefs <- ggplot(thr.coefs.draws.sum, aes(coef, str_to_title(Order), xmin = .lower, xmax = .upper)) +
  geom_errorbarh(colour = "#b2182b", size = 2, height = 0)+
  geom_errorbarh(data = thr.coefs.draws.sum7, colour = "#d6604d", size = 2, height = 0) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Parameter values") +
  ylab("Order") +
  theme_minimal(base_size = 10)

newdat_sus <- data.frame(AOH_area.ha.log.z = median(all_nodecline2$AOH_area.ha.log.z ),
                         perc_inKBA = seq(from = 0, to = 100, length.out = 101))

sussp.draws <- add_epred_draws(object = Allsp_NoDec_mod, newdata = newdat_sus, re_formula = NA)

sussp.draws.sum <- sussp.draws %>%
  group_by(perc_inKBA) %>% median_hdci(.epred, .width = .9)
sussp.draws.sum7 <- sussp.draws %>% 
  group_by(perc_inKBA) %>% median_hdci(.epred, .width = .7)

nodecline_plt <- ggplot(all_nodecline2, aes(perc_inKBA, no_decline)) +
  geom_point(alpha = .2) +
  geom_ribbon(data = sussp.draws.sum, aes(x = perc_inKBA, y = .epred, ymin = .lower, ymax = .upper),
              fill = "#2166ac") +
  geom_ribbon(data = sussp.draws.sum7, aes(x = perc_inKBA, y = .epred, ymin = .lower, ymax = .upper),
              fill = "#4393c3") +
  geom_line(data = sussp.draws.sum, aes(perc_inKBA, .epred), linewidth = 1) +
  xlab("% of AOH in KBA") +
  ylab("Probabilty of use not negatively \n impacting extinction risk") +
  theme_minimal(base_size = 10)

sus.coefs.draws <- as.data.frame(coef(Allsp_NoDec_mod, summary = FALSE)$Order[,,"perc_inKBA"]) %>%
  pivot_longer(everything(), names_to = "Order", values_to = "coef")

sus.coefs.draws.sum <- sus.coefs.draws %>% 
  filter(!Order %in% c("AVIIFORMES", "LEPTOSOMIFORMES", "OPISTHOCOMIFORMES",
                       "PHOENICOPTERIFORMES", "EURYPYGIFORMES", "CARIAMIFORMES",
                       "MESITORNITHIFORMES", "COLIIFORMES")) %>% 
  group_by(Order) %>% median_hdci(coef, .width = .9)
sus.coefs.draws.sum7 <- sus.coefs.draws%>% 
  filter(!Order %in% c("AVIIFORMES", "LEPTOSOMIFORMES", "OPISTHOCOMIFORMES",
                       "PHOENICOPTERIFORMES", "EURYPYGIFORMES", "CARIAMIFORMES",
                       "MESITORNITHIFORMES", "COLIIFORMES"))  %>% 
  group_by(Order) %>% median_hdci(coef, .width = .7)

nodecline_coefs <- ggplot(sus.coefs.draws.sum, aes(coef, str_to_title(Order), xmin = .lower, xmax = .upper)) +
  geom_errorbarh(colour = "#2166ac", size = 2, height = 0)+
  geom_errorbarh(data = sus.coefs.draws.sum7, colour = "#4393c3", size = 2, height = 0) +
  geom_point() +
  ylab("Order") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Parameter values") +
  theme_minimal(base_size = 10)



## arrangement
library(ggpubr)
sus_thr_plt <- ggarrange(thrsp_thr_by_trade_plt, nodecline_plt, 
                         thrsp_thr_by_trade_coefs, nodecline_coefs,
                         heights = c(1, 1.5), labels = c("A.", "B.", "C.", "D."))

ggsave(path = "Outputs/draft_figs", sus_thr_plt, 
       filename = "Sus_thr_perc.png",  bg = "white",
       device = "png", width = 18, height = 18, units = "cm")

# additional plot

Thr_hist <- ggplot(thr_sp_thr_by_trade2, aes(perc_inKBA)) +
  geom_histogram(fill = "tomato") +
  facet_wrap(~Trade_a_threat, ncol = 1) +
  ylab("Count") + xlab("% of AOH in KBA") +
  theme_minimal()

NoDec_hist <- ggplot(all_nodecline2, aes(perc_inKBA)) +
  geom_histogram(fill = "dodgerblue") +
  facet_wrap(~no_decline, ncol = 1) +
  ylab("Count") + xlab("% of AOH in KBA") +
  theme_minimal()

hist_plt_arr <- ggarrange(Thr_hist, NoDec_hist, labels = c("A.", "B."),
                          ncol = 1)

ggsave(path = "Outputs/draft_figs", hist_plt_arr, 
       filename = "Thr_nodec_hists.png",  bg = "white",
       device = "png", width = 12, height = 15, units = "cm")


