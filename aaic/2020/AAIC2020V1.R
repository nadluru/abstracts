#' ---
#' title: "R code for the analysis conducted in 'Geodesic path differences in the Alzheimer's disease connectome project (ADCP)'"
#' author: "Nagesh Adluru"
#' output:
#'    html_document:
#'      toc: true
#'      toc_float: true
#' ---
#'
#' # Initialization
# Loading the libraries ====
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(parallel)
library(scales)
library(grid)
library(gridExtra)
library(earth)
library(broom)

library(expm)
library(parallel)
library(magrittr)
library(tidyr)
library(DataCombine)
library(ggridges)
library(GGally)
library(lemon)
library(forcats)
library(stringr)
library(hexbin)
library(latex2exp)
library(ggExtra)
library(scales)
library(ggrepel)
library(tidyquant)
library(data.table)
library(tibbletime)
library(ggpubr)
library(RColorBrewer)
library(lme4)
library(blme)
library(purrr)
library(ggsci)
library(scales)
library(dplyr)
library(tidyr)
library(magrittr)
library(stringr)
library(latex2exp)
library(forcats)
library(broom)
library(purrr)
library(scales)
library(ggpol)
library(ggrepel)
library(lemon)
library(ggsignif)
library(gmodels)
library(margins)
library(lmPerm)
library(permuco)
library(permute)
library(modelr)
library(cowplot)
library(igraph)
library(ggraph)
library(seriation)

# tract count data basic trial 4 grouping into lobes fixed (poster production) =====
basicdemo = read.csv(paste0(csvroot, 'BasicDemographics3.csv'),
                     na.strings = c('', ' '))
basicdemo$Subject.ID.Number %<>% as.factor
basicdemo$Sex %<>% trimws %>% as.factor


connectomeids = read.table(paste0(csvroot, 'ConnectomeIDsList.txt'), header = F) %>% as.matrix
basicdemo %<>% filter(Subject.ID.Number %in% connectomeids)

basicdemo %<>% mutate(Connectomes = paste0('/Users/nadluru/ADCPConnectomes/MeanLengthConnectomesV3/meanlength_', 
                                           plyr::mapvalues(Consensus.Diagnosis, 'Healthy Norm', 'HC'), '_', 
                                           Subject.ID.Number, '_noweights.csv'))
basicdemo %<>% mutate(Connectomes2 = paste0('/Users/nadluru/ADCPConnectomes/WeightedTractCountConnectomes/connectome_', 
                                            plyr::mapvalues(Consensus.Diagnosis, 'Healthy Norm', 'HC'), '_', 
                                            Subject.ID.Number, '.csv'))
basicdemo %>% pull(Connectomes) %>% write.table(paste0(csvroot, 'ConnectomeListMacV3.txt'), row.names = F, col.names = F)
basicdemo %>% pull(Connectomes2) %>% write.table(paste0(csvroot, 'ConnectomeListMacV3_WTC.txt'), row.names = F, col.names = F)

nodes = read.csv('/Users/nadluru/IIT2019/LUT_GM_Desikan_COG_Mass.csv', as.is = T, header = T)
hemispheres = c('Left', 'Right')
lobes = nodes$Lobe %>% unique
fromLobes = map(6:2, ~ map((.x - 1):1, ~ paste0(lobes[[.y]], "_", lobes[[.x]]), .y = .x)) %>% unlist
toLobes = map(1:5, ~ map((.x + 1):6, ~ paste0(lobes[[.y]], "_", lobes[[.x]]), .y = .x)) %>% unlist
pb = progress_estimated(68)
distancedf = pmap(list(basicdemo$Connectomes,
                       basicdemo$Connectomes2,
                       basicdemo$Subject.ID.Number), ~ {
                         pb$tick()$print()
                         ml = read.csv(..1,
                                       sep = ",",
                                       header = F) %>% as.matrix
                         ml = ml + t(ml)
                         wtc = read.csv(..2,
                                        sep = ",",
                                        header = F) %>% as.matrix
                         wtc = wtc + t(wtc)
                         #dm = ml / wtc
                         net = graph_from_adjacency_matrix(wtc^(-1),
                                                           weighted = T,
                                                           diag = F,
                                                           mode = "undirected")
                         distancedf =  graph_from_adjacency_matrix(
                           distances(net) %>% as.matrix,
                           weighted = T,
                           diag = F,
                           mode = "undirected"
                         ) %>%
                           igraph::as_data_frame("edges") %>%
                           left_join(nodes %>% select(Vertex, Name, Hemisphere, Lobe),
                                     by = c("from" = "Vertex")) %>%
                           left_join(nodes %>% select(Vertex, Name, Hemisphere, Lobe),
                                     by = c("to" = "Vertex")) %>%
                           mutate(
                             HemHem = plyr::mapvalues(
                               paste0(Hemisphere.x, "_", Hemisphere.y),
                               from = "Right_Left",
                               to = "Left_Right"
                             ),
                             LobeLobe = plyr::mapvalues(paste0(Lobe.x, "_", Lobe.y),
                                                        from = fromLobes,
                                                        to = toLobes)
                           ) %>%
                           pivot_longer(c(HemHem, LobeLobe),
                                        names_to = 'RegionalClass',
                                        values_to = 'RegionalName') %>%
                           group_by(RegionalName) %>%
                           summarise(MeanDistance = mean(weight, na.rm = T)) %>%
                           mutate(ID = ..3)
                       }) %>%
  bind_rows %>%
  left_join(basicdemo, by = c('ID' = 'Subject.ID.Number'))

distancedf %>% write.csv(paste0(csvroot, 'DistancesRegionalV5.csv'),
                         row.names = F)

# Linear models trial 4 grouping into lobes fixed (poster production) ======
distancedf = read.csv(paste0(csvroot, 'DistancesRegionalV5.csv'))
mdls = distancedf %>%
  left_join(distancedf %>% 
              filter(Consensus.Diagnosis  %>% str_detect('Healthy') &
                       is.finite(MeanDistance)) %>% 
              group_by(RegionalName) %>% 
              summarise(MD = mean(MeanDistance, na.rm = T))) %>%
  mutate(MeanDistance = MeanDistance / MD) %>%
  filter(!str_detect(RegionalName, "Cerebell.*") &
           !is.infinite(MeanDistance) &
           !(RegionalName %in% map_chr(1:6, ~paste0(lobes[[.x]], "_", lobes[[.x]]))) &
           !(RegionalName %in% map2_chr(c(1,1,2),c(1,2,2),~paste0(hemispheres[[.x]], "_", hemispheres[[.y]])))) %>%
  group_by(RegionalName) %>%
  group_map(~ {
    mdl = lm(MeanDistance ~ Consensus.Diagnosis + Age + Sex,
             data = .x) %>% tidy %>% cbind(RegionalName = .y)
  }) %>%
  bind_rows

p = mdls %>% filter(str_detect(term, 'Consensus')) %>%
  mutate(term = gsub("Consensus.Diagnosis", "", term),
         term = plyr::mapvalues(term, from = c('Healthy Norm',
                                               'MCI'),
                                to = c('CU-AD',
                                       'MCI-AD')),
         RegionalName = fct_relevel(RegionalName,
                                    'Left_Left',
                                    'Left_Right',
                                    'Right_Right')) %>%
  ggplot(aes(x = RegionalName, y = estimate)) +
  geom_point() +
  geom_text(aes(label = ifelse(p.value < 0.05, "***", ""),
                y = estimate + std.error), vjust = 0) +
  geom_errorbar(aes(ymin = estimate - std.error,
                    ymax = estimate + std.error),
                width = 0.5) +
  geom_hline(yintercept = 0,
             linetype = 2) +
  facet_rep_grid(. ~ term, scales = "free") +
  gtheme +
  labs(x = "",
       y = TeX("$\\Delta$ \\[RGPL\\]")) +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        legend.title = element_blank(),
        panel.spacing = unit(0.000001, "lines"),
        axis.ticks.length = unit(0.2, "cm")) +
  guides(color = guide_legend(label.position = 'right'))
p
pdf(paste0(figroot, 'MainEffects_Poster', '.pdf'),
    width = 8.75,
    height = 4.75)
print(p)
dev.off()

visdf = distancedf %>%
  left_join(distancedf %>% 
              filter(Consensus.Diagnosis %>% str_detect('Healthy') &
                       is.finite(MeanDistance)) %>% 
              group_by(RegionalName) %>% 
              summarise(MD = mean(MeanDistance, na.rm = T))) %>%
  mutate(MeanDistance = MeanDistance / MD) %>%
  filter(!str_detect(RegionalName, "Cerebell.*") &
           !is.infinite(MeanDistance) &
           !(RegionalName %in% map_chr(1:6, ~paste0(lobes[[.x]], "_", lobes[[.x]]))) &
           !(RegionalName %in% map2_chr(c(1,1,2),c(1,2,2),~paste0(hemispheres[[.x]], "_", hemispheres[[.y]])))) %>%
  mutate(Consensus.Diagnosis = plyr::mapvalues(
    Consensus.Diagnosis,
    from = 'Healthy Norm',
    to = 'CU'
  ))
p = visdf %>%
  filter(MeanDistance <= 1.5) %>% 
  ggplot(aes(x = MeanDistance, 
             color = Consensus.Diagnosis)) + 
  geom_density() + 
  facet_rep_wrap(RegionalName ~ ., 
                 scales = "free_y",
                 ncol = 5) + 
  scale_x_continuous(breaks = seq(0.7,
                                  1.5,
                                  by = 0.1)) +
  gtheme +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        strip.text.x = element_text(size = 11),
        legend.title = element_blank(),
        axis.ticks.length = unit(0.2, "cm")) +
  labs(x = 'Relative geodesic path length (RGPL)',
       y = 'Sample density') +
  geom_vline(data = visdf %>% 
               group_by(RegionalName, 
                        Consensus.Diagnosis) %>% 
               summarise(MD2 = mean(MeanDistance, na.rm = T)), 
             aes(xintercept = MD2, 
                 color = Consensus.Diagnosis), 
             linetype = 2)
p
pdf(paste0(figroot, 'Distributions_RGND_Poster', '.pdf'),
    width = 10.0,
    height = 4.75)
print(p)
dev.off()

