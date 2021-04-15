#' ---
#' title: "R code for the analysis conducted in 'Analysis of the individual rates of change of diffusion tensor imaging measures from an accelerated longitudinal autism study'"
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
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)
library(scales)
library(lemon)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(blme)
library(purrr)
library(broom)
library(broom.mixed)
library(effects)
library(AICcmodavg)
library(wesanderson)
library(tidyverse)
library(Hotelling)

# Initializing variables ====
rm(list = ls(all = TRUE))
csvroot = 'CSVs/'
figroot = 'Figures/'

# ggplot theme =====
txtsize = 14
gtheme = theme(
  # legend
  legend.key = element_rect(colour = "black"),
  legend.title = element_blank(),
  legend.text = element_text(size = txtsize),
  legend.position = "top",
  legend.background = element_blank(),
  
  # text and axis
  strip.text.x = element_text(size = txtsize),
  strip.text.y = element_text(size = txtsize),
  axis.text = element_text(colour = "black",
                           size = txtsize),
  plot.title = element_text(size = txtsize, hjust = 0.5),
  axis.title = element_text(size = txtsize),
  axis.line = element_line(),
  strip.background = element_blank(),
  
  # panel
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                  colour = "gray"),
  panel.border = element_blank(),
  
  # ticks
  axis.ticks.length = unit(0.25, "cm")
)

#' # Loading data and basic filters with visits, age and regions
imgnmeta = read.csv(paste0(csvroot, 'MedianDTIDemoLong_TP1Corrected.csv')) %>% 
  filter(!is.na(GroupLabel) & !is.na(Age_Edit)) %>% 
  filter(!(ROIName %in% c('lT', 'rSCP', 'rT', 'F'))) %>% 
  filter(Age_Edit >= 10) %>% 
  group_by(ID, ROIName, DTIMeasureName) %>% 
  arrange(Age_Edit, .by_group = T) %>% 
  mutate(Start.Age = first(Age_Edit), TotalVisits = n()) %>% 
  ungroup() %>% filter(TotalVisits > 2) %>% 
  mutate(ID = as.factor(ID), 
         TotalVisits = as.factor(paste0('# visits: ', TotalVisits))) %>%
  filter(ROIName %in% 
           (read.csv(paste0(csvroot, 'JHU_labels_aa.csv')) %>% 
              filter(Include == 'Yes') %>% pull(Abbreviations)))

imgnmeta %<>% 
  anti_join(imgnmeta %>% group_by(ROIName) %>% 
              group_map(~{
                .x %>% filter((DTIMeasureName %in% c('MD', 'RD') & 
                                 DTIMeasureMedian > 0.85) | 
                                (DTIMeasureName == 'FA' & 
                                   DTIMeasureMedian < 0.25)) %>% 
                  select(ID, ROIName) %>% distinct
              }, .keep = T) %>% bind_rows)

#' # Accelerated longitudinal design (Figure 1)
dfSampleInfo = imgnmeta %>% 
  filter(ROIName == 'BCC' & DTIMeasureName == 'FA') %>% 
  group_by(TotalVisits, GroupLabel) %>% 
  summarise(n = n(), sub = length(unique(ID)), 
            .groups = 'drop') %>% 
  mutate(x = 5, y = 100)

#+ fig.width=7.0, fig.height=7.0, warning=F
p = imgnmeta %>% 
  filter(ROIName == 'BCC' & DTIMeasureName == 'FA') %>% 
  ggplot(aes(x = Age_Edit, y = fct_reorder(ID, Start.Age))) + 
  geom_line(aes(group = ID), alpha = 0.2) + 
  geom_point(size = 1, shape = 21, fill = NA) + 
  geom_text(data = dfSampleInfo, 
            aes(x = Inf, y = -Inf, 
                label = paste0('# samples: ', n)), 
            size = txtsize/3, hjust = 1, vjust = -1) + 
  geom_text(data = dfSampleInfo, 
            aes(x = Inf, y = -Inf, 
                label = paste0('# subjects: ', sub)), 
            size = txtsize/3, hjust = 1, vjust = -2.5) + 
  facet_rep_grid(GroupLabel~TotalVisits, scales = 'free_y') + 
  labs(x = 'Age [y]', y = 'Subjects ordered by age at visit 1') + 
  gtheme + 
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.x = element_line(size = 0.2, 
                                          linetype = 'solid', 
                                          color = 'gray'), 
        axis.ticks.length.y = unit(0.125, 'cm')) + 
  scale_x_continuous(breaks = pretty_breaks(n = 10))
p
ggsave(paste0(figroot, 'ALD_Final.pdf'), 
       p, 
       width = 7.0, 
       height = 7.0)

#' # Analysis
#' ## Subject level slopes and means from mixed models
slopedf = imgnmeta %>% 
  group_by(ROIName, DTIMeasureName) %>% 
  group_map(~{
    (lmer(DTIMeasureMedian ~ Age_Edit + (Age_Edit|ID), 
          data = .x) %>% coef)$ID %>% 
      rownames_to_column('ID') %>% 
      cbind(.y)
    }) %>% bind_rows %>% 
  left_join(imgnmeta %>% select(ID, GroupLabel) %>% distinct)

meandf = imgnmeta %>% 
  group_by(ROIName, DTIMeasureName) %>% group_map(~{
    (lmer(DTIMeasureMedian ~ (1|ID), 
          data = .x) %>% coef)$ID %>% 
      rownames_to_column('ID') %>% 
      cbind(.y)
    }) %>% bind_rows %>% 
  rename(meanDTI = `(Intercept)`) %>% 
  left_join(imgnmeta %>% select(ID,GroupLabel) %>% distinct)

slopemeandf = slopedf %>% inner_join(meandf)

#' ## Bhattacharyya Coefficient
nbins_Mean = 30
nbins_Slope = 30
r = raster::raster(xmn = min(slopemeandf$meanDTI), 
                   xmx = max(slopemeandf$meanDTI), 
                   ymn = min(slopemeandf$Age_Edit), 
                   ymx = max(slopemeandf$Age_Edit), 
                   ncol = nbins_Mean, nrow = nbins_Slope)
BCdf = slopemeandf %>% 
  group_by(ROIName, DTIMeasureName) %>% 
  group_map(~{
  asd = raster::rasterize(.x %>% filter(GroupLabel == 'ASD') %>% 
                            select(meanDTI, Age_Edit), r, fun = 'count')
  asdcount = .x %>% filter(GroupLabel == 'ASD') %>% nrow
  tdccount = .x %>% filter(GroupLabel == 'TDC') %>% nrow
  tdc = raster::rasterize(.x %>% filter(GroupLabel == 'TDC') %>% 
                            select(meanDTI, Age_Edit), r, fun = 'count')
  
cbind(bc = sum(sqrt((raster::values(tdc))/tdccount*(raster::values(asd))/asdcount), na.rm = T), .y)
}) %>% bind_rows

slopemeanbcdf = slopemeandf %>% inner_join(BCdf)

#' ## Univariate t tests
dfpval = slopemeandf %>% 
  rename(Slope = Age_Edit) %>% 
  inner_join(imgnmeta %>% select(ID, Age_Edit) %>% 
               group_by(ID) %>% 
               summarise(meanAge = mean(Age_Edit), 
                         .groups = 'drop')) %>% 
  group_by(ROIName, DTIMeasureName) %>% 
  group_map(~{
    lm(Slope ~ meanAge + GroupLabel, data = .x) %>% 
      tidy() %>% 
      filter(term == 'GroupLabelTDC') %>% 
      select(p.value) %>% 
      cbind(.y, DV='UniSlope') %>% 
      rbind(lm(meanDTI ~ meanAge + GroupLabel, data = .x) %>% 
              tidy() %>% filter(term == 'GroupLabelTDC') %>% 
              select(p.value) %>% cbind(.y, DV='UniMean'))
    }) %>% bind_rows %>% 
  mutate(padj = p.adjust(p.value, method = 'BH')) %>% 
  filter(padj <= 0.05)

#' ## Hotelling T2 test
dfbivariatepval = slopemeanbcdf %>% 
  group_by(ROIName, DTIMeasureName) %>% 
  group_map(~{
    cbind(pval = hotelling.test(.~GroupLabel, .x %>% 
                                  select(meanDTI, Age_Edit,
                                         GroupLabel))$pval, 
          .y, bcvis = .x$bc[[1]])}) %>% bind_rows %>% 
  mutate(padj = p.adjust(pval, method = 'BH')) %>% 
  filter(padj <= 0.05)

#' # Slopes vs. age (Figure 2)
#+ fig.width=9.45, fig.height=5.95, warning=F
p = dfpval %>% filter(DV == 'UniSlope') %>% 
  right_join(slopemeandf) %>% na.omit %>% 
  rename(Slope = Age_Edit) %>% 
  inner_join(imgnmeta %>% select(ID, Age_Edit) %>% 
               group_by(ID) %>% 
               summarise(meanAge = mean(Age_Edit), 
                         .groups = 'drop')) %>% 
  ggplot(aes(x = meanAge, 
             y = Slope, 
             color = GroupLabel)) + 
  geom_point(shape = 21, fill = NA) + 
  facet_rep_wrap(DTIMeasureName ~ fct_relevel(ROIName, 
                                              c('GCC', 
                                                'BCC', 
                                                'SCC', 
                                                'lALIC', 
                                                'rALIC', 
                                                'rPLIC')), 
                 scales = 'free', ncol = 4) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) + 
  labs(x = 'Mean age [y]', 
       y = 'Slope (dDTI/dAge)') + 
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5)) + 
  scale_color_brewer(palette = 'Set1') + 
  theme(panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) + 
  gtheme
p
ggsave(paste0(figroot, 'SlopeVsAge_Final.pdf'), 
       p, 
       width = 9.45, 
       height = 5.95)
#' # Means vs. age (Figure 3)
#+ fig.width=6.8, fig.height=7.2, warning=F
p = dfpval %>% 
  filter(DV == 'UniMean') %>% 
  select(ROIName, DTIMeasureName) %>% 
  distinct %>% arrange(DTIMeasureName) %>% 
  mutate(RD = paste0(ROIName, '-', DTIMeasureName)) %>% 
  filter(RD %in% c('GCC-FA', 
                   'BCC-FA', 
                   'SCC-FA', 
                   'GCC-MD', 
                   'lSLF-MD', 
                   'SCC-MD', 
                   'GCC-RD', 
                   'BCC-RD', 
                   'SCC-RD')) %>% 
  right_join(slopemeandf) %>% na.omit %>% 
  rename(Slope = Age_Edit) %>% 
  inner_join(imgnmeta %>% select(ID, Age_Edit) %>% 
               group_by(ID) %>% 
               summarise(meanAge = mean(Age_Edit), 
                         .groups = 'drop')) %>% 
  ggplot(aes(x = meanAge, 
             y = meanDTI, 
             color = GroupLabel)) + 
  geom_point(shape = 21, fill = NA) + 
  facet_rep_wrap(DTIMeasureName ~ fct_relevel(ROIName, 
                                              c('GCC', 
                                                'BCC', 
                                                'SCC', 
                                                'lSLF')), 
                 scales = 'free', ncol = 3) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) + 
  labs(x = 'Mean age [y]', 
       y = 'Mean (over the age range) DTI') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_color_brewer(palette = 'Set1') + 
  theme(panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) + 
  gtheme
p
ggsave(paste0(figroot, 'MeanVsAge_Final.pdf'), 
       p, 
       width = 6.8, 
       height = 7.2)

#' # Slopes vs. means plot (Figure 4)
#+ fig.width=12.85, fig.height=5.95, warning=F
p = dfbivariatepval %>% 
  right_join(slopemeanbcdf) %>% na.omit %>% 
  ggplot(aes(x = meanDTI, y = Age_Edit, color = GroupLabel)) + 
  stat_ellipse(type = 'norm', level = 0.9545) + 
  geom_point(shape = 21, fill = NA) + 
  geom_text(data = dfbivariatepval, 
            aes(x = Inf, 
                y = -Inf, 
                label = paste0('BC: ', round(bcvis, 3))), 
            hjust = 1, vjust = -1, inherit.aes = F) + 
  facet_rep_wrap(DTIMeasureName ~ fct_relevel(ROIName, 
                                              c('GCC', 
                                                'BCC', 
                                                'SCC', 
                                                'lALIC', 
                                                'rALIC', 
                                                'lSLF')), 
                 scales = 'free', ncol = 6) + 
  labs(x = 'Mean (over the age range) DTI', 
       y = 'Slope (dDTI/dAge)') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
  scale_color_brewer(palette = 'Set1') + 
  gtheme
p
ggsave(paste0(figroot, 'SlopeMean_Final.pdf'), 
       p, 
       width = 12.85, 
       height = 5.95)