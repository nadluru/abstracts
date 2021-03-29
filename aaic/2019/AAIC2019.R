#' ---
#' title: "R code for the analysis conducted in 'Difference Analysis of DWI and Cortical Thickness in the Alzheimer's Disease Connectome Project (ADCP)'"
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
library(DataCombine)
library(ggplot2)
library(ggridges)
library(GGally)
library(lemon)

# Variables ====
rm(list = ls(all = TRUE))
csvroot = 'CSVs/'
figroot = 'Figures/'

#' # Data 
#' ## Basic demographics
basicdemo = read.csv(paste0(csvroot, 'BasicDemographics3.csv'),
                     na.strings = c('', ' '))
basicdemo$Subject.ID.Number %<>% as.factor
basicdemo$Sex %<>% trimws %>% as.factor

#' ## JHU DTI
jhudtibasicdemo = read.csv(paste0(csvroot, 'DTIMeasures_JHU_BasicDemo.csv'))
library(plyr)
jhudtibasicdemo$Consensus.Diagnosis %<>% revalue(c("Healthy Norm" = "CU"))
detach('package:plyr')
jhuboxplotinfo = jhudtibasicdemo %>%
  filter(MeasureName != 'ad') %>%
  mutate(AgeGroup = cut_interval(Age, 2)) %>%
  group_by(WhiteMatter, 
           MeasureName, 
           Sex, 
           AgeGroup, 
           Consensus.Diagnosis) %>%
  summarise(RMean = mean(Mean),
            RMedian = median(Mean),
            RIQR = IQR(Mean))

RMean_AD = jhuboxplotinfo %>% 
  filter(Consensus.Diagnosis == 'AD') %>% 
  spread(Consensus.Diagnosis, RMean) %>%
  mutate(RMean = `AD`)
RMean_MCI = jhuboxplotinfo %>% 
  filter(Consensus.Diagnosis == 'MCI') %>% 
  spread(Consensus.Diagnosis, RMean) %>%
  mutate(RMean = `MCI`)
RMean_HN = jhuboxplotinfo %>% 
  filter(Consensus.Diagnosis == 'CU') %>% 
  spread(Consensus.Diagnosis, RMean) %>% 
  mutate(RMean = `CU`)

DTIRMeanDiffs = RMean_AD
DTIRMeanDiffs$AD.minus.MCI = RMean_AD$RMean - RMean_MCI$RMean
DTIRMeanDiffs$AD.minus.CU = RMean_AD$RMean - RMean_HN$RMean
DTIRMeanDiffs$MCI.minus.CU = RMean_MCI$RMean - RMean_HN$RMean
DTIRMeanDiffs %<>% 
  gather(Group.Difference, 
         RMeanDifference, 
         AD.minus.MCI:MCI.minus.CU)

#' ## Freesurfer
fsthicknessdemo = read.csv(paste0(csvroot, 'FreeSurferThicknessDemo.csv')) %>%
  filter(Cortex != 'MeanThickness')
library(plyr)
fsthicknessdemo$Consensus.Diagnosis %<>% revalue(c("Healthy Norm" = "CU"))
detach('package:plyr')
fsboxplotinfo = fsthicknessdemo %>%
  mutate(AgeGroup = cut_interval(Age, 2)) %>%
  group_by(Cortex,
           Sex, 
           AgeGroup, 
           Consensus.Diagnosis) %>%
  summarise(RMean = mean(MeanThickness),
            RMedian = median(MeanThickness),
            RIQR = IQR(MeanThickness))

RMean_AD = fsboxplotinfo %>% 
  filter(Consensus.Diagnosis == 'AD') %>% 
  spread(Consensus.Diagnosis, RMean) %>%
  mutate(RMean = `AD`)
RMean_MCI = fsboxplotinfo %>% 
  filter(Consensus.Diagnosis == 'MCI') %>% 
  spread(Consensus.Diagnosis, RMean) %>%
  mutate(RMean = `MCI`)
RMean_HN = fsboxplotinfo %>% 
  filter(Consensus.Diagnosis == 'CU') %>% 
  spread(Consensus.Diagnosis, RMean) %>% 
  mutate(RMean = `CU`)

FSRMeanDiffs = RMean_AD
FSRMeanDiffs$AD.minus.MCI = RMean_AD$RMean - RMean_MCI$RMean
FSRMeanDiffs$AD.minus.CU = RMean_AD$RMean - RMean_HN$RMean
FSRMeanDiffs$MCI.minus.CU = RMean_MCI$RMean - RMean_HN$RMean
FSRMeanDiffs %<>% 
  gather(Group.Difference, 
         RMeanDifference, 
         AD.minus.MCI:MCI.minus.CU)

#' ## JHU + Freesurfer
jhufsboxplotinfo = merge(jhuboxplotinfo, 
                         fsboxplotinfo, 
                         by.x = c('Sex', 
                                  'AgeGroup', 
                                  'Consensus.Diagnosis'), 
                         by.y = c('Sex', 
                                  'AgeGroup', 
                                  'Consensus.Diagnosis'), 
                         all = TRUE)

DTIFSRMeanDiffs = merge(DTIRMeanDiffs, 
                        FSRMeanDiffs, 
                        by.x = c('Sex', 
                                 'AgeGroup', 
                                 'Group.Difference'), 
                        by.y = c('Sex', 
                                 'AgeGroup', 
                                 'Group.Difference'), 
                        all = TRUE)

#' # Plots 
#' ## Basic demographics
#+ fig.width=8.0, fig.height=5.0, warning=F
p = basicdemo %>% 
  mutate(Age.Group = cut_interval(Age, 2)) %>% 
  ggplot(aes(x = Sex, 
         fill = Consensus.Diagnosis)) + 
  geom_bar() + 
  geom_text(aes(label = ..count..),
            stat = 'count',
            position = position_stack(0.5),
            size = 8) +
  facet_rep_wrap(Age.Group~.) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        plot.title = element_text(size = 20, hjust = 0.5),
        panel.border = element_blank(), 
        axis.line = element_line()) +
  coord_capped_cart(bottom = 'both',
                    left = 'both') +
  labs(y = 'Sample size')
p
pdf(paste0(figroot, 'BasicDemo', '.pdf'),
    width = 8, height = 5)
print(p)
dev.off()

#+ fig.width=9.0, fig.height=4.5, warning=F
p = basicdemo %>% 
  mutate(Age.Group = cut_interval(Age, 2)) %>% 
  ggplot(aes(x = Age,
             y = Age.Group,
             color = Consensus.Diagnosis)) + 
  geom_density_ridges(aes(point_color = Consensus.Diagnosis),
                      jittered_points = TRUE,
                      point_alpha = 1,
                      alpha = 0.5, scale = 1.0) +
  facet_rep_wrap(Sex~.) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        plot.title = element_text(size = 20, hjust = 0.5),
        panel.border = element_blank(), 
        axis.line = element_line()) +
  coord_capped_cart(bottom = 'both',
                    left = 'both') +
  labs(x = 'Age [years]',
    y = 'Age group')
p
pdf(paste0(figroot, 'BasicDemoAge', '.pdf'),
    width = 9, height = 4.5)
print(p)
dev.off()

#+ fig.width=6.0, fig.height=5.0, warning=F
p = jhudtibasicdemo %>% select(SubjectID, Sex, Age, Consensus.Diagnosis) %>% unique %>% 
  mutate(Age.Group = cut_interval(Age, 2)) %>% 
  ggplot(aes(x = Sex, 
             fill = Consensus.Diagnosis)) + 
  geom_bar() + 
  geom_text(aes(label = ..count..),
            stat = 'count',
            position = position_stack(0.5),
            size = 8) +
  facet_rep_wrap(Age.Group~.) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.position = "bottom",
        panel.border = element_blank(), 
        axis.line = element_line()) +
  coord_capped_cart(bottom = 'both',
                    left = 'both') +
  labs(y = 'Sample size')
p
pdf(paste0(figroot, 'ImgBasicDemo', '.pdf'),
    width = 6, height = 5)
print(p)
dev.off()

#+ fig.width=6.0, fig.height=5.0, warning=F
p = jhudtibasicdemo %>% select(SubjectID, Sex, Age, Consensus.Diagnosis) %>% unique %>% 
  mutate(Age.Group = cut_interval(Age, 2)) %>% 
  ggplot(aes(x = Age,
             y = Age.Group,
             color = Consensus.Diagnosis)) + 
  geom_density_ridges(aes(point_color = Consensus.Diagnosis),
                      jittered_points = TRUE,
                      point_alpha = 1,
                      alpha = 0.5, scale = 1.0) +
  facet_rep_wrap(Sex~.) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.position = "bottom",
        panel.border = element_blank(), 
        axis.line = element_line()) +
  coord_capped_cart(bottom = 'both',
                    left = 'both') +
  labs(x = 'Age [years]',
       y = 'Age group')
p
pdf(paste0(figroot, 'ImgBasicDemoAge', '.pdf'),
    width = 6, height = 5)
print(p)
dev.off()

#' ## JHU DTI
#+ fig.width=15.0, fig.height=6.0, warning=F
wmplot = jhuboxplotinfo %>% 
  filter(RMean <= 1.1) %>%
  ggplot(aes(x = RMean, 
             y = AgeGroup,
             color = Consensus.Diagnosis)) + 
  geom_density_ridges(aes(point_color = Consensus.Diagnosis),
                      jittered_points = TRUE,
                      point_alpha = 1,
                      alpha = 0.5, scale = 1.0) +
  facet_rep_wrap(Sex~MeasureName, nrow = 2,
                 scales = "free") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        panel.border = element_blank(), 
        axis.line = element_line()) +
  coord_capped_cart(bottom = 'both',
                    left = 'both') +
  labs(x = 'DTI mean',
       y = 'Age group')
wmplot
pdf(paste0(figroot, 'DTIMean', '.pdf'),
    width = 15, height = 6)
print(wmplot)
dev.off()

#+ fig.width=16.5, fig.height=6.0, warning=F
wmplot = DTIRMeanDiffs %>% 
  ggplot(aes(x = RMeanDifference, 
             y = AgeGroup,
             color = Group.Difference)) + 
  geom_density_ridges(aes(point_color = Group.Difference),
                      jittered_points = TRUE,
                      point_alpha = 1,
                      alpha = 0.5, 
                      scale = 1.0) +
  geom_vline(xintercept = 0,
             linetype = 'dashed') +
  facet_rep_wrap(Sex~MeasureName, nrow = 2,
                 scales = "free") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        panel.border = element_blank(), 
        axis.line = element_line()) +
  coord_capped_cart(bottom = 'both',
                    left = 'both') +
  labs(x = 'DTI mean difference',
       y = 'Age group')
wmplot
pdf(paste0(figroot, 'DTIMeanDiffs', '.pdf'),
    width = 16.5, height = 6)
print(wmplot)
dev.off()

#+ fig.width=15.0, fig.height=6.0, warning=F
wmplot = jhuboxplotinfo %>% 
  ggplot(aes(x = RIQR/RMean, 
             y = AgeGroup,
             color = Consensus.Diagnosis)) + 
  geom_density_ridges(aes(point_color = Consensus.Diagnosis),
                      jittered_points = TRUE,
                      point_alpha = 1,
                      alpha = 0.5, scale = 1.0) +
  facet_rep_wrap(Sex~MeasureName, nrow = 2,
                 scales = "free") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        panel.border = element_blank(), 
        axis.line = element_line()) +
  coord_capped_cart(bottom = 'both',
                    left = 'both') +
  labs(x = 'Normalized IQR of DTI',
       y = 'Age group')
wmplot
pdf(paste0(figroot, 'DTIIQR', '.pdf'),
    width = 15, height = 6)
print(wmplot)
dev.off()

#+ fig.width=4.0, fig.height=6.0, warning=F
wmplot = jhuboxplotinfo %>% 
  filter(RMean <= 1.1 & MeasureName == 'fa') %>%
  ggplot(aes(x = RMean, 
             y = AgeGroup,
             color = Consensus.Diagnosis)) + 
  geom_density_ridges(aes(point_color = Consensus.Diagnosis),
                      jittered_points = TRUE,
                      point_alpha = 1,
                      alpha = 0.5, scale = 1.0) +
  facet_rep_wrap(Sex~MeasureName, nrow = 2) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        panel.border = element_blank(), 
        axis.line = element_line()) +
  guides(point_color = guide_legend(title.position = "top",
                                    title.hjust = 0.5,
                                    nrow = 2,
                                    byrow = TRUE)) +
  coord_capped_cart(bottom = 'both',
                    left = 'both') +
  labs(x = 'DTI mean',
       y = 'Age group')
wmplot
pdf(paste0(figroot, 'DTIMean_FA', '.pdf'),
    width = 4, height = 6)
print(wmplot)
dev.off()

#+ fig.width=4.0, fig.height=6.0, warning=F
wmplot = DTIRMeanDiffs %>% 
  filter(MeasureName == 'fa') %>%
  ggplot(aes(x = RMeanDifference, 
             y = AgeGroup,
             color = Group.Difference)) + 
  geom_density_ridges(aes(point_color = Group.Difference),
                      jittered_points = TRUE,
                      point_alpha = 1,
                      alpha = 0.5, 
                      scale = 1.0) +
  geom_vline(xintercept = 0,
             linetype = 'dashed') +
  facet_rep_wrap(Sex~MeasureName, nrow = 2) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        legend.position = "bottom",
        axis.line = element_line()) +
  guides(point_color = guide_legend(title.position = "top",
                                    title.hjust = 0.5,
                                    nrow = 2,
                                    byrow = TRUE)) +
  coord_capped_cart(bottom = 'both',
                    left = 'both') +
  labs(x = 'DTI mean difference',
       y = 'Age group')
wmplot
pdf(paste0(figroot, 'DTIMeanDiffs_FA', '.pdf'),
    width = 4, height = 6)
print(wmplot)
dev.off()

#' ## Freesurfer
#+ fig.width=4.0, fig.height=6.0, warning=F
gmplot = fsboxplotinfo %>%
  ggplot(aes(x = RMean, 
             y = AgeGroup,
             color = Consensus.Diagnosis)) + 
  geom_density_ridges(aes(point_color = Consensus.Diagnosis),
                      jittered_points = TRUE,
                      point_alpha = 1,
                      alpha = 0.5, 
                      scale = 1.0) +
  facet_rep_wrap(Sex~., nrow = 2,
                 scales = "free") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        legend.position = "bottom",
        axis.line = element_line()) +
  guides(point_color = guide_legend(title.position = "top",
                                    title.hjust = 0.5,
                                    nrow = 2,
                                    byrow = TRUE)) +
  coord_capped_cart(bottom = 'both',
                    left = 'both') +
  labs(x = 'Cortical thickness mean',
       y = 'Age group')
gmplot
pdf(paste0(figroot, 'CorticalThicknessMean', '.pdf'),
    width = 4, height = 6)
print(gmplot)
dev.off()

#+ fig.width=4.0, fig.height=6.0, warning=F
gmplot = FSRMeanDiffs %>% 
  ggplot(aes(x = RMeanDifference, 
             y = AgeGroup,
             color = Group.Difference)) + 
  geom_density_ridges(aes(point_color = Group.Difference),
                      jittered_points = TRUE,
                      point_alpha = 1,
                      alpha = 0.5, 
                      scale = 1.0) +
  geom_vline(xintercept = 0,
             linetype = 'dashed') +
  facet_rep_wrap(Sex~., nrow = 2,
                 scales = "free") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        legend.position = "bottom",
        axis.line = element_line()) +
  guides(point_color = guide_legend(title.position = "top",
                                    title.hjust = 0.5,
                                    nrow = 2,
                                    byrow = TRUE)) +
  coord_capped_cart(bottom = 'both',
                    left = 'both') +
  labs(x = 'Cortical thickness mean difference',
       y = 'Age group')
gmplot
pdf(paste0(figroot, 'CorticalThicknessMeanDiff', '.pdf'),
    width = 4, height = 6)
print(gmplot)
dev.off()

#+ fig.width=8.0, fig.height=6.0, warning=F
gmplot = fsboxplotinfo %>%
  ggplot(aes(x = RIQR/RMean, 
             y = AgeGroup,
             color = Consensus.Diagnosis)) + 
  geom_density_ridges(aes(point_color = Consensus.Diagnosis),
                      jittered_points = TRUE,
                      point_alpha = 1,
                      alpha = 0.5, 
                      scale = 1.0) +
  facet_rep_wrap(Sex~., nrow = 2,
                 scales = "free") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        panel.border = element_blank(), 
        axis.line = element_line()) +
  coord_capped_cart(bottom = 'both',
                    left = 'both') +
  labs(x = 'Normalized IQR of cortical thickness',
       y = 'Age group')
gmplot
pdf(paste0(figroot, 'CorticalThicknessIQR', '.pdf'),
        width = 8, height = 6)
print(gmplot)
dev.off()

#' ## JHU DTI + Freesurfer
p = jhufsboxplotinfo %>%
  filter(RMean.x <= 1.1) %>%
  group_by(Sex) %>% 
  do(
    plots = ggplot(data = ., 
                   aes(x = RMean.y, 
                       y = RMean.x, 
                       color = Consensus.Diagnosis)) + 
      geom_density2d() +
      scale_fill_gradient(low = "green", high = "red") +
      scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
      facet_rep_grid(MeasureName~AgeGroup, 
                     scales = "free") +
      theme(axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            strip.text = element_text(size = 20),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            plot.title = element_text(size = 20, hjust = 0.5),
            panel.border = element_blank(), 
            axis.line = element_line()) +
      coord_capped_cart(bottom = 'both',
                        left = 'both') +
      labs(x = 'Cortical thickness mean',
           y = 'DTI mean')
  ) %>%
  rowwise() %>%
  do(plotsanno = .$plots + 
       labs(title = paste(.$Sex)),
     filename = paste('DTIFSMean', .$Sex, 
                      sep = '_') %>% trimws
  )
#+ fig.width=10.0, fig.height=7.0, warning=F
p$plotsanno
p %>% 
  rowwise() %>%
  do(x = pdf(paste0(figroot, gsub("/", "", .$filename), '.pdf'),
             width = 10, height = 7),
     y = print(.$plotsanno),
     z = dev.off()
  )

p = DTIFSRMeanDiffs %>% 
  group_by(Sex) %>% 
  do(
    plots = ggplot(data = ., 
                   aes(x = RMeanDifference.y, 
                       y = RMeanDifference.x, 
                       color = Group.Difference)) + 
      geom_density2d() + 
      geom_vline(xintercept = 0,
                 linetype = 'dashed') +
      geom_hline(yintercept = 0,
                 linetype = 'dashed') +
      facet_rep_grid(MeasureName~AgeGroup, 
                     scales = "free_y") +
      theme(axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            strip.text = element_text(size = 20),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            plot.title = element_text(size = 20, hjust = 0.5),
            panel.border = element_blank(), 
            axis.line = element_line()) +
      coord_capped_cart(bottom = 'both',
                        left = 'both') +
      labs(x = 'Cortical thickness mean difference',
           y = 'DTI mean difference')
  ) %>%
  rowwise() %>%
  do(plotsanno = .$plots + 
       labs(title = paste(.$Sex)),
     filename = paste('DTIFSMeanDiff', .$Sex, 
                      sep = '_') %>% trimws
  )

#+ fig.width=10.0, fig.height=7.0, warning=F
p$plotsanno
p %>% 
  rowwise() %>%
  do(x = pdf(paste0(figroot, gsub("/", "", .$filename), '.pdf'),
             width = 10, height = 7),
     y = print(.$plotsanno),
     z = dev.off()
  )

p = jhufsboxplotinfo %>% 
  group_by(Sex) %>% 
  do(
    plots = ggplot(data = ., 
                   aes(x = RIQR.y/RMean.y, 
                       y = RIQR.x/RMean.x, 
                       color = Consensus.Diagnosis)) + 
      geom_density2d() + 
      facet_rep_grid(MeasureName~AgeGroup, 
                     scales = "free") +
      theme(axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            strip.text = element_text(size = 20),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            plot.title = element_text(size = 20, hjust = 0.5),
            panel.border = element_blank(), 
            axis.line = element_line()) +
      coord_capped_cart(bottom = 'both',
                        left = 'both') +
      labs(x = 'Normalized IQR of cortical thickness',
           y = 'Normalized IQR of DTI')
  ) %>%
  rowwise() %>%
  do(plotsanno = .$plots + 
       labs(title = paste(.$Sex)),
     filename = paste('DTIFSIQR', .$Sex, 
                      sep = '_') %>% trimws
  )

#+ fig.width=10.0, fig.height=7.0, warning=F
p$plotsanno
p %>% 
  rowwise() %>%
  do(x = pdf(paste0(figroot, gsub("/", "", .$filename), '.pdf'),
             width = 10, height = 7),
     y = print(.$plotsanno),
     z = dev.off()
  )