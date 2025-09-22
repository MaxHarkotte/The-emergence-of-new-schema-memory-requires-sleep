# Harkotte et al (2025) - The emergence of new schema memory requires sleep
# Analyses of behavioral statistics, including summary statistics, plotting, inferential statistics
# last modified: August 2025
# maximilian.harkotte@gmail.com

rm(list = ls()) # clear workspace
cat("\014") # clear console

# 00 - Load toolboxes -----------------------------------------------------
library(psych)
library(ggplot2)
library(svglite)
library(ggpubr)
library(tidyverse)
library(lme4)
library(tibble)
library(cocor)
library(ggsignif)
library(rstatix) 

# Add aesthetics for later plotting
limits = aes(ymax = mean + (se), ymin = mean - (se)) # (1.96*se) for confidence intervals
dodge = position_dodge(width = 0.8)

# 01 - Data paths ---------------------------------------------------------
dataPath <- "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/" 
# dataPath <- "/Volumes/born_animal/Max/02_Schema_memory_formation/02_Schema_during_sleep/"

setwd(dataPath)

# 02 - Import data --------------------------------------------------------
test_wide <-
  read.csv2(
    "02-VideoScorings/02-Tables/01-DataSummary/00-TestClean.csv",
    header = TRUE,
    sep = ";",
    stringsAsFactors = TRUE
  )

enc_wide <-
  read.csv2(
    "02-VideoScorings/02-Tables/01-DataSummary/00-EncodingClean.csv",
    header = TRUE,
    sep = ";",
    stringsAsFactors = TRUE
  )

test_rearing <-
  read.csv2(
    "02-VideoScorings/02-Tables/01-DataSummary/00-TestClean_Rearing.csv",
    header = TRUE,
    sep = ";",
    stringsAsFactors = TRUE
  )

sleep <-
  read.csv2(
    "02-VideoScorings/02-Tables/01-DataSummary/00-Sleep.csv",
    header = TRUE,
    sep = ",",
    stringsAsFactors = TRUE
  )

# 03 - DRs main experiment  ------------------------------------------
test_schema <- subset(test_wide, test_wide$Type == "Schema")

Cum_DR_schema <-
  subset(
    test_schema,
    select = c(
      "Animal",
      "Type",
      "Retention",
      "Cum_DiRa_min_1",
      "Cum_DiRa_min_2",
      "Cum_DiRa_min_3",
      "Cum_DiRa_min_4",
      "Cum_DiRa_min_5"
    )
  )

Cum_DR_schema <-
  pivot_longer(
    Cum_DR_schema,
    cols = 4:8 ,
    names_to = "Minute",
    values_to = "DR"
  )

Cum_DR_schema$Minute <- as.factor(Cum_DR_schema$Minute)

Cum_DR_sum_schema = describeBy(
  Cum_DR_schema$DR,
  list(Cum_DR_schema$Retention, Cum_DR_schema$Minute),
  mat = TRUE,
  digits = 2
)

# Overall statistics
basic_all_hlm <-
  lmer(DR ~ Minute + Retention + (1 | Animal),
       data = Cum_DR_schema,
       REML = FALSE)

basic_ret_hlm <-
  lmer(DR ~ Minute + (1 | Animal),
       data = Cum_DR_schema,
       REML = FALSE)

basic_min_hlm <-
  lmer(DR ~ Retention + (1 | Animal),
       data = Cum_DR_schema,
       REML = FALSE)

interaction_all_hlm <-
  lmer(DR ~ Minute  * Retention + (1 | Animal),
       data = Cum_DR_schema,
       REML = FALSE)


anova(basic_ret_hlm, basic_all_hlm) # test main effect retention
anova(basic_min_hlm, basic_all_hlm) # test main effect min
anova(interaction_all_hlm, basic_all_hlm) # test interaction effect with minute

# Post-hoc Tests
stat.test1 <- Cum_DR_schema %>%
  group_by(Retention, Minute) %>%
  t_test(DR ~ 1, mu = 0) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('####', '###', '##', '#' ,  'ns')) %>%
  mutate(group1 = Retention, group2 = Minute) %>%
  add_xy_position(x = "Minute",
                  group = "Retention",
                  step.increase = 0) %>%
  mutate(y.position = 1.1, xmax = xmax -0.2)

stat.test2 <- Cum_DR_schema %>%
  group_by(Minute) %>%
  t_test(DR ~ Retention) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('****', '***', '**', '*' , 'ns'))  %>%
  add_xy_position(x = "Minute",
                  dodge = 0.8,
                  step.increase = 0)%>%
  mutate(y.position = 0.5)

# Plotting
DR_Cum_plot_schema <-
  ggplot(data = Cum_DR_sum_schema, aes(x = group2, y = mean, fill = group1)) +
  geom_bar(
    stat = 'identity',
    position = dodge,
    width = .8,
    colour = "black",
    linewidth = .2
  ) +
  geom_errorbar(limits, position = dodge, width = 0.3,
                linewidth = .2) +
  geom_jitter(
    data = Cum_DR_schema,
    aes(x = Minute, y = DR, fill = Retention),
    shape = 21,  # Hollow point with fill
    size = 0.8,
    alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)
  ) +
  stat_pvalue_manual(
    stat.test1,
    label = "{p.signif}",
    x = "xmax",
    remove.bracket = TRUE,
    hide.ns = TRUE,
    size = 2) +
  stat_pvalue_manual(
    stat.test2,
    label = "{p.signif}",
    tip.length = 0.01,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  ) +
  theme_classic() +
  scale_y_continuous(name = "Discrimination ratio",
                     breaks = seq(-1, 1, 0.5),
                     limits = c(-1, 1.3)) +
  scale_x_discrete(
    name = "Minutes",
    labels = c("1", "2", "3", "4", "5"),
    limits = c(
      "Cum_DiRa_min_1",
      "Cum_DiRa_min_2",
      "Cum_DiRa_min_3",
      "Cum_DiRa_min_4",
      "Cum_DiRa_min_5"
    )
  ) +
  scale_fill_manual(
    name = "Retention condition",
    labels = c(paste0("Sleep (N = 10)"),
               paste0("Sleep deprivation (N = 11)")),
    limits = c("Sleep", "Wake"),
    values = c("Sleep" = "#595959", "Wake" = "#e6e6e6ff")
  ) +
  geom_hline(colour = "black",
             yintercept = 0,
             size = .1)


DR_Cum_plot_schema <- DR_Cum_plot_schema +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = "none")          
DR_Cum_plot_schema

# Save figure 
ggsave(
  file = "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/04-Figures/schema_DR.svg",
  plot = DR_Cum_plot_schema,
  width = 55,
  height = 35,
  units = "mm"
)


# 04 - Control Parameters test - main experiment --------------------------
# Total exploration time 
test_expl_time <- subset(test_wide, test_wide$Type == "Schema")

test_expl_time <-
  subset(
    test_expl_time,
    select = c(
      "Animal",
      "Type",
      "Retention",
      "Total_exp_time"
    )
  )


test_expl_time_sum = describeBy(
  test_expl_time$Total_exp_time,
  list(test_expl_time$Retention),
  mat = TRUE,
  digits = 2
)

# Comparisons

stat.test1 <- test_expl_time %>%
  t_test(Total_exp_time ~ Retention) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('****', '***', '**', '*' , 'ns'))  %>%
  add_xy_position(x = "Retention",
                  step.increase = 0)%>%
  mutate(y.position = y.position + 8)

test_expl_time_plot <-
  ggplot(data = test_expl_time_sum, aes(x = group1, y = mean, fill = group1)) +
  geom_bar(
    stat = 'identity',
    position = dodge,
    width = .8,
    colour = "black",
    linewidth = .2
  ) +
  geom_errorbar(limits, position = dodge, width = 0.3,
                linewidth = .2) +
  geom_point(
    data = test_expl_time,
    aes(x = Retention, y = Total_exp_time, fill = Retention),
    shape = 21,  
    size = 0.8,
    alpha = 0.6
  ) +
  stat_pvalue_manual(
    stat.test1,
    label = "{p.signif}",
    tip.length = 0.01,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  ) +
  theme_classic() +
  scale_y_continuous(name = "Exploration time (s)") +
  scale_x_discrete(
    name = NULL, 
    labels = NULL,
    limits = c("Sleep", "Wake")
  ) +
  scale_fill_manual(
    values = c("Sleep" = "#595959", "Wake" = "#e6e6e6ff")
  ) 


test_expl_time_plot <- test_expl_time_plot +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = "none")          

test_expl_time_plot

# Disctance traveled 

test_dist <- subset(test_wide, test_wide$Type == "Schema")

test_dist <-
  subset(
    test_dist,
    select = c(
      "Animal",
      "Type",
      "Retention",
      "Cum_dist_min_5"
    )
  )


test_dist_sum = describeBy(
  test_dist$Cum_dist_min_5,
  list(test_dist$Retention),
  mat = TRUE,
  digits = 2
)

# Comparisons

stat.test1 <- test_dist %>%
  t_test(Cum_dist_min_5 ~ Retention) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('****', '***', '**', '*' , 'ns'))  %>%
  add_xy_position(x = "Retention",
                  step.increase = 0)%>%
  mutate(y.position = y.position + 8)

test_dist_plot <-
  ggplot(data = test_dist_sum, aes(x = group1, y = mean, fill = group1)) +
  geom_bar(
    stat = 'identity',
    position = dodge,
    width = .8,
    colour = "black",
    linewidth = .2
  ) +
  geom_errorbar(limits, position = dodge, width = 0.3,
                linewidth = .2) +
  geom_point(
    data = test_dist,
    aes(x = Retention, y = Cum_dist_min_5, fill = Retention),
    shape = 21,  
    size = 0.8,
    alpha = 0.6
  ) +
  stat_pvalue_manual(
    stat.test1,
    label = "{p.signif}",
    tip.length = 0.01,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  ) +
  theme_classic() +
  scale_y_continuous(name = "Distance (m)") +
  scale_x_discrete(
    name = NULL, 
    labels = NULL,
    limits = c("Sleep", "Wake")
  ) +
  scale_fill_manual(
    values = c("Sleep" = "#595959", "Wake" = "#e6e6e6ff")
  ) 


test_dist_plot<- test_dist_plot +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = "none")          

test_dist_plot


# Rearing
Cum_test_rear <-
  subset(
    test_rearing,
    select = c(
      "Animal",
      "Type",
      "Retention",
      "Cum_rear_exp_5"
    )
  )

Cum_test_rear_sum = describeBy(
  Cum_test_rear$Cum_rear_exp_5,
  list(Cum_test_rear$Retention),
  mat = TRUE,
  digits = 2
)

# Comparisons
anova_test_rear <- aov(Cum_rear_exp_5 ~ Retention , data = Cum_test_rear)
summary(anova_test_rear)

stat.test1 <- Cum_test_rear %>%
  t_test(Cum_rear_exp_5 ~ Retention) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('****', '***', '**', '*' , 'ns'))  %>%
  add_xy_position(x = "Retention",
                  step.increase = 0)%>%
  mutate(y.position = y.position + 8)

Cum_test_rear_plot <-
  ggplot(data = Cum_test_rear_sum, aes(x = group1, y = mean, fill = group1)) +
  geom_bar(
    stat = 'identity',
    position = dodge,
    width = .8,
    colour = "black",
    linewidth = .2
  ) +
  geom_errorbar(limits, position = dodge, width = 0.3,
                linewidth = .2) +
  geom_point(
    data = Cum_test_rear,
    aes(x = Retention, y = Cum_rear_exp_5, fill = Retention),
    shape = 21,  
    size = 0.8,
    alpha = 0.6,
    color = "black"
  ) +
  stat_pvalue_manual(
    stat.test1,
    label = "{p.signif}",
    tip.length = 0.01,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  ) +
  theme_classic() +
  scale_y_continuous(name = "Rearing time (s)") +
  scale_x_discrete(
    name = NULL,
    labels = NULL,
    limits = c("Sleep", "Wake")
  ) +
  scale_fill_manual(
    values = c("Sleep" = "#595959", "Wake" = "#e6e6e6ff")
  ) 

Cum_test_rear_plot <- Cum_test_rear_plot +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = "none")          

Cum_test_rear_plot

# save figure 
# ggsave(
#   file = "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/04-Figures/test_main_expl.svg",
#   plot = test_expl_time_plot,
#   width = 16,
#   height = 28,
#   units = "mm"
# )
# 
# ggsave(
#   file = "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/04-Figures/test_main_dist.svg",
#   plot = test_dist_plot,
#   width = 14,
#   height = 28,
#   units = "mm"
# )
# 
# ggsave(
#   file = "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/04-Figures/test_main_rear.svg",
#   plot = Cum_test_rear_plot,
#   width = 15,
#   height = 28,
#   units = "mm"
# )
# 



# 05 - Comparing DRs with encoding - main experiment ----------------------------------------
enc_schema <- subset(enc_wide, enc_wide$Type == "Schema")
enc_schema$Sampling <- as.factor(as.numeric(enc_schema$Sampling))

# DRs for first and last encoding episode (sleep condition)
enc_trials <- subset(enc_schema, enc_schema$Sampling == c("1", "8"))
enc_trials <- subset(enc_trials, enc_trials$Retention == "Sleep")

enc_trials <-
  subset(
    enc_trials,
    select = c(
      "Animal",
      "Type",
      "Sampling",
      "Cum_DiRa_min_2",
      "Cum_DiRa_min_5"
    )
  )

enc_trials <-
  pivot_longer(
    enc_trials,
    cols = 4:5 ,
    names_to = "Minute",
    values_to = "DR"
  )

enc_trials$Minute <- as.factor(enc_trials$Minute)
enc_trials$Sampling <- as.factor(as.numeric(enc_trials$Sampling))

enc_trials_sum = describeBy(
  enc_trials$DR,
  list(enc_trials$Minute, enc_trials$Sampling),
  mat = TRUE,
  digits = 2
)

# Comparisons
stat.test1 <- enc_trials %>%
  group_by(Minute, Sampling) %>%
  t_test(DR ~ 1, mu = 0) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('####', '###', '##', '#' ,  'ns')) %>%
  mutate(group1 = Minute, group2 = Sampling) %>%
  add_xy_position(x = "Sampling",
                  group = "Minute",
                  step.increase = 0) %>%
  mutate(y.position = 1.1, xmax = xmax -0.2)

stat.test2 <- enc_trials %>%
  group_by(Sampling) %>%
  t_test(DR ~ Minute) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('****', '***', '**', '*' , 'ns'))  %>%
  add_xy_position(x = "Sampling",
                  dodge = 0.8,
                  step.increase = 0)%>%
  mutate(y.position = 0.5)

enc_trials_plot<-
  ggplot(data = enc_trials_sum, aes(x = group2, y = mean, fill = group1)) +
  geom_bar(
    stat = 'identity',
    position = dodge,
    width = .8,
    colour = "black",
    linewidth = .2
  ) +
  geom_errorbar(limits, position = dodge, width = 0.3,
                linewidth = .2) +
  geom_jitter(
    data = enc_trials,
    aes(x = Sampling, y = DR, fill = Minute),
    shape = 21,  # Hollow point with fill
    size = 0.8,
    alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)
  ) +
  stat_pvalue_manual(
    stat.test1,
    label = "{p.signif}",
    x = "xmax",
    remove.bracket = TRUE,
    hide.ns = TRUE,
    size = 2) +
  
  theme_classic() +
  scale_y_continuous(name = "Discrimination ratio",
                     breaks = seq(-1, 1, 0.5),
                     limits = c(-1, 1.3)) +
  scale_x_discrete(
    name = "Encoding Episode",
    labels = c("First", "Last"),
    limits = c(
      "1",
      "8"
    )
  ) +
  scale_fill_manual(
    name = "Minute",
    labels = c("2","5"),
    limits = c("Cum_DiRa_min_2", "Cum_DiRa_min_5"),
    values = c("Cum_DiRa_min_2" = "#D40202", "Cum_DiRa_min_5" = "black")
  ) +
  geom_hline(colour = "black",
             yintercept = 0,
             size = .1)

enc_trials_plot <- enc_trials_plot +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = "none")          
enc_trials_plot

# save plot
ggsave(
  file = "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/04-Figures/enc_plot.svg",
  plot = enc_trials_plot,
  width = 32,
  height = 33,
  units = "mm"
)

# Correlation between last encoding DR and test DR (sleep condition)
recall <- subset(test_schema, test_schema$Retention == "Sleep")
last_enc <- subset(enc_schema, enc_schema$Sampling == "8")
last_enc <- subset(last_enc, last_enc$Type == "Schema")
last_enc <- subset(last_enc, last_enc$Retention == "Sleep")

DR_corr <-
  merge(recall, last_enc, by.x = "Animal", by.y = "Animal")

DR_corr_5 <- subset(DR_corr, select = c('Animal', 'Cum_DiRa_min_5.x', 'Cum_DiRa_min_5.y'))
DR_corr_2 <- subset(DR_corr, select = c('Animal', 'Cum_DiRa_min_2.x', 'Cum_DiRa_min_2.y'))


DR_corr_5  <- DR_corr_5  %>%
  mutate(Minute = "5")%>%
  rename(DR_Test = Cum_DiRa_min_5.x,
         DR_Encoding= Cum_DiRa_min_5.y)

DR_corr_2  <- DR_corr_2  %>%
  mutate(Minute = "2")%>%
  rename(DR_Test = Cum_DiRa_min_2.x,
         DR_Encoding= Cum_DiRa_min_2.y)

DR_corr <- rbind(DR_corr_2, DR_corr_5)

ggscatter(
  DR_corr,
  x = "DR_Encoding",
  y = "DR_Test",
  color = "Minute",
  shape = 21,
  size = 1,
  add = "reg.line",
  add.params = list(aes(color = Minute), fill = "lightgray", size = 0.5),
  conf.int = TRUE,
  cor.coef = TRUE,
  facet.by = "Minute",  
  cor.coef.method = "spearman",   
  xlab = "DR last episode",
  ylab = "DR Retrieval",
  palette = c("black", "#D40202")
)

DR_corr_plot <- ggscatter(
  DR_corr,
  x = "DR_Encoding",
  y = "DR_Test",
  fill = "Minute",
  shape = 21,
  size = 0.8,
  add = "reg.line",
  add.params = list(color = "Minute", fill = "lightgray", size = 0.5),
  conf.int = TRUE,
  cor.coef = FALSE,
  cor.coef.method = "spearman",   
  xlab = "DR last episode",
  ylab = "DR Retrieval",
  palette = c("black", "#D40202")
)

DR_corr_plot <- DR_corr_plot +
  theme_classic(base_size = 7) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7),
   legend.position = "none"  # <-- removes the legend
  )

DR_corr_plot

# Save plot
ggsave(
  file = "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/04-Figures/corr_plot.svg",
  plot = DR_corr_plot,
  width = 32,
  height = 33,
  units = "mm"
)

# 06 - DRs control experiment  --------------------------------------
test_hab <- subset(test_wide, test_wide$Type == "Interference")

Cum_DR_hab <-
  subset(
    test_hab,
    select = c(
      "Animal",
      "Type",
      "Retention",
      "Cum_DiRa_min_1",
      "Cum_DiRa_min_2",
      "Cum_DiRa_min_3",
      "Cum_DiRa_min_4",
      "Cum_DiRa_min_5"
    )
  )

Cum_DR_hab <-
  pivot_longer(
    Cum_DR_hab,
    cols = 4:8 ,
    names_to = "Minute",
    values_to = "DR"
  )

Cum_DR_hab$Minute <- as.factor(Cum_DR_hab$Minute)

Cum_DR_sum_hab = describeBy(
  Cum_DR_hab$DR,
  list(Cum_DR_hab$Retention, Cum_DR_hab$Minute),
  mat = TRUE,
  digits = 2
)

# Overall statistics
basic_all_hlm <-
  lmer(DR ~ Minute + Retention + (1 | Animal),
       data = Cum_DR_hab,
       REML = FALSE)

basic_ret_hlm <-
  lmer(DR ~ Minute + (1 | Animal),
       data = Cum_DR_hab,
       REML = FALSE)

basic_min_hlm <-
  lmer(DR ~ Retention + (1 | Animal),
       data = Cum_DR_hab,
       REML = FALSE)

interaction_all_hlm <-
  lmer(DR ~ Minute  * Retention + (1 | Animal),
       data = Cum_DR_hab,
       REML = FALSE)


anova(basic_ret_hlm, basic_all_hlm) # test main effect retention
anova(basic_min_hlm, basic_all_hlm) # test main effect min
anova(interaction_all_hlm, basic_all_hlm) # test interaction effect with minute

# Post-hoc Tests
stat.test1 <- Cum_DR_hab %>%
  group_by(Retention, Minute) %>%
  t_test(DR ~ 1, mu = 0) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('####', '###', '##', '#' ,  'ns')) %>%
  mutate(group1 = Retention, group2 = Minute) %>%
  add_xy_position(x = "Minute",
                  group = "Retention",
                  step.increase = 0) %>%
  mutate(y.position = 1.1, xmax = xmax -0.2)

stat.test2 <- Cum_DR_hab %>%
  group_by(Minute) %>%
  t_test(DR ~ Retention) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('****', '***', '**', '*' , 'ns'))  %>%
  add_xy_position(x = "Minute",
                  dodge = 0.8,
                  step.increase = 0)%>%
  mutate(y.position = 0.5)

# Plotting
DR_Cum_plot_hab <-
  ggplot(data = Cum_DR_sum_hab, aes(x = group2, y = mean, fill = group1)) +
  geom_bar(
    stat = 'identity',
    position = dodge,
    width = .8,
    colour = "black",
    linewidth = .2
  ) +
  geom_errorbar(limits, position = dodge, width = 0.3,
                linewidth = .2) +
  geom_jitter(
    data = Cum_DR_hab,
    aes(x = Minute, y = DR, fill = Retention),
    shape = 21,  # Hollow point with fill
    size = 0.8,
    alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)
  ) +
  stat_pvalue_manual(
    stat.test1,
    label = "{p.signif}",
    x = "xmax",
    remove.bracket = TRUE,
    hide.ns = TRUE,
    size = 2) +
  stat_pvalue_manual(
    stat.test2,
    label = "{p.signif}",
    tip.length = 0.01,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  ) +
  theme_classic() +
  scale_y_continuous(name = "Discrimination ratio",
                     breaks = seq(-1, 1, 0.5),
                     limits = c(-1, 1.3)) +
  scale_x_discrete(
    name = "Minutes",
    labels = c("1", "2", "3", "4", "5"),
    limits = c(
      "Cum_DiRa_min_1",
      "Cum_DiRa_min_2",
      "Cum_DiRa_min_3",
      "Cum_DiRa_min_4",
      "Cum_DiRa_min_5"
    )
  ) +
  scale_fill_manual(
    name = "Retention condition",
    labels = c(paste0("Sleep (N = 10)"),
               paste0("Sleep deprivation (N = 11)")),
    limits = c("Sleep", "Wake"),
    values = c("Sleep" = "#595959", "Wake" = "#d7dbdd")
  ) +
  geom_hline(colour = "black",
             yintercept = 0,
             size = .1)


DR_Cum_plot_hab <- DR_Cum_plot_hab +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = "none")          
DR_Cum_plot_hab

# Save figure 
ggsave(
  file = "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/04-Figures/interference_DR.svg",
  plot = DR_Cum_plot_hab,
  width = 37,
  height = 28,
  units = "mm"
)


# 07 - Control Parameters test - control experiment --------------------------
# Total exploration time 
test_expl_time <- subset(test_wide, test_wide$Type == "Interference")

test_expl_time <-
  subset(
    test_expl_time,
    select = c(
      "Animal",
      "Type",
      "Retention",
      "Total_exp_time"
    )
  )


test_expl_time_sum = describeBy(
  test_expl_time$Total_exp_time,
  list(test_expl_time$Retention),
  mat = TRUE,
  digits = 2
)

# Comparisons

stat.test1 <- test_expl_time %>%
  t_test(Total_exp_time ~ Retention) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('****', '***', '**', '*' , 'ns'))  %>%
  add_xy_position(x = "Retention",
                  step.increase = 0)%>%
  mutate(y.position = y.position + 8)

test_expl_time_plot <-
  ggplot(data = test_expl_time_sum, aes(x = group1, y = mean, fill = group1)) +
  geom_bar(
    stat = 'identity',
    position = dodge,
    width = .8,
    colour = "black",
    linewidth = .2
  ) +
  geom_errorbar(limits, position = dodge, width = 0.3,
                linewidth = .2) +
  geom_point(
    data = test_expl_time,
    aes(x = Retention, y = Total_exp_time, fill = Retention),
    shape = 21,  
    size = 0.8,
    alpha = 0.6
  ) +
  stat_pvalue_manual(
    stat.test1,
    label = "{p.signif}",
    tip.length = 0.01,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  ) +
  theme_classic() +
  scale_y_continuous(name = "Exploration time (s)") +
  scale_x_discrete(
    name = NULL, 
    labels = NULL,
    limits = c("Sleep", "Wake")
  ) +
  scale_fill_manual(
    values = c("Sleep" = "#595959", "Wake" = "#d7dbdd")
  ) 


test_expl_time_plot <- test_expl_time_plot +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = "none")          

test_expl_time_plot

# Disctance traveled 

test_dist <- subset(test_wide, test_wide$Type == "Interference")

test_dist <-
  subset(
    test_dist,
    select = c(
      "Animal",
      "Type",
      "Retention",
      "Cum_dist_min_5"
    )
  )


test_dist_sum = describeBy(
  test_dist$Cum_dist_min_5,
  list(test_dist$Retention),
  mat = TRUE,
  digits = 2
)

# Comparisons

stat.test1 <- test_dist %>%
  t_test(Cum_dist_min_5 ~ Retention) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('****', '***', '**', '*' , 'ns'))  %>%
  add_xy_position(x = "Retention",
                  step.increase = 0)%>%
  mutate(y.position = y.position + 8)

test_dist_plot <-
  ggplot(data = test_dist_sum, aes(x = group1, y = mean, fill = group1)) +
  geom_bar(
    stat = 'identity',
    position = dodge,
    width = .8,
    colour = "black",
    linewidth = .2
  ) +
  geom_errorbar(limits, position = dodge, width = 0.3,
                linewidth = .2) +
  geom_point(
    data = test_dist,
    aes(x = Retention, y = Cum_dist_min_5, fill = Retention),
    shape = 21,  
    size = 0.8,
    alpha = 0.6
  ) +
  stat_pvalue_manual(
    stat.test1,
    label = "{p.signif}",
    tip.length = 0.01,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  ) +
  theme_classic() +
  scale_y_continuous(name = "Distance (m)") +
  scale_x_discrete(
    name = NULL, 
    labels = NULL,
    limits = c("Sleep", "Wake")
  ) +
  scale_fill_manual(
    values = c("Sleep" = "#595959", "Wake" = "#d7dbdd")
  ) 


test_dist_plot<- test_dist_plot +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = "none")          

test_dist_plot

# save figure 
ggsave(
  file = "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/04-Figures/test_ctrl_expl.svg",
  plot = test_expl_time_plot,
  width = 16,
  height = 28,
  units = "mm"
)

ggsave(
  file = "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/04-Figures/test_ctrl_dist.svg",
  plot = test_dist_plot,
  width = 14,
  height = 28,
  units = "mm"
)

# 08 - Encoding comparison across experiments -----------------------------------
enc_sleep <- subset(enc_wide, enc_wide$Retention == "Sleep")

Cum_enc_sleep <-
  subset(
    enc_sleep,
    select = c(
      "Animal",
      "Type",
      "Retention",
      "Sampling",
      "Total_exp_time"
    )
  )

Cum_enc_sleep$Retention = as.factor(as.character(Cum_enc_sleep$Retention))

Cum_enc_sleep_sum = describeBy(
  Cum_enc_sleep$Total_exp_time,
  list(Cum_enc_sleep$Type, Cum_enc_sleep$Retention, Cum_enc_sleep$Sampling),
  mat = TRUE,
  digits = 2
)

# Overall statistics
basic_all_hlm <-
  lmer(Total_exp_time ~ Sampling + Type + (1 | Animal),
       data = Cum_enc_sleep,
       REML = FALSE)

basic_type_hlm <-
  lmer(Total_exp_time ~ Sampling + (1 | Animal),
       data = Cum_enc_sleep,
       REML = FALSE)

basic_sampling_hlm <-
  lmer(Total_exp_time ~ Type + (1 | Animal),
       data = Cum_enc_sleep,
       REML = FALSE)

interaction_all_hlm <-
  lmer(Total_exp_time ~ Sampling  * Type + (1 | Animal),
       data = Cum_enc_sleep,
       REML = FALSE)


anova(basic_type_hlm, basic_all_hlm) # test main effect experiment type
anova(basic_sampling_hlm, basic_all_hlm) # test main effect min
anova(interaction_all_hlm, basic_all_hlm) # test interaction effect with minute


# Post-hoc Comparisons
stat.test2 <- Cum_enc_sleep %>%
  group_by(Sampling) %>%
  t_test(Total_exp_time ~ Type) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('****', '***', '**', '*' , 'ns'))  %>%
  add_xy_position(x = "Sampling",
                  dodge = 0.8,
                  step.increase = 0)%>%
  mutate(y.position = 0.5)

enc_sleep_exp_plot<-
  ggplot(data = Cum_enc_sleep_sum, aes(x = group3, y = mean, fill = group1)) +
  geom_bar(
    stat = 'identity',
    position = dodge,
    width = .8,
    colour = "black",
    linewidth = .2
  ) +
  geom_errorbar(limits, position = dodge, width = 0.3,
                linewidth = .2) +
  geom_jitter(
    data = Cum_enc_sleep,
    aes(x = Sampling, y = Total_exp_time, fill = Type),
    shape = 21,  # Hollow point with fill
    size = 0.8,
    alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)
  ) +
  stat_pvalue_manual(
    stat.test2,
    label = "{p.signif}",
    x = "xmax",
    remove.bracket = TRUE,
    hide.ns = TRUE,
    size = 2) +
  theme_classic() +
  scale_y_continuous(name = "Exploration time (s)") +
  scale_x_discrete(
    name = "Sleep Groups - Encoding Episode"
  ) +
  scale_fill_manual(
    name = "Habituation",
    labels = c("Empty","Objects"),
    limits = c("Schema", "Interference"),
    values = c("Schema" = "#D40202", "Interference" = "#FFFFFF")
  ) 

enc_sleep_exp_plot <- enc_sleep_exp_plot +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = "none")          
enc_sleep_exp_plot


# save figure 
ggsave(
  file = "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/04-Figures/encoding_comp.svg",
  plot = enc_sleep_exp_plot,
  width = 45,
  height = 34,
  units = "mm"
)


# 09 - Comparison of post-encoding sleep ----------------------------------
sleep$SleepDuration <- as.numeric(as.character(sleep$SleepDuration))
sleep <- reorder_levels(sleep, "Type", order = c("Schema", "Interference"))


Sleep_sum = describeBy(
  sleep$SleepDuration,
  list(sleep$Type),
  mat = TRUE,
  digits = 2
)

# Comparisons

stat.test1 <- sleep %>%
  t_test(SleepDuration ~ Type) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('****', '***', '**', '*' , 'ns'))  %>%
  add_xy_position(x = "Type",
                  step.increase = 0)%>%
  mutate(y.position = y.position + 8)

sleep_plot <-
  ggplot(data = Sleep_sum, aes(x = group1, y = mean, fill = group1)) +
  geom_bar(
    stat = 'identity',
    position = dodge,
    width = .8,
    colour = "black",
    linewidth = .2
  ) +
  geom_errorbar(limits, position = dodge, width = 0.3,
                linewidth = .2) +
  geom_point(
    data = sleep,
    aes(x = Type, y = SleepDuration, fill = Type),
    shape = 21,  
    size = 0.8,
    alpha = 0.6
  ) +
  stat_pvalue_manual(
    stat.test1,
    label = "{p.signif}",
    tip.length = 0.01,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  ) +
  theme_classic() +
  scale_y_continuous(name = "Post-encoding sleep (min)") +
  scale_x_discrete(
    name = NULL, 
    labels = NULL,
    limits = c("Schema", "Interference")
  ) +
  scale_fill_manual(
    values = c("Interference" = "#D40202", "Schema" = "#ffffff")
  ) 


sleep_plot <- sleep_plot +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = "none")          

sleep_plot

# save figure 
ggsave(
  file = "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/04-Figures/sleep_comp.svg",
  plot = sleep_plot,
  width = 17,
  height = 30,
  units = "mm"
)

# 10 - Total exploration time comparison during encoding ------------------
Cum_enc_exp <-
  subset(
    enc_wide,
    select = c(
      "Animal",
      "Type",
      "Retention",
      "Sampling",
      "Total_exp_time"
    )
  )

Cum_enc_exp <- subset(Cum_enc_exp, Cum_enc_exp$Type == "Schema")

Cum_enc_exp_sum = describeBy(
  Cum_enc_exp$Total_exp_time,
  list(Cum_enc_exp$Retention, Cum_enc_exp$Sampling),
  mat = TRUE,
  digits = 2
)

# Comparisons
# Overall statistics
basic_all_hlm <-
  lmer(Total_exp_time ~ Sampling + Retention + (1 | Animal),
       data = Cum_enc_exp,
       REML = FALSE)

basic_ret_hlm <-
  lmer(Total_exp_time ~ Sampling + (1 | Animal),
       data = Cum_enc_exp,
       REML = FALSE)

basic_min_hlm <-
  lmer(Total_exp_time ~ Retention + (1 | Animal),
       data = Cum_enc_exp,
       REML = FALSE)

interaction_all_hlm <-
  lmer(Total_exp_time ~ Sampling  * Retention + (1 | Animal),
       data = Cum_enc_exp,
       REML = FALSE)


anova(basic_ret_hlm, basic_all_hlm) # test main effect retention
anova(basic_min_hlm, basic_all_hlm) # test main effect min
anova(interaction_all_hlm, basic_all_hlm) # test interaction effect with minute

# Mean exploration time
Cum_enc_exp__total_sum = describeBy(
  Cum_enc_exp$Total_exp_time,
  list(Cum_enc_exp$Retention),
  mat = TRUE,
  digits = 2
)

Cum_enc_exp__total_sum


# Plotting
stat.test1 <- Cum_enc_exp %>%
  group_by(Sampling) %>%
  t_test(Total_exp_time ~ Retention) %>%
  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols =  c('****', '***', '**', '*' , 'ns'))  %>%
  add_xy_position(x = "Sampling",
                  step.increase = 0)%>%
  mutate(y.position = y.position + 4)

Cum_enc_exp_plot <-
  ggplot(data = Cum_enc_exp_sum, aes(x = group2, y = mean, fill = group1)) +
  geom_bar(
    stat = 'identity',
    position = dodge,
    width = .8,
    colour = "black"
  ) +
  geom_errorbar(limits, position = dodge, width = 0.3) +
  geom_jitter(
    data = Cum_enc_exp,
    aes(x = Sampling, y = Total_exp_time, fill = Retention),
    shape = 21,  
    size = 1.2,
    alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)
  ) +
  stat_pvalue_manual(
    stat.test1,
    label = "{p.signif}",
    tip.length = 0.01,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  ) +
  theme_classic() +
  scale_y_continuous(name = "Total exploration time (s)") +
  scale_x_discrete(
    name = "Encoding Episode"
  ) +
  scale_fill_manual(
    name = "Post-encoding condition",
    labels = c("Sleep","Sleep deprivation (SD)"),
    limits = c("Sleep", "Wake"),
    values = c("Sleep" = "#595959", "Wake" = "#d7dbdd")
  ) 

Cum_enc_exp_plot <- Cum_enc_exp_plot +   theme_classic(base_size = 8) +  # set global font and size
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    strip.text = element_text(size = 8),
    legend.position = "none")          

Cum_enc_exp_plot
