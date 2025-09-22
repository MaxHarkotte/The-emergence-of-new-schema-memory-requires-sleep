# Harkotte et al (2025) - Spatial schema memory formation in rats is dependent on post-encoding sleep 
# Analyses of cfos counts across brain regions, including correlational network analyses
# last modified: August 2025
# maximilian.harkotte@gmail.com

rm(list = ls()) # clear workspace
cat("\014") # clear console

# 00 - Load packages -------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(car)
library(performance)
library(rstatix)
library(igraph)
library(scales)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(lme4)

# 01 - Source file ---------------------------------------------------------
dataPath <- "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/05-Histology"
# dataPath <- "/Volumes/born_animal/Max/02_Schema_memory_formation/02_Schema_during_sleep/05-Histology"
setwd(dataPath)

# 02 - Read in data --------------------------------------------------------
cfos <- read.csv2("cfos.csv", header = TRUE, sep = ",")

# 03 - Data wrangling ------------------------------------------------------
cfos_long <- cfos %>%
  pivot_longer(
    cols = -c(ID, condition),
    names_to = "BrainRegion", 
    values_to = "CellCount"
  )

# Compute mean cell count per brain region for home cage control
baseline_means <- cfos %>%
  filter(condition == "home_cage") %>%
  summarise(across(-c(ID, condition), ~ mean(.x, na.rm = TRUE)))

# Normalize by home cage
adjusted_cfos <- cfos %>%
  filter(condition != "home_cage") %>%
  mutate(across(-c(ID, condition), ~ . / baseline_means[[cur_column()]]))

adjusted_cfos[ ,c('nucleus_reuniens', 'reticular_nucleus')] <- list(NULL) # ignore thalamic data

adjusted_cfos_long <- adjusted_cfos %>%
  pivot_longer(
    cols = -c(ID, condition),
    names_to = "BrainRegion",
    values_to = "CellCount"
  )

# Convert factors
cfos_long <- cfos_long %>% mutate(ID = as.factor(ID), condition = as.factor(condition))
adjusted_cfos_long <- adjusted_cfos_long %>% mutate(ID = as.factor(ID), condition = as.factor(condition))

# 04 - mPFC --------------------------------------------------------------------
mPFC_adjusted <- adjusted_cfos_long %>%
  filter(BrainRegion %in% c("prelimbic", "infralimbic", "anterior_CG1", "medial_CG1", "medial_CG2")) 

mPFC_adjusted$BrainRegion = factor(mPFC_adjusted$BrainRegion, levels= c("prelimbic", "infralimbic", "anterior_CG1", "medial_CG1", "medial_CG2"))

mPFC_adjusted$BrainRegion <- factor(mPFC_adjusted$BrainRegion)

# Overall statistics
basic_all_hlm <-
  lmer(CellCount ~ condition + BrainRegion + (1 | ID),
       data = mPFC_adjusted,
       REML = FALSE)

basic_reg_hlm <-
  lmer(CellCount ~ condition + (1 | ID),
       data = mPFC_adjusted,
       REML = FALSE)

basic_cond_hlm <-
  lmer(CellCount ~ BrainRegion + (1 | ID),
       data = mPFC_adjusted,
       REML = FALSE)

interaction_all_hlm <-
  lmer(CellCount ~ condition  * BrainRegion + (1 | ID),
       data = mPFC_adjusted,
       REML = FALSE)


anova(basic_reg_hlm, basic_all_hlm) # test main effect region
anova(basic_cond_hlm, basic_all_hlm) # test main effect post-encoding condition
anova(interaction_all_hlm, basic_all_hlm) # test interaction effect 

# Post-hoc Tests
stat.test1 <- mPFC_adjusted %>%
  group_by(condition, BrainRegion) %>%
  t_test(CellCount ~ 1, mu = 1) %>%
  add_significance("p", symbols =  c('####', '###', '##', '#', 'ns')) %>%
  mutate(group1 = condition, group2 = BrainRegion) %>%
  add_xy_position(x = "BrainRegion",
                  group = "condition",
                  step.increase = 0) %>%
  mutate(y.position = y.position + 0.2, xmax = xmax -0.2)

stat.test2 <- mPFC_adjusted %>%
  group_by(BrainRegion) %>%
  t_test(CellCount ~ condition) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "BrainRegion",
                  dodge = 0.8,
                  step.increase = 0)

stat.test3 <- mPFC_adjusted %>%
  group_by(condition) %>%
  t_test(CellCount ~ BrainRegion) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "BrainRegion",
                  group = "condition",
                  step.increase = 0)

# Plotting
bxp_mPFC_adj <- ggboxplot(
  mPFC_adjusted,
  x = "BrainRegion",
  y = "CellCount",
  fill = "condition",
  width = 0.6,
  lwd = 0.3,  
  outlier.size = 0.1, 
  ylim = c(0,3.8)
) + 
  stat_pvalue_manual(
  stat.test1,
  label = "{p.signif}",
  x = "xmax",
  remove.bracket = TRUE,
  hide.ns = TRUE,
  size = 2
) +
  labs(x = "Brain Region", y = "c-Fos+ cell count (%)") +
  geom_hline(colour = "black", yintercept = 1, linetype="dashed", linewidth = 0.2) +
  stat_pvalue_manual(
    stat.test2,
    label = "{p.adj.signif}",
    tip.length = 0.01,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  ) +
  stat_pvalue_manual(
    stat.test3,
    label = "{p.adj.signif}",
    tip.length = 0.01,
    step.increase = 0.07,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  )+
  scale_x_discrete(
    name = "Medial Prefrontal Cortex",
    labels = c("PL", "IL", "aCG1", "mCG1", "mCG2"),
    limits = c(
      "prelimbic",
      "infralimbic",
      "anterior_CG1",
      "medial_CG1",
      "medial_CG2"
      
    )
  )+
  scale_fill_manual(
    name = "Retention condition",
    labels = c("Sleep (N = 10)", "Wake (N = 10)"),
    limits = c("sleep", "sleep_deprivation"), 
    values = c("sleep" = "#595959", "sleep_deprivation" = "#d7dbdd")
  )  

bxp_mPFC_adj  <- bxp_mPFC_adj +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7)
  )
bxp_mPFC_adj 

# 05 - Hippocampus -------------------------------------------------------------
hippocampus_adjusted <- adjusted_cfos_long %>%
  filter(grepl("^CA", BrainRegion) | BrainRegion == "dentate")

hippocampus_adjusted$BrainRegion <- factor(
  hippocampus_adjusted$BrainRegion)

# Overall statistics
basic_all_hlm <-
  lmer(CellCount ~ condition + BrainRegion + (1 | ID),
       data = hippocampus_adjusted,
       REML = FALSE)

basic_reg_hlm <-
  lmer(CellCount ~ condition + (1 | ID),
       data = hippocampus_adjusted,
       REML = FALSE)

basic_cond_hlm <-
  lmer(CellCount ~ BrainRegion + (1 | ID),
       data = hippocampus_adjusted,
       REML = FALSE)

interaction_all_hlm <-
  lmer(CellCount ~ condition  * BrainRegion + (1 | ID),
       data = hippocampus_adjusted,
       REML = FALSE)


anova(basic_reg_hlm, basic_all_hlm) # test main effect region
anova(basic_cond_hlm, basic_all_hlm) # test main effect post-encoding condition
anova(interaction_all_hlm, basic_all_hlm) # test interaction effect 

# Post-hoc Tests
stat.test1 <- hippocampus_adjusted %>%
  group_by(condition, BrainRegion) %>%
  t_test(CellCount ~ 1, mu = 1) %>%
  add_significance("p", symbols =  c('####', '###', '##', '#', 'ns')) %>%
  mutate(group1 = condition, group2 = BrainRegion) %>%
  add_xy_position(x = "BrainRegion",
                  group = "condition",
                  step.increase = 0) %>%
  mutate(y.position = y.position + 0.2)

stat.test2 <- hippocampus_adjusted %>%
  group_by(BrainRegion) %>%
  t_test(CellCount ~ condition) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "BrainRegion",
                  dodge = 0.8,
                  step.increase = 0)

stat.test3 <- hippocampus_adjusted %>%
  group_by(condition) %>%
  t_test(CellCount ~ BrainRegion) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "BrainRegion",
                  group = "condition",
                  step.increase = 0)

# Plotting
bxp_hippocampus_adj <- ggboxplot(
  hippocampus_adjusted,
  x = "BrainRegion",
  y = "CellCount",
  fill = "condition",
  width = 0.6,
  lwd = 0.3,  
  outlier.size = 0.1, 
  ylim = c(0,3.8)
) + stat_pvalue_manual(
  stat.test1,
  label = "{p.signif}",
  x = "xmax",
  remove.bracket = TRUE,
  hide.ns = TRUE,
  size = 2
) +
  labs(x = "Brain Region", y = "c-Fos+ cell count (%)") +
  theme_classic(base_size = 18) +
  geom_hline(colour = "black", yintercept = 1, linetype="dashed", linewidth = 0.2) +
  stat_pvalue_manual(
    stat.test2,
    label = "{p.adj.signif}",
    tip.length = 0.01,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  ) +
  stat_pvalue_manual(
    stat.test3,
    label = "{p.adj.signif}",
    tip.length = 0.01,
    step.increase = 0.07,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  )+
  scale_x_discrete(
    name = "Hippocampus",
    labels = c("CA1dist", "CA1prox", "CA3", "DG"),
    limits = c(
      "CA1dist",
      "CA1prox",
      "CA3",
      "dentate"
    )
  )+
  scale_fill_manual(
    name = "Retention condition",
    labels = c("Sleep (N = 10)", "Wake (N = 10)"),
    limits = c("sleep", "sleep_deprivation"), 
    values = c("sleep" = "#595959", "sleep_deprivation" = "#d7dbdd")
  )     

bxp_hippocampus_adj  <- bxp_hippocampus_adj +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7)
  )

bxp_hippocampus_adj  

mPFC_HC <- ggarrange(
  bxp_mPFC_adj,
  bxp_hippocampus_adj,
  widths = c(1.1, 0.9),
  ncol = 2,
  nrow = 1, 
  common.legend = TRUE
)

mPFC_HC

# 06 - Cortex --------------------------------------------------------
cortex_adjusted <- adjusted_cfos_long %>%
  filter(BrainRegion %in% c("parietal", "perirhinal", "entorhinal"))

cortex_adjusted$BrainRegion <- factor(
  cortex_adjusted$BrainRegion)

cortex_adjusted$BrainRegion = factor(cortex_adjusted$BrainRegion, levels= c("parietal", "entorhinal", "perirhinal"))

# ANOVA
# Overall statistics
basic_all_hlm <-
  lmer(CellCount ~ condition + BrainRegion + (1 | ID),
       data = cortex_adjusted,
       REML = FALSE)

basic_reg_hlm <-
  lmer(CellCount ~ condition + (1 | ID),
       data = cortex_adjusted,
       REML = FALSE)

basic_cond_hlm <-
  lmer(CellCount ~ BrainRegion + (1 | ID),
       data = cortex_adjusted,
       REML = FALSE)

interaction_all_hlm <-
  lmer(CellCount ~ condition  * BrainRegion + (1 | ID),
       data = cortex_adjusted,
       REML = FALSE)


anova(basic_reg_hlm, basic_all_hlm) # test main effect region
anova(basic_cond_hlm, basic_all_hlm) # test main effect post-encoding condition
anova(interaction_all_hlm, basic_all_hlm) # test interaction effect 

# Post-hoc Tests
stat.test1 <- cortex_adjusted %>%
  group_by(condition, BrainRegion) %>%
  t_test(CellCount ~ 1, mu = 1) %>%
  add_significance("p", symbols =  c('####', '###', '##', '#', 'ns')) %>%
  mutate(group1 = condition, group2 = BrainRegion) %>%
  add_xy_position(x = "BrainRegion",
                  group = "condition",
                  step.increase = 0) %>%
  mutate(y.position = y.position + 0.2)

stat.test2 <- cortex_adjusted %>%
  group_by(BrainRegion) %>%
  t_test(CellCount ~ condition) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "BrainRegion",
                  dodge = 0.8,
                  step.increase = 0)

stat.test3 <- cortex_adjusted %>%
  group_by(condition) %>%
  t_test(CellCount ~ BrainRegion) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "BrainRegion",
                  group = "condition",
                  step.increase = 0)


# Plotting
bxp_cortex_adj <- ggboxplot(
  cortex_adjusted,
  x = "BrainRegion",
  y = "CellCount",
  fill = "condition",
  width = 0.6,
  lwd = 0.3,  
  outlier.size = 0.1, 
  ylim = c(0,6.5)
) + stat_pvalue_manual(
  stat.test1,
  label = "{p.signif}",
  x = "xmax",
  remove.bracket = TRUE,
  hide.ns = TRUE,
  size = 2
) +
  labs(x = "Brain Region", y = "c-Fos+ cell count (%)") +
  geom_hline(colour = "black", yintercept = 1, linetype="dashed", linewidth = 0.2) +
  stat_pvalue_manual(
    stat.test2,
    label = "{p.adj.signif}",
    tip.length = 0.01,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  ) +
  stat_pvalue_manual(
    stat.test3,
    label = "{p.adj.signif}",
    tip.length = 0.01,
    step.increase = 0.07,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  )+
  scale_x_discrete(
    name = "Posterior Cortex",
    labels = c("PAR", "ENT", "PRH"),
    limits = c(
      "parietal",
      "entorhinal",
      "perirhinal"
    )
  )+
  scale_fill_manual(
    name = "Retention condition",
    labels = c("Sleep (N = 10)", "Wake (N = 10)"),
    limits = c("sleep", "sleep_deprivation"), 
    values = c("sleep" = "#595959", "sleep_deprivation" = "#d7dbdd"),
    guide = "none"
  )     

bxp_cortex_adj <- bxp_cortex_adj +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7)
  )

bxp_cortex_adj

# 07 - Retrosplenial cortex ----------------------------------------------------
retro_adjusted <- adjusted_cfos_long %>%
  filter(
    BrainRegion %in% c(
      "anterior_retrosplenial_agranular",
      "anterior_retrosplenial_granular",
      "medial_retrosplenial_agranular",
      "medial_retrosplenial_granular",
      "posterior_retrosplenial_agranular",
      "posterior_retrosplenial_granular"
    )
  )

retro_adjusted$BrainRegion <- factor(
  retro_adjusted$BrainRegion,
  levels = c(
    "anterior_retrosplenial_granular",
    "medial_retrosplenial_granular",
    "posterior_retrosplenial_granular",
    "anterior_retrosplenial_agranular",
    "medial_retrosplenial_agranular",
    "posterior_retrosplenial_agranular"
  )
)
# ANOVA
basic_all_hlm <-
  lmer(CellCount ~ condition + BrainRegion + (1 | ID),
       data = retro_adjusted,
       REML = FALSE)

basic_reg_hlm <-
  lmer(CellCount ~ condition + (1 | ID),
       data = retro_adjusted,
       REML = FALSE)

basic_cond_hlm <-
  lmer(CellCount ~ BrainRegion + (1 | ID),
       data = retro_adjusted,
       REML = FALSE)

interaction_all_hlm <-
  lmer(CellCount ~ condition  * BrainRegion + (1 | ID),
       data = retro_adjusted,
       REML = FALSE)


anova(basic_reg_hlm, basic_all_hlm) # test main effect region
anova(basic_cond_hlm, basic_all_hlm) # test main effect post-encoding condition
anova(interaction_all_hlm, basic_all_hlm) # test interaction effect

# Post-hoc Tests
stat.test1 <- retro_adjusted %>%
  group_by(condition, BrainRegion) %>%
  t_test(CellCount ~ 1, mu = 1) %>%
  add_significance("p", symbols =  c('####', '###', '##', '#', 'ns')) %>%
  mutate(group1 = condition, group2 = BrainRegion) %>%
  add_xy_position(x = "BrainRegion",
                  group = "condition",
                  step.increase = 0) %>%
  mutate(y.position = y.position + 0.2)

stat.test2 <- retro_adjusted %>%
  group_by(BrainRegion) %>%
  t_test(CellCount ~ condition) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "BrainRegion",
                  dodge = 0.8,
                  step.increase = 0)

stat.test3 <- retro_adjusted %>%
  group_by(condition) %>%
  t_test(CellCount ~ BrainRegion) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "BrainRegion",
                  group = "condition",
                  step.increase = 0)


# Plotting
bxp_retro_adj <- ggboxplot(
  retro_adjusted,
  x = "BrainRegion",
  y = "CellCount",
  fill = "condition",
  width = 0.6,
  lwd = 0.3,  
  outlier.size = 0.1, 
  ylim = c(0,6.5)
) + stat_pvalue_manual(
  stat.test1,
  label = "{p.signif}",
  x = "xmax",
  remove.bracket = TRUE,
  hide.ns = TRUE,
  size = 2
) +
  labs(x = "Brain Region", y = "c-Fos+ cell count (%)") +
  geom_hline(colour = "black", yintercept = 1, linetype="dashed", linewidth = 0.2) +
  stat_pvalue_manual(
    stat.test2,
    label = "{p.adj.signif}",
    tip.length = 0.01,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  ) +
  stat_pvalue_manual(
    stat.test3,
    label = "{p.adj.signif}",
    tip.length = 0.01,
    step.increase = 0.07,
    bracket.nudge.y = 0.8,
    hide.ns = TRUE,
    size = 2.5
  )+
  scale_x_discrete(
    name = "Cortex - Retrosplenial granular/agranular",
    labels = c("aRSG", "mRSG", "pRSG", "aRSA", "mRSA", "pRSA"),
    limits = c(
      "anterior_retrosplenial_granular",
      "medial_retrosplenial_granular",
      "posterior_retrosplenial_granular",
      "anterior_retrosplenial_agranular",
      "medial_retrosplenial_agranular",
      "posterior_retrosplenial_agranular"
    )
  )+
  scale_fill_manual(
    name = "Retention condition",
    labels = c("Sleep (N = 10)", "Wake (N = 10)"),
    limits = c("sleep", "sleep_deprivation"), 
    values = c("sleep" = "#595959", "sleep_deprivation" = "#d7dbdd"),
    guide = "none"
  )     

bxp_retro_adj <- bxp_retro_adj +   theme_classic(base_size = 7) +  # set global font and size
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7)
  )

bxp_retro_adj 

cortex <- ggarrange(
  bxp_retro_adj,
  bxp_cortex_adj,
  widths = c(1.3, 0.7),
  ncol = 2,
  nrow = 1
)

cortex 

Fig3_lower_panel <- ggarrange(
  mPFC_HC,
  cortex,
  ncol = 1,
  nrow = 2
)

Fig3_lower_panel

# Save plot 
ggsave(
  file = "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/04-Figures/Fig3_lower_panel.svg",
  plot = Fig3_lower_panel,
  width = 165,
  height = 110,
  units = "mm"
)

# 08 - Correlational network analysis ----------------------------------------

# Threshold 
n <- 10
df <- n-2
alpha<-0.05

t_critical <- qt(1 - alpha/2, df)

r_critical <- t_critical / sqrt(t_critical^2 + df)
print(r_critical)


categories <- c(
  rep("mPFC", 5), rep("Cortex", 9), 
  rep("Hippocampus", 4)
)

colors <- c("mPFC" = "#AE282C", "Cortex" = "#F0C571", "Hippocampus" = "#2066a8")

# SLEEP DEPRIVATION 
cfos_wake <- subset(adjusted_cfos, adjusted_cfos$condition == "sleep_deprivation")
subset_data <- cfos_wake[3:20]

# Brain Regions and Colors (cfos-icv)
brain_regions <- colnames(subset_data)

# Correlation matrix
correlation_matrix <- cor(subset_data, method = "pearson", use = "pairwise.complete.obs")

# Substitute values in the main diagonal with NA
correlation_matrix[lower.tri(correlation_matrix)] <- NA

diag(correlation_matrix) <- NA

# Extract only significant correlations based on the threshold
significant_correlations <- sum(abs(correlation_matrix) > r_critical, na.rm = TRUE)

# Binary matrix for significant connection
binary_matrix <- abs(correlation_matrix) >= r_critical
edges <- which(binary_matrix, arr.ind = TRUE)

# Dataframe of connections
edges_df <- data.frame(
  from = rownames(correlation_matrix)[edges[, 1]],
  to = colnames(correlation_matrix)[edges[, 2]],
  weight = correlation_matrix[edges]
)

# Graph with all the nodes, even the ones without sign.correlations
brain_graph <- graph_from_data_frame(edges_df, directed = FALSE, vertices = data.frame(name = brain_regions))

# Convert to tidygraph
graph_tbl <- as_tbl_graph(brain_graph)

# Add x/y layout coordinates
layout <- layout_in_circle(brain_graph)
layout_df <- as.data.frame(layout)
colnames(layout_df) <- c("x", "y")

mean_activity <- colMeans(subset_data, na.rm = TRUE)  # named vector: names = region names

# Also pass category, color, and mean activity explicitly
vertex_df <- data.frame(
  name = brain_regions,
  category = categories,
  activity = mean_activity,
  x = layout_df$x,
  y = layout_df$y
)

vertex_df <- vertex_df %>%
  mutate(
    color = colors[category]
  )

vertex_df$category <- factor(vertex_df$category, levels = names(colors))

# Join vertex attributes to graph_tbl
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(vertex_df, by = "name")

# Abbreviated region labels
labels_abbr <- c("PL", "IL", "aCG1", "mCG1", "mCG2", 
                 "aRSA", "aRSG", "mRSA", "mRSG", "pRSA", "pRSG", 
                 "PRH", "ENT", "PAR", "CA1dist", "CA1prox", "CA3", "DG")

graph_tbl <- graph_tbl %>% 
  mutate(label = labels_abbr)

# Plot with legends
wake_network <- ggraph(graph_tbl, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(width = abs(weight)), color = "black", show.legend = FALSE) +
  scale_edge_width(range = c(0.2, 1.2)) +
  geom_node_point(aes(size = activity, fill = category), shape = 21, color = "black") +
  geom_node_text(aes(x = 1.15 * x, y = 1.15 * y, label = label), size = 2) +
  scale_fill_manual(values = colors, name = "Brain area",
                    guide = guide_legend(override.aes = list(size = 6, shape = 21, fill = colors, stroke = 0.5))) +
  scale_size_continuous(limits = c(0, 3),breaks = c(0., 1, 2), range = c(1, 6), name = "Relative cFos activity [%]") +
  theme_void() +
  ggtitle("Wake")+
  guides(size = "none") +
  theme( text = element_text(size = 7),
         legend.text = element_text(size = 7),
         legend.title = element_text(size = 7),
         plot.title = element_text(hjust = 0.5, face = "bold"))

wake_network 

net_wake  = brain_graph

# SLEEP
cfos_sleep <- subset(adjusted_cfos, adjusted_cfos$condition == "sleep")
subset_data <- cfos_sleep[3:20]

# Correlation matrix
correlation_matrix <- cor(subset_data, method = "pearson", use = "pairwise.complete.obs")

# Substitute values in the main diagonal with NA
correlation_matrix[lower.tri(correlation_matrix)] <- NA

diag(correlation_matrix) <- NA

# Extract only significant correlations based on the threshold
significant_correlations <- sum(abs(correlation_matrix) > r_critical, na.rm = TRUE)

# Binary matrix for significant connection
binary_matrix <- abs(correlation_matrix) >= r_critical
edges <- which(binary_matrix, arr.ind = TRUE)

# Dataframe of connections
edges_df <- data.frame(
  from = rownames(correlation_matrix)[edges[, 1]],
  to = colnames(correlation_matrix)[edges[, 2]],
  weight = correlation_matrix[edges]
)

# Graph with all the nodes, even the ones without sign.correlations
brain_graph <- graph_from_data_frame(edges_df, directed = FALSE, vertices = data.frame(name = brain_regions))

# Convert to tidygraph
graph_tbl <- as_tbl_graph(brain_graph)

# Add x/y layout coordinates
layout <- layout_in_circle(brain_graph)
layout_df <- as.data.frame(layout)
colnames(layout_df) <- c("x", "y")

mean_activity <- colMeans(subset_data, na.rm = TRUE)  # named vector: names = region names

# Also pass category, color, and mean activity explicitly
vertex_df <- data.frame(
  name = brain_regions,
  category = categories,
  activity = mean_activity,
  x = layout_df$x,
  y = layout_df$y
)

vertex_df <- vertex_df %>%
  mutate(
    color = colors[category]
  )

vertex_df$category <- factor(vertex_df$category, levels = names(colors))

# Join vertex attributes to graph_tbl
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(vertex_df, by = "name")

# Abbreviated region labels
labels_abbr <- c("PL", "IL", "aCG1", "mCG1", "mCG2", 
                 "aRSA", "aRSG", "mRSA", "mRSG", "pRSA", "pRSG", 
                 "PRH", "ENT", "PAR", "CA1dist", "CA1prox", "CA3", "DG")

graph_tbl <- graph_tbl %>% 
  mutate(label = labels_abbr)

# Plot with legends
sleep_network <- ggraph(graph_tbl, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(width = abs(weight)), color = "black", show.legend = FALSE) +
  scale_edge_width(range = c(0.2, 1.2)) +
  geom_node_point(aes(size = activity, fill = category), shape = 21, color = "black") +
  geom_node_text(aes(x = 1.15 * x, y = 1.15 * y, label = label), size = 2) +
  scale_fill_manual(values = colors, name = "Brain area",
                    guide = guide_legend(override.aes = list(size = 6, shape = 21, fill = colors, stroke = 0.5))) +
  scale_size_continuous(limits = c(0, 3),breaks = c(0.1,  1, 2),range = c(1, 6), name = "Relative cFos activity [%]") +
  theme_void() +
  ggtitle("Sleep")+
  guides(
    size = "none"
  )+
  theme( text = element_text(size = 7),
         legend.text = element_text(size = 7),
         legend.title = element_text(size = 7),
         plot.title = element_text(hjust = 0.5, face = "bold"))

sleep_network

net_sleep  = brain_graph

networks <- ggarrange(
  sleep_network,
  wake_network,
  ncol = 2,
  nrow = 1, 
  common.legend = TRUE,
  legend = "right"
)

networks

# Save plot
ggsave(
  file = "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/04-Figures/Fig3_networks.svg",
  plot = networks,
  width = 160,
  height = 70,
  units = "mm"
)


# 09 - Network statistics ------------------------------------------------------

#observed difference in edge density
obs_edge_diff = edge_density(net_sleep)-edge_density(net_wake)

#get graph nodes
nodes = unique(c(V(net_sleep)$name,V(net_wake)$name)) #20 nodes

#combine edges 
edges_sleep = as_data_frame(net_sleep, what = "edges")[,1:2] #86 edges
edges_wake  = as_data_frame(net_wake, what = "edges")[, 1:2] #24 edges
all_edges = rbind(edges_sleep,edges_wake)

#Permutatio test
set.seed(43)
n_perm     = 1000
perm_edge_diffs = numeric(n_perm)

for(i in 1:n_perm){
  #randomly assign edges to graphs 
  assign_to_s <- sample(c(TRUE,FALSE), size = nrow(all_edges), replace = TRUE)
  s_edges = all_edges[assign_to_s, ]
  w_edges = all_edges[!assign_to_s, ]
  
  #rebuild graphs 
  s_perm = graph_from_data_frame(s_edges, directed = FALSE, vertices = nodes)
  w_perm = graph_from_data_frame(w_edges, directed = FALSE, vertices = nodes)
  
  #compute permutation difference
  perm_edge_diffs[i] = edge_density(s_perm) - edge_density(w_perm)
}

#compute p-value
p_val = mean(abs(perm_edge_diffs >= obs_edge_diff))
cat('Observed edge density difference:', round(obs_edge_diff,3),'\n')
cat('Permutation p-value: < 0.001') 


#plot the results 

perm_df <- data.frame(EdgeDiff = perm_edge_diffs)

perm_res <- ggplot(perm_df, aes(x = EdgeDiff)) +
  geom_histogram(bins = 15, fill = "gray80", color = "black") +
  geom_vline(xintercept = obs_edge_diff, color = "red", size = 0.3) +  # observed difference
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", size = 0.3) +  # permutation mean
  labs(
    x = "Edge Density Difference",
    y = "Count"
  ) +
  theme_classic(base_size = 6)+     # set global font and size
  theme(
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    strip.text = element_text(size = 6)
  )

perm_res 

# Save plot 
# ggsave(
#   file = "/Volumes/born_animal/Max/02_Schema_memory_formation/02_Schema_during_sleep//04-Figures/perm_res.svg",
#   plot = perm_res ,
#   width = 55,
#   height = 38,
#   units = "mm"
# )


# 10 - Correlational network analysis with DR -----------------------------

# install.packages("corrplot")  # if needed
library(corrplot)

dataPath <- "Z:/Max/02_Schema_memory_formation/02_Schema_during_sleep/" 
# dataPath <- "/Volumes/born_animal/Max/02_Schema_memory_formation/02_Schema_during_sleep/"

setwd(dataPath)

test_wide <-
  read.csv2(
    "02-VideoScorings/02-Tables/01-DataSummary/00-TestClean.csv",
    header = TRUE,
    sep = ";",
    stringsAsFactors = TRUE
  )

test_wide <- subset(test_wide, test_wide$Type == "Schema")
test_wide <-   subset(
  test_wide,
  select = c(
    "Animal",
    "Retention",
    "Cum_DiRa_min_1",
    "Cum_DiRa_min_2",
    "Cum_DiRa_min_3",
    "Cum_DiRa_min_4",
    "Cum_DiRa_min_5"
  )
)

names(test_wide)[names(test_wide) == "Animal"] <- "ID"

# Threshold 
n <- 10
df <- n-2
alpha<-0.05

t_critical <- qt(1 - alpha/2, df)

r_critical <- t_critical / sqrt(t_critical^2 + df)
print(r_critical)


categories <- c(
  rep("mPFC", 5), rep("Cortex", 9), 
  rep("Hippocampus", 4), 
  "Memory"
)

colors <- c("mPFC" = "#AE282C", "Cortex" = "#F0C571", "Hippocampus" = "#2066a8", "Memory" = "black")

# SLEEP DEPRIVATION 
cfos_wake <- subset(adjusted_cfos, adjusted_cfos$condition == "sleep_deprivation")
subset_data <- merge(cfos_wake, test_wide, by = "ID")

subset_data <-   subset(
  subset_data,
  select = c(
    "prelimbic",
    "infralimbic",
    "anterior_CG1",
    "medial_CG1",
    "medial_CG2",
    "anterior_retrosplenial_agranular",
    "anterior_retrosplenial_granular",
    "medial_retrosplenial_agranular",
    "medial_retrosplenial_granular",
    "posterior_retrosplenial_agranular",
    "posterior_retrosplenial_granular",
    "perirhinal",
    "entorhinal",
    "parietal",
    "CA1dist",
    "CA1prox",
    "CA3",
    "dentate",
    "Cum_DiRa_min_2"
  )
)

names(subset_data)[names(subset_data) == "Cum_DiRa_min_2"] <- "Memory"

# Brain Regions and Colors (cfos-icv)
brain_regions <- colnames(subset_data)

# Correlation matrix
correlation_matrix <- cor(subset_data, method = "pearson", use = "pairwise.complete.obs")

# Substitute values in the main diagonal with NA
correlation_matrix[lower.tri(correlation_matrix)] <- NA

diag(correlation_matrix) <- NA



corrplot(correlation_matrix, method = "color", col = colorRampPalette(c("blue","white","red"))(200),
         tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")



# SLEEP
cfos_sleep <- subset(adjusted_cfos, adjusted_cfos$condition == "sleep")

subset_data <- merge(cfos_sleep, test_wide, by = "ID")

subset_data <-   subset(
  subset_data,
  select = c(
    "prelimbic",
    "infralimbic",
    "anterior_CG1",
    "medial_CG1",
    "medial_CG2",
    "anterior_retrosplenial_agranular",
    "anterior_retrosplenial_granular",
    "medial_retrosplenial_agranular",
    "medial_retrosplenial_granular",
    "posterior_retrosplenial_agranular",
    "posterior_retrosplenial_granular",
    "perirhinal",
    "entorhinal",
    "parietal",
    "CA1dist",
    "CA1prox",
    "CA3",
    "dentate",
    "Cum_DiRa_min_2"
  )
)

names(subset_data)[names(subset_data) == "Cum_DiRa_min_2"] <- "Memory"
brain_regions <- colnames(subset_data)


# Correlation matrix
correlation_matrix <- cor(subset_data, method = "pearson", use = "pairwise.complete.obs")



# Substitute values in the main diagonal with NA
correlation_matrix[lower.tri(correlation_matrix)] <- NA

diag(correlation_matrix) <- NA


corrplot(correlation_matrix, method = "color", col = colorRampPalette(c("blue","white","red"))(200),
         tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")


# 11 - cFos ~ Memory performance ------------------------------------------
all_data <- merge(adjusted_cfos, test_wide, by = "ID")

sleep_data <- subset(all_data, all_data$condition == "sleep")

sleep_data$mean_network <- rowMeans(sleep_data[, c("prelimbic", 
                                                   "infralimbic",
                                                   "anterior_CG1",
                                                   "medial_CG1",
                                                   "medial_CG2",
                                                   "anterior_retrosplenial_agranular",
                                                   "anterior_retrosplenial_granular",
                                                   "medial_retrosplenial_agranular",
                                                   "medial_retrosplenial_granular",
                                                   "posterior_retrosplenial_agranular",
                                                   "posterior_retrosplenial_granular",
                                                   "perirhinal",
                                                   "parietal",
                                                   "CA1dist",
                                                   "CA1prox",
                                                   "CA3",
                                                   "dentate")], na.rm = TRUE)

wake_data <- subset(all_data, all_data$condition == "sleep_deprivation")

wake_data$mean_network <- rowMeans(wake_data[, c("prelimbic", 
                                                   "infralimbic",
                                                   "anterior_CG1",
                                                   "medial_CG1",
                                                   "medial_CG2",
                                                   "anterior_retrosplenial_granular",
                                                   "medial_retrosplenial_agranular",
                                                   "medial_retrosplenial_granular",
                                                   "posterior_retrosplenial_agranular",
                                                   "posterior_retrosplenial_granular",
                                                   "perirhinal",
                                                   "parietal",
                                                   "entorhinal",
                                                   "CA1prox",
                                                   "CA3",
                                                   "dentate")], na.rm = TRUE)



# Correlation between last encoding DR and test DR (sleep condition)

sleep_data <- subset(sleep_data, select = c('ID','condition', 'mean_network', 'Cum_DiRa_min_2'))
wake_data <- subset(wake_data, select = c('ID','condition' , 'mean_network', 'Cum_DiRa_min_2'))

corr_data <- rbind(sleep_data, wake_data )

# Compute correlation and p-value per condition
corr_values <- corr_data %>%
  group_by(condition) %>%
  summarise(
    cor_test = list(cor.test(mean_network, Cum_DiRa_min_2, method = "spearman")),
    .groups = "drop"
  ) %>%
  mutate(
    rho = sapply(cor_test, function(x) x$estimate),
    pval = sapply(cor_test, function(x) x$p.value),
    label = paste0("rho = ", round(rho, 2), ", p = ", signif(pval, 2))
  )

# Precompute y positions for labels to avoid overlaps
# Example: stagger from max y-value
max_y <- max(corr_data$Cum_DiRa_min_2)
corr_values <- corr_values %>%
  arrange(condition) %>%
  mutate(ypos = max_y - 0.1 * (row_number() - 1) * diff(range(corr_data$Cum_DiRa_min_2)))

# Base scatter plot
DR_corr_plot <- ggscatter(
  corr_data,
  x = "mean_network",
  y = "Cum_DiRa_min_2",
  color = "condition",
  fill = "condition",
  shape = 21,
  size = 2,
  add = "reg.line",
  add.params = list(color = "condition", fill = "lightgray", size = 0.5),
  conf.int = TRUE,
  xlab = "Network activity",
  ylab = "DR at Test",
  palette = c("black", "#D40202"),
  group = "condition"
)

# Add labels using precomputed y positions
DR_corr_plot +
  geom_text(
    data = corr_values,
    aes(
      x = max(corr_data$mean_network),
      y = ypos,
      label = label,
      color = condition
    ),
    hjust = 1
  )



# Stepwise - multiple regression
all_data <- merge(adjusted_cfos, test_wide, by = "ID")

all_data$mean_HC <- rowMeans(all_data[, c("CA1dist",
                                          "CA1prox",
                                          "CA3",
                                          'dentate')], na.rm = TRUE)

all_data$mean_CG <- rowMeans(all_data[, c("anterior_CG1",
                                          "medial_CG1",
                                          "medial_CG2")], na.rm = TRUE)

all_data$mean_retro <- rowMeans(all_data[, c("anterior_retrosplenial_agranular",
                                             "anterior_retrosplenial_granular",
                                             "medial_retrosplenial_agranular",
                                             "medial_retrosplenial_granular",
                                             "posterior_retrosplenial_agranular",
                                             "posterior_retrosplenial_granular")], na.rm = TRUE)

all_data <- subset(all_data, all_data$condition == "sleep")

models <- list(
  HC = lm(Cum_DiRa_min_2 ~ mean_HC, data = all_data),
  HC_CG = lm(Cum_DiRa_min_2 ~ mean_HC + mean_CG, data = all_data),
  HC_CG_retro = lm(Cum_DiRa_min_2 ~ mean_HC + mean_CG + mean_retro, data = all_data)
)

# Extract R2 and adj R2
model_stats <- data.frame(
  Model = names(models),
  R2 = sapply(models, function(m) summary(m)$r.squared),
  Adj_R2 = sapply(models, function(m) summary(m)$adj.r.squared),
  AIC = sapply(models, AIC)
)

print(model_stats)


ggplot(model_stats, aes(x = Model)) +
  geom_col(aes(y = R2), fill = "steelblue", alpha = 0.7) +
  geom_point(aes(y = Adj_R2), color = "red", size = 3) +
  geom_text(aes(y = R2, label = round(R2, 2)), vjust = -0.5) +
  geom_text(aes(y = Adj_R2, label = round(Adj_R2, 2)), vjust = 1.5, color = "red") +
  ylab("Variance explained") +
  ggtitle("Raw R² (bars) and Adjusted R² (red points)") +
  theme_minimal()


models <- list(
  model_CA3 = lm(Cum_DiRa_min_2 ~ CA3, data = all_data),
  model_CA3_aCG1 = lm(Cum_DiRa_min_2 ~ CA3 + anterior_CG1, data = all_data),
  model_CA3_aCG1_aRSG = lm(Cum_DiRa_min_2 ~ CA3 + anterior_CG1 + anterior_retrosplenial_granular, data = all_data),
  model_CA3_aCG1_aRSG = lm(Cum_DiRa_min_2 ~ CA3 + anterior_CG1 + anterior_retrosplenial_granular, data = all_data)
)

# Extract R2 and adj R2
model_stats <- data.frame(
  Model = names(models),
  R2 = sapply(models, function(m) summary(m)$r.squared),
  Adj_R2 = sapply(models, function(m) summary(m)$adj.r.squared),
  AIC = sapply(models, AIC)
)

print(model_stats)

ggplot(model_stats, aes(x = Model)) +
  geom_col(aes(y = R2), fill = "steelblue", alpha = 0.7) +
  geom_point(aes(y = Adj_R2), color = "red", size = 3) +
  geom_text(aes(y = R2, label = round(R2, 2)), vjust = -0.5) +
  geom_text(aes(y = Adj_R2, label = round(Adj_R2, 2)), vjust = 1.5, color = "red") +
  ylab("Variance explained") +
  ggtitle("Raw R² (bars) and Adjusted R² (red points)") +
  theme_minimal()

