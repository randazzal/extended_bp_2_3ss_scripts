library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggeasy)
library(readxl)
library(writexl)

mic <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/detailed_list.txt", sep = "\t")
avg <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/avg_regroup.txt", sep = "\t")
avg <- avg[ ,c("EVENT", "tissue")]

# exon A 3ss
nn_micro <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/conservation/nn_micro_A.csv", sep = ",")
nn_micro[,1] <- NULL
neu_micro <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/conservation/neu_micro_A.csv", sep = ",")
neu_micro[,1] <- NULL
micro_pc <- cbind(neu_micro, nn_micro)
rm(neu_micro)
rm(nn_micro)
nn_avg <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/conservation/nn_avg_A.csv", sep = ",")
nn_avg[,1] <- NULL
neu_avg <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/conservation/neu_avg_A.csv", sep = ",")
neu_avg[,1] <- NULL
avg_pc <- cbind(neu_avg, nn_avg)
rm(neu_avg)
rm(nn_avg)

##### separate plus and minus strand for each size class
micro_plus <- read_excel(path = "/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/intron_length/no_neuro_lengths.xlsx", sheet = "micro")
micro_minus <- micro_plus %>%
  filter(micro_plus$strand == "-")
micro_plus <- micro_plus %>%
  filter(micro_plus$strand == "+")
avg_plus <- read_excel(path = "/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/intron_length/no_neuro_lengths.xlsx", sheet = "avg")
avg_minus <- avg_plus %>%
  filter(avg_plus$strand == "-")
avg_plus <- avg_plus %>%
  filter(avg_plus$strand == "+")

# fix directionality
plus_micro <- micro_pc[ ,colnames(micro_pc) %in% micro_plus$EVENT]
plus_micro$position <- -450:450
minus_micro <- micro_pc[ ,colnames(micro_pc) %in% micro_minus$EVENT]
minus_micro$position <- 450:-450
micro_score <- merge(plus_micro, minus_micro, by="position")
micro_score$position <- NULL
rownames(micro_score) <- -450:450
rm(micro_minus)
rm(micro_plus)
rm(plus_micro)
rm(minus_micro)
rm(micro_pc)

plus_avg <- avg_pc[ ,colnames(avg_pc) %in% avg_plus$EVENT]
plus_avg$position <- -450:450
minus_avg <- avg_pc[ ,colnames(avg_pc) %in% avg_minus$EVENT]
minus_avg$position <- 450:-450
avg_score <- merge(plus_avg, minus_avg, by="position")
avg_score$position <- NULL
rownames(avg_score) <- -450:450
rm(avg_minus)
rm(avg_plus)
rm(plus_avg)
rm(minus_avg)
rm(avg_pc)

# separate by tissue expression 
rare_list <- mic[mic$tissue == "rare", ]
rare_micro <- micro_score[ ,colnames(micro_score) %in% rare_list$EVENT]
common_list <- mic[mic$tissue == "common", ]
common_micro <- micro_score[ ,colnames(micro_score) %in% common_list$EVENT]
higher_list <- mic[mic$tissue == "fdr_g", ]
higher_micro <- micro_score[ ,colnames(micro_score) %in% higher_list$EVENT]
lower_list <- mic[mic$tissue == "fdr_l", ]
lower_micro <- micro_score[ ,colnames(micro_score) %in% lower_list$EVENT]
other_list <- mic[mic$tissue == "other", ]
other_micro <- micro_score[ ,colnames(micro_score) %in% other_list$EVENT]
neuro_list <- rbind(higher_list, lower_list)
neuro_micro <- micro_score[ ,colnames(micro_score) %in% neuro_list$EVENT]
nn_list <- rbind(rare_list, common_list, other_list)
nn_micro <- micro_score[ ,colnames(micro_score) %in% nn_list$EVENT]

rare_list <- avg[avg$tissue == "rare", ]
rare_avg <- avg_score[ ,colnames(avg_score) %in% rare_list$EVENT]
common_list <- avg[avg$tissue == "common", ]
common_avg <- avg_score[ ,colnames(avg_score) %in% common_list$EVENT]
higher_list <- avg[avg$tissue == "fdr_g", ]
higher_avg <- avg_score[ ,colnames(avg_score) %in% higher_list$EVENT]
lower_list <- avg[avg$tissue == "fdr_l", ]
lower_avg <- avg_score[ ,colnames(avg_score) %in% lower_list$EVENT]
other_list <- avg[avg$tissue == "other", ]
other_avg <- avg_score[ ,colnames(avg_score) %in% other_list$EVENT]
neuro_list <- rbind(higher_list, lower_list)
neuro_avg <- avg_score[ ,colnames(avg_score) %in% neuro_list$EVENT]
nn_list <- rbind(rare_list, common_list, other_list)
nn_avg <- avg_score[ ,colnames(avg_score) %in% nn_list$EVENT]

#neural
lim_avg <- neuro_avg %>%
  mutate(position = -450:450)
lim_avg <- lim_avg[lim_avg$position >= -250 & lim_avg$position <= -1, ]

lim_micro <- neuro_micro %>%
  mutate(position = -450:450)
lim_micro <- lim_micro[lim_micro$position >= -250 & lim_micro$position <= -1, ]

matching_cols <- grep("^HsaEX", colnames(lim_avg), value = TRUE)
narrow_avg <- lim_avg %>%
  pivot_longer(cols = all_of(matching_cols),
               names_to = "event",
               values_to = "value") 
narrow_avg$group <- "Neural Midexon"

matching_cols <- grep("^HsaEX", colnames(lim_micro), value = TRUE)
narrow_micro <- lim_micro %>%
  pivot_longer(cols = all_of(matching_cols),
               names_to = "event",
               values_to = "value")
narrow_micro$group <- "Neural Microexon"
less <- rbind(narrow_avg, narrow_micro)

medians_avg_neu <- narrow_avg %>%
  filter(position >= -150 & position <= -1) %>%
  group_by(position, group) %>%
  summarise(
    med = median(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

medians_micro_neu <- narrow_micro %>%
  filter(position >= -150 & position <= -1) %>%
  group_by(position, group) %>%
  summarise(
    med = median(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

medians_avg_neu <- medians_avg_neu %>% mutate(type = "Neural Midexon")
medians_micro_neu <- medians_micro_neu %>% mutate(type = "Neural Microexon")

#non_neural
lim_avg <- nn_avg %>%
  mutate(position = -450:450)
lim_avg <- lim_avg[lim_avg$position >= -250 & lim_avg$position <= -1, ]

lim_micro <- nn_micro %>%
  mutate(position = -450:450)
lim_micro <- lim_micro[lim_micro$position >= -250 & lim_micro$position <= -1, ]

matching_cols <- grep("^HsaEX", colnames(lim_avg), value = TRUE)
narrow_avg <- lim_avg %>%
  pivot_longer(cols = all_of(matching_cols),
               names_to = "event",
               values_to = "value") 
narrow_avg$group <- "Non-Neural Midexon"

matching_cols <- grep("^HsaEX", colnames(lim_micro), value = TRUE)
narrow_micro <- lim_micro %>%
  pivot_longer(cols = all_of(matching_cols),
               names_to = "event",
               values_to = "value")
narrow_micro$group <- "Non-Neural Microexon"
combo <- rbind(narrow_avg, narrow_micro)

medians_avg_nn <- narrow_avg %>%
  filter(position >= -150 & position <= -1) %>%
  group_by(position, group) %>%
  summarise(
    med = median(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

medians_micro_nn <- narrow_micro %>%
  filter(position >= -150 & position <= -1) %>%
  group_by(position, group) %>%
  summarise(
    med = median(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

medians_avg_nn <- medians_avg_nn %>% mutate(type = "Non-Neural Midexon")
medians_micro_nn <- medians_micro_nn %>% mutate(type = "Non-Neural Microexon")

ggplot() +
  # Microexon ribbon
  geom_ribbon(
    data = medians_micro_neu,
    aes(x = position, ymin = mean - sd, ymax = mean + sd, fill = type),
    alpha = 0.25
  ) +
  geom_ribbon(
    data = medians_micro_nn,
    aes(x = position, ymin = mean - sd, ymax = mean + sd, fill = type),
    alpha = 0.25
  ) +
  # Avg exon ribbon
  geom_ribbon(
    data = medians_avg_neu,
    aes(x = position, ymin = mean - sd, ymax = mean + sd, fill = type),
    alpha = 0.25
  ) +
  geom_ribbon(
    data = medians_avg_nn,
    aes(x = position, ymin = mean - sd, ymax = mean + sd, fill = type),
    alpha = 0.25
  ) +
  # Microexon mean line
  geom_line(
    data = medians_micro_neu,
    aes(x = position, y = mean, color = type),
    linewidth = 1
  ) +
  geom_line(
    data = medians_micro_nn,
    aes(x = position, y = mean, color = type),
    linewidth = 1
  ) +
  # Avg exon mean line
  geom_line(
    data = medians_avg_neu,
    aes(x = position, y = mean, color = type),
    linewidth = 1
  ) +
  geom_line(
    data = medians_avg_nn,
    aes(x = position, y = mean, color = type),
    linewidth = 1
  ) +
  # Manual colors — guaranteed stable
  scale_color_manual(values = c(
    "Neural Midexon"   = "#5e3c80",
    "Neural Microexon" = "#fdb863",
    "Non-Neural Midexon"   = "#b2abd9",
    "Non-Neural Microexon" = "#e66115"
  )) +
  scale_fill_manual(values = c(
    "Neural Midexon"   = "#5e3c80",
    "Neural Microexon" = "#fdb863",
    "Non-Neural Midexon"   = "#b2abd9",
    "Non-Neural Microexon" = "#e66115"
  )) +
  
  scale_x_continuous(breaks = c(-150, -100, -50)) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 1))

mega <- rbind(less, combo)
compare_groups <- function(df, g1, g2) {
  df %>%
    filter(group %in% c(g1, g2)) %>%
    group_by(position) %>%
    summarise(p_value = wilcox.test(value ~ group)$p.value, .groups = "drop")
}
p_nmicro_navg <- compare_groups(mega, "Neural Microexon", "Neural Midexon")
p_nmicro_nonm <- compare_groups(mega, "Neural Microexon", "Non-Neural Microexon")
p_nmicro_nonavg <- compare_groups(mega, "Neural Microexon", "Non-Neural Midexon")

p_nmicro_navg$group <- c("neumicro neuavg")
p_nmicro_nonm$group <- c("neumicro nnmicro")
p_nmicro_nonavg$group <- c("neumicro nnavg")
merged <- rbind(p_nmicro_navg, p_nmicro_nonm, p_nmicro_nonavg)

# Merge results
p_values_df <- purrr::reduce(list(
  p_nmicro_navg, #p_value.x
  p_nmicro_nonm, #p_value.y
  p_nmicro_nonavg #p_value
), dplyr::full_join, by = "position")
p_values_df <- p_values_df %>%
  mutate(across(starts_with("p_"), ~ p.adjust(.x, method = "fdr")))
p_values_df %>%
  pivot_longer(cols = starts_with("p_"), names_to = "comparison", values_to = "p_value") %>%
  ggplot(aes(x = position, y = p_value, color = comparison)) +
  geom_point() +
  geom_hline(yintercept = 0.001, linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "pink") +
  theme_minimal() +
  ylim(0,1) +
  xlim(-150,0) +
  labs(title = "Significance of Microexon vs. Other Groups at Each Position",
       x = "Position",
       y = "p-value",
       color = "Comparison")
