library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(ggeasy)

micro_noUGC <- read.delim("/data2/lackey_lab/randazza/microexons/motifs/UGC_motif_chick/DIY/hsa_vastDB_micro_noUGC.txt", sep = "\t")
micro_UGC <- read.delim("/data2/lackey_lab/randazza/microexons/motifs/UGC_motif_chick/DIY/hsa_vastDB_micro_UGC.txt", sep = "\t")
mid_noUGC <- read.delim("/data2/lackey_lab/randazza/microexons/motifs/UGC_motif_chick/DIY/hsa_vastDB_mid_noUGC.txt", sep = "\t")
mid_UGC <- read.delim("/data2/lackey_lab/randazza/microexons/motifs/UGC_motif_chick/DIY/hsa_vastDB_mid_UGC.txt", sep = "\t")

psi_values <- read.delim("/data2/lackey_lab/randazza/microexons/devo/human/keio/all_psi_values.txt", sep = "\t")

mega_midUGC <- merge(psi_values, mid_UGC, by.x = "name", by.y = "EVENT.y")
mega_microUGC <- merge(psi_values, micro_UGC, by.x = "name", by.y = "EVENT.y")
mega_mid_noUGC <- merge(psi_values, mid_noUGC, by.x = "name", by.y = "EVENT.y")
mega_micro_noUGC <- merge(psi_values, micro_noUGC, by.x = "name", by.y = "EVENT.y")

add_group_column <- function(df) {
  df %>%
    mutate(group = ifelse(tissue %in% c("fdr_g", "fdr_l"),
                          "neural",
                          "non-neural"))
}
mega_microUGC <- add_group_column(mega_microUGC)
mega_microUGC$class <- "UGC"
mega_micro_noUGC <- add_group_column(mega_micro_noUGC)
mega_micro_noUGC$class <- "No UGC"
mega_midUGC <- add_group_column(mega_midUGC)
mega_midUGC$class <- "UGC"
mega_mid_noUGC <- add_group_column(mega_mid_noUGC)
mega_mid_noUGC$class <- "No UGC"

mega_micro <- rbind(mega_microUGC, mega_micro_noUGC)
mega_micro$platform <- mega_micro$ss_dist + mega_micro$length
mega_mic_neu <- mega_micro %>%
  filter(group == "neural")
mega_mic_neu_lim <- mega_mic_neu[mega_mic_neu$IncLevel_mean_neuron < 0.98 & mega_mic_neu$IncLevel_mean_neuron > 0.02
                                 & mega_mic_neu$IncLevel_mean_npc > 0.02 & mega_mic_neu$IncLevel_mean_npc < 0.98, ]
mega_mid <- rbind(mega_midUGC, mega_mid_noUGC)
mega_mid$platform <- mega_mid$ss_dist + mega_mid$length
mega_mid_neu <- mega_mid %>%
  filter(group == "neural")
mega_mid_neu_lim <- mega_mid_neu[mega_mid_neu$IncLevel_mean_neuron < 0.98 & mega_mid_neu$IncLevel_mean_neuron > 0.02
                                 & mega_mid_neu$IncLevel_mean_npc > 0.02 & mega_mid_neu$IncLevel_mean_npc < 0.98, ]

temp <- mega_mic_neu[, c("class", "ppt_len", "ss_dist", "TGC_count", "platform")]
temp2 <- mega_mid_neu[, c("class", "ppt_len", "ss_dist", "TGC_count", "platform")]
temp$group <- "Microexon"
temp$set <- paste(temp$group, temp$class, sep = "_")
temp2$group <- "Midexon"
temp2$set <- paste(temp2$group, temp2$class, sep = "_")
temp <- rbind(temp, temp2)
temp %>%
  ggplot(aes(x = set, y = ss_dist)) +
  geom_boxplot() +
  theme_minimal()
test <- wilcox.test(temp[temp$set == "Microexon_No UGC",]$ss_dist,
                    temp[temp$set == "Microexon_UGC",]$ss_dist) #2.39e-6
test <- wilcox.test(temp[temp$set == "Midexon_No UGC",]$ss_dist,
                    temp[temp$set == "Midexon_UGC",]$ss_dist) #2.91e-13
test <- wilcox.test(temp[temp$set == "Microexon_UGC",]$ss_dist,
                    temp[temp$set == "Midexon_UGC",]$ss_dist) #0.0989
test <- wilcox.test(temp[temp$set == "Microexon_No UGC",]$ss_dist,
                    temp[temp$set == "Midexon_No UGC",]$ss_dist) #0.615
temp %>%
  ggplot(aes(x = set, y = ppt_len)) +
  geom_boxplot() +
  theme_minimal()
test <- wilcox.test(temp[temp$set == "Microexon_No UGC",]$ppt_len,
                    temp[temp$set == "Microexon_UGC",]$ppt_len) #0.00504
test <- wilcox.test(temp[temp$set == "Midexon_No UGC",]$ppt_len,
                    temp[temp$set == "Midexon_UGC",]$ppt_len) #0.0453
test <- wilcox.test(temp[temp$set == "Microexon_UGC",]$ppt_len,
                    temp[temp$set == "Midexon_UGC",]$ppt_len) #1.1e-16
test <- wilcox.test(temp[temp$set == "Microexon_No UGC",]$ppt_len,
                    temp[temp$set == "Midexon_No UGC",]$ppt_len) #0.001
temp %>%
  ggplot(aes(x = TGC_count, y = platform, color = set)) +
  geom_point() +
  theme_minimal()