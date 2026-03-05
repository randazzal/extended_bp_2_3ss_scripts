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
mega_mid <- rbind(mega_midUGC, mega_mid_noUGC)

################################ Figures 4E and 4F
micro_trim <- mega_micro[, c("IncLevel_mean_npc", "IncLevel_mean_neuron", "TGC_count", "group", "class", "tissue")]
micro_trim_long <- micro_trim %>%
  pivot_longer(cols = c("IncLevel_mean_npc", "IncLevel_mean_neuron"),
               values_to = "PSI", names_to = "stage")
micro_trim_long$stage <- factor(micro_trim_long$stage, levels = c("IncLevel_mean_npc", "IncLevel_mean_neuron"))
micro_trim_long %>%
  #filter(group == "neural") %>%
  filter(tissue == "fdr_l") %>%
  ggplot(aes(x = stage, y = PSI, fill = class)) +
  geom_boxplot(position = "dodge") +
  theme_minimal()
test <- wilcox.test(micro_trim_long[micro_trim_long$class == "UGC" & micro_trim_long$group == "neural" & micro_trim_long$stage == "IncLevel_mean_npc",]$PSI,
                    micro_trim_long[micro_trim_long$class == "UGC" & micro_trim_long$group == "neural" & micro_trim_long$stage == "IncLevel_mean_neuron",]$PSI) #5.49e-44
test <- wilcox.test(micro_trim_long[micro_trim_long$class == "No UGC" & micro_trim_long$group == "neural" & micro_trim_long$stage == "IncLevel_mean_npc",]$PSI,
                    micro_trim_long[micro_trim_long$class == "UGC" & micro_trim_long$group == "neural" & micro_trim_long$stage == "IncLevel_mean_npc",]$PSI) #0.00194
test <- wilcox.test(micro_trim_long[micro_trim_long$class == "No UGC" & micro_trim_long$group == "neural" & micro_trim_long$stage == "IncLevel_mean_npc",]$PSI,
                    micro_trim_long[micro_trim_long$class == "No UGC" & micro_trim_long$group == "neural" & micro_trim_long$stage == "IncLevel_mean_neuron",]$PSI) #1.05e-6
test <- wilcox.test(micro_trim_long[micro_trim_long$class == "No UGC" & micro_trim_long$group == "neural" & micro_trim_long$stage == "IncLevel_mean_neuron",]$PSI,
                    micro_trim_long[micro_trim_long$class == "UGC" & micro_trim_long$group == "neural" & micro_trim_long$stage == "IncLevel_mean_neuron",]$PSI) #0.00579

mid_trim <- mega_mid[, c("IncLevel_mean_npc", "IncLevel_mean_neuron", "TGC_count", "group", "class")]
mid_trim_long <- mid_trim %>%
  pivot_longer(cols = c("IncLevel_mean_npc", "IncLevel_mean_neuron"),
               values_to = "PSI", names_to = "stage")
mid_trim_long$stage <- factor(mid_trim_long$stage, levels = c("IncLevel_mean_npc", "IncLevel_mean_neuron"))
mid_trim_long %>%
  filter(group == "neural") %>%
  ggplot(aes(x = stage, y = PSI, fill = class)) +
  geom_boxplot(position = "dodge") +
  theme_minimal()
test <- wilcox.test(mid_trim_long[mid_trim_long$class == "UGC" & mid_trim_long$group == "neural" & mid_trim_long$stage == "IncLevel_mean_npc",]$PSI,
                    mid_trim_long[mid_trim_long$class == "UGC" & mid_trim_long$group == "neural" & mid_trim_long$stage == "IncLevel_mean_neuron",]$PSI) #0.908
test <- wilcox.test(mid_trim_long[mid_trim_long$class == "No UGC" & mid_trim_long$group == "neural" & mid_trim_long$stage == "IncLevel_mean_npc",]$PSI,
                    mid_trim_long[mid_trim_long$class == "UGC" & mid_trim_long$group == "neural" & mid_trim_long$stage == "IncLevel_mean_npc",]$PSI) #0.344
test <- wilcox.test(mid_trim_long[mid_trim_long$class == "No UGC" & mid_trim_long$group == "neural" & mid_trim_long$stage == "IncLevel_mean_npc",]$PSI,
                    mid_trim_long[mid_trim_long$class == "No UGC" & mid_trim_long$group == "neural" & mid_trim_long$stage == "IncLevel_mean_neuron",]$PSI) #0.144
test <- wilcox.test(mid_trim_long[mid_trim_long$class == "No UGC" & mid_trim_long$group == "neural" & mid_trim_long$stage == "IncLevel_mean_neuron",]$PSI,
                    mid_trim_long[mid_trim_long$class == "UGC" & mid_trim_long$group == "neural" & mid_trim_long$stage == "IncLevel_mean_neuron",]$PSI) #0.665
