library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggeasy)
library(readxl)
library(writexl)

final <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/detailed_list.txt", sep = "\t")
noneuro_micro <- read_excel(path = "/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/intron_length/no_neuro_lengths.xlsx", sheet = "micro")
colnames(noneuro_micro) <- c("EVENT", "GENE", "chromosome", "C15SS", "C23SS", "strand", "exonA_3ss", "exonA_5ss", "up", "down")
noneuro_micro$up_length <- noneuro_micro$up - 2
noneuro_micro$down_length <- noneuro_micro$down - 2
micro <- noneuro_micro %>%
  filter(EVENT %in% final$EVENT)
micro <- merge(micro, final, by = "EVENT")

final <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/avg_regroup.txt", sep = "\t")
final <- final[ ,c("EVENT", "tissue")]
noneuro_avg <- read_excel(path = "/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/intron_length/no_neuro_lengths.xlsx", sheet = "avg")
colnames(noneuro_avg) <- c("EVENT", "GENE", "chromosome", "C15SS", "C23SS", "strand", "exonA_3ss", "exonA_5ss", "up", "down")
noneuro_avg$up_length <- noneuro_avg$up - 2
noneuro_avg$down_length <- noneuro_avg$down - 2
avg <- noneuro_avg %>%
  filter(EVENT %in% final$EVENT)
avg <- merge(avg, final, by = "EVENT")

#upstream intron
m_up <- micro[ ,c("EVENT", "up_length", "tissue")] %>%
  mutate(group = "micro")
a_up <- avg[ ,c("EVENT", "up_length", "tissue")] %>%
  mutate(group = "average")
avg_neu <- a_up %>%
  filter(tissue != "common" & tissue != "rare" & tissue != "other") %>%
  mutate(group = "Human Neural Midexon")
mic_neu <- m_up %>%
  filter(tissue != "common" & tissue != "rare" & tissue != "other") %>%
  mutate(group = "Human Neural Microexon")
avg_nn <- a_up %>%
  filter(tissue != "fdr_l" & tissue != "fdr_g") %>%
  mutate(group = "Human Non-Neural Midexon")
mic_nn <- m_up %>%
  filter(tissue != "fdr_l" & tissue != "fdr_g") %>%
  mutate(group = "Human Non-Neural Microexon")

combo <- rbind(avg_neu, avg_nn, mic_neu, mic_nn)
combo %>%
  ggplot(aes(x = group, y = up_length, fill = group)) +
  geom_boxplot(fill = "#E6E6FA") +
  theme_minimal() +
  labs(title = "Upstream Intron Lengths", 
       x = "Exon Class", 
       y = "Length (nts)") +
  easy_center_title() +
  ylim(0,15000)

combo$log_length <- log10(combo$up_length)
combo %>%
  ggplot(aes(x = group, y = log_length, fill = group)) +
  geom_boxplot(fill = "#E6E6FA") +
  theme_minimal() +
  labs(title = "Upstream Intron Lengths", 
       x = "Exon Class", 
       y = "Length (nts) (log10)") +
  easy_center_title() +
  ylim(0,6)

test <- wilcox.test(mic_neu$up_length, avg_neu$up_length) #0.0206
test <- wilcox.test(mic_neu$up_length, avg_nn$up_length) #0.00663
test <- wilcox.test(mic_neu$up_length, mic_nn$up_length) #0.0561
test <- wilcox.test(mic_nn$up_length, avg_nn$up_length) #1.05e-10
test <- wilcox.test(mic_nn$up_length, avg_neu$up_length) #1.07e-07
test <- wilcox.test(avg_neu$up_length, avg_nn$up_length) #0.786

#downstream intron
m_down <- micro[ ,c("EVENT", "down_length", "tissue")]
a_down <- avg[ ,c("EVENT", "down_length", "tissue")] 
avg_neu <- a_down %>%
  filter(tissue != "common" & tissue != "rare" & tissue != "other") %>%
  mutate(group = "Human Neural Midexon")
mic_neu <- m_down %>%
  filter(tissue != "common" & tissue != "rare" & tissue != "other") %>%
  mutate(group = "Human Neural Microexon")
avg_nn <- a_down %>%
  filter(tissue != "fdr_l" & tissue != "fdr_g") %>%
  mutate(group = "Human Non-Neural Midexon")
mic_nn <- m_down %>%
  filter(tissue != "fdr_l" & tissue != "fdr_g") %>%
  mutate(group = "Human Non-Neural Microexon")

combo <- rbind(avg_neu, avg_nn, mic_neu, mic_nn)
combo %>%
  ggplot(aes(x = group, y = down_length, fill = group)) +
  geom_boxplot(fill = "#E6E6FA") +
  theme_minimal() +
  labs(title = "Downstream Intron Lengths", 
       x = "Exon Class", 
       y = "Length (nts)") +
  easy_center_title() +
  ylim(0,15000)

combo$log_length <- log10(combo$down_length)
combo %>%
  ggplot(aes(x = group, y = log_length, fill = group)) +
  geom_boxplot(fill = "#E6E6FA") +
  theme_minimal() +
  labs(title = "Upstream Intron Lengths", 
       x = "Exon Class", 
       y = "Length (nts) (log10)") +
  easy_center_title() +
  ylim(0,6)

test <- wilcox.test(mic_neu$down_length, avg_neu$down_length) #0.000361
test <- wilcox.test(mic_neu$down_length, avg_nn$down_length) #1.47e-05
test <- wilcox.test(mic_neu$down_length, mic_nn$down_length) #0.313
test <- wilcox.test(mic_nn$down_length, avg_nn$down_length) #5.85e-12
test <- wilcox.test(avg_neu$down_length, mic_nn$down_length) #5.95e-08
test <- wilcox.test(avg_neu$down_length, avg_nn$down_length) #0.558
