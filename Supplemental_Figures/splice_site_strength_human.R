setwd("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/SS_strength/")
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggeasy)
library(readxl)
library(writexl)

final <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/detailed_list.txt", sep = "\t")
limited <- read.delim("micro_A3ss.bed", sep = "\t", header = FALSE)
final <- final %>%
  filter(EVENT %in% limited$V4)
final <- final[order(final$EVENT), ]
rm(limited)

micro3 <- read.delim("micro_A3ss_score.txt", sep = "\t", header = FALSE)
colnames(micro3) <- c("sequence", "3ss_score")
micro5 <- read.delim("micro_A5ss_score.txt", sep = "\t", header = FALSE)
colnames(micro5) <- c("sequence5", "5ss_score")
micro <- cbind(final, micro3, micro5)
micro$exon <- "A"
microC1 <- read.delim("micro_C15ss_score.txt", sep = "\t", header = FALSE)
colnames(microC1) <- c("sequencec15", "5ss_score")
microC2 <- read.delim("micro_C23ss_score.txt", sep = "\t", header = FALSE)
colnames(microC2) <- c("sequencec23", "3ss_score")

micro_A3 <- micro[,c("EVENT", "3ss_score", "exon", "tissue")]
micro_A5 <- micro[,c("EVENT", "5ss_score", "exon", "tissue")]

micro_C1 <- cbind(final, microC1)
micro_C1$exon <- "C"
microC1 <- micro_C1[,c("EVENT", "5ss_score", "exon", "tissue")]
micro_C2 <- cbind(final, microC2)
micro_C2$exon <- "C"
microC2 <- micro_C2[,c("EVENT", "3ss_score", "exon", "tissue")]
rm(micro_C1)
rm(micro_C2)
rm(micro3)
rm(micro5)

micro_neu <- micro[micro$tissue != "common" & micro$tissue != "rare" & micro$tissue != "other", ]
micro_neu$group <- "Neural Microexon"
micro_neu <- micro_neu[, c("EVENT", "tissue", "sequence", "3ss_score", "sequence5",
                           "5ss_score", "exon", "group")]
micro_nn <- micro[micro$tissue != "fdr_g" & micro$tissue != "fdr_l", ]
micro_nn$group <- "Non-Neural Microexon"
micro_nn <- micro_nn[, c("EVENT", "tissue", "sequence", "3ss_score", "sequence5",
                           "5ss_score", "exon", "group")]
micro_neu_C1 <- microC1[microC1$tissue != "common" & microC1$tissue != "rare" & microC1$tissue != "other", ]
micro_neu_C1$group <- "Neural Microexon"
micro_nn_C1 <- microC1[microC1$tissue != "fdr_g" & microC1$tissue != "fdr_l", ]
micro_nn_C1$group <- "Non-Neural Microexon"
micro_neu_C2 <- microC2[microC2$tissue != "common" & microC2$tissue != "rare" & microC2$tissue != "other", ]
micro_neu_C2$group <- "Neural Microexon"
micro_nn_C2 <- microC2[microC2$tissue != "fdr_g" & microC2$tissue != "fdr_l", ]
micro_nn_C2$group <- "Non-Neural Microexon"

final <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/avg_regroup.txt", sep = "\t")
final <- final[ ,c("EVENT", "tissue")]
limited <- read.delim("avg_A3ss.bed", sep = "\t", header = FALSE)
final <- final %>%
  filter(EVENT %in% limited$V4)
final <- final[order(final$EVENT), ]
rm(limited)

avg3 <- read.delim("avg_A3ss_score.txt", sep = "\t", header = FALSE)
colnames(avg3) <- c("sequence", "3ss_score")
avg5 <- read.delim("avg_A5ss_score.txt", sep = "\t", header = FALSE)
colnames(avg5) <- c("sequence5", "5ss_score")
avg <- cbind(final, avg3, avg5)
avg$exon <- "A"
avgC1 <- read.delim("avg_C15ss_score.txt", sep = "\t", header = FALSE)
colnames(avgC1) <- c("sequencec15", "5ss_score")
avgC2 <- read.delim("avg_C23ss_score.txt", sep = "\t", header = FALSE)
colnames(avgC2) <- c("sequencec23", "3ss_score")

avg_A3 <- avg[,c("EVENT", "3ss_score", "exon", "tissue")]
avg_A5 <- avg[,c("EVENT", "5ss_score", "exon", "tissue")]

avg_C1 <- cbind(final, avgC1)
avg_C1$exon <- "C"
avgC1 <- avg_C1[,c("EVENT", "5ss_score", "exon", "tissue")]
avg_C2 <- cbind(final, avgC2)
avg_C2$exon <- "C"
avgC2 <- avg_C2[,c("EVENT", "3ss_score", "exon", "tissue")]
rm(avg_C1)
rm(avg_C2)
rm(avg3)
rm(avg5)

avg_neu <- avg[avg$tissue != "common" & avg$tissue != "rare" & avg$tissue != "other", ]
avg_neu$group <- "Neural Midexon"
avg_nn <- avg[avg$tissue != "fdr_g" & avg$tissue != "fdr_l", ]
avg_nn$group <- "Non-Neural Midexon"
avg_neu_C1 <- avgC1[avgC1$tissue != "common" & avgC1$tissue != "rare" & avgC1$tissue != "other", ]
avg_neu_C1$group <- "Neural Midexon"
avg_nn_C1 <- avgC1[avgC1$tissue != "fdr_g" & avgC1$tissue != "fdr_l", ]
avg_nn_C1$group <- "Non-Neural Midexon"
avg_neu_C2 <- avgC2[avgC2$tissue != "common" & avgC2$tissue != "rare" & avgC2$tissue != "other", ]
avg_neu_C2$group <- "Neural Midexon"
avg_nn_C2 <- avgC2[avgC2$tissue != "fdr_g" & avgC2$tissue != "fdr_l", ]
avg_nn_C2$group <- "Non-Neural Midexon"

combo <- rbind(avg_neu, avg_nn, micro_neu, micro_nn)
combo %>%
  ggplot(aes(x = group, y = `3ss_score`, fill = group)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c("Neural Midexon" = "#5e3c99", "Non-Neural Microexon" = "#e66101", "Non-Neural Midexon" = "#b2abd2", "Neural Microexon" = "#fdb863")) +
  theme(legend.position = "none") +
  ylim(-35,15)
combo %>%
  ggplot(aes(x = group, y = `5ss_score`, fill = group)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c("Neural Midexon" = "#5e3c99", "Non-Neural Microexon" = "#e66101", "Non-Neural Midexon" = "#b2abd2", "Neural Microexon" = "#fdb863")) +
  theme(legend.position = "none") +
  ylim(-35,15)
test <- wilcox.test(micro_neu$`3ss_score`, avg_neu$`3ss_score`) #0.00171
test <- wilcox.test(micro_neu$`3ss_score`, avg_nn$`3ss_score`) #0.0015
test <- wilcox.test(micro_neu$`3ss_score`, micro_nn$`3ss_score`) #0.158
test <- wilcox.test(micro_nn$`3ss_score`, avg_nn$`3ss_score`) #6.31e-15
test <- wilcox.test(micro_nn$`3ss_score`, avg_neu$`3ss_score`) #1.02e-11
test <- wilcox.test(avg_nn$`3ss_score`, avg_neu$`3ss_score`) #0.556

test <- wilcox.test(micro_neu$`5ss_score`, avg_neu$`5ss_score`) #0.000661
test <- wilcox.test(micro_neu$`5ss_score`, avg_nn$`5ss_score`) #1.87e-06
test <- wilcox.test(micro_neu$`5ss_score`, micro_nn$`5ss_score`) #0.94
test <- wilcox.test(micro_nn$`5ss_score`, avg_nn$`5ss_score`) #9.25e-12
test <- wilcox.test(micro_nn$`5ss_score`, avg_neu$`5ss_score`) #9.85e-06
test <- wilcox.test(avg_nn$`5ss_score`, avg_neu$`5ss_score`) #0.02

combo <- rbind(avg_neu_C1, avg_nn_C1, micro_neu_C1, micro_nn_C1)
combo %>%
  ggplot(aes(x = group, y = `5ss_score`, fill = group)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c("Neural Midexon" = "#5e3c99", "Non-Neural Microexon" = "#e66101", "Non-Neural Midexon" = "#b2abd2", "Neural Microexon" = "#fdb863")) +
  theme(legend.position = "none") +
  ylim(-35,15)
test <- wilcox.test(micro_neu_C1$`5ss_score`, avg_neu_C1$`5ss_score`) #0.329
test <- wilcox.test(micro_neu_C1$`5ss_score`, avg_nn_C1$`5ss_score`) #0.122
test <- wilcox.test(micro_neu_C1$`5ss_score`, micro_nn_C1$`5ss_score`) #0.0507
test <- wilcox.test(micro_nn_C1$`5ss_score`, avg_nn_C1$`5ss_score`) #0.267
test <- wilcox.test(micro_nn_C1$`5ss_score`, avg_neu_C1$`5ss_score`) #0.143
test <- wilcox.test(avg_neu_C1$`5ss_score`, avg_nn_C1$`5ss_score`) #0.391

combo <- rbind(avg_neu_C2, avg_nn_C2, micro_neu_C2, micro_nn_C2)
combo %>%
  ggplot(aes(x = group, y = `3ss_score`, fill = group)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c("Neural Midexon" = "#5e3c99", "Non-Neural Microexon" = "#e66101", "Non-Neural Midexon" = "#b2abd2", "Neural Microexon" = "#fdb863")) +
  theme(legend.position = "none") +
  ylim(-35,15)
test <- wilcox.test(micro_neu_C2$`3ss_score`, avg_neu_C2$`3ss_score`) #0.323
test <- wilcox.test(micro_neu_C2$`3ss_score`, avg_nn_C2$`3ss_score`) #0.0898
test <- wilcox.test(micro_neu_C2$`3ss_score`, micro_nn_C2$`3ss_score`) #0.466
test <- wilcox.test(avg_neu_C2$`3ss_score`, micro_nn_C2$`3ss_score`) #0.0193
test <- wilcox.test(avg_nn_C2$`3ss_score`, micro_nn_C2$`3ss_score`) #0.000217
test <- wilcox.test(avg_nn_C2$`3ss_score`, avg_neu_C2$`3ss_score`) #0.231
