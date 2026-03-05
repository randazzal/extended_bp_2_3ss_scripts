library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(ggeasy)

neu_avg <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/SVM-BPfinder/neu_avg_intron_results.txt", sep = "\t")
nn_avg <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/SVM-BPfinder/nn_avg_intron_results.txt", sep = "\t")
neu_micro <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/SVM-BPfinder/neu_micro_intron_results.txt", sep = "\t")
nn_micro <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/SVM-BPfinder/nn_micro_intron_results.txt", sep = "\t")
nnC2_avg <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/SVM-BPfinder/neu_avgC2_results.txt", sep = "\t")
neuC2_avg <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/SVM-BPfinder/nn_avgC2_results.txt", sep = "\t")
nnC2_micro <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/SVM-BPfinder/neu_microC2_results.txt", sep = "\t")
neuC2_micro <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/SVM-BPfinder/nn_microC2_results.txt", sep = "\t")
controls <- read.delim("/data2/lackey_lab/randazza/microexons/Austin_data/SVM-BPfinder/AAH_control_results.txt", sep = "\t")

control_max <- controls %>%
  group_by(seq_id) %>%
  slice_max(order_by = svm_scr, n = 1) %>%
  ungroup()
control_max$group <- "AAH Controls"
neuavg_max <- neu_avg %>%
  group_by(seq_id) %>%
  slice_max(order_by = svm_scr, n = 1) %>%
  ungroup() #1042
neuavg_max$group <- "Neural Midexons"
nnavg_max <- nn_avg %>%
  group_by(seq_id) %>%
  slice_max(order_by = svm_scr, n = 1) %>%
  ungroup() #21599
nnavg_max$group <- "Non-Neural Midexons"
nnmicro_max <- nn_micro %>%
  group_by(seq_id) %>%
  slice_max(order_by = svm_scr, n = 1) %>%
  ungroup() #543
nnmicro_max$group <- "Non-Neural Microexons"
neumicro_max <- neu_micro %>%
  group_by(seq_id) %>%
  slice_max(order_by = svm_scr, n = 1) %>%
  ungroup() #270
neumicro_max$group <- "Neural Microexons"
neuavgC2_max <- neuC2_avg %>%
  group_by(seq_id) %>%
  slice_max(order_by = svm_scr, n = 1) %>%
  ungroup() #1042
neuavgC2_max$group <- "Neural Midexons C2"
nnavgC2_max <- nnC2_avg %>%
  group_by(seq_id) %>%
  slice_max(order_by = svm_scr, n = 1) %>%
  ungroup() #21599
nnavgC2_max$group <- "Non-Neural Midexons C2"
nnmicroC2_max <- nnC2_micro %>%
  group_by(seq_id) %>%
  slice_max(order_by = svm_scr, n = 1) %>%
  ungroup() #543
nnmicroC2_max$group <- "Non-Neural Microexons C2"
neumicroC2_max <- neuC2_micro %>%
  group_by(seq_id) %>%
  slice_max(order_by = svm_scr, n = 1) %>%
  ungroup() #270
neumicroC2_max$group <- "Neural Microexons C2"

combo <- rbind(neuavg_max, nnavg_max, nnmicro_max, neumicro_max)
mega <- rbind(neuavg_max, nnavg_max, nnmicro_max, neumicro_max, neuavgC2_max, nnavgC2_max, nnmicroC2_max, neumicroC2_max, control_max)

combo$group <- factor(combo$group, levels = c("Neural Microexons",
                                              "Neural Midexons",
                                              "Non-Neural Microexons",
                                              "Non-Neural Midexons"))
combo %>%
  ggplot(aes(x = group, y = ss_dist, fill = group)) +
  geom_boxplot(position = "dodge") +
  theme_minimal() +
  labs(title = "Branchpoint Distance from 3'ss",
       subtitle = "(highest scoring branchpoint)",
       x = "Exon Class", y = "Distance") +
  ggeasy::easy_center_title()+
  scale_fill_manual(values = c("Neural Midexons" = "#5e3c99",
                               "Non-Neural Microexons" = "#e66101",
                               "Non-Neural Midexons" = "#b2abd2",
                               "Neural Microexons" = "#fdb863")) +
  theme(legend.position = "none") +
  ylim(0,100)
test <- wilcox.test(neumicro_max$ss_dist, neuavg_max$ss_dist)#4.96e-07
test <- wilcox.test(neumicro_max$ss_dist, nnavg_max$ss_dist)#1.81e-11
test <- wilcox.test(neumicro_max$ss_dist, nnmicro_max$ss_dist)#0.0113
test <- wilcox.test(neuavg_max$ss_dist, nnmicro_max$ss_dist)#0.00246
test <- wilcox.test(nnavg_max$ss_dist, nnmicro_max$ss_dist)#7.99e-08
test <- wilcox.test(neuavg_max$ss_dist, nnavg_max$ss_dist)#0.0307

combo %>%
  ggplot(aes(x = group, y = ppt_len, fill = group)) +
  geom_boxplot(position = "dodge") +
  theme_minimal() +
  labs(title = "Polypyrimidine Tract Length",
       subtitle = "(highest scoring branchpoint)",
       x = "Exon Class", y = "Length") +
  ggeasy::easy_center_title()+
  scale_fill_manual(values = c("Neural Midexons" = "#5e3c99",
                               "Non-Neural Microexons" = "#e66101",
                               "Non-Neural Midexons" = "#b2abd2",
                               "Neural Microexons" = "#fdb863")) +
  theme(legend.position = "none") +
  ylim(0,100)
test <- wilcox.test(neumicro_max$ppt_len, neuavg_max$ppt_len)#8.97e-23
test <- wilcox.test(neumicro_max$ppt_len, nnavg_max$ppt_len)#1.98e-36
test <- wilcox.test(neumicro_max$ppt_len, nnmicro_max$ppt_len)#3.63e-09
test <- wilcox.test(neuavg_max$ppt_len, nnmicro_max$ppt_len)#1.3e-07
test <- wilcox.test(nnavg_max$ppt_len, nnmicro_max$ppt_len)#4.94e-19
test <- wilcox.test(neuavg_max$ppt_len, nnavg_max$ppt_len)#0.00114

combo$log <- log2(combo$ppt_off)
combo %>%
  ggplot(aes(x = group, y = log, fill = group)) +
  geom_boxplot(position = "dodge") +
  theme_minimal() +
  labs(title = "Distance between BP and PPT",
       subtitle = "(highest scoring branchpoint)",
       x = "Exon Class", y = "log2(Distance)") +
  ggeasy::easy_center_title()+
  scale_fill_manual(values = c("Neural Midexons" = "#5e3c99",
                               "Non-Neural Microexons" = "#e66101",
                               "Non-Neural Midexons" = "#b2abd2",
                               "Neural Microexons" = "#fdb863")) +
  theme(legend.position = "none") +
  ylim(0,7)
test <- wilcox.test(neumicro_max$ppt_off, neuavg_max$ppt_off)#3.56e-06
test <- wilcox.test(neumicro_max$ppt_off, nnavg_max$ppt_off)#4.13e-07
test <- wilcox.test(neumicro_max$ppt_off, nnmicro_max$ppt_off)#0.00146
test <- wilcox.test(neuavg_max$ppt_off, nnmicro_max$ppt_off)#0.125
test <- wilcox.test(nnavg_max$ppt_off, nnmicro_max$ppt_off)#0.0791
test <- wilcox.test(neuavg_max$ppt_off, nnavg_max$ppt_off)#0.885

mega$group <- factor(mega$group, levels = c("Neural Microexons",
                                            "Neural Midexons",
                                            "Non-Neural Microexons",
                                            "Non-Neural Midexons",
                                            "Neural Microexons C2",
                                            "Neural Midexons C2",
                                            "Non-Neural Microexons C2",
                                            "Non-Neural Midexons C2",
                                            "AAH Controls"))
mega %>%
  ggplot(aes(x = group, y = ss_dist, fill = group)) +
  geom_boxplot(position = "dodge") +
  theme_minimal() +
  labs(title = "Branchpoint Distance from 3'ss",
       subtitle = "(highest scoring branchpoint)",
       x = "Exon Class", y = "Distance") +
  ggeasy::easy_center_title() +
  theme(legend.position = "none") +
  ylim(0,100)
test <- wilcox.test(neumicro_max$ss_dist, neuavgC2_max$ss_dist)#5.59e-16
test <- wilcox.test(neumicro_max$ss_dist, nnavgC2_max$ss_dist)#1.93e-14
test <- wilcox.test(neumicro_max$ss_dist, nnmicroC2_max$ss_dist)#1.58e-07
test <- wilcox.test(neumicro_max$ss_dist, neumicroC2_max$ss_dist)#2.76e-09
test <- wilcox.test(neumicro_max$ss_dist, control_max$ss_dist)#0.000602

mega %>%
  ggplot(aes(x = group, y = ppt_len, fill = group)) +
  geom_boxplot(position = "dodge") +
  theme_minimal() +
  labs(title = "Polypyrimidine Tract Length",
       subtitle = "(highest scoring branchpoint)",
       x = "Exon Class", y = "Length") +
  ggeasy::easy_center_title()+
  theme(legend.position = "none") +
  ylim(0,100)
test <- wilcox.test(neumicro_max$ppt_len, neuavgC2_max$ppt_len)#5.08e-32
test <- wilcox.test(neumicro_max$ppt_len, nnavgC2_max$ppt_len)#2.74e-24
test <- wilcox.test(neumicro_max$ppt_len, nnmicroC2_max$ppt_len)#4.2e-15
test <- wilcox.test(neumicro_max$ppt_len, neumicroC2_max$ppt_len)#2.45e-22
test <- wilcox.test(neumicro_max$ppt_len, control_max$ppt_len)#0.0112

mega$log <- log2(mega$ppt_off)
mega %>%
  ggplot(aes(x = group, y = log, fill = group)) +
  geom_boxplot(position = "dodge") +
  theme_minimal() +
  labs(title = "Distance between BP and PPT",
       subtitle = "(highest scoring branchpoint)",
       x = "Exon Class", y = "log2(Distance)") +
  ggeasy::easy_center_title()+
  theme(legend.position = "none") +
  ylim(0,7)
test <- wilcox.test(neumicro_max$ppt_off, neuavgC2_max$ppt_off)#9.26e-05
test <- wilcox.test(neumicro_max$ppt_off, nnavgC2_max$ppt_off)#8.47e-05
test <- wilcox.test(neumicro_max$ppt_off, nnmicroC2_max$ppt_off)#0.178
test <- wilcox.test(neumicro_max$ppt_off, neumicroC2_max$ppt_off)#0.00753
test <- wilcox.test(neumicro_max$ppt_off, control_max$ppt_off)#0.26
