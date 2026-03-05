library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(ggeasy)

### chicken mid vs micro exon predicted pairing probabilities (450 is start of exon)
mega_mic <- read.delim("/data2/lackey_lab/randazza/microexons/devo/chicken/fixed_beds/microA3_bp/mega_gga_micro.tsv", sep = "\t")
mega_mic$Position <- as.integer(mega_mic$Position)

micro_noUGC <- read.delim("/data2/lackey_lab/randazza/microexons/motifs/UGC_motif_chick/DIY/gga_vastDB_micro_noUGC.txt", sep = "\t")
micro_UGC <- read.delim("/data2/lackey_lab/randazza/microexons/motifs/UGC_motif_chick/DIY/gga_vastDB_micro_UGC.txt", sep = "\t")
mid_noUGC <- read.delim("/data2/lackey_lab/randazza/microexons/motifs/UGC_motif_chick/DIY/gga_vastDB_mid_noUGC.txt", sep = "\t")
mid_UGC <- read.delim("/data2/lackey_lab/randazza/microexons/motifs/UGC_motif_chick/DIY/gga_vastDB_mid_UGC.txt", sep = "\t")

#### filter for neural
micro_noUGC <- micro_noUGC %>%
  filter(tissue %in% c("fdr_l", "fdr_g"))
micro_UGC <- micro_UGC %>%
  filter(tissue %in% c("fdr_l", "fdr_g"))
mid_noUGC <- mid_noUGC %>%
  filter(tissue %in% c("fdr_l", "fdr_g"))
mid_UGC <- mid_UGC %>%
  filter(tissue %in% c("fdr_l", "fdr_g"))

## divide into upstream, bp 2 3ss, exon, downstream

mega_micro_noUGC <- mega_mic %>%
  left_join(micro_noUGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    branchpoint_pos = 450 - ss_dist
  ) %>%
  filter(
    Position >= branchpoint_pos,
    Position < 450
  )
mega_micro_UGC <- mega_mic %>%
  left_join(micro_UGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    branchpoint_pos = 450 - ss_dist
  ) %>%
  filter(
    Position >= branchpoint_pos,
    Position < 450
  )
exon_micro_noUGC <- mega_mic %>%
  left_join(micro_noUGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    exon_pos = 450 + LENGTH - 1
  ) %>%
  filter(
    Position >= 450,
    Position <= exon_pos
  )
exon_micro_UGC <- mega_mic %>%
  left_join(micro_UGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    exon_pos = 450 + LENGTH - 1
  ) %>%
  filter(
    Position >= 450,
    Position <= exon_pos
  )

up_micro_noUGC <- mega_mic %>%
  left_join(micro_noUGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    branchpoint_pos = 450 - ss_dist
  ) %>%
  filter(
    Position >= branchpoint_pos - 100,
    Position < branchpoint_pos
  )
up_micro_UGC <- mega_mic %>%
  left_join(micro_UGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    branchpoint_pos = 450 - ss_dist
  ) %>%
  filter(
    Position >= branchpoint_pos - 100,
    Position < branchpoint_pos
  )

down_micro_noUGC <- mega_mic %>%
  left_join(micro_noUGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    exon_pos = 450 + LENGTH - 1
  ) %>%
  filter(
    Position > exon_pos,
    Position <= exon_pos + 100
  )
down_micro_UGC <- mega_mic %>%
  left_join(micro_UGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    exon_pos = 450 + LENGTH - 1
  ) %>%
  filter(
    Position > exon_pos,
    Position <= exon_pos + 100
  )

rm(mega_mic)
rm(micro_noUGC)
rm(micro_UGC)

mega_mic <- read.delim("/data2/lackey_lab/randazza/microexons/devo/chicken/fixed_beds/avgA3_bp/mega_gga_mid.tsv", sep = "\t")
mega_mic$Position <- as.integer(mega_mic$Position)
mega_mid_noUGC <- mega_mic %>%
  left_join(mid_noUGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    branchpoint_pos = 450 - ss_dist
  ) %>%
  filter(
    Position >= branchpoint_pos,
    Position < 450
  )
mega_mid_UGC <- mega_mic %>%
  left_join(mid_UGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    branchpoint_pos = 450 - ss_dist
  ) %>%
  filter(
    Position >= branchpoint_pos,
    Position < 450
  )

exon_mid_noUGC <- mega_mic %>%
  left_join(mid_noUGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    exon_pos = 450 + LENGTH - 1
  ) %>%
  filter(
    Position >= 450,
    Position <= exon_pos
  )
exon_mid_UGC <- mega_mic %>%
  left_join(mid_UGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    exon_pos = 450 + LENGTH - 1
  ) %>%
  filter(
    Position >= 450,
    Position <= exon_pos
  )

up_mid_noUGC <- mega_mic %>%
  left_join(mid_noUGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    branchpoint_pos = 450 - ss_dist
  ) %>%
  filter(
    Position >= branchpoint_pos - 100,
    Position < branchpoint_pos
  )
up_mid_UGC <- mega_mic %>%
  left_join(mid_UGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    branchpoint_pos = 450 - ss_dist
  ) %>%
  filter(
    Position >= branchpoint_pos - 100,
    Position < branchpoint_pos
  )

down_mid_noUGC <- mega_mic %>%
  left_join(mid_noUGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    exon_pos = 450 + LENGTH - 1
  ) %>%
  filter(
    Position > exon_pos,
    Position <= exon_pos + 100
  )
down_mid_UGC <- mega_mic %>%
  left_join(mid_UGC, by = c("DP_File" = "EVENT.x")) %>%
  mutate(
    exon_pos = 450 + LENGTH - 1
  ) %>%
  filter(
    Position > exon_pos,
    Position <= exon_pos + 100
  )
rm(mega_mic)
rm(mid_noUGC)
rm(mid_UGC)

down_micro_noUGC <- down_micro_noUGC[, c("Position", "SumPairingProbability...", "DP_File")]
down_micro_UGC <- down_micro_UGC[, c("Position", "SumPairingProbability...", "DP_File")]
down_mid_noUGC <- down_mid_noUGC[, c("Position", "SumPairingProbability...", "DP_File")]
down_mid_UGC <- down_mid_UGC[, c("Position", "SumPairingProbability...", "DP_File")]
down_micro_noUGC$region <- "Downstream"
down_micro_UGC$region <- "Downstream"
down_mid_noUGC$region <- "Downstream"
down_mid_UGC$region <- "Downstream"

up_micro_noUGC <- up_micro_noUGC[, c("Position", "SumPairingProbability...", "DP_File")]
up_micro_UGC <- up_micro_UGC[, c("Position", "SumPairingProbability...", "DP_File")]
up_mid_noUGC <- up_mid_noUGC[, c("Position", "SumPairingProbability...", "DP_File")]
up_mid_UGC <- up_mid_UGC[, c("Position", "SumPairingProbability...", "DP_File")]
up_micro_noUGC$region <- "Upstream"
up_micro_UGC$region <- "Upstream"
up_mid_noUGC$region <- "Upstream"
up_mid_UGC$region <- "Upstream"

exon_micro_noUGC <- exon_micro_noUGC[, c("Position", "SumPairingProbability...", "DP_File")]
exon_micro_UGC <- exon_micro_UGC[, c("Position", "SumPairingProbability...", "DP_File")]
exon_mid_noUGC <- exon_mid_noUGC[, c("Position", "SumPairingProbability...", "DP_File")]
exon_mid_UGC <- exon_mid_UGC[, c("Position", "SumPairingProbability...", "DP_File")]
exon_micro_noUGC$region <- "Exon"
exon_micro_UGC$region <- "Exon"
exon_mid_noUGC$region <- "Exon"
exon_mid_UGC$region <- "Exon"

mega_micro_noUGC <- mega_micro_noUGC[, c("Position", "SumPairingProbability...", "DP_File")]
mega_micro_UGC <- mega_micro_UGC[, c("Position", "SumPairingProbability...", "DP_File")]
mega_mid_noUGC <- mega_mid_noUGC[, c("Position", "SumPairingProbability...", "DP_File")]
mega_mid_UGC <- mega_mid_UGC[, c("Position", "SumPairingProbability...", "DP_File")]
mega_micro_noUGC$region <- "BP2SS"
mega_micro_UGC$region <- "BP2SS"
mega_mid_noUGC$region <- "BP2SS"
mega_mid_UGC$region <- "BP2SS"

mic_UGC <- rbind(mega_micro_UGC, down_micro_UGC, up_micro_UGC, exon_micro_UGC)
mic_noUGC <- rbind(mega_micro_noUGC, down_micro_noUGC, up_micro_noUGC, exon_micro_noUGC)
mid_noUGC <- rbind(mega_mid_noUGC, down_mid_noUGC, up_mid_noUGC, exon_mid_noUGC)
mid_UGC <- rbind(mega_mid_UGC, down_mid_UGC, up_mid_UGC, exon_mid_UGC)

mic_noUGC$group <- "Micro No UGC"
mic_noUGC$class <- "Microexon"
mic_UGC$group <- "Micro UGC"
mic_UGC$class <- "Microexon"
mid_noUGC$group <- "Mid No UGC"
mid_noUGC$class <- "Midexon"
mid_UGC$group <- "Mid UGC"
mid_UGC$class <- "Midexon"

mega_mega <- rbind(mic_noUGC, mic_UGC, mid_noUGC, mid_UGC)
mega_mega$SumPairingProbability... <- as.integer(mega_mega$SumPairingProbability...)

mean_prob_by_event <- mega_mega %>%
  group_by(DP_File, region) %>%
  summarise(
    mean_SumProbability = mean(SumPairingProbability..., na.rm = TRUE)
  ) %>%
  ungroup()

mean_prob_by_event <- merge(mean_prob_by_event, mega_mega, by = c("DP_File", "region"))
mean_prob_by_event$region <- factor(mean_prob_by_event$region, levels = c("Upstream", "BP2SS", "Exon", "Downstream"))

mean_prob_by_event %>%
  ggplot(aes(x = group, y = mean_SumProbability)) +
  geom_boxplot() +
  theme_minimal()
mean_prob_by_event %>%
  ggplot(aes(x = region, y = mean_SumProbability, fill = group)) +
  geom_boxplot() +
  theme_minimal()
mean_prob_by_event %>%
  ggplot(aes(x = region, y = SumPairingProbability..., fill = group)) +
  geom_boxplot() +
  theme_minimal()
mean_prob_by_event %>%
  ggplot(aes(x = region, y = mean_SumProbability, fill = class)) +
  geom_boxplot() +
  theme_minimal()
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "Upstream",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "Upstream",]$mean_SumProbability,) #0.0982 n 0.223
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "Upstream",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "Upstream",]$mean_SumProbability,) #1.42e-141 n 0.152
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "Upstream",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "Upstream",]$mean_SumProbability,) #1.47e-142 n 5.3e-23
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "Upstream",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "Upstream",]$mean_SumProbability,) #0 n 2.46e-16

test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "Upstream",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "Upstream",]$SumPairingProbability...,) #0.607
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "Upstream",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "Upstream",]$SumPairingProbability...,) #0.599
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "Upstream",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "Upstream",]$SumPairingProbability...,) #0.013
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "Upstream",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "Upstream",]$SumPairingProbability...,) #0.0188

test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "Downstream",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "Downstream",]$mean_SumProbability,) #1.58e-51 n 0.0658
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "Downstream",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "Downstream",]$mean_SumProbability,) #1.68e-54 n 2.96e-9
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "Downstream",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "Downstream",]$mean_SumProbability,) #4.25e-49 n 1.65e-31
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "Downstream",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "Downstream",]$mean_SumProbability,) #0 n 5.18e-67

test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "Downstream",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "Downstream",]$SumPairingProbability...,) #0.197
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "Downstream",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "Downstream",]$SumPairingProbability...,) #0.0605
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "Downstream",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "Downstream",]$SumPairingProbability...,) #0.00559
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "Downstream",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "Downstream",]$SumPairingProbability...,) #3.66e-5

test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "Exon",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "Exon",]$mean_SumProbability,) #0.0694 n 0.0282
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "Exon",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "Exon",]$mean_SumProbability,) #0 n 2.97e-135
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "Exon",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "Exon",]$mean_SumProbability,) #0.258 n 0.422
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "Exon",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "Exon",]$mean_SumProbability,) #9.48e-9 n 0.383

test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "Exon",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "Exon",]$SumPairingProbability...,) #0.222
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "Exon",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "Exon",]$SumPairingProbability...,) #0.000177
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "Exon",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "Exon",]$SumPairingProbability...,) #0.317
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "Exon",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "Exon",]$SumPairingProbability...,) #0.263

test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "BP2SS",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "BP2SS",]$mean_SumProbability,) #0.000607 n 0.262
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "BP2SS",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "BP2SS",]$mean_SumProbability,) #2.47e-36 n 0.0766
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "BP2SS",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "BP2SS",]$mean_SumProbability,) #8.16e-47 n 3.66e-37
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "BP2SS",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "BP2SS",]$mean_SumProbability,) #0 n 2.55e-198

test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "BP2SS",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "BP2SS",]$SumPairingProbability...,) #0.962
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "BP2SS",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "BP2SS",]$SumPairingProbability...,) #0.188
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC" & mean_prob_by_event$region == "BP2SS",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC" & mean_prob_by_event$region == "BP2SS",]$SumPairingProbability...,) #1.92e-5
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro UGC" & mean_prob_by_event$region == "BP2SS",]$SumPairingProbability...,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC" & mean_prob_by_event$region == "BP2SS",]$SumPairingProbability...,) #1.31e-20

test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Micro UGC",]$mean_SumProbability,) #0.000607
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC",]$mean_SumProbability,) #2.47e-36
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC",]$mean_SumProbability,) #8.16e-47
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro No UGC",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC",]$mean_SumProbability,) #1.13e-89
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro UGC",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid No UGC",]$mean_SumProbability,) #3.91e-260
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$group == "Micro UGC",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$group == "Mid UGC",]$mean_SumProbability,) #0

test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$class == "Microexon" & mean_prob_by_event$region == "BP2SS",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$class == "Midexon" & mean_prob_by_event$region == "BP2SS",]$mean_SumProbability,) #0 n 6.74e-238
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$class == "Microexon" & mean_prob_by_event$region == "Upstream",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$class == "Midexon" & mean_prob_by_event$region == "Upstream",]$mean_SumProbability,) #0 n 1.03e-37
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$class == "Microexon" & mean_prob_by_event$region == "Downstream",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$class == "Midexon" & mean_prob_by_event$region == "Downstream",]$mean_SumProbability,) #0 n 9.76e-85
test <- wilcox.test(mean_prob_by_event[mean_prob_by_event$class == "Microexon" & mean_prob_by_event$region == "Exon",]$mean_SumProbability,
                    mean_prob_by_event[mean_prob_by_event$class == "Midexon" & mean_prob_by_event$region == "Exon",]$mean_SumProbability,) #4.76e-07 n 0.867