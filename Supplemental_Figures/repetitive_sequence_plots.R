library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)

setwd("/data2/lackey_lab/randazza/microexons/repeats/")
hsa_micro_up <- read.delim("hsa_micro_up_parsed.tsv")
hsa_micro_up <- hsa_micro_up %>%
  separate(data, into = c("start", "end", "size", "copy_number", "consensus_size", "matches%", "indels%", "score",
                          "As", "Cs", "Gs", "Ts", "entropy", "motif", "aligned_repeat", "consensus"), 
           sep = " ", remove = TRUE)
count <- unique(hsa_micro_up$intron) #336

hsa_micro_down <- read.delim("hsa_micro_down_parsed.tsv")
hsa_micro_down <- hsa_micro_down %>%
  separate(data, into = c("start", "end", "size", "copy_number", "consensus_size", "matches%", "indels%", "score",
                          "As", "Cs", "Gs", "Ts", "entropy", "motif", "aligned_repeat", "consensus"), 
           sep = " ", remove = TRUE)
count <- unique(hsa_micro_down$intron) #305

hsa_avg_up <- read.delim("hsa_avg_up_parsed.tsv")
hsa_avg_up <- hsa_avg_up %>%
  separate(data, into = c("start", "end", "size", "copy_number", "consensus_size", "matches%", "indels%", "score",
                          "As", "Cs", "Gs", "Ts", "entropy", "motif", "aligned_repeat", "consensus"), 
           sep = " ", remove = TRUE)
count <- unique(hsa_avg_up$intron) #10702

hsa_avg_down <- read.delim("hsa_avg_down_parsed.tsv")
hsa_avg_down <- hsa_avg_down %>%
  separate(data, into = c("start", "end", "size", "copy_number", "consensus_size", "matches%", "indels%", "score",
                          "As", "Cs", "Gs", "Ts", "entropy", "motif", "aligned_repeat", "consensus"), 
           sep = " ", remove = TRUE)
count <- unique(hsa_avg_down$intron) #10900

gga_micro_up <- read.delim("gga_micro_up_parsed.tsv")
gga_micro_up <- gga_micro_up %>%
  separate(data, into = c("start", "end", "size", "copy_number", "consensus_size", "matches%", "indels%", "score",
                          "As", "Cs", "Gs", "Ts", "entropy", "motif", "aligned_repeat", "consensus"), 
           sep = " ", remove = TRUE)
count <- unique(gga_micro_up$intron) #126

gga_micro_down <- read.delim("gga_micro_down_parsed.tsv")
gga_micro_down <- gga_micro_down %>%
  separate(data, into = c("start", "end", "size", "copy_number", "consensus_size", "matches%", "indels%", "score",
                          "As", "Cs", "Gs", "Ts", "entropy", "motif", "aligned_repeat", "consensus"), 
           sep = " ", remove = TRUE)
count <- unique(gga_micro_down$intron) #106

gga_avg_up <- read.delim("gga_avg_up_parsed.tsv")
gga_avg_up <- gga_avg_up %>%
  separate(data, into = c("start", "end", "size", "copy_number", "consensus_size", "matches%", "indels%", "score",
                          "As", "Cs", "Gs", "Ts", "entropy", "motif", "aligned_repeat", "consensus"), 
           sep = " ", remove = TRUE)
count <- unique(gga_avg_up$intron) #1876

gga_avg_down <- read.delim("gga_avg_down_parsed.tsv")
gga_avg_down <- gga_avg_down %>%
  separate(data, into = c("start", "end", "size", "copy_number", "consensus_size", "matches%", "indels%", "score",
                          "As", "Cs", "Gs", "Ts", "entropy", "motif", "aligned_repeat", "consensus"), 
           sep = " ", remove = TRUE)
count <- unique(gga_avg_down$intron) #1739

###plotting
convert_and_add_total <- function(df) {
  df$size <- as.numeric(df$size)
  df$copy_number <- as.numeric(df$copy_number)
  
  df$total <- df$size * df$copy_number
  
  return(df)
}
hsa_micro_down <- convert_and_add_total(hsa_micro_down)
hsa_micro_up <- convert_and_add_total(hsa_micro_up)
hsa_avg_down <- convert_and_add_total(hsa_avg_down)
hsa_avg_up <- convert_and_add_total(hsa_avg_up)
gga_micro_down <- convert_and_add_total(gga_micro_down)
gga_micro_up <- convert_and_add_total(gga_micro_up)
gga_avg_down <- convert_and_add_total(gga_avg_down)
gga_avg_up <- convert_and_add_total(gga_avg_up)

sum_by_intron <- function(df) {
  df %>%
    group_by(intron) %>%
    summarize(total_sum = sum(total, na.rm = TRUE))
}
hsa_micro_down_sum <- sum_by_intron(hsa_micro_down)
hsa_micro_up_sum <- sum_by_intron(hsa_micro_up)
hsa_avg_down_sum <- sum_by_intron(hsa_avg_down)
hsa_avg_up_sum <- sum_by_intron(hsa_avg_up)
gga_micro_down_sum <- sum_by_intron(gga_micro_down)
gga_micro_up_sum <- sum_by_intron(gga_micro_up)
gga_avg_down_sum <- sum_by_intron(gga_avg_down)
gga_avg_up_sum <- sum_by_intron(gga_avg_up)

hsa_micro_down_sum$species <- "Hsa Microexon"
hsa_micro_down_sum$group <- "micro"
hsa_micro_up_sum$species <- "Hsa Microexon"
hsa_micro_up_sum$group <- "micro"

gga_micro_down_sum$species <- "Gga Microexon"
gga_micro_down_sum$group <- "micro"
gga_micro_up_sum$species <- "Gga Microexon"
gga_micro_up_sum$group <- "micro"

hsa_avg_down_sum$species <- "Hsa Midexon"
hsa_avg_down_sum$group <- "mid"
hsa_avg_up_sum$species <- "Hsa Midexon"
hsa_avg_up_sum$group <- "mid"

gga_avg_down_sum$species <- "Gga Midexon"
gga_avg_down_sum$group <- "mid"
gga_avg_up_sum$species <- "Gga Midexon"
gga_avg_up_sum$group <- "mid"

combo_up <- rbind(hsa_micro_up_sum, hsa_avg_up_sum, gga_micro_up_sum, gga_avg_up_sum)
combo_up$species <- factor(combo_up$species, levels = c("Hsa Microexon","Hsa Midexon", "Gga Microexon", "Gga Midexon"))
combo_down <- rbind(hsa_micro_down_sum, hsa_avg_down_sum, gga_micro_down_sum, gga_avg_down_sum)
combo_down$species <- factor(combo_down$species, levels = c("Hsa Microexon","Hsa Midexon", "Gga Microexon", "Gga Midexon"))

combo_up %>%
  ggplot(aes(x = species, y = total_sum, fill = group)) +
  geom_boxplot() #+
  #ylim(0,2000)
combo_up$log <- log10(combo_up$total_sum)
combo_up %>%
  ggplot(aes(x = species, y = log, fill = group)) +
  geom_boxplot() +
  theme_minimal()
test <- wilcox.test(hsa_micro_up_sum$total_sum, hsa_avg_up_sum$total_sum) #0.00109
test <- wilcox.test(hsa_micro_up_sum$total_sum, gga_avg_up_sum$total_sum) #6.1e-18
test <- wilcox.test(hsa_micro_up_sum$total_sum, gga_micro_up_sum$total_sum) #2.09e-07
test <- wilcox.test(gga_micro_up_sum$total_sum, gga_avg_up_sum$total_sum) #0.693

combo_down %>%
  ggplot(aes(x = species, y = total_sum, fill = group)) +
  geom_boxplot() #+
  #ylim(0,2000)
combo_up$log <- log10(combo_up$total_sum)
combo_up %>%
  ggplot(aes(x = species, y = log, fill = group)) +
  geom_boxplot()
test <- wilcox.test(hsa_micro_down_sum$total_sum, hsa_avg_down_sum$total_sum) #0.000129
test <- wilcox.test(hsa_micro_down_sum$total_sum, gga_avg_down_sum$total_sum) #5.08e-12
test <- wilcox.test(hsa_micro_down_sum$total_sum, gga_micro_down_sum$total_sum) #4.64e-06
test <- wilcox.test(gga_micro_down_sum$total_sum, gga_avg_down_sum$total_sum) #0.202

#split into neural and non-neural
hsa_micro_up_sum$EVENT <- sub("^@", "", sub("::.*", "", hsa_micro_up_sum$intron))
hsa_micro_down_sum$EVENT <- sub("^@", "", sub("::.*", "", hsa_micro_down_sum$intron))
hsa_avg_up_sum$EVENT <- sub("^@", "", sub("::.*", "", hsa_avg_up_sum$intron))
hsa_avg_down_sum$EVENT <- sub("^@", "", sub("::.*", "", hsa_avg_down_sum$intron))
gga_micro_up_sum$EVENT <- sub("^@", "", sub("::.*", "", gga_micro_up_sum$intron))
gga_micro_down_sum$EVENT <- sub("^@", "", sub("::.*", "", gga_micro_down_sum$intron))
gga_avg_up_sum$EVENT <- sub("^@", "", sub("::.*", "", gga_avg_up_sum$intron))
gga_avg_down_sum$EVENT <- sub("^@", "", sub("::.*", "", gga_avg_down_sum$intron))

micro_hsa <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/detailed_list.txt", sep = "\t")
hsa_neuro_micro <- micro_hsa %>%
  filter(tissue %in% c("fdr_g", "fdr_l")) %>%
  pull(EVENT)
hsa_nn_micro <- micro_hsa %>%
  filter(tissue %in% c("common", "rare", "other")) %>%
  pull(EVENT)
avg_hsa <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/avg_regroup.txt", sep = "\t")
hsa_neuro_avg <- avg_hsa %>%
  filter(tissue %in% c("fdr_g", "fdr_l")) %>%
  pull(EVENT)
hsa_nn_avg <- avg_hsa %>%
  filter(tissue %in% c("common", "rare", "other")) %>%
  pull(EVENT)
micro_gga <- read.delim("/data2/lackey_lab/randazza/microexons/devo/chicken/gga_classes.txt", sep = "\t")
gga_neuro_micro <- micro_gga %>%
  filter(tissue %in% c("fdr_g", "fdr_l")) %>%
  pull(EVENT)
gga_nn_micro <- micro_gga %>%
  filter(tissue %in% c("common", "rare", "other")) %>%
  pull(EVENT)
avg_gga <- read.delim("/data2/lackey_lab/randazza/microexons/devo/chicken/avg_control/gga_avg_classes.txt", sep = "\t")
gga_neuro_avg <- avg_gga %>%
  filter(tissue %in% c("fdr_g", "fdr_l")) %>%
  pull(EVENT)
gga_nn_avg <- avg_gga %>%
  filter(tissue %in% c("common", "rare", "other")) %>%
  pull(EVENT)

hsa_micro_up_sum_neuro <- hsa_micro_up_sum %>%
  filter(EVENT %in% hsa_neuro_micro) %>%
  mutate(species = "Hsa Neural Microexon")
hsa_micro_up_sum_nn <- hsa_micro_up_sum %>%
  filter(EVENT %in% hsa_nn_micro) %>%
  mutate(species = "Hsa Non-Neural Microexon")
hsa_micro_down_sum_neuro <- hsa_micro_down_sum %>%
  filter(EVENT %in% hsa_neuro_micro) %>%
  mutate(species = "Hsa Neural Microexon")
hsa_micro_down_sum_nn <- hsa_micro_down_sum %>%
  filter(EVENT %in% hsa_nn_micro) %>%
  mutate(species = "Hsa Non-Neural Microexon")

hsa_avg_up_sum_neuro <- hsa_avg_up_sum %>%
  filter(EVENT %in% hsa_neuro_avg) %>%
  mutate(species = "Hsa Neural Midexon")
hsa_avg_up_sum_nn <- hsa_avg_up_sum %>%
  filter(EVENT %in% hsa_nn_avg) %>%
  mutate(species = "Hsa Non-Neural Midexon")
hsa_avg_down_sum_neuro <- hsa_avg_down_sum %>%
  filter(EVENT %in% hsa_neuro_avg) %>%
  mutate(species = "Hsa Neural Midexon")
hsa_avg_down_sum_nn <- hsa_avg_down_sum %>%
  filter(EVENT %in% hsa_nn_avg) %>%
  mutate(species = "Hsa Non-Neural Midexon")

gga_micro_up_sum_neuro <- gga_micro_up_sum %>%
  filter(EVENT %in% gga_neuro_micro) %>%
  mutate(species = "Gga Neural Microexon")
gga_micro_up_sum_nn <- gga_micro_up_sum %>%
  filter(EVENT %in% gga_nn_micro) %>%
  mutate(species = "Gga Non-Neural Microexon")
gga_micro_down_sum_neuro <- gga_micro_down_sum %>%
  filter(EVENT %in% gga_neuro_micro) %>%
  mutate(species = "Gga Neural Microexon")
gga_micro_down_sum_nn <- gga_micro_down_sum %>%
  filter(EVENT %in% gga_nn_micro) %>%
  mutate(species = "Gga Non-Neural Microexon")

gga_avg_up_sum_neuro <- gga_avg_up_sum %>%
  filter(EVENT %in% gga_neuro_avg) %>%
  mutate(species = "Gga Neural Midexon")
gga_avg_up_sum_nn <- gga_avg_up_sum %>%
  filter(EVENT %in% gga_nn_avg) %>%
  mutate(species = "Gga Non-Neural Midexon")
gga_avg_down_sum_neuro <- gga_avg_down_sum %>%
  filter(EVENT %in% gga_neuro_avg) %>%
  mutate(species = "Gga Neural Midexon")
gga_avg_down_sum_nn <- gga_avg_down_sum %>%
  filter(EVENT %in% gga_nn_avg) %>%
  mutate(species = "Gga Non-Neural Midexon")

combo_up <- rbind(hsa_micro_up_sum_neuro, hsa_avg_up_sum_neuro, gga_micro_up_sum_neuro, gga_avg_up_sum_neuro,
                  hsa_micro_up_sum_nn, hsa_avg_up_sum_nn, gga_micro_up_sum_nn, gga_avg_up_sum_nn)
combo_up$species <- factor(combo_up$species, levels = c("Hsa Neural Microexon","Gga Neural Microexon",
                                                        "Hsa Neural Midexon", "Gga Neural Midexon",
                                                        "Hsa Non-Neural Microexon","Gga Non-Neural Microexon",
                                                        "Hsa Non-Neural Midexon", "Gga Non-Neural Midexon"))
combo_down <- rbind(hsa_micro_down_sum_neuro, hsa_avg_down_sum_neuro, gga_micro_down_sum_neuro, gga_avg_down_sum_neuro,
                    hsa_micro_down_sum_nn, hsa_avg_down_sum_nn, gga_micro_down_sum_nn, gga_avg_down_sum_nn)
combo_down$species <- factor(combo_down$species, levels = c("Hsa Neural Microexon","Gga Neural Microexon",
                                                            "Hsa Neural Midexon", "Gga Neural Midexon",
                                                            "Hsa Non-Neural Microexon","Gga Non-Neural Microexon",
                                                            "Hsa Non-Neural Midexon", "Gga Non-Neural Midexon"))

combo_up %>%
  ggplot(aes(x = species, y = total_sum, fill = group)) +
  geom_boxplot() +
  #ylim(0,2000) +
  theme(legend.position = "none") +
  theme_minimal()
combo_up$log <- log10(combo_up$total_sum)
combo_up %>%
  ggplot(aes(x = species, y = log, fill = group)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  theme_minimal()
test <- wilcox.test(hsa_micro_up_sum_neuro$total_sum, gga_micro_up_sum_neuro$total_sum) #241.7481 98.00385 0.0071
test <- wilcox.test(hsa_avg_up_sum_neuro$total_sum, gga_avg_up_sum_neuro$total_sum) #614.0179 150.745 1.04e-21
test <- wilcox.test(hsa_micro_up_sum_nn$total_sum, gga_micro_up_sum_nn$total_sum) #581.2061 180.1 2.55e-05
test <- wilcox.test(hsa_avg_up_sum_nn$total_sum, gga_avg_up_sum_nn$total_sum) #609.8172 143.8317 5.44e-133
test <- wilcox.test(hsa_micro_up_sum_neuro$total_sum, hsa_avg_up_sum_neuro$total_sum) #0.0018
test <- wilcox.test(hsa_micro_up_sum_neuro$total_sum, hsa_avg_up_sum_nn$total_sum) #0.000317
test <- wilcox.test(hsa_micro_up_sum_neuro$total_sum, hsa_micro_up_sum_nn$total_sum) #0.0348
test <- wilcox.test(hsa_micro_up_sum_nn$total_sum, hsa_avg_up_sum_nn$total_sum) #0.122
test <- wilcox.test(gga_micro_up_sum_neuro$total_sum, gga_avg_up_sum_neuro$total_sum) #0.852
test <- wilcox.test(gga_micro_up_sum_neuro$total_sum, gga_avg_up_sum_nn$total_sum) #0.495
test <- wilcox.test(gga_micro_up_sum_neuro$total_sum, gga_micro_up_sum_nn$total_sum) #0.641
test <- wilcox.test(gga_micro_up_sum_nn$total_sum, gga_avg_up_sum_nn$total_sum) #0.915

combo_down %>%
  ggplot(aes(x = species, y = total_sum, fill = group)) +
  geom_boxplot() +
  ylim(0,2000) +
  theme(legend.position = "none") +
  theme_minimal()
combo_down$log <- log10(combo_down$total_sum)
combo_down %>%
  ggplot(aes(x = species, y = log, fill = group)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  theme_minimal()
test <- wilcox.test(hsa_micro_down_sum_neuro$total_sum, gga_micro_down_sum_neuro$total_sum) #338.6698 147.0347 0.0106
test <- wilcox.test(hsa_avg_down_sum_neuro$total_sum, gga_avg_down_sum_neuro$total_sum) #506.1063 235.3121 1.84e-12
test <- wilcox.test(hsa_micro_down_sum_nn$total_sum, gga_micro_down_sum_nn$total_sum) #419.1469 136.7807 0.000193
test <- wilcox.test(hsa_avg_down_sum_nn$total_sum, gga_avg_down_sum_nn$total_sum) #590.8773 139.3824 1.02e-138
test <- wilcox.test(hsa_micro_down_sum_neuro$total_sum, hsa_avg_down_sum_neuro$total_sum) #0.0128
test <- wilcox.test(hsa_micro_down_sum_neuro$total_sum, hsa_avg_down_sum_nn$total_sum) #0.000865
test <- wilcox.test(hsa_micro_down_sum_neuro$total_sum, hsa_micro_down_sum_nn$total_sum) #0.178
test <- wilcox.test(hsa_micro_down_sum_nn$total_sum, hsa_avg_down_sum_nn$total_sum) #0.0142
test <- wilcox.test(gga_micro_down_sum_neuro$total_sum, gga_avg_down_sum_neuro$total_sum) #0.0859
test <- wilcox.test(gga_micro_down_sum_neuro$total_sum, gga_avg_down_sum_nn$total_sum) #0.292
test <- wilcox.test(gga_micro_down_sum_neuro$total_sum, gga_micro_down_sum_nn$total_sum) #0.581
test <- wilcox.test(gga_micro_down_sum_nn$total_sum, gga_avg_down_sum_nn$total_sum) #0.597

##repeatmasker overlaps
hsa_micro_up <- read.delim("bed_files/hsa_micro_up_repeatmask_overlap.txt", sep = "\t", header = FALSE)
hsa_micro_down <- read.delim("bed_files/hsa_micro_down_repeatmask_overlap.txt", sep = "\t", header = FALSE)
hsa_avg_down <- read.delim("bed_files/hsa_avg_down_repeatmask_overlap.txt", sep = "\t", header = FALSE)
hsa_avg_up <- read.delim("bed_files/hsa_avg_up_repeatmask_overlap.txt", sep = "\t", header = FALSE)
gga_micro_up <- read.delim("bed_files/gga_micro_up_repeatmask_overlap.txt", sep = "\t", header = FALSE)
gga_micro_down <- read.delim("bed_files/gga_micro_down_repeatmask_overlap.txt", sep = "\t", header = FALSE)
gga_avg_down <- read.delim("bed_files/gga_avg_down_repeatmask_overlap.txt", sep = "\t", header = FALSE)
gga_avg_up <- read.delim("bed_files/gga_avg_up_repeatmask_overlap.txt", sep = "\t", header = FALSE)

sum_by_intron <- function(df) {
  df %>%
    group_by(V4) %>%
    summarize(total_sum = sum(V16, na.rm = TRUE))
}
hsa_micro_down_sum <- sum_by_intron(hsa_micro_down)
hsa_micro_up_sum <- sum_by_intron(hsa_micro_up)
hsa_avg_down_sum <- sum_by_intron(hsa_avg_down)
hsa_avg_up_sum <- sum_by_intron(hsa_avg_up)
gga_micro_down_sum <- sum_by_intron(gga_micro_down)
gga_micro_up_sum <- sum_by_intron(gga_micro_up)
gga_avg_down_sum <- sum_by_intron(gga_avg_down)
gga_avg_up_sum <- sum_by_intron(gga_avg_up)

micro_hsa <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/detailed_list.txt", sep = "\t")
hsa_neuro_micro <- micro_hsa %>%
  filter(tissue %in% c("fdr_g", "fdr_l")) %>%
  pull(EVENT)
hsa_nn_micro <- micro_hsa %>%
  filter(tissue %in% c("common", "rare", "other")) %>%
  pull(EVENT)
avg_hsa <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/avg_regroup.txt", sep = "\t")
hsa_neuro_avg <- avg_hsa %>%
  filter(tissue %in% c("fdr_g", "fdr_l")) %>%
  pull(EVENT)
hsa_nn_avg <- avg_hsa %>%
  filter(tissue %in% c("common", "rare", "other")) %>%
  pull(EVENT)
micro_gga <- read.delim("/data2/lackey_lab/randazza/microexons/devo/chicken/gga_classes.txt", sep = "\t")
gga_neuro_micro <- micro_gga %>%
  filter(tissue %in% c("fdr_g", "fdr_l")) %>%
  pull(EVENT)
gga_nn_micro <- micro_gga %>%
  filter(tissue %in% c("common", "rare", "other")) %>%
  pull(EVENT)
avg_gga <- read.delim("/data2/lackey_lab/randazza/microexons/devo/chicken/avg_control/gga_avg_classes.txt", sep = "\t")
gga_neuro_avg <- avg_gga %>%
  filter(tissue %in% c("fdr_g", "fdr_l")) %>%
  pull(EVENT)
gga_nn_avg <- avg_gga %>%
  filter(tissue %in% c("common", "rare", "other")) %>%
  pull(EVENT)

hsa_micro_up_sum_neuro <- hsa_micro_up_sum %>%
  filter(V4 %in% hsa_neuro_micro) %>%
  mutate(species = "Hsa Neural Microexon")
hsa_micro_up_sum_nn <- hsa_micro_up_sum %>%
  filter(V4 %in% hsa_nn_micro) %>%
  mutate(species = "Hsa Non-Neural Microexon")
hsa_micro_down_sum_neuro <- hsa_micro_down_sum %>%
  filter(V4 %in% hsa_neuro_micro) %>%
  mutate(species = "Hsa Neural Microexon")
hsa_micro_down_sum_nn <- hsa_micro_down_sum %>%
  filter(V4 %in% hsa_nn_micro) %>%
  mutate(species = "Hsa Non-Neural Microexon")

hsa_avg_up_sum_neuro <- hsa_avg_up_sum %>%
  filter(V4 %in% hsa_neuro_avg) %>%
  mutate(species = "Hsa Neural Midexon")
hsa_avg_up_sum_nn <- hsa_avg_up_sum %>%
  filter(V4 %in% hsa_nn_avg) %>%
  mutate(species = "Hsa Non-Neural Midexon")
hsa_avg_down_sum_neuro <- hsa_avg_down_sum %>%
  filter(V4 %in% hsa_neuro_avg) %>%
  mutate(species = "Hsa Neural Midexon")
hsa_avg_down_sum_nn <- hsa_avg_down_sum %>%
  filter(V4 %in% hsa_nn_avg) %>%
  mutate(species = "Hsa Non-Neural Midexon")

gga_micro_up_sum_neuro <- gga_micro_up_sum %>%
  filter(V4 %in% gga_neuro_micro) %>%
  mutate(species = "Gga Neural Microexon")
gga_micro_up_sum_nn <- gga_micro_up_sum %>%
  filter(V4 %in% gga_nn_micro) %>%
  mutate(species = "Gga Non-Neural Microexon")
gga_micro_down_sum_neuro <- gga_micro_down_sum %>%
  filter(V4 %in% gga_neuro_micro) %>%
  mutate(species = "Gga Neural Microexon")
gga_micro_down_sum_nn <- gga_micro_down_sum %>%
  filter(V4 %in% gga_nn_micro) %>%
  mutate(species = "Gga Non-Neural Microexon")

gga_avg_up_sum_neuro <- gga_avg_up_sum %>%
  filter(V4 %in% gga_neuro_avg) %>%
  mutate(species = "Gga Neural Midexon")
gga_avg_up_sum_nn <- gga_avg_up_sum %>%
  filter(V4 %in% gga_nn_avg) %>%
  mutate(species = "Gga Non-Neural Midexon")
gga_avg_down_sum_neuro <- gga_avg_down_sum %>%
  filter(V4 %in% gga_neuro_avg) %>%
  mutate(species = "Gga Neural Midexon")
gga_avg_down_sum_nn <- gga_avg_down_sum %>%
  filter(V4 %in% gga_nn_avg) %>%
  mutate(species = "Gga Non-Neural Midexon")

combo_up <- rbind(hsa_micro_up_sum_neuro, hsa_avg_up_sum_neuro, gga_micro_up_sum_neuro, gga_avg_up_sum_neuro,
                  hsa_micro_up_sum_nn, hsa_avg_up_sum_nn, gga_micro_up_sum_nn, gga_avg_up_sum_nn)
combo_up$species <- factor(combo_up$species, levels = c("Hsa Neural Microexon","Gga Neural Microexon",
                                                        "Hsa Neural Midexon", "Gga Neural Midexon",
                                                        "Hsa Non-Neural Microexon","Gga Non-Neural Microexon",
                                                        "Hsa Non-Neural Midexon", "Gga Non-Neural Midexon"))
combo_down <- rbind(hsa_micro_down_sum_neuro, hsa_avg_down_sum_neuro, gga_micro_down_sum_neuro, gga_avg_down_sum_neuro,
                    hsa_micro_down_sum_nn, hsa_avg_down_sum_nn, gga_micro_down_sum_nn, gga_avg_down_sum_nn)
combo_down$species <- factor(combo_down$species, levels = c("Hsa Neural Microexon","Gga Neural Microexon",
                                                            "Hsa Neural Midexon", "Gga Neural Midexon",
                                                            "Hsa Non-Neural Microexon","Gga Non-Neural Microexon",
                                                            "Hsa Non-Neural Midexon", "Gga Non-Neural Midexon"))

combo_up %>%
  ggplot(aes(x = species, y = total_sum, fill = species)) +
  geom_boxplot() +
  #ylim(0,2000) +
  theme(legend.position = "none") +
  theme_minimal()
combo_up$log <- log10(combo_up$total_sum)
combo_up %>%
  ggplot(aes(x = species, y = log, fill = species)) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  ylim(0,5)
test <- wilcox.test(hsa_micro_up_sum_neuro$total_sum, gga_micro_up_sum_neuro$total_sum) #1008.784 347.1705 5.71e-09
test <- wilcox.test(hsa_avg_up_sum_neuro$total_sum, gga_avg_up_sum_neuro$total_sum) #2357.439 423.8791 3.46e-57
test <- wilcox.test(hsa_micro_up_sum_nn$total_sum, gga_micro_up_sum_nn$total_sum) #1380.358 470.8696 1.63e-10
test <- wilcox.test(hsa_avg_up_sum_nn$total_sum, gga_avg_up_sum_nn$total_sum) #2196.025 555.9946 0
test <- wilcox.test(hsa_micro_up_sum_neuro$total_sum, hsa_avg_up_sum_neuro$total_sum) #0.000165
test <- wilcox.test(hsa_micro_up_sum_neuro$total_sum, hsa_avg_up_sum_nn$total_sum) #1.69e-07
test <- wilcox.test(hsa_micro_up_sum_neuro$total_sum, hsa_micro_up_sum_nn$total_sum) #0.2
test <- wilcox.test(hsa_micro_up_sum_nn$total_sum, hsa_avg_up_sum_nn$total_sum) #4.41e-07
test <- wilcox.test(gga_micro_up_sum_neuro$total_sum, gga_avg_up_sum_neuro$total_sum) #0.665
test <- wilcox.test(gga_micro_up_sum_neuro$total_sum, gga_avg_up_sum_nn$total_sum) #0.0128
test <- wilcox.test(gga_micro_up_sum_neuro$total_sum, gga_micro_up_sum_nn$total_sum) #0.208
test <- wilcox.test(gga_micro_up_sum_nn$total_sum, gga_avg_up_sum_nn$total_sum) #0.307

combo_down %>%
  ggplot(aes(x = species, y = total_sum, fill = species)) +
  geom_boxplot() +
  #ylim(0,2000) +
  theme(legend.position = "none")
combo_down$log <- log10(combo_down$total_sum)
combo_down %>%
  ggplot(aes(x = species, y = log, fill = species)) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  ylim(0,5)
test <- wilcox.test(hsa_micro_down_sum_neuro$total_sum, gga_micro_down_sum_neuro$total_sum) #1105.924 394.8442 5.19e-11
test <- wilcox.test(hsa_avg_down_sum_neuro$total_sum, gga_avg_down_sum_neuro$total_sum) #2302.146 530.0851 6.59e-49
test <- wilcox.test(hsa_micro_down_sum_nn$total_sum, gga_micro_down_sum_nn$total_sum) #1738.199 340.6667 3.84e-20
test <- wilcox.test(hsa_avg_down_sum_nn$total_sum, gga_avg_down_sum_nn$total_sum) #2155.025 515.1775 0
test <- wilcox.test(hsa_micro_down_sum_neuro$total_sum, hsa_avg_down_sum_neuro$total_sum) #3.83e-05
test <- wilcox.test(hsa_micro_down_sum_neuro$total_sum, hsa_avg_down_sum_nn$total_sum) #2.19e-07
test <- wilcox.test(hsa_micro_down_sum_neuro$total_sum, hsa_micro_down_sum_nn$total_sum) #0.0249
test <- wilcox.test(hsa_micro_down_sum_nn$total_sum, hsa_avg_down_sum_nn$total_sum) #0.000967
test <- wilcox.test(gga_micro_down_sum_neuro$total_sum, gga_avg_down_sum_neuro$total_sum) #0.037
test <- wilcox.test(gga_micro_down_sum_neuro$total_sum, gga_avg_down_sum_nn$total_sum) #0.00204
test <- wilcox.test(gga_micro_down_sum_neuro$total_sum, gga_micro_down_sum_nn$total_sum) #0.855
test <- wilcox.test(gga_micro_down_sum_nn$total_sum, gga_avg_down_sum_nn$total_sum) #7.17e-05
