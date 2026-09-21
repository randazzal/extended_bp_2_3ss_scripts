library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(readxl)

##compare expression of known splicing factors
galgal6_genes <- read.delim("/data2/lackey_lab/DownloadedSequenceData/randazza/new_chick/gene_id_name.txt", sep = "\t", header = FALSE)

human_sets <- read_xlsx("/data2/lackey_lab/DownloadedSequenceData/randazza/chick_seq/human_splicing_genes.xlsx", sheet = "Sheet1")
##Splicing related genes found with GSEA gene sets and then overlapped with chicken genes
alt_splicing <- human_sets$GOBP_alt_splicing_spliceosome %>%
  na.omit()
ss_recon <- human_sets$GOBP_ss_recognition %>%
  na.omit()
mrna_cis <- human_sets$GOBP_mrna_cis_splicing %>%
  na.omit()
alt_reg <- human_sets$GOBP_reg_of_alt_splicing %>%
  na.omit()
mrna_react <- human_sets$reactome_mrna_splicing %>%
  na.omit()
reg_mrna <- human_sets$GOBP_reg_mrna_splicing_spliceosome %>%
  na.omit()
reg_rna <- human_sets$GOBP_reg_rna_splicing %>%
  na.omit()
endo <- human_sets$GOBP_endonucleolytic %>%
  na.omit()
trans <- human_sets$GOBP_transesterification %>%
  na.omit()
rna_splicing <- human_sets$GOBP_rna_splicing %>%
  na.omit()

all_splicing_genes <- unique(c(
  alt_splicing,
  ss_recon,
  endo,
  alt_reg,
  mrna_react,
  mrna_cis,
  reg_mrna,
  reg_rna,
  rna_splicing,
  trans
))

chick_splicing <- galgal6_genes %>%
  filter(tolower(V2) %in% tolower(all_splicing_genes))
names <- data.frame(
  V1 = c("ENSGALG00000064136", "ENSGALG00000066499"),
  V2 = c("NOVA1", "SRRM4")
)
chick_splicing <- rbind(chick_splicing, names)
write.table(chick_splicing, "/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/chick_splicing_genes.txt", sep = "\t")

gois <- read.delim("/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/chick_splicing_genes.txt", sep = "\t")
tpms <- read.delim("/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/all_tpms.txt", sep = "\t")
tpms <- merge(tpms, gois, by.x = "IDs", by.y = "V1")
tpms <- tpms %>%
  filter(!IDs %in% c("ENSGALG00000009850"))
missing <- gois %>%
  filter(!V1 %in% tpms$IDs)

brain_15 <- tpms[,c("IDs","V2","Female_15B1","Female_15B2","Female_15B3","Male_15B1","Male_15B2","Male_15B3","Male_15B4")]
female_cols <- grep("^Female_", names(brain_15), value = TRUE)
male_cols   <- grep("^Male_", names(brain_15), value = TRUE)

brain_15$Female_mean_15 <- rowMeans(brain_15[, female_cols], na.rm = TRUE)
brain_15$Male_mean_15   <- rowMeans(brain_15[, male_cols], na.rm = TRUE)

all_cols <- c(female_cols, male_cols)
brain_15$Overall_mean_15 <- rowMeans(brain_15[, all_cols], na.rm = TRUE)

heart_15 <- tpms[,c("IDs","V2","Female_15H1","Female_15H2","Female_15H3","Male_15H1","Male_15H2","Male_15H3","Female_15H4")]
female_cols <- grep("^Female_", names(heart_15), value = TRUE)
male_cols   <- grep("^Male_", names(heart_15), value = TRUE)

heart_15$Female_mean_15 <- rowMeans(heart_15[, female_cols], na.rm = TRUE)
heart_15$Male_mean_15   <- rowMeans(heart_15[, male_cols], na.rm = TRUE)

all_cols <- c(female_cols, male_cols)
heart_15$Overall_mean_15 <- rowMeans(heart_15[, all_cols], na.rm = TRUE)

brain_27 <- tpms[,c("IDs","V2","Female_5B1","Female_5B2","Female_5B3","Male_5B1","Male_5B2","Male_5B3")]
female_cols <- grep("^Female_", names(brain_27), value = TRUE)
male_cols   <- grep("^Male_", names(brain_27), value = TRUE)

brain_27$Female_mean_27 <- rowMeans(brain_27[, female_cols], na.rm = TRUE)
brain_27$Male_mean_27   <- rowMeans(brain_27[, male_cols], na.rm = TRUE)

all_cols <- c(female_cols, male_cols)
brain_27$Overall_mean_27 <- rowMeans(brain_27[, all_cols], na.rm = TRUE)

heart_27 <- tpms[,c("IDs","V2","Female_5H1","Female_5H2","Female_5H3","Male_5H1","Male_5H2","Male_5H3")]
female_cols <- grep("^Female_", names(heart_27), value = TRUE)
male_cols   <- grep("^Male_", names(heart_27), value = TRUE)

heart_27$Female_mean_27 <- rowMeans(heart_27[, female_cols], na.rm = TRUE)
heart_27$Male_mean_27   <- rowMeans(heart_27[, male_cols], na.rm = TRUE)

all_cols <- c(female_cols, male_cols)
heart_27$Overall_mean_27 <- rowMeans(heart_27[, all_cols], na.rm = TRUE)

brain_15 <- brain_15[,c("IDs", "V2", "Female_mean_15", "Male_mean_15", "Overall_mean_15")]
brain_27 <- brain_27[,c("IDs", "V2", "Female_mean_27", "Male_mean_27", "Overall_mean_27")]

heart_15 <- heart_15[,c("IDs", "V2", "Female_mean_15", "Male_mean_15", "Overall_mean_15")]
heart_27 <- heart_27[,c("IDs", "V2", "Female_mean_27", "Male_mean_27", "Overall_mean_27")]

combo_brain <- merge(brain_15, brain_27, by = c("IDs", "V2"))
combo_heart <- merge(heart_15, heart_27, by = c("IDs", "V2"))

combo_brain$female_score <- log(combo_brain$Female_mean_27/combo_brain$Female_mean_15)
combo_brain$male_score <- log(combo_brain$Male_mean_27/combo_brain$Male_mean_15)
combo_brain$overall_score <- log(combo_brain$Overall_mean_27/combo_brain$Overall_mean_15)

combo_heart$female_score <- log(combo_heart$Female_mean_27/combo_heart$Female_mean_15)
combo_heart$male_score <- log(combo_heart$Male_mean_27/combo_heart$Male_mean_15)
combo_heart$overall_score <- log(combo_heart$Overall_mean_27/combo_heart$Overall_mean_15)

plot_df <- combo_heart %>%
  arrange(overall_score) %>%
  mutate(cumulative_fraction = row_number() / n())

targets <- plot_df %>%
  filter(V2 %in% c("NOVA1", "SRRM4"))

ggplot(plot_df, aes(x = overall_score, y = cumulative_fraction)) +
  geom_point(color = "grey40", size = 1.8, alpha = 0.7) +
  geom_point(
    data = targets,
    aes(color = V2),
    size = 3
  ) +
  geom_text(
    data = targets,
    aes(label = V2, color = V2),
    nudge_y = 0.04,
    check_overlap = TRUE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  labs(
    x = "HH27–HH15 expression change (log2 fold change)",
    y = "Cumulative fraction of splicing factors"
  ) +
  theme_classic()

plot_df <- combo_brain %>%
  arrange(overall_score) %>%
  mutate(cumulative_fraction = row_number() / n())

targets <- plot_df %>%
  filter(V2 %in% c("NOVA1", "SRRM4"))

ggplot(plot_df, aes(x = overall_score, y = cumulative_fraction)) +
  geom_point(color = "grey40", size = 1.8, alpha = 0.7) +
  geom_point(
    data = targets,
    aes(color = V2),
    size = 3
  ) +
  geom_text(
    data = targets,
    aes(label = V2, color = V2),
    nudge_y = 0.04,
    check_overlap = TRUE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  labs(
    x = "HH27–HH15 expression change (log2 fold change)",
    y = "Cumulative fraction of splicing factors"
  ) +
  theme_classic()
