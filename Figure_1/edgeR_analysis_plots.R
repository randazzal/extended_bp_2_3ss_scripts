library(edgeR)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(ggrepel)

summary <- read.table("/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/counts.txt.summary", sep = "\t")
BioCC_input = read.table(file = "/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/counts.txt", sep = "\t", header = T)
#change row names to gene id
rownames(BioCC_input) <- BioCC_input$Geneid
#delete unnecessary columns
BioCC_input[,1:6] <- NULL
colnames(BioCC_input) <- c("x5B1", "x5B2", "x5B3", "x5B4", "x5B5", "x5B6", "x5H1", "x5H2", "x5H3", "x5H4",
                           "x5H5", "x5H6", "x7B1", "x7B2", "x7B3", "x7B4", "x7B5", "x7B6", "x7H1", "x7H3", "x7H4",
                           "x7H7", "x7H9", "x7H10", "x9B10",
                           "x9B11", "x9B13", "x9B5", "x9B6", "x9H4", "x9H5", "x9H13", "x9H14", "x9H15", "x9H9",
                           "S11B12", "S11B13", "S11B14", "S11B15",
                           "S11B1", "S11B4", "S11B19", "S11H12", "S11H13", "S11H14", "S11H1", "S11H4", "S11H5",
                           "S11H19", "S15B1", "S15B4", "S15B5", "S15B6", "S15B7", "S15B8", "S15B9", "S15H1", "S15H4",
                           "S15H5", "S15H6", "S15H7", "S15H8", "S15H9", "9B15", "9H22", "11H20")

BioCC_input <- BioCC_input[, c("x5B1", "x5B2", "x5B3", "x5B4", "x5B5", "x5B6", "x5H1", "x5H2", "x5H3", "x5H4",
                               "x5H5", "x5H6", "x7B1", "x7B2", "x7B3", "x7B4", "x7B5", "x7B6", "x7H1", "x7H3", "x7H4",
                               "x7H7", "x7H9", "x7H10", "x9B10", 
                               "x9B11", "x9B13", "x9B5", "x9B6", "9B15", "x9H4", "x9H5", "x9H13", "x9H14", "x9H15", "x9H9", "9H22",
                               "S11B12", "S11B13", "S11B14", "S11B15",
                               "S11B1", "S11B4", "S11B19", "S11H12", "S11H13", "S11H14", "S11H1", "S11H4", "S11H5",
                               "S11H19", "11H20", "S15B1", "S15B4", "S15B5", "S15B6", "S15B7", "S15B8", "S15B9", "S15H1", "S15H4",
                               "S15H5", "S15H6", "S15H7", "S15H8", "S15H9")]
#removes genes that were not expressed in at least half of samples
#if value is equal to zero then row is summed, rows greater than 3 are kept
BioCC_input_zero_filt <- BioCC_input[rowSums(BioCC_input == 0) <= 30,]
#eliminate those with low expression based on global median, for loop
Medians = numeric()
for (i in 1:nrow(BioCC_input_zero_filt)){
  Medians[i] <- median(as.numeric(BioCC_input_zero_filt[i,]))
}
BioCC_final <- BioCC_input_zero_filt[Medians>=10,]
rpk <- BioCC_final
#assign samples to groups
Groups <- as.factor(c("A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B", "C", "C", "C", "C", "C",
                      "C", "D", "D", "D", "D", "D", "D", "E", "E", "E", "E", "E", "E", "F", "F", "F", "F", "F", "F", "F", "G", "G", "G", "G", "G", "G",
                      "G", "H", "H", "H", "H", "H", "H", "H", "H", "I", "I", "I", "I", "I", "I", "I", "J", "J", "J", "J", "J",
                      "J", "J"))

#create DGEList, counts = numeric matrix of read counts, group = factor giving the conditions of each sample
rpk.norm.g <- DGEList(counts = rpk, group = Groups)
#calculate normalization factor for TMM normalization (accounts for RNA composition differences)
rpk.norm.g <- calcNormFactors(rpk.norm.g, method = "TMM")
#calculate counts per million (for further analysis)
cpm_data <- cpm(rpk.norm.g)
plotMDS(rpk.norm.g, labels=Groups)
###TPM calculations
eff.lib.size <- rpk.norm.g$samples$lib.size * rpk.norm.g$samples$norm.factors
gene_lengths <- read.delim("/data2/lackey_lab/DownloadedSequenceData/randazza/new_chick/chicken_gene_lengths.txt", sep = "\t")
gene.length.kb <- gene_lengths/1000
sum(rownames(rpk.norm.g$counts) %in% rownames(gene.length.kb))
gene_len_vec <- setNames(gene.length.kb$x, rownames(gene.length.kb))
common_genes <- intersect(names(gene_len_vec), rownames(rpk.norm.g$counts))
gene_len_aligned <- gene_len_vec[common_genes]
rpk.norm.g$counts <- rpk.norm.g$counts[common_genes, ]
rpk_data <- rpk.norm.g$counts / gene_len_aligned
tpm.tmm <- t(t(rpk_data) / (eff.lib.size / 1e6))

colnames(tpm.tmm) <- c("Male_5B1", "Female_5B1", "Female_5B2", "Male_5B2", "Female_5B3", "Male_5B3",
                       "Male_5H1", "Female_5H1", "Female_5H2", "Female_5H3", "Male_5H2", "Male_5H3",
                       "Female_7B1", "Male_7B1", "Female_7B2", "Female_7B3", "Male_7B2", "Male_7B3",
                       "Female_7H1", "Female_7H2", "Male_7H1", "Female_7H3", "Male_7H2", "Male_7H3",
                       "Female_9B1", "Male_9B1", "Male_9B2", "Female_9B2", "Male_9B3", "Female_9B3", 
                       "Female_9H1", "Female_9H2", "Female_9H3", "Male_9H1", "Male_9H2", "Female_9H4", "Male_9H3", 
                       "Male_11B1", "Male_11B2", "Female_11B1", "Female_11B2", "Male_11B3", "Male_11B4", "Female_11B3",
                       "Male_11H1", "Male_11H2", "Female_11H1", "Male_11H3", "Male_11H4", "Male_11H5", "Female_11H2", "Female_11H3",
                       "Female_15B1", "Female_15B2", "Male_15B1", "Male_15B2", "Male_15B3", "Female_15B3", "Male_15B4",
                       "Female_15H1", "Female_15H2", "Female_15H3", "Male_15H1", "Male_15H2", "Female_15H4", "Male_15H3")
tpm_data <- as.data.frame(tpm.tmm)
tpm_data$IDs <- rownames(tpm_data)
write.table(tpm_data, "/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/all_tpms.txt", sep = "\t")

gois_micro <- c("ENSGALG00000004045", "ENSGALG00000047813", "ENSGALG00000016419", "ENSGALG00000007167",
          "ENSGALG00000033512", "ENSGALG00000035156", "ENSGALG00000010967", "ENSGALG00000003969",
          "ENSGALG00000017075", "ENSGALG00000036677",
          "ENSGALG00000032366", "ENSGALG00000037423", "ENSGALG00000015511")

limited <- tpm_data %>%
  filter(IDs %in% gois_micro)
limited_long <- limited %>%
  pivot_longer(cols = c("Male_5B1", "Female_5B1", "Female_5B2", "Male_5B2", "Female_5B3", "Male_5B3",
                        "Male_5H1", "Female_5H1", "Female_5H2", "Female_5H3", "Male_5H2", "Male_5H3",
                        "Female_7B1", "Male_7B1", "Female_7B2", "Female_7B3", "Male_7B2", "Male_7B3",
                        "Female_7H1", "Female_7H2", "Male_7H1", "Female_7H3", "Male_7H2", "Male_7H3",
                        "Female_9B1", "Male_9B1", "Male_9B2", "Female_9B2", "Male_9B3",
                        "Female_9H1", "Female_9H2", "Female_9H3", "Male_9H1", "Male_9H2", "Female_9H4",
                        "Male_11B1", "Male_11B2", "Female_11B1", "Female_11B2", "Male_11B3", "Male_11B4", "Female_11B3",
                        "Male_11H1", "Male_11H2", "Female_11H1", "Male_11H3", "Male_11H4", "Male_11H5", "Female_11H2",
                        "Female_15B1", "Female_15B2", "Male_15B1", "Male_15B2", "Male_15B3", "Female_15B3", "Male_15B4",
                        "Female_15H1", "Female_15H2", "Female_15H3", "Male_15H1", "Male_15H2", "Female_15H4", "Male_15H3",
                        "Female_9B3", "Male_9H3", "Female_11H3"),
               names_to = "Sample", values_to = "TPMs")
tpm_long <- limited_long %>%
  mutate(
    Sex = str_extract(Sample, "Male|Female"),
    Stage = str_extract(Sample, "(?<=_)[0-9]+"),
    Tissue = sub(".*_[0-9]+([A-Z])[0-9]+$", "\\1", Sample),
    Replicate = str_extract(Sample, "[0-9]+$")
  )
names_micro <- data.frame(
  IDs = c("ENSGALG00000004045", "ENSGALG00000047813", "ENSGALG00000016419", "ENSGALG00000007167",
          "ENSGALG00000033512", "ENSGALG00000035156", "ENSGALG00000010967", "ENSGALG00000003969",
          "ENSGALG00000017075", "ENSGALG00000036677",
          "ENSGALG00000032366", "ENSGALG00000037423", "ENSGALG00000015511"),
  Gene = c("AGAP1", "APBB2", "ASAP2", "CLEC16A", "CPEB4", "DCTN1", "DOCK7", "EMC1", "FRY", "ITSN1",
           "MEF2A", "PTPRK", "ROBO1")
)

tpm_long <- merge(tpm_long, names_micro, 
                  by.x = "IDs", by.y = "IDs", 
                  all.x = TRUE)
tpm_summary <- tpm_long %>%
  group_by(IDs, Sex, Tissue, Stage) %>%
  summarise(mean_TPM = mean(TPMs, na.rm = TRUE), .groups = "drop")

ggplot() +
  geom_col(
    data = tpm_summary,
    aes(x = Stage, y = mean_TPM, fill = Sex),
    position = position_dodge(width = 0.8),
    color = "black"
  ) +
  geom_jitter(
    data = tpm_long,
    aes(x = Stage, y = TPMs, color = Sex),
    width = 0.15,
    size = 2,
    alpha = 0.7
  ) +
  facet_grid(Gene ~ Tissue, scales = "free_y") +
  labs(
    x = "Developmental Stage",
    y = "TPM",
    title = "Gene Expression by Sex, Tissue, and Stage"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

tpm_long %>%
  filter(Tissue == "B") %>%
  ggplot(aes(x = Stage, y = TPMs, fill = Sex)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Sex), position = position_dodge(width = 0.8), size = 1.5, alpha = 0.7) +
  facet_wrap(~Gene, scales = "free_y") +
  theme_bw() +
  labs(title = "Expression in Brain", y = "TPM", x = "Stage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tpm_long %>%
  filter(Tissue == "H") %>%
  ggplot(aes(x = Stage, y = TPMs, fill = Sex)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Sex), position = position_dodge(width = 0.8), size = 1.5, alpha = 0.7) +
  facet_wrap(~Gene) +
  theme_bw() +
  labs(title = "Expression in Heart", y = "TPM", x = "Stage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,70)

## RBP expression levels
gois_rbp <- c("ENSGALG00000064136", "ENSGALG00000066499", "ENSGALG00000012540", "ENSGALG00000032889", "ENSGALG00000034325", "ENSGALG00000011305", "ENSGALG00000009288", "ENSGALG00000001962")

limited <- tpm_data %>%
  filter(IDs %in% gois_rbp)
limited_long <- limited %>%
  pivot_longer(cols = c("Male_5B1", "Female_5B1", "Female_5B2", "Male_5B2", "Female_5B3", "Male_5B3",
                        "Male_5H1", "Female_5H1", "Female_5H2", "Female_5H3", "Male_5H2", "Male_5H3",
                        "Female_7B1", "Male_7B1", "Female_7B2", "Female_7B3", "Male_7B2", "Male_7B3",
                        "Female_7H1", "Female_7H2", "Male_7H1", "Female_7H3", "Male_7H2", "Male_7H3",
                        "Female_9B1", "Male_9B1", "Male_9B2", "Female_9B2", "Male_9B3",
                        "Female_9H1", "Female_9H2", "Female_9H3", "Male_9H1", "Male_9H2", "Female_9H4",
                        "Male_11B1", "Male_11B2", "Female_11B1", "Female_11B2", "Male_11B3", "Male_11B4", "Female_11B3",
                        "Male_11H1", "Male_11H2", "Female_11H1", "Male_11H3", "Male_11H4", "Male_11H5", "Female_11H2",
                        "Female_15B1", "Female_15B2", "Male_15B1", "Male_15B2", "Male_15B3", "Female_15B3", "Male_15B4",
                        "Female_15H1", "Female_15H2", "Female_15H3", "Male_15H1", "Male_15H2", "Female_15H4", "Male_15H3",
                        "Female_9B3", "Male_9H3", "Female_11H3"),
               names_to = "Sample", values_to = "TPMs")
tpm_long <- limited_long %>%
  mutate(
    Sex = str_extract(Sample, "Male|Female"),
    Stage = str_extract(Sample, "(?<=_)[0-9]+"),
    Tissue = sub(".*_[0-9]+([A-Z])[0-9]+$", "\\1", Sample),
    Replicate = str_extract(Sample, "[0-9]+$")
  )

names_rbp <- data.frame(
  IDs = c("ENSGALG00000064136", "ENSGALG00000066499", "ENSGALG00000012540", "ENSGALG00000032889", "ENSGALG00000034325", "ENSGALG00000011305", "ENSGALG00000009288", "ENSGALG00000001962"),
  Gene = c("NOVA1", "SRRM4", "RBFOX2", "RBFOX3", "RBFOX1", "SRSF11", "RNPS1", "PTBP1")
)

tpm_long <- merge(tpm_long, names_rbp, 
                  by.x = "IDs", by.y = "IDs", 
                  all.x = TRUE)
tpm_summary <- tpm_long %>%
  group_by(IDs, Sex, Tissue, Stage) %>%
  summarise(mean_TPM = mean(TPMs, na.rm = TRUE), .groups = "drop")

ggplot() +
  geom_col(
    data = tpm_summary,
    aes(x = Stage, y = mean_TPM, fill = Sex),
    position = position_dodge(width = 0.8),
    color = "black"
  ) +
  geom_jitter(
    data = tpm_long,
    aes(x = Stage, y = TPMs, color = Sex),
    width = 0.15,
    size = 2,
    alpha = 0.7
  ) +
  facet_grid(Gene ~ Tissue, scales = "free_y") +
  labs(
    x = "Developmental Stage",
    y = "TPM",
    title = "Gene Expression by Sex, Tissue, and Stage"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

tpm_long %>%
  filter(Tissue == "B") %>%
  ggplot(aes(x = Stage, y = TPMs, fill = Sex)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Sex), position = position_dodge(width = 0.8), size = 1.5, alpha = 0.7) +
  facet_wrap(~Gene, scales = "free_y") +
  theme_bw() +
  labs(title = "Expression in Brain", y = "TPM", x = "Stage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tpm_long %>%
  filter(Tissue == "H") %>%
  ggplot(aes(x = Stage, y = TPMs, fill = Sex)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Sex), position = position_dodge(width = 0.8), size = 1.5, alpha = 0.7) +
  facet_wrap(~Gene) +
  theme_bw() +
  labs(title = "Expression in Heart", y = "TPM", x = "Stage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,70)